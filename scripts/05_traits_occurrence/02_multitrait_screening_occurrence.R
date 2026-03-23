
##%#############################################################%##
#                                                               #
####  glmmTMB Occurrence (ALL traits, 2-way only) + AIC + plots  ####
#                                                               #
##%#############################################################%##
# 思路框架 / Workflow (CN + EN)
#
# 0) Packages / 包管理
#    - Load all required packages; install if missing.
#    - Set OMP threads for TMB if helpful.
#
# 1) Data / 数据
#    - Read myoccdata.csv
#    - Occur: 0/1
#    - UI2: factor with fixed order
#    - Random effect IDs: SS, SSBS, Best_guess_binomial
#
# 2) Trait correlation / 性状相关性
#    - Correlation matrix + p-values (Hmisc::rcorr)
#    - Professional heatmaps:
#        (i) ggcorrplot (p<0.05 blank)
#        (ii) corrplot (classic)
#    - Export: PNG + PDF (+ PPTX for ggplot-based figures)
#
# 3) FULL model (no parallel) / 全模型拟合（不并行）
#    - Fixed effects: two-way interactions only
#         UI2 * StdTmeanAnomalyRS
#       + UI2 * (all traits)
#       + StdTmeanAnomalyRS * (all traits)
#    - Random effects:
#       (1|SS) + (1|SSBS) + (1|Best_guess_binomial)
#    - Use robust optimizer controls for convergence.
#
# 4) Plotting rules (ABS probability) / 作图固定规则（绝对概率）
#    A) UI2 × Trait (ABS):
#       - Fix temperature at scenarios (0, 1, median)
#       - Fix ALL OTHER traits at MEDIAN
#    B) UI2 × Temperature (ABS):
#       - Fix ALL traits at MEDIAN
#    C) Temperature × Trait (ABS):
#       - Vary temperature (x-axis) and vary the focal trait (3 levels: 25/50/75%)
#       - Fix ALL OTHER traits at MEDIAN
#       - Show panels by UI2 (facet UI2)
#
# 5) Percent change %Δp (optional but kept) / 相对变化 %Δp（保留）
#    - Baseline for %Δp:
#       UI2 = PV, Temp = 0, ALL traits = MEDIAN (and focal trait = MEDIAN)
#    - %Δp = (p - p0)/p0 * 100
#
# 6) AIC backward selection (parallel) / AIC逐步筛选（并行）
#    - Parallelize “drop-one-term refits” using future.apply
#    - Keep hierarchical principle:
#        * interactions can be dropped
#        * main effects cannot be dropped if still used in any interaction
#
# 7) BEST model plots / 最优模型重复绘图
#    - Repeat (4A/4B/4C) using best model
#
# 8) Output format / 输出
#    - For all ggplots: print first, then save PNG + PDF + PPTX
##%#############################################################%##

rm(list = ls())

## ============================================================
## 0) Packages / 包管理（缺失则安装）
## ============================================================
pkgs <- c(
  "glmmTMB", "TMB",
  "dplyr", "tidyr", "stringr",
  "ggplot2",
  "MASS",
  "future", "future.apply",
  "Matrix",
  "Hmisc",
  "ggcorrplot",
  "corrplot", "RColorBrewer",
  "officer", "rvg"
)

to_install <- pkgs[!vapply(pkgs, requireNamespace, FUN.VALUE = logical(1), quietly = TRUE)]
if (length(to_install) > 0) install.packages(to_install, dependencies = TRUE)

suppressPackageStartupMessages({
  library(glmmTMB)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(MASS)
  library(future)
  library(future.apply)
  library(Matrix)
  library(Hmisc)
  library(ggcorrplot)
  library(corrplot)
  library(RColorBrewer)
  library(officer)
  library(rvg)
})

have_topptx <- exists("topptx")

## Threads / 线程（可选）
n_cores <- max(1, parallel::detectCores(logical = TRUE) - 1)
Sys.setenv(OMP_NUM_THREADS = n_cores)

## ============================================================
## 1) Global settings / 全局设置
## ============================================================
out_root <- "glmmTMB_occurrence_alltraits_AIC_2way_ABS_PCT_medianControls"
dir.create(out_root, showWarnings = FALSE, recursive = TRUE)

dir_corr       <- file.path(out_root, "00_TraitCorrelation"); dir.create(dir_corr, showWarnings = FALSE, recursive = TRUE)
dir_models     <- file.path(out_root, "99_Models");           dir.create(dir_models, showWarnings = FALSE, recursive = TRUE)

dir_full_plots <- file.path(out_root, "01_Plots_FULLmodel");  dir.create(dir_full_plots, showWarnings = FALSE, recursive = TRUE)
dir_best_plots <- file.path(out_root, "02_Plots_BESTmodel");  dir.create(dir_best_plots, showWarnings = FALSE, recursive = TRUE)

UI2_levels <- c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High","Urban")
pal <- c("#009E73", "#0072B2", "#E69F00", "#D55E00", "#CC79A7")
names(pal) <- UI2_levels

BASE_UI2   <- "Primary vegetation"
BASE_TEMP0 <- 0

TEMP_SCENARIOS  <- c(0, 1)   # always include
ADD_TEMP_MEDIAN <- TRUE

TITLE_SIZE <- 14
BASE_SIZE  <- 13

## ============================================================
## 2) Read data / 读取数据
## ============================================================
myoccdata <- readRDS("myoccdata.rds")
myoccdata$Occur <- as.integer(myoccdata$Occur)
myoccdata$UI2   <- factor(myoccdata$UI2, levels = UI2_levels)

stopifnot(all(c("SS","SSBS","Best_guess_binomial") %in% names(myoccdata)))
stopifnot("StdTmeanAnomalyRS" %in% names(myoccdata))

## Traits (exclude TL/DC; no synonym correction) / 性状变量（不含TL/DC，不做同物异名）
trait_vars_cont <- c(
  "RS.rs", "HB.rs", "TR.rs", "HWI.rs", "GL.rs", "CS.rs",
  "Tmin_position", "Tmax_position"
)
miss_traits <- setdiff(trait_vars_cont, names(myoccdata))
if (length(miss_traits) > 0) stop("Missing trait columns: ", paste(miss_traits, collapse = ", "))

## Trait labels / 性状标签（用于图标题）
trait_labels <- c(
  "RS.rs" = "Range size (RS.rs; standardised log RS)",
  "HB.rs" = "Habitat breadth (HB.rs; standardised sqrt HB)",
  "TR.rs" = "Thermal tolerance breadth (TR.rs; standardised TR)",
  "HWI.rs"= "Dispersal ability (HWI.rs; standardised HWI)",
  "GL.rs" = "Generation length (GL.rs; standardised log GL)",
  "CS.rs" = "Clutch size (CS.rs; standardised log CS)",
  "Tmin_position" = "Tmin_position (0–1)",
  "Tmax_position" = "Tmax_position (0–1)"
)

## Medians for control / 中位数控制值（核心：你要求所有“其它变量”用中位数固定）
med_temp   <- median(myoccdata$StdTmeanAnomalyRS, na.rm = TRUE)
trait_med  <- lapply(trait_vars_cont, function(v) median(myoccdata[[v]], na.rm = TRUE))
names(trait_med) <- trait_vars_cont

## ============================================================
## 3) Trait correlation / 性状相关性（含p值）
## ============================================================
trait_mat <- myoccdata %>%
  dplyr::select(all_of(trait_vars_cont)) %>%
  mutate(across(everything(), as.numeric)) %>%
  as.matrix()

rc <- Hmisc::rcorr(trait_mat, type = "pearson")
cor_mat <- rc$r
p_mat   <- rc$P

write.csv(cor_mat, file.path(dir_corr, "trait_correlation_r.csv"))
write.csv(p_mat,   file.path(dir_corr, "trait_correlation_p.csv"))

## helper: print then export PNG/PDF/PPTX for ggplot objects / 先看再导出
save_and_show_plot <- function(p, png_file, pdf_file = NULL, pptx_file = NULL,
                               width = 10, height = 7, dpi = 300) {
  print(p)
  ggsave(png_file, p, width = width, height = height, dpi = dpi)
  if (!is.null(pdf_file)) ggsave(pdf_file, p, width = width, height = height)
  
  if (!is.null(pptx_file)) {
    if (exists("topptx")) {
      topptx(p, pptx_file)
    } else if (requireNamespace("officer", quietly = TRUE) && requireNamespace("rvg", quietly = TRUE)) {
      doc <- officer::read_pptx()
      doc <- officer::add_slide(doc, layout = "Title and Content", master = "Office Theme")
      doc <- officer::ph_with(doc, rvg::dml(ggobj = p), location = officer::ph_location_fullsize())
      print(doc, target = pptx_file)
    } else {
      message("PPTX export skipped: ", pptx_file)
    }
  }
  invisible(TRUE)
}

## ggcorrplot (p<0.05 blank) / 美观热图（不显著留白）
p_corr <- ggcorrplot::ggcorrplot(
  cor_mat,
  method = "square",
  type = "lower",
  lab = TRUE,
  p.mat = p_mat,
  sig.level = 0.05,
  insig = "blank"
) +
  theme_classic(base_size = BASE_SIZE) +
  theme(plot.title = element_text(size = TITLE_SIZE, face = "bold")) +
  labs(
    title = "Trait correlation (Pearson r) with significance (p<0.05)\n性状相关性（Pearson r）及显著性（p<0.05）"
  )

save_and_show_plot(
  p_corr,
  png_file  = file.path(dir_corr, "Corr_ggcorrplot.png"),
  pdf_file  = file.path(dir_corr, "Corr_ggcorrplot.pdf"),
  pptx_file = file.path(dir_corr, "Corr_ggcorrplot.pptx"),
  width = 10.2, height = 7.8
)

## corrplot (classic) / 经典相关性图（先看后保存PNG/PDF）
corr_cols <- colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(200)

## view first / 先查看
corrplot::corrplot(cor_mat, method = "color", type = "lower",
                   col = corr_cols, tl.col = "black", tl.srt = 45,
                   p.mat = p_mat, sig.level = 0.05, insig = "blank",
                   addCoef.col = "black", number.cex = 0.7)

## then save / 再保存
png(file.path(dir_corr, "Corr_corrplot.png"), width = 1600, height = 1200, res = 200)
corrplot::corrplot(cor_mat, method = "color", type = "lower",
                   col = corr_cols, tl.col = "black", tl.srt = 45,
                   p.mat = p_mat, sig.level = 0.05, insig = "blank",
                   addCoef.col = "black", number.cex = 0.7)
dev.off()

pdf(file.path(dir_corr, "Corr_corrplot.pdf"), width = 10.2, height = 7.8)
corrplot::corrplot(cor_mat, method = "color", type = "lower",
                   col = corr_cols, tl.col = "black", tl.srt = 45,
                   p.mat = p_mat, sig.level = 0.05, insig = "blank",
                   addCoef.col = "black", number.cex = 0.7)
dev.off()

message("Trait correlation outputs saved in: ", dir_corr)

## ============================================================
## 4) FULL model fit (no parallel) / 拟合全模型（不并行）
## ============================================================
traits_term <- paste(trait_vars_cont, collapse = " + ")

full_formula <- as.formula(paste0(
  "Occur ~ ",
  "UI2 * StdTmeanAnomalyRS",
  " + UI2 * (", traits_term, ")",
  " + StdTmeanAnomalyRS * (", traits_term, ")",
  " + (1|SS) + (1|SSBS) + (1|Best_guess_binomial)"
))
message("\nFULL formula:\n", deparse(full_formula))

## robust control for convergence / 更稳健的优化控制（提高收敛概率）
ctrl_robust <- glmmTMBControl(
  optCtrl = list(iter.max = 1e5, eval.max = 1e5),
  optimizer = optim,
  optArgs = list(method = "BFGS")
)

m_full <- glmmTMB(
  formula = full_formula,
  family  = binomial(link = "logit"),
  data    = myoccdata,
  control = ctrl_robust
)

saveRDS(m_full, file.path(dir_models, "model_FULL_glmmTMB.rds"))

sink(file.path(dir_models, "FULL_model_summary.txt"))
cat("==== FULL model summary ====\n")
print(summary(m_full))
sink()

## ============================================================
## 5) Prediction engine (MVN draws) / 预测引擎（参数抽样近似不确定性）
## ============================================================
inv_logit <- function(x) 1/(1+exp(-x))

get_beta_V <- function(mod) {
  b <- fixef(mod)$cond
  V <- as.matrix(vcov(mod)$cond)
  list(beta = b, V = V)
}

## draw-based prediction on fixed effects / 固定效应参数抽样预测
pred_draws_glmmTMB <- function(mod, newdata, ndraw = 2000, seed = 1) {
  set.seed(seed)
  bv <- get_beta_V(mod)
  beta_hat <- bv$beta
  V <- bv$V
  
  X <- model.matrix(delete.response(terms(mod)), newdata)
  need <- colnames(X)
  
  beta_hat <- beta_hat[need]
  V <- V[need, need, drop = FALSE]
  
  B <- MASS::mvrnorm(n = ndraw, mu = beta_hat, Sigma = V)  # ndraw x p
  eta <- B %*% t(X)                                        # ndraw x n
  inv_logit(eta)
}

summ_draws <- function(mat) {
  data.frame(
    mean  = apply(mat, 2, mean),
    low95 = apply(mat, 2, quantile, probs = 0.025),
    high95= apply(mat, 2, quantile, probs = 0.975)
  )
}

pct_change_draws <- function(p, p0) {
  p0v <- as.vector(p0)
  p0v[p0v < 1e-9] <- 1e-9
  sweep(p, 1, p0v, FUN = function(x, b) (x - b)/b * 100)
}

## Save a ggplot to PNG+PDF+PPTX after viewing / 先看再导出三格式
save_plot3 <- function(p, outdir, tag, width = 12.5, height = 5.6, dpi = 300) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  save_and_show_plot(
    p,
    png_file  = file.path(outdir, paste0(tag, ".png")),
    pdf_file  = file.path(outdir, paste0(tag, ".pdf")),
    pptx_file = file.path(outdir, paste0(tag, ".pptx")),
    width = width, height = height, dpi = dpi
  )
}

## ============================================================
## 6) Plot builders with your MEDIAN control rules
## ============================================================

## 6A) UI2 × Trait (Temp fixed; other traits median) / UI2×性状（温度固定；其它性状中位数）
plot_UI2_by_trait_fixedTemp_ABS_PCT <- function(
    mod, data,
    trait_var, trait_label,
    temp_var = "StdTmeanAnomalyRS",
    temp_values = TEMP_SCENARIOS,
    add_temp_median = TRUE,
    n = 220,
    ndraw = 2000,
    seed = 1
) {
  ## temperature scenarios / 温度情景
  tv <- temp_values
  tmed <- median(data[[temp_var]], na.rm = TRUE)
  if (add_temp_median) tv <- sort(unique(c(tv, tmed)))
  
  tlabels <- sapply(tv, function(v) {
    if (isTRUE(all.equal(v, 0))) return("0")
    if (isTRUE(all.equal(v, 1))) return("1")
    if (isTRUE(all.equal(v, tmed))) return("the median")
    format(signif(v, 3), scientific = FALSE)
  }, USE.NAMES = FALSE)
  
  ## trait seq / 性状取值序列
  xseq <- seq(
    quantile(data[[trait_var]], 0.02, na.rm = TRUE),
    quantile(data[[trait_var]], 0.98, na.rm = TRUE),
    length.out = n
  )
  
  ## baseline for %Δp: PV + temp=0 + ALL traits median / %Δp基线（全中位数）
  nd0 <- data.frame(UI2 = factor(BASE_UI2, levels = levels(data$UI2)))
  nd0[[temp_var]] <- BASE_TEMP0
  for (v in trait_vars_cont) nd0[[v]] <- trait_med[[v]]
  nd0[[trait_var]] <- trait_med[[trait_var]]
  p0 <- pred_draws_glmmTMB(mod, nd0, ndraw = ndraw, seed = seed)  # draws x 1
  
  out_abs <- list()
  out_pct <- list()
  
  for (i in seq_along(tv)) {
    tmp0 <- tv[i]
    tlab <- tlabels[i]
    
    nd <- expand.grid(UI2 = levels(data$UI2), xx = xseq)
    nd$UI2 <- factor(nd$UI2, levels = levels(data$UI2))
    nd[[trait_var]] <- nd$xx
    nd$xx <- NULL
    nd[[temp_var]] <- tmp0
    
    ## other traits fixed at MEDIAN / 其它性状固定中位数
    for (v in trait_vars_cont) {
      if (!v %in% names(nd)) nd[[v]] <- trait_med[[v]]
      if (v != trait_var) nd[[v]] <- trait_med[[v]]
    }
    
    p_draw <- pred_draws_glmmTMB(mod, nd, ndraw = ndraw, seed = seed + i)
    s_abs  <- summ_draws(p_draw)
    
    p_pct  <- pct_change_draws(p_draw, p0)
    s_pct  <- summ_draws(p_pct)
    
    out_abs[[i]] <- bind_cols(nd, s_abs) %>% mutate(temp_label = tlab)
    out_pct[[i]] <- bind_cols(nd, s_pct) %>% mutate(temp_label = tlab)
  }
  
  dfA <- bind_rows(out_abs)
  dfP <- bind_rows(out_pct)
  
  p_abs <- ggplot(dfA, aes(x = .data[[trait_var]], group = UI2)) +
    geom_ribbon(aes(ymin = low95, ymax = high95, fill = UI2), alpha = 0.22, colour = NA) +
    geom_line(aes(y = mean, color = UI2), linewidth = 1.05) +
    facet_wrap(~temp_label, nrow = 1) +
    scale_fill_manual(values = pal, name = "Land-use type") +
    scale_color_manual(values = pal, name = "Land-use type") +
    theme_classic(base_size = BASE_SIZE) +
    theme(plot.title = element_text(size = TITLE_SIZE, face = "bold")) +
    labs(
      title = paste0("ABS: UI2 × ", trait_label, " (Temp fixed; other traits = median)\n",
                     "绝对概率：UI2×性状（温度固定；其它性状=中位数）"),
      x = trait_label,
      y = "Probability of occurrence"
    )
  
  p_pct <- ggplot(dfP, aes(x = .data[[trait_var]], group = UI2)) +
    geom_ribbon(aes(ymin = low95, ymax = high95, fill = UI2), alpha = 0.22, colour = NA) +
    geom_line(aes(y = mean, color = UI2), linewidth = 1.05) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5) +
    facet_wrap(~temp_label, nrow = 1) +
    scale_fill_manual(values = pal, name = "Land-use type") +
    scale_color_manual(values = pal, name = "Land-use type") +
    theme_classic(base_size = BASE_SIZE) +
    theme(plot.title = element_text(size = TITLE_SIZE, face = "bold")) +
    labs(
      title = paste0("%Δp: UI2 × ", trait_label, " (Temp fixed; baseline = PV@0 & all medians)\n",
                     "相对变化：UI2×性状（温度固定；基线=PV@0且全中位数）"),
      x = trait_label,
      y = "Change in probability of occurrence (%)"
    )
  
  list(abs_plot = p_abs, pct_plot = p_pct, abs_df = dfA, pct_df = dfP)
}

## 6B) UI2 × Temperature (traits median) / UI2×温度（性状全中位数）
plot_UI2_by_temp_ABS_PCT <- function(
    mod, data,
    temp_var = "StdTmeanAnomalyRS",
    temp_label = "Standardised mean temperature anomaly",
    n = 260,
    ndraw = 2000,
    seed = 1
) {
  ## temp sequence / 温度序列
  xseq <- seq(
    quantile(data[[temp_var]], 0.02, na.rm = TRUE),
    quantile(data[[temp_var]], 0.98, na.rm = TRUE),
    length.out = n
  )
  
  ## newdata grid / 预测网格：UI2 × temp；traits fixed at medians
  nd <- expand.grid(UI2 = levels(data$UI2), tmp = xseq)
  names(nd)[names(nd) == "tmp"] <- temp_var
  nd$UI2 <- factor(nd$UI2, levels = levels(data$UI2))
  for (v in trait_vars_cont) nd[[v]] <- trait_med[[v]]
  
  ## baseline for %Δp / %Δp基线：PV + temp=0 + traits median
  nd0 <- data.frame(UI2 = factor(BASE_UI2, levels = levels(data$UI2)))
  nd0[[temp_var]] <- BASE_TEMP0
  for (v in trait_vars_cont) nd0[[v]] <- trait_med[[v]]
  p0 <- pred_draws_glmmTMB(mod, nd0, ndraw = ndraw, seed = seed)  # draws x 1
  
  p_draw <- pred_draws_glmmTMB(mod, nd, ndraw = ndraw, seed = seed + 10)
  s_abs  <- summ_draws(p_draw)
  
  p_pct  <- pct_change_draws(p_draw, p0)
  s_pct  <- summ_draws(p_pct)
  
  dfA <- bind_cols(nd, s_abs)
  dfP <- bind_cols(nd, s_pct)
  
  p_abs <- ggplot(dfA, aes(x = .data[[temp_var]], group = UI2)) +
    geom_ribbon(aes(ymin = low95, ymax = high95, fill = UI2), alpha = 0.22, colour = NA) +
    geom_line(aes(y = mean, color = UI2), linewidth = 1.10) +
    scale_fill_manual(values = pal, name = "Land-use type") +
    scale_color_manual(values = pal, name = "Land-use type") +
    theme_classic(base_size = BASE_SIZE) +
    theme(plot.title = element_text(size = TITLE_SIZE, face = "bold")) +
    labs(
      title = paste0("ABS: UI2 × Temperature (all traits = median)\n",
                     "绝对概率：UI2×温度（所有性状=中位数）"),
      x = temp_label,
      y = "Probability of occurrence"
    )
  
  p_pct <- ggplot(dfP, aes(x = .data[[temp_var]], group = UI2)) +
    geom_ribbon(aes(ymin = low95, ymax = high95, fill = UI2), alpha = 0.22, colour = NA) +
    geom_line(aes(y = mean, color = UI2), linewidth = 1.10) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5) +
    scale_fill_manual(values = pal, name = "Land-use type") +
    scale_color_manual(values = pal, name = "Land-use type") +
    theme_classic(base_size = BASE_SIZE) +
    theme(plot.title = element_text(size = TITLE_SIZE, face = "bold")) +
    labs(
      title = paste0("%Δp: UI2 × Temperature (baseline = PV@0 & all medians)\n",
                     "相对变化：UI2×温度（基线=PV@0且全中位数）"),
      x = temp_label,
      y = "Change in probability of occurrence (%)"
    )
  
  list(abs_plot = p_abs, pct_plot = p_pct, abs_df = dfA, pct_df = dfP)
}

## 6C) Temperature × Trait (other traits median; trait levels 25/50/75; facet by UI2)
## 温度×性状（其它性状中位数；性状取3水平；按UI2分面）
plot_temp_by_trait3levels_facetUI2_ABS_PCT <- function(
    mod, data,
    trait_var, trait_label,
    temp_var = "StdTmeanAnomalyRS",
    temp_label = "Standardised mean temperature anomaly",
    n = 220,
    ndraw = 2000,
    seed = 1
) {
  ## temp sequence / 温度序列
  tseq <- seq(
    quantile(data[[temp_var]], 0.02, na.rm = TRUE),
    quantile(data[[temp_var]], 0.98, na.rm = TRUE),
    length.out = n
  )
  
  ## focal trait 3 levels / 性状三水平（25/50/75%）
  qv <- as.numeric(quantile(data[[trait_var]], probs = c(0.25, 0.5, 0.75), na.rm = TRUE, type = 7))
  lvl_lab <- c("Low (25%)", "Median (50%)", "High (75%)")
  
  ## build grid: UI2 × temp × traitLevel / 构建网格
  nd <- expand.grid(
    UI2 = levels(data$UI2),
    tmp = tseq,
    L = seq_along(qv)
  )
  names(nd)[names(nd) == "tmp"] <- temp_var
  nd$UI2 <- factor(nd$UI2, levels = levels(data$UI2))
  nd[[trait_var]] <- qv[nd$L]
  nd$TraitLevel <- factor(nd$L, levels = 1:3, labels = lvl_lab)
  nd$L <- NULL
  
  ## other traits fixed at median / 其它性状固定中位数
  for (v in trait_vars_cont) {
    if (!v %in% names(nd)) nd[[v]] <- trait_med[[v]]
    if (v != trait_var) nd[[v]] <- trait_med[[v]]
  }
  
  ## baseline for %Δp: PV + temp=0 + ALL traits median (focal trait median) / %Δp基线
  nd0 <- data.frame(UI2 = factor(BASE_UI2, levels = levels(data$UI2)))
  nd0[[temp_var]] <- BASE_TEMP0
  for (v in trait_vars_cont) nd0[[v]] <- trait_med[[v]]
  nd0[[trait_var]] <- trait_med[[trait_var]]
  p0 <- pred_draws_glmmTMB(mod, nd0, ndraw = ndraw, seed = seed)
  
  ## predict draws / 预测
  p_draw <- pred_draws_glmmTMB(mod, nd, ndraw = ndraw, seed = seed + 20)
  s_abs  <- summ_draws(p_draw)
  
  p_pct  <- pct_change_draws(p_draw, p0)
  s_pct  <- summ_draws(p_pct)
  
  dfA <- bind_cols(nd, s_abs)
  dfP <- bind_cols(nd, s_pct)
  
  p_abs <- ggplot(dfA, aes(x = .data[[temp_var]], group = TraitLevel)) +
    geom_ribbon(aes(ymin = low95, ymax = high95, fill = TraitLevel), alpha = 0.18, colour = NA) +
    geom_line(aes(y = mean, color = TraitLevel), linewidth = 1.05) +
    facet_wrap(~UI2, nrow = 1) +
    theme_classic(base_size = BASE_SIZE) +
    theme(plot.title = element_text(size = TITLE_SIZE, face = "bold"),
          legend.position = "right") +
    labs(
      title = paste0("ABS: Temperature × ", trait_label, " (other traits = median; facet by UI2)\n",
                     "绝对概率：温度×性状（其它性状=中位数；按UI2分面）"),
      x = temp_label,
      y = "Probability of occurrence",
      color = trait_label,
      fill  = trait_label
    )
  
  p_pct <- ggplot(dfP, aes(x = .data[[temp_var]], group = TraitLevel)) +
    geom_ribbon(aes(ymin = low95, ymax = high95, fill = TraitLevel), alpha = 0.18, colour = NA) +
    geom_line(aes(y = mean, color = TraitLevel), linewidth = 1.05) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5) +
    facet_wrap(~UI2, nrow = 1) +
    theme_classic(base_size = BASE_SIZE) +
    theme(plot.title = element_text(size = TITLE_SIZE, face = "bold"),
          legend.position = "right") +
    labs(
      title = paste0("%Δp: Temperature × ", trait_label, " (baseline = PV@0 & all medians; facet by UI2)\n",
                     "相对变化：温度×性状（基线=PV@0且全中位数；按UI2分面）"),
      x = temp_label,
      y = "Change in probability of occurrence (%)",
      color = trait_label,
      fill  = trait_label
    )
  
  list(abs_plot = p_abs, pct_plot = p_pct, abs_df = dfA, pct_df = dfP)
}

## ============================================================
## 7) FULL model plotting (A/B/C) / 全模型三类两重交互作图
## ============================================================

## 7B) UI2 × Temperature (traits median)
res_full_UI2xTemp <- plot_UI2_by_temp_ABS_PCT(
  mod = m_full, data = myoccdata,
  temp_var = "StdTmeanAnomalyRS",
  temp_label = "Standardised mean temperature anomaly",
  n = 260, ndraw = 2000, seed = 1000
)
save_plot3(res_full_UI2xTemp$abs_plot, dir_full_plots, "FULL_UI2xTemp_ABS", width = 11.5, height = 5.6)
save_plot3(res_full_UI2xTemp$pct_plot, dir_full_plots, "FULL_UI2xTemp_PCT", width = 11.5, height = 5.6)

## 7A + 7C) For each trait: UI2×Trait (temp fixed) and Temp×Trait (3 levels; facet UI2)
all_full_trait_plots <- list()

for (v in trait_vars_cont) {
  message("FULL model plotting trait: ", v)
  
  ## UI2 × Trait (temp fixed; other traits median)
  r1 <- plot_UI2_by_trait_fixedTemp_ABS_PCT(
    mod = m_full, data = myoccdata,
    trait_var = v, trait_label = trait_labels[[v]],
    n = 220, ndraw = 2000, seed = 1100
  )
  save_plot3(r1$abs_plot, dir_full_plots, paste0("FULL_UI2x", v, "_TempFixed_ABS"), width = 12.8, height = 5.6)
  save_plot3(r1$pct_plot, dir_full_plots, paste0("FULL_UI2x", v, "_TempFixed_PCT"), width = 12.8, height = 5.6)
  
  ## Temp × Trait (3 trait levels; other traits median; facet UI2)
  r2 <- plot_temp_by_trait3levels_facetUI2_ABS_PCT(
    mod = m_full, data = myoccdata,
    trait_var = v, trait_label = trait_labels[[v]],
    temp_var = "StdTmeanAnomalyRS",
    temp_label = "Standardised mean temperature anomaly",
    n = 220, ndraw = 2000, seed = 1200
  )
  save_plot3(r2$abs_plot, dir_full_plots, paste0("FULL_TempX", v, "_Trait3levels_FacetUI2_ABS"), width = 14.0, height = 5.2)
  save_plot3(r2$pct_plot, dir_full_plots, paste0("FULL_TempX", v, "_Trait3levels_FacetUI2_PCT"), width = 14.0, height = 5.2)
  
  all_full_trait_plots[[v]] <- list(UI2xTrait = r1, TempXTrait = r2)
}

message("DONE: FULL model plots saved in: ", dir_full_plots)

## ============================================================
## 8) AIC backward selection (parallel, A2方案: worker readRDS)
##    AIC逐步筛选（并行；worker 内部读取模型RDS；AIC阈值=2或5）
## ============================================================

## ---- 推荐：避免过度并行（多进程时让每个进程单线程）----
Sys.setenv(OMP_NUM_THREADS = 1)

## ---- 0) 保存 full model 到磁盘（只保存一次）----
dir.create(dir_models, showWarnings = FALSE, recursive = TRUE)

m_full_path <- file.path(dir_models, "AIC_m_full.rds")
saveRDS(m_full, m_full_path)

## 当前迭代模型的“工作文件”（每步更新后覆盖写入）
m_curr_path <- file.path(dir_models, "AIC_m_curr.rds")
saveRDS(m_full, m_curr_path)

## ---- 1) 工具函数 / utilities ----
get_fixed_terms <- function(mod) {
  tt <- terms(mod)
  attr(tt, "term.labels")
}

has_interaction <- function(term) grepl(":", term)

## hierarchical principle / 层级原则
## - 允许删除交互项
## - 删除主效应时，要求当前模型中不包含任何涉及该主效应的交互项
hierarchical_ok_to_drop <- function(term, current_terms) {
  if (has_interaction(term)) return(TRUE)
  involved <- any(grepl(paste0("(^|:)", term, "(:|$)"), current_terms) & has_interaction(current_terms))
  !involved
}

## ---- 2) worker里执行：读取模型 -> update drop term -> AIC ----
## IMPORTANT:
## - 这个函数不要引用外部大对象（避免FUN巨大）
## - 只靠参数输入：m_path、term、ctrl_robust
fit_drop_one_term_from_rds <- function(trm, m_path, ctrl_robust) {
  mod <- readRDS(m_path)
  
  m2 <- try(
    update(mod, as.formula(paste(". ~ . -", trm)), control = ctrl_robust),
    silent = TRUE
  )
  
  if (inherits(m2, "try-error")) {
    return(data.frame(term = trm, AIC = NA_real_))
  } else {
    return(data.frame(term = trm, AIC = AIC(m2)))
  }
}

## ---- 3) 并行 drop1 AIC table（A2：workers读m_curr_path） ----
drop1_AIC_table_parallel_A2 <- function(m_path, workers = n_cores, verbose = FALSE) {
  
  ## 主进程读取一次，拿到terms和AIC0
  mod0 <- readRDS(m_path)
  current_terms <- get_fixed_terms(mod0)
  aic0 <- AIC(mod0)
  
  cand <- current_terms[sapply(current_terms, hierarchical_ok_to_drop, current_terms = current_terms)]
  if (length(cand) == 0) return(NULL)
  
  ## plan 控制（自包含）
  old_plan <- future::plan()
  future::plan(future::multisession, workers = workers)
  on.exit(future::plan(old_plan), add = TRUE)
  
  ## 并行评估所有候选drop
  res_list <- future.apply::future_lapply(
    cand,
    FUN = fit_drop_one_term_from_rds,
    m_path = m_path,
    ctrl_robust = ctrl_robust
  )
  
  tab <- dplyr::bind_rows(res_list)
  tab$deltaAIC <- tab$AIC - aic0
  tab <- tab %>% dplyr::arrange(deltaAIC)
  
  if (isTRUE(verbose)) {
    message("drop1 candidates = ", length(cand), " | current AIC = ", round(aic0, 2))
  }
  
  tab
}

## ---- 4) Stepwise backward (A2) with AIC threshold ----
## aic_threshold: 2 or 5
## - accept drop only if deltaAIC <= -aic_threshold
stepwise_AIC_backward_parallel_A2 <- function(
    m_start_path,
    m_work_path,
    max_steps = 200,
    workers = n_cores,
    verbose = TRUE,
    aic_threshold = 2,   # <<< 设为2或5
    save_tables = TRUE
) {
  
  stopifnot(is.numeric(aic_threshold), length(aic_threshold) == 1, aic_threshold > 0)
  
  ## 初始化工作模型文件
  file.copy(m_start_path, m_work_path, overwrite = TRUE)
  
  ## 可选：保存每步drop1表
  tab_dir <- file.path(dir_models, "AIC_drop1_tables")
  if (isTRUE(save_tables)) dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)
  
  for (k in 1:max_steps) {
    
    tab <- drop1_AIC_table_parallel_A2(m_path = m_work_path, workers = workers)
    
    if (isTRUE(save_tables) && !is.null(tab)) {
      write.csv(tab, file.path(tab_dir, paste0("drop1_step_", k, ".csv")), row.names = FALSE)
    }
    
    if (is.null(tab) || all(!is.finite(tab$AIC))) {
      if (verbose) message("No valid candidate models. Stop.")
      break
    }
    
    best <- tab[1, , drop = FALSE]
    
    ## 当前模型 AIC（主进程读取一次）
    m_curr <- readRDS(m_work_path)
    aic_curr <- AIC(m_curr)
    
    if (verbose) {
      message("\n[Step ", k, "] Current AIC = ", round(aic_curr, 2))
      message("  Best drop: ", best$term,
              "  -> AIC = ", round(best$AIC, 2),
              "  (ΔAIC = ", round(best$deltaAIC, 2), ")",
              " | threshold = ", aic_threshold)
    }
    
    ## 接受规则：AIC至少下降aic_threshold
    if (is.finite(best$deltaAIC) && best$deltaAIC <= -aic_threshold) {
      
      m_new <- update(
        m_curr,
        as.formula(paste(". ~ . -", best$term)),
        control = ctrl_robust
      )
      
      ## 保存每一步模型（可追溯）
      step_path <- file.path(dir_models, paste0("AIC_step_", k, "_model.rds"))
      saveRDS(m_new, step_path)
      
      ## 覆盖写回工作模型文件（下一轮workers就读这个）
      saveRDS(m_new, m_work_path)
      
    } else {
      if (verbose) message("  No further AIC improvement beyond threshold. Stop.")
      break
    }
  }
  
  ## 返回最终模型对象
  readRDS(m_work_path)
}

## ---- 5) RUN selection / 执行AIC筛选 ----
## 你要阈值=2：aic_threshold = 2
## 你要阈值=5：aic_threshold = 5
AIC_THRESHOLD <- 5   # <<< 改成 5 即可

m_best <- stepwise_AIC_backward_parallel_A2(
  m_start_path = m_full_path,
  m_work_path  = m_curr_path,
  max_steps = 200,
  workers = n_cores,
  verbose = TRUE,
  aic_threshold = AIC_THRESHOLD,
  save_tables = TRUE
)

## 保存最优模型
saveRDS(m_best, file.path(dir_models, paste0("model_BEST_AIC_glmmTMB_thr", AIC_THRESHOLD, ".rds")))

## 输出摘要
sink(file.path(dir_models, paste0("BEST_model_summary_thr", AIC_THRESHOLD, ".txt")))
cat("==== BEST (AIC) model summary ====\n")
cat("AIC improvement threshold =", AIC_THRESHOLD, "\n\n")
print(summary(m_best))
sink()

## export fixed effects table / 导出固定效应表
fixef_tab <- as.data.frame(summary(m_best)$coefficients$cond)
fixef_tab$term <- rownames(fixef_tab)
fixef_tab <- fixef_tab %>% dplyr::relocate(term)
sig_tab <- fixef_tab %>% dplyr::filter(`Pr(>|z|)` < 0.05)

write.csv(fixef_tab, file.path(dir_models, paste0("bestModel_fixedEffects_all_thr", AIC_THRESHOLD, ".csv")), row.names = FALSE)
write.csv(sig_tab,   file.path(dir_models, paste0("bestModel_fixedEffects_significant_p<0.05_thr", AIC_THRESHOLD, ".csv")), row.names = FALSE)

message("AIC selection done (threshold=", AIC_THRESHOLD, "). BEST model saved in: ", dir_models)


## ============================================================
## 9) BEST model plotting (A/B/C) / 最优模型重复三类作图
##    （此处逻辑不变：继续用 m_best）
## ============================================================

## 9B) UI2 × Temperature (traits median)
res_best_UI2xTemp <- plot_UI2_by_temp_ABS_PCT(
  mod = m_best, data = myoccdata,
  temp_var = "StdTmeanAnomalyRS",
  temp_label = "Standardised mean temperature anomaly",
  n = 260, ndraw = 2000, seed = 2000
)
save_plot3(res_best_UI2xTemp$abs_plot, dir_best_plots, paste0("BEST_UI2xTemp_ABS_thr", AIC_THRESHOLD), width = 11.5, height = 5.6)
save_plot3(res_best_UI2xTemp$pct_plot, dir_best_plots, paste0("BEST_UI2xTemp_PCT_thr", AIC_THRESHOLD), width = 11.5, height = 5.6)

## 9A + 9C) per trait
all_best_trait_plots <- list()

for (v in trait_vars_cont) {
  message("BEST model plotting trait: ", v)
  
  r1 <- plot_UI2_by_trait_fixedTemp_ABS_PCT(
    mod = m_best, data = myoccdata,
    trait_var = v, trait_label = trait_labels[[v]],
    n = 220, ndraw = 2000, seed = 2100
  )
  save_plot3(r1$abs_plot, dir_best_plots, paste0("BEST_UI2x", v, "_TempFixed_ABS_thr", AIC_THRESHOLD), width = 12.8, height = 5.6)
  save_plot3(r1$pct_plot, dir_best_plots, paste0("BEST_UI2x", v, "_TempFixed_PCT_thr", AIC_THRESHOLD), width = 12.8, height = 5.6)
  
  r2 <- plot_temp_by_trait3levels_facetUI2_ABS_PCT(
    mod = m_best, data = myoccdata,
    trait_var = v, trait_label = trait_labels[[v]],
    temp_var = "StdTmeanAnomalyRS",
    temp_label = "Standardised mean temperature anomaly",
    n = 220, ndraw = 2000, seed = 2200
  )
  save_plot3(r2$abs_plot, dir_best_plots, paste0("BEST_TempX", v, "_Trait3levels_FacetUI2_ABS_thr", AIC_THRESHOLD), width = 14.0, height = 5.2)
  save_plot3(r2$pct_plot, dir_best_plots, paste0("BEST_TempX", v, "_Trait3levels_FacetUI2_PCT_thr", AIC_THRESHOLD), width = 14.0, height = 5.2)
  
  all_best_trait_plots[[v]] <- list(UI2xTrait = r1, TempXTrait = r2)
}

message("\nALL DONE ✅")
message("Trait correlation: ", dir_corr)
message("FULL plots:       ", dir_full_plots)
message("BEST plots:       ", dir_best_plots)
message("Models:           ", dir_models)

## ============================================================
## Build FULL model with 3-way interactions
## 在 m_full 的基础上：加入 UI2 × 温度 × 每个性状 的三重交互
## ============================================================

## 你前面已经定义过（示例）：
## trait_vars_cont <- c("RS.rs","HB.rs","TR.rs","HWI.rs","GL.rs","CS.rs","Tmin_position","Tmax_position")

stopifnot(all(trait_vars_cont %in% names(myoccdata)))
stopifnot(all(c("UI2","StdTmeanAnomalyRS","Occur","SS","SSBS","Best_guess_binomial") %in% names(myoccdata)))

## 1) 固定效应：UI2 * Temp * (all traits)
trait_block <- paste(trait_vars_cont, collapse = " + ")

## 2) 随机效应：按你之前的设置（可按需要调整）
rand_block <- "(1|SS) + (1|SSBS) + (1|Best_guess_binomial)"

## 3) 三重交互 full formula
full_formula_3way <- as.formula(
  paste0("Occur ~ UI2 * StdTmeanAnomalyRS * (", trait_block, ") + ", rand_block)
)

cat("\n===== FULL formula (3-way) =====\n")
print(full_formula_3way)

## 4) 拟合（更稳健的控制参数用你已有 ctrl_robust）
##    注意：你说“m_full基础上添加”，可以选择 update 或直接重拟合
##    - update 更方便，但如果公式变化很大，重拟合更干净
m_full_3way <- glmmTMB(
  formula = full_formula_3way,
  family  = binomial(link = "logit"),
  data    = myoccdata,
  control = ctrl_robust
)

saveRDS(m_full_3way, file.path(dir_models, "model_FULL_3way_glmmTMB.rds"))

sink(file.path(dir_models, "FULL_3way_model_summary.txt"))
cat("==== FULL 3-way model summary ====\n")
print(summary(m_full_3way))
sink()

message("DONE: m_full_3way fitted and saved.")

m_full_3way <- update(m_full, formula = full_formula_3way, control = ctrl_robust)

## ============================================================
## 8) AIC backward selection (parallel, A2方案: worker readRDS)
##    UPDATED RULES / 更新规则：
##    - 可删除：所有交互项 + 所有性状主效应（trait main effects）
##    - 不可删除：UI2 主效应、温度主效应(StdTmeanAnomalyRS)
##    - AIC阈值：2（只有ΔAIC <= -2 才接受）
##    - 并行：worker 内部 readRDS(m_path)，避免传大对象
## ============================================================

Sys.setenv(OMP_NUM_THREADS = 1)

## ---------- 8.0 Save full model and working model ----------
dir.create(dir_models, showWarnings = FALSE, recursive = TRUE)

m_full_path <- file.path(dir_models, "AIC_m_full.rds")
saveRDS(m_full, m_full_path)

m_curr_path <- file.path(dir_models, "AIC_m_curr.rds")
saveRDS(m_full, m_curr_path)

AIC_THRESHOLD <- 2

## ---------- 8.1 Utilities ----------
has_interaction <- function(term) grepl(":", term)

get_fixed_terms <- function(mod) {
  tt <- terms(mod)
  attr(tt, "term.labels")
}

## Identify protected terms (never drop) / 保护项：永不删除
## NOTE: these are main effects names as they appear in term.labels
PROTECT_MAIN <- c("UI2", "StdTmeanAnomalyRS")

## Identify trait main effects allowed to drop / 可删的性状主效应
## IMPORTANT: you must define trait_vars_cont earlier in your script.
## Example: trait_vars_cont <- c("RS.rs","HB.rs","TR.rs","HWI.rs","GL.rs","CS.rs","Tmin_position","Tmax_position")
is_trait_main_effect <- function(term, trait_vars) {
  term %in% trait_vars
}

## Decide which terms are eligible to drop / 判定可删项
## Rules:
##  - if term is protected main effect -> NOT droppable
##  - else if term is interaction -> droppable
##  - else if term is trait main effect -> droppable
##  - else -> NOT droppable (conservative)
ok_to_drop_term <- function(term, trait_vars = trait_vars_cont) {
  if (term %in% PROTECT_MAIN) return(FALSE)
  if (has_interaction(term)) return(TRUE)
  if (is_trait_main_effect(term, trait_vars)) return(TRUE)
  FALSE
}

## ---------- 8.2 Worker function (readRDS -> update -> AIC) ----------
fit_drop_one_term_from_rds <- function(trm, m_path, ctrl_robust) {
  mod <- readRDS(m_path)
  
  m2 <- try(
    update(mod, as.formula(paste(". ~ . -", trm)), control = ctrl_robust),
    silent = TRUE
  )
  
  if (inherits(m2, "try-error")) {
    return(data.frame(term = trm, AIC = NA_real_))
  } else {
    return(data.frame(term = trm, AIC = AIC(m2)))
  }
}

## ---------- 8.3 Drop1 AIC table (parallel; A2) ----------
drop1_AIC_table_parallel_A2_rule <- function(m_path, workers = n_cores, verbose = FALSE) {
  
  ## master reads current model once
  mod0 <- readRDS(m_path)
  current_terms <- get_fixed_terms(mod0)
  aic0 <- AIC(mod0)
  
  cand <- current_terms[sapply(current_terms, ok_to_drop_term)]
  if (length(cand) == 0) return(NULL)
  
  old_plan <- future::plan()
  future::plan(future::multisession, workers = workers)
  on.exit(future::plan(old_plan), add = TRUE)
  
  res_list <- future.apply::future_lapply(
    cand,
    FUN = fit_drop_one_term_from_rds,
    m_path = m_path,
    ctrl_robust = ctrl_robust
  )
  
  tab <- dplyr::bind_rows(res_list)
  tab$deltaAIC <- tab$AIC - aic0
  tab <- tab %>% dplyr::arrange(deltaAIC)
  
  if (isTRUE(verbose)) {
    message("Candidates = ", length(cand), " | AIC0 = ", round(aic0, 2))
  }
  
  tab
}

## ---------- 8.4 Stepwise backward (parallel; A2) ----------
stepwise_AIC_backward_parallel_A2_rule <- function(
    m_start_path,
    m_work_path,
    max_steps = 200,
    workers = n_cores,
    verbose = TRUE,
    aic_threshold = 2,
    save_tables = TRUE
) {
  
  stopifnot(is.numeric(aic_threshold), length(aic_threshold) == 1, aic_threshold > 0)
  
  file.copy(m_start_path, m_work_path, overwrite = TRUE)
  
  tab_dir <- file.path(dir_models, "AIC_drop1_tables")
  if (isTRUE(save_tables)) dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)
  
  for (k in 1:max_steps) {
    
    tab <- drop1_AIC_table_parallel_A2_rule(m_path = m_work_path, workers = workers)
    
    if (isTRUE(save_tables) && !is.null(tab)) {
      write.csv(tab, file.path(tab_dir, paste0("drop1_step_", k, ".csv")), row.names = FALSE)
    }
    
    if (is.null(tab) || all(!is.finite(tab$AIC))) {
      if (verbose) message("No valid candidate models. Stop.")
      break
    }
    
    best <- tab[1, , drop = FALSE]
    
    ## current AIC
    m_curr <- readRDS(m_work_path)
    aic_curr <- AIC(m_curr)
    
    if (verbose) {
      message("\n[Step ", k, "] Current AIC = ", round(aic_curr, 2))
      message("  Best drop: ", best$term,
              " -> AIC = ", round(best$AIC, 2),
              " (ΔAIC = ", round(best$deltaAIC, 2), ")",
              " | threshold = ", aic_threshold)
    }
    
    ## accept only if AIC drops by at least threshold
    if (is.finite(best$deltaAIC) && best$deltaAIC <= -aic_threshold) {
      
      m_new <- update(
        m_curr,
        as.formula(paste(". ~ . -", best$term)),
        control = ctrl_robust
      )
      
      step_path <- file.path(dir_models, paste0("AIC_step_", k, "_model.rds"))
      saveRDS(m_new, step_path)
      
      saveRDS(m_new, m_work_path)
      
    } else {
      if (verbose) message("  No further AIC improvement beyond threshold. Stop.")
      break
    }
  }
  
  readRDS(m_work_path)
}

## ---------- 8.5 RUN selection ----------
m_best <- stepwise_AIC_backward_parallel_A2_rule(
  m_start_path = m_full_path,
  m_work_path  = m_curr_path,
  max_steps = 200,
  workers = n_cores,
  verbose = TRUE,
  aic_threshold = AIC_THRESHOLD,
  save_tables = TRUE
)

saveRDS(m_best, file.path(dir_models, "model_BEST_AIC_glmmTMB_thr2_relaxedDrop.rds"))

sink(file.path(dir_models, "BEST_model_summary_thr2_relaxedDrop.txt"))
cat("==== BEST (AIC) model summary ====\n")
cat("Rule: drop interactions + trait main effects; keep UI2 & StdTmeanAnomalyRS main effects\n")
cat("AIC threshold = ", AIC_THRESHOLD, "\n\n")
print(summary(m_best))
sink()

## export fixed effects table / 导出固定效应表
fixef_tab <- as.data.frame(summary(m_best)$coefficients$cond)
fixef_tab$term <- rownames(fixef_tab)
fixef_tab <- fixef_tab %>% dplyr::relocate(term)
sig_tab <- fixef_tab %>% dplyr::filter(`Pr(>|z|)` < 0.05)

write.csv(fixef_tab, file.path(dir_models, "bestModel_fixedEffects_all_thr2_relaxedDrop.csv"), row.names = FALSE)
write.csv(sig_tab,   file.path(dir_models, "bestModel_fixedEffects_significant_p<0.05_thr2_relaxedDrop.csv"), row.names = FALSE)

message("AIC selection done (threshold=2, relaxed drops). BEST model saved in: ", dir_models)


## ============================================================
## 9) BEST model plotting (A/B/C) / 最优模型重复三类作图
##    （你原来的作图代码可直接沿用）
## ============================================================

## 9B) UI2 × Temperature (traits median)
res_best_UI2xTemp <- plot_UI2_by_temp_ABS_PCT(
  mod = m_best, data = myoccdata,
  temp_var = "StdTmeanAnomalyRS",
  temp_label = "Standardised mean temperature anomaly",
  n = 260, ndraw = 2000, seed = 2000
)
save_plot3(res_best_UI2xTemp$abs_plot, dir_best_plots, "BEST_UI2xTemp_ABS_thr2_relaxedDrop", width = 11.5, height = 5.6)
save_plot3(res_best_UI2xTemp$pct_plot, dir_best_plots, "BEST_UI2xTemp_PCT_thr2_relaxedDrop", width = 11.5, height = 5.6)

## 9A + 9C) per trait
all_best_trait_plots <- list()

for (v in trait_vars_cont) {
  message("BEST model plotting trait: ", v)
  
  r1 <- plot_UI2_by_trait_fixedTemp_ABS_PCT(
    mod = m_best, data = myoccdata,
    trait_var = v, trait_label = trait_labels[[v]],
    n = 220, ndraw = 2000, seed = 2100
  )
  save_plot3(r1$abs_plot, dir_best_plots, paste0("BEST_UI2x", v, "_TempFixed_ABS_thr2_relaxedDrop"), width = 12.8, height = 5.6)
  save_plot3(r1$pct_plot, dir_best_plots, paste0("BEST_UI2x", v, "_TempFixed_PCT_thr2_relaxedDrop"), width = 12.8, height = 5.6)
  
  r2 <- plot_temp_by_trait3levels_facetUI2_ABS_PCT(
    mod = m_best, data = myoccdata,
    trait_var = v, trait_label = trait_labels[[v]],
    temp_var = "StdTmeanAnomalyRS",
    temp_label = "Standardised mean temperature anomaly",
    n = 220, ndraw = 2000, seed = 2200
  )
  save_plot3(r2$abs_plot, dir_best_plots, paste0("BEST_TempX", v, "_Trait3levels_FacetUI2_ABS_thr2_relaxedDrop"), width = 14.0, height = 5.2)
  save_plot3(r2$pct_plot, dir_best_plots, paste0("BEST_TempX", v, "_Trait3levels_FacetUI2_PCT_thr2_relaxedDrop"), width = 14.0, height = 5.2)
  
  all_best_trait_plots[[v]] <- list(UI2xTrait = r1, TempXTrait = r2)
}

message("\nALL DONE ✅")
message("Trait correlation: ", dir_corr)
message("FULL plots:       ", dir_full_plots)
message("BEST plots:       ", dir_best_plots)
message("Models:           ", dir_models)





## ============================================================
## Clean trait display names for legends / 性状图例显示名
## ============================================================
trait_display_name <- c(
  RS.rs = "Range size",
  HB.rs = "Habitat breadth",
  TR.rs = "Thermal tolerance breadth",
  HWI.rs = "Dispersal ability",
  GL.rs = "Generation length",
  CS.rs = "Clutch size",
  Tmin_position = "Thermal minimum position",
  Tmax_position = "Thermal maximum position"
)
plot_temp_by_trait_with_legend <- function(
    mod, data,
    trait_var,
    temp_var = "StdTmeanAnomalyRS",
    n = 220,
    ndraw = 2000,
    seed = 1,
    
    ## quantiles for legend
    qs = c(0.10, 0.50, 0.90),
    
    ## plotting params (可手调)
    base_size = 13,
    title_size = 13,
    legend_title_size = 13,
    legend_text_size = 12,
    line_width = 1.2,
    ribbon_alpha = 0.22
) {
  
  set.seed(seed)
  
  ## ---- quantile values of the focal trait ----
  q_vals <- quantile(data[[trait_var]], probs = qs, na.rm = TRUE, type = 7)
  names(q_vals) <- c("Low (10%)", "Median (50%)", "High (90%)")
  
  ## ---- temperature sequence ----
  temp_seq <- seq(
    quantile(data[[temp_var]], 0.02, na.rm = TRUE),
    quantile(data[[temp_var]], 0.98, na.rm = TRUE),
    length.out = n
  )
  
  ## ---- baseline (for ABS only, still needed for predictions) ----
  nd <- expand.grid(
    UI2 = levels(data$UI2),
    temp = temp_seq,
    Q = names(q_vals)
  )
  nd$UI2 <- factor(nd$UI2, levels = levels(data$UI2))
  nd[[temp_var]] <- nd$temp
  nd$temp <- NULL
  
  ## ---- set focal trait by quantile, others at median ----
  for (v in trait_vars_cont) {
    if (v == trait_var) {
      nd[[v]] <- q_vals[nd$Q]
    } else {
      nd[[v]] <- median(data[[v]], na.rm = TRUE)
    }
  }
  
  ## ---- predict (fixed effects only) ----
  p_draw <- pred_draws_glmmTMB(
    mod = mod,
    newdata = nd,
    ndraw = ndraw,
    seed = seed
  )
  
  summ <- summ_draws(p_draw)
  df <- cbind(nd, summ)
  
  ## ---- nicer labels for legend ----
  df$Q <- factor(
    df$Q,
    levels = c("Low (10%)", "Median (50%)", "High (90%)")
  )
  
  ## ---- plot ----
  p <- ggplot(df, aes(x = .data[[temp_var]], colour = Q, fill = Q)) +
    geom_ribbon(
      aes(ymin = low95, ymax = high95),
      alpha = ribbon_alpha,
      colour = NA
    ) +
    geom_line(
      aes(y = mean),
      linewidth = line_width
    ) +
    facet_wrap(~UI2, nrow = 1) +
    scale_colour_manual(
      values = c("Low (10%)" = "#E41A1C",
                 "Median (50%)" = "#4DAF4A",
                 "High (90%)" = "#377EB8"),
      name = trait_display_name[[trait_var]]
    ) +
    scale_fill_manual(
      values = c("Low (10%)" = "#E41A1C",
                 "Median (50%)" = "#4DAF4A",
                 "High (90%)" = "#377EB8"),
      name = trait_display_name[[trait_var]]
    ) +
    theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(size = title_size, face = "bold"),
      legend.title = element_text(size = legend_title_size),
      legend.text  = element_text(size = legend_text_size),
      legend.position = "top"
    ) +
    labs(
      title = paste0(
        "ABS: Temperature anomoly × ", trait_display_name[[trait_var]],
        "\n(other traits fixed at median)"
      ),
      x = "Standardised mean temperature anomaly",
      y = "Probability of occurrence"
    )
  
  p
}

p_CS_temp <- plot_temp_by_trait_with_legend(
  mod  = m_full,
  data = myoccdata,
  trait_var = "CS.rs",
  
  ## 手动微调参数（你可以继续改）
  base_size = 13,
  title_size = 16,
  legend_title_size = 14,
  legend_text_size = 12,
  line_width = 1.25,
  ribbon_alpha = 0.25
)

## 先查看
print(p_CS_temp)

## 再导出（三格式）
save_plot3(
  p_CS_temp,
  outdir = dir_best_plots,
  tag = "BEST_Temperature_x_ClutchSize_ABS_Legend_10_50_90",
  width = 13,
  height = 5.8
)

for (v in trait_vars_cont) {
  message("Plotting Temperature × ", v)
  
  p <- plot_temp_by_trait_with_legend(
    mod  = m_full,
    data = myoccdata,
    trait_var = v
  )
  
  print(p)
  
  save_plot3(
    p,
    outdir = dir_best_plots,
    tag = paste0("BEST_Temperature_x_", v, "_ABS_Legend_10_50_90"),
    width = 13,
    height = 5.8
  )
}


##%#############################################################%##
#  ADD-ON: Make 3-way plot like your example
#  Row = trait level (Low/Med/High)
#  Col = UI2 (land-use)
#  Color = UI2 (different colors for land-use)
#  Ribbon = 95% CI
#
#  CN: 生成你给的那种三重交互图：
#      每行=性状三档（10/50/90%），每列=土地利用类型(UI2)，
#      颜色=土地利用类型（每列一个颜色），带95%CI阴影
##%#############################################################%##

## ============================================================
## 0) A clean plotting function: Temp × UI2 × Trait(3-level)
## ============================================================
plot_threeway_temp_UI2_traitRows_colorUI2 <- function(
    mod, data,
    trait_var,
    temp_var = "StdTmeanAnomalyRS",
    qs = c(0.10, 0.50, 0.90),
    q_labels = c("Low (10%)", "Median (50%)", "High (90%)"),
    n = 220,
    ndraw = 2000,
    seed = 1,
    
    ## labels / 可手动调
    title_text = NULL,
    xlab = "Standardised mean temperature anomaly (StdTmeanAnomalyRS)",
    ylab = "Predicted probability of occurrence",
    base_size = 13,
    title_size = 13,
    strip_size = 11,
    line_width = 1.15,
    ribbon_alpha = 0.22
) {
  stopifnot(length(qs) == 3, length(q_labels) == 3)
  
  ## ---- trait 3-level values / 性状三档数值 ----
  q_vals <- as.numeric(quantile(data[[trait_var]], probs = qs, na.rm = TRUE, type = 7))
  names(q_vals) <- q_labels
  
  ## ---- temperature sequence / 温度序列 ----
  temp_seq <- seq(
    quantile(data[[temp_var]], 0.02, na.rm = TRUE),
    quantile(data[[temp_var]], 0.98, na.rm = TRUE),
    length.out = n
  )
  
  ## ---- build grid: UI2 × temp × Q / 构建预测网格 ----
  nd <- expand.grid(
    UI2 = levels(data$UI2),
    Q   = factor(q_labels, levels = q_labels),
    tmp = temp_seq
  )
  names(nd)[names(nd) == "tmp"] <- temp_var
  nd$UI2 <- factor(nd$UI2, levels = levels(data$UI2))
  
  ## ---- set focal trait by Q, others fixed at median / 其它性状中位数固定 ----
  for (v in trait_vars_cont) {
    if (v == trait_var) {
      nd[[v]] <- q_vals[as.character(nd$Q)]
    } else {
      nd[[v]] <- median(data[[v]], na.rm = TRUE)
    }
  }
  
  ## ---- predict draws (fixed effects) / 后验抽样预测 ----
  p_draw <- pred_draws_glmmTMB(mod, nd, ndraw = ndraw, seed = seed)
  summ   <- summ_draws(p_draw)
  df     <- cbind(nd, summ)
  
  ## ---- title default ----
  if (is.null(title_text)) {
    title_text <- paste0(
      "Three-way interaction: UI2 × standardised warming anomaly × ",
      trait_display_name[[trait_var]],
      " (95% CI)"
    )
  }
  
  ## ---- make plot ----
  p <- ggplot(df, aes(x = .data[[temp_var]])) +
    geom_ribbon(
      aes(ymin = low95, ymax = high95, fill = UI2),
      alpha = ribbon_alpha,
      colour = NA
    ) +
    geom_line(
      aes(y = mean, colour = UI2),
      linewidth = line_width
    ) +
    facet_grid(Q ~ UI2, scales = "fixed") +   # rows=Q, cols=UI2
    scale_colour_manual(values = pal, name = "Land-use type") +
    scale_fill_manual(values = pal, name = "Land-use type") +
    theme_classic(base_size = base_size) +
    theme(
      legend.position = "none",  # 像示例图一样：每列就是一种UI2颜色，所以可隐藏图例
      plot.title = element_text(size = title_size, face = "bold"),
      strip.text = element_text(size = strip_size, face = "bold"),
      strip.background = element_rect(fill = "white", color = "black")
    ) +
    labs(
      title = title_text,
      x = xlab,
      y = ylab
    )
  
  p
}


plot_threeway_temp_UI2_traitRows_colorUI2 <- function(
    mod, data,
    trait_var,
    temp_var = "StdTmeanAnomalyRS",
    qs = c(0.10, 0.50, 0.90),
    q_labels = c("Low (10%)", "Median (50%)", "High (90%)"),
    n = 220,
    ndraw = 2000,
    seed = 1,
    
    ## labels / 可手调
    title_text = NULL,
    xlab = "Standardised mean temperature anomaly (StdTmeanAnomalyRS)",
    ylab = "Predicted probability of occurrence",
    base_size = 13,
    title_size = 14,
    strip_size = 10,        # ← 行标签字体稍小
    line_width = 1.15,
    ribbon_alpha = 0.22
) {
  
  ## ---- trait 3-level values ----
  q_vals <- as.numeric(quantile(data[[trait_var]], probs = qs, na.rm = TRUE, type = 7))
  names(q_vals) <- q_labels
  
  ## ---- temperature sequence ----
  temp_seq <- seq(
    quantile(data[[temp_var]], 0.02, na.rm = TRUE),
    quantile(data[[temp_var]], 0.98, na.rm = TRUE),
    length.out = n
  )
  
  ## ---- prediction grid ----
  nd <- expand.grid(
    UI2 = levels(data$UI2),
    Q   = factor(q_labels, levels = q_labels),
    tmp = temp_seq
  )
  names(nd)[names(nd) == "tmp"] <- temp_var
  nd$UI2 <- factor(nd$UI2, levels = levels(data$UI2))
  
  ## ---- set traits ----
  for (v in trait_vars_cont) {
    if (v == trait_var) {
      nd[[v]] <- q_vals[as.character(nd$Q)]
    } else {
      nd[[v]] <- median(data[[v]], na.rm = TRUE)
    }
  }
  
  ## ---- predict ----
  p_draw <- pred_draws_glmmTMB(mod, nd, ndraw = ndraw, seed = seed)
  summ   <- summ_draws(p_draw)
  df     <- cbind(nd, summ)
  
  ## ============================================================
  ## NEW: build row strip labels / 构造“性状 + 分位”的行标签
  ## ============================================================
  trait_name <- trait_display_name[[trait_var]]
  
  row_lab_map <- c(
    "Low (10%)"     = paste0(trait_name, " (Lower 10%)"),
    "Median (50%)"  = paste0(trait_name, " (Median)"),
    "High (90%)"    = paste0(trait_name, " (Upper 90%)")
  )
  
  ## ---- default title ----
  if (is.null(title_text)) {
    title_text <- paste0(
      "Three-way interaction: UI2 × standardised warming anomaly × ",
      trait_name
    )
  }
  
  ## ---- plot ----
  p <- ggplot(df, aes(x = .data[[temp_var]])) +
    geom_ribbon(
      aes(ymin = low95, ymax = high95, fill = UI2),
      alpha = ribbon_alpha,
      colour = NA
    ) +
    geom_line(
      aes(y = mean, colour = UI2),
      linewidth = line_width
    ) +
    facet_grid(
      Q ~ UI2,
      labeller = labeller(Q = row_lab_map)   # ← 行标签替换在这里
    ) +
    scale_colour_manual(values = pal) +
    scale_fill_manual(values = pal) +
    theme_classic(base_size = base_size) +
    theme(
      legend.position = "none",
      plot.title = element_text(size = title_size, face = "bold"),
      strip.text.y = element_text(size = strip_size),   # ← 行标签字体
      strip.text.x = element_text(size = strip_size + 1, face = "bold"),
      strip.background = element_rect(fill = "white", colour = "black")
    ) +
    labs(
      title = title_text,
      x = xlab,
      y = ylab
    )
  
  p
}

## ============================================================
## 1) Example: Clutch size (CS.rs) / 示例：Clutch size
## ============================================================
p_three_CS <- plot_threeway_temp_UI2_traitRows_colorUI2(
  mod  = m_full,          # 或 m_full
  data = myoccdata,
  trait_var = "CS.rs",
  
  ## Low/Median/High = 10/50/90
  qs = c(0.10, 0.50, 0.90),
  q_labels = c("Low (10%)", "Median (50%)", "High (90%)"),
  
  ## 手动调参（可按你审美继续改）
  base_size = 13,
  title_size = 14,
  strip_size = 11,
  line_width = 1.2,
  ribbon_alpha = 0.22
)

## 先查看
print(p_three_CS)

## 再导出（三格式）
save_plot3(
  p_three_CS,
  outdir = dir_best_plots,
  tag = "BEST_ThreeWay_UI2_x_Temp_x_ClutchSize_rowsTrait_colsUI2_colorUI2",
  width = 12.5,
  height = 6.8
)

## ============================================================
## 2) Run for ALL traits / 对所有性状批量输出（可选）
## ============================================================
for (v in trait_vars_cont) {
  message("3-way plotting (rows=trait levels, cols=UI2): ", v)
  
  p <- plot_threeway_temp_UI2_traitRows_colorUI2(
    mod  = m_full,
    data = myoccdata,
    trait_var = v,
    qs = c(0.10, 0.50, 0.90),
    q_labels = c("Low (10%)", "Median (50%)", "High (90%)")
  )
  
  print(p)
  
  save_plot3(
    p,
    outdir = dir_best_plots,
    tag = paste0("BEST_ThreeWay_UI2_x_Temp_x_", v, "_rowsTrait_colsUI2_colorUI2"),
    width = 12.5,
    height = 6.8
  )
}

AIC(m_full)
