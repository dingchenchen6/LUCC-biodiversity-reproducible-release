##%######################################################%##
#                                                          #
####   Run climate × LUI models (lme4 + brms) + plots   ####
#                                                          #
##%######################################################%##
# 目的 / Purpose:
#   - 复现 Tim/Charlie 的 LUI × climate anomaly 模型
#   - 同时拟合 lme4（频率学派）与 brms（贝叶斯）
#   - 完整输出：数据分布图、模型诊断、tab_model 表、预测曲线、模型对比图
#
# 输出 / Outputs:
#   - outDir/ 下面保存：模型对象(rds)、诊断图(pdf/png/pptx)、表格(html/csv)、预测数据(rds)
#
# 注意 / Notes:
#   - 本脚本假设你的 Functions.R 已包含：
#       StdCenterPredictor(), BackTransformCentreredPredictor(), RescaleAbundance()
#   - 若没有 BackTransformCentreredPredictor()，请用你自己的反标准化函数替代
#
# 作者 / Author: (整理版 by ChatGPT)
##%######################################################%##

rm(list = ls())

## =========================
## 0) Directories
## =========================
predictsDataDir <- "5_PREDICTSMatchPropNatHab/"
outDir <- "6_RunLUClimateModels_brms/"
if(!dir.exists(outDir)) dir.create(outDir, recursive = TRUE)

log_file <- file.path(outDir, "log_full_pipeline.txt")
sink(log_file, split = TRUE)

t.start <- Sys.time()
cat("Start time: ", as.character(t.start), "\n\n")

## =========================
## 1) Libraries
## =========================
pkg_needed <- c(
  "devtools","StatisticalModels","predictsFunctions",
  "sjPlot","ggplot2","cowplot",
  "lme4","MASS","broom.mixed","performance",
  "DHARMa",
  "brms","bayesplot","cmdstanr","loo",
  "officer","rvg"
)

to_install <- pkg_needed[!sapply(pkg_needed, requireNamespace, quietly = TRUE)]
if(length(to_install) > 0){
  cat("Installing missing packages:\n", paste(to_install, collapse = ", "), "\n")
  install.packages(to_install, dependencies = TRUE)
}

library(StatisticalModels)
library(predictsFunctions)
library(sjPlot)
library(ggplot2)
library(cowplot)
library(lme4)
library(MASS)
library(broom.mixed)
library(performance)
library(DHARMa)
library(brms)
library(bayesplot)
library(cmdstanr)
library(loo)
library(officer)
library(rvg)

## Functions.R (你的自定义函数)
source("Functions.R")

## =========================
## 2) Output subfolders
## =========================
plotDir  <- file.path(outDir, "plots_all")
modelDir <- file.path(outDir, "models_all")
dataDir  <- file.path(outDir, "data_all")
tabDir   <- file.path(outDir, "tables_all")

dir.create(plotDir,  showWarnings = FALSE, recursive = TRUE)
dir.create(modelDir, showWarnings = FALSE, recursive = TRUE)
dir.create(dataDir,  showWarnings = FALSE, recursive = TRUE)
dir.create(tabDir,   showWarnings = FALSE, recursive = TRUE)

## =========================
## 3) Helpers: save plots in 3 formats (pdf/png/pptx)
## =========================
save_plot_3formats <- function(p, prefix,
                               width_mm = 183, height_mm = 100,
                               dpi_png = 600){
  stopifnot(inherits(p, "ggplot") || inherits(p, "gtable") || inherits(p, "patchwork"))
  pdf_path <- file.path(plotDir, paste0(prefix, ".pdf"))
  png_path <- file.path(plotDir, paste0(prefix, ".png"))
  ppt_path <- file.path(plotDir, paste0(prefix, ".pptx"))
  
  # PDF
  ggsave(pdf_path, plot = p, width = width_mm, height = height_mm, units = "mm", dpi = 300)
  
  # PNG
  ggsave(png_path, plot = p, width = width_mm, height = height_mm, units = "mm", dpi = dpi_png)
  
  # PPTX (vector)
  doc <- read_pptx()
  doc <- add_slide(doc, layout = "Title and Content", master = "Office Theme")
  doc <- ph_with(doc, value = prefix, location = ph_location_type(type = "title"))
  doc <- ph_with(doc, dml(ggobj = p), location = ph_location_fullsize())
  print(doc, target = ppt_path)
  
  invisible(list(pdf = pdf_path, png = png_path, pptx = ppt_path))
}

## =========================
## 4) Theme & palette (风格保持不变 / keep style)
## =========================
ui2_cols <- c("#009E73", "#0072B2", "#E69F00", "#D55E00", "pink")
ui2_levs <- c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High","Urban")

theme_pub <- theme_bw() +
  theme(
    aspect.ratio = 1,
    title = element_text(size = 8, face = "bold"),
    axis.text = element_text(size = 7),
    axis.title = element_text(size = 7),
    legend.background = element_blank(),
    legend.text = element_text(size = 6),
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.2),
    panel.border = element_rect(linewidth = 0.2),
    axis.ticks = element_line(linewidth = 0.2)
  )

## =========================
## 5) Read data & preprocess
## =========================
cat("\n[Step] Read predictsSites...\n")
predictsSites <- readRDS(file.path(predictsDataDir, "PREDICTSSitesWithClimateAndNatHab1.rds"))
predictsSites <- predictsSites@data

# UI2 factor & reference level
predictsSites$UI2 <- factor(predictsSites$UI2)
if("Primary vegetation" %in% levels(predictsSites$UI2)){
  predictsSites$UI2 <- relevel(predictsSites$UI2, ref = "Primary vegetation")
}
predictsSites <- droplevels(predictsSites)

# Standardise climate anomalies (中心化+标准化 / center+scale)
predictsSites$StdTmeanAnomalyRS <- StdCenterPredictor(predictsSites$StdTmeanAnomaly)
predictsSites$StdTmaxAnomalyRS  <- StdCenterPredictor(predictsSites$StdTmaxAnomaly)

# Combine natural habitat (PV+SV)
predictsSites$NH_1000  <- predictsSites$PV_1000  + predictsSites$SV_1000
predictsSites$NH_3000  <- predictsSites$PV_3000  + predictsSites$SV_3000
predictsSites$NH_5000  <- predictsSites$PV_5000  + predictsSites$SV_5000
predictsSites$NH_10000 <- predictsSites$PV_10000 + predictsSites$SV_10000

# Standardise NH
predictsSites$NH_1000.rs  <- StdCenterPredictor(predictsSites$NH_1000)
predictsSites$NH_3000.rs  <- StdCenterPredictor(predictsSites$NH_3000)
predictsSites$NH_5000.rs  <- StdCenterPredictor(predictsSites$NH_5000)
predictsSites$NH_10000.rs <- StdCenterPredictor(predictsSites$NH_10000)

# Rescale abundance and log values
# 注意 / Important:
#   RescaleAbundance() 通常会生成 LogAbund，常见定义是 log(Total_abundance + 0.01)
predictsSites <- RescaleAbundance(predictsSites)
predictsSites <- droplevels(predictsSites)

# Remove NA climate thresholds
predictsSites <- predictsSites[!is.na(predictsSites$avg_temp), ]

# Force UI2 level order consistent (如果某些 level 不存在也没关系)
predictsSites$UI2 <- factor(predictsSites$UI2, levels = intersect(ui2_levs, levels(predictsSites$UI2)))
predictsSites <- droplevels(predictsSites)

# Save cleaned data (文件名统一修正)
saveRDS(predictsSites, file = file.path(dataDir, "PREDICTSSiteData_clean.rds"))

cat("[Info] N sites:", nrow(predictsSites), "\n")
cat("[Info] UI2 levels:", paste(levels(predictsSites$UI2), collapse = ", "), "\n")
cat("[Info] Unique SSBS:", length(unique(predictsSites$SSBS)), "\n\n")

## =========================
## 6) Data distribution checks (显示图像 + 保存三种格式)
## =========================
cat("\n[Step] Data distribution plots...\n")

# 6.1 LogAbund distribution overall
p_dist_ab1 <- ggplot(predictsSites, aes(x = LogAbund)) +
  geom_histogram(bins = 60) +
  labs(x = "LogAbund", y = "Count", title = "Distribution of LogAbund (overall)") +
  theme_pub
save_plot_3formats(p_dist_ab1, "dist_LogAbund_overall", width_mm = 160, height_mm = 120)

# 6.2 LogAbund by UI2
p_dist_ab2 <- ggplot(predictsSites, aes(x = LogAbund, fill = UI2)) +
  geom_density(alpha = 0.25) +
  scale_fill_manual(values = ui2_cols) +
  labs(x = "LogAbund", y = "Density", title = "LogAbund density by UI2") +
  theme_pub +
  theme(legend.position = "right")
save_plot_3formats(p_dist_ab2, "dist_LogAbund_by_UI2", width_mm = 180, height_mm = 120)

# 6.3 Species richness distribution
p_dist_sr1 <- ggplot(predictsSites, aes(x = Species_richness)) +
  geom_histogram(bins = 80) +
  labs(x = "Species richness", y = "Count", title = "Distribution of species richness") +
  theme_pub
save_plot_3formats(p_dist_sr1, "dist_SpeciesRichness_overall", width_mm = 160, height_mm = 120)

# 6.4 Climate anomaly distributions
p_dist_clim <- ggplot(predictsSites, aes(x = StdTmeanAnomalyRS, fill = UI2)) +
  geom_density(alpha = 0.25) +
  scale_fill_manual(values = ui2_cols) +
  labs(x = "StdTmeanAnomalyRS", y = "Density", title = "StdTmeanAnomalyRS by UI2") +
  theme_pub +
  theme(legend.position = "right")
save_plot_3formats(p_dist_clim, "dist_StdTmeanAnomalyRS_by_UI2", width_mm = 180, height_mm = 120)

## =========================
## 7) Build model datasets
## =========================
cat("\n[Step] Build model datasets...\n")

ab_mean <- subset(predictsSites, !is.na(LogAbund) & !is.na(StdTmeanAnomalyRS))
sr_mean <- subset(predictsSites, !is.na(Species_richness) & !is.na(StdTmeanAnomalyRS))

ab_max  <- subset(predictsSites, !is.na(LogAbund) & !is.na(StdTmaxAnomalyRS))
sr_max  <- subset(predictsSites, !is.na(Species_richness) & !is.na(StdTmaxAnomalyRS))

saveRDS(ab_mean, file.path(dataDir, "ab_mean.rds"))
saveRDS(sr_mean, file.path(dataDir, "sr_mean.rds"))
saveRDS(ab_max,  file.path(dataDir, "ab_max.rds"))
saveRDS(sr_max,  file.path(dataDir, "sr_max.rds"))

cat("[Info] ab_mean:", nrow(ab_mean), " sr_mean:", nrow(sr_mean), "\n")
cat("[Info] ab_max :", nrow(ab_max),  " sr_max :", nrow(sr_max),  "\n\n")

## =========================
## 8) lme4 models (ordinary linear terms; not poly)
## =========================
cat("\n[Step] Fit lme4 models...\n")

# 8.1 Abundance (mean anomaly) - Gaussian LMM
# 中文：LogAbund ~ UI2 * StdTmeanAnomalyRS + (1|SS) + (1|SSB)
m_lme4_ab_mean <- lmer(
  LogAbund ~ UI2 * StdTmeanAnomalyRS + (1|SS) + (1|SSB),
  data = ab_mean,
  REML = TRUE
)

# 8.2 Richness (mean anomaly) - Negative binomial GLMM
# 中文：用 glmer.nb 比 poisson 更能处理过度离散
m_lme4_sr_mean_nb <- glmer.nb(
  Species_richness ~ UI2 * StdTmeanAnomalyRS + (1|SS) + (1|SSB) + (1|SSBS),
  data = sr_mean,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

# 8.3 Abundance (max anomaly)
m_lme4_ab_max <- lmer(
  LogAbund ~ UI2 * StdTmaxAnomalyRS + (1|SS) + (1|SSB),
  data = ab_max,
  REML = TRUE
)

# 8.4 Richness (max anomaly) - NB
m_lme4_sr_max_nb <- glmer.nb(
  Species_richness ~ UI2 * StdTmaxAnomalyRS + (1|SS) + (1|SSB) + (1|SSBS),
  data = sr_max,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

saveRDS(m_lme4_ab_mean,   file.path(modelDir, "lme4_ab_mean.rds"))
saveRDS(m_lme4_sr_mean_nb,file.path(modelDir, "lme4_sr_mean_nb.rds"))
saveRDS(m_lme4_ab_max,    file.path(modelDir, "lme4_ab_max.rds"))
saveRDS(m_lme4_sr_max_nb, file.path(modelDir, "lme4_sr_max_nb.rds"))

## 8.5 tab_model outputs (HTML)
tab_model(m_lme4_ab_mean, transform = NULL,
          file = file.path(tabDir, "lme4_Abundance_MeanAnom.html"))
tab_model(m_lme4_sr_mean_nb, transform = NULL,
          file = file.path(tabDir, "lme4_Richness_MeanAnom_NB.html"))
tab_model(m_lme4_ab_max, transform = NULL,
          file = file.path(tabDir, "lme4_Abundance_MaxAnom.html"))
tab_model(m_lme4_sr_max_nb, transform = NULL,
          file = file.path(tabDir, "lme4_Richness_MaxAnom_NB.html"))

## =========================
## 9) lme4 diagnostics (DHARMa + performance)
## =========================
cat("\n[Step] lme4 diagnostics...\n")

save_dharma <- function(model, prefix){
  # DHARMa 对 LMM 也可用，但主要对 GLMM 更关键
  simres <- DHARMa::simulateResiduals(model, n = 1000)
  png(file.path(plotDir, paste0(prefix, "_DHARMa_residual.png")), 2400, 1800, res = 300)
  plot(simres)
  dev.off()
  
  # 额外检验：dispersion / zero-inflation（对计数模型有意义）
  if(inherits(model, "glmerMod")){
    zi <- DHARMa::testZeroInflation(simres)
    dp <- DHARMa::testDispersion(simres)
    out <- DHARMa::testOutliers(simres)
    
    writeLines(c(
      paste0("Model: ", prefix),
      paste0("ZeroInflation test p: ", zi$p.value),
      paste0("Dispersion test p: ", dp$p.value),
      paste0("Outliers test p: ", out$p.value)
    ), con = file.path(plotDir, paste0(prefix, "_DHARMa_tests.txt")))
  }
  
  # performance 包：R2 / ICC / check_collinearity 等
  perf_txt <- capture.output({
    cat("\n== performance::model_performance ==\n")
    print(performance::model_performance(model))
    cat("\n== performance::r2 ==\n")
    print(performance::r2(model))
    cat("\n== performance::icc ==\n")
    print(performance::icc(model))
    cat("\n== performance::check_collinearity ==\n")
    print(performance::check_collinearity(model))
  })
  writeLines(perf_txt, con = file.path(plotDir, paste0(prefix, "_performance.txt")))
  
  invisible(TRUE)
}

save_dharma(m_lme4_ab_mean,    "lme4_ab_mean")
save_dharma(m_lme4_sr_mean_nb, "lme4_sr_mean_nb")
save_dharma(m_lme4_ab_max,     "lme4_ab_max")
save_dharma(m_lme4_sr_max_nb,  "lme4_sr_max_nb")

## =========================
## 10) brms backend setup
## =========================
cat("\n[Step] Setup cmdstanr...\n")
if (!cmdstanr::cmdstan_version(TRUE)) {
  cat("CmdStan not found -> installing...\n")
  cmdstanr::install_cmdstan()
}

## =========================
## 11) brms priors
## =========================
# 中文：这里给相对温和的先验，避免过强约束；你也可用你已有 prior_student/prior_nb
b_prior_gaus <- c(
  prior(normal(0, 1), class = "b"),
  prior(student_t(3, 0, 2.5), class = "sd"),
  prior(student_t(3, 0, 2.5), class = "sigma")
)

b_prior_nb <- c(
  prior(normal(0, 2), class = "b"),
  prior(student_t(3, 0, 2.5), class = "sd"),
  prior(gamma(2, 0.2), class = "shape")
)

## =========================
## 12) brms models (ordinary linear terms)
## =========================
cat("\n[Step] Fit brms models...\n")

# 12.1 Abundance mean anomaly (Gaussian)
b_brms_ab_mean <- brm(
  LogAbund ~ UI2 * StdTmeanAnomalyRS + (1|SS) + (1|SSB),
  data    = ab_mean,
  family  = gaussian(),
  prior   = b_prior_gaus,
  chains  = 4, iter = 6000, warmup = 3000,
  seed    = 123,
  backend = "cmdstanr",
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  save_pars = save_pars(all = TRUE)
)
pp_check(b_brms_ab_mean)
# 12.2 Richness mean anomaly (NegBin)
b_brms_sr_mean <- brm(
  Species_richness ~ UI2 * StdTmeanAnomalyRS + (1|SS) + (1|SSB) + (1|SSBS),
  data    = sr_mean,
  family  = negbinomial(),
  prior   = b_prior_nb,
  chains  = 4, iter = 8000, warmup = 4000,
  seed    = 123,
  backend = "cmdstanr",
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  save_pars = save_pars(all = TRUE)
)
pp_check(b_brms_sr_mean)
# 12.3 Abundance max anomaly
b_brms_ab_max <- brm(
  LogAbund ~ UI2 * StdTmaxAnomalyRS + (1|SS) + (1|SSB),
  data    = ab_max,
  family  = gaussian(),
  prior   = b_prior_gaus,
  chains  = 4, iter = 6000, warmup = 3000,
  seed    = 123,
  backend = "cmdstanr",
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  save_pars = save_pars(all = TRUE)
)

# 12.4 Richness max anomaly
b_brms_sr_max <- brm(
  Species_richness ~ UI2 * StdTmaxAnomalyRS + (1|SS) + (1|SSB) + (1|SSBS),
  data    = sr_max,
  family  = negbinomial(),
  prior   = b_prior_nb,
  chains  = 4, iter = 8000, warmup = 4000,
  seed    = 123,
  backend = "cmdstanr",
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  save_pars = save_pars(all = TRUE)
)

saveRDS(b_brms_ab_mean, file.path(modelDir, "brms_ab_mean.rds"))
saveRDS(b_brms_sr_mean, file.path(modelDir, "brms_sr_mean.rds"))
saveRDS(b_brms_ab_max,  file.path(modelDir, "brms_ab_max.rds"))
saveRDS(b_brms_sr_max,  file.path(modelDir, "brms_sr_max.rds"))

## =========================
## 13) brms diagnostics (pp_check + NUTS + loo + bayes_R2)
## =========================
cat("\n[Step] brms diagnostics...\n")

save_brms_diag <- function(fit, prefix){
  # 13.1 summary text
  writeLines(capture.output(print(summary(fit))),
             con = file.path(plotDir, paste0(prefix, "_summary.txt")))
  
  # 13.2 NUTS params
  np <- nuts_params(fit)
  n_div <- sum(np$Parameter == "divergent__" & np$Value == 1)
  bfmi_vals <- bayesplot::bfmi(np)
  max_td <- max(np$Value[np$Parameter == "treedepth__"], na.rm = TRUE)
  writeLines(c(
    paste0("Model: ", prefix),
    paste0("Divergent transitions: ", n_div),
    paste0("BFMI: ", paste(round(bfmi_vals, 3), collapse = ", ")),
    paste0("Max treedepth observed: ", max_td)
  ), con = file.path(plotDir, paste0(prefix, "_NUTS.txt")))
  
  # 13.3 trace + dens (前10个参数)
  pars <- head(parnames(fit), 10)
  png(file.path(plotDir, paste0(prefix, "_trace.png")), 2400, 1800, res = 300)
  print(mcmc_trace(as.array(fit), pars = pars))
  dev.off()
  png(file.path(plotDir, paste0(prefix, "_dens.png")), 2400, 1800, res = 300)
  print(mcmc_dens(as.array(fit), pars = pars))
  dev.off()
  
  # 13.4 PPC
  p_pp <- pp_check(fit, ndraws = 200)
  save_plot_3formats(p_pp, paste0(prefix, "_ppcheck"), width_mm = 180, height_mm = 120)
  
  # 13.5 bayes_R2
  r2 <- brms::bayes_R2(fit)
  write.csv(data.frame(R2 = as.numeric(r2)),
            file.path(tabDir, paste0(prefix, "_bayesR2_draws.csv")),
            row.names = FALSE)
  
  # 13.6 loo
  loo_res <- loo::loo(fit)
  saveRDS(loo_res, file.path(modelDir, paste0(prefix, "_loo.rds")))
  writeLines(capture.output(print(loo_res)),
             con = file.path(plotDir, paste0(prefix, "_loo.txt")))
  
  invisible(TRUE)
}

save_brms_diag(b_brms_ab_mean, "brms_ab_mean")
save_brms_diag(b_brms_sr_mean, "brms_sr_mean")
save_brms_diag(b_brms_ab_max,  "brms_ab_max")
save_brms_diag(b_brms_sr_max,  "brms_sr_max")



## =========================
## 14) Predictions + % change curves (REWRITE / 重写版)
## =========================
cat("\n[Step 14] Predictions & effect plots (rewritten, stable)...\n")

## ---------------------------------------------------------
## 14.0 统一的工具函数（稳定版）
## ---------------------------------------------------------

# (1) 构建预测网格 / Build prediction grid
# - xRS: 标准化后的 anomaly（如 StdTmeanAnomalyRS）
# - xOrigAll: 用于反标准化的原始 anomaly 向量（如 predictsSites$StdTmeanAnomaly）
build_grid2 <- function(df, xRS, xOrigName, xOrigAll,
                        grid_n = 300,
                        ui2_levels = levels(df$UI2)) {
  nd <- expand.grid(
    UI2 = factor(ui2_levels, levels = ui2_levels),
    xRS = seq(min(df[[xRS]], na.rm = TRUE),
              max(df[[xRS]], na.rm = TRUE),
              length.out = grid_n)
  )
  names(nd)[names(nd) == "xRS"] <- xRS
  
  # 反标准化 x 轴（保证和你原图一致）
  nd[[xOrigName]] <- BackTransformCentreredPredictor(
    transformedX = nd[[xRS]],
    originalX    = xOrigAll
  )
  
  # 参考行：Primary vegetation 且 x 最接近 0
  refRow <- which(
    nd$UI2 == "Primary vegetation" &
      abs(nd[[xOrigName]]) == min(abs(nd[[xOrigName]][nd$UI2 == "Primary vegetation"]), na.rm = TRUE)
  )[1]
  
  # 每个 UI2 的 2.5%~97.5% 截断区间（避免外推极端）
  Qlist <- lapply(ui2_levels, function(u) {
    quantile(df[[xRS]][df$UI2 == u], probs = c(0.025, 0.975), na.rm = TRUE)
  })
  names(Qlist) <- ui2_levels
  
  list(nd = nd, refRow = refRow, Qlist = Qlist)
}

# (2) quantile trimming：把超出区间的 Pred 置 NA（不会连线）
trim_quantile2 <- function(nd, xRS, Qlist) {
  for (u in names(Qlist)) {
    idx <- nd$UI2 == u & (nd[[xRS]] < Qlist[[u]][1] | nd[[xRS]] > Qlist[[u]][2])
    nd$PredMedian[idx] <- NA
    nd$PredLower[idx]  <- NA
    nd$PredUpper[idx]  <- NA
  }
  nd
}

# (3) brms：在 link(η) 尺度计算相对变化（强烈推荐）
#    rel = exp(eta - eta_ref)  -> 稳定，不受 -0.01 影响
pred_percent_brms_linkratio <- function(fit, nd, refRow) {
  eta <- posterior_linpred(fit, newdata = nd, re_formula = NA)  # draws x N
  rel <- exp(eta - eta[, refRow])                               # draws x N
  
  nd$PredMedian <- apply(rel, 2, median, na.rm = TRUE) * 100 - 100
  nd$PredLower  <- apply(rel, 2, quantile, probs = 0.025, na.rm = TRUE) * 100 - 100
  nd$PredUpper  <- apply(rel, 2, quantile, probs = 0.975, na.rm = TRUE) * 100 - 100
  nd
}

# (4) lme4：固定效应参数抽样（vcov）得到 η draws，然后 link ratio
#    - 对 lmer (Gaussian identity)：η = Xβ
#    - 对 glmer.nb / poisson (log link)：η = Xβ；ratio仍用 exp(η-η_ref)
pred_percent_lme4_linkratio_vcov <- function(fit, nd, refRow, nsim = 1000) {
  # 固定效应设计矩阵（仅 fixed effects）
  X <- model.matrix(delete.response(terms(fit)), nd)
  
  beta_hat <- lme4::fixef(fit)
  V <- as.matrix(vcov(fit))
  
  # 防止数值问题（偶尔 vcov 非正定）
  # 如遇到报错，可考虑 Matrix::nearPD
  draws_beta <- MASS::mvrnorm(n = nsim, mu = beta_hat, Sigma = V)
  
  eta <- draws_beta %*% t(X)     # nsim x N
  rel <- exp(eta - eta[, refRow])# nsim x N
  
  nd$PredMedian <- apply(rel, 2, median, na.rm = TRUE) * 100 - 100
  nd$PredLower  <- apply(rel, 2, quantile, probs = 0.025, na.rm = TRUE) * 100 - 100
  nd$PredUpper  <- apply(rel, 2, quantile, probs = 0.975, na.rm = TRUE) * 100 - 100
  nd
}

# (可选) 更严格的 lme4 不确定性：bootMer（会慢很多）
# pred_percent_lme4_linkratio_bootMer <- function(fit, nd, refRow, nsim = 500){
#   X <- model.matrix(delete.response(terms(fit)), nd)
#   FUN <- function(ff) {
#     eta <- drop(X %*% fixef(ff))
#     rel <- exp(eta - eta[refRow])
#     rel
#   }
#   bb <- lme4::bootMer(fit, FUN = FUN, nsim = nsim, use.u = FALSE, type = "parametric")
#   rel_mat <- bb$t   # nsim x N
#   nd$PredMedian <- apply(rel_mat, 2, median) * 100 - 100
#   nd$PredLower  <- apply(rel_mat, 2, quantile, 0.025) * 100 - 100
#   nd$PredUpper  <- apply(rel_mat, 2, quantile, 0.975) * 100 - 100
#   nd
# }

# (5) 统一作图：先排序，避免“锯齿连线”
make_effect_plot2 <- function(df, xname,
                              ylab_txt, xlab_txt,
                              x_breaks, x_limits, y_breaks, y_limits,
                              legend_pos = c(0.2,0.8),
                              show_legend = TRUE,
                              panel_title = NULL) {
  df <- df[order(df$UI2, df[[xname]]), ]  # <<<<<< 关键：排序
  
  p <- ggplot(df, aes(x = .data[[xname]], y = PredMedian, group = UI2)) +
    geom_line(aes(col = UI2), linewidth = 0.75, na.rm = TRUE) +
    geom_ribbon(aes(ymin = PredLower, ymax = PredUpper, fill = UI2),
                alpha = 0.2, na.rm = TRUE) +
    geom_hline(yintercept = 0, lty = "dashed", linewidth = 0.2) +
    scale_fill_manual(values = ui2_cols) +
    scale_colour_manual(values = ui2_cols) +
    scale_x_continuous(breaks = x_breaks, limits = x_limits) +
    scale_y_continuous(breaks = y_breaks, limits = y_limits) +
    labs(x = xlab_txt, y = ylab_txt) +
    theme_pub +
    theme(legend.position = if (show_legend) legend_pos else "none")
  
  if (!is.null(panel_title)) p <- p + ggtitle(panel_title)
  p
}

## ---------------------------------------------------------
## 14.1 Mean anomaly：Abundance（brms vs lme4）
## ---------------------------------------------------------

g_ab_mean2 <- build_grid2(
  df = ab_mean,
  xRS = "StdTmeanAnomalyRS",
  xOrigName = "StdTmeanAnomaly",
  xOrigAll = predictsSites$StdTmeanAnomaly,
  grid_n = 300,
  ui2_levels = levels(ab_mean$UI2)
)

nd_ab_mean_brms2 <- pred_percent_brms_linkratio(b_brms_ab_mean, g_ab_mean2$nd, g_ab_mean2$refRow)
nd_ab_mean_brms2 <- trim_quantile2(nd_ab_mean_brms2, "StdTmeanAnomalyRS", g_ab_mean2$Qlist)
nd_ab_mean_brms2$Model <- "brms"

nd_ab_mean_lme4_2 <- pred_percent_lme4_linkratio_vcov(m_lme4_ab_mean, g_ab_mean2$nd, g_ab_mean2$refRow, nsim = 1000)
nd_ab_mean_lme4_2 <- trim_quantile2(nd_ab_mean_lme4_2, "StdTmeanAnomalyRS", g_ab_mean2$Qlist)
nd_ab_mean_lme4_2$Model <- "lme4"

saveRDS(nd_ab_mean_brms2, file.path(dataDir, "pred_ab_mean_brms_STABLE.rds"))
saveRDS(nd_ab_mean_lme4_2, file.path(dataDir, "pred_ab_mean_lme4_STABLE.rds"))

p_ab_mean_brms2 <- make_effect_plot2(
  nd_ab_mean_brms2, "StdTmeanAnomaly",
  "Change in total abundance (%)", "Standardised Temperature Anomaly",
  x_breaks = c(0,0.5,1,1.5,2), x_limits = c(0,2),
  y_breaks = c(-50,-25,0,25,50,75,100,125,150), y_limits = c(-50,150),
  show_legend = TRUE, panel_title = "a (brms)"
)

p_ab_mean_lme4_2 <- make_effect_plot2(
  nd_ab_mean_lme4_2, "StdTmeanAnomaly",
  "Change in total abundance (%)", "Standardised Temperature Anomaly",
  x_breaks = c(0,0.5,1,1.5,2), x_limits = c(0,2),
  y_breaks = c(-50,-25,0,25,50,75,100,125,150), y_limits = c(-50,150),
  show_legend = FALSE, panel_title = "a (lme4)"
)

p_ab_mean_compare2 <- cowplot::plot_grid(p_ab_mean_brms2, p_ab_mean_lme4_2, ncol = 2)
save_plot_3formats(p_ab_mean_compare2, "Compare_MeanAnom_Abundance_brms_vs_lme4_STABLE",
                   width_mm = 183, height_mm = 100)


## ---------------------------------------------------------
## 14.2 Mean anomaly：Richness（brms vs lme4）
## ---------------------------------------------------------

g_sr_mean2 <- build_grid2(
  df = sr_mean,
  xRS = "StdTmeanAnomalyRS",
  xOrigName = "StdTmeanAnomaly",
  xOrigAll = predictsSites$StdTmeanAnomaly,
  grid_n = 300,
  ui2_levels = levels(sr_mean$UI2)
)

nd_sr_mean_brms2 <- pred_percent_brms_linkratio(b_brms_sr_mean, g_sr_mean2$nd, g_sr_mean2$refRow)
nd_sr_mean_brms2 <- trim_quantile2(nd_sr_mean_brms2, "StdTmeanAnomalyRS", g_sr_mean2$Qlist)
nd_sr_mean_brms2$Model <- "brms"

nd_sr_mean_lme4_2 <- pred_percent_lme4_linkratio_vcov(m_lme4_sr_mean_nb, g_sr_mean2$nd, g_sr_mean2$refRow, nsim = 1000)
nd_sr_mean_lme4_2 <- trim_quantile2(nd_sr_mean_lme4_2, "StdTmeanAnomalyRS", g_sr_mean2$Qlist)
nd_sr_mean_lme4_2$Model <- "lme4"

saveRDS(nd_sr_mean_brms2, file.path(dataDir, "pred_sr_mean_brms_STABLE.rds"))
saveRDS(nd_sr_mean_lme4_2, file.path(dataDir, "pred_sr_mean_lme4_STABLE.rds"))

p_sr_mean_brms2 <- make_effect_plot2(
  nd_sr_mean_brms2, "StdTmeanAnomaly",
  "Change in species richness (%)", "Standardised Temperature Anomaly",
  x_breaks = c(0,0.5,1,1.5,2), x_limits = c(0,2),
  y_breaks = c(-75,-50,-25,0,25,50,75,100), y_limits = c(-75,100),
  show_legend = FALSE, panel_title = "b (brms)"
)

p_sr_mean_lme4_2 <- make_effect_plot2(
  nd_sr_mean_lme4_2, "StdTmeanAnomaly",
  "Change in species richness (%)", "Standardised Temperature Anomaly",
  x_breaks = c(0,0.5,1,1.5,2), x_limits = c(0,2),
  y_breaks = c(-75,-50,-25,0,25,50,75,100), y_limits = c(-75,100),
  show_legend = FALSE, panel_title = "b (lme4)"
)

p_sr_mean_compare2 <- cowplot::plot_grid(p_sr_mean_brms2, p_sr_mean_lme4_2, ncol = 2)
save_plot_3formats(p_sr_mean_compare2, "Compare_MeanAnom_Richness_brms_vs_lme4_STABLE",
                   width_mm = 183, height_mm = 100)

# Figure2 风格（Mean anomaly）
p_fig2_brms2 <- cowplot::plot_grid(p_ab_mean_brms2, p_sr_mean_brms2, ncol = 2)
save_plot_3formats(p_fig2_brms2, "Figure2_MeanAnom_Abun_Rich_brms_STABLE",
                   width_mm = 183, height_mm = 100)

p_fig2_lme4_2 <- cowplot::plot_grid(p_ab_mean_lme4_2, p_sr_mean_lme4_2, ncol = 2)
save_plot_3formats(p_fig2_lme4_2, "Figure2_MeanAnom_Abun_Rich_lme4_STABLE",
                   width_mm = 183, height_mm = 100)


## ---------------------------------------------------------
## 14.3 Max anomaly：Abundance（brms vs lme4）
## ---------------------------------------------------------

g_ab_max2 <- build_grid2(
  df = ab_max,
  xRS = "StdTmaxAnomalyRS",
  xOrigName = "StdTmaxAnomaly",
  xOrigAll = predictsSites$StdTmaxAnomaly,
  grid_n = 150,
  ui2_levels = levels(ab_max$UI2)
)

nd_ab_max_brms2 <- pred_percent_brms_linkratio(b_brms_ab_max, g_ab_max2$nd, g_ab_max2$refRow)
nd_ab_max_brms2 <- trim_quantile2(nd_ab_max_brms2, "StdTmaxAnomalyRS", g_ab_max2$Qlist)
nd_ab_max_brms2$Model <- "brms"

nd_ab_max_lme4_2 <- pred_percent_lme4_linkratio_vcov(m_lme4_ab_max, g_ab_max2$nd, g_ab_max2$refRow, nsim = 1000)
nd_ab_max_lme4_2 <- trim_quantile2(nd_ab_max_lme4_2, "StdTmaxAnomalyRS", g_ab_max2$Qlist)
nd_ab_max_lme4_2$Model <- "lme4"

saveRDS(nd_ab_max_brms2, file.path(dataDir, "pred_ab_max_brms_STABLE.rds"))
saveRDS(nd_ab_max_lme4_2, file.path(dataDir, "pred_ab_max_lme4_STABLE.rds"))

p_ab_max_brms2 <- make_effect_plot2(
  nd_ab_max_brms2, "StdTmaxAnomaly",
  "Change in total abundance (%)", "Maximum Temperature Anomaly",
  x_breaks = c(-1,0,1,2,3,4,5), x_limits = c(-1,5),
  y_breaks = c(-60,-40,-20,0,20,40,60), y_limits = c(-65,60),
  show_legend = TRUE, panel_title = "a (brms)"
)

p_ab_max_lme4_2 <- make_effect_plot2(
  nd_ab_max_lme4_2, "StdTmaxAnomaly",
  "Change in total abundance (%)", "Maximum Temperature Anomaly",
  x_breaks = c(-1,0,1,2,3,4,5), x_limits = c(-1,5),
  y_breaks = c(-60,-40,-20,0,20,40,60), y_limits = c(-65,60),
  show_legend = FALSE, panel_title = "a (lme4)"
)

p_ab_max_compare2 <- cowplot::plot_grid(p_ab_max_brms2, p_ab_max_lme4_2, ncol = 2)
save_plot_3formats(p_ab_max_compare2, "Compare_MaxAnom_Abundance_brms_vs_lme4_STABLE",
                   width_mm = 183, height_mm = 100)


## ---------------------------------------------------------
## 14.4 Max anomaly：Richness（brms vs lme4）
## ---------------------------------------------------------

g_sr_max2 <- build_grid2(
  df = sr_max,
  xRS = "StdTmaxAnomalyRS",
  xOrigName = "StdTmaxAnomaly",
  xOrigAll = predictsSites$StdTmaxAnomaly,
  grid_n = 150,
  ui2_levels = levels(sr_max$UI2)
)

nd_sr_max_brms2 <- pred_percent_brms_linkratio(b_brms_sr_max, g_sr_max2$nd, g_sr_max2$refRow)
nd_sr_max_brms2 <- trim_quantile2(nd_sr_max_brms2, "StdTmaxAnomalyRS", g_sr_max2$Qlist)
nd_sr_max_brms2$Model <- "brms"

nd_sr_max_lme4_2 <- pred_percent_lme4_linkratio_vcov(m_lme4_sr_max_nb, g_sr_max2$nd, g_sr_max2$refRow, nsim = 1000)
nd_sr_max_lme4_2 <- trim_quantile2(nd_sr_max_lme4_2, "StdTmaxAnomalyRS", g_sr_max2$Qlist)
nd_sr_max_lme4_2$Model <- "lme4"

saveRDS(nd_sr_max_brms2, file.path(dataDir, "pred_sr_max_brms_STABLE.rds"))
saveRDS(nd_sr_max_lme4_2, file.path(dataDir, "pred_sr_max_lme4_STABLE.rds"))

p_sr_max_brms2 <- make_effect_plot2(
  nd_sr_max_brms2, "StdTmaxAnomaly",
  "Change in species richness (%)", "Maximum Temperature Anomaly",
  x_breaks = c(-1,0,1,2,3,4,5), x_limits = c(-1,5),
  y_breaks = c(-60,-40,-20,0,20,40,60), y_limits = c(-65,60),
  show_legend = FALSE, panel_title = "b (brms)"
)

p_sr_max_lme4_2 <- make_effect_plot2(
  nd_sr_max_lme4_2, "StdTmaxAnomaly",
  "Change in species richness (%)", "Maximum Temperature Anomaly",
  x_breaks = c(-1,0,1,2,3,4,5), x_limits = c(-1,5),
  y_breaks = c(-60,-40,-20,0,20,40,60), y_limits = c(-65,60),
  show_legend = FALSE, panel_title = "b (lme4)"
)

p_sr_max_compare2 <- cowplot::plot_grid(p_sr_max_brms2, p_sr_max_lme4_2, ncol = 2)
save_plot_3formats(p_sr_max_compare2, "Compare_MaxAnom_Richness_brms_vs_lme4_STABLE",
                   width_mm = 183, height_mm = 100)

# ExtendedData 风格（Max anomaly）
p_ext2_brms2 <- cowplot::plot_grid(p_ab_max_brms2, p_sr_max_brms2, ncol = 2)
save_plot_3formats(p_ext2_brms2, "ExtendedData2_MaxAnom_Abun_Rich_brms_STABLE",
                   width_mm = 183, height_mm = 100)

p_ext2_lme4_2 <- cowplot::plot_grid(p_ab_max_lme4_2, p_sr_max_lme4_2, ncol = 2)
save_plot_3formats(p_ext2_lme4_2, "ExtendedData2_MaxAnom_Abun_Rich_lme4_STABLE",
                   width_mm = 183, height_mm = 100)


## =========================
## 15) Coefficient comparison (forest plot: lme4 vs brms) [optional but useful]
## =========================
cat("\n[Step 15] Coefficient forest plots (lme4 vs brms)...\n")

tidy_lme4_fixed <- function(model, model_name){
  tt <- broom.mixed::tidy(model, effects = "fixed", conf.int = TRUE)
  tt$model <- model_name
  tt[, c("term","estimate","conf.low","conf.high","model")]
}

tidy_brms_fixed <- function(fit, model_name){
  ps <- as.data.frame(posterior_summary(fit, probs = c(0.025, 0.975)))
  ps$term <- rownames(ps)
  rownames(ps) <- NULL
  ps <- ps[grepl("^b_", ps$term), ]
  ps$term <- gsub("^b_", "", ps$term)
  out <- data.frame(
    term = ps$term,
    estimate = ps$Estimate,
    conf.low = ps$Q2.5,
    conf.high = ps$Q97.5,
    model = model_name
  )
  out
}

coef_plot2 <- function(df, title_txt){
  df <- df[df$term != "Intercept", ]
  df$term <- factor(df$term, levels = rev(unique(df$term)))
  
  ggplot(df, aes(x = estimate, y = term, xmin = conf.low, xmax = conf.high, shape = model)) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3) +
    geom_pointrange(position = position_dodge(width = 0.6), linewidth = 0.3) +
    labs(x = "Effect size (link scale)", y = NULL, title = title_txt) +
    theme_pub +
    theme(legend.position = "bottom", aspect.ratio = NA)
}

# Mean anomaly abundance
df_coef_ab_mean2 <- rbind(
  tidy_lme4_fixed(m_lme4_ab_mean, "lme4"),
  tidy_brms_fixed(b_brms_ab_mean, "brms")
)
p_coef_ab_mean2 <- coef_plot2(df_coef_ab_mean2, "Fixed effects: Abundance (mean anomaly)")
save_plot_3formats(p_coef_ab_mean2, "Coef_Abundance_MeanAnom_lme4_vs_brms_STABLE",
                   width_mm = 180, height_mm = 150)

# Mean anomaly richness
df_coef_sr_mean2 <- rbind(
  tidy_lme4_fixed(m_lme4_sr_mean_nb, "lme4"),
  tidy_brms_fixed(b_brms_sr_mean, "brms")
)
p_coef_sr_mean2 <- coef_plot2(df_coef_sr_mean2, "Fixed effects: Richness (mean anomaly)")
save_plot_3formats(p_coef_sr_mean2, "Coef_Richness_MeanAnom_lme4_vs_brms_STABLE",
                   width_mm = 180, height_mm = 160)


## =========================
## 16) Finish
## =========================
cat("\n[Step 16] Finished Step14–Step15 rewrite.\n")
cat("Key change: all percent-change computed on link scale via exp(eta - eta_ref).\n")
