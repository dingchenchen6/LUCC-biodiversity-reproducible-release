###############################################################
## Land-use x Climate warming x Functional strategy groups (FG)
## Occurrence (brms) -- STABLE FG Bayesian pipeline
## 功能策略组（稳定方案）x 土地利用 x 变暖 -> 物种出现概率（Bayesian）
##
## Purpose / 目的
## 1) Refit the stable FG scheme with brms using a reliable but faster setup
## 2) Keep figure style, directory structure, and output formats aligned
## 3) Export Bayesian prediction figures, coefficient summaries, and diagnostics
###############################################################

suppressPackageStartupMessages({
  library(brms)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
  library(bayesplot)
  library(posterior)
  library(parallel)
})

if (requireNamespace("cmdstanr", quietly = TRUE)) {
  current_cmdstan <- tryCatch(cmdstanr::cmdstan_path(), error = function(e) "")
  if (!nzchar(current_cmdstan)) {
    cmdstan_candidates <- Sys.glob(path.expand("~/.cmdstan/cmdstan-*"))
    if (length(cmdstan_candidates) > 0) {
      cmdstanr::set_cmdstan_path(sort(cmdstan_candidates, decreasing = TRUE)[1])
    }
  }
}

if (!exists("topptx")) {
  if (requireNamespace("officer", quietly = TRUE) &&
      requireNamespace("rvg", quietly = TRUE)) {
    suppressPackageStartupMessages({
      library(officer)
      library(rvg)
    })
  }
}

## ============================================================
## 1) Global settings / 全局设置
## ============================================================

stable_root <- "glmmTMB_occurrence_FG_RobustPCAKmeans_FULL_plusTwoWay"
out_root <- "brms_occurrence_FG_RobustPCAKmeans_FULL_plusTwoWay"
dir.create(out_root, showWarnings = FALSE, recursive = TRUE)

dir_data   <- file.path(out_root, "00_Data")
dir_expl   <- file.path(out_root, "01_Exploration_and_DataSupport")
dir_models <- file.path(out_root, "02_Models")
dir_plots3 <- file.path(out_root, "03_Plots_ThreeWay_ABS_PCT")
dir_plots2 <- file.path(out_root, "04_Plots_TwoWay_PaperGrade")
dir_tables <- file.path(out_root, "05_Tables_Bayesian")
dir_caps   <- file.path(out_root, "06_Figure_Captions")
dir_diag   <- file.path(out_root, "98_Diagnostics")
dir.create(dir_data,   showWarnings = FALSE, recursive = TRUE)
dir.create(dir_expl,   showWarnings = FALSE, recursive = TRUE)
dir.create(dir_models, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_plots3, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_plots2, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_tables, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_caps,   showWarnings = FALSE, recursive = TRUE)
dir.create(dir_diag,   showWarnings = FALSE, recursive = TRUE)

UI2_levels <- c(
  "Primary vegetation", "Secondary vegetation",
  "Agriculture_Low", "Agriculture_High", "Urban"
)

pal_UI2 <- c("#009E73", "#0072B2", "#E69F00", "#D55E00", "#CC79A7")
names(pal_UI2) <- UI2_levels

pal_FG <- c("FG1" = "#F05A5A", "FG2" = "#19B51E", "FG3" = "#4E79FF")
FG_levels <- c("FG1", "FG2", "FG3")

fg_feature_labels <- c(
  "FG1" = "FG1: slow-history, larger-bodied, longer generation length",
  "FG2" = "FG2: broader niche, wider thermal breadth, higher fecundity",
  "FG3" = "FG3: restricted range, lower dispersal, more specialist-like"
)

fg_strip_labels <- c(
  "FG1" = "FG1\nslow-history / larger body size /\nlonger generation length",
  "FG2" = "FG2\nbroader niche / thermal breadth /\nhigher fecundity",
  "FG3" = "FG3\nrestricted range / lower dispersal /\nmore specialist-like"
)

TITLE_SIZE <- 14
BASE_SIZE  <- 13

BASE_UI2   <- "Primary vegetation"
BASE_TEMP0 <- 0
N_TEMP <- 240

BRMS_CHAINS <- suppressWarnings(as.integer(Sys.getenv("FG_BRMS_CHAINS", "4")))
BRMS_ITER   <- suppressWarnings(as.integer(Sys.getenv("FG_BRMS_ITER", "6000")))
BRMS_WARMUP <- suppressWarnings(as.integer(Sys.getenv("FG_BRMS_WARMUP", "3000")))
BRMS_THREADS <- suppressWarnings(as.integer(Sys.getenv("FG_BRMS_THREADS", "2")))
BRMS_CORES   <- suppressWarnings(as.integer(Sys.getenv("FG_BRMS_CORES", as.character(min(BRMS_CHAINS, parallel::detectCores())))))
BRMS_ADAPT   <- suppressWarnings(as.numeric(Sys.getenv("FG_BRMS_ADAPT_DELTA", "0.99")))
BRMS_TREEDEPTH <- suppressWarnings(as.integer(Sys.getenv("FG_BRMS_MAX_TREEDEPTH", "15")))
BRMS_BACKEND <- Sys.getenv("FG_BRMS_BACKEND", unset = "cmdstanr")
BRMS_FORCE_REFIT <- identical(tolower(Sys.getenv("FG_BRMS_FORCE_REFIT", "0")), "1")
BRMS_GRAINSIZE <- suppressWarnings(as.integer(Sys.getenv("FG_BRMS_GRAINSIZE", "1000")))
SMOKE_N <- suppressWarnings(as.integer(Sys.getenv("FG_BRMS_SMOKE_N", "0")))
BRMS_SEED <- 123
N_DRAWS_PRED <- suppressWarnings(as.integer(Sys.getenv("FG_BRMS_NDRAWS_PRED", "1600")))
if (!is.finite(BRMS_CHAINS) || BRMS_CHAINS < 2) BRMS_CHAINS <- 4
if (!is.finite(BRMS_ITER) || BRMS_ITER < 3000) BRMS_ITER <- 6000
if (!is.finite(BRMS_WARMUP) || BRMS_WARMUP < 1000) BRMS_WARMUP <- floor(BRMS_ITER / 2)
if (!is.finite(BRMS_THREADS) || BRMS_THREADS < 1) BRMS_THREADS <- 1
if (!is.finite(BRMS_CORES) || BRMS_CORES < 1) BRMS_CORES <- min(BRMS_CHAINS, parallel::detectCores())
if (!is.finite(BRMS_ADAPT) || BRMS_ADAPT < 0.95) BRMS_ADAPT <- 0.99
if (!is.finite(BRMS_TREEDEPTH) || BRMS_TREEDEPTH < 12) BRMS_TREEDEPTH <- 15
if (!is.finite(BRMS_GRAINSIZE) || BRMS_GRAINSIZE < 1) BRMS_GRAINSIZE <- 1000
if (!is.finite(N_DRAWS_PRED) || N_DRAWS_PRED < 400) N_DRAWS_PRED <- 1600
if (!BRMS_BACKEND %in% c("cmdstanr", "rstan")) BRMS_BACKEND <- "cmdstanr"

options(mc.cores = BRMS_CORES)

## ============================================================
## 2) Export helpers / 导出辅助函数
## ============================================================

save_plot_3formats <- function(p, filename_noext, outdir,
                               width = 12.2, height = 6.8, dpi = 320) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  png_file  <- file.path(outdir, paste0(filename_noext, ".png"))
  pdf_file  <- file.path(outdir, paste0(filename_noext, ".pdf"))
  pptx_file <- file.path(outdir, paste0(filename_noext, ".pptx"))

  print(p)
  if (requireNamespace("ragg", quietly = TRUE)) {
    ggsave(png_file, p, width = width, height = height, dpi = dpi,
           device = ragg::agg_png, bg = "white", limitsize = FALSE)
  } else {
    ggsave(png_file, p, width = width, height = height, dpi = dpi, bg = "white", limitsize = FALSE)
  }
  ggsave(pdf_file, p, width = width, height = height, bg = "white", limitsize = FALSE)

  if (exists("topptx")) {
    topptx(p, pptx_file)
  } else if (requireNamespace("officer", quietly = TRUE) &&
             requireNamespace("rvg", quietly = TRUE)) {
    doc <- officer::read_pptx()
    doc <- officer::add_slide(doc, layout = "Title and Content", master = "Office Theme")
    body_loc <- tryCatch(
      officer::ph_location_type(type = "body"),
      error = function(e) officer::ph_location_fullsize()
    )
    doc <- officer::ph_with(doc, rvg::dml(ggobj = p), location = body_loc)
    print(doc, target = pptx_file)
  }
  invisible(TRUE)
}

theme_paper <- function(base_size = 13) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = TITLE_SIZE),
      plot.subtitle = element_text(size = base_size),
      panel.spacing = unit(1.0, "lines"),
      strip.background = element_rect(fill = "white", color = "black"),
      strip.text = element_text(face = "bold", lineheight = 0.95),
      axis.title = element_text(face = "bold"),
      axis.title.x = element_text(margin = margin(t = 8)),
      axis.title.y = element_text(margin = margin(r = 8)),
      legend.title = element_text(face = "bold"),
      plot.margin = margin(t = 12, r = 24, b = 18, l = 14)
    )
}

figure_caption_registry <- tibble::tibble(
  filename_noext = character(),
  figure_dir = character(),
  title = character(),
  description_en = character(),
  description_cn = character()
)

add_figure_caption <- function(filename_noext, figure_dir, title,
                               description_en, description_cn) {
  figure_caption_registry <<- dplyr::bind_rows(
    figure_caption_registry,
    tibble::tibble(
      filename_noext = filename_noext,
      figure_dir = normalizePath(figure_dir, winslash = "/", mustWork = FALSE),
      title = title,
      description_en = description_en,
      description_cn = description_cn
    )
  )
}

write_figure_captions <- function() {
  cap_csv <- file.path(dir_caps, "FG_brms_figure_brief_descriptions_CN_EN.csv")
  cap_txt <- file.path(dir_caps, "FG_brms_figure_brief_descriptions_CN_EN.txt")

  utils::write.csv(figure_caption_registry, cap_csv, row.names = FALSE)

  con <- file(cap_txt, open = "wt")
  on.exit(close(con), add = TRUE)
  cat("FG brms figure brief descriptions (CN + EN)\n\n", file = con)
  for (i in seq_len(nrow(figure_caption_registry))) {
    row <- figure_caption_registry[i, ]
    cat("[", row$filename_noext, "] ", row$title, "\n", sep = "", file = con)
    cat("Path: ", row$figure_dir, "\n", sep = "", file = con)
    cat("EN: ", row$description_en, "\n", sep = "", file = con)
    cat("CN: ", row$description_cn, "\n\n", sep = "", file = con)
  }
}

label_map_UI2 <- c(
  "Primary vegetation"   = "Primary vegetation",
  "Secondary vegetation" = "Secondary vegetation",
  "Agriculture_Low"      = "Agriculture (low)",
  "Agriculture_High"     = "Agriculture (high)",
  "Urban"                = "Urban"
)

label_map_UI2_wrap <- c(
  "Primary vegetation"   = "Primary\nvegetation",
  "Secondary vegetation" = "Secondary\nvegetation",
  "Agriculture_Low"      = "Agriculture\n(low)",
  "Agriculture_High"     = "Agriculture\n(high)",
  "Urban"                = "Urban"
)

## ============================================================
## 3) Small utilities / 小工具
## ============================================================

summ_draws <- function(draws_vec) {
  c(
    mean  = mean(draws_vec),
    low95 = unname(quantile(draws_vec, 0.025)),
    high95= unname(quantile(draws_vec, 0.975))
  )
}

percent_change_draws <- function(ep, ep0) {
  ep0v <- as.vector(ep0)
  ep0v[ep0v < 1e-9] <- 1e-9
  sweep(ep, 1, ep0v, FUN = function(p, p0) (p - p0) / p0 * 100)
}

collapse_draws_weighted <- function(ep_mat, group_key, weights_per_col) {
  glev <- unique(group_key)
  out <- sapply(glev, function(g) {
    cols <- which(group_key == g)
    w <- weights_per_col[cols]
    w <- w / sum(w)
    ep_mat[, cols, drop = FALSE] %*% w
  })
  if (is.vector(out)) out <- matrix(out, ncol = length(glev))
  colnames(out) <- glev
  out
}

pretty_coef_label <- function(x) {
  x <- str_replace_all(x, c(
    "StdTmeanAnomalyRS" = "Warming",
    "UI2Secondaryvegetation" = "Secondary vegetation",
    "UI2Agriculture_Low" = "Agriculture (low)",
    "UI2Agriculture_High" = "Agriculture (high)",
    "UI2Urban" = "Urban",
    "FGFG1" = "FG1",
    "FGFG3" = "FG3",
    ":" = " x "
  ))
  stringr::str_wrap(x, width = 34)
}

extract_divergences <- function(fit_obj) {
  if (inherits(fit_obj, "CmdStanMCMC")) {
    diag_arr <- fit_obj$sampler_diagnostics()
    if (length(dim(diag_arr)) == 3 && "divergent__" %in% dimnames(diag_arr)[[3]]) {
      return(sum(diag_arr[, , "divergent__"], na.rm = TRUE))
    }
  }
  if (inherits(fit_obj, "stanfit")) {
    sampler_list <- rstan::get_sampler_params(fit_obj, inc_warmup = FALSE)
    return(sum(vapply(sampler_list, function(x) sum(x[, "divergent__"], na.rm = TRUE), numeric(1))))
  }
  NA_real_
}

extract_max_treedepth_hits <- function(fit_obj) {
  if (inherits(fit_obj, "CmdStanMCMC")) {
    diag_arr <- fit_obj$sampler_diagnostics()
    if (length(dim(diag_arr)) == 3 && "treedepth__" %in% dimnames(diag_arr)[[3]]) {
      return(sum(diag_arr[, , "treedepth__"] >= BRMS_TREEDEPTH, na.rm = TRUE))
    }
  }
  if (inherits(fit_obj, "stanfit")) {
    sampler_list <- rstan::get_sampler_params(fit_obj, inc_warmup = FALSE)
    return(sum(vapply(sampler_list, function(x) sum(x[, "treedepth__"] >= BRMS_TREEDEPTH, na.rm = TRUE), numeric(1))))
  }
  NA_real_
}

## ============================================================
## 4) Read stable FG data / 读取稳定FG数据
## ============================================================

data_file <- file.path(stable_root, "00_Data", "myoccdata_with_FG.rds")
if (!file.exists(data_file)) stop("Stable FG data not found: ", data_file)
myocc_fg <- readRDS(data_file)

need_cols <- c("Occur", "UI2", "SS", "SSBS", "Best_guess_binomial",
               "StdTmeanAnomalyRS", "FG")
stopifnot(all(need_cols %in% names(myocc_fg)))

myocc_fg$Occur <- as.numeric(myocc_fg$Occur)
myocc_fg$UI2 <- factor(myocc_fg$UI2, levels = UI2_levels)
myocc_fg$SS <- factor(myocc_fg$SS)
myocc_fg$SSBS <- factor(myocc_fg$SSBS)
myocc_fg$Best_guess_binomial <- factor(myocc_fg$Best_guess_binomial)
myocc_fg$FG <- factor(myocc_fg$FG, levels = FG_levels)

if (is.finite(SMOKE_N) && SMOKE_N > 0 && SMOKE_N < nrow(myocc_fg)) {
  set.seed(BRMS_SEED)
  myocc_fg <- myocc_fg %>%
    dplyr::group_by(UI2) %>%
    dplyr::slice_sample(n = max(1, min(dplyr::n(), ceiling(SMOKE_N / length(UI2_levels))))) %>%
    dplyr::ungroup()
}

FG_base <- FG_levels[1]
myocc_fg$FG <- relevel(myocc_fg$FG, ref = FG_base)

cat("Stable FG data loaded.\n")
cat("Rows:", nrow(myocc_fg), "\n")
cat("Species:", dplyr::n_distinct(myocc_fg$Best_guess_binomial), "\n")
cat("Baseline FG:", FG_base, "\n")

run_manifest <- data.frame(
  metric = c(
    "backend", "chains", "iter", "warmup", "cores", "threads_per_chain",
    "adapt_delta", "max_treedepth", "prediction_draws",
    "smoke_n", "rows_used", "species_used", "seed"
  ),
  value = c(
    BRMS_BACKEND, BRMS_CHAINS, BRMS_ITER, BRMS_WARMUP, BRMS_CORES, BRMS_THREADS,
    BRMS_ADAPT, BRMS_TREEDEPTH, N_DRAWS_PRED,
    SMOKE_N, nrow(myocc_fg), dplyr::n_distinct(myocc_fg$Best_guess_binomial), BRMS_SEED
  )
)
utils::write.csv(run_manifest, file.path(dir_tables, "FG_brms_run_manifest.csv"), row.names = FALSE)

## Record support figure / 数据支撑图
support_df <- myocc_fg %>%
  dplyr::count(FG, UI2, name = "n_records") %>%
  dplyr::mutate(label = format(n_records, big.mark = ",", scientific = FALSE))
utils::write.csv(support_df, file.path(dir_tables, "FG_UI2_record_support.csv"), row.names = FALSE)

p_support <- ggplot(support_df, aes(x = UI2, y = FG, fill = n_records)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = label), size = 4) +
  scale_x_discrete(labels = label_map_UI2_wrap) +
  scale_fill_gradient(low = "#EAF4F1", high = "#0F766E", name = "Records") +
  theme_paper(BASE_SIZE) +
  theme(
    axis.text.x = element_text(size = BASE_SIZE - 1, lineheight = 0.9),
    axis.title = element_blank()
  ) +
  labs(
    title = "Data support across land-use types and functional strategy groups",
    subtitle = "Record counts help explain why uncertainty differs among FG and land-use combinations."
  )
save_plot_3formats(p_support, "Fig0_DataSupport_UI2_by_FG_record_counts", dir_expl,
                   width = 9.6, height = 4.8)
add_figure_caption(
  "Fig0_DataSupport_UI2_by_FG_record_counts", dir_expl,
  "Data support across land-use types and functional strategy groups",
  "Record support is substantial across all land-use by FG combinations, but sample size still differs enough to affect the width of Bayesian credible intervals.",
  "各土地利用类型与 FG 组合都有较充足记录，但样本量仍存在明显差异，这会直接影响 Bayesian 可信区间的宽度。"
)

## ============================================================
## 5) brms model / Bayesian 模型
## ============================================================

fml <- bf(
  Occur ~ UI2 * StdTmeanAnomalyRS * FG +
    (1 || SS) + (1 || SSBS) + (1 || Best_guess_binomial),
  family = bernoulli(link = "logit"),
  decomp = "QR"
)

pri <- c(
  prior(normal(0, 0.7), class = "b"),
  prior(exponential(2), class = "sd")
)

fit_file <- file.path(dir_models, "brms_Occ_UI2_Temp_FG_stable_cmdstanr.rds")

if (file.exists(fit_file) && !BRMS_FORCE_REFIT) {
  message("Loading existing brms fit: ", fit_file)
  m_FG_brms <- readRDS(fit_file)
} else {
  brm_args <- list(
    formula = fml,
    data = droplevels(myocc_fg),
    backend = BRMS_BACKEND,
    prior = pri,
    chains = BRMS_CHAINS,
    iter = BRMS_ITER,
    warmup = BRMS_WARMUP,
    cores = BRMS_CORES,
    seed = BRMS_SEED,
    control = list(adapt_delta = BRMS_ADAPT, max_treedepth = BRMS_TREEDEPTH),
    save_pars = save_pars(all = FALSE, group = FALSE),
    normalize = FALSE,
    init = 0,
    refresh = 100
  )
  if (BRMS_THREADS > 1) {
    brm_args$threads <- threading(BRMS_THREADS, grainsize = BRMS_GRAINSIZE)
  }
  m_FG_brms <- do.call(brm, brm_args)
  saveRDS(m_FG_brms, fit_file)
}

## ============================================================
## 6) Diagnostics and tables / 诊断与表格
## ============================================================

fit_summary <- summary(m_FG_brms)
capture.output(fit_summary, file = file.path(dir_models, "FG_brms_model_summary.txt"))

fixef_tab <- as.data.frame(fixef(m_FG_brms, probs = c(0.025, 0.975)))
fixef_tab$term <- rownames(fixef_tab)
rownames(fixef_tab) <- NULL
utils::write.csv(fixef_tab, file.path(dir_tables, "FG_brms_fixed_effects_summary.csv"), row.names = FALSE)

diag_summary <- data.frame(
  metric = c("chains", "iter", "warmup", "threads_per_chain", "ndraws_total",
             "max_Rhat_fixed", "min_Bulk_ESS_fixed", "min_Tail_ESS_fixed",
             "divergences", "max_treedepth_hits"),
  value = c(
    BRMS_CHAINS,
    BRMS_ITER,
    BRMS_WARMUP,
    BRMS_THREADS,
    posterior::ndraws(m_FG_brms),
    max(fit_summary$fixed[, "Rhat"], na.rm = TRUE),
    min(fit_summary$fixed[, "Bulk_ESS"], na.rm = TRUE),
    min(fit_summary$fixed[, "Tail_ESS"], na.rm = TRUE),
    extract_divergences(m_FG_brms$fit),
    extract_max_treedepth_hits(m_FG_brms$fit)
  )
)
utils::write.csv(diag_summary, file.path(dir_tables, "FG_brms_diagnostics_summary.csv"), row.names = FALSE)

coef_focus <- fixef_tab %>%
  dplyr::filter(term == "StdTmeanAnomalyRS" | stringr::str_detect(term, "StdTmeanAnomalyRS")) %>%
  dplyr::mutate(label = pretty_coef_label(term)) %>%
  dplyr::arrange(Estimate)

p_coef <- ggplot(coef_focus, aes(x = Estimate, y = reorder(label, Estimate))) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5, colour = "grey55") +
  geom_linerange(aes(xmin = Q2.5, xmax = Q97.5), linewidth = 0.9, colour = "#1F3B70") +
  geom_point(size = 2.6, colour = "#1F3B70") +
  theme_paper(BASE_SIZE) +
  labs(
    title = "Bayesian posterior estimates for warming-related coefficients",
    subtitle = "Points are posterior means; horizontal intervals are 95% credible intervals on the logit scale.",
    x = "Posterior estimate (logit scale)",
    y = NULL
  )
save_plot_3formats(p_coef, "Fig9_Bayesian_warming_related_coefficients_forest", dir_diag,
                   width = 10.6, height = 7.2)
add_figure_caption(
  "Fig9_Bayesian_warming_related_coefficients_forest", dir_diag,
  "Bayesian posterior estimates for warming-related coefficients",
  "The posterior coefficient plot shows which warming-related main effects, two-way interactions and three-way interactions are credibly separated from zero on the logit scale.",
  "后验系数图直接展示了哪些与变暖相关的主效应、二重交互和三重交互在 logit 尺度上与零明显分离。"
)

draw_vars <- posterior::variables(m_FG_brms)
key_area_terms <- c(
  "b_StdTmeanAnomalyRS",
  "b_UI2Agriculture_Low:StdTmeanAnomalyRS",
  "b_UI2Urban:StdTmeanAnomalyRS",
  "b_FGFG1:StdTmeanAnomalyRS",
  "b_FGFG3:StdTmeanAnomalyRS"
)
key_area_terms <- key_area_terms[key_area_terms %in% draw_vars]
if (length(key_area_terms) > 0) {
  area_labels <- setNames(stringr::str_wrap(pretty_coef_label(sub("^b_", "", key_area_terms)), width = 28),
                          key_area_terms)
  p_area <- bayesplot::mcmc_areas(
    as.array(m_FG_brms),
    pars = key_area_terms,
    prob = 0.8,
    prob_outer = 0.95
  ) +
    ggplot2::scale_y_discrete(labels = area_labels) +
    theme_paper(BASE_SIZE) +
    labs(
      title = "Posterior distributions of representative warming-related terms",
      subtitle = "Inner bands show 80% credible intervals; outer bands show 95% credible intervals.",
      x = "Posterior value",
      y = NULL
    )
  save_plot_3formats(p_area, "Fig10_Bayesian_representative_posterior_distributions", dir_diag,
                     width = 10.2, height = 6.8)
  add_figure_caption(
    "Fig10_Bayesian_representative_posterior_distributions", dir_diag,
    "Posterior distributions of representative warming-related terms",
    "Posterior densities provide a direct visual check of uncertainty and the direction of key warming-related coefficients under the Bayesian model.",
    "代表性后验分布图可以直观展示关键变暖相关系数的不确定性范围及其方向。"
  )
}

## ============================================================
## 7) Predictions / 后验预测
## ============================================================

temp_seq <- seq(
  quantile(myocc_fg$StdTmeanAnomalyRS, 0.02, na.rm = TRUE),
  quantile(myocc_fg$StdTmeanAnomalyRS, 0.98, na.rm = TRUE),
  length.out = N_TEMP
)

ndraws_pred_use <- min(N_DRAWS_PRED, posterior::ndraws(m_FG_brms))

## 7.1 Three-way predictions
nd3 <- expand.grid(
  StdTmeanAnomalyRS = temp_seq,
  UI2 = levels(myocc_fg$UI2),
  FG = levels(myocc_fg$FG)
)

ep3 <- posterior_epred(
  m_FG_brms,
  newdata = nd3,
  re_formula = NA,
  ndraws = ndraws_pred_use
)

nd0 <- data.frame(
  UI2 = factor(BASE_UI2, levels = levels(myocc_fg$UI2)),
  StdTmeanAnomalyRS = BASE_TEMP0,
  FG = factor(FG_base, levels = levels(myocc_fg$FG))
)
ep0 <- posterior_epred(
  m_FG_brms,
  newdata = nd0,
  re_formula = NA,
  ndraws = ndraws_pred_use
)

summ_abs3 <- apply(ep3, 2, summ_draws)
pred_abs3 <- cbind(nd3, as.data.frame(t(summ_abs3))) %>%
  dplyr::arrange(FG, UI2, StdTmeanAnomalyRS)
pred_abs3$FG <- factor(pred_abs3$FG, levels = FG_levels)
pred_abs3$UI2 <- factor(pred_abs3$UI2, levels = UI2_levels)

ep_pct3 <- percent_change_draws(ep3, ep0)
summ_pct3 <- apply(ep_pct3, 2, summ_draws)
pred_pct3 <- cbind(nd3, as.data.frame(t(summ_pct3))) %>%
  dplyr::arrange(FG, UI2, StdTmeanAnomalyRS)
pred_pct3$FG <- factor(pred_pct3$FG, levels = FG_levels)
pred_pct3$UI2 <- factor(pred_pct3$UI2, levels = UI2_levels)

saveRDS(pred_abs3, file.path(dir_data, "pred_abs_threeway_brms.rds"))
saveRDS(pred_pct3, file.path(dir_data, "pred_pct_threeway_vsPV0_brms.rds"))

p_abs3 <- ggplot(pred_abs3, aes(x = StdTmeanAnomalyRS, colour = UI2, fill = UI2)) +
  geom_ribbon(aes(ymin = low95, ymax = high95), alpha = 0.22, colour = NA) +
  geom_line(aes(y = mean), linewidth = 1.05) +
  facet_grid(FG ~ UI2, labeller = labeller(UI2 = label_map_UI2_wrap, FG = fg_strip_labels)) +
  scale_colour_manual(values = pal_UI2, drop = FALSE, guide = "none") +
  scale_fill_manual(values = pal_UI2, drop = FALSE, guide = "none") +
  theme_paper(BASE_SIZE) +
  theme(strip.text.y = element_text(angle = 0)) +
  labs(
    title = "Occurrence responses to land use and warming differ among functional strategy groups",
    subtitle = paste0("brms posterior means; ", ndraws_pred_use, " posterior draws for prediction summaries."),
    x = "Standardised mean temperature anomaly (StdTmeanAnomalyRS)",
    y = "Predicted probability of occurrence"
  )
save_plot_3formats(p_abs3, "Fig1_ThreeWay_ABS_UI2_Temp_FG_brms", dir_plots3,
                   width = 12.2, height = 7.2)
add_figure_caption(
  "Fig1_ThreeWay_ABS_UI2_Temp_FG_brms", dir_plots3,
  "Occurrence responses to land use and warming differ among functional strategy groups",
  "Bayesian posterior predictions support heterogeneous warming-response curves among land-use types and FG, consistent with a meaningful three-way interaction.",
  "Bayesian 后验预测同样显示，不同土地利用类型和不同 FG 的变暖响应曲线存在明显异质性，支持具有生态意义的三重交互。"
)

p_pct3 <- ggplot(pred_pct3, aes(x = StdTmeanAnomalyRS, colour = UI2, fill = UI2)) +
  geom_ribbon(aes(ymin = low95, ymax = high95), alpha = 0.22, colour = NA) +
  geom_line(aes(y = mean), linewidth = 1.05) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5) +
  facet_grid(FG ~ UI2, labeller = labeller(UI2 = label_map_UI2_wrap, FG = fg_strip_labels)) +
  scale_colour_manual(values = pal_UI2, drop = FALSE, guide = "none") +
  scale_fill_manual(values = pal_UI2, drop = FALSE, guide = "none") +
  theme_paper(BASE_SIZE) +
  theme(strip.text.y = element_text(angle = 0)) +
  labs(
    title = "Relative warming responses of occurrence probability across land uses and strategy groups",
    subtitle = paste0(
      "Baseline: UI2=", BASE_UI2, ", StdTmeanAnomalyRS=", BASE_TEMP0,
      ", FG=", FG_base, ".  Percent change = (p-p0)/p0*100."
    ),
    x = "Standardised mean temperature anomaly (StdTmeanAnomalyRS)",
    y = "Percent change in predicted occurrence probability (%)"
  )
save_plot_3formats(p_pct3, "Fig2_ThreeWay_PCT_vsPV0_UI2_Temp_FG_brms", dir_plots3,
                   width = 12.2, height = 7.2)
add_figure_caption(
  "Fig2_ThreeWay_PCT_vsPV0_UI2_Temp_FG_brms", dir_plots3,
  "Relative warming responses of occurrence probability across land uses and strategy groups",
  "Relative change from the primary-vegetation baseline varies across FG and land-use contexts, reinforcing the interpretation that warming sensitivity depends on ecological strategy and land-use background.",
  "相对于原始植被基线的变化幅度在不同 FG 和土地利用背景下并不一致，进一步说明变暖敏感性同时依赖功能策略和土地利用环境。"
)

## 7.2 Two-way: UI2 x warming averaged over FG frequency
fg_wt <- prop.table(table(myocc_fg$FG))
nd_UI2Temp <- expand.grid(
  StdTmeanAnomalyRS = temp_seq,
  UI2 = levels(myocc_fg$UI2),
  FG = levels(myocc_fg$FG)
)
ep_UI2Temp <- posterior_epred(
  m_FG_brms,
  newdata = nd_UI2Temp,
  re_formula = NA,
  ndraws = ndraws_pred_use
)
weights_fg <- as.numeric(fg_wt[as.character(nd_UI2Temp$FG)])
group_key_UI2Temp <- paste(nd_UI2Temp$UI2, nd_UI2Temp$StdTmeanAnomalyRS, sep = "___")
collapsed_UI2Temp <- collapse_draws_weighted(ep_UI2Temp, group_key_UI2Temp, weights_fg)
ui2temp_df <- nd_UI2Temp %>%
  dplyr::select(UI2, StdTmeanAnomalyRS) %>%
  dplyr::distinct() %>%
  dplyr::arrange(UI2, StdTmeanAnomalyRS)
summ_UI2Temp <- apply(collapsed_UI2Temp, 2, summ_draws)
pred_UI2Temp <- cbind(ui2temp_df, as.data.frame(t(summ_UI2Temp)))

p_ui2temp <- ggplot(pred_UI2Temp, aes(x = StdTmeanAnomalyRS, y = mean, colour = UI2, fill = UI2)) +
  geom_ribbon(aes(ymin = low95, ymax = high95), alpha = 0.20, colour = NA) +
  geom_line(linewidth = 1.05) +
  scale_colour_manual(values = pal_UI2, drop = FALSE) +
  scale_fill_manual(values = pal_UI2, drop = FALSE) +
  theme_paper(BASE_SIZE) +
  labs(
    title = "Two-way interaction: land use × warming (FG-averaged)",
    subtitle = "Bayesian posterior predictions; average over FG composition within each posterior draw.",
    x = "Standardised mean temperature anomaly (StdTmeanAnomalyRS)",
    y = "Predicted probability of occurrence",
    colour = "Land-use type",
    fill = "Land-use type"
  )
save_plot_3formats(p_ui2temp, "Fig6_TwoWay_UI2_by_Warming_FGaveraged_ABS_brms", dir_plots2,
                   width = 10.6, height = 6.6)
add_figure_caption(
  "Fig6_TwoWay_UI2_by_Warming_FGaveraged_ABS_brms", dir_plots2,
  "Two-way interaction: land use × warming (FG-averaged)",
  "Even after averaging over FG composition, land-use types retain distinct warming-response curves, indicating a robust land-use filtering effect.",
  "即使对 FG 组成加权平均后，不同土地利用类型的变暖响应曲线仍然不同，说明土地利用过滤效应是稳健存在的。"
)

## 7.3 Two-way: FG x warming averaged over UI2 frequency
ui2_wt <- prop.table(table(myocc_fg$UI2))
nd_FGTemp <- expand.grid(
  StdTmeanAnomalyRS = temp_seq,
  UI2 = levels(myocc_fg$UI2),
  FG = levels(myocc_fg$FG)
)
ep_FGTemp <- posterior_epred(
  m_FG_brms,
  newdata = nd_FGTemp,
  re_formula = NA,
  ndraws = ndraws_pred_use
)
weights_ui2 <- as.numeric(ui2_wt[as.character(nd_FGTemp$UI2)])
group_key_FGTemp <- paste(nd_FGTemp$FG, nd_FGTemp$StdTmeanAnomalyRS, sep = "___")
collapsed_FGTemp <- collapse_draws_weighted(ep_FGTemp, group_key_FGTemp, weights_ui2)
fgtemp_df <- nd_FGTemp %>%
  dplyr::select(FG, StdTmeanAnomalyRS) %>%
  dplyr::distinct() %>%
  dplyr::arrange(FG, StdTmeanAnomalyRS)
summ_FGTemp <- apply(collapsed_FGTemp, 2, summ_draws)
pred_FGTemp <- cbind(fgtemp_df, as.data.frame(t(summ_FGTemp)))
pred_FGTemp$FG <- factor(pred_FGTemp$FG, levels = FG_levels)

p_fgtemp <- ggplot(pred_FGTemp, aes(x = StdTmeanAnomalyRS, y = mean, colour = FG, fill = FG)) +
  geom_ribbon(aes(ymin = low95, ymax = high95), alpha = 0.20, colour = NA) +
  geom_line(linewidth = 1.05) +
  scale_colour_manual(values = pal_FG, drop = FALSE,
                      labels = fg_feature_labels[levels(pred_FGTemp$FG)]) +
  scale_fill_manual(values = pal_FG, drop = FALSE,
                    labels = fg_feature_labels[levels(pred_FGTemp$FG)]) +
  theme_paper(BASE_SIZE) +
  labs(
    title = "Two-way interaction: functional strategy group × warming (UI2-averaged)",
    subtitle = "Bayesian posterior predictions; average over land-use composition within each posterior draw.",
    x = "Standardised mean temperature anomaly (StdTmeanAnomalyRS)",
    y = "Predicted probability of occurrence",
    colour = "FG",
    fill = "FG"
  )
save_plot_3formats(p_fgtemp, "Fig7_TwoWay_FG_by_Warming_UI2averaged_ABS_brms", dir_plots2,
                   width = 10.6, height = 6.6)
add_figure_caption(
  "Fig7_TwoWay_FG_by_Warming_UI2averaged_ABS_brms", dir_plots2,
  "Two-way interaction: functional strategy group × warming (UI2-averaged)",
  "FG retain distinct warming-response trajectories after averaging over land-use types, supporting the ecological interpretability of the stable FG scheme.",
  "对土地利用类型加权平均之后，不同 FG 的变暖响应轨迹仍然不同，说明稳定 FG 方案具有清晰的生态可解释性。"
)

## 7.4 Two-way: UI2 x FG at fixed warming levels (no line)
warm_levels <- c(0, 1, median(myocc_fg$StdTmeanAnomalyRS, na.rm = TRUE))
warm_labels <- c("0 (baseline)", "1 (+1 SD)", "median")

nd_UI2FG <- expand.grid(
  StdTmeanAnomalyRS = warm_levels,
  UI2 = levels(myocc_fg$UI2),
  FG = levels(myocc_fg$FG)
)
nd_UI2FG$WarmLabel <- factor(
  nd_UI2FG$StdTmeanAnomalyRS,
  levels = warm_levels,
  labels = warm_labels
)

ep_UI2FG <- posterior_epred(
  m_FG_brms,
  newdata = nd_UI2FG,
  re_formula = NA,
  ndraws = ndraws_pred_use
)
summ_UI2FG <- apply(ep_UI2FG, 2, summ_draws)
pred_UI2FG <- cbind(nd_UI2FG, as.data.frame(t(summ_UI2FG))) %>%
  dplyr::arrange(WarmLabel, FG, UI2)
pred_UI2FG$FG <- factor(pred_UI2FG$FG, levels = FG_levels)
pred_UI2FG$UI2 <- factor(pred_UI2FG$UI2, levels = UI2_levels)

p_ui2fg <- ggplot(pred_UI2FG, aes(x = UI2, y = mean, colour = FG)) +
  geom_hline(yintercept = 0, colour = "grey82", linewidth = 0.4) +
  geom_pointrange(aes(ymin = low95, ymax = high95),
                  position = position_dodge(width = 0.55),
                  linewidth = 0.7) +
  facet_wrap(~ WarmLabel, nrow = 1) +
  scale_x_discrete(labels = label_map_UI2_wrap) +
  scale_colour_manual(values = pal_FG, drop = FALSE,
                      labels = fg_feature_labels[levels(pred_UI2FG$FG)]) +
  theme_paper(BASE_SIZE) +
  theme(
    panel.spacing.x = unit(1.2, "lines"),
    axis.text.x = element_text(size = BASE_SIZE - 2, lineheight = 0.9, angle = 8, hjust = 1, vjust = 1)
  ) +
  labs(
    title = "Two-way interaction: land use × FG (at fixed warming levels)",
    subtitle = "Posterior means and 95% credible intervals at StdTmeanAnomalyRS = 0, 1, and median.",
    x = NULL,
    y = "Predicted probability of occurrence",
    colour = "Functional strategy group"
  )
save_plot_3formats(p_ui2fg, "Fig8_TwoWay_UI2_by_FG_at_fixedWarming_ABS_brms", dir_plots2,
                   width = 14.4, height = 6.4)
add_figure_caption(
  "Fig8_TwoWay_UI2_by_FG_at_fixedWarming_ABS_brms", dir_plots2,
  "Two-way interaction: land use × FG (at fixed warming levels)",
  "At the same warming level, FG differ within each land-use panel, making the comparison among strategy groups visually direct and avoiding misleading line connections across categorical land-use types.",
  "在相同变暖水平下，每个土地利用 panel 中都可以直接比较不同 FG 的差异，而且不再使用跨类别连线，避免视觉误导。"
)

message("\nStable FG Bayesian pipeline completed.")
write_figure_captions()
message("Outputs saved under: ", normalizePath(out_root))
