###############################################################
## Community diversity main rewrite
## 群落多样性主线重整版：LMM + Negative binomial richness
###############################################################

script_dir <- normalizePath(dirname(sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grepl("^--file=", commandArgs(trailingOnly = FALSE))][1])))
source(file.path(script_dir, "00_shared_helpers.R"))
suppressPackageStartupMessages({
  library(lme4)
  library(glmmTMB)
  library(MASS)
  library(DHARMa)
})

project_root <- get_project_root(script_dir)
out_root <- file.path(project_root, "results", "community", "diversity_main_nb")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_root, "01_DataSupport"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_root, "02_Models"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_root, "03_Plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_root, "04_Tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_root, "05_Figure_Captions"), recursive = TRUE, showWarnings = FALSE)

captions <- list()

data_candidates <- c(
  file.path(project_root, "data", "derived_public", "predicts_sites_used.rds")
)
data_file <- data_candidates[file.exists(data_candidates)][1]
if (is.na(data_file)) stop("No processed site-level PREDICTS data found.")

dat <- readRDS(data_file)
dat$UI2 <- factor(dat$UI2, levels = ui2_levels)

support_mean <- dat %>% filter(!is.na(StdTmeanAnomalyRS), !is.na(UI2))
support_max <- dat %>% filter(!is.na(StdTmaxAnomalyRS), !is.na(UI2))

write.csv(
  support_mean %>% count(UI2, name = "n_sites_mean"),
  file.path(out_root, "04_Tables", "support_mean_by_UI2.csv"),
  row.names = FALSE
)
write.csv(
  support_max %>% count(UI2, name = "n_sites_max"),
  file.path(out_root, "04_Tables", "support_max_by_UI2.csv"),
  row.names = FALSE
)

p_support_mean <- build_ui2_support_plot(
  support_mean,
  "StdTmeanAnomalyRS",
  "Standardised mean temperature anomaly (StdTmeanAnomalyRS)",
  "Site support for mean temperature anomaly models"
)
save_plot_3formats(
  p_support_mean,
  "Fig1_Diversity_support_mean_anomaly",
  file.path(out_root, "01_DataSupport"),
  width_mm = 183,
  height_mm = 105,
  slide_title = "Diversity mean-anomaly support"
)
captions <- record_caption(
  captions, "Fig1", "Fig1_Diversity_support_mean_anomaly",
  "均温异常模型的站点样本支撑分布",
  "Site support for the mean-anomaly diversity models",
  "不同土地利用类型在均温异常上的样本支撑并不完全一致，因此模型曲线应结合经验范围一起解释。",
  "The empirical support differs across land-use types, so model curves should be interpreted together with the observed support."
)

p_support_max <- build_ui2_support_plot(
  support_max,
  "StdTmaxAnomalyRS",
  "Standardised maximum temperature anomaly (StdTmaxAnomalyRS)",
  "Site support for maximum temperature anomaly models"
)
save_plot_3formats(
  p_support_max,
  "Fig2_Diversity_support_max_anomaly",
  file.path(out_root, "01_DataSupport"),
  width_mm = 183,
  height_mm = 105,
  slide_title = "Diversity max-anomaly support"
)
captions <- record_caption(
  captions, "Fig2", "Fig2_Diversity_support_max_anomaly",
  "极端高温异常模型的站点样本支撑分布",
  "Site support for the maximum-anomaly diversity models",
  "极端温度异常的经验支撑范围更宽但也更不均匀，这会影响不同土地利用曲线末端的不确定性。",
  "Support for maximum anomalies is broader but also more uneven, which affects uncertainty near the ends of the curves."
)

fit_or_load <- function(path, expr) {
  if (file.exists(path)) return(readRDS(path))
  mod <- eval.parent(substitute(expr))
  saveRDS(mod, path)
  mod
}

dat_ab_mean <- dat %>% filter(!is.na(LogAbund), !is.na(StdTmeanAnomalyRS))
dat_ab_max <- dat %>% filter(!is.na(LogAbund), !is.na(StdTmaxAnomalyRS))
dat_rich_mean <- dat %>% filter(!is.na(Species_richness), !is.na(StdTmeanAnomalyRS))
dat_rich_max <- dat %>% filter(!is.na(Species_richness), !is.na(StdTmaxAnomalyRS))

m_ab_mean <- fit_or_load(
  file.path(out_root, "02_Models", "m_abundance_mean_lmer.rds"),
  lmer(LogAbund ~ UI2 * StdTmeanAnomalyRS + (1 | SS) + (1 | SSB), data = dat_ab_mean, REML = FALSE)
)
m_ab_max <- fit_or_load(
  file.path(out_root, "02_Models", "m_abundance_max_lmer.rds"),
  lmer(LogAbund ~ UI2 * StdTmaxAnomalyRS + (1 | SS) + (1 | SSB), data = dat_ab_max, REML = FALSE)
)
m_rich_mean <- fit_or_load(
  file.path(out_root, "02_Models", "m_richness_mean_nb.rds"),
  glmmTMB(
    Species_richness ~ UI2 * StdTmeanAnomalyRS + (1 | SS) + (1 | SSB) + (1 | SSBS),
    data = dat_rich_mean,
    family = nbinom2(link = "log")
  )
)
m_rich_max <- fit_or_load(
  file.path(out_root, "02_Models", "m_richness_max_nb.rds"),
  glmmTMB(
    Species_richness ~ UI2 * StdTmaxAnomalyRS + (1 | SS) + (1 | SSB) + (1 | SSBS),
    data = dat_rich_max,
    family = nbinom2(link = "log")
  )
)

capture.output(summary(m_ab_mean), file = file.path(out_root, "04_Tables", "summary_abundance_mean_lmer.txt"))
capture.output(summary(m_ab_max), file = file.path(out_root, "04_Tables", "summary_abundance_max_lmer.txt"))
capture.output(summary(m_rich_mean), file = file.path(out_root, "04_Tables", "summary_richness_mean_nb.txt"))
capture.output(summary(m_rich_max), file = file.path(out_root, "04_Tables", "summary_richness_max_nb.txt"))

diag_mean <- simulateResiduals(m_rich_mean, n = 500)
diag_max <- simulateResiduals(m_rich_max, n = 500)
sink(file.path(out_root, "04_Tables", "nb_richness_dharma_tests.txt"))
cat("Mean anomaly richness NB\n")
print(testDispersion(diag_mean))
print(testZeroInflation(diag_mean))
print(testOutliers(diag_mean, type = "bootstrap"))
cat("\nMax anomaly richness NB\n")
print(testDispersion(diag_max))
print(testZeroInflation(diag_max))
print(testOutliers(diag_max, type = "bootstrap"))
sink()

coef_div <- bind_rows(
  coef_table_lmer(m_ab_mean) %>% mutate(model = "Abundance_mean"),
  coef_table_lmer(m_ab_max) %>% mutate(model = "Abundance_max"),
  coef_table_glmmtmb(m_rich_mean) %>% mutate(model = "Richness_mean_NB"),
  coef_table_glmmtmb(m_rich_max) %>% mutate(model = "Richness_max_NB")
)
write.csv(coef_div, file.path(out_root, "04_Tables", "diversity_fixed_effects_summary.csv"), row.names = FALSE)

coef_forest_df <- coef_div %>%
  filter(grepl("StdT", term)) %>%
  mutate(term_label = paste(model, term, sep = " | "))

p_coef <- plot_fixed_effect_forest(
  coef_forest_df,
  term_col = "term_label",
  title = "Warming-related coefficients in the rewritten diversity models",
  subtitle = "Abundance is modelled with LMM and richness with negative-binomial mixed models.",
  xlab = "Fixed-effect estimate"
)
save_plot_3formats(
  p_coef,
  "Fig3_Diversity_warming_coefficients",
  file.path(out_root, "03_Plots"),
  width_mm = 183,
  height_mm = 140,
  slide_title = "Diversity warming coefficients"
)
captions <- record_caption(
  captions, "Fig3", "Fig3_Diversity_warming_coefficients",
  "重写后的多样性模型中，变暖相关系数在 abundance 和 richness 之间并不完全一致",
  "Warming-related coefficients differ between abundance and richness in the rewritten diversity models",
  "把 abundance 和 richness 放在同一系数森林图中，有助于直接比较两类多样性指标对均温和极端高温异常的响应方向。",
  "Placing abundance and richness on the same coefficient forest makes it easier to compare their responses to mean and maximum warming anomalies."
)

predict_lmer_draws <- function(mod, nd, n_draw = 2000, seed = 123) {
  set.seed(seed)
  beta_hat <- lme4::fixef(mod)
  V <- as.matrix(vcov(mod))
  X <- model.matrix(delete.response(terms(mod)), nd)
  beta_draw <- MASS::mvrnorm(n = n_draw, mu = beta_hat, Sigma = V)
  X %*% t(beta_draw)
}

predict_glmmtmb_draws <- function(mod, nd, n_draw = 2000, seed = 123) {
  set.seed(seed)
  beta_hat <- glmmTMB::fixef(mod)$cond
  V <- as.matrix(vcov(mod)$cond)
  X <- model.matrix(lme4::nobars(formula(mod)[-2]), nd)
  beta_draw <- MASS::mvrnorm(n = n_draw, mu = beta_hat, Sigma = V)
  X %*% t(beta_draw)
}

make_percent_change <- function(model, dat_model, x_var, response_type = c("abundance", "richness"),
                                x_lab, fig_stub, fig_title, caption_cn, caption_en) {
  response_type <- match.arg(response_type)
  x_seq <- seq(min(dat_model[[x_var]], na.rm = TRUE), max(dat_model[[x_var]], na.rm = TRUE), length.out = 220)
  nd <- expand.grid(
    tmp = x_seq,
    UI2 = factor(ui2_levels, levels = ui2_levels)
  )
  names(nd)[1] <- x_var
  nd$SS <- dat_model$SS[1]
  nd$SSB <- dat_model$SSB[1]
  if ("SSBS" %in% names(dat_model)) nd$SSBS <- dat_model$SSBS[1]

  if (inherits(model, "lmerMod")) {
    eta_draw <- predict_lmer_draws(model, nd)
    mu_draw <- exp(eta_draw)
  } else {
    eta_draw <- predict_glmmtmb_draws(model, nd)
    mu_draw <- exp(eta_draw)
  }

  pv_rows <- which(nd$UI2 == "Primary vegetation")
  ref_row <- pv_rows[which.min(abs(nd[[x_var]][pv_rows]))]
  mu_rel <- sweep(mu_draw, 2, mu_draw[ref_row, ], "/")

  q_list <- lapply(ui2_levels, function(u) {
    quantile(dat_model[[x_var]][dat_model$UI2 == u], probs = c(0.025, 0.975), na.rm = TRUE)
  })
  names(q_list) <- ui2_levels
  for (u in ui2_levels) {
    q <- q_list[[u]]
    idx <- nd$UI2 == u & (nd[[x_var]] < q[1] | nd[[x_var]] > q[2])
    mu_rel[idx, ] <- NA
  }

  out <- nd %>%
    mutate(
      PredMedian = apply(mu_rel, 1, median, na.rm = TRUE) * 100 - 100,
      PredLower = apply(mu_rel, 1, quantile, probs = 0.025, na.rm = TRUE) * 100 - 100,
      PredUpper = apply(mu_rel, 1, quantile, probs = 0.975, na.rm = TRUE) * 100 - 100
    )
  out_plot <- out %>%
    filter(is.finite(PredMedian), is.finite(PredLower), is.finite(PredUpper))

  p <- ggplot(out_plot, aes(x = .data[[x_var]], y = PredMedian, colour = UI2, fill = UI2)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
    geom_ribbon(aes(ymin = PredLower, ymax = PredUpper), alpha = 0.18, linewidth = 0) +
    geom_line(linewidth = 0.7) +
    scale_colour_manual(values = ui2_cols, labels = label_map_ui2) +
    scale_fill_manual(values = ui2_cols, labels = label_map_ui2) +
    labs(
      title = fig_title,
      subtitle = "Percent change is expressed relative to Primary vegetation at anomaly = 0.",
      x = x_lab,
      y = if (response_type == "abundance") "Change in total abundance (%)" else "Change in species richness (%)"
    ) +
    theme_paper() +
    theme(legend.position = "bottom")

  save_plot_3formats(
    p,
    fig_stub,
    file.path(out_root, "03_Plots"),
    width_mm = 183,
    height_mm = 110,
    slide_title = fig_title
  )

  write.csv(out, file.path(out_root, "04_Tables", paste0(fig_stub, "_predictions.csv")), row.names = FALSE)

  captions <<- record_caption(
    captions,
    gsub("_.*$", "", fig_stub),
    fig_stub,
    caption_cn,
    fig_title,
    "各土地利用类型的曲线均以原生植被在异常值为 0 处作为统一基线，因此更适合横向比较。",
    "All land-use curves are expressed relative to the same Primary-vegetation baseline at anomaly = 0, making cross-land-use comparisons clearer."
  )
}

make_percent_change(
  m_ab_mean, dat_ab_mean, "StdTmeanAnomalyRS", "abundance",
  "Standardised mean temperature anomaly (StdTmeanAnomalyRS)",
  "Fig4_Diversity_abundance_mean_percent_change",
  "Abundance responses to mean warming anomaly across land-use types",
  "不同土地利用类型下，总丰度对均温异常的响应方向和幅度并不一致",
  "Abundance responses to mean warming anomaly across land-use types"
)
make_percent_change(
  m_rich_mean, dat_rich_mean, "StdTmeanAnomalyRS", "richness",
  "Standardised mean temperature anomaly (StdTmeanAnomalyRS)",
  "Fig5_Diversity_richness_mean_percent_change",
  "Richness responses to mean warming anomaly across land-use types",
  "采用 NB 重写后，物种丰富度对均温异常的结果不再依赖存在过度离散问题的 Poisson 主模型",
  "Richness responses to mean warming anomaly across land-use types"
)
make_percent_change(
  m_ab_max, dat_ab_max, "StdTmaxAnomalyRS", "abundance",
  "Standardised maximum temperature anomaly (StdTmaxAnomalyRS)",
  "Fig6_Diversity_abundance_max_percent_change",
  "Abundance responses to maximum warming anomaly across land-use types",
  "极端高温异常下，总丰度对土地利用背景的敏感性更加明显",
  "Abundance responses to maximum warming anomaly across land-use types"
)
make_percent_change(
  m_rich_max, dat_rich_max, "StdTmaxAnomalyRS", "richness",
  "Standardised maximum temperature anomaly (StdTmaxAnomalyRS)",
  "Fig7_Diversity_richness_max_percent_change",
  "Richness responses to maximum warming anomaly across land-use types",
  "极端高温异常下的 richness-NB 主图适合作为替代原 Poisson richness 结果的主版本",
  "Richness responses to maximum warming anomaly across land-use types"
)

sink(file.path(out_root, "04_Tables", "analysis_notes.txt"))
cat("Canonical source script:\n")
cat("scripts/03_main_models/01_run_lu_climate_models_main.R\n\n")
cat("This rewrite replaces Poisson richness as the default main model with negative-binomial mixed models.\n")
cat("Important note: DHARMa still indicates residual dispersion departures, but the severe Poisson overdispersion problem is no longer the main issue.\n")
sink()

write_caption_registry(
  captions,
  file.path(out_root, "05_Figure_Captions"),
  "diversity_main_nb_figure_descriptions_CN_EN"
)

message("Finished: ", out_root)
