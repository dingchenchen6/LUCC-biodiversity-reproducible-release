###############################################################
## Functional community structure: main Bray-ATA rewrite
## 功能群落结构主线重整版
###############################################################

script_dir <- normalizePath(dirname(sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grepl("^--file=", commandArgs(trailingOnly = FALSE))][1])))
source(file.path(script_dir, "00_shared_helpers.R"))
suppressPackageStartupMessages({
  library(lme4)
  library(tidyr)
})

project_root <- get_project_root(script_dir)
base_dir <- file.path(project_root, "results", "community", "functional_main_bray")
out_root <- file.path(project_root, "results", "community", "functional_main_bray")

dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_root, "01_DataSupport"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_root, "02_Models"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_root, "03_Plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_root, "04_Tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_root, "05_Figure_Captions"), recursive = TRUE, showWarnings = FALSE)

captions <- list()

cs_data <- readRDS(file.path(project_root, "data", "derived_public", "functional_pairs_with_covariates.rds"))
pred_abs <- readRDS(file.path(project_root, "results", "community", "functional_main_bray", "source_data", "lmer_pred_abs_D_all.rds"))
pred_diff0 <- readRDS(file.path(project_root, "results", "community", "functional_main_bray", "source_data", "lmer_pred_diff_vsPV0_fixed_D_all.rds"))
ata <- readRDS(file.path(project_root, "results", "community", "functional_main_bray", "source_data", "site_trait_matrix_ATA01.rds"))

pred_abs$UI2 <- factor(pred_abs$UI2, levels = ui2_levels)
pred_diff0$UI2 <- factor(pred_diff0$UI2, levels = ui2_levels)

presence_diag <- data.frame(
  stat = c("min_presence_fraction", "median_presence_fraction", "max_presence_fraction", "overall_zero_fraction"),
  value = c(
    min(rowMeans(ata > 0, na.rm = TRUE), na.rm = TRUE),
    median(rowMeans(ata > 0, na.rm = TRUE), na.rm = TRUE),
    max(rowMeans(ata > 0, na.rm = TRUE), na.rm = TRUE),
    mean(as.matrix(ata) == 0, na.rm = TRUE)
  )
)
write.csv(presence_diag, file.path(out_root, "04_Tables", "functional_presence_diagnostic.csv"), row.names = FALSE)

p_support <- build_ui2_support_plot(
  cs_data,
  x_var = "pair_clim_mean",
  x_lab = "Pair-level mean temperature anomaly (pair_clim_mean)",
  title = "Functional pair support across the climate gradient",
  subtitle = "Support is shared with the compositional analysis, but interpretation here is in ATA-based functional space."
)
save_plot_3formats(
  p_support,
  "Fig1_Functional_pair_climate_support",
  file.path(out_root, "01_DataSupport"),
  width_mm = 183,
  height_mm = 105,
  slide_title = "Functional pair climate support"
)
captions <- record_caption(
  captions, "Fig1", "Fig1_Functional_pair_climate_support",
  "功能群落结构配对样本在气候梯度上的支撑与组成结构分析一致",
  "Empirical climate support for functional site pairs",
  "功能群落结构分析与组成结构分析共享同一套站点 pair，因此气候支撑差异同样需要在结果解释中保留。",
  "The functional analysis uses the same site-pair support as the compositional analysis, so the same climate-support constraints apply."
)

main_metrics <- c("Bray-Curtis", "Jaccard", "Sorensen")
main_components <- c("Total", "Turnover", "Gradient/Nestedness")
pred_abs_main <- pred_abs %>%
  filter(metric %in% main_metrics, component %in% main_components) %>%
  mutate(
    metric = factor(metric, levels = main_metrics),
    component = factor(component, levels = main_components)
  )

p_abs <- ggplot(pred_abs_main, aes(x = pair_clim_mean, y = PredMedian, colour = UI2, fill = UI2)) +
  geom_ribbon(aes(ymin = PredLower, ymax = PredUpper), alpha = 0.18, linewidth = 0) +
  geom_line(linewidth = 0.7) +
  scale_colour_manual(values = ui2_cols, labels = label_map_ui2) +
  scale_fill_manual(values = ui2_cols, labels = label_map_ui2) +
  facet_grid(component ~ metric, labeller = labeller(metric = label_map_metric, component = label_map_component)) +
  labs(
    title = "Functional dissimilarity in ATA space varies across land-use types and warming levels",
    subtitle = "Nine-panel summary across Bray-Curtis, Jaccard, and Sorensen total, turnover, and gradient or nestedness components.",
    x = "Mean temperature anomaly of the site pair",
    y = "Predicted functional dissimilarity (D)"
  ) +
  theme_paper() +
  coord_cartesian(clip = "off") +
  theme(legend.position = "bottom")

save_plot_3formats(
  p_abs,
  "Fig2_Functional_Bray_absolute_predictions",
  file.path(out_root, "03_Plots"),
  width_mm = 260,
  height_mm = 185,
  slide_title = "Functional multi-metric predictions"
)
captions <- record_caption(
  captions, "Fig2", "Fig2_Functional_Bray_absolute_predictions",
  "ATA 空间中的功能结构九宫格主图显示三类距离指标在不同土地利用和变暖背景下并不完全一致",
  "Functional responses in ATA space differ across land-use types and also vary among Bray-Curtis, Jaccard, and Sorensen metrics",
  "把三类功能距离指标和 total、turnover、gradient/nestedness 三个组分并列，可以直接比较不同功能距离定义是否支持同样的生态结论。",
  "Showing the three functional dissimilarity metrics alongside their total, turnover, and gradient or nestedness components makes it much easier to evaluate whether different distance definitions support the same ecological conclusion."
)

pred_diff0_main <- pred_diff0 %>%
  filter(metric %in% main_metrics, component %in% main_components) %>%
  mutate(
    metric = factor(metric, levels = main_metrics),
    component = factor(component, levels = main_components)
  )

p_diff0 <- ggplot(pred_diff0_main, aes(x = pair_clim_mean, y = PredMedian, colour = UI2, fill = UI2)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
  geom_ribbon(aes(ymin = PredLower, ymax = PredUpper), alpha = 0.18, linewidth = 0) +
  geom_line(linewidth = 0.7) +
  scale_colour_manual(values = ui2_cols, labels = label_map_ui2) +
  scale_fill_manual(values = ui2_cols, labels = label_map_ui2) +
  facet_grid(component ~ metric, labeller = labeller(metric = label_map_metric, component = label_map_component)) +
  labs(
    title = "Functional dissimilarity relative to Primary vegetation at climate mean = 0",
    subtitle = "Positive values indicate higher functional dissimilarity than the fixed PV baseline.",
    x = "Mean temperature anomaly of the site pair",
    y = "Difference in functional dissimilarity relative to PV@0"
  ) +
  theme_paper() +
  coord_cartesian(clip = "off") +
  theme(legend.position = "bottom")

save_plot_3formats(
  p_diff0,
  "Fig3_Functional_Bray_diff_vsPV0",
  file.path(out_root, "03_Plots"),
  width_mm = 260,
  height_mm = 185,
  slide_title = "Functional multi-metric contrasts vs PV@0"
)
captions <- record_caption(
  captions, "Fig3", "Fig3_Functional_Bray_diff_vsPV0",
  "固定原生植被基线下，不同土地利用在三类功能距离指标中的偏离方向和幅度并不完全相同",
  "Against a fixed Primary-vegetation baseline, the direction and magnitude of functional shifts are not identical across the three functional dissimilarity metrics",
  "固定基线让不同功能距离指标可以在同一参考系下比较，更容易判断土地利用是否稳定地推动功能同质化。",
  "A fixed baseline places the different functional dissimilarity metrics on the same reference scale, making it easier to judge whether land use consistently promotes functional homogenization."
)

coef_csv <- file.path(out_root, "04_Tables", "functional_bray_main_fixed_effects.csv")
fallback_model_root <- file.path(dirname(project_root), "FunctionalHomogenization_UI2_full_rewrite_DISSIMILARITY", "models_lmer_logitD")

build_func_model_path <- function(metric, component) {
  metric_code <- c("Bray-Curtis" = "bray", "Jaccard" = "jac", "Sorensen" = "sor")[metric]
  component_code <- dplyr::case_when(
    component == "Total" ~ "total",
    component == "Turnover" ~ "turnover",
    metric == "Bray-Curtis" & component == "Gradient/Nestedness" ~ "gradient",
    TRUE ~ "nestedness"
  )
  file.path(fallback_model_root, paste0("mFUNC_", metric_code, "_", component_code, ".rds"))
}

model_specs <- expand.grid(
  metric = main_metrics,
  component = main_components,
  stringsAsFactors = FALSE
) %>%
  dplyr::mutate(
    model = paste(metric, component, sep = " | "),
    path = vapply(seq_len(n()), function(i) build_func_model_path(metric[i], component[i]), character(1))
  )

if (all(file.exists(model_specs$path))) {
  coef_tabs <- lapply(seq_len(nrow(model_specs)), function(i) {
    tb <- coef_table_lmer(readRDS(model_specs$path[i]), keep_pattern = "UI2|pair_clim_mean")
    tb$model <- model_specs$model[i]
    tb$metric <- model_specs$metric[i]
    tb$component <- model_specs$component[i]
    tb
  })
  coef_tab <- bind_rows(coef_tabs)
  write.csv(coef_tab, coef_csv, row.names = FALSE)
} else if (file.exists(coef_csv)) {
  coef_tab <- read.csv(coef_csv, check.names = FALSE)
} else {
  coef_tab <- NULL
  warning("No source_models directory or precomputed fixed-effect table found; skipping coefficient forest plot.")
}

if (!is.null(coef_tab)) {
  coef_forest_df <- coef_tab %>%
    filter(grepl("pair_clim_mean", term)) %>%
    mutate(
      metric = factor(metric, levels = main_metrics),
      component = factor(component, levels = main_components),
      term_pretty = factor(pretty_climate_term(term),
                           levels = rev(c(
                             "Climate main effect",
                             "Secondary vegetation x climate",
                             "Agriculture (low) x climate",
                             "Agriculture (high) x climate",
                             "Urban x climate"
                           )))
    )

  p_coef <- ggplot(coef_forest_df, aes(x = Estimate, y = term_pretty)) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3) +
    geom_errorbar(aes(xmin = Estimate - 1.96 * `Std. Error`,
                      xmax = Estimate + 1.96 * `Std. Error`),
                  orientation = "y", width = 0.12, linewidth = 0.35) +
    geom_point(size = 1.8, colour = "black") +
    facet_grid(component ~ metric, labeller = labeller(metric = label_map_metric, component = label_map_component)) +
    labs(
      title = "Warming-related coefficients across functional dissimilarity metrics",
      subtitle = "Displayed terms are the climate main effect and land-use-specific climate interactions.",
      x = "Fixed-effect estimate on the logit(D) scale",
      y = NULL
    ) +
    theme_paper(base_size = 10) +
    coord_cartesian(clip = "off")

  save_plot_3formats(
    p_coef,
    "Fig4_Functional_Bray_warming_coefficients",
    file.path(out_root, "03_Plots"),
    width_mm = 255,
    height_mm = 185,
    slide_title = "Functional multi-metric warming coefficients"
  )
  captions <- record_caption(
    captions, "Fig4", "Fig4_Functional_Bray_warming_coefficients",
    "三类功能距离指标及其三个组分的气候相关系数并不完全一致，说明 warming 对功能结构的影响取决于所用距离定义",
    "Warming-related coefficients are not identical across Bray-Curtis, Jaccard, and Sorensen functional metrics or among their components",
    "把 nine models 的气候主效应和土地利用特异气候交互放在同一张多 panel 图里，可以直接比较不同功能距离定义对 warming-sensitivity 的一致性。",
    "Displaying the climate main effect and land-use-specific climate interactions from all nine functional models in one facetted figure makes it easier to compare how consistently different distance definitions identify warming sensitivity."
  )
}

sink(file.path(out_root, "04_Tables", "analysis_notes.txt"))
cat("Canonical source script:\n")
cat("scripts/04_community/03_functional_main_rewrite.R and scripts/03_main_models/01_run_lu_climate_models_main.R\n\n")
cat("This rewrite now shows Bray-Curtis, Jaccard, and Sorensen in parallel.\n")
cat("Each metric is displayed for total, turnover, and gradient or nestedness components in combined multi-panel figures.\n")
cat("Public-release note:\n")
cat("- In the lightweight GitHub release, source_models/*.rds may be omitted to reduce repository size.\n")
cat("- When those model objects are absent, the coefficient forest is rebuilt from the precomputed fixed-effect table instead.\n")
sink()

write_caption_registry(
  captions,
  file.path(out_root, "05_Figure_Captions"),
  "functional_main_bray_figure_descriptions_CN_EN"
)

message("Finished: ", out_root)
