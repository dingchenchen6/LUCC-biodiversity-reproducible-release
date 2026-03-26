###############################################################
## Community compositional structure: main Bray rewrite
## 群落组成结构主线重整版
###############################################################

script_dir <- normalizePath(dirname(sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grepl("^--file=", commandArgs(trailingOnly = FALSE))][1])))
source(file.path(script_dir, "00_shared_helpers.R"))
suppressPackageStartupMessages({
  library(lme4)
  library(tidyr)
})

project_root <- get_project_root(script_dir)
base_dir <- file.path(project_root, "results", "community", "compositional_main_bray")
out_root <- file.path(project_root, "results", "community", "compositional_main_bray")

dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_root, "01_DataSupport"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_root, "02_Models"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_root, "03_Plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_root, "04_Tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_root, "05_Figure_Captions"), recursive = TRUE, showWarnings = FALSE)

captions <- list()

cs_data <- readRDS(file.path(project_root, "data", "derived_public", "compositional_pairs_with_covariates.rds"))
pred_abs <- readRDS(file.path(project_root, "results", "community", "compositional_main_bray", "source_data", "lmer_pred_abs_D_all.rds"))
pred_diff0 <- readRDS(file.path(project_root, "results", "community", "compositional_main_bray", "source_data", "lmer_pred_diff_vsPV0_fixed_D_all.rds"))

pred_abs$UI2 <- factor(pred_abs$UI2, levels = ui2_levels)
pred_diff0$UI2 <- factor(pred_diff0$UI2, levels = ui2_levels)

pair_support <- cs_data %>%
  count(UI2, name = "n_pairs") %>%
  mutate(prop = n_pairs / sum(n_pairs))
write.csv(pair_support, file.path(out_root, "04_Tables", "pair_support_by_UI2.csv"), row.names = FALSE)

climate_support <- cs_data %>%
  group_by(UI2) %>%
  summarise(
    n_pairs = n(),
    clim_min = min(pair_clim_mean, na.rm = TRUE),
    clim_q25 = quantile(pair_clim_mean, 0.25, na.rm = TRUE),
    clim_median = median(pair_clim_mean, na.rm = TRUE),
    clim_q75 = quantile(pair_clim_mean, 0.75, na.rm = TRUE),
    clim_max = max(pair_clim_mean, na.rm = TRUE),
    .groups = "drop"
  )
write.csv(climate_support, file.path(out_root, "04_Tables", "pair_climate_support_by_UI2.csv"), row.names = FALSE)

p_support <- build_ui2_support_plot(
  cs_data,
  x_var = "pair_clim_mean",
  x_lab = "Pair-level mean temperature anomaly (pair_clim_mean)",
  title = "Compositional pair support across the climate gradient",
  subtitle = "Urban has the narrowest support; predictions should remain clipped to the empirical range."
)
save_plot_3formats(
  p_support,
  "Fig1_Compositional_pair_climate_support",
  file.path(out_root, "01_DataSupport"),
  width_mm = 183,
  height_mm = 105,
  slide_title = "Compositional pair climate support"
)
captions <- record_caption(
  captions, "Fig1", "Fig1_Compositional_pair_climate_support",
  "群落组成结构配对样本在气候梯度上的支撑分布",
  "Empirical climate support for compositional site pairs",
  "不同土地利用类型的 pair_clim_mean 支撑范围并不一致，Urban 最窄，因此主图预测必须限制在经验支撑区间内。",
  "Climate support differs by land-use type, with Urban showing the narrowest range; this justifies clipping predictions to the empirical support."
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
    title = "Compositional dissimilarity responds differently to warming across land-use types",
    subtitle = "Nine-panel summary across Bray-Curtis, Jaccard, and Sorensen total, turnover, and gradient or nestedness components.",
    x = "Mean temperature anomaly of the site pair",
    y = "Predicted compositional dissimilarity (D)"
  ) +
  theme_paper() +
  coord_cartesian(clip = "off") +
  theme(legend.position = "bottom")

save_plot_3formats(
  p_abs,
  "Fig2_Compositional_Bray_absolute_predictions",
  file.path(out_root, "03_Plots"),
  width_mm = 260,
  height_mm = 185,
  slide_title = "Compositional multi-metric predictions"
)
captions <- record_caption(
  captions, "Fig2", "Fig2_Compositional_Bray_absolute_predictions",
  "群落组成结构九宫格主图显示 Bray-Curtis、Jaccard 和 Sorensen 三类指标在不同土地利用类型下对变暖的响应并不完全一致",
  "Compositional responses to warming differ across land-use types and also vary among Bray-Curtis, Jaccard, and Sorensen metrics",
  "把三类距离指标和 total、turnover、gradient/nestedness 三个组分合并到同一张多 panel 图里，可以直接比较不同指标对同一生态问题给出的响应是否一致。",
  "Combining the three dissimilarity metrics and their total, turnover, and gradient or nestedness components into one multi-panel figure makes it much easier to compare whether different metrics tell a consistent ecological story."
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
    title = "Compositional change relative to Primary vegetation at climate mean = 0",
    subtitle = "Positive values indicate higher dissimilarity than the fixed PV baseline.",
    x = "Mean temperature anomaly of the site pair",
    y = "Difference in dissimilarity relative to PV@0"
  ) +
  theme_paper() +
  coord_cartesian(clip = "off") +
  theme(legend.position = "bottom")

save_plot_3formats(
  p_diff0,
  "Fig3_Compositional_Bray_diff_vsPV0",
  file.path(out_root, "03_Plots"),
  width_mm = 260,
  height_mm = 185,
  slide_title = "Compositional multi-metric contrasts vs PV@0"
)
captions <- record_caption(
  captions, "Fig3", "Fig3_Compositional_Bray_diff_vsPV0",
  "相对于固定的原生植被基线，不同土地利用在三类距离指标下都表现出不同方向和强度的群落结构偏离",
  "Relative to the fixed Primary-vegetation baseline, modified land uses show different magnitudes and directions of structural change across the three dissimilarity metrics",
  "固定 PV@0 基线让三类指标可以在同一框架下比较，也更适合展示土地利用是否一致地推动同质化或异质化。",
  "The fixed PV@0 baseline provides a shared frame of reference across all three metrics, making it easier to compare whether land use consistently promotes homogenization or differentiation."
)

coef_csv <- file.path(out_root, "04_Tables", "bray_main_fixed_effects.csv")
fallback_model_root <- file.path(dirname(project_root), "Homogenization_UI2_full_rewrite_DISSIMILARITY", "models_lmer_logitD")

build_comp_model_path <- function(metric, component) {
  metric_code <- c("Bray-Curtis" = "bray", "Jaccard" = "jac", "Sorensen" = "sor")[metric]
  component_code <- dplyr::case_when(
    component == "Total" ~ "total",
    component == "Turnover" ~ "turnover",
    metric == "Bray-Curtis" & component == "Gradient/Nestedness" ~ "gradient",
    TRUE ~ "nestedness"
  )
  file.path(fallback_model_root, paste0("m_", metric_code, "_", component_code, ".rds"))
}

model_specs <- expand.grid(
  metric = main_metrics,
  component = main_components,
  stringsAsFactors = FALSE
) %>%
  dplyr::mutate(
    model = paste(metric, component, sep = " | "),
    path = vapply(seq_len(n()), function(i) build_comp_model_path(metric[i], component[i]), character(1))
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
      title = "Warming-related coefficients across compositional dissimilarity metrics",
      subtitle = "Displayed terms are the climate main effect and land-use-specific climate interactions.",
      x = "Fixed-effect estimate on the logit(D) scale",
      y = NULL
    ) +
    theme_paper(base_size = 10) +
    coord_cartesian(clip = "off")

  save_plot_3formats(
    p_coef,
    "Fig4_Compositional_Bray_warming_coefficients",
    file.path(out_root, "03_Plots"),
    width_mm = 255,
    height_mm = 185,
    slide_title = "Compositional multi-metric warming coefficients"
  )
  captions <- record_caption(
    captions, "Fig4", "Fig4_Compositional_Bray_warming_coefficients",
    "三类距离指标及其三个组分的气候相关系数并不完全一致，说明不同距离定义会强调群落结构变化的不同方面",
    "Warming-related coefficients are not identical across Bray-Curtis, Jaccard, and Sorensen metrics or among their components",
    "把 nine models 的气候主效应和土地利用特异气候交互放在同一张多 panel 图里，可以直接比较不同距离定义对 warming-sensitivity 的一致性。",
    "Putting the climate main effect and land-use-specific climate interactions from all nine models into one facetted figure makes it much easier to compare how consistently the different dissimilarity metrics identify warming sensitivity."
  )
}

sink(file.path(out_root, "04_Tables", "analysis_notes.txt"))
cat("Canonical source script:\n")
cat("scripts/04_community/02_compositional_main_rewrite.R and scripts/03_main_models/01_run_lu_climate_models_main.R\n\n")
cat("This rewrite now shows Bray-Curtis, Jaccard, and Sorensen in parallel.\n")
cat("Each metric is displayed for total, turnover, and gradient or nestedness components in combined multi-panel figures.\n\n")
cat("Public-release note:\n")
cat("- In the lightweight GitHub release, source_models/*.rds may be omitted to reduce repository size.\n")
cat("- When those model objects are absent, the coefficient forest is rebuilt from the precomputed fixed-effect table instead.\n\n")
cat("Support warning:\n")
cat("- Urban climate support is narrow.\n")
cat("- Predictions should be interpreted only within the empirical pair_clim_mean range of each land-use type.\n")
sink()

write_caption_registry(
  captions,
  file.path(out_root, "05_Figure_Captions"),
  "compositional_main_bray_figure_descriptions_CN_EN"
)

message("Finished: ", out_root)
