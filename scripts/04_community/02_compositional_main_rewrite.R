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

main_components <- c("Total", "Turnover", "Gradient/Nestedness")
pred_abs_main <- pred_abs %>%
  filter(metric == "Bray-Curtis", component %in% main_components) %>%
  mutate(component = factor(component, levels = main_components))

p_abs <- ggplot(pred_abs_main, aes(x = pair_clim_mean, y = PredMedian, colour = UI2, fill = UI2)) +
  geom_ribbon(aes(ymin = PredLower, ymax = PredUpper), alpha = 0.18, linewidth = 0) +
  geom_line(linewidth = 0.7) +
  scale_colour_manual(values = ui2_cols, labels = label_map_ui2) +
  scale_fill_manual(values = ui2_cols, labels = label_map_ui2) +
  facet_wrap(~component, ncol = 3) +
  labs(
    title = "Compositional dissimilarity responds differently to warming across land-use types",
    subtitle = "Main results shown only for Bray-Curtis total, turnover, and gradient components.",
    x = "Pair-level mean temperature anomaly (pair_clim_mean)",
    y = "Predicted compositional dissimilarity (D)"
  ) +
  theme_paper() +
  theme(legend.position = "bottom")

save_plot_3formats(
  p_abs,
  "Fig2_Compositional_Bray_absolute_predictions",
  file.path(out_root, "03_Plots"),
  width_mm = 220,
  height_mm = 120,
  slide_title = "Compositional Bray main predictions"
)
captions <- record_caption(
  captions, "Fig2", "Fig2_Compositional_Bray_absolute_predictions",
  "群落组成 Bray 主结果显示不同土地利用类型对变暖的响应轨迹不同",
  "Bray-based compositional responses to warming differ across land-use types",
  "Bray total、turnover 和 gradient 的主图能够覆盖组成结构变化的整体、替代和梯度三层含义，比把 9 个指标并列展示更集中、更易解释。",
  "Bray total, turnover, and gradient together capture overall dissimilarity, replacement, and gradient components, providing a clearer main narrative than showing all nine metrics in parallel."
)

pred_diff0_main <- pred_diff0 %>%
  filter(metric == "Bray-Curtis", component %in% main_components) %>%
  mutate(component = factor(component, levels = main_components))

p_diff0 <- ggplot(pred_diff0_main, aes(x = pair_clim_mean, y = PredMedian, colour = UI2, fill = UI2)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
  geom_ribbon(aes(ymin = PredLower, ymax = PredUpper), alpha = 0.18, linewidth = 0) +
  geom_line(linewidth = 0.7) +
  scale_colour_manual(values = ui2_cols, labels = label_map_ui2) +
  scale_fill_manual(values = ui2_cols, labels = label_map_ui2) +
  facet_wrap(~component, ncol = 3) +
  labs(
    title = "Compositional homogenization relative to Primary vegetation at climate mean = 0",
    subtitle = "Positive values indicate higher dissimilarity than the fixed PV baseline.",
    x = "Pair-level mean temperature anomaly (pair_clim_mean)",
    y = "Difference in dissimilarity relative to PV@0"
  ) +
  theme_paper() +
  theme(legend.position = "bottom")

save_plot_3formats(
  p_diff0,
  "Fig3_Compositional_Bray_diff_vsPV0",
  file.path(out_root, "03_Plots"),
  width_mm = 220,
  height_mm = 120,
  slide_title = "Compositional Bray contrasts vs PV@0"
)
captions <- record_caption(
  captions, "Fig3", "Fig3_Compositional_Bray_diff_vsPV0",
  "相对于固定的原生植被基线，农业与次生植被在部分 Bray 组成维度上表现出更高的差异化",
  "Relative to the fixed Primary-vegetation baseline, some modified land uses show elevated Bray dissimilarity",
  "固定 PV@0 基线更适合答辩和写作，因为它避免了“同气候条件下的 PV 对照”被误解为真实配对观测值。",
  "The fixed PV@0 baseline is easier to explain because it avoids confusing same-climate PV contrasts with observed pairwise comparisons."
)

coef_csv <- file.path(out_root, "04_Tables", "bray_main_fixed_effects.csv")
model_paths <- c(
  Bray_Total = file.path(project_root, "results", "community", "compositional_main_bray", "source_models", "m_bray_total.rds"),
  Bray_Turnover = file.path(project_root, "results", "community", "compositional_main_bray", "source_models", "m_bray_turnover.rds"),
  Bray_Gradient = file.path(project_root, "results", "community", "compositional_main_bray", "source_models", "m_bray_gradient.rds")
)

if (all(file.exists(model_paths))) {
  mods <- lapply(model_paths, readRDS)
  coef_tabs <- lapply(names(mods), function(nm) {
    tb <- coef_table_lmer(mods[[nm]], keep_pattern = "UI2|pair_clim_mean")
    tb$model <- nm
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
    mutate(term_label = paste(model, term, sep = " | "))

  p_coef <- plot_fixed_effect_forest(
    coef_forest_df,
    term_col = "term_label",
    title = "Warming-related coefficients in the compositional Bray models",
    subtitle = "Displayed terms are the climate main effect and land-use-specific climate interactions.",
    xlab = "Fixed-effect estimate on the logit(D) scale"
  )

  save_plot_3formats(
    p_coef,
    "Fig4_Compositional_Bray_warming_coefficients",
    file.path(out_root, "03_Plots"),
    width_mm = 183,
    height_mm = 120,
    slide_title = "Compositional Bray warming coefficients"
  )
  captions <- record_caption(
    captions, "Fig4", "Fig4_Compositional_Bray_warming_coefficients",
    "组成结构 Bray 主线的气候相关系数显示不同组分对变暖的敏感性并不一致",
    "Warming-related coefficients differ across the compositional Bray components",
    "这一图把主效应和土地利用特异的气候交互放在同一坐标系中，更适合快速判断哪一部分是 warming-sensitive。",
    "Placing the climate main effect and land-use-specific climate interactions on the same axis makes it easier to see which Bray component is most warming-sensitive."
  )
}

sink(file.path(out_root, "04_Tables", "analysis_notes.txt"))
cat("Canonical source script:\n")
cat("scripts/04_community/02_compositional_main_rewrite.R and scripts/03_main_models/01_run_lu_climate_models_main.R\n\n")
cat("This rewrite keeps Bray total / turnover / gradient as the main narrative.\n")
cat("Jaccard and Sorensen remain supplementary and are not re-plotted here.\n\n")
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
