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
    title = "Functional Bray dissimilarity in ATA space varies across land-use types and warming levels",
    subtitle = "Main results shown only for the continuous ATA Bray family.",
    x = "Pair-level mean temperature anomaly (pair_clim_mean)",
    y = "Predicted functional dissimilarity (D)"
  ) +
  theme_paper() +
  theme(legend.position = "bottom")

save_plot_3formats(
  p_abs,
  "Fig2_Functional_Bray_absolute_predictions",
  file.path(out_root, "03_Plots"),
  width_mm = 220,
  height_mm = 120,
  slide_title = "Functional Bray main predictions"
)
captions <- record_caption(
  captions, "Fig2", "Fig2_Functional_Bray_absolute_predictions",
  "ATA 空间中的功能 Bray 主图显示农业与城市化并不只是改变组成，也改变功能差异结构",
  "Functional Bray responses in ATA space differ across land-use types",
  "连续 ATA Bray 指标保留了功能性状加权后的梯度信息，比把 ATA 二值化后的 presence/absence 指标更稳定、更可解释。",
  "Continuous ATA Bray metrics retain trait-weighted gradient information and are more stable and interpretable than binarized ATA presence/absence metrics."
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
    title = "Functional dissimilarity relative to Primary vegetation at climate mean = 0",
    subtitle = "Positive values indicate higher functional dissimilarity than the fixed PV baseline.",
    x = "Pair-level mean temperature anomaly (pair_clim_mean)",
    y = "Difference in functional dissimilarity relative to PV@0"
  ) +
  theme_paper() +
  theme(legend.position = "bottom")

save_plot_3formats(
  p_diff0,
  "Fig3_Functional_Bray_diff_vsPV0",
  file.path(out_root, "03_Plots"),
  width_mm = 220,
  height_mm = 120,
  slide_title = "Functional Bray contrasts vs PV@0"
)
captions <- record_caption(
  captions, "Fig3", "Fig3_Functional_Bray_diff_vsPV0",
  "固定原生植被基线下可以更直观看到不同土地利用在功能结构上的偏离方向和强度",
  "A fixed Primary-vegetation baseline clarifies the direction and magnitude of functional shifts",
  "这一表达方式更适合答辩，因为它把“相对原生植被偏离多少”讲得比同气候对照更直接。",
  "This contrast is especially useful in presentations because it makes the deviation from Primary vegetation more directly interpretable."
)

coef_csv <- file.path(out_root, "04_Tables", "functional_bray_main_fixed_effects.csv")
model_paths <- c(
  Functional_Bray_Total = file.path(project_root, "results", "community", "functional_main_bray", "source_models", "mFUNC_bray_total.rds"),
  Functional_Bray_Turnover = file.path(project_root, "results", "community", "functional_main_bray", "source_models", "mFUNC_bray_turnover.rds"),
  Functional_Bray_Gradient = file.path(project_root, "results", "community", "functional_main_bray", "source_models", "mFUNC_bray_gradient.rds")
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
    title = "Warming-related coefficients in the functional Bray models",
    subtitle = "Climate interactions are strongest in the total and gradient components of ATA-space dissimilarity.",
    xlab = "Fixed-effect estimate on the logit(D) scale"
  )

  save_plot_3formats(
    p_coef,
    "Fig4_Functional_Bray_warming_coefficients",
    file.path(out_root, "03_Plots"),
    width_mm = 183,
    height_mm = 120,
    slide_title = "Functional Bray warming coefficients"
  )
  captions <- record_caption(
    captions, "Fig4", "Fig4_Functional_Bray_warming_coefficients",
    "功能 Bray 模型中与变暖相关的项主要集中在 total 和 gradient 组分上",
    "Warming-related coefficients are strongest in the total and gradient components of functional Bray dissimilarity",
    "这说明 warming 对功能结构的影响并不只是简单替代，更涉及功能梯度或差异幅度的改变。",
    "This suggests that warming affects functional structure not only through replacement, but also through changes in gradient-like functional separation."
  )
}

sink(file.path(out_root, "04_Tables", "analysis_notes.txt"))
cat("Canonical source script:\n")
cat("scripts/04_community/03_functional_main_rewrite.R and scripts/03_main_models/01_run_lu_climate_models_main.R\n\n")
cat("This rewrite keeps Bray total / turnover / gradient in ATA space as the main narrative.\n")
cat("Functional Jaccard/Sorensen are intentionally excluded from the main plots because ATA>0 is nearly degenerate in the current data.\n")
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
