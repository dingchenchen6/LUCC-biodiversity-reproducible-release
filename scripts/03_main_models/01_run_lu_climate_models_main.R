##%######################################################%##
#                                                          #
####   Run the climate/LUI models + diagnostics + plots   ####
#                                                          #
##%######################################################%##
# EN: Full end-to-end script based on our conversation:
#   (1) Read & preprocess PREDICTS site data
#   (2) Fit GLMERSelect models (as in your original script)
#   (3) Save model tables + stats
#   (4) Make manuscript-style figures (Mean anomaly + Max anomaly)
#   (5) Model checking/diagnostics:
#       - Gaussian LMM (abundance): residual vs fitted + QQ (ggplot)
#       - Poisson GLMM (richness): DHARMa diagnostics (base plot) + tests
#   (6) brms (cmdstanr) models at the end (NO poly), plus diagnostics + plots
#
# 中文：基于以上对话的“从头到尾”完整脚本：
#   (1) 读取并预处理 PREDICTS 站点数据
#   (2) 用 GLMERSelect 拟合模型（保留你原始流程）
#   (3) 输出模型表格与统计量
#   (4) 生成论文风格图（平均/最大温度变异）
#   (5) 模型诊断：
#       - 高斯 LMM（丰度）：残差-拟合值 + QQ 图（ggplot）
#       - 计数 GLMM（丰富度）：DHARMa 诊断图（base plot）+ 检验结果
#   (6) 代码最后用 brms（cmdstanr）建模与出图（不使用 poly）
#
# NOTE / 注意：
# - rstandard() 不适用于 lmerMod/merMod，所以 LMM 标准化残差手动计算：resid/sigma
# - DHARMa 的 QQ/Residual 图可能不是 ggplot，所以用“设备捕获”方式导出三种格式
#
# OUTPUT / 输出：
# - 所有主要图：PDF + PNG + PPTX（三种格式）
# - 诊断图：PDF + PNG + PPTX（三种格式）
# - brms 诊断：trace/density/pp_check + LOO + bayes_R2 等输出到文件
#
# ========================================================= #

rm(list = ls())

## =========================
## 0) Directories & logging
## =========================
# EN: adjust to your paths
# 中文：按你的工程目录修改
predictsDataDir <- "5_PREDICTSMatchPropNatHab/"
outDir          <- "6_RunLUClimateModels_new/"

if(!dir.exists(outDir)) dir.create(outDir, recursive = TRUE)

plotDir <- file.path(outDir, "plots_3formats")
diagDir <- file.path(outDir, "diagnostics_3formats")
brmsDir <- file.path(outDir, "brms")
brmsModelDir <- file.path(brmsDir, "models")
brmsPlotDir  <- file.path(brmsDir, "plots_3formats")
brmsDiagDir  <- file.path(brmsDir, "diagnostics")

dir.create(plotDir,      showWarnings = FALSE, recursive = TRUE)
dir.create(diagDir,      showWarnings = FALSE, recursive = TRUE)
dir.create(brmsDir,      showWarnings = FALSE, recursive = TRUE)
dir.create(brmsModelDir, showWarnings = FALSE, recursive = TRUE)
dir.create(brmsPlotDir,  showWarnings = FALSE, recursive = TRUE)
dir.create(brmsDiagDir,  showWarnings = FALSE, recursive = TRUE)

log_file <- file.path(outDir, "log.txt")
sink(log_file, split = TRUE)
t.start <- Sys.time()
cat("\n[START] ", as.character(t.start), "\n")

## =========================
## 1) Libraries
## =========================
# EN: core packages
# 中文：核心包
suppressPackageStartupMessages({
  library(devtools)
  
  library(StatisticalModels)
  library(predictsFunctions)
  source("Functions.R")
  
  library(lme4)
  library(sjPlot)
  library(cowplot)
  library(ggplot2)
  
  # diagnostics
  library(DHARMa)
  library(performance)
  
  # exporting
  library(officer)
  library(rvg)
  
  # brms
  library(brms)
  library(cmdstanr)
  library(bayesplot)
  library(loo)
})

## =========================
## 2) Unified palette / levels
## =========================
# EN: keep the same land-use order + colors (as your original)
# 中文：保持你原始配色与顺序
ui2_levs <- c("Primary vegetation",
              "Secondary vegetation",
              "Agriculture_Low",
              "Agriculture_High",
              "Urban")

ui2_cols <- c("#009E73", "#0072B2", "#E69F00", "#D55E00", "pink")

## =========================
## 3) Helper: export 3 formats (ggplot)
## =========================
# EN: export ggplot -> PDF + PNG + PPTX
# 中文：ggplot 三格式导出
export_plot_3 <- function(p, prefix, out_dir,
                          width_mm = 160, height_mm = 120,
                          dpi_pdf = 300, dpi_png = 600,
                          slide_title = NULL) {
  stopifnot(inherits(p, "ggplot"))
  if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  pdf_path <- file.path(out_dir, paste0(prefix, ".pdf"))
  png_path <- file.path(out_dir, paste0(prefix, ".png"))
  ppt_path <- file.path(out_dir, paste0(prefix, ".pptx"))
  
  ggsave(pdf_path, p, width = width_mm, height = height_mm, units = "mm", dpi = dpi_pdf)
  ggsave(png_path, p, width = width_mm, height = height_mm, units = "mm", dpi = dpi_png)
  
  doc <- read_pptx()
  doc <- add_slide(doc, layout = "Title and Content", master = "Office Theme")
  doc <- ph_with(doc, value = (if(is.null(slide_title)) prefix else slide_title),
                 location = ph_location_type(type = "title"))
  doc <- ph_with(doc, dml(ggobj = p), location = ph_location_fullsize())
  print(doc, target = ppt_path)
  
  invisible(list(pdf = pdf_path, png = png_path, ppt = ppt_path))
}

## =========================
## 4) Helper: export 3 formats (base plot / non-ggplot)  ✅ for DHARMa
## =========================
# EN: DHARMa plots sometimes are NOT ggplot, so we capture devices
# 中文：DHARMa 图可能不是 ggplot，用设备捕获三格式导出
export_baseplot_3 <- function(plot_fun, prefix, out_dir,
                              width_mm = 160, height_mm = 120,
                              dpi_png = 600) {
  stopifnot(is.function(plot_fun))
  if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  pdf_path <- file.path(out_dir, paste0(prefix, ".pdf"))
  png_path <- file.path(out_dir, paste0(prefix, ".png"))
  ppt_path <- file.path(out_dir, paste0(prefix, ".pptx"))
  
  # PDF
  grDevices::pdf(pdf_path, width = width_mm/25.4, height = height_mm/25.4)
  plot_fun()
  grDevices::dev.off()
  
  # PNG
  grDevices::png(png_path, width = width_mm/25.4, height = height_mm/25.4,
                 units = "in", res = dpi_png)
  plot_fun()
  grDevices::dev.off()
  
  # PPTX (place PNG)
  doc <- read_pptx()
  doc <- add_slide(doc, layout = "Title and Content", master = "Office Theme")
  doc <- ph_with(doc, value = prefix, location = ph_location_type(type = "title"))
  doc <- ph_with(doc, external_img(png_path), location = ph_location_fullsize())
  print(doc, target = ppt_path)
  
  invisible(list(pdf = pdf_path, png = png_path, ppt = ppt_path))
}

## =========================
## 5) Helper: LMM diagnostics (Gaussian lmer)  ✅ FIX rstandard() issue
## =========================
# EN: lmerMod has no rstandard(); use resid/sigma
# 中文：lmerMod 没有 rstandard()，改用 resid/sigma 得到标准化残差
lmm_diag_plots <- function(mod) {
  res_raw <- resid(mod)
  sig     <- sigma(mod)
  stdres  <- res_raw / sig
  
  df <- data.frame(
    fitted = fitted(mod),
    resid  = res_raw,
    stdres = stdres
  )
  
  p_rvfit <- ggplot(df, aes(fitted, resid)) +
    geom_point(size = 0.8, alpha = 0.6) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
    theme_bw() +
    xlab("Fitted values / 拟合值") +
    ylab("Residuals / 残差") +
    ggtitle("Gaussian LMM: Residuals vs fitted / 残差-拟合值")
  
  p_qq <- ggplot(df, aes(sample = stdres)) +
    stat_qq(size = 0.8, alpha = 0.6) +
    stat_qq_line(linewidth = 0.4) +
    theme_bw() +
    xlab("Theoretical quantiles / 理论分位数") +
    ylab("Standardised residuals / 标准化残差") +
    ggtitle("Gaussian LMM: Normal Q-Q / 正态QQ图")
  
  list(rvfit = p_rvfit, qq = p_qq)
}

## =========================
## 6) Helper: DHARMa plots + tests (for count GLMM)
## =========================
# EN: use bootstrap in testOutliers for integer distributions
# 中文：对计数分布 outlier 检验建议用 bootstrap（更稳）
dharma_plots <- function(mod, nsim = 500) {
  sim <- DHARMa::simulateResiduals(mod, n = nsim)
  
  out_test  <- DHARMa::testOutliers(sim, type = "bootstrap")
  disp_test <- DHARMa::testDispersion(sim)
  zi_test   <- DHARMa::testZeroInflation(sim)
  
  plot_fun_qq <- function() {
    DHARMa::plotQQunif(sim)
    title("DHARMa QQ (uniform) / QQ图（均匀分布）")
  }
  plot_fun_res <- function() {
    DHARMa::plotResiduals(sim)
    title("DHARMa residual diagnostics / 残差诊断图")
  }
  
  list(sim = sim,
       plot_fun = list(qq = plot_fun_qq, res = plot_fun_res),
       tests = list(out = out_test, disp = disp_test, zi = zi_test))
}

## =========================================================
## Part A) Read & preprocess data
## =========================================================

cat("\n[DATA] reading PREDICTS site data ...\n")
predictsSites <- readRDS(file.path(predictsDataDir, "PREDICTSSitesWithClimateAndNatHab1.rds"))
predictsSites <- predictsSites@data

# UI2 factor + reference level
predictsSites$UI2 <- factor(predictsSites$UI2)
predictsSites$UI2 <- relevel(predictsSites$UI2, ref = "Primary vegetation")
predictsSites$UI2 <- factor(predictsSites$UI2, levels = ui2_levs)

# Standardise climate anomaly predictors
predictsSites$StdTmeanAnomalyRS <- StdCenterPredictor(predictsSites$StdTmeanAnomaly)
predictsSites$StdTmaxAnomalyRS  <- StdCenterPredictor(predictsSites$StdTmaxAnomaly)

# Natural habitat (NH): combine PV + SV
predictsSites$NH_1000  <- predictsSites$PV_1000  + predictsSites$SV_1000
predictsSites$NH_3000  <- predictsSites$PV_3000  + predictsSites$SV_3000
predictsSites$NH_5000  <- predictsSites$PV_5000  + predictsSites$SV_5000
predictsSites$NH_10000 <- predictsSites$PV_10000 + predictsSites$SV_10000

predictsSites$NH_1000.rs  <- StdCenterPredictor(predictsSites$NH_1000)
predictsSites$NH_3000.rs  <- StdCenterPredictor(predictsSites$NH_3000)
predictsSites$NH_5000.rs  <- StdCenterPredictor(predictsSites$NH_5000)
predictsSites$NH_10000.rs <- StdCenterPredictor(predictsSites$NH_10000)

# Abundance transforms in predictsFunctions
predictsSites <- RescaleAbundance(predictsSites)

# Drop unused factor levels
predictsSites <- droplevels(predictsSites)

# Remove NA climate values (threshold filters)
predictsSites <- predictsSites[!is.na(predictsSites$avg_temp), ]

# Quick correlation checks (optional)
cat("\n[CHECK] correlations ...\n")
cat("cor(avg_temp, TmeanAnomaly) = ", cor(predictsSites$avg_temp, predictsSites$TmeanAnomaly, use="complete.obs"), "\n")
cat("cor(avg_temp, StdTmeanAnomaly) = ", cor(predictsSites$avg_temp, predictsSites$StdTmeanAnomaly, use="complete.obs"), "\n")
cat("cor(TmeanAnomaly, StdTmeanAnomaly) = ", cor(predictsSites$TmeanAnomaly, predictsSites$StdTmeanAnomaly, use="complete.obs"), "\n")

# Save processed dataset
saveRDS(predictsSites, file = file.path(outDir, "PREDICTSSiteData.rds"))

cat("\n[INFO] N SSBS levels: ", length(unique(predictsSites$SSBS)), "\n")

## =========================================================
## Part B) GLMERSelect models (same logic as your script)
## =========================================================

## 1) Abundance ~ mean anomaly (Gaussian on LogAbund)
cat("\n[MODEL] GLMERSelect: Abundance ~ mean anomaly ...\n")
model_data <- predictsSites[!is.na(predictsSites$LogAbund), ]
model_data <- model_data[!is.na(model_data$StdTmeanAnomalyRS), ]

MeanAnomalyModelAbund <- GLMERSelect(
  modelData = model_data,
  responseVar = "LogAbund",
  fitFamily = "gaussian",
  fixedFactors = "UI2",
  fixedTerms = list(StdTmeanAnomalyRS = 1),
  randomStruct = "(1|SS)+(1|SSB)",
  fixedInteractions = c("UI2:poly(StdTmeanAnomalyRS,1)"),
  saveVars = c("Species_richness", "Total_abundance", "SSBS", "NH_3000")
)
save(MeanAnomalyModelAbund, file = file.path(outDir, "MeanAnomalyModelAbund.rdata"))

## 2) Richness ~ mean anomaly (Poisson)
cat("\n[MODEL] GLMERSelect: Richness ~ mean anomaly ...\n")
model_data <- predictsSites[!is.na(predictsSites$StdTmeanAnomalyRS), ]

MeanAnomalyModelRich <- GLMERSelect(
  modelData = model_data,
  responseVar = "Species_richness",
  fitFamily = "poisson",
  fixedFactors = "UI2",
  fixedTerms = list(StdTmeanAnomalyRS = 1),
  randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",
  fixedInteractions = c("UI2:poly(StdTmeanAnomalyRS,1)"),
  saveVars = c("SSBS")
)
save(MeanAnomalyModelRich, file = file.path(outDir, "MeanAnomalyModelRich.rdata"))

## 3) Abundance ~ max anomaly (Gaussian)
cat("\n[MODEL] GLMERSelect: Abundance ~ max anomaly ...\n")
model_data <- predictsSites[!is.na(predictsSites$LogAbund), ]
model_data <- model_data[!is.na(model_data$StdTmaxAnomalyRS), ]

MaxAnomalyModelAbund <- GLMERSelect(
  modelData = model_data,
  responseVar = "LogAbund",
  fitFamily = "gaussian",
  fixedFactors = "UI2",
  fixedTerms = list(StdTmaxAnomalyRS = 1),
  randomStruct = "(1|SS)+(1|SSB)",
  fixedInteractions = c("UI2:poly(StdTmaxAnomalyRS,1)"),
  saveVars = c("Species_richness", "Total_abundance", "SSBS", "NH_3000")
)
save(MaxAnomalyModelAbund, file = file.path(outDir, "MaxAnomalyModelAbund.rdata"))

## 4) Richness ~ max anomaly (Poisson)
cat("\n[MODEL] GLMERSelect: Richness ~ max anomaly ...\n")
model_data <- predictsSites[!is.na(predictsSites$StdTmaxAnomalyRS), ]

MaxAnomalyModelRich <- GLMERSelect(
  modelData = model_data,
  responseVar = "Species_richness",
  fitFamily = "poisson",
  fixedFactors = "UI2",
  fixedTerms = list(StdTmaxAnomalyRS = 1),
  randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",
  fixedInteractions = c("UI2:poly(StdTmaxAnomalyRS,1)"),
  saveVars = c("SSBS")
)
save(MaxAnomalyModelRich, file = file.path(outDir, "MaxAnomalyModelRich.rdata"))

## =========================================================
## Part C) Output tables + stats (sjPlot + GLMERSelect stats)
## =========================================================
cat("\n[OUTPUT] model tables + stats ...\n")

tab_model(MeanAnomalyModelAbund$model, transform = NULL,
          file = file.path(outDir, "AbunMeanAnom_output_table.html"))
tab_model(MeanAnomalyModelRich$model, transform = NULL,
          file = file.path(outDir, "RichMeanAnom_output_table.html"))
tab_model(MaxAnomalyModelAbund$model, transform = NULL,
          file = file.path(outDir, "AbunMaxAnom_output_table.html"))
tab_model(MaxAnomalyModelRich$model, transform = NULL,
          file = file.path(outDir, "RichMaxAnom_output_table.html"))

# Save GLMERSelect stats tables
MeanAbunStats <- as.data.frame(MeanAnomalyModelAbund$stats)
MeanRichStats <- as.data.frame(MeanAnomalyModelRich$stats)
MaxAbunStats  <- as.data.frame(MaxAnomalyModelAbund$stats)
MaxRichStats  <- as.data.frame(MaxAnomalyModelRich$stats)

checksig <- function(x) ifelse(x <= 0.05, "Yes", "No")

MeanAbunStats$significant <- sapply(MeanAbunStats$P, checksig)
MeanRichStats$significant <- sapply(MeanRichStats$P, checksig)
MaxAbunStats$significant  <- sapply(MaxAbunStats$P, checksig)
MaxRichStats$significant  <- sapply(MaxRichStats$P, checksig)

write.csv(MeanAbunStats, file = file.path(outDir, "MeanAnomAbun_Stats.csv"), row.names = FALSE)
write.csv(MeanRichStats, file = file.path(outDir, "MeanAnomRich_Stats.csv"), row.names = FALSE)
write.csv(MaxAbunStats,  file = file.path(outDir, "MaxAnomAbun_Stats.csv"),  row.names = FALSE)
write.csv(MaxRichStats,  file = file.path(outDir, "MaxAnomRich_Stats.csv"),  row.names = FALSE)

## =========================================================
## Part D) Manuscript figures (GLMERSelect)  + export 3 formats
## =========================================================
cat("\n[PLOTS] manuscript-style figures (GLMERSelect) ...\n")

exclQuantiles <- c(0.025, 0.975)

## ---------- Figure 2: Mean anomaly (Abundance + Richness) ----------
# EN: build prediction grid (same as your original)
# 中文：构建预测网格（沿用你原始逻辑）

# ---- Abundance (mean) ----
nd_mean_ab <- expand.grid(
  StdTmeanAnomalyRS = seq(min(MeanAnomalyModelAbund$data$StdTmeanAnomalyRS, na.rm=TRUE),
                          max(MeanAnomalyModelAbund$data$StdTmeanAnomalyRS, na.rm=TRUE),
                          length.out = 300),
  UI2 = factor(ui2_levs, levels = ui2_levs)
)

nd_mean_ab$StdTmeanAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd_mean_ab$StdTmeanAnomalyRS,
  originalX    = predictsSites$StdTmeanAnomaly
)

nd_mean_ab$LogAbund <- 0
nd_mean_ab$Species_richness <- 0

refRow <- which(
  nd_mean_ab$UI2 == "Primary vegetation" &
    nd_mean_ab$StdTmeanAnomaly == min(abs(nd_mean_ab$StdTmeanAnomaly))
)[1]

# quantile trimming per UI2
Qlist_mean_ab <- lapply(ui2_levs, function(u) {
  quantile(MeanAnomalyModelAbund$data$StdTmeanAnomalyRS[MeanAnomalyModelAbund$data$UI2 == u],
           probs = exclQuantiles, na.rm = TRUE)
})
names(Qlist_mean_ab) <- ui2_levs

a.preds <- PredictGLMERRandIter(model = MeanAnomalyModelAbund$model, data = nd_mean_ab)
a.preds <- exp(a.preds) - 0.01
a.preds <- sweep(a.preds, 2, a.preds[refRow, ], "/")

# trim
for(u in ui2_levs){
  Q <- Qlist_mean_ab[[u]]
  idx <- nd_mean_ab$UI2 == u & (nd_mean_ab$StdTmeanAnomalyRS < Q[1] | nd_mean_ab$StdTmeanAnomalyRS > Q[2])
  a.preds[idx, ] <- NA
}

nd_mean_ab$PredMedian <- (apply(a.preds, 1, median, na.rm=TRUE) * 100) - 100
nd_mean_ab$PredLower  <- (apply(a.preds, 1, quantile, probs=0.025, na.rm=TRUE) * 100) - 100
nd_mean_ab$PredUpper  <- (apply(a.preds, 1, quantile, probs=0.975, na.rm=TRUE) * 100) - 100

p_mean_ab <- ggplot(nd_mean_ab, aes(x = StdTmeanAnomaly, y = PredMedian)) +
  geom_line(aes(col = UI2), linewidth = 0.75) +
  geom_ribbon(aes(ymin = PredLower, ymax = PredUpper, fill = UI2), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", linewidth = 0.2) +
  scale_fill_manual(values = ui2_cols) +
  scale_colour_manual(values = ui2_cols) +
  theme_bw() +
  scale_x_continuous(breaks = c(0,0.5,1,1.5,2), limits = c(0,2)) +
  scale_y_continuous(breaks = c(-50,-25,0,25,50,75,100,125,150), limits = c(-50,150)) +
  ylab("Change in total abundance (%)") +
  xlab("Standardised Temperature Anomaly") +
  theme(aspect.ratio = 1,
        title = element_text(size = 8, face = "bold"),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 7),
        legend.position = c(0.2, 0.8),
        legend.background = element_blank(),
        legend.text = element_text(size = 6),
        legend.title = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.2),
        panel.border = element_rect(linewidth = 0.2),
        axis.ticks = element_line(linewidth = 0.2)) +
  ggtitle("a")

# ---- Richness (mean) ----
nd_mean_rich <- expand.grid(
  StdTmeanAnomalyRS = seq(min(MeanAnomalyModelRich$data$StdTmeanAnomalyRS, na.rm=TRUE),
                          max(MeanAnomalyModelRich$data$StdTmeanAnomalyRS, na.rm=TRUE),
                          length.out = 300),
  UI2 = factor(ui2_levs, levels = ui2_levs)
)

nd_mean_rich$StdTmeanAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd_mean_rich$StdTmeanAnomalyRS,
  originalX    = predictsSites$StdTmeanAnomaly
)

nd_mean_rich$LogAbund <- 0
nd_mean_rich$Species_richness <- 0

refRow2 <- which(
  nd_mean_rich$UI2 == "Primary vegetation" &
    nd_mean_rich$StdTmeanAnomaly == min(abs(nd_mean_rich$StdTmeanAnomaly))
)[1]

Qlist_mean_rich <- lapply(ui2_levs, function(u) {
  quantile(MeanAnomalyModelRich$data$StdTmeanAnomalyRS[MeanAnomalyModelRich$data$UI2 == u],
           probs = exclQuantiles, na.rm = TRUE)
})
names(Qlist_mean_rich) <- ui2_levs

s.preds <- PredictGLMERRandIter(model = MeanAnomalyModelRich$model, data = nd_mean_rich, nIters = 10000)
s.preds <- exp(s.preds)
s.preds <- sweep(s.preds, 2, s.preds[refRow2, ], "/")

for(u in ui2_levs){
  Q <- Qlist_mean_rich[[u]]
  idx <- nd_mean_rich$UI2 == u & (nd_mean_rich$StdTmeanAnomalyRS < Q[1] | nd_mean_rich$StdTmeanAnomalyRS > Q[2])
  s.preds[idx, ] <- NA
}

nd_mean_rich$PredMedian <- (apply(s.preds, 1, median, na.rm=TRUE) * 100) - 100
nd_mean_rich$PredLower  <- (apply(s.preds, 1, quantile, probs=0.025, na.rm=TRUE) * 100) - 100
nd_mean_rich$PredUpper  <- (apply(s.preds, 1, quantile, probs=0.975, na.rm=TRUE) * 100) - 100

p_mean_rich <- ggplot(nd_mean_rich, aes(x = StdTmeanAnomaly, y = PredMedian)) +
  geom_line(aes(col = UI2), linewidth = 0.75) +
  geom_ribbon(aes(ymin = PredLower, ymax = PredUpper, fill = UI2), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", linewidth = 0.2) +
  scale_fill_manual(values = ui2_cols) +
  scale_colour_manual(values = ui2_cols) +
  theme_bw() +
  scale_x_continuous(breaks = c(0,0.5,1,1.5,2), limits = c(0,2)) +
  scale_y_continuous(breaks = c(-75,-50,-25,0,25,50,75,100), limits = c(-75,100)) +
  ylab("Change in species richness (%)") +
  xlab("Standardised Temperature Anomaly") +
  theme(aspect.ratio = 1,
        title = element_text(size = 8, face = "bold"),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 7),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.2),
        panel.border = element_rect(linewidth = 0.2),
        axis.ticks = element_line(linewidth = 0.2)) +
  ggtitle("b")

p_mean_combo <- cowplot::plot_grid(p_mean_ab, p_mean_rich, ncol = 2)
export_plot_3(p_mean_combo, "Figure2_MeanAnom_Abun_Rich_GLMMselect",
              plotDir, width_mm = 183, height_mm = 100,
              slide_title = "Figure 2 (GLMERSelect): Mean anomaly × UI2")

## ---------- Extended Data: Max anomaly (Abundance + Richness) ----------
# ---- Abundance (max) ----
nd_max_ab <- expand.grid(
  StdTmaxAnomalyRS = seq(min(MaxAnomalyModelAbund$data$StdTmaxAnomalyRS, na.rm=TRUE),
                         max(MaxAnomalyModelAbund$data$StdTmaxAnomalyRS, na.rm=TRUE),
                         length.out = 150),
  UI2 = factor(ui2_levs, levels = ui2_levs)
)

nd_max_ab$StdTmaxAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd_max_ab$StdTmaxAnomalyRS,
  originalX    = predictsSites$StdTmaxAnomaly
)

nd_max_ab$LogAbund <- 0
nd_max_ab$Species_richness <- 0

refRow3 <- which(
  nd_max_ab$UI2 == "Primary vegetation" &
    nd_max_ab$StdTmaxAnomaly == min(abs(nd_max_ab$StdTmaxAnomaly))
)[1]

Qlist_max_ab <- lapply(ui2_levs, function(u) {
  quantile(MaxAnomalyModelAbund$data$StdTmaxAnomalyRS[MaxAnomalyModelAbund$data$UI2 == u],
           probs = exclQuantiles, na.rm = TRUE)
})
names(Qlist_max_ab) <- ui2_levs

a2.preds <- PredictGLMERRandIter(model = MaxAnomalyModelAbund$model, data = nd_max_ab, nIters = 10000)
a2.preds <- exp(a2.preds) - 0.01
a2.preds <- sweep(a2.preds, 2, a2.preds[refRow3, ], "/")

for(u in ui2_levs){
  Q <- Qlist_max_ab[[u]]
  idx <- nd_max_ab$UI2 == u & (nd_max_ab$StdTmaxAnomalyRS < Q[1] | nd_max_ab$StdTmaxAnomalyRS > Q[2])
  a2.preds[idx, ] <- NA
}

nd_max_ab$PredMedian <- (apply(a2.preds, 1, median, na.rm=TRUE) * 100) - 100
nd_max_ab$PredLower  <- (apply(a2.preds, 1, quantile, probs=0.025, na.rm=TRUE) * 100) - 100
nd_max_ab$PredUpper  <- (apply(a2.preds, 1, quantile, probs=0.975, na.rm=TRUE) * 100) - 100

p_max_ab <- ggplot(nd_max_ab, aes(x = StdTmaxAnomaly, y = PredMedian)) +
  geom_line(aes(col = UI2), linewidth = 0.75) +
  geom_ribbon(aes(ymin = PredLower, ymax = PredUpper, fill = UI2), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", linewidth = 0.2) +
  scale_fill_manual(values = ui2_cols) +
  scale_colour_manual(values = ui2_cols) +
  theme_bw() +
  scale_x_continuous(breaks = c(-1,0,1,2,3,4,5), limits = c(-1,5)) +
  scale_y_continuous(breaks = c(-60,-40,-20,0,20,40,60), limits = c(-65,60)) +
  ylab("Change in total abundance (%)") +
  xlab("Maximum Temperature Anomaly") +
  theme(aspect.ratio = 1,
        title = element_text(size = 8, face = "bold"),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 7),
        legend.position = c(0.2, 0.8),
        legend.background = element_blank(),
        legend.text = element_text(size = 6),
        legend.title = element_blank()) +
  ggtitle("a")

# ---- Richness (max) ----
nd_max_rich <- expand.grid(
  StdTmaxAnomalyRS = seq(min(MaxAnomalyModelRich$data$StdTmaxAnomalyRS, na.rm=TRUE),
                         max(MaxAnomalyModelRich$data$StdTmaxAnomalyRS, na.rm=TRUE),
                         length.out = 150),
  UI2 = factor(ui2_levs, levels = ui2_levs)
)

nd_max_rich$StdTmaxAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd_max_rich$StdTmaxAnomalyRS,
  originalX    = predictsSites$StdTmaxAnomaly
)

nd_max_rich$LogAbund <- 0
nd_max_rich$Species_richness <- 0

refRow4 <- which(
  nd_max_rich$UI2 == "Primary vegetation" &
    nd_max_rich$StdTmaxAnomaly == min(abs(nd_max_rich$StdTmaxAnomaly))
)[1]

Qlist_max_rich <- lapply(ui2_levs, function(u) {
  quantile(MaxAnomalyModelRich$data$StdTmaxAnomalyRS[MaxAnomalyModelRich$data$UI2 == u],
           probs = exclQuantiles, na.rm = TRUE)
})
names(Qlist_max_rich) <- ui2_levs

s2.preds <- PredictGLMERRandIter(model = MaxAnomalyModelRich$model, data = nd_max_rich, nIters = 10000)
s2.preds <- exp(s2.preds)
s2.preds <- sweep(s2.preds, 2, s2.preds[refRow4, ], "/")

for(u in ui2_levs){
  Q <- Qlist_max_rich[[u]]
  idx <- nd_max_rich$UI2 == u & (nd_max_rich$StdTmaxAnomalyRS < Q[1] | nd_max_rich$StdTmaxAnomalyRS > Q[2])
  s2.preds[idx, ] <- NA
}

nd_max_rich$PredMedian <- (apply(s2.preds, 1, median, na.rm=TRUE) * 100) - 100
nd_max_rich$PredLower  <- (apply(s2.preds, 1, quantile, probs=0.025, na.rm=TRUE) * 100) - 100
nd_max_rich$PredUpper  <- (apply(s2.preds, 1, quantile, probs=0.975, na.rm=TRUE) * 100) - 100

##%=========================================================%##
##  Max anomaly × Species richness  (GLMERSelect model)
##  Fix: robust refRow + safe trimming
##%=========================================================%##

# ---- 0) prediction grid ----
nd_max_rich <- expand.grid(
  StdTmaxAnomalyRS = seq(
    from = min(MaxAnomalyModelRich$data$StdTmaxAnomalyRS, na.rm = TRUE),
    to   = max(MaxAnomalyModelRich$data$StdTmaxAnomalyRS, na.rm = TRUE),
    length.out = 150
  ),
  UI2 = factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High","Urban"),
               levels = levels(MaxAnomalyModelRich$data$UI2))
)

# back-transform x for plotting (原始尺度用于画图)
nd_max_rich$StdTmaxAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd_max_rich$StdTmaxAnomalyRS,
  originalX    = predictsSites$StdTmaxAnomaly
)

# dummy cols (not used but keep same structure)
nd_max_rich$LogAbund <- 0
nd_max_rich$Species_richness <- 0

# ---- 1) robust reference row: PV within PV rows, closest to 0 anomaly ----
pv_idx <- which(nd_max_rich$UI2 == "Primary vegetation")

# IMPORTANT: use which.min(abs()) instead of "== min(abs())"
refRow <- pv_idx[ which.min(abs(nd_max_rich$StdTmaxAnomaly[pv_idx])) ]

# sanity check
stopifnot(length(refRow) == 1, is.finite(refRow))

# ---- 2) get prediction draws (nrow(nd) × nIters) ----
# 说明 / Note:
# PredictGLMERRandIter returns a matrix: rows=newdata, cols=iterations
s.preds.tmax <- PredictGLMERRandIter(
  model = MaxAnomalyModelRich$model,
  data  = nd_max_rich,
  nIters = 10000
)

# back-transform to richness scale (poisson log-link in GLMERSelect richness models)
s.preds.tmax <- exp(s.preds.tmax)

# ---- 3) convert to relative change vs reference (PV at ~0 anomaly) ----
# Divide EACH iteration column by its reference-row value
den <- s.preds.tmax[refRow, ]  # length = nIters

# if denominator has any NA/0, fix safely
bad_den <- (!is.finite(den)) | den <= 0
if(any(bad_den)){
  warning("Reference-row predictions contain NA/<=0 in some iterations. Replacing with small positive value.")
  den[bad_den] <- min(den[!bad_den], na.rm = TRUE)
  if(!is.finite(den[1])) den[bad_den] <- 1e-6
}

s.preds.tmax_rel <- sweep(s.preds.tmax, MARGIN = 2, STATS = den, FUN = "/")

# ---- 4) SAFE trimming by UI2 quantiles (optional but recommended) ----
exclQuantiles <- c(0.025, 0.975)

for(u in levels(nd_max_rich$UI2)){
  
  x_u <- MaxAnomalyModelRich$data$StdTmaxAnomalyRS[MaxAnomalyModelRich$data$UI2 == u]
  x_u <- x_u[is.finite(x_u)]
  
  # 如果该 land-use 样本太少，跳过裁剪（否则容易全 NA）
  if(length(x_u) < 30){
    message("[trim] Skip ", u, " because n<30 (n=", length(x_u), ")")
    next
  }
  
  Q <- quantile(x_u, probs = exclQuantiles, na.rm = TRUE)
  
  # Q 无效也跳过
  if(any(!is.finite(Q))){
    message("[trim] Skip ", u, " because Q invalid: ", paste(Q, collapse=", "))
    next
  }
  
  idx <- nd_max_rich$UI2 == u &
    (nd_max_rich$StdTmaxAnomalyRS < Q[1] | nd_max_rich$StdTmaxAnomalyRS > Q[2])
  
  # 只把超出范围的行设 NA（不会全杀）
  s.preds.tmax_rel[idx, ] <- NA
}

# ---- 5) summarise to PredMedian/Lower/Upper (% change) ----
nd_max_rich$PredMedian <- (apply(s.preds.tmax_rel, 1, median,   na.rm = TRUE) * 100) - 100
nd_max_rich$PredLower  <- (apply(s.preds.tmax_rel, 1, quantile, probs = 0.025, na.rm = TRUE) * 100) - 100
nd_max_rich$PredUpper  <- (apply(s.preds.tmax_rel, 1, quantile, probs = 0.975, na.rm = TRUE) * 100) - 100

# ---- 6) final sanity checks ----
cat("\n[check] NA counts:\n")
print(colSums(is.na(nd_max_rich[, c("PredMedian","PredLower","PredUpper")])))

cat("\n[check] non-NA points per UI2:\n")
print(with(nd_max_rich, tapply(!is.na(PredMedian), UI2, sum)))



p_max_rich <- ggplot(nd_max_rich, aes(x = StdTmaxAnomaly, y = PredMedian)) +
  geom_line(aes(col = UI2), linewidth = 0.75) +
  geom_ribbon(aes(ymin = PredLower, ymax = PredUpper, fill = UI2), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", linewidth = 0.2) +
  scale_fill_manual(values = ui2_cols) +
  scale_colour_manual(values = ui2_cols) +
  theme_bw() +
  scale_x_continuous(breaks = c(-1,0,1,2,3,4,5), limits = c(-1,5)) +
  scale_y_continuous(breaks = c(-60,-40,-20,0,20,40,60), limits = c(-65,60)) +
  ylab("Change in species richness (%)") +
  xlab("Maximum Temperature Anomaly") +
  theme(aspect.ratio = 1,
        title = element_text(size = 8, face = "bold"),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 7),
        legend.position = "none") +
  ggtitle("b")

p_max_combo <- cowplot::plot_grid(p_max_ab, p_max_rich, ncol = 2)
export_plot_3(p_max_combo, "ExtendedData_MaxAnom_Abun_Rich_GLMMselect",
              plotDir, width_mm = 183, height_mm = 100,
              slide_title = "Extended Data (GLMERSelect): Max anomaly × UI2")

## =========================================================
## Part E) Diagnostics for GLMERSelect models (3 formats)
## =========================================================
cat("\n[DIAG] diagnostics for GLMERSelect models ...\n")

# 1) Gaussian abundance LMM diagnostics (mean + max)
diag_ab_mean <- lmm_diag_plots(MeanAnomalyModelAbund$model)
export_plot_3(diag_ab_mean$rvfit, "Diag_GLMMselect_AbunMean_RVFit", diagDir, 160, 120)
export_plot_3(diag_ab_mean$qq,    "Diag_GLMMselect_AbunMean_QQ",    diagDir, 160, 120)

diag_ab_max <- lmm_diag_plots(MaxAnomalyModelAbund$model)
export_plot_3(diag_ab_max$rvfit,  "Diag_GLMMselect_AbunMax_RVFit",  diagDir, 160, 120)
export_plot_3(diag_ab_max$qq,     "Diag_GLMMselect_AbunMax_QQ",     diagDir, 160, 120)

# 2) Richness Poisson GLMM diagnostics via DHARMa (mean + max)
d_mean <- dharma_plots(MeanAnomalyModelRich$model, nsim = 500)
export_baseplot_3(d_mean$plot_fun$qq,  "Diag_GLMMselect_RichMean_DHARMa_QQ",  diagDir, 160, 120)
export_baseplot_3(d_mean$plot_fun$res, "Diag_GLMMselect_RichMean_DHARMa_Res", diagDir, 160, 120)

writeLines(c(
  "DHARMa tests / DHARMa检验（Mean anomaly Richness）",
  paste0("Outliers (bootstrap) p = ", d_mean$tests$out$p.value),
  paste0("Dispersion p = ", d_mean$tests$disp$p.value),
  paste0("Zero-inflation p = ", d_mean$tests$zi$p.value)
), con = file.path(diagDir, "Diag_GLMMselect_RichMean_DHARMa_tests.txt"))

d_max <- dharma_plots(MaxAnomalyModelRich$model, nsim = 500)
export_baseplot_3(d_max$plot_fun$qq,  "Diag_GLMMselect_RichMax_DHARMa_QQ",  diagDir, 160, 120)
export_baseplot_3(d_max$plot_fun$res, "Diag_GLMMselect_RichMax_DHARMa_Res", diagDir, 160, 120)

writeLines(c(
  "DHARMa tests / DHARMa检验（Max anomaly Richness）",
  paste0("Outliers (bootstrap) p = ", d_max$tests$out$p.value),
  paste0("Dispersion p = ", d_max$tests$disp$p.value),
  paste0("Zero-inflation p = ", d_max$tests$zi$p.value)
), con = file.path(diagDir, "Diag_GLMMselect_RichMax_DHARMa_tests.txt"))

## =========================================================
## Part F) brms (cmdstanr) models + diagnostics + plots (NO poly)
## =========================================================
cat("\n[BRMS] setting cmdstanr backend ...\n")

# EN: ensure CmdStan exists
# 中文：确保 CmdStan 可用
if (is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE))) {
  cat("[BRMS] CmdStan not found; installing...\n")
  cmdstanr::install_cmdstan()
}

# Priors (simple, robust defaults)
prior_gaus <- c(
  prior(normal(0, 2), class = "b"),
  prior(student_t(3, 0, 2.5), class = "sd"),
  prior(student_t(3, 0, 2.5), class = "sigma")
)
prior_nb <- c(
  prior(normal(0, 2), class = "b"),
  prior(student_t(3, 0, 2.5), class = "sd"),
  prior(gamma(2, 0.1), class = "shape")
)

# Data subsets
brms_ab_mean   <- predictsSites[!is.na(predictsSites$LogAbund) & !is.na(predictsSites$StdTmeanAnomalyRS), ]
brms_rich_mean <- predictsSites[!is.na(predictsSites$StdTmeanAnomalyRS), ]
brms_ab_max    <- predictsSites[!is.na(predictsSites$LogAbund) & !is.na(predictsSites$StdTmaxAnomalyRS), ]
brms_rich_max  <- predictsSites[!is.na(predictsSites$StdTmaxAnomalyRS), ]

for(df in list(brms_ab_mean, brms_rich_mean, brms_ab_max, brms_rich_max)){
  # just to force evaluation; no-op
}

brms_ab_mean$UI2   <- factor(brms_ab_mean$UI2,   levels = ui2_levs)
brms_rich_mean$UI2 <- factor(brms_rich_mean$UI2, levels = ui2_levs)
brms_ab_max$UI2    <- factor(brms_ab_max$UI2,    levels = ui2_levs)
brms_rich_max$UI2  <- factor(brms_rich_max$UI2,  levels = ui2_levs)

saveRDS(brms_ab_mean,   file = file.path(brmsDir, "data_ab_mean.rds"))
saveRDS(brms_rich_mean, file = file.path(brmsDir, "data_rich_mean.rds"))
saveRDS(brms_ab_max,    file = file.path(brmsDir, "data_ab_max.rds"))
saveRDS(brms_rich_max,  file = file.path(brmsDir, "data_rich_max.rds"))

## ---------- brms models (NO poly) ----------
cat("\n[BRMS] fitting models (NO poly) ...\n")

# A1) Abundance ~ UI2 * StdTmeanAnomalyRS  (Gaussian on LogAbund)
b_ab_mean <- brm(
  LogAbund ~ UI2 * StdTmeanAnomalyRS + (1|SS) + (1|SSB),
  data    = brms_ab_mean,
  family  = gaussian(),
  prior   = prior_gaus,
  chains  = 4,
  iter    = 6000,
  warmup  = 3000,
  seed    = 123,
  backend = "cmdstanr",
  control = list(adapt_delta = 0.97, max_treedepth = 13),
  save_pars = save_pars(all = TRUE)
)
saveRDS(b_ab_mean, file = file.path(brmsModelDir, "brms_abund_meanAnom_nopoly.rds"))

# A2) Richness ~ UI2 * StdTmeanAnomalyRS  (NegBin)
b_rich_mean <- brm(
  Species_richness ~ UI2 * StdTmeanAnomalyRS + (1|SS) + (1|SSB) + (1|SSBS),
  data    = brms_rich_mean,
  family  = negbinomial(),
  prior   = prior_nb,
  chains  = 4,
  iter    = 8000,
  warmup  = 4000,
  seed    = 123,
  backend = "cmdstanr",
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  save_pars = save_pars(all = TRUE)
)
saveRDS(b_rich_mean, file = file.path(brmsModelDir, "brms_rich_meanAnom_nopoly.rds"))

# B1) Abundance ~ UI2 * StdTmaxAnomalyRS  (Gaussian on LogAbund)
b_ab_max <- brm(
  LogAbund ~ UI2 * StdTmaxAnomalyRS + (1|SS) + (1|SSB),
  data    = brms_ab_max,
  family  = gaussian(),
  prior   = prior_gaus,
  chains  = 4,
  iter    = 6000,
  warmup  = 3000,
  seed    = 123,
  backend = "cmdstanr",
  control = list(adapt_delta = 0.97, max_treedepth = 13),
  save_pars = save_pars(all = TRUE)
)
saveRDS(b_ab_max, file = file.path(brmsModelDir, "brms_abund_maxAnom_nopoly.rds"))

# B2) Richness ~ UI2 * StdTmaxAnomalyRS  (NegBin)
b_rich_max <- brm(
  Species_richness ~ UI2 * StdTmaxAnomalyRS + (1|SS) + (1|SSB) + (1|SSBS),
  data    = brms_rich_max,
  family  = negbinomial(),
  prior   = prior_nb,
  chains  = 4,
  iter    = 8000,
  warmup  = 4000,
  seed    = 123,
  backend = "cmdstanr",
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  save_pars = save_pars(all = TRUE)
)
saveRDS(b_rich_max, file = file.path(brmsModelDir, "brms_rich_maxAnom_nopoly.rds"))

## ---------- brms diagnostics ----------
cat("\n[BRMS] saving diagnostics ...\n")

save_brms_diag <- function(fit, prefix, out_dir, pp_ndraws = 200, do_loo = TRUE) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # posterior summary
  post_sum <- as.data.frame(posterior_summary(fit, probs = c(0.025, 0.975)))
  post_sum$param <- rownames(post_sum)
  rownames(post_sum) <- NULL
  write.csv(post_sum, file.path(out_dir, paste0(prefix, "_posterior_summary.csv")), row.names = FALSE)
  
  # diagnostics text
  np <- nuts_params(fit)
  n_div <- sum(np$Parameter == "divergent__" & np$Value == 1)
  bfmi_vals <- bayesplot::bfmi(np)
  max_td <- max(np$Value[np$Parameter == "treedepth__"], na.rm = TRUE)
  n_tdmax <- sum(np$Parameter == "treedepth__" & np$Value >= max_td, na.rm = TRUE)
  
  writeLines(c(
    paste0("Model: ", prefix),
    "Backend: cmdstanr",
    paste0("Divergent transitions: ", n_div),
    paste0("BFMI (per chain): ", paste(round(bfmi_vals, 3), collapse = ", ")),
    paste0("Max treedepth observed: ", max_td),
    paste0("Max treedepth reached (count): ", n_tdmax),
    "",
    "---- summary(fit) ----",
    capture.output(print(summary(fit)))
  ), con = file.path(out_dir, paste0(prefix, "_diagnostics.txt")))
  
  # trace + density (top parameters)
  pnames <- parnames(fit)
  key_pars <- pnames[grepl("^b_", pnames) | grepl("^sd_", pnames) | grepl("^sigma$", pnames) | grepl("^shape$", pnames)]
  key_pars <- head(key_pars, 12)
  
  if(length(key_pars) > 0){
    png(file.path(out_dir, paste0(prefix, "_trace.png")), 2400, 1800, res = 300)
    print(bayesplot::mcmc_trace(as.array(fit), pars = key_pars))
    dev.off()
    
    png(file.path(out_dir, paste0(prefix, "_dens.png")), 2400, 1800, res = 300)
    print(bayesplot::mcmc_dens(as.array(fit), pars = key_pars))
    dev.off()
  }
  
  # pp_check
  png(file.path(out_dir, paste0(prefix, "_ppcheck.png")), 2400, 1600, res = 300)
  print(pp_check(fit, ndraws = pp_ndraws))
  dev.off()
  
  # bayes_R2
  r2 <- tryCatch(brms::bayes_R2(fit), error = function(e) NULL)
  if(!is.null(r2)){
    write.csv(data.frame(R2 = as.numeric(r2)),
              file.path(out_dir, paste0(prefix, "_bayesR2_draws.csv")), row.names = FALSE)
    r2_sum <- quantile(as.numeric(r2), probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
    writeLines(paste0("bayes_R2 (2.5%, 50%, 97.5%): ", paste(round(r2_sum, 4), collapse = ", ")),
               con = file.path(out_dir, paste0(prefix, "_bayesR2_summary.txt")))
  }
  
  # LOO
  if(do_loo){
    loo_res <- tryCatch(loo::loo(fit), error = function(e) NULL)
    if(!is.null(loo_res)){
      saveRDS(loo_res, file.path(out_dir, paste0(prefix, "_loo.rds")))
      writeLines(capture.output(print(loo_res)),
                 con = file.path(out_dir, paste0(prefix, "_loo.txt")))
      pk <- loo::pareto_k_table(loo_res)
      write.csv(as.data.frame(pk),
                file.path(out_dir, paste0(prefix, "_pareto_k_table.csv")), row.names = FALSE)
    } else {
      writeLines("LOO failed. Check save_pars(all=TRUE) and model.",
                 con = file.path(out_dir, paste0(prefix, "_loo_failed.txt")))
    }
  }
  invisible(TRUE)
}

save_brms_diag(b_ab_mean,   "brms_abund_meanAnom_nopoly", brmsDiagDir)
save_brms_diag(b_rich_mean, "brms_rich_meanAnom_nopoly",  brmsDiagDir)
save_brms_diag(b_ab_max,    "brms_abund_maxAnom_nopoly",  brmsDiagDir)
save_brms_diag(b_rich_max,  "brms_rich_maxAnom_nopoly",   brmsDiagDir)

## ---------- brms predictions + plots (same style as GLMERSelect figures) ----------
cat("\n[BRMS] predictions + manuscript-style plots ...\n")

# EN: helper to compute percent change relative to PV at closest-to-zero anomaly
# 中文：计算相对 PV@最接近0 的百分比变化
brms_pred_percent_change <- function(fit, nd, refRow,
                                     type = c("abund_log", "count"),
                                     probs = c(0.025, 0.975)) {
  type <- match.arg(type)
  
  ep <- posterior_epred(fit, newdata = nd, re_formula = NA)  # draws x rows (response scale)
  
  if(type == "abund_log"){
    # LogAbund -> abundance
    ep <- exp(ep) - 0.01
  }
  
  # relative to reference row (same draw)
  rel <- sweep(ep, 1, ep[, refRow], "/")
  
  nd$PredMedian <- (apply(rel, 2, median, na.rm=TRUE) * 100) - 100
  nd$PredLower  <- (apply(rel, 2, quantile, probs[1], na.rm=TRUE) * 100) - 100
  nd$PredUpper  <- (apply(rel, 2, quantile, probs[2], na.rm=TRUE) * 100) - 100
  nd
}

# EN: trim outside data quantiles within each UI2
# 中文：按各 UI2 的数据分位数裁剪极端范围
trim_by_ui2 <- function(df, ui2_levels, Qlist, xRS) {
  for(u in ui2_levels){
    Q <- Qlist[[u]]
    idx <- df$UI2 == u & (df[[xRS]] < Q[1] | df[[xRS]] > Q[2])
    df$PredMedian[idx] <- NA
    df$PredLower[idx]  <- NA
    df$PredUpper[idx]  <- NA
  }
  df
}

# ========== Mean anomaly (brms) ==========
nd_bmean_ab <- expand.grid(
  StdTmeanAnomalyRS = seq(min(brms_ab_mean$StdTmeanAnomalyRS, na.rm=TRUE),
                          max(brms_ab_mean$StdTmeanAnomalyRS, na.rm=TRUE),
                          length.out = 300),
  UI2 = factor(ui2_levs, levels = ui2_levs)
)
nd_bmean_ab$StdTmeanAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd_bmean_ab$StdTmeanAnomalyRS,
  originalX    = predictsSites$StdTmeanAnomaly
)

ref_bmean <- which(
  nd_bmean_ab$UI2 == "Primary vegetation" &
    nd_bmean_ab$StdTmeanAnomaly == min(abs(nd_bmean_ab$StdTmeanAnomaly))
)[1]

Q_bmean <- lapply(ui2_levs, function(u){
  quantile(brms_ab_mean$StdTmeanAnomalyRS[brms_ab_mean$UI2 == u],
           probs = exclQuantiles, na.rm = TRUE)
}); names(Q_bmean) <- ui2_levs

nd_bmean_ab2 <- brms_pred_percent_change(b_ab_mean, nd_bmean_ab, ref_bmean, type="abund_log")
nd_bmean_ab2 <- trim_by_ui2(nd_bmean_ab2, ui2_levs, Q_bmean, "StdTmeanAnomalyRS")

p_bmean_ab <- ggplot(nd_bmean_ab2, aes(x = StdTmeanAnomaly, y = PredMedian)) +
  geom_line(aes(col = UI2), linewidth = 0.75) +
  geom_ribbon(aes(ymin = PredLower, ymax = PredUpper, fill = UI2), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", linewidth = 0.2) +
  scale_fill_manual(values = ui2_cols) +
  scale_colour_manual(values = ui2_cols) +
  theme_bw() +
  scale_x_continuous(breaks = c(0,0.5,1,1.5,2), limits = c(0,2)) +
  scale_y_continuous(breaks = c(-50,-25,0,25,50,75,100,125,150), limits = c(-50,150)) +
  ylab("Change in total abundance (%)") +
  xlab("Standardised Temperature Anomaly") +
  theme(aspect.ratio = 1,
        title = element_text(size = 8, face = "bold"),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 7),
        legend.position = c(0.2, 0.8),
        legend.background = element_blank(),
        legend.text = element_text(size = 6),
        legend.title = element_blank()) +
  ggtitle("a (brms)")

nd_bmean_rich <- expand.grid(
  StdTmeanAnomalyRS = seq(min(brms_rich_mean$StdTmeanAnomalyRS, na.rm=TRUE),
                          max(brms_rich_mean$StdTmeanAnomalyRS, na.rm=TRUE),
                          length.out = 300),
  UI2 = factor(ui2_levs, levels = ui2_levs)
)
nd_bmean_rich$StdTmeanAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd_bmean_rich$StdTmeanAnomalyRS,
  originalX    = predictsSites$StdTmeanAnomaly
)
ref_bmean2 <- which(
  nd_bmean_rich$UI2 == "Primary vegetation" &
    nd_bmean_rich$StdTmeanAnomaly == min(abs(nd_bmean_rich$StdTmeanAnomaly))
)[1]
Q_bmean2 <- lapply(ui2_levs, function(u){
  quantile(brms_rich_mean$StdTmeanAnomalyRS[brms_rich_mean$UI2 == u],
           probs = exclQuantiles, na.rm = TRUE)
}); names(Q_bmean2) <- ui2_levs

nd_bmean_rich2 <- brms_pred_percent_change(b_rich_mean, nd_bmean_rich, ref_bmean2, type="count")
nd_bmean_rich2 <- trim_by_ui2(nd_bmean_rich2, ui2_levs, Q_bmean2, "StdTmeanAnomalyRS")

p_bmean_rich <- ggplot(nd_bmean_rich2, aes(x = StdTmeanAnomaly, y = PredMedian)) +
  geom_line(aes(col = UI2), linewidth = 0.75) +
  geom_ribbon(aes(ymin = PredLower, ymax = PredUpper, fill = UI2), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", linewidth = 0.2) +
  scale_fill_manual(values = ui2_cols) +
  scale_colour_manual(values = ui2_cols) +
  theme_bw() +
  scale_x_continuous(breaks = c(0,0.5,1,1.5,2), limits = c(0,2)) +
  scale_y_continuous(breaks = c(-75,-50,-25,0,25,50,75,100), limits = c(-75,100)) +
  ylab("Change in species richness (%)") +
  xlab("Standardised Temperature Anomaly") +
  theme(aspect.ratio = 1,
        title = element_text(size = 8, face = "bold"),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 7),
        legend.position = "none") +
  ggtitle("b (brms)")

p_bmean_combo <- cowplot::plot_grid(p_bmean_ab, p_bmean_rich, ncol = 2)
export_plot_3(p_bmean_combo, "Figure2_MeanAnom_Abun_Rich_brms_nopoly",
              brmsPlotDir, width_mm = 183, height_mm = 100,
              slide_title = "Figure 2 (brms, no poly): Mean anomaly × UI2")

# ========== Max anomaly (brms) ==========
nd_bmax_ab <- expand.grid(
  StdTmaxAnomalyRS = seq(min(brms_ab_max$StdTmaxAnomalyRS, na.rm=TRUE),
                         max(brms_ab_max$StdTmaxAnomalyRS, na.rm=TRUE),
                         length.out = 150),
  UI2 = factor(ui2_levs, levels = ui2_levs)
)
nd_bmax_ab$StdTmaxAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd_bmax_ab$StdTmaxAnomalyRS,
  originalX    = predictsSites$StdTmaxAnomaly
)
ref_bmax <- which(
  nd_bmax_ab$UI2 == "Primary vegetation" &
    nd_bmax_ab$StdTmaxAnomaly == min(abs(nd_bmax_ab$StdTmaxAnomaly))
)[1]
Q_bmax <- lapply(ui2_levs, function(u){
  quantile(brms_ab_max$StdTmaxAnomalyRS[brms_ab_max$UI2 == u],
           probs = exclQuantiles, na.rm = TRUE)
}); names(Q_bmax) <- ui2_levs

nd_bmax_ab2 <- brms_pred_percent_change(b_ab_max, nd_bmax_ab, ref_bmax, type="abund_log")
nd_bmax_ab2 <- trim_by_ui2(nd_bmax_ab2, ui2_levs, Q_bmax, "StdTmaxAnomalyRS")

p_bmax_ab <- ggplot(nd_bmax_ab2, aes(x = StdTmaxAnomaly, y = PredMedian)) +
  geom_line(aes(col = UI2), linewidth = 0.75) +
  geom_ribbon(aes(ymin = PredLower, ymax = PredUpper, fill = UI2), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", linewidth = 0.2) +
  scale_fill_manual(values = ui2_cols) +
  scale_colour_manual(values = ui2_cols) +
  theme_bw() +
  scale_x_continuous(breaks = c(-1,0,1,2,3,4,5), limits = c(-1,5)) +
  scale_y_continuous(breaks = c(-60,-40,-20,0,20,40,60), limits = c(-65,60)) +
  ylab("Change in total abundance (%)") +
  xlab("Maximum Temperature Anomaly") +
  theme(aspect.ratio = 1,
        title = element_text(size = 8, face = "bold"),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 7),
        legend.position = c(0.2, 0.8),
        legend.background = element_blank(),
        legend.text = element_text(size = 6),
        legend.title = element_blank()) +
  ggtitle("a (brms)")

nd_bmax_rich <- expand.grid(
  StdTmaxAnomalyRS = seq(min(brms_rich_max$StdTmaxAnomalyRS, na.rm=TRUE),
                         max(brms_rich_max$StdTmaxAnomalyRS, na.rm=TRUE),
                         length.out = 150),
  UI2 = factor(ui2_levs, levels = ui2_levs)
)
nd_bmax_rich$StdTmaxAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd_bmax_rich$StdTmaxAnomalyRS,
  originalX    = predictsSites$StdTmaxAnomaly
)
ref_bmax2 <- which(
  nd_bmax_rich$UI2 == "Primary vegetation" &
    nd_bmax_rich$StdTmaxAnomaly == min(abs(nd_bmax_rich$StdTmaxAnomaly))
)[1]
Q_bmax2 <- lapply(ui2_levs, function(u){
  quantile(brms_rich_max$StdTmaxAnomalyRS[brms_rich_max$UI2 == u],
           probs = exclQuantiles, na.rm = TRUE)
}); names(Q_bmax2) <- ui2_levs

nd_bmax_rich2 <- brms_pred_percent_change(b_rich_max, nd_bmax_rich, ref_bmax2, type="count")
nd_bmax_rich2 <- trim_by_ui2(nd_bmax_rich2, ui2_levs, Q_bmax2, "StdTmaxAnomalyRS")

p_bmax_rich <- ggplot(nd_bmax_rich2, aes(x = StdTmaxAnomaly, y = PredMedian)) +
  geom_line(aes(col = UI2), linewidth = 0.75) +
  geom_ribbon(aes(ymin = PredLower, ymax = PredUpper, fill = UI2), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", linewidth = 0.2) +
  scale_fill_manual(values = ui2_cols) +
  scale_colour_manual(values = ui2_cols) +
  theme_bw() +
  scale_x_continuous(breaks = c(-1,0,1,2,3,4,5), limits = c(-1,5)) +
  scale_y_continuous(breaks = c(-60,-40,-20,0,20,40,60), limits = c(-65,60)) +
  ylab("Change in species richness (%)") +
  xlab("Maximum Temperature Anomaly") +
  theme(aspect.ratio = 1,
        title = element_text(size = 8, face = "bold"),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 7),
        legend.position = "none") +
  ggtitle("b (brms)")

p_bmax_combo <- cowplot::plot_grid(p_bmax_ab, p_bmax_rich, ncol = 2)
export_plot_3(p_bmax_combo, "ExtendedData_MaxAnom_Abun_Rich_brms_nopoly",
              brmsPlotDir, width_mm = 183, height_mm = 100,
              slide_title = "Extended Data (brms, no poly): Max anomaly × UI2")

## =========================================================
## Finish
## =========================================================
t.end <- Sys.time()
cat("\n[FINISHED] ", as.character(t.end), "\n")
cat("[TIME] elapsed: ", round(difftime(t.end, t.start, units="mins"), 2), " minutes\n")

cat("\n[OUTPUT PATHS]\n")
cat("GLMERSelect plots: ", plotDir, "\n")
cat("GLMERSelect diagnostics: ", diagDir, "\n")
cat("brms models: ", brmsModelDir, "\n")
cat("brms plots: ", brmsPlotDir, "\n")
cat("brms diagnostics: ", brmsDiagDir, "\n")

sink()
