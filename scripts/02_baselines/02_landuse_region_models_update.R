
####=========================================================####
####  Land use × Region effects on biodiversity (PREDICTS)  ####
####  土地利用 × 地区 对多样性影响的分析完整脚本              ####
####  全数据 + 剔除小样本 两套模型 & 森林图                 ####
####=========================================================####

rm(list = ls())

##----------------------------------------------------------##
## 0. 目录 & R 包 / Directories & packages
##----------------------------------------------------------##

inDir  <- "1_PreparePREDICTSData/"
outDir <- "2_LUI_RegionModels/"
if (!dir.exists(outDir)) dir.create(outDir, recursive = TRUE)

## StatisticalModels: 这里只用 PredictGLMERRandIter()
library(StatisticalModels)
library(dplyr)
library(ggplot2)
library(forcats)
library(DHARMa)             # 模型诊断 / model diagnostics
library(tibble)
library(brms)               # 贝叶斯 GLMM (可选，可注释掉)
library(lme4)               # glmer / lmer 混合效应模型

options(mc.cores = parallel::detectCores())

##----------------------------------------------------------##
## 1. 读取数据 / Read data
##----------------------------------------------------------##

Sites    <- readRDS(file.path(inDir, "PREDICTSSiteData2.rds"))
Predicts <- readRDS(file.path(inDir, "PREDICTSDatabase1.rds"))

##----------------------------------------------------------##
## 2. 合并地区信息 + 重构 UN_region（美洲拆分为南/北美）
##    Join region info & split Americas into N/S America
##----------------------------------------------------------##

Sites <- Sites %>%
  left_join(
    Predicts %>%
      dplyr::select(SSBS, UN_region, UN_subregion, Realm, Ecoregion) %>%
      distinct(SSBS, .keep_all = TRUE),
    by = "SSBS"
  )

# 使用 UN_subregion 把 Americas 拆成 North America & South America
# Use UN_subregion to split Americas into North America & South America
Sites <- Sites %>%
  mutate(
    UN_region = case_when(
      UN_region == "Americas" & UN_subregion == "North America" ~ "North America",
      UN_region == "Americas" & UN_subregion %in% c("Caribbean",
                                                    "Central America",
                                                    "South America") ~ "South America",
      TRUE ~ as.character(UN_region)
    )
  )

# UI2: 五类土地利用/利用强度 / Land-use intensity categories
ui2_levels <- c("Primary vegetation",
                "Secondary vegetation",
                "Agriculture_Low",
                "Agriculture_High",
                "Urban")

Sites$UI2       <- factor(Sites$UI2, levels = ui2_levels)
Sites$UN_region <- factor(
  Sites$UN_region,
  levels = c("Africa", "Europe", "Asia",
             "North America", "South America", "Oceania")
)

cat("Raw data size:", nrow(Sites), "\n")
cat("Raw UI2 × UN_region table:\n")
print(table(Sites$UI2, Sites$UN_region))

##----------------------------------------------------------##
## 3. 清洗异常数据 / Clean invalid biodiversity data
##----------------------------------------------------------##
##  - 移除 NA / Inf / -Inf / Simpson=0
##  - Remove NA / Inf / -Inf / Simpson=0 (non-biological artefacts)
##----------------------------------------------------------##

Sites_clean <- Sites %>%
  filter(
    !is.na(Species_richness),
    !is.na(LogAbund),
    !is.na(Simpson_diversity)
  ) %>%
  filter(
    is.finite(LogAbund),
    LogAbund >= 0,
    Simpson_diversity > 0,
    is.finite(Simpson_diversity)
  )

cat("Data size after cleaning:", nrow(Sites_clean), "\n")

##----------------------------------------------------------##
## 4. 查看数据分布特征 / Explore data distributions (EDA)
##----------------------------------------------------------##
## 使用清洗后的全数据做探索性分析
## Use Sites_clean for overall EDA
##----------------------------------------------------------##

cat("\n=== Summary (cleaned data): Species_richness / LogAbund / Simpson_diversity ===\n")
print(summary(Sites_clean[, c("Species_richness","LogAbund","Simpson_diversity")]))

# 直方图 / Histograms
p_sr_hist <- ggplot(Sites_clean, aes(Species_richness)) +
  geom_histogram(bins = 40, colour = "black", fill = "grey80") +
  theme_bw() +
  labs(title = "Species richness", x = "Species richness", y = "Count")

p_ab_hist <- ggplot(Sites_clean, aes(LogAbund)) +
  geom_histogram(bins = 40, colour = "black", fill = "grey80") +
  theme_bw() +
  labs(title = "LogAbund (log(Abund+1))", x = "LogAbund", y = "Count")

p_simp_hist <- ggplot(Sites_clean, aes(Simpson_diversity)) +
  geom_histogram(bins = 40, colour = "black", fill = "grey80") +
  theme_bw() +
  labs(title = "Inverse Simpson diversity", x = "Simpson_diversity", y = "Count")

print(p_sr_hist); print(p_ab_hist); print(p_simp_hist)

# 按 UI2 看箱线图 / Boxplots by land-use
p_sr_ui2 <- ggplot(Sites_clean, aes(UI2, Species_richness)) +
  geom_boxplot(outlier.alpha = 0.4) +
  theme_bw() +
  labs(title = "Species richness by land-use", x = "UI2", y = "Species richness") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_ab_ui2 <- ggplot(Sites_clean, aes(UI2, LogAbund)) +
  geom_boxplot(outlier.alpha = 0.4) +
  theme_bw() +
  labs(title = "LogAbund by land-use", x = "UI2", y = "LogAbund") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_simp_ui2 <- ggplot(Sites_clean, aes(UI2, Simpson_diversity)) +
  geom_boxplot(outlier.alpha = 0.4) +
  theme_bw() +
  labs(title = "Simpson diversity by land-use", x = "UI2", y = "Simpson_diversity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p_sr_ui2)
print(p_ab_ui2)
print(p_simp_ui2)

##----------------------------------------------------------##
## 5. UI2 × UN_region 样本量 & 低样本组合 / Sample sizes
##----------------------------------------------------------##
## 这里用 n ≤ 40 识别“小样本组合”，可按需要调整阈值
## Identify low-sample combinations (n ≤ 40) for later filtering
##----------------------------------------------------------##

sr_n <- Sites_clean %>%
  count(UI2, UN_region, name = "n")

cat("\n=== UI2 × UN_region sample sizes (cleaned) ===\n")
print(sr_n)

low_groups <- sr_n %>% filter(n <= 40)

cat("\n=== Low-sample UI2 × UN_region combinations (n ≤ 40) ===\n")
print(low_groups)

##----------------------------------------------------------##
## 6. 构建两套数据：全数据 vs 剔除小样本
##    Build two datasets: full vs filtered (n > 40)
##----------------------------------------------------------##

# 6.1 全数据（仅基础清洗）
# Full data (cleaned only)
Sites_full <- Sites_clean

# 6.2 剔除低样本组合（n ≤ 40）
# Filtered data: remove low-sample UI2 × UN_region combinations
Sites_filt <- Sites_clean %>%
  anti_join(low_groups, by = c("UI2", "UN_region"))

cat("\nData size (full, cleaned only):", nrow(Sites_full), "\n")
cat("Data size (filtered, n > 40 combos):", nrow(Sites_filt), "\n")

cat("\nFull UI2 × UN_region table:\n")
print(table(Sites_full$UI2, Sites_full$UN_region))
cat("\nFiltered UI2 × UN_region table:\n")
print(table(Sites_filt$UI2, Sites_filt$UN_region))

##----------------------------------------------------------##
## 7. 构建模型数据（全数据 & 过滤数据）
##    Prepare model data (full & filtered)
##----------------------------------------------------------##

prep_model_data <- function(Sites_obj) {
  md_sr <- na.omit(Sites_obj[, c("Species_richness", "UI2", "SS", "SSB", "SSBS", "UN_region")])
  md_ab <- na.omit(Sites_obj[, c("LogAbund",        "UI2", "SS", "SSB", "SSBS", "UN_region")])
  md_sd <- na.omit(Sites_obj[, c("Simpson_diversity","UI2", "SS", "SSB", "SSBS", "UN_region")])
  md_sd$LogSimpson <- log(md_sd$Simpson_diversity + 1)
  
  md_sr$UI2 <- factor(md_sr$UI2, levels = ui2_levels)
  md_ab$UI2 <- factor(md_ab$UI2, levels = ui2_levels)
  md_sd$UI2 <- factor(md_sd$UI2, levels = ui2_levels)
  
  md_sr$UN_region <- droplevels(factor(md_sr$UN_region))
  md_ab$UN_region <- droplevels(factor(md_ab$UN_region))
  md_sd$UN_region <- droplevels(factor(md_sd$UN_region))
  
  list(sr = md_sr, ab = md_ab, sd = md_sd)
}

md_full <- prep_model_data(Sites_full)
md_filt <- prep_model_data(Sites_filt)

model_full_sr <- md_full$sr
model_full_ab <- md_full$ab
model_full_sd <- md_full$sd

model_filt_sr <- md_filt$sr
model_filt_ab <- md_filt$ab
model_filt_sd <- md_filt$sd

##----------------------------------------------------------##
## 8. lme4 模型：全数据 & 过滤数据
##    lme4 models: full & filtered
##----------------------------------------------------------##
## 这里显式拟合两类模型：
##  - “Global” 模型：~ UI2        （只看土地利用主效应）
##  - “Region” 模型：~ UI2*UN_region（土地利用 × 地区交互）
##----------------------------------------------------------##

cat("\n=== Fitting lme4 models: FULL data ===\n")

##---------- FULL data: Species richness ----------##
# Global model (UI2 only)
sr_full_UI2 <- glmer(
  Species_richness ~ UI2 + (1|SS) + (1|SSB) + (1|SSBS),
  data   = model_full_sr,
  family = poisson(link = "log")
)

# Region-interaction model (UI2 × UN_region)
sr_full_UI2reg <- glmer(
  Species_richness ~ UI2 * UN_region + (1|SS) + (1|SSB) + (1|SSBS),
  data   = model_full_sr,
  family = poisson(link = "log")
)

##---------- FULL data: Total abundance (LogAbund) ----------##
# Global model (UI2 only) – Gaussian LMM on LogAbund
ab_full_UI2 <- lmer(
  LogAbund ~ UI2 + (1|SS) + (1|SSB),
  data = model_full_ab,
  REML = FALSE
)

# Region-interaction model
ab_full_UI2reg <- lmer(
  LogAbund ~ UI2 * UN_region + (1|SS) + (1|SSB),
  data = model_full_ab,
  REML = FALSE
)

##---------- FULL data: Simpson diversity (LogSimpson) ----------##
sd_full_UI2 <- lmer(
  LogSimpson ~ UI2 + (1|SS) + (1|SSB),
  data = model_full_sd,
  REML = FALSE
)

sd_full_UI2reg <- lmer(
  LogSimpson ~ UI2 * UN_region + (1|SS) + (1|SSB),
  data = model_full_sd,
  REML = FALSE
)

cat("\n=== Fitting lme4 models: FILTERED data (n > 40) ===\n")

##---------- FILTERED data: Species richness ----------##
sr_filt_UI2 <- glmer(
  Species_richness ~ UI2 + (1|SS) + (1|SSB) + (1|SSBS),
  data   = model_filt_sr,
  family = poisson(link = "log")
)

sr_filt_UI2reg <- glmer(
  Species_richness ~ UI2 * UN_region + (1|SS) + (1|SSB) + (1|SSBS),
  data   = model_filt_sr,
  family = poisson(link = "log")
)

##---------- FILTERED data: Total abundance ----------##
ab_filt_UI2 <- lmer(
  LogAbund ~ UI2 + (1|SS) + (1|SSB),
  data = model_filt_ab,
  REML = FALSE
)

ab_filt_UI2reg <- lmer(
  LogAbund ~ UI2 * UN_region + (1|SS) + (1|SSB),
  data = model_filt_ab,
  REML = FALSE
)

##---------- FILTERED data: Simpson diversity ----------##
sd_filt_UI2 <- lmer(
  LogSimpson ~ UI2 + (1|SS) + (1|SSB),
  data = model_filt_sd,
  REML = FALSE
)

sd_filt_UI2reg <- lmer(
  LogSimpson ~ UI2 * UN_region + (1|SS) + (1|SSB),
  data = model_filt_sd,
  REML = FALSE
)

cat("\n=== Model summary (filtered SR: Species_richness ~ UI2 * UN_region) ===\n")
print(summary(sr_filt_UI2reg))

##----------------------------------------------------------##
## 9. 过度离散 + DHARMa 诊断（以过滤数据模型为主）
##    Overdispersion & DHARMa diagnostics (filtered models)
##----------------------------------------------------------##
## 这里我们对三个模型分别做诊断：
##  1) sr_filt_UI2reg : Poisson（物种丰富度）
##  2) ab_filt_UI2reg : Gaussian（总丰度 LogAbund）
##  3) sd_filt_UI2reg : Gaussian（LogSimpson）
##
## 对每个模型：
##  - 使用 DHARMa::simulateResiduals() 生成模拟残差
##  - plot() 画出标准 4-panel 诊断图，并保存为 PNG
##  - 使用 testUniformity(), testDispersion(), testZeroInflation()
##    检查残差的分布、离散度和零膨胀（主要针对 Poisson 模型）
##
## 判读标准（Interpretation guidelines）：
##  - testUniformity: p > 0.05 → 残差分布接近均匀，说明拟合尚可；
##                    p < 0.05 → 模型拟合有系统偏差（需要检查模型结构）。
##  - testDispersion: p > 0.05 → 无显著过度/不足离散；
##                    p < 0.05 → 过度离散（Poisson 常见），
##                                 需考虑加随机效应、改 family、加入解释变量等。
##  - testZeroInflation: p > 0.05 → 零值比例与模型预期一致；
##                       p < 0.05 → 零膨胀或零值不足，考虑零膨胀模型或改 family。
##  - 图形上：
##      * Residuals vs predicted：无明显系统趋势或漏斗形；
##      * QQ 图：点大致沿 1:1 线分布；
##      * Residuals vs predictors：不应随预测变量呈明确函数关系；
##      * 直方图 / ECDF：均匀分布为佳（对 DHARMa 残差）。
##----------------------------------------------------------##

## 9.1 自定义 Poisson 过度离散检验函数（仅对 SR 模型使用）
overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp  <- residuals(model, type = "pearson")
  Pearson.chisq <- sum(rp^2)
  ratio <- Pearson.chisq / rdf
  pval  <- pchisq(Pearson.chisq, df = rdf, lower.tail = FALSE)
  c(chisq = Pearson.chisq, ratio = ratio, rdf = rdf, p = pval)
}

cat("\n=== Overdispersion test (Poisson SR, filtered data) ===\n")
print(overdisp_fun(sr_filt_UI2reg))
## 解读：
##  - ratio ~ 1 且 p > 0.05 → 无明显过度离散（OK）
##  - ratio >> 1 且 p < 0.05 → 存在过度离散（需要考虑模型调整）


## 9.2 DHARMa 诊断：Species richness 模型（Poisson）
set.seed(123)
sim_sr_f <- simulateResiduals(sr_filt_UI2reg, n = 1000)

## 屏幕上查看
plot(sim_sr_f, main = "SR ~ UI2 × Region (filtered, Poisson)")

## 保存为 PNG（4-panel）
png(file.path(outDir, "DHARMa_SR_filtered.png"),
    width = 7, height = 7, units = "in", res = 300)
plot(sim_sr_f, main = "SR ~ UI2 × Region (filtered, Poisson)")
dev.off()

cat("\n=== DHARMa tests for SR (Poisson, filtered data) ===\n")
cat("Uniformity test:\n")
print(testUniformity(sim_sr_f))
cat("\nDispersion test:\n")
print(testDispersion(sim_sr_f))
cat("\nZero-inflation test:\n")
print(testZeroInflation(sim_sr_f))
## 若三个检验的 p 值均 > 0.05，且诊断图无明显模式，一般认为拟合较好；
## 若任一检验 p < 0.05，则需回到模型结构检查是否缺少关键变量、
## 随机效应结构是否合理，或 Poisson family 是否合适等。


## 9.3 DHARMa 诊断：Total abundance 模型（Gaussian, LogAbund）
set.seed(123)
sim_ab_f <- simulateResiduals(ab_filt_UI2reg, n = 1000)

plot(sim_ab_f)

png(file.path(outDir, "DHARMa_Abundance_filtered.png"),
    width = 7, height = 7, units = "in", res = 300)
plot(sim_ab_f)
dev.off()

cat("\n=== DHARMa tests for Abundance (Gaussian, filtered data) ===\n")
cat("Uniformity test:\n")
print(testUniformity(sim_ab_f))
cat("\nDispersion test:\n")
print(testDispersion(sim_ab_f))
## 对 Gaussian 模型，主要关注：
##  - Uniformity：p > 0.05 → 残差近似对称、无明显偏态
##  - Dispersion：p > 0.05 → 方差结构合理
##  - 诊断图中残差 vs 预测值应随机散布，无明显趋势或漏斗形；
##    QQ 图点大致沿直线分布，如偏离严重说明正态性或方差齐性存在问题。


## 9.4 DHARMa 诊断：Simpson 多样性 模型（Gaussian, LogSimpson）
set.seed(123)
sim_sd_f <- simulateResiduals(sd_filt_UI2reg, n = 1000)

plot(sim_sd_f, main = "LogSimpson ~ UI2 × Region (filtered, Gaussian)")

png(file.path(outDir, "DHARMa_Simpson_filtered.png"),
    width = 7, height = 7, units = "in", res = 300)
plot(sim_sd_f, main = "LogSimpson ~ UI2 × Region (filtered, Gaussian)")
dev.off()

cat("\n=== DHARMa tests for Simpson diversity (Gaussian, filtered data) ===\n")
cat("Uniformity test:\n")
print(testUniformity(sim_sd_f))
cat("\nDispersion test:\n")
print(testDispersion(sim_sd_f))
## 判读标准与 Abundance 模型相似：
##  - 若 p 值普遍 > 0.05 且图形上残差无明显结构，即认为模型残差行为合理；
##  - 若 QQ 图偏离严重或 Uniformity/Dispersion p < 0.05，
##    需考虑变换、增加/修改解释变量或随机效应结构。

##----------------------------------------------------------##
## 10. 辅助函数：逆 link & Δ% 相对 Primary vegetation
##     Helper functions for predictions and % change
##----------------------------------------------------------##

# eta -> response
#  - pois:   exp(eta)
#  - log1p:  exp(eta) - 1
#  - id:     identity
inv_link <- function(eta_mat, type = c("pois", "log1p", "id")) {
  type <- match.arg(type)
  switch(
    type,
    "pois"  = exp(eta_mat),
    "log1p" = exp(eta_mat) - 1,
    "id"    = eta_mat
  )
}

# 将不同 UI2 的预测值转为相对 Primary vegetation 的百分比变化
# Convert predictions to % change relative to Primary vegetation
rel_to_primary <- function(pred_mat) {
  rel <- sweep(pred_mat, 2, pred_mat[1, ], "/")
  data.frame(
    median = (apply(rel, 1, median)   * 100) - 100,
    lower  = (apply(rel, 1, quantile, 0.025) * 100) - 100,
    upper  = (apply(rel, 1, quantile, 0.975) * 100) - 100
  )
}

# 对每个地区，计算 UI2 相对该地区 Primary vegetation 的%变化
# For each region, compute % change relative to that region’s Primary vegetation
region_delta_from_mat <- function(pred_mat, newdata, regions, metric_label) {
  out_list <- vector("list", length(regions))
  names(out_list) <- regions
  
  for (r in regions) {
    idx_r  <- which(newdata$UN_region == r)
    idx_pv <- idx_r[newdata$UI2[idx_r] == "Primary vegetation"]
    preds_r <- pred_mat[idx_r, , drop = FALSE]
    rel     <- sweep(preds_r, 2, pred_mat[idx_pv, ], "/")
    med <- (apply(rel, 1, median)   * 100) - 100
    lwr <- (apply(rel, 1, quantile, 0.025) * 100) - 100
    upr <- (apply(rel, 1, quantile, 0.975) * 100) - 100
    
    out_list[[r]] <- data.frame(
      UN_region = r,
      UI2       = newdata$UI2[idx_r],
      median    = med,
      lower     = lwr,
      upper     = upr,
      metric    = metric_label
    )
  }
  
  bind_rows(out_list)
}

##----------------------------------------------------------##
## 11. 计算 global & region Δ%（全数据 / 过滤数据）
##     Compute global & region % changes (full / filtered)
##----------------------------------------------------------##
## 这里的“global”使用 UI2-only 模型；“region”使用 UI2×UN_region 模型
## "global": UI2-only models; "region": UI2 × UN_region models
##----------------------------------------------------------##

build_effects <- function(sr_global, ab_global, sd_global,
                          sr_reg,    ab_reg,    sd_reg,
                          regions) {
  
  ##---------- Global-level (UI2-only models) ----------##
  nd_global <- data.frame(
    UI2               = factor(ui2_levels, levels = ui2_levels),
    Species_richness  = 0,
    LogAbund          = 0,
    LogSimpson        = 0,
    Simpson_diversity = 0
  )
  
  # SR global
  sr_eta_g  <- PredictGLMERRandIter(sr_global, nd_global)
  sr_resp_g <- inv_link(sr_eta_g, "pois")
  sr_delta  <- rel_to_primary(sr_resp_g)
  rich_global <- cbind(
    UI2    = nd_global$UI2,
    sr_delta,
    metric = "Δ species richness (%)"
  )
  
  # Abundance global
  ab_eta_g  <- PredictGLMERRandIter(ab_global, nd_global)
  ab_resp_g <- inv_link(ab_eta_g, "log1p")
  ab_delta  <- rel_to_primary(ab_resp_g)
  abun_global <- cbind(
    UI2    = nd_global$UI2,
    ab_delta,
    metric = "Δ total abundance (%)"
  )
  
  # Simpson global
  sd_eta_g  <- PredictGLMERRandIter(sd_global, nd_global)
  sd_resp_g <- inv_link(sd_eta_g, "log1p")
  sd_delta  <- rel_to_primary(sd_resp_g)
  simp_global <- cbind(
    UI2    = nd_global$UI2,
    sd_delta,
    metric = "Δ Simpson’s diversity (%)"
  )
  
  global_all <- bind_rows(rich_global, abun_global, simp_global)
  global_all[global_all$UI2 == "Primary vegetation", c("lower", "upper")] <- NA
  
  ##---------- Region-level (UI2 × UN_region models) ----------##
  nd_reg <- expand.grid(
    UI2       = factor(ui2_levels, levels = ui2_levels),
    UN_region = regions
  )
  nd_reg$Species_richness  <- 0
  nd_reg$LogAbund          <- 0
  nd_reg$LogSimpson        <- 0
  nd_reg$Simpson_diversity <- 0
  
  # SR region
  sr_eta_reg  <- PredictGLMERRandIter(sr_reg, nd_reg)
  sr_resp_reg <- inv_link(sr_eta_reg, "pois")
  rich_region <- region_delta_from_mat(
    pred_mat     = sr_resp_reg,
    newdata      = nd_reg,
    regions      = regions,
    metric_label = "Δ species richness (%)"
  )
  
  # Abundance region
  ab_eta_reg  <- PredictGLMERRandIter(ab_reg, nd_reg)
  ab_resp_reg <- inv_link(ab_eta_reg, "log1p")
  abun_region <- region_delta_from_mat(
    pred_mat     = ab_resp_reg,
    newdata      = nd_reg,
    regions      = regions,
    metric_label = "Δ total abundance (%)"
  )
  
  # Simpson region
  sd_eta_reg  <- PredictGLMERRandIter(sd_reg, nd_reg)
  sd_resp_reg <- inv_link(sd_eta_reg, "log1p")
  simp_region <- region_delta_from_mat(
    pred_mat     = sd_resp_reg,
    newdata      = nd_reg,
    regions      = regions,
    metric_label = "Δ Simpson’s diversity (%)"
  )
  
  list(
    global = global_all,
    region = bind_rows(rich_region, abun_region, simp_region)
  )
}

regions_full <- levels(model_full_sr$UN_region)
regions_filt <- levels(model_filt_sr$UN_region)

cat("\n=== Building effects (FULL data) ===\n")
eff_full <- build_effects(
  sr_global = sr_full_UI2,
  ab_global = ab_full_UI2,
  sd_global = sd_full_UI2,
  sr_reg    = sr_full_UI2reg,
  ab_reg    = ab_full_UI2reg,
  sd_reg    = sd_full_UI2reg,
  regions   = regions_full
)

cat("\n=== Building effects (FILTERED data) ===\n")
eff_filt <- build_effects(
  sr_global = sr_filt_UI2,
  ab_global = ab_filt_UI2,
  sd_global = sd_filt_UI2,
  sr_reg    = sr_filt_UI2reg,
  ab_reg    = ab_filt_UI2reg,
  sd_reg    = sd_filt_UI2reg,
  regions   = regions_filt
)

global_full <- eff_full$global
region_full <- eff_full$region
global_filt <- eff_filt$global
region_filt <- eff_filt$region

##----------------------------------------------------------##
## 12. 准备森林图数据（非原生土地利用）/ Plot data
##----------------------------------------------------------##
## 目标图形：
##   - 按 UI2 分列分面（4 列）
##   - 按指标分行（3 行：Abundance – Richness – Simpson）
##   - x 轴为 Region，点颜色和形状编码 Region
##   - 灰色横带 + 虚线 = global UI2-only 模型的效应
##   - 灰色柱子（副坐标轴） = 样本量
##----------------------------------------------------------##

nonprim_levels <- c("Secondary vegetation",
                    "Agriculture_Low",
                    "Agriculture_High",
                    "Urban")

ui2_labels <- c("Secondary vegetation",
                "Agriculture (low)",
                "Agriculture (high)",
                "Urban")

# 指标顺序 / metric order in facets (top→bottom)
metric_levels <- c("Δ total abundance (%)",
                   "Δ species richness (%)",
                   "Δ Simpson’s diversity (%)")

prep_plot_df <- function(global_all, region_all) {
  
  # region-level: 只保留非原生土地利用
  region_all2 <- region_all %>%
    filter(UI2 %in% nonprim_levels) %>%
    mutate(
      UI2      = factor(UI2, levels = nonprim_levels),
      UI2_plot = factor(UI2, levels = nonprim_levels, labels = ui2_labels),
      metric   = factor(metric, levels = metric_levels)
    )
  
  # global-level：同样只保留非原生土地利用，并添加 UI2_plot
  global_all2 <- global_all %>%
    filter(UI2 %in% nonprim_levels) %>%
    mutate(
      UI2      = factor(UI2, levels = nonprim_levels),
      UI2_plot = factor(UI2, levels = nonprim_levels, labels = ui2_labels),
      metric   = factor(metric, levels = metric_levels)
    )
  
  list(global = global_all2, region = region_all2)
}

plot_full  <- prep_plot_df(global_full, global_region <- region_full)
plot_filt  <- prep_plot_df(global_filt, global_region <- region_filt)

global_full_plot <- plot_full$global
region_full_plot <- plot_full$region
global_filt_plot <- plot_filt$global
region_filt_plot <- plot_filt$region

##----------------------------------------------------------##
## 13. 森林图函数：灰色柱子表示样本量 & global 背景带
##     Forest plot: grey bars for sample size & global band
##----------------------------------------------------------##

# 地区颜色和点型（可按期刊配色微调）
region_cols <- c(
  "Africa"        = "#1b9e77",
  "Europe"        = "#d95f02",
  "Asia"          = "#7570b3",
  "North America" = "#e7298a",
  "South America" = "#66a61e",
  "Oceania"       = "#e6ab02"
)

region_shapes <- c(
  "Africa"        = 16,
  "Europe"        = 17,
  "Asia"          = 15,
  "North America" = 3,
  "South America" = 4,
  "Oceania"       = 8
)

##----------------------------------------------------------##
## 13. 森林图函数：灰色柱子从底端代表样本量 + 多地区彩色点
##     Forest plot: grey bars = sample size, coloured points = regions
##----------------------------------------------------------##
## 输入：
##  - region_plot: 由 build_effects() + prep_region_plot() 得到的 region_xxx_plot
##  - global_plot: 由 build_effects() + prep_region_plot() 得到的 global_xxx_plot
##  - Sites_obj  : 对应的 Sites_full 或 Sites_filt
##  - tag        : 图标题中使用的标签（如 "FULL (all cleaned sites)"）
##----------------------------------------------------------##

make_forest_plot <- function(region_plot, global_plot, Sites_obj, tag) {
  
  ##------------------------##
  ## 1. 合并样本量信息 n    ##
  ##------------------------##
  # 每个土地利用 × 地区的样本量
  # sample size for each UI2 × UN_region combination
  n_df <- Sites_obj %>%
    dplyr::filter(UI2 %in% nonprim_levels) %>%
    dplyr::count(UN_region, UI2, name = "n")
  
  # region_plot 已经只包含非原生土地利用；这里合并 n
  region_all2 <- region_plot %>%
    dplyr::left_join(n_df, by = c("UN_region", "UI2")) %>%
    dplyr::filter(!is.na(n)) %>%                     # 去掉没有数据的组合
    dplyr::mutate(
      # 统一地区顺序 / make sure region order is consistent
      UN_region = factor(
        UN_region,
        levels = c("Africa", "Asia", "Europe",
                   "North America", "Oceania", "South America")
      ),
      # UI2 和 UI2_plot 顺序
      UI2      = factor(UI2,      levels = nonprim_levels),
      UI2_plot = factor(UI2_plot, levels = ui2_labels),
      # 指标顺序：Abundance（顶行）→ Richness → Simpson
      metric   = factor(metric,   levels = metric_levels),
      # 为 geom_rect 准备数值型 x 位置
      x_num    = as.numeric(UN_region)
    )
  
  ## global 也确保 factor 顺序一致
  global_all2 <- global_plot %>%
    dplyr::filter(UI2 %in% nonprim_levels) %>%
    dplyr::mutate(
      UI2      = factor(UI2,      levels = nonprim_levels),
      UI2_plot = factor(UI2_plot, levels = ui2_labels),
      metric   = factor(metric,   levels = metric_levels)
    )
  
  ##------------------------##
  ## 2. y 轴范围 & n 的缩放  ##
  ##------------------------##
  # 左轴：多样性变化（%）
  y_min <- -100
  y_max <-  100
  effect_breaks <- seq(y_min, y_max, by = 25)  # 可改 20 或 50
  
  # 把样本量 n 线性映射到整个 y 轴高度：
  # n = 0   → y = y_min（面板最底端）
  # n = max → y = y_max（面板最顶端）
  max_n   <- max(region_all2$n, na.rm = TRUE)
  scale_n <- (y_max - y_min) / max_n
  
  # 右侧 sample size 轴刻度（4 个即可，避免太挤）
  n_breaks_raw <- c(0,
                    max_n / 3,
                    2 * max_n / 3,
                    max_n)
  n_breaks <- round(n_breaks_raw, -1)  # 四舍五入到 10 的倍数
  
  ##------------------------##
  ## 3. 开始绘图             ##
  ##------------------------##
  p <- ggplot(region_all2, aes(x = UN_region)) +
    
    ## 灰色样本量柱子：从 y_min 开始往上，长度 ∝ n
    ## Grey bars for sample size: start at bottom (y_min), height ∝ n
    geom_rect(
      data = region_all2,
      aes(xmin = x_num - 0.35,
          xmax = x_num + 0.35,
          ymin = y_min,
          ymax = y_min + n * scale_n),
      inherit.aes = FALSE,
      fill   = "grey90",
      colour = NA
    ) +
    
    ## global 背景带 & global 中位数虚线
    geom_rect(
      data        = global_all2,
      aes(xmin = -Inf, xmax = Inf, ymin = lower, ymax = upper),
      inherit.aes = FALSE,
      fill        = "grey85",
      colour      = NA
    ) +
    geom_hline(
      data        = global_all2,
      aes(yintercept = median),
      inherit.aes = FALSE,
      linetype    = "dotted",
      colour      = "grey40"
    ) +
    
    ## 多样性变化的 0 线
    geom_hline(yintercept = 0, linetype = "dashed", size = 0.3) +
    
    ## 各地区 effect 的置信区间和点
    geom_errorbar(
      aes(ymin = lower, ymax = upper, colour = UN_region),
      width = 0, size = 0.35
    ) +
    geom_point(
      aes(y = median, colour = UN_region, shape = UN_region),
      size = 2
    ) +
    
    ## 行：指标；列：土地利用 / use intensity
    facet_grid(metric ~ UI2_plot, scales = "fixed") +
    
    ## 地区颜色 & 点型（在前面定义的 region_cols / region_shapes）
    scale_colour_manual(
      name   = "Region / 地区",
      values = region_cols,
      drop   = FALSE
    ) +
    scale_shape_manual(
      name   = "Region / 地区",
      values = region_shapes,
      drop   = FALSE
    ) +
    
    ## 左轴：相对原生植被的变化（%）
    ## 右轴：样本量 sample size (n)，使用反变换恢复 n
    scale_y_continuous(
      limits   = c(y_min, y_max),
      breaks   = effect_breaks,
      name     = "Change relative to primary vegetation (%)\n相对原生植被的变化 (%)",
      sec.axis = sec_axis(
        trans  = ~ (. - y_min) / scale_n,   # y → n
        breaks = n_breaks,
        name   = "Sample size (n)"
      )
    ) +
    
    labs(
      x     = NULL,
      title = paste0("Land-use × Region effects (", tag, ")")
    ) +
    
    theme_bw() +
    theme(
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.text.x        = element_text(angle = 45, hjust = 1, size = 7),
      axis.text.y        = element_text(size = 7),
      strip.text         = element_text(size = 9),
      strip.background   = element_rect(fill = "white"),
      legend.position    = "top",
      legend.title       = element_text(size = 9, face = "bold"),
      legend.text        = element_text(size = 8),
      aspect.ratio       = 0.8,              # 每个面板“瘦高一点”
      panel.spacing.y    = unit(1.2, "lines")
    )
  
  p
}

# 全数据森林图
p_forest_full <- make_forest_plot(
  region_plot = region_full_plot,
  global_plot = global_full_plot,
  Sites_obj   = Sites_full,
  tag         = "FULL (all cleaned sites)"
)

# 过滤数据森林图
p_forest_filt <- make_forest_plot(
  region_plot = region_filt_plot,
  global_plot = global_filt_plot,
  Sites_obj   = Sites_filt,
  tag         = "FILTERED (n > 40 per UI2 × Region)"
)


print(p_forest_full)

ggsave(file.path(outDir, "LUI_region_forest_FULL_likeFigureB.png"),
       p_forest_full, width = 9, height = 8, units = "in", dpi = 600)

print(p_forest_filt)

ggsave(file.path(outDir, "LUI_region_forest_FILTERED_likeFigureB.png"),
       p_forest_filt, width = 9, height = 8, units = "in", dpi = 600)

##----------------------------------------------------------##
## 14. 诊断：Global vs 低样本组合（3 个多样性指数）
##     Diagnostics: global vs low-sample combinations
##----------------------------------------------------------##
## 思路：
##  - Global：所有站点整体分布
##  - Low-sample：所有 n ≤ 阈值 的 UI2 × UN_region 组合中的站点
##  - 逐个指数画箱线图，对比 Global vs 低样本组合
##----------------------------------------------------------##

####-------------------- 14.1 Species richness --------------------####

low_data_diag_sr <- Sites_clean %>%
  semi_join(low_groups, by = c("UI2", "UN_region")) %>%
  dplyr::select(UI2, UN_region, Species_richness) %>%
  mutate(group = paste(UI2, UN_region, sep = " | "))

global_diag_sr <- Sites_clean %>%
  mutate(group = "Global") %>%
  dplyr::select(group, Species_richness)

compare_sr <- bind_rows(global_diag_sr, low_data_diag_sr)

n_df_sr <- compare_sr %>%
  group_by(group) %>%
  summarise(
    n     = n(),
    y_pos = max(Species_richness, na.rm = TRUE) * 1.05,
    .groups = "drop"
  )

p_sr_diag <- ggplot(compare_sr,
                    aes(x = group, y = Species_richness, fill = group)) +
  geom_boxplot(outlier.alpha = 0.5) +
  geom_text(
    data = n_df_sr,
    aes(x = group, y = y_pos, label = paste0("n = ", n)),
    inherit.aes = FALSE,
    vjust = 0, size = 3
  ) +
  theme_bw() +
  labs(x = "Group / 组别",
       y = "Species richness",
       title = "Species richness: Global vs low-sample UI2 × Region groups") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

print(p_sr_diag)
ggsave(file.path(outDir, "SR_lowGroups_vs_global.png"),
       p_sr_diag, width = 8, height = 5, units = "in", dpi = 600)


####-------------------- 14.2 Total abundance (LogAbund) --------------------####

low_data_diag_ab <- Sites_clean %>%
  semi_join(low_groups, by = c("UI2", "UN_region")) %>%
  dplyr::select(UI2, UN_region, LogAbund) %>%
  mutate(group = paste(UI2, UN_region, sep = " | "))

global_diag_ab <- Sites_clean %>%
  mutate(group = "Global") %>%
  dplyr::select(group, LogAbund)

compare_ab <- bind_rows(global_diag_ab, low_data_diag_ab)

n_df_ab <- compare_ab %>%
  group_by(group) %>%
  summarise(
    n     = n(),
    y_pos = max(LogAbund, na.rm = TRUE) * 1.05,
    .groups = "drop"
  )

p_ab_diag <- ggplot(compare_ab,
                    aes(x = group, y = LogAbund, fill = group)) +
  geom_boxplot(outlier.alpha = 0.5) +
  geom_text(
    data = n_df_ab,
    aes(x = group, y = y_pos, label = paste0("n = ", n)),
    inherit.aes = FALSE,
    vjust = 0, size = 3
  ) +
  theme_bw() +
  labs(x = "Group / 组别",
       y = "Total abundance (log scale)\n总丰度（对数尺度）",
       title = "Total abundance: Global vs low-sample UI2 × Region groups") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

print(p_ab_diag)
ggsave(file.path(outDir, "Abundance_lowGroups_vs_global.png"),
       p_ab_diag, width = 8, height = 5, units = "in", dpi = 600)


####-------------------- 14.3 Simpson diversity --------------------####

low_data_diag_sd <- Sites_clean %>%
  semi_join(low_groups, by = c("UI2", "UN_region")) %>%
  dplyr::select(UI2, UN_region, Simpson_diversity) %>%
  mutate(group = paste(UI2, UN_region, sep = " | "))

global_diag_sd <- Sites_clean %>%
  mutate(group = "Global") %>%
  dplyr::select(group, Simpson_diversity)

compare_simpson <- bind_rows(global_diag_sd, low_data_diag_sd)

n_df_simpson <- compare_simpson %>%
  group_by(group) %>%
  summarise(
    n     = n(),
    y_pos = max(Simpson_diversity, na.rm = TRUE) * 1.05,
    .groups = "drop"
  )

p_simpson_diag <- ggplot(compare_simpson,
                         aes(x = group, y = Simpson_diversity, fill = group)) +
  geom_boxplot(outlier.alpha = 0.5) +
  geom_text(
    data = n_df_simpson,
    aes(x = group, y = y_pos, label = paste0("n = ", n)),
    inherit.aes = FALSE,
    vjust = 0, size = 3
  ) +
  theme_bw() +
  labs(x = "Group / 组别",
       y = "Simpson diversity",
       title = "Simpson diversity: Global vs low-sample UI2 × Region groups") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

print(p_simpson_diag)
ggsave(file.path(outDir, "Simpson_lowGroups_vs_global.png"),
       p_simpson_diag, width = 8, height = 5, units = "in", dpi = 600)

##==========================================================##
## 15. brms 贝叶斯模型（全数据 + 过滤数据）
##     brms Bayesian models (full & filtered data)
##==========================================================##
## 说明 / Notes:
##  - 本节使用 brms 对三个多样性指标拟合层级贝叶斯模型；
##  - 同时拟合：
##       * 全数据 full（仅基础清洗）
##       * 过滤数据 filtered（剔除 UI2 × Region 样本量 ≤ 阈值的组合）
##  - 物种丰富度：同时拟合 Poisson 与 Negative Binomial，以比较拟合质量；
##  - 丰度与 Simpson：此处先给出 Gaussian(log) 版本，与前面 GLMER 对应；
##    更复杂的分布（Gamma / NB）已在后续步骤（16–18 节）展开。
##----------------------------------------------------------##

##---------------- 15.0 并行与后端设置 / Parallel & backend ----------------##

## 使用所有可用核心进行链间并行（强烈推荐）
## Use all available cores for parallel chains
## 并行设置：使用总核心数 - 1
n_cores <- max(parallel::detectCores() - 1, 1)
options(mc.cores = n_cores)
cat("Using", n_cores, "cores for parallel sampling.\n")

## 可选：使用 cmdstanr 后端加速（通常显著快于 rstan）
## Optional: use cmdstanr backend for faster sampling
## 需要先安装 cmdstanr 并编译 CmdStan，详见官方文档
## Uncomment if cmdstanr is installed and configured:

install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
library(cmdstanr)
cmdstanr::install_cmdstan()
#cmdstanr::set_cmdstan_path()  # 只需运行一次
backend_brms <- "cmdstanr"
## 若不使用 cmdstanr，则保持默认后端：
backend_brms <- "rstan"

##---------------- 15.1 统一先验 / Priors ----------------##

## 物种丰富度（Poisson 或 NB）的先验 / Priors for SR models
b_prior_sr_pois <- c(
  prior(normal(0, 2), class = "b"),          # 固定效应 / fixed effects
  prior(student_t(3, 0, 2.5), class = "sd")  # 随机效应标准差 / random-effect SD
)

## Negative Binomial（多一个 shape 参数）
b_prior_sr_nb <- c(
  prior(normal(0, 2), class = "b"),
  prior(student_t(3, 0, 2.5), class = "sd"),
  prior(gamma(2, 0.1), class = "shape")      # 负二项的形状参数 / NB shape
)

## 丰度 & Simpson（Gaussian on log）的先验
b_prior_gaus <- c(
  prior(normal(0, 2), class = "b"),
  prior(student_t(3, 0, 2.5), class = "sd"),
  prior(student_t(3, 0, 2.5), class = "sigma")
)

## 控制参数 / Sampler control
brm_ctrl <- list(adapt_delta = 0.95, max_treedepth = 12)

##==========================================================##
## 15.2 物种丰富度模型 / Species richness models
##==========================================================##
## 模型结构 / Model structure:
##   Species_richness ~ UI2 * UN_region + (1|SS) + (1|SSB) + (1|SSBS)
## 分别对：
##   - 全数据 model_full_sr
##   - 过滤数据 model_filt_sr
## 拟合 Poisson 和 Negative Binomial 模型。
## 之后可通过 LOO/AIC、pp_check 等比较哪个分布更合适。
##==========================================================##

##----- 15.2.1 全数据：Poisson SR 模型 / FULL data: Poisson SR -----##

b_sr_full_pois <- brm(
  Species_richness ~ UI2 * UN_region + (1|SS) + (1|SSB) + (1|SSBS),
  data    = model_full_sr,
  family  = poisson(),
  prior   = b_prior_sr_pois,
  chains  = 4, iter = 4000, warmup = 2000,
  control = brm_ctrl,
  seed    = 123,
  backend = backend_brms
)

##----- 15.2.2 过滤数据：Poisson SR 模型 / FILTERED: Poisson SR -----##

b_sr_filt_pois <- brm(
  Species_richness ~ UI2 * UN_region + (1|SS) + (1|SSB) + (1|SSBS),
  data    = model_filt_sr,
  family  = poisson(),
  prior   = b_prior_sr_pois,
  chains  = 4, iter = 4000, warmup = 2000,
  control = brm_ctrl,
  seed    = 124,
  backend = backend_brms
)

##----- 15.2.3 过滤数据：Negative Binomial SR 模型 / FILTERED: NB SR -----##

b_sr_filt_nb <- brm(
  Species_richness ~ UI2 * UN_region + (1|SS) + (1|SSB) + (1|SSBS),
  data    = model_filt_sr,
  family  = negbinomial(),       # 关键：负二项分布 / Negative binomial
  prior   = b_prior_sr_nb,
  chains  = 4, iter = 4000, warmup = 2000,
  control = brm_ctrl,
  seed    = 125,
  backend = backend_brms
)
b_prior_sr_nb <- c(
  prior(normal(0, 2), class = "b"),                 # fixed effects
  prior(student_t(3, 0, 2.5), class = "sd"),        # random-effect SDs
  prior(normal(0, 1), class = "shape", lb = 0)      # shape > 0，收紧一些
)

b_sr_filt_nb <- brm(
  Species_richness ~ UI2 * UN_region + (1|SS) + (1|SSB) + (1|SSBS),
  data    = model_filt_sr,
  family  = negbinomial(link = "log"),
  prior   = b_prior_sr_nb,
  chains  = 4,
  iter    = 6000,
  warmup  = 3000,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  backend = "cmdstanr",
  threads = threading(4),
  seed    = 123
)

##----- 15.2.4 SR 模型摘要与诊断 / Summaries & diagnostics -----##

cat("\n=== brms SR models (full & filtered) ===\n")
print(summary(b_sr_full_pois), digits = 2)
print(summary(b_sr_filt_pois), digits = 2)
print(summary(b_sr_filt_nb),   digits = 2)

## Rhat 接近 1、有效样本量大，说明收敛良好；
## 若有 Rhat > 1.01 或 n_eff 很小，则需增加迭代次数或调整 adapt_delta。

## 轨迹图 / Trace plots（重点检查 NB 过滤模型）
plot(b_sr_filt_nb)

## 后验预测检验 / Posterior predictive checks
pp_check(b_sr_filt_pois, nsamples = 200)
pp_check(b_sr_filt_nb,   nsamples = 200)

## 判读：
##  - pp_check 曲线/柱状图若与观测分布相近，说明模型拟合较好；
##  - 若 Poisson 与 NB 对比中，NB 在尾部或方差结构上匹配更好，
##    则可在文中说明“NB better captures count variability than Poisson”.

##（可选）使用 loo 比较 Poisson vs NB（过滤数据）
## Optional: compare Poisson vs NB via LOO
# loo_sr_pois <- loo(b_sr_filt_pois)
# loo_sr_nb   <- loo(b_sr_filt_nb)
# loo_compare(loo_sr_pois, loo_sr_nb)


##==========================================================##
## 15.3 总丰度模型 / Total abundance models (LogAbund)
##==========================================================##
## 模型结构：
##   LogAbund ~ UI2 * UN_region + (1|SS) + (1|SSB)
## 这里先用 Gaussian family 拟合 log(Abund+1)，
## 与前面的 GLMER 分析保持一致，便于对比。
## 更适合的计数模型 / Gamma 模型已在后续步骤中给出。
##==========================================================##

##----- 15.3.1 全数据：LogAbund ~ Gaussian -----##

b_ab_full <- brm(
  LogAbund ~ UI2 * UN_region + (1|SS) + (1|SSB),
  data    = model_full_ab,
  family  = gaussian(),
  prior   = b_prior_gaus,
  chains  = 4, iter = 4000, warmup = 2000,
  control = brm_ctrl,
  seed    = 126,
  backend = backend_brms
)

##----- 15.3.2 过滤数据：LogAbund ~ Gaussian -----##

b_ab_filt <- brm(
  LogAbund ~ UI2 * UN_region + (1|SS) + (1|SSB),
  data    = model_filt_ab,
  family  = gaussian(),
  prior   = b_prior_gaus,
  chains  = 4, iter = 4000, warmup = 2000,
  control = brm_ctrl,
  seed    = 127,
  backend = backend_brms
)

cat("\n=== brms Abundance models (full & filtered) ===\n")
print(summary(b_ab_full), digits = 2)
print(summary(b_ab_filt), digits = 2)

## 轨迹图与 pp_check（可重点看过滤版）
plot(b_ab_filt)
pp_check(b_ab_filt, nsamples = 200)

## 判读：
##  - 轨迹图应呈“毛毛虫”状、链间无明显系统差异；
##  - pp_check 若显示预测分布在尾部与观测差异大，
##    则可在后续步骤中采用更合适的分布（如 NB 或 Gamma）。


##==========================================================##
## 15.4 Simpson 多样性模型 / Simpson diversity models (LogSimpson)
##==========================================================##
## 模型结构：
##   LogSimpson ~ UI2 * UN_region + (1|SS) + (1|SSB)
## 与丰度类似，先用 Gaussian 拟合 log-transformed 指数，
## 与前面的 GLMER 保持一致。更自然的 Gamma 模型
## 已在第 18 节给出。
##==========================================================##

##----- 15.4.1 全数据：LogSimpson ~ Gaussian -----##

b_sd_full <- brm(
  LogSimpson ~ UI2 * UN_region + (1|SS) + (1|SSB),
  data    = model_full_sd,
  family  = gaussian(),
  prior   = b_prior_gaus,
  chains  = 4, iter = 4000, warmup = 2000,
  control = brm_ctrl,
  seed    = 128,
  backend = backend_brms
)

##----- 15.4.2 过滤数据：LogSimpson ~ Gaussian -----##

b_sd_filt <- brm(
  LogSimpson ~ UI2 * UN_region + (1|SS) + (1|SSB),
  data    = model_filt_sd,
  family  = gaussian(),
  prior   = b_prior_gaus,
  chains  = 4, iter = 4000, warmup = 2000,
  control = brm_ctrl,
  seed    = 129,
  backend = backend_brms
)

cat("\n=== brms Simpson models (full & filtered) ===\n")
print(summary(b_sd_full), digits = 2)
print(summary(b_sd_filt), digits = 2)

## 轨迹与预测检验
plot(b_sd_filt)
pp_check(b_sd_filt, nsamples = 200)

## 判读：
##  - 若 Gaussian(LogSimpson) 在 pp_check 中表现不佳，
##    可在后续步骤采用 Gamma(Simpson_diversity) 等更合适的分布；
##  - 这里的 Gaussian 模型主要用于与传统 GLMM 结果对比。

##==========================================================##
## 15.5 边际效应图（例：过滤数据模型）
##     Marginal effects (UI2 × UN_region, filtered)
##==========================================================##

ce_sr_nb <- conditional_effects(b_sr_filt_nb, effects = "UI2:UN_region")
ce_ab_f  <- conditional_effects(b_ab_filt,   effects = "UI2:UN_region")
ce_sd_f  <- conditional_effects(b_sd_filt,   effects = "UI2:UN_region")

plot(ce_sr_nb)
plot(ce_ab_f)
plot(ce_sd_f)

## 说明：
##  - 这些边际效应图展示了在过滤数据上 UI2 × Region 对三个指标的影响；
##  - 在正式文章中可以选取 NB SR + 改进后的丰度/Simpson 模型的边际效应图，
##    作为主结果展示（当前 Gaussian 版可作为 Supplementary）。


##----------------------------------------------------------##
## 16. 改进的物种丰富度模型：负二项与零膨胀（glmmTMB）
##     Improved SR models: Negative Binomial & Zero-inflated (glmmTMB)
##----------------------------------------------------------##
## 背景：
##  - 第 9 节 Poisson 模型（sr_filt_UI2reg）的 DHARMa 诊断显示：
##      * 非均匀残差（uniformity p << 0.05）
##      * 显著离散性问题（underdispersion）
##      * 零值比例不匹配（zero-inflation p << 0.05）
##  - 说明 Poisson family（方差 = 均值）假设不合适，需要更灵活的计数分布。
##
## 这里使用 glmmTMB 拟合两类模型：
##  1) 负二项（nbinom2）：允许方差 > 均值，是生态计数数据最常用的替代方案；
##  2) 零膨胀负二项（nbinom2 + ziformula）：若零值问题依然存在，可进一步优化。
##
## 并用 DHARMa 对两个模型进行诊断，并与原 Poisson 模型在 AIC 上进行比较。
##----------------------------------------------------------##

## 如未加载 glmmTMB，先加载

library(glmmTMB)

##-----------------------------##
## 16.1 负二项模型（nbinom2）  ##
##-----------------------------##

cat("\n=== Fitting Negative Binomial model for Species richness (filtered data) ===\n")

sr_nb <- glmmTMB::glmmTMB(
  Species_richness ~ UI2 * UN_region +
    (1 | SS) + (1 | SSB) + (1 | SSBS),
  data   = model_filt_sr,
  family = glmmTMB::nbinom2()   # 方差 = mu + mu^2/shape
)

summary(sr_nb)

## AIC 对比：Poisson vs NB
cat("\n=== AIC comparison: Poisson vs Negative Binomial (SR, filtered data) ===\n")
aic_comp_sr <- AIC(sr_filt_UI2reg, sr_nb)
print(aic_comp_sr)
## 解读：
##  - AIC 更小的模型通常更优；
##  - 若 nb 模型 AIC 明显小于 Poisson，说明负二项更适合物种丰富度数据。 


## DHARMa 诊断：负二项 SR 模型
set.seed(123)
sim_sr_nb <- DHARMa::simulateResiduals(sr_nb, n = 1000)

## 屏幕查看
plot(sim_sr_nb, main = "DHARMa NB: SR ~ UI2 × Region (filtered, nbinom2)")

## 保存为 PNG
png(file.path(outDir, "DHARMa_SR_NB_filtered.png"),
    width = 7, height = 7, units = "in", res = 300)
plot(sim_sr_nb, main = "DHARMa NB: SR ~ UI2 × Region (filtered, nbinom2)")
dev.off()

cat("\n=== DHARMa tests for SR (Negative Binomial, filtered data) ===\n")
cat("Uniformity test:\n")
print(DHARMa::testUniformity(sim_sr_nb))
cat("\nDispersion test:\n")
print(DHARMa::testDispersion(sim_sr_nb))
cat("\nZero-inflation test:\n")
print(DHARMa::testZeroInflation(sim_sr_nb))

## 判读：
##  - 若 Uniformity / Dispersion / Zero-inflation 的 p 值大多 > 0.05，
##    且诊断图中残差无明显模式，则 NB 模型在拟合和残差行为上优于 Poisson；
##  - 若仍有明显问题，可进一步尝试零膨胀模型或 COM-Poisson 等更灵活分布。


##--------------------------------------##
## 16.2 零膨胀负二项模型（可选，若需要） ##
##--------------------------------------##
## 若从 16.1 中看到 Zero-inflation p < 0.05（零值问题仍存在），
## 则可以拟合零膨胀负二项模型进一步改进。
## 若不需要，可将本小节整段注释掉。

cat("\n=== Fitting Zero-inflated Negative Binomial model for Species richness (filtered data) ===\n")

sr_nb_zi <- glmmTMB::glmmTMB(
  Species_richness ~ UI2 * UN_region +
    (1 | SS) + (1 | SSB) + (1 | SSBS),
  ziformula = ~ 1,             # 最简单的零膨胀结构：全局常数
  data      = model_filt_sr,
  family    = glmmTMB::nbinom2()
)

summary(sr_nb_zi)

## AIC 对比：Poisson vs NB vs ZI-NB
cat("\n=== AIC comparison: Poisson vs NB vs ZI-NB (SR, filtered data) ===\n")
aic_comp_sr_zi <- AIC(sr_filt_UI2reg, sr_nb, sr_nb_zi)
print(aic_comp_sr_zi)
## 一般情况：
##  - 若 ZI-NB AIC 显著更小 → 零膨胀结构有助于改善拟合；
##  - 若 NB 和 ZI-NB AIC 接近，可优先选择更简单的 NB 模型。


## DHARMa 诊断：ZI-NB SR 模型
set.seed(123)
sim_sr_nb_zi <- DHARMa::simulateResiduals(sr_nb_zi, n = 1000)

plot(sim_sr_nb_zi, main = "DHARMa ZI-NB: SR ~ UI2 × Region (filtered)")

png(file.path(outDir, "DHARMa_SR_ZINB_filtered.png"),
    width = 7, height = 7, units = "in", res = 300)
plot(sim_sr_nb_zi, main = "DHARMa ZI-NB: SR ~ UI2 × Region (filtered)")
dev.off()

cat("\n=== DHARMa tests for SR (Zero-inflated NB, filtered data) ===\n")
cat("Uniformity test:\n")
print(DHARMa::testUniformity(sim_sr_nb_zi))
cat("\nDispersion test:\n")
print(DHARMa::testDispersion(sim_sr_nb_zi))
cat("\nZero-inflation test:\n")
print(DHARMa::testZeroInflation(sim_sr_nb_zi))

## 判读建议：
##  - 若 ZI-NB 模型在 Uniformity / Dispersion / Zero-inflation 检验中
##    显示 p 值大多 > 0.05，且 AIC 最低，则可将其视为物种丰富度的“主推模型”；
##  - 若 NB 模型已表现良好（诊断和 AIC 都不错），可优先使用 NB 模型，
##    并在文中说明：对零膨胀结构的拓展并未显著改善模型拟合。

##----------------------------------------------------------##
## 17. 改进 Total abundance 模型：负二项计数模型（glmmTMB）
##     Improved Abundance model: Negative Binomial (glmmTMB)
##----------------------------------------------------------##
## 背景：
##  - 第 9 节中，Gaussian 模型（ab_filt_UI2reg: LogAbund）在 DHARMa 诊断中
##    显示残差分布严重偏离均匀、显著 underdispersion。
##  - 这表明 Gaussian family 不适合丰度数据。
##
## 思路：
##  - 利用 LogAbund = log(Abund + 1) 的关系，反推原始计数：
##        Abund_raw = exp(LogAbund) - 1
##  - 使用 glmmTMB 对 Abund_raw 拟合负二项 GLMM（nbinom2），
##    并用 DHARMa 检验残差行为，比较 AIC。
##----------------------------------------------------------##

## 17.1 构造丰度计数变量 Abund_raw
model_filt_ab2 <- model_filt_ab %>%
  dplyr::mutate(
    Abund_raw = exp(LogAbund) - 1   # 精确反变换：Abund = exp(log(Abund+1)) - 1
  )

## 可选：若你不希望出现非整数（理论上不会），可加 round()
## model_filt_ab2$Abund_raw <- round(model_filt_ab2$Abund_raw)

## 17.2 拟合负二项 GLMM（Abund_raw ~ UI2 * UN_region）
cat("\n=== Fitting Negative Binomial model for Total abundance (filtered data) ===\n")

ab_nb <- glmmTMB::glmmTMB(
  Abund_raw ~ UI2 * UN_region + (1 | SS) + (1 | SSB),
  data   = model_filt_ab2,
  family = glmmTMB::nbinom2()
)

summary(ab_nb)

## 17.3 与原 Gaussian 模型比较 AIC（ab_filt_UI2reg）
cat("\n=== AIC comparison: Gaussian(LogAbund) vs NB(Abund_raw) ===\n")
aic_comp_ab <- AIC(ab_filt_UI2reg, ab_nb)
print(aic_comp_ab)
## 解读：
##  - 若 ab_nb 的 AIC 显著低于 ab_filt_UI2reg，说明负二项计数模型更适合丰度数据；
##  - 若两者差异不大，可在文中说明两种建模方式结果一致，但计数模型残差表现更好。

## 17.4 DHARMa 诊断：Abund_raw 负二项模型
set.seed(123)
sim_ab_nb <- DHARMa::simulateResiduals(ab_nb, n = 1000)

## 屏幕查看诊断图（4-panel）
plot(sim_ab_nb, main = "DHARMa NB: Abund_raw ~ UI2 × Region (filtered, nbinom2)")
plot(sim_ab_nb)
## 保存为 PNG
png(file.path(outDir, "DHARMa_Abundance_NB_filtered.png"),
    width = 7, height = 7, units = "in", res = 300)
plot(sim_ab_nb, main = "DHARMa NB: Abund_raw ~ UI2 × Region (filtered, nbinom2)")
dev.off()

cat("\n=== DHARMa tests for Total abundance (Negative Binomial, filtered data) ===\n")
cat("Uniformity test:\n")
print(DHARMa::testUniformity(sim_ab_nb))
cat("\nDispersion test:\n")
print(DHARMa::testDispersion(sim_ab_nb))
cat("\nZero-inflation test:\n")
print(DHARMa::testZeroInflation(sim_ab_nb))

## 判读标准：
##  - Uniformity: p > 0.05 → 残差近似均匀，模型整体拟合合理；
##  - Dispersion: p > 0.05 → 无明显过度/不足离散；
##  - Zero-inflation: p > 0.05 → 零值比例与模型预期一致；
##  - 图形上 Residuals vs Predicted 随机散布，无明显趋势/漏斗形，
##    QQ 图点基本沿 1:1 线分布，即可认为负二项模型在丰度上表现良好。

####===================  End of script  ===================####










####=========================================================####
####  Land use × Region effects on biodiversity (PREDICTS)  ####
####  土地利用 × 地区 对多样性影响的分析完整脚本              ####
####  全数据 + 剔除小样本 两套模型 & 诊断 + 森林图 + brms    ####
####  图像三种格式：PNG / PDF / PPTX                         ####
####=========================================================####

rm(list = ls())

##----------------------------------------------------------##
## 0. 目录 & R 包 / Directories & packages
##----------------------------------------------------------##

inDir  <- "1_PreparePREDICTSData/"
outDir <- "2_LUI_RegionModels_update/"

if (!dir.exists(outDir)) dir.create(outDir, recursive = TRUE)

## 子文件夹 / Sub-folders
dir_models_glmer   <- file.path(outDir, "models_glmer")
dir_models_brms    <- file.path(outDir, "models_brms")
dir_models_glmmTMB <- file.path(outDir, "models_glmmTMB")

dir_plots_eda      <- file.path(outDir, "plots_eda")
dir_plots_dharma   <- file.path(outDir, "plots_dharma")
dir_plots_forest   <- file.path(outDir, "plots_forest")
dir_plots_lowN     <- file.path(outDir, "plots_lowN")
dir_plots_brms     <- file.path(outDir, "plots_brms")

dir_data_effects   <- file.path(outDir, "data_effects")
dir_data_diag      <- file.path(outDir, "data_diag")

dirs_to_make <- c(
  dir_models_glmer, dir_models_brms, dir_models_glmmTMB,
  dir_plots_eda, dir_plots_dharma, dir_plots_forest,
  dir_plots_lowN, dir_plots_brms,
  dir_data_effects, dir_data_diag
)

for (d in dirs_to_make) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

## R packages
library(StatisticalModels)  # PredictGLMERRandIter()
library(dplyr)
library(ggplot2)
library(forcats)
library(DHARMa)
library(tibble)
library(brms)
library(lme4)
library(glmmTMB)
library(export)            # graph2ppt

options(mc.cores = parallel::detectCores()-1)

##----------------------------------------------------------##
## 1. 读取数据 / Read data
##----------------------------------------------------------##

Sites    <- readRDS(file.path(inDir, "PREDICTSSiteData2.rds"))
Predicts <- readRDS(file.path(inDir, "PREDICTSDatabase1.rds"))

##----------------------------------------------------------##
## 2. 合并地区信息 & 重构 UN_region（美洲拆南/北） 
##    Join region info & split Americas into N/S America
##----------------------------------------------------------##

Sites <- Sites %>%
  dplyr::left_join(
    Predicts %>%
      dplyr::select(SSBS, UN_region, UN_subregion, Realm, Ecoregion) %>%
      dplyr::distinct(SSBS, .keep_all = TRUE),
    by = "SSBS"
  )

Sites <- Sites %>%
  dplyr::mutate(
    UN_region = dplyr::case_when(
      UN_region == "Americas" & UN_subregion == "North America" ~ "North America",
      UN_region == "Americas" & UN_subregion %in%
        c("Caribbean","Central America","South America")        ~ "South America",
      TRUE ~ as.character(UN_region)
    )
  )

ui2_levels <- c("Primary vegetation",
                "Secondary vegetation",
                "Agriculture_Low",
                "Agriculture_High",
                "Urban")

Sites$UI2 <- factor(Sites$UI2, levels = ui2_levels)
Sites$UN_region <- factor(
  Sites$UN_region,
  levels = c("Africa", "Europe", "Asia",
             "North America", "South America", "Oceania")
)

cat("Raw data size:", nrow(Sites), "\n")
cat("Raw UI2 × UN_region table:\n")
print(table(Sites$UI2, Sites$UN_region))

##----------------------------------------------------------##
## 3. 清洗异常数据 / Clean invalid biodiversity data
##----------------------------------------------------------##

Sites_clean <- Sites %>%
  dplyr::filter(
    !is.na(Species_richness),
    !is.na(LogAbund),
    !is.na(Simpson_diversity)
  ) %>%
  dplyr::filter(
    is.finite(LogAbund),
    LogAbund >= 0,                      # 排除逻辑上不合理值
    Simpson_diversity > 0,
    is.finite(Simpson_diversity)
  )

cat("Data size after cleaning:", nrow(Sites_clean), "\n")

##----------------------------------------------------------##
## 4. 数据分布特征 / Data distribution (EDA)
##----------------------------------------------------------##

cat("\n=== Summary (cleaned data): SR / LogAbund / Simpson_diversity ===\n")
print(summary(Sites_clean[, c("Species_richness","LogAbund","Simpson_diversity")]))

## 4.1 全局直方图 / Global histograms
p_sr_hist <- ggplot(Sites_clean, aes(Species_richness)) +
  geom_histogram(bins = 40, colour = "black", fill = "grey80") +
  theme_bw() +
  labs(title = "Species richness", x = "Species richness", y = "Count")

p_ab_hist <- ggplot(Sites_clean, aes(LogAbund)) +
  geom_histogram(bins = 40, colour = "black", fill = "grey80") +
  theme_bw() +
  labs(title = "LogAbund (log(Abund+1))", x = "LogAbund", y = "Count")

p_simp_hist <- ggplot(Sites_clean, aes(Simpson_diversity)) +
  geom_histogram(bins = 40, colour = "black", fill = "grey80") +
  theme_bw() +
  labs(title = "Inverse Simpson diversity", x = "Simpson_diversity", y = "Count")

print(p_sr_hist); print(p_ab_hist); print(p_simp_hist)

## 保存三种格式 / Save in PNG + PDF + PPTX
ggsave(file.path(dir_plots_eda, "Hist_SR.png"),   p_sr_hist, width = 6, height = 4, dpi = 400)
ggsave(file.path(dir_plots_eda, "Hist_SR.pdf"),   p_sr_hist, width = 6, height = 4)
export::graph2ppt(file = file.path(dir_plots_eda, "Hist_SR.pptx"),
                  x = p_sr_hist, width = 6, height = 4)

ggsave(file.path(dir_plots_eda, "Hist_LogAbund.png"),  p_ab_hist, width = 6, height = 4, dpi = 400)
ggsave(file.path(dir_plots_eda, "Hist_LogAbund.pdf"),  p_ab_hist, width = 6, height = 4)
export::graph2ppt(file = file.path(dir_plots_eda, "Hist_LogAbund.pptx"),
                  x = p_ab_hist, width = 6, height = 4)

ggsave(file.path(dir_plots_eda, "Hist_Simpson.png"),   p_simp_hist, width = 6, height = 4, dpi = 400)
ggsave(file.path(dir_plots_eda, "Hist_Simpson.pdf"),   p_simp_hist, width = 6, height = 4)
export::graph2ppt(file = file.path(dir_plots_eda, "Hist_Simpson.pptx"),
                  x = p_simp_hist, width = 6, height = 4)

## 4.2 按土地利用类型 / By land-use (看不同类型的响应范围)
p_sr_ui2 <- ggplot(Sites_clean, aes(UI2, Species_richness)) +
  geom_boxplot(outlier.alpha = 0.4) +
  theme_bw() +
  labs(title = "Species richness by land-use",
       x = "UI2 (land-use type)", y = "Species richness") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_ab_ui2 <- ggplot(Sites_clean, aes(UI2, LogAbund)) +
  geom_boxplot(outlier.alpha = 0.4) +
  theme_bw() +
  labs(title = "LogAbund by land-use",
       x = "UI2 (land-use type)", y = "LogAbund") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_simp_ui2 <- ggplot(Sites_clean, aes(UI2, Simpson_diversity)) +
  geom_boxplot(outlier.alpha = 0.4) +
  theme_bw() +
  labs(title = "Simpson diversity by land-use",
       x = "UI2 (land-use type)", y = "Simpson_diversity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p_sr_ui2); print(p_ab_ui2); print(p_simp_ui2)

ggsave(file.path(dir_plots_eda, "Box_SR_UI2.png"),   p_sr_ui2, width = 7, height = 4.5, dpi = 400)
ggsave(file.path(dir_plots_eda, "Box_SR_UI2.pdf"),   p_sr_ui2, width = 7, height = 4.5)
export::graph2ppt(file = file.path(dir_plots_eda, "Box_SR_UI2.pptx"),
                  x = p_sr_ui2, width = 7, height = 4.5)

ggsave(file.path(dir_plots_eda, "Box_LogAbund_UI2.png"),  p_ab_ui2, width = 7, height = 4.5, dpi = 400)
ggsave(file.path(dir_plots_eda, "Box_LogAbund_UI2.pdf"),  p_ab_ui2, width = 7, height = 4.5)
export::graph2ppt(file = file.path(dir_plots_eda, "Box_LogAbund_UI2.pptx"),
                  x = p_ab_ui2, width = 7, height = 4.5)

ggsave(file.path(dir_plots_eda, "Box_Simpson_UI2.png"),   p_simp_ui2, width = 7, height = 4.5, dpi = 400)
ggsave(file.path(dir_plots_eda, "Box_Simpson_UI2.pdf"),   p_simp_ui2, width = 7, height = 4.5)
export::graph2ppt(file = file.path(dir_plots_eda, "Box_Simpson_UI2.pptx"),
                  x = p_simp_ui2, width = 7, height = 4.5)

## 4.3 UI2 × Region 基本统计（方便检查某些组合是否极端）
summary_by_group <- Sites_clean %>%
  dplyr::group_by(UI2, UN_region) %>%
  dplyr::summarise(
    n      = dplyr::n(),
    SR_med = stats::median(Species_richness),
    SR_max = max(Species_richness),
    Ab_med = stats::median(LogAbund),
    Ab_max = max(LogAbund),
    Sim_med = stats::median(Simpson_diversity),
    Sim_max = max(Simpson_diversity),
    .groups = "drop"
  )

cat("\n=== UI2 × UN_region summaries (cleaned) ===\n")
print(summary_by_group)

write.csv(summary_by_group,
          file.path(dir_data_diag, "Summary_by_UI2_region.csv"),
          row.names = FALSE)

##----------------------------------------------------------##
## 5. 样本量 & 低样本组合 / Sample sizes & low-N combos
##----------------------------------------------------------##

sr_n <- Sites_clean %>%
  dplyr::count(UI2, UN_region, name = "n")

cat("\n=== UI2 × UN_region sample sizes (cleaned) ===\n")
print(sr_n)

low_groups <- sr_n %>% dplyr::filter(n <= 40)

cat("\n=== Low-sample UI2 × UN_region combinations (n ≤ 40) ===\n")
print(low_groups)

write.csv(sr_n,
          file.path(dir_data_diag, "UI2_region_sample_sizes.csv"),
          row.names = FALSE)
write.csv(low_groups,
          file.path(dir_data_diag, "UI2_region_lowSampleGroups_n40.csv"),
          row.names = FALSE)

##----------------------------------------------------------##
## 6. 全数据 vs 剔除小样本数据集 / Full vs filtered datasets
##----------------------------------------------------------##

Sites_full <- Sites_clean
Sites_filt <- Sites_clean %>%
  dplyr::anti_join(low_groups, by = c("UI2","UN_region"))

cat("\nData size (full):", nrow(Sites_full), "\n")
cat("Data size (filtered n>40):", nrow(Sites_filt), "\n")

cat("\nFull UI2 × UN_region table:\n")
print(table(Sites_full$UI2, Sites_full$UN_region))
cat("\nFiltered UI2 × UN_region table:\n")
print(table(Sites_filt$UI2, Sites_filt$UN_region))

##----------------------------------------------------------##
## 7. 模型数据 / Model data frames
##----------------------------------------------------------##

prep_model_data <- function(Sites_obj) {
  md_sr <- stats::na.omit(Sites_obj[, c("Species_richness","UI2","SS","SSB","SSBS","UN_region")])
  md_ab <- stats::na.omit(Sites_obj[, c("LogAbund","UI2","SS","SSB","SSBS","UN_region")])
  md_sd <- stats::na.omit(Sites_obj[, c("Simpson_diversity","UI2","SS","SSB","SSBS","UN_region")])
  md_sd$LogSimpson <- log(md_sd$Simpson_diversity + 1)
  
  md_sr$UI2 <- factor(md_sr$UI2, levels = ui2_levels)
  md_ab$UI2 <- factor(md_ab$UI2, levels = ui2_levels)
  md_sd$UI2 <- factor(md_sd$UI2, levels = ui2_levels)
  
  md_sr$UN_region <- droplevels(factor(md_sr$UN_region))
  md_ab$UN_region <- droplevels(factor(md_ab$UN_region))
  md_sd$UN_region <- droplevels(factor(md_sd$UN_region))
  
  list(sr = md_sr, ab = md_ab, sd = md_sd)
}

md_full <- prep_model_data(Sites_full)
md_filt <- prep_model_data(Sites_filt)

model_full_sr <- md_full$sr
model_full_ab <- md_full$ab
model_full_sd <- md_full$sd

model_filt_sr <- md_filt$sr
model_filt_ab <- md_filt$ab
model_filt_sd <- md_filt$sd

##----------------------------------------------------------##
## 8. lme4 模型：全数据 & 过滤数据 / lme4 models
##----------------------------------------------------------##

cat("\n=== Fitting lme4 models: FULL data ===\n")

sr_full_UI2 <- glmer(
  Species_richness ~ UI2 + (1|SS) + (1|SSB) + (1|SSBS),
  data   = model_full_sr,
  family = poisson(link = "log")
)

sr_full_UI2reg <- glmer(
  Species_richness ~ UI2 * UN_region + (1|SS) + (1|SSB) + (1|SSBS),
  data   = model_full_sr,
  family = poisson(link = "log")
)

ab_full_UI2 <- lmer(
  LogAbund ~ UI2 + (1|SS) + (1|SSB),
  data = model_full_ab,
  REML = FALSE
)

ab_full_UI2reg <- lmer(
  LogAbund ~ UI2 * UN_region + (1|SS) + (1|SSB),
  data = model_full_ab,
  REML = FALSE
)

sd_full_UI2 <- lmer(
  LogSimpson ~ UI2 + (1|SS) + (1|SSB),
  data = model_full_sd,
  REML = FALSE
)

sd_full_UI2reg <- lmer(
  LogSimpson ~ UI2 * UN_region + (1|SS) + (1|SSB),
  data = model_full_sd,
  REML = FALSE
)

cat("\n=== Fitting lme4 models: FILTERED data ===\n")

sr_filt_UI2 <- glmer(
  Species_richness ~ UI2 + (1|SS) + (1|SSB) + (1|SSBS),
  data   = model_filt_sr,
  family = poisson(link = "log")
)

sr_filt_UI2reg <- glmer(
  Species_richness ~ UI2 * UN_region + (1|SS) + (1|SSB) + (1|SSBS),
  data   = model_filt_sr,
  family = poisson(link = "log")
)

ab_filt_UI2 <- lmer(
  LogAbund ~ UI2 + (1|SS) + (1|SSB),
  data = model_filt_ab,
  REML = FALSE
)

ab_filt_UI2reg <- lmer(
  LogAbund ~ UI2 * UN_region + (1|SS) + (1|SSB),
  data = model_filt_ab,
  REML = FALSE
)

sd_filt_UI2 <- lmer(
  LogSimpson ~ UI2 + (1|SS) + (1|SSB),
  data = model_filt_sd,
  REML = FALSE
)

sd_filt_UI2reg <- lmer(
  LogSimpson ~ UI2 * UN_region + (1|SS) + (1|SSB),
  data = model_filt_sd,
  REML = FALSE
)

cat("\n=== Model summary (filtered SR: Species_richness ~ UI2 * UN_region) ===\n")
print(summary(sr_filt_UI2reg))

## 保存 lme4 模型 / Save lme4 models
saveRDS(sr_full_UI2,    file.path(dir_models_glmer, "sr_full_UI2.rds"))
saveRDS(sr_full_UI2reg, file.path(dir_models_glmer, "sr_full_UI2reg.rds"))
saveRDS(sr_filt_UI2,    file.path(dir_models_glmer, "sr_filt_UI2.rds"))
saveRDS(sr_filt_UI2reg, file.path(dir_models_glmer, "sr_filt_UI2reg.rds"))

saveRDS(ab_full_UI2,    file.path(dir_models_glmer, "ab_full_UI2.rds"))
saveRDS(ab_full_UI2reg, file.path(dir_models_glmer, "ab_full_UI2reg.rds"))
saveRDS(ab_filt_UI2,    file.path(dir_models_glmer, "ab_filt_UI2.rds"))
saveRDS(ab_filt_UI2reg, file.path(dir_models_glmer, "ab_filt_UI2reg.rds"))

saveRDS(sd_full_UI2,    file.path(dir_models_glmer, "sd_full_UI2.rds"))
saveRDS(sd_full_UI2reg, file.path(dir_models_glmer, "sd_full_UI2reg.rds"))
saveRDS(sd_filt_UI2,    file.path(dir_models_glmer, "sd_filt_UI2.rds"))
saveRDS(sd_filt_UI2reg, file.path(dir_models_glmer, "sd_filt_UI2reg.rds"))

##----------------------------------------------------------##
## 9. 过度离散 & DHARMa 诊断：全数据 vs 过滤数据
##    Overdispersion & diagnostics (full vs filtered)
##----------------------------------------------------------##
## 诊断逻辑（写在注释中方便稿件）：
##  - SR：Poisson 模型常见过度离散或零膨胀 → 看 Pearson χ²/df、DHARMa
##  - 丰度、Simpson 的 Gaussian(Log) 确认残差是否接近正态、方差齐性
##  - 全数据 vs 过滤数据：若过滤后残差更好、异常点更少，说明剔除
##    极低样本组合是合理的 sensitivity check
##----------------------------------------------------------##

overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp  <- residuals(model, type = "pearson")
  Pearson.chisq <- sum(rp^2)
  ratio <- Pearson.chisq / rdf
  pval  <- stats::pchisq(Pearson.chisq, df = rdf, lower.tail = FALSE)
  c(chisq = Pearson.chisq, ratio = ratio, rdf = rdf, p = pval)
}

## 9.1 Poisson SR: full vs filtered
od_full  <- overdisp_fun(sr_full_UI2reg)
od_filt  <- overdisp_fun(sr_filt_UI2reg)

cat("\n=== Overdispersion test (Poisson SR, full vs filtered) ===\n")
print(rbind(full = od_full, filtered = od_filt))
write.table(rbind(full = od_full, filtered = od_filt),
            file.path(dir_data_diag, "Overdispersion_SR_full_vs_filtered.txt"),
            col.names = NA, quote = FALSE)

## 判读：ratio >> 1 且 p<0.05 = 过度离散；ratio<<1 且 p<0.05 = under-dispersion

## 9.2 DHARMa: SR (filtered) —— 主推图
set.seed(123)
sim_sr_f <- DHARMa::simulateResiduals(sr_filt_UI2reg, n = 1000)

## 屏幕显示一次（方便交互查看）/ Show on screen
plot(sim_sr_f, main = "SR ~ UI2 × Region (filtered, Poisson)")

## 保存 PNG & PDF（必须 dev.off()）
png(file.path(dir_plots_dharma, "DHARMa_SR_filtered.png"),
    width = 7, height = 7, units = "in", res = 300)
plot(sim_sr_f, main = "SR ~ UI2 × Region (filtered, Poisson)")
dev.off()

pdf(file.path(dir_plots_dharma, "DHARMa_SR_filtered.pdf"),
    width = 7, height = 7)
plot(sim_sr_f, main = "SR ~ UI2 × Region (filtered, Poisson)")
dev.off()

## 再画一次并导出 PPTX（graph2ppt 捕获当前图形窗口）
plot(sim_sr_f, main = "SR ~ UI2 × Region (filtered, Poisson)")
export::graph2ppt(file = file.path(dir_plots_dharma, "DHARMa_SR_filtered.pptx"),
                  width = 7, height = 7)

## 数值检验结果写入文件 + 注释判断依据
sink(file.path(dir_data_diag, "DHARMa_tests_SR_Poisson_filtered.txt"))
cat("Uniformity test (should be p>0.05 for good fit):\n")
print(DHARMa::testUniformity(sim_sr_f))
cat("\nDispersion test (p<0.05 → over/under-dispersion):\n")
print(DHARMa::testDispersion(sim_sr_f))
cat("\nZero-inflation test (p<0.05 → zero inflation/deflation):\n")
print(DHARMa::testZeroInflation(sim_sr_f))
sink()

## 9.3 DHARMa: Abundance Gaussian (filtered)
set.seed(123)
sim_ab_f <- DHARMa::simulateResiduals(ab_filt_UI2reg, n = 1000)
plot(sim_ab_f, main = "LogAbund ~ UI2 × Region (filtered, Gaussian)")

png(file.path(dir_plots_dharma, "DHARMa_Abundance_filtered.png"),
    width = 7, height = 7, units = "in", res = 300)
plot(sim_ab_f, main = "LogAbund ~ UI2 × Region (filtered, Gaussian)")
dev.off()

pdf(file.path(dir_plots_dharma, "DHARMa_Abundance_filtered.pdf"),
    width = 7, height = 7)
plot(sim_ab_f, main = "LogAbund ~ UI2 × Region (filtered, Gaussian)")
dev.off()

plot(sim_ab_f, main = "LogAbund ~ UI2 × Region (filtered, Gaussian)")
export::graph2ppt(file = file.path(dir_plots_dharma, "DHARMa_Abundance_filtered.pptx"),
                  width = 7, height = 7)

sink(file.path(dir_data_diag, "DHARMa_tests_Abundance_Gaussian_filtered.txt"))
cat("Uniformity test (residual distribution):\n")
print(DHARMa::testUniformity(sim_ab_f))
cat("\nDispersion test:\n")
print(DHARMa::testDispersion(sim_ab_f))
sink()

## 9.4 DHARMa: Simpson Gaussian (filtered)
set.seed(123)
sim_sd_f <- DHARMa::simulateResiduals(sd_filt_UI2reg, n = 1000)
plot(sim_sd_f, main = "LogSimpson ~ UI2 × Region (filtered, Gaussian)")

png(file.path(dir_plots_dharma, "DHARMa_Simpson_filtered.png"),
    width = 7, height = 7, units = "in", res = 300)
plot(sim_sd_f, main = "LogSimpson ~ UI2 × Region (filtered, Gaussian)")
dev.off()

pdf(file.path(dir_plots_dharma, "DHARMa_Simpson_filtered.pdf"),
    width = 7, height = 7)
plot(sim_sd_f, main = "LogSimpson ~ UI2 × Region (filtered, Gaussian)")
dev.off()

plot(sim_sd_f, main = "LogSimpson ~ UI2 × Region (filtered, Gaussian)")
export::graph2ppt(file = file.path(dir_plots_dharma, "DHARMa_Simpson_filtered.pptx"),
                  width = 7, height = 7)

#sink(file.path(dir_data_diag, "DHARMa_tests_Simpson_Gaussian_filtered.txt"))
cat("Uniformity test:\n")
print(DHARMa::testUniformity(sim_sd_f))
cat("\nDispersion test:\n")
print(DHARMa::testDispersion(sim_sd_f))
#sink()

##----------------------------------------------------------##
## 10. 逆 link & Δ% 相对 Primary vegetation
##----------------------------------------------------------##

inv_link <- function(eta_mat, type = c("pois","log1p","id")) {
  type <- match.arg(type)
  switch(
    type,
    "pois"  = exp(eta_mat),
    "log1p" = exp(eta_mat) - 1,
    "id"    = eta_mat
  )
}

rel_to_primary <- function(pred_mat) {
  rel <- sweep(pred_mat, 2, pred_mat[1,], "/")
  data.frame(
    median = (apply(rel, 1, stats::median)   * 100) - 100,
    lower  = (apply(rel, 1, stats::quantile, 0.025) * 100) - 100,
    upper  = (apply(rel, 1, stats::quantile, 0.975) * 100) - 100
  )
}

region_delta_from_mat <- function(pred_mat, newdata, regions, metric_label) {
  out_list <- vector("list", length(regions))
  names(out_list) <- regions
  
  for (r in regions) {
    idx_r  <- which(newdata$UN_region == r)
    idx_pv <- idx_r[newdata$UI2[idx_r] == "Primary vegetation"]
    preds_r <- pred_mat[idx_r, , drop = FALSE]
    rel     <- sweep(preds_r, 2, pred_mat[idx_pv, ], "/")
    med <- (apply(rel, 1, stats::median)   * 100) - 100
    lwr <- (apply(rel, 1, stats::quantile, 0.025) * 100) - 100
    upr <- (apply(rel, 1, stats::quantile, 0.975) * 100) - 100
    
    out_list[[r]] <- data.frame(
      UN_region = r,
      UI2       = newdata$UI2[idx_r],
      median    = med,
      lower     = lwr,
      upper     = upr,
      metric    = metric_label
    )
  }
  dplyr::bind_rows(out_list)
}

##----------------------------------------------------------##
## 11. global & region Δ%：全数据 vs 过滤数据
##----------------------------------------------------------##

build_effects <- function(sr_global, ab_global, sd_global,
                          sr_reg,    ab_reg,    sd_reg,
                          regions) {
  nd_global <- data.frame(
    UI2               = factor(ui2_levels, levels = ui2_levels),
    Species_richness  = 0,
    LogAbund          = 0,
    LogSimpson        = 0,
    Simpson_diversity = 0
  )
  
  sr_eta_g  <- PredictGLMERRandIter(sr_global, nd_global)
  sr_resp_g <- inv_link(sr_eta_g, "pois")
  sr_delta  <- rel_to_primary(sr_resp_g)
  rich_global <- cbind(
    UI2    = nd_global$UI2,
    sr_delta,
    metric = "Δ species richness (%)"
  )
  
  ab_eta_g  <- PredictGLMERRandIter(ab_global, nd_global)
  ab_resp_g <- inv_link(ab_eta_g, "log1p")
  ab_delta  <- rel_to_primary(ab_resp_g)
  abun_global <- cbind(
    UI2    = nd_global$UI2,
    ab_delta,
    metric = "Δ total abundance (%)"
  )
  
  sd_eta_g  <- PredictGLMERRandIter(sd_global, nd_global)
  sd_resp_g <- inv_link(sd_eta_g, "log1p")
  sd_delta  <- rel_to_primary(sd_resp_g)
  simp_global <- cbind(
    UI2    = nd_global$UI2,
    sd_delta,
    metric = "Δ Simpson’s diversity (%)"
  )
  
  global_all <- dplyr::bind_rows(rich_global, abun_global, simp_global)
  global_all[global_all$UI2 == "Primary vegetation", c("lower","upper")] <- NA
  
  nd_reg <- expand.grid(
    UI2       = factor(ui2_levels, levels = ui2_levels),
    UN_region = regions
  )
  nd_reg$Species_richness  <- 0
  nd_reg$LogAbund          <- 0
  nd_reg$LogSimpson        <- 0
  nd_reg$Simpson_diversity <- 0
  
  sr_eta_reg  <- PredictGLMERRandIter(sr_reg, nd_reg)
  sr_resp_reg <- inv_link(sr_eta_reg, "pois")
  rich_region <- region_delta_from_mat(
    pred_mat     = sr_resp_reg,
    newdata      = nd_reg,
    regions      = regions,
    metric_label = "Δ species richness (%)"
  )
  
  ab_eta_reg  <- PredictGLMERRandIter(ab_reg, nd_reg)
  ab_resp_reg <- inv_link(ab_eta_reg, "log1p")
  abun_region <- region_delta_from_mat(
    pred_mat     = ab_resp_reg,
    newdata      = nd_reg,
    regions      = regions,
    metric_label = "Δ total abundance (%)"
  )
  
  sd_eta_reg  <- PredictGLMERRandIter(sd_reg, nd_reg)
  sd_resp_reg <- inv_link(sd_eta_reg, "log1p")
  simp_region <- region_delta_from_mat(
    pred_mat     = sd_resp_reg,
    newdata      = nd_reg,
    regions      = regions,
    metric_label = "Δ Simpson’s diversity (%)"
  )
  
  list(
    global = global_all,
    region = dplyr::bind_rows(rich_region, abun_region, simp_region)
  )
}

regions_full <- levels(model_full_sr$UN_region)
regions_filt <- levels(model_filt_sr$UN_region)

eff_full <- build_effects(
  sr_global = sr_full_UI2,
  ab_global = ab_full_UI2,
  sd_global = sd_full_UI2,
  sr_reg    = sr_full_UI2reg,
  ab_reg    = ab_full_UI2reg,
  sd_reg    = sd_full_UI2reg,
  regions   = regions_full
)

eff_filt <- build_effects(
  sr_global = sr_filt_UI2,
  ab_global = ab_filt_UI2,
  sd_global = sd_filt_UI2,
  sr_reg    = sr_filt_UI2reg,
  ab_reg    = ab_filt_UI2reg,
  sd_reg    = sd_filt_UI2reg,
  regions   = regions_filt
)

global_full <- eff_full$global
region_full <- eff_full$region
global_filt <- eff_filt$global
region_filt <- eff_filt$region

write.csv(global_full,
          file.path(dir_data_effects, "Global_effects_FULL.csv"),
          row.names = FALSE)
write.csv(region_full,
          file.path(dir_data_effects, "Region_effects_FULL.csv"),
          row.names = FALSE)
write.csv(global_filt,
          file.path(dir_data_effects, "Global_effects_FILTERED.csv"),
          row.names = FALSE)
write.csv(region_filt,
          file.path(dir_data_effects, "Region_effects_FILTERED.csv"),
          row.names = FALSE)

##----------------------------------------------------------##
## 12. 森林图数据，关注异常预测（置信区间过宽）
##----------------------------------------------------------##

nonprim_levels <- c("Secondary vegetation",
                    "Agriculture_Low",
                    "Agriculture_High",
                    "Urban")

ui2_labels <- c("Secondary vegetation",
                "Agriculture (low)",
                "Agriculture (high)",
                "Urban")

metric_levels <- c("Δ total abundance (%)",
                   "Δ species richness (%)",
                   "Δ Simpson’s diversity (%)")

prep_plot_df <- function(global_all, region_all) {
  region_all2 <- region_all %>%
    dplyr::filter(UI2 %in% nonprim_levels) %>%
    dplyr::mutate(
      UI2      = factor(UI2, levels = nonprim_levels),
      UI2_plot = factor(UI2, levels = nonprim_levels, labels = ui2_labels),
      metric   = factor(metric, levels = metric_levels),
      CI_width = upper - lower           # 置信区间宽度，方便后面检查异常
    )
  
  global_all2 <- global_all %>%
    dplyr::filter(UI2 %in% nonprim_levels) %>%
    dplyr::mutate(
      UI2      = factor(UI2, levels = nonprim_levels),
      UI2_plot = factor(UI2, levels = nonprim_levels, labels = ui2_labels),
      metric   = factor(metric, levels = metric_levels)
    )
  
  list(global = global_all2, region = region_all2)
}

plot_full  <- prep_plot_df(global_full, region_full)
plot_filt  <- prep_plot_df(global_filt, region_filt)

global_full_plot <- plot_full$global
region_full_plot <- plot_full$region
global_filt_plot <- plot_filt$global
region_filt_plot <- plot_filt$region

## 12.1 检查哪些组合 CI 特别宽（通常对应样本量少或异质性高）
ci_flag_full <- region_full_plot %>%
  dplyr::group_by(metric) %>%
  dplyr::mutate(
    CI_q90 = stats::quantile(CI_width, 0.9),
    wide_CI = CI_width > CI_q90
  ) %>%
  dplyr::ungroup() %>%
  dplyr::filter(wide_CI)

write.csv(ci_flag_full,
          file.path(dir_data_diag, "WideCI_combinations_FULL.csv"),
          row.names = FALSE)
#sink()
ci_flag_filt <- region_filt_plot %>%
  dplyr::group_by(metric) %>%
  dplyr::mutate(
    CI_q90 = stats::quantile(CI_width, 0.9),
    wide_CI = CI_width > CI_q90
  ) %>%
  dplyr::ungroup() %>%
  dplyr::filter(wide_CI)

write.csv(ci_flag_filt,
          file.path(dir_data_diag, "WideCI_combinations_FILTERED.csv"),
          row.names = FALSE)

## 解释：稿件中可以说 “we flagged the 10% widest CIs as potentially unstable
## predictions, which mostly occurred in low-sample land-use–region combinations”.

##----------------------------------------------------------##
## 13. 森林图（投稿版版式，适合主文 Figure）
##----------------------------------------------------------##

region_cols <- c(
  "Africa"        = "#1b9e77",
  "Europe"        = "#d95f02",
  "Asia"          = "#7570b3",
  "North America" = "#e7298a",
  "South America" = "#66a61e",
  "Oceania"       = "#e6ab02"
)

region_shapes <- c(
  "Africa"        = 16,
  "Europe"        = 17,
  "Asia"          = 15,
  "North America" = 3,
  "South America" = 4,
  "Oceania"       = 8
)

make_forest_plot <- function(region_plot, global_plot, Sites_obj, tag) {
  n_df <- Sites_obj %>%
    dplyr::filter(UI2 %in% nonprim_levels) %>%
    dplyr::count(UN_region, UI2, name = "n")
  
  region_all2 <- region_plot %>%
    dplyr::left_join(n_df, by = c("UN_region","UI2")) %>%
    dplyr::filter(!is.na(n)) %>%
    dplyr::mutate(
      UN_region = factor(
        UN_region,
        levels = c("Africa","Asia","Europe",
                   "North America","Oceania","South America")
      ),
      UI2      = factor(UI2,      levels = nonprim_levels),
      UI2_plot = factor(UI2_plot, levels = ui2_labels),
      metric   = factor(metric,   levels = metric_levels),
      x_num    = as.numeric(UN_region)
    )
  
  global_all2 <- global_plot %>%
    dplyr::filter(UI2 %in% nonprim_levels) %>%
    dplyr::mutate(
      UI2      = factor(UI2,      levels = nonprim_levels),
      UI2_plot = factor(UI2_plot, levels = ui2_labels),
      metric   = factor(metric,   levels = metric_levels)
    )
  
  y_min <- -100
  y_max <-  100
  effect_breaks <- seq(y_min, y_max, by = 25)
  
  max_n   <- max(region_all2$n, na.rm = TRUE)
  scale_n <- (y_max - y_min) / max_n
  
  n_breaks_raw <- c(0, max_n/3, 2*max_n/3, max_n)
  n_breaks <- round(n_breaks_raw, -1)
  
  p <- ggplot(region_all2, aes(x = UN_region)) +
    
    ## 样本量灰柱 / grey bars for sample size
    geom_rect(
      data = region_all2,
      aes(xmin = x_num - 0.35,
          xmax = x_num + 0.35,
          ymin = y_min,
          ymax = y_min + n * scale_n),
      inherit.aes = FALSE,
      fill = "grey95",
      colour = NA
    ) +
    
    ## global UI2-only 背景带 + 中位数虚线
    geom_rect(
      data        = global_all2,
      aes(xmin = -Inf, xmax = Inf, ymin = lower, ymax = upper),
      inherit.aes = FALSE,
      fill        = "grey88",
      colour      = NA
    ) +
    geom_hline(
      data        = global_all2,
      aes(yintercept = median),
      inherit.aes = FALSE,
      linetype    = "dotted",
      colour      = "grey40"
    ) +
    
    geom_hline(yintercept = 0, linetype = "dashed", size = 0.3) +
    
    geom_errorbar(
      aes(ymin = lower, ymax = upper, colour = UN_region),
      width = 0, size = 0.4
    ) +
    geom_point(
      aes(y = median, colour = UN_region, shape = UN_region),
      size = 1.8, stroke = 0.2
    ) +
    
    facet_grid(metric ~ UI2_plot, scales = "fixed") +
    
    scale_colour_manual(
      name   = "Region",
      values = region_cols,
      drop   = FALSE
    ) +
    scale_shape_manual(
      name   = "Region",
      values = region_shapes,
      drop   = FALSE
    ) +
    
    scale_y_continuous(
      limits   = c(y_min, y_max),
      breaks   = effect_breaks,
      name     = "Change relative to primary vegetation (%)",
      sec.axis = sec_axis(
        trans  = ~ (. - y_min) / scale_n,
        breaks = n_breaks,
        name   = "Sample size (n)"
      )
    ) +
    
    labs(
      x     = NULL,
      title = tag
    ) +
    
    theme_bw(base_size = 9) +
    theme(
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(size = 0.2, colour = "grey90"),
      axis.text.x        = element_text(angle = 45, hjust = 1, size = 7),
      axis.text.y        = element_text(size = 7),
      strip.text         = element_text(size = 8, face = "bold"),
      strip.background   = element_rect(fill = "white"),
      legend.position    = "top",
      legend.title       = element_text(size = 8),
      legend.text        = element_text(size = 7),
      aspect.ratio       = 0.8,
      panel.spacing.y    = unit(1.2, "lines"),
      plot.title         = element_text(size = 9, face = "bold", hjust = 0)
    )
  
  p
}

p_forest_full <- make_forest_plot(
  region_plot = region_full_plot,
  global_plot = global_full_plot,
  Sites_obj   = Sites_full,
  tag         = "A  Full data (all sites)"
)

p_forest_filt <- make_forest_plot(
  region_plot = region_filt_plot,
  global_plot = global_filt_plot,
  Sites_obj   = Sites_filt,
  tag         = "B  Filtered data (n > 40 per UI2 × region)"
)

print(p_forest_full)
print(p_forest_filt)

## 保存为投稿图版：PNG + PDF + PPTX
ggsave(file.path(dir_plots_forest, "Figure_LUI_region_FULL.png"),
       p_forest_full, width = 9, height = 8, units = "in", dpi = 600)
ggsave(file.path(dir_plots_forest, "Figure_LUI_region_FULL.pdf"),
       p_forest_full, width = 9, height = 8)
export::graph2ppt(file = file.path(dir_plots_forest, "Figure_LUI_region_FULL.pptx"),
                  x = p_forest_full, width = 9, height = 8)

ggsave(file.path(dir_plots_forest, "Figure_LUI_region_FILTERED.png"),
       p_forest_filt, width = 9, height = 8, units = "in", dpi = 600)
ggsave(file.path(dir_plots_forest, "Figure_LUI_region_FILTERED.pdf"),
       p_forest_filt, width = 9, height = 8)
export::graph2ppt(file = file.path(dir_plots_forest, "Figure_LUI_region_FILTERED.pptx"),
                  x = p_forest_filt, width = 9, height = 8)

##==========================================================##
## 14. 低样本组合 vs Global 箱线图（诊断 + 说明过滤合理性）
##==========================================================##

### 14.1 SR
low_data_diag_sr <- Sites_clean %>%
  dplyr::semi_join(low_groups, by = c("UI2","UN_region")) %>%
  dplyr::select(UI2, UN_region, Species_richness) %>%
  dplyr::mutate(group = paste(UI2, UN_region, sep = " | "))

global_diag_sr <- Sites_clean %>%
  dplyr::mutate(group = "Global") %>%
  dplyr::select(group, Species_richness)

compare_sr <- dplyr::bind_rows(global_diag_sr, low_data_diag_sr)

n_df_sr <- compare_sr %>%
  dplyr::group_by(group) %>%
  dplyr::summarise(
    n     = dplyr::n(),
    y_pos = max(Species_richness, na.rm = TRUE) * 1.05,
    .groups = "drop"
  )

p_sr_diag <- ggplot(compare_sr,
                    aes(x = group, y = Species_richness, fill = group)) +
  geom_boxplot(outlier.alpha = 0.5) +
  geom_text(
    data = n_df_sr,
    aes(x = group, y = y_pos, label = paste0("n = ", n)),
    inherit.aes = FALSE,
    vjust = 0, size = 3
  ) +
  theme_bw() +
  labs(x = "Group",
       y = "Species richness",
       title = "Species richness: Global vs low-sample UI2 × region groups") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

print(p_sr_diag)

ggsave(file.path(dir_plots_lowN, "SR_lowGroups_vs_global.png"),
       p_sr_diag, width = 8, height = 5, units = "in", dpi = 600)
ggsave(file.path(dir_plots_lowN, "SR_lowGroups_vs_global.pdf"),
       p_sr_diag, width = 8, height = 5)
export::graph2ppt(file = file.path(dir_plots_lowN, "SR_lowGroups_vs_global.pptx"),
                  x = p_sr_diag, width = 8, height = 5)

### 14.2 Abundance
low_data_diag_ab <- Sites_clean %>%
  dplyr::semi_join(low_groups, by = c("UI2","UN_region")) %>%
  dplyr::select(UI2, UN_region, LogAbund) %>%
  dplyr::mutate(group = paste(UI2, UN_region, sep = " | "))

global_diag_ab <- Sites_clean %>%
  dplyr::mutate(group = "Global") %>%
  dplyr::select(group, LogAbund)

compare_ab <- dplyr::bind_rows(global_diag_ab, low_data_diag_ab)

n_df_ab <- compare_ab %>%
  dplyr::group_by(group) %>%
  dplyr::summarise(
    n     = dplyr::n(),
    y_pos = max(LogAbund, na.rm = TRUE) * 1.05,
    .groups = "drop"
  )

p_ab_diag <- ggplot(compare_ab,
                    aes(x = group, y = LogAbund, fill = group)) +
  geom_boxplot(outlier.alpha = 0.5) +
  geom_text(
    data = n_df_ab,
    aes(x = group, y = y_pos, label = paste0("n = ", n)),
    inherit.aes = FALSE,
    vjust = 0, size = 3
  ) +
  theme_bw() +
  labs(x = "Group",
       y = "Total abundance (log scale)",
       title = "Total abundance: Global vs low-sample UI2 × region groups") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

print(p_ab_diag)

ggsave(file.path(dir_plots_lowN, "Abundance_lowGroups_vs_global.png"),
       p_ab_diag, width = 8, height = 5, units = "in", dpi = 600)
ggsave(file.path(dir_plots_lowN, "Abundance_lowGroups_vs_global.pdf"),
       p_ab_diag, width = 8, height = 5)
export::graph2ppt(file = file.path(dir_plots_lowN, "Abundance_lowGroups_vs_global.pptx"),
                  x = p_ab_diag, width = 8, height = 5)

### 14.3 Simpson
low_data_diag_sd <- Sites_clean %>%
  dplyr::semi_join(low_groups, by = c("UI2","UN_region")) %>%
  dplyr::select(UI2, UN_region, Simpson_diversity) %>%
  dplyr::mutate(group = paste(UI2, UN_region, sep = " | "))

global_diag_sd <- Sites_clean %>%
  dplyr::mutate(group = "Global") %>%
  dplyr::select(group, Simpson_diversity)

compare_simpson <- dplyr::bind_rows(global_diag_sd, low_data_diag_sd)

n_df_simpson <- compare_simpson %>%
  dplyr::group_by(group) %>%
  dplyr::summarise(
    n     = dplyr::n(),
    y_pos = max(Simpson_diversity, na.rm = TRUE) * 1.05,
    .groups = "drop"
  )

p_simpson_diag <- ggplot(compare_simpson,
                         aes(x = group, y = Simpson_diversity, fill = group)) +
  geom_boxplot(outlier.alpha = 0.5) +
  geom_text(
    data = n_df_simpson,
    aes(x = group, y = y_pos, label = paste0("n = ", n)),
    inherit.aes = FALSE,
    vjust = 0, size = 3
  ) +
  theme_bw() +
  labs(x = "Group",
       y = "Simpson diversity",
       title = "Simpson diversity: Global vs low-sample UI2 × region groups") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

print(p_simpson_diag)

ggsave(file.path(dir_plots_lowN, "Simpson_lowGroups_vs_global.png"),
       p_simpson_diag, width = 8, height = 5, units = "in", dpi = 600)
ggsave(file.path(dir_plots_lowN, "Simpson_lowGroups_vs_global.pdf"),
       p_simpson_diag, width = 8, height = 5)
export::graph2ppt(file = file.path(dir_plots_lowN, "Simpson_lowGroups_vs_global.pptx"),
                  x = p_simpson_diag, width = 8, height = 5)

##==========================================================##
## 15. brms 模型：全数据 & 过滤数据 + 投稿级边际效应图
##==========================================================##

backend_brms <- "rstan"
n_cores <- max(parallel::detectCores() - 1, 1)
options(mc.cores = n_cores)
cat("Using", n_cores, "cores for brms.\n")

b_prior_sr_pois <- c(
  brms::prior(normal(0, 2), class = "b"),
  brms::prior(student_t(3, 0, 2.5), class = "sd")
)

b_prior_sr_nb <- c(
  brms::prior(normal(0, 2), class = "b"),
  brms::prior(student_t(3, 0, 2.5), class = "sd"),
  brms::prior(normal(0, 1), class = "shape", lb = 0)
)
####先验调整######
b_prior_sr_nb <- c(
  prior(normal(0, 2), class = "b"),
  prior(student_t(3, 0, 2.5), class = "sd"),
  prior(gamma(2, 0.1), class = "shape")
)

b_prior_gaus <- c(
  brms::prior(normal(0, 2), class = "b"),
  brms::prior(student_t(3, 0, 2.5), class = "sd"),
  brms::prior(student_t(3, 0, 2.5), class = "sigma")
)

brm_ctrl <- list(adapt_delta = 0.95, max_treedepth = 12)

## 15.1 Species richness
b_sr_full_pois <- brm(
  Species_richness ~ UI2 * UN_region + (1|SS) + (1|SSB) + (1|SSBS),
  data    = model_full_sr,
  family  = poisson(),
  prior   = b_prior_sr_pois,
  chains  = 4, iter = 4000, warmup = 2000,
  control = brm_ctrl,
  seed    = 2001,
  backend = backend_brms
)
b_sr_full_nb <- brm(
  Species_richness ~ UI2 * UN_region + (1|SS) + (1|SSB) + (1|SSBS),
  data    = model_full_sr,
  family  = negbinomial(),          # <- 关键：NB
  prior   = b_prior_sr_nb,          # <- 应包含 shape 先验
  chains  = 4, iter = 6000, warmup = 3000,
  cores   = 4,                      # <- 建议显式写
  backend = backend_brms,           # "cmdstanr"
  threads = threading(1),           # <- 先用 1 最稳
  init    = 0,                      # <- 避免慢链
  control = brm_ctrl,
  seed    = 2001
)

b_sr_filt_nb <- brm(
  Species_richness ~ UI2 * UN_region + (1|SS) + (1|SSB) + (1|SSBS),
  data    = model_filt_sr,
  family  = negbinomial(link = "log"),
  prior   = b_prior_sr_nb,
  chains  = 4,
  iter    = 6000,
  warmup  = 3000,
  control = list(adapt_delta = 0.95, max_treedepth = 12),
  backend = backend_brms,
  seed    = 2002
)

cat("\n=== brms SR models (full Poisson vs filtered NB) ===\n")
print(summary(b_sr_full_nb), digits = 2)
print(summary(b_sr_filt_nb),   digits = 2)

saveRDS(b_sr_full_nb, file.path(dir_models_brms, "b_sr_full_nb.rds"))
saveRDS(b_sr_filt_nb,   file.path(dir_models_brms, "b_sr_filt_nb.rds"))

## pp_check：比较 full Poisson vs filtered NB
p_pp_sr_full <- pp_check(b_sr_full_nb, nsamples = 200)
p_pp_sr_filt <- pp_check(b_sr_filt_nb,   nsamples = 200)

ggsave(file.path(dir_plots_brms, "ppcheck_SR_full_Poisson.png"),
       p_pp_sr_full, width = 6, height = 4, dpi = 400)
ggsave(file.path(dir_plots_brms, "ppcheck_SR_full_Poisson.pdf"),
       p_pp_sr_full, width = 6, height = 4)
export::graph2ppt(file = file.path(dir_plots_brms, "ppcheck_SR_full_Poisson.pptx"),
                  x = p_pp_sr_full, width = 6, height = 4)

ggsave(file.path(dir_plots_brms, "ppcheck_SR_filtered_NB.png"),
       p_pp_sr_filt, width = 6, height = 4, dpi = 400)
ggsave(file.path(dir_plots_brms, "ppcheck_SR_filtered_NB.pdf"),
       p_pp_sr_filt, width = 6, height = 4)
export::graph2ppt(file = file.path(dir_plots_brms, "ppcheck_SR_filtered_NB.pptx"),
                  x = p_pp_sr_filt, width = 6, height = 4)

## 15.2 边际效应图（投稿版）：filtered NB 模型
ce_sr_nb <- conditional_effects(b_sr_filt_nb, effects = "UI2:UN_region")

p_ce_sr_nb <- plot(ce_sr_nb, plot = FALSE)[[1]] +
  theme_bw(base_size = 9) +
  labs(
    x = "Land-use type",
    y = "Predicted species richness\n(negative binomial, filtered data)",
    title = "C  Species richness responses (brms NB, filtered data)"
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(size = 0.2, colour = "grey90"),
    axis.text.x      = element_text(angle = 45, hjust = 1),
    strip.text       = element_text(size = 8, face = "bold")
  )

print(p_ce_sr_nb)

ggsave(file.path(dir_plots_brms, "Figure_SR_brms_NB_filtered.png"),
       p_ce_sr_nb, width = 7.5, height = 6, units = "in", dpi = 600)
ggsave(file.path(dir_plots_brms, "Figure_SR_brms_NB_filtered.pdf"),
       p_ce_sr_nb, width = 7.5, height = 6)
export::graph2ppt(file = file.path(dir_plots_brms, "Figure_SR_brms_NB_filtered.pptx"),
                  x = p_ce_sr_nb, width = 7.5, height = 6)

## 这里你可以把森林图 (A–B) + 边际效应图 (C) 组合成 Figure 2/3：
##  A: 全数据森林图，展示潜在异常宽 CI（但说明仅用于对比）
##  B: 过滤数据森林图（主结果）
##  C: brms NB 边际效应（展示更平滑/稳健的 region × UI2 响应）

##==========================================================##
## 16. glmmTMB SR NB / ZI-NB + 诊断（与 lme4 / brms 对比）
##==========================================================##

sr_nb <- glmmTMB::glmmTMB(
  Species_richness ~ UI2 * UN_region +
    (1 | SS) + (1 | SSB) + (1 | SSBS),
  data   = model_filt_sr,
  family = glmmTMB::nbinom2()
)

sr_nb_zi <- glmmTMB::glmmTMB(
  Species_richness ~ UI2 * UN_region +
    (1 | SS) + (1 | SSB) + (1 | SSBS),
  ziformula = ~ 1,
  data      = model_filt_sr,
  family    = glmmTMB::nbinom2()
)

saveRDS(sr_nb,    file.path(dir_models_glmmTMB, "sr_nb_filtered.rds"))
saveRDS(sr_nb_zi, file.path(dir_models_glmmTMB, "sr_nb_ZINB_filtered.rds"))

aic_comp_sr <- AIC(sr_filt_UI2reg, sr_nb, sr_nb_zi)
cat("\n=== AIC: Poisson vs NB vs ZI-NB (SR, filtered) ===\n")
print(aic_comp_sr)
write.csv(as.data.frame(aic_comp_sr),
          file.path(dir_data_diag, "AIC_SR_Pois_vs_NB_vs_ZINB.csv"))

## 诊断：若 NB / ZI-NB 的 DHARMa 更接近均匀且 AIC 明显更低，
## 在稿件中可说明 “moving from Poisson to NB (and optionally ZI-NB)
##  greatly improved residual behaviour and model fit”.

set.seed(123)
sim_sr_nb    <- DHARMa::simulateResiduals(sr_nb,    n = 1000)
set.seed(123)
sim_sr_nb_zi <- DHARMa::simulateResiduals(sr_nb_zi, n = 1000)

## NB
plot(sim_sr_nb, main = "NB: SR ~ UI2 × Region (filtered)")
png(file.path(dir_plots_dharma, "DHARMa_SR_NB_filtered.png"),
    width = 7, height = 7, units = "in", res = 300)
plot(sim_sr_nb, main = "NB: SR ~ UI2 × Region (filtered)")
dev.off()
pdf(file.path(dir_plots_dharma, "DHARMa_SR_NB_filtered.pdf"),
    width = 7, height = 7)
plot(sim_sr_nb, main = "NB: SR ~ UI2 × Region (filtered)")
dev.off()
plot(sim_sr_nb, main = "NB: SR ~ UI2 × Region (filtered)")
export::graph2ppt(file = file.path(dir_plots_dharma, "DHARMa_SR_NB_filtered.pptx"),
                  width = 7, height = 7)

sink(file.path(dir_data_diag, "DHARMa_tests_SR_NB_filtered.txt"))
cat("Uniformity test:\n");      print(DHARMa::testUniformity(sim_sr_nb))
cat("\nDispersion test:\n");    print(DHARMa::testDispersion(sim_sr_nb))
cat("\nZero-inflation test:\n");print(DHARMa::testZeroInflation(sim_sr_nb))
sink()

## ZI-NB
plot(sim_sr_nb_zi, main = "ZI-NB: SR ~ UI2 × Region (filtered)")
png(file.path(dir_plots_dharma, "DHARMa_SR_ZINB_filtered.png"),
    width = 7, height = 7, units = "in", res = 300)
plot(sim_sr_nb_zi, main = "ZI-NB: SR ~ UI2 × Region (filtered)")
dev.off()
pdf(file.path(dir_plots_dharma, "DHARMa_SR_ZINB_filtered.pdf"),
    width = 7, height = 7)
plot(sim_sr_nb_zi, main = "ZI-NB: SR ~ UI2 × Region (filtered)")
dev.off()
plot(sim_sr_nb_zi, main = "ZI-NB: SR ~ UI2 × Region (filtered)")
export::graph2ppt(file = file.path(dir_plots_dharma, "DHARMa_SR_ZINB_filtered.pptx"),
                  width = 7, height = 7)

sink(file.path(dir_data_diag, "DHARMa_tests_SR_ZINB_filtered.txt"))
cat("Uniformity test:\n");      print(DHARMa::testUniformity(sim_sr_nb_zi))
cat("\nDispersion test:\n");    print(DHARMa::testDispersion(sim_sr_nb_zi))
cat("\nZero-inflation test:\n");print(DHARMa::testZeroInflation(sim_sr_nb_zi))
sink()

##==========================================================##
## 17. glmmTMB Abundance (NB on counts) + 诊断
##==========================================================##

model_filt_ab2 <- model_filt_ab %>%
  dplyr::mutate(
    Abund_raw = exp(LogAbund) - 1
  )

ab_nb <- glmmTMB::glmmTMB(
  Abund_raw ~ UI2 * UN_region + (1 | SS) + (1 | SSB),
  data   = model_filt_ab2,
  family = glmmTMB::nbinom2()
)

saveRDS(ab_nb, file.path(dir_models_glmmTMB, "ab_nb_filtered.rds"))

aic_comp_ab <- AIC(ab_filt_UI2reg, ab_nbcat("\n=== AIC: Gaussian(LogAbund) vs NB(Abund_raw) ===\n"))
print(aic_comp_ab)
write.csv(as.data.frame(aic_comp_ab),
          file.path(dir_data_diag, "AIC_Abundance_Gaussian_vs_NB.csv"))

set.seed(123)
sim_ab_nb <- DHARMa::simulateResiduals(ab_nb, n = 1000)

plot(sim_ab_nb, main = "NB: Abund_raw ~ UI2 × Region (filtered)")
png(file.path(dir_plots_dharma, "DHARMa_Abundance_NB_filtered.png"),
    width = 7, height = 7, units = "in", res = 300)
plot(sim_ab_nb, main = "NB: Abund_raw ~ UI2 × Region (filtered)")
dev.off()
pdf(file.path(dir_plots_dharma, "DHARMa_Abundance_NB_filtered.pdf"),
    width = 7, height = 7)
plot(sim_ab_nb, main = "NB: Abund_raw ~ UI2 × Region (filtered)")
dev.off()
plot(sim_ab_nb, main = "NB: Abund_raw ~ UI2 × Region (filtered)")
export::graph2ppt(file = file.path(dir_plots_dharma, "DHARMa_Abundance_NB_filtered.pptx"),
                  width = 7, height = 7)

sink(file.path(dir_data_diag, "DHARMa_tests_Abundance_NB_filtered.txt"))
cat("Uniformity test:\n");      print(DHARMa::testUniformity(sim_ab_nb))
cat("\nDispersion test:\n");    print(DHARMa::testDispersion(sim_ab_nb))
cat("\nZero-inflation test:\n");print(DHARMa::testZeroInflation(sim_ab_nb))
sink()

####===================  End of script  ===================####


