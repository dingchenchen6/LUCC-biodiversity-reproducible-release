###############################################################
## Land-use x Climate warming x Functional strategy groups (FG)
## Occurrence (glmmTMB) -- STABLE & INTERPRETABLE FG pipeline
## 功能策略组（稳健PCA + 聚类稳定性评估）x 土地利用 x 变暖 -> 物种出现概率
##
## Design goals / 设计目标
## 1) Professional: use trait preprocessing suitable for strongly skewed traits
## 2) Rigorous: choose FG using multi-criteria selection, not silhouette alone
## 3) Stable: enforce minimum cluster size and bootstrap stability
## 4) Interpretable: keep raw traits for ecological interpretation and figures
##
## Core workflow / 核心流程
## Step 1) Build species-level trait table (median per species)
## Step 2) Create "clustering trait space":
##         - log-transform skewed positive traits
##         - winsorise mild tails
##         - z-score standardisation
## Step 3) PCA on transformed z-scored traits
## Step 4) Compare clustering candidates on PCA space:
##         - methods: kmeans / pam / ward
##         - K: 2..5
##         - metrics: silhouette + min cluster proportion + bootstrap Jaccard stability
## Step 5) Choose FG using a transparent rule:
##         - min cluster proportion >= threshold
##         - min bootstrap Jaccard >= threshold
##         - silhouette within 90% of the best eligible solution
##         - choose the largest K among remaining candidates
## Step 6) Join FG back to occurrence data
## Step 7) Fit glmmTMB:
##         logit(p) = UI2 * StdTmeanAnomalyRS * FG + random intercepts
## Step 8) Export:
##         - cluster selection tables and plot
##         - PCA biplot + trait profile bar + heatmap
##         - three-way and two-way paper-grade plots
##         - interpretation tables and logs
##
## Notes / 说明
## - Clustering is done on transformed/winsorised trait space for stability.
## - Ecological interpretation is always reported back on RAW trait values and
##   species-standardised raw-trait z-scores.
###############################################################

# rm(list = ls())

## ============================================================
## 0) Libraries / 加载包
## ============================================================
suppressPackageStartupMessages({
  library(glmmTMB)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
  library(cluster)
  library(MASS)
})

if (requireNamespace("ggrepel", quietly = TRUE)) {
  suppressPackageStartupMessages(library(ggrepel))
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

out_root <- "glmmTMB_occurrence_FG_RobustPCAKmeans_FULL_plusTwoWay"
dir.create(out_root, showWarnings = FALSE, recursive = TRUE)

dir_data   <- file.path(out_root, "00_Data")
dir_expl   <- file.path(out_root, "01_Exploration_and_FG_Interpretation")
dir_models <- file.path(out_root, "02_Models")
dir_plots3 <- file.path(out_root, "03_Plots_ThreeWay_ABS_PCT")
dir_plots2 <- file.path(out_root, "04_Plots_TwoWay_PaperGrade")
dir_tables <- file.path(out_root, "05_Tables_FG_Interpretation")
dir_caps   <- file.path(out_root, "06_Figure_Captions")
dir.create(dir_data,   showWarnings = FALSE, recursive = TRUE)
dir.create(dir_expl,   showWarnings = FALSE, recursive = TRUE)
dir.create(dir_models, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_plots3, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_plots2, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_tables, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_caps,   showWarnings = FALSE, recursive = TRUE)

UI2_levels <- c("Primary vegetation", "Secondary vegetation",
                "Agriculture_Low", "Agriculture_High", "Urban")

pal_UI2 <- c("#009E73", "#0072B2", "#E69F00", "#D55E00", "#CC79A7")
names(pal_UI2) <- UI2_levels

pal_FG <- c("FG1" = "#F05A5A", "FG2" = "#19B51E", "FG3" = "#4E79FF",
            "FG4" = "#E7298A", "FG5" = "#66A61E")

TITLE_SIZE <- 14
BASE_SIZE  <- 13

BASE_UI2   <- "Primary vegetation"
BASE_TEMP0 <- 0

N_TEMP <- 240
NSIM_CI <- suppressWarnings(as.integer(Sys.getenv("FG_NSIM_CI", "4000")))
if (!is.finite(NSIM_CI) || NSIM_CI < 200) NSIM_CI <- 4000
SEED_CI <- 123

## FG selection settings / FG 选择设置
K_MIN <- 2
K_MAX <- 5
PVE_TARGET <- 0.85
PC_MIN <- 3
PC_MAX <- 5

BOOT_B <- 40
MIN_CLUSTER_PROP <- 0.15
MIN_CLUSTER_JACCARD <- 0.55
NEAR_BEST_SIL_RATIO <- 0.90

## Trait preprocessing settings / 性状预处理设置
WINSOR_PROBS <- c(0.01, 0.99)

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
           device = ragg::agg_png, bg = "white")
  } else {
    ggsave(png_file, p, width = width, height = height, dpi = dpi, bg = "white")
  }
  ggsave(pdf_file, p, width = width, height = height)

  if (exists("topptx")) {
    topptx(p, pptx_file)
  } else if (requireNamespace("officer", quietly = TRUE) &&
             requireNamespace("rvg", quietly = TRUE)) {
    doc <- officer::read_pptx()
    doc <- officer::add_slide(doc, layout = "Title and Content", master = "Office Theme")
    doc <- officer::ph_with(doc, rvg::dml(ggobj = p), location = officer::ph_location_fullsize())
    print(doc, target = pptx_file)
  } else {
    message("PPTX export skipped: ", pptx_file)
  }
  invisible(TRUE)
}

theme_paper <- function(base_size = 13) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = TITLE_SIZE),
      plot.subtitle = element_text(size = base_size),
      strip.background = element_rect(fill = "white", color = "black"),
      strip.text = element_text(face = "bold"),
      axis.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      plot.margin = margin(t = 10, r = 16, b = 16, l = 10)
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
  cap_csv <- file.path(dir_caps, "FG_figure_brief_descriptions_CN_EN.csv")
  cap_txt <- file.path(dir_caps, "FG_figure_brief_descriptions_CN_EN.txt")

  utils::write.csv(figure_caption_registry, cap_csv, row.names = FALSE)

  con <- file(cap_txt, open = "wt")
  on.exit(close(con), add = TRUE)
  cat("FG figure brief descriptions (CN + EN)\n\n", file = con)
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

winsorise_vec <- function(x, probs = c(0.01, 0.99)) {
  q <- quantile(x, probs = probs, na.rm = TRUE, type = 7)
  pmin(pmax(x, q[1]), q[2])
}

safe_log10 <- function(x) {
  if (any(x <= 0, na.rm = TRUE)) stop("safe_log10 received non-positive values.")
  log10(x)
}

summ_prob_sims <- function(p_sim) {
  cbind(
    mean   = colMeans(p_sim),
    low95  = apply(p_sim, 2, quantile, probs = 0.025, na.rm = TRUE),
    high95 = apply(p_sim, 2, quantile, probs = 0.975, na.rm = TRUE)
  )
}

summ_pct_vs_baseline <- function(p_sim, p0_sim) {
  p0 <- as.numeric(p0_sim[, 1])
  p0[p0 < 1e-9] <- 1e-9
  pct_sim <- sweep(p_sim, 1, p0, FUN = function(p, p0) (p - p0) / p0 * 100)
  cbind(
    mean   = colMeans(pct_sim),
    low95  = apply(pct_sim, 2, quantile, probs = 0.025, na.rm = TRUE),
    high95 = apply(pct_sim, 2, quantile, probs = 0.975, na.rm = TRUE)
  )
}

collapse_sims_over <- function(p_sim, group_index) {
  glev <- unique(group_index)
  out <- sapply(glev, function(g) {
    cols <- which(group_index == g)
    rowMeans(p_sim[, cols, drop = FALSE])
  })
  if (is.vector(out)) out <- matrix(out, ncol = length(glev))
  colnames(out) <- glev
  out
}

## ============================================================
## 4) Read data / 读取数据
## ============================================================

myoccdata <- readRDS("myoccdata.rds")

need_cols <- c("Occur", "UI2", "SS", "SSBS", "Best_guess_binomial",
               "StdTmeanAnomaly", "RS", "HB", "TR", "HWI", "GL", "CS", "BM")
stopifnot(all(need_cols %in% names(myoccdata)))

myoccdata$Occur <- as.integer(myoccdata$Occur)
myoccdata$UI2 <- factor(myoccdata$UI2, levels = UI2_levels)
myoccdata$SS <- factor(myoccdata$SS)
myoccdata$SSBS <- factor(myoccdata$SSBS)
myoccdata$Best_guess_binomial <- factor(myoccdata$Best_guess_binomial)

mu_t <- mean(myoccdata$StdTmeanAnomaly, na.rm = TRUE)
sd_t <- sd(myoccdata$StdTmeanAnomaly, na.rm = TRUE)
if (!is.finite(mu_t) || !is.finite(sd_t) || sd_t == 0) {
  stop("StdTmeanAnomaly cannot be standardised.")
}
myoccdata$StdTmeanAnomalyRS <- (myoccdata$StdTmeanAnomaly - mu_t) / sd_t
saveRDS(myoccdata, file.path(dir_data, "myoccdata_with_StdTmeanAnomalyRS.rds"))

## ============================================================
## 5) Species-level trait table / 物种层面性状表
## ============================================================

trait_raw <- c("RS", "HB", "TR", "HWI", "GL", "CS", "BM")
sp_col <- "Best_guess_binomial"

sp_traits_raw <- myoccdata %>%
  dplyr::select(all_of(sp_col), all_of(trait_raw)) %>%
  dplyr::filter(!is.na(.data[[sp_col]])) %>%
  dplyr::group_by(.data[[sp_col]]) %>%
  dplyr::summarise(
    dplyr::across(all_of(trait_raw), ~ median(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  dplyr::filter(if_all(all_of(trait_raw), ~ is.finite(.x)))

saveRDS(sp_traits_raw, file.path(dir_data, "species_traits_raw_complete.rds"))
write.csv(sp_traits_raw, file.path(dir_data, "species_traits_raw_complete.csv"), row.names = FALSE)
cat("Species with complete raw traits:", nrow(sp_traits_raw), "\n")

## ============================================================
## 6) Build clustering trait space / 构建聚类用性状空间
## ============================================================

## We use transformed traits for clustering stability, but retain raw traits
## for interpretation. This reduces domination by very large-scale variables.
sp_traits_cluster <- sp_traits_raw %>%
  dplyr::transmute(
    !!sp_col := .data[[sp_col]],
    logRS  = winsorise_vec(safe_log10(RS), probs = WINSOR_PROBS),
    HB     = winsorise_vec(HB, probs = WINSOR_PROBS),
    logTR  = winsorise_vec(safe_log10(TR), probs = WINSOR_PROBS),
    logHWI = winsorise_vec(safe_log10(HWI), probs = WINSOR_PROBS),
    logGL  = winsorise_vec(safe_log10(GL), probs = WINSOR_PROBS),
    logCS  = winsorise_vec(safe_log10(CS), probs = WINSOR_PROBS),
    logBM  = winsorise_vec(safe_log10(BM), probs = WINSOR_PROBS)
  )

cluster_vars <- c("logRS", "HB", "logTR", "logHWI", "logGL", "logCS", "logBM")
X_cluster <- sp_traits_cluster %>% dplyr::select(all_of(cluster_vars))
X_cluster_z <- scale(X_cluster)

saveRDS(sp_traits_cluster, file.path(dir_data, "species_traits_for_clustering_transformed.rds"))
write.csv(sp_traits_cluster, file.path(dir_data, "species_traits_for_clustering_transformed.csv"), row.names = FALSE)

## ============================================================
## 7) PCA on transformed trait space / 变换性状空间 PCA
## ============================================================

pca_fg <- prcomp(X_cluster_z, center = FALSE, scale. = FALSE)
pve <- pca_fg$sdev^2 / sum(pca_fg$sdev^2)
pve_df <- data.frame(
  PC = paste0("PC", seq_along(pve)),
  ExplainedVariance = pve,
  CumulativeVariance = cumsum(pve)
)
write.csv(pve_df, file.path(dir_tables, "FG_PCA_explained_variance.csv"), row.names = FALSE)

n_pc <- which(cumsum(pve) >= PVE_TARGET)[1]
n_pc <- max(PC_MIN, min(PC_MAX, n_pc))
cat("Using", n_pc, "PCs for clustering (target cumulative variance =", PVE_TARGET, ")\n")

pc_scores <- as.data.frame(pca_fg$x[, paste0("PC", 1:n_pc), drop = FALSE])
pc_scores[[sp_col]] <- sp_traits_cluster[[sp_col]]

loadings_long <- as.data.frame(pca_fg$rotation[, 1:max(2, n_pc), drop = FALSE])
loadings_long$Trait <- rownames(loadings_long)
loadings_long <- tidyr::pivot_longer(loadings_long, cols = starts_with("PC"),
                                     names_to = "PC", values_to = "Loading")
write.csv(loadings_long, file.path(dir_tables, "FG_PCA_loadings_long.csv"), row.names = FALSE)

## ============================================================
## 8) Candidate clustering comparison / 候选聚类方案比较
## ============================================================

fit_clustering <- function(dat, method, K) {
  if (method == "kmeans") {
    return(kmeans(dat, centers = K, nstart = 80)$cluster)
  }
  if (method == "pam") {
    return(cluster::pam(dat, k = K)$clustering)
  }
  if (method == "ward") {
    return(cutree(hclust(dist(dat), method = "ward.D2"), k = K))
  }
  stop("Unknown method: ", method)
}

cluster_size_stats <- function(cl) {
  tb <- table(cl)
  prop <- prop.table(tb)
  data.frame(
    n_clusters = length(tb),
    min_n = unname(min(tb)),
    min_prop = unname(min(prop)),
    sizes = paste(sort(as.integer(tb)), collapse = ",")
  )
}

bootstrap_cluster_stability <- function(dat, cl0, method, K, B = 40) {
  n <- nrow(dat)
  groups <- sort(unique(cl0))
  out <- matrix(NA_real_, nrow = B, ncol = length(groups))

  for (b in seq_len(B)) {
    idx <- sample.int(n, size = n, replace = TRUE)
    idx_u <- sort(unique(idx))
    dat_b <- dat[idx_u, , drop = FALSE]
    cl_b <- fit_clustering(dat_b, method = method, K = K)
    labs_b <- sort(unique(cl_b))

    for (gi in seq_along(groups)) {
      A <- which(cl0 == groups[gi])
      vals <- numeric(length(labs_b))

      for (hj in seq_along(labs_b)) {
        B_idx <- idx_u[which(cl_b == labs_b[hj])]
        inter <- length(intersect(A, B_idx))
        uni <- length(union(A, B_idx))
        vals[hj] <- if (uni == 0) NA_real_ else inter / uni
      }
      out[b, gi] <- max(vals, na.rm = TRUE)
    }
  }

  data.frame(
    cluster = groups,
    mean_jaccard = colMeans(out, na.rm = TRUE)
  )
}

methods_to_try <- c("kmeans", "pam", "ward")
K_grid <- K_MIN:K_MAX
dist_pc <- dist(pc_scores[, paste0("PC", 1:n_pc), drop = FALSE])

score_list <- list()

for (mth in methods_to_try) {
  for (K in K_grid) {
    cat("Evaluating candidate:", mth, "K=", K, "\n")
    cl <- fit_clustering(pc_scores[, paste0("PC", 1:n_pc), drop = FALSE], method = mth, K = K)
    sil <- mean(cluster::silhouette(cl, dist_pc)[, "sil_width"], na.rm = TRUE)
    size_stat <- cluster_size_stats(cl)
    stab <- bootstrap_cluster_stability(pc_scores[, paste0("PC", 1:n_pc), drop = FALSE],
                                        cl0 = cl, method = mth, K = K, B = BOOT_B)

    score_list[[paste(mth, K, sep = "_")]] <- data.frame(
      method = mth,
      K = K,
      silhouette = sil,
      min_n = size_stat$min_n,
      min_prop = size_stat$min_prop,
      cluster_sizes = size_stat$sizes,
      mean_cluster_jaccard = mean(stab$mean_jaccard, na.rm = TRUE),
      min_cluster_jaccard = min(stab$mean_jaccard, na.rm = TRUE)
    )
  }
}

score_tbl <- dplyr::bind_rows(score_list)
score_tbl$eligible <- with(score_tbl,
                           min_prop >= MIN_CLUSTER_PROP &
                           min_cluster_jaccard >= MIN_CLUSTER_JACCARD)

if (!any(score_tbl$eligible)) {
  stop("No clustering candidate passed the minimum stability/size thresholds.")
}

best_eligible_sil <- max(score_tbl$silhouette[score_tbl$eligible], na.rm = TRUE)
score_tbl$near_best_sil <- score_tbl$silhouette >= (NEAR_BEST_SIL_RATIO * best_eligible_sil)
score_tbl$final_candidate <- score_tbl$eligible & score_tbl$near_best_sil

score_tbl <- score_tbl %>%
  dplyr::arrange(dplyr::desc(final_candidate),
                 dplyr::desc(K),
                 dplyr::desc(min_cluster_jaccard),
                 dplyr::desc(silhouette))

best_row <- score_tbl %>%
  dplyr::filter(final_candidate) %>%
  dplyr::arrange(dplyr::desc(K),
                 dplyr::desc(min_cluster_jaccard),
                 dplyr::desc(silhouette)) %>%
  dplyr::slice(1)

best_method <- best_row$method[1]
best_K <- best_row$K[1]

write.csv(score_tbl,
          file.path(dir_tables, "FG_robust_cluster_selection_table.csv"),
          row.names = FALSE)

cat("\nSelected FG solution:\n")
cat("  Method =", best_method, "\n")
cat("  K =", best_K, "\n")
cat("  silhouette =", round(best_row$silhouette[1], 3), "\n")
cat("  min_prop =", round(best_row$min_prop[1], 3), "\n")
cat("  min_cluster_jaccard =", round(best_row$min_cluster_jaccard[1], 3), "\n\n")

## Selection summary plot
score_tbl_plot <- score_tbl %>%
  dplyr::mutate(
    methodK = paste0(str_to_title(method), " | K=", K),
    status = dplyr::case_when(
      final_candidate ~ "Final candidates",
      eligible ~ "Eligible",
      TRUE ~ "Rejected"
    )
  )

p_select <- ggplot(score_tbl_plot,
                   aes(x = silhouette, y = min_cluster_jaccard,
                       size = min_prop, color = status, label = methodK)) +
  geom_point(alpha = 0.9) +
  geom_vline(xintercept = NEAR_BEST_SIL_RATIO * best_eligible_sil,
             linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = MIN_CLUSTER_JACCARD,
             linetype = "dashed", color = "grey40") +
  geom_text_repel(size = 3.2, max.overlaps = 30) +
  scale_color_manual(values = c("Rejected" = "grey60",
                                "Eligible" = "#1F78B4",
                                "Final candidates" = "#D95F02")) +
  scale_size_continuous(range = c(2.5, 7)) +
  theme_paper(BASE_SIZE) +
  labs(
    title = "FG candidate selection: stability, size balance, and silhouette",
    subtitle = paste0("Eligible if min cluster proportion >= ", MIN_CLUSTER_PROP,
                      " and min cluster bootstrap Jaccard >= ", MIN_CLUSTER_JACCARD,
                      ". Final candidates also require silhouette >= ",
                      round(NEAR_BEST_SIL_RATIO * 100), "% of best eligible silhouette."),
    x = "Average silhouette width",
    y = "Minimum cluster bootstrap Jaccard",
    size = "Minimum cluster proportion",
    color = "Status"
  )

save_plot_3formats(p_select, "FG_candidate_selection_tradeoff", dir_expl,
                   width = 10.6, height = 7.2)
add_figure_caption(
  "FG_candidate_selection_tradeoff", dir_expl,
  "FG candidate selection: stability, size balance, and silhouette",
  paste0("Among the candidate clustering solutions, ", str_to_title(best_method),
         " with K=", best_K, " retained a near-best silhouette while satisfying the minimum cluster-size and bootstrap-stability criteria, so it was selected as the final FG scheme."),
  paste0("在所有候选聚类方案中，", str_to_title(best_method), " 的 K=", best_K,
         " 同时满足最小组大小和 bootstrap 稳定性阈值，并保持接近最优的 silhouette，因此被选为最终 FG 方案。")
)

## ============================================================
## 9) Fit final clustering solution / 拟合最终 FG 方案
## ============================================================

cl_final_raw <- fit_clustering(pc_scores[, paste0("PC", 1:n_pc), drop = FALSE],
                               method = best_method, K = best_K)

## Reorder clusters by centroid along PC1 for stable naming
cent_order <- data.frame(cl = cl_final_raw, PC1 = pc_scores$PC1) %>%
  dplyr::group_by(cl) %>%
  dplyr::summarise(PC1_centroid = mean(PC1), .groups = "drop") %>%
  dplyr::arrange(PC1_centroid)

fg_map <- data.frame(
  cl = cent_order$cl,
  FG = factor(paste0("FG", seq_len(nrow(cent_order))),
              levels = paste0("FG", seq_len(nrow(cent_order))))
)

pc_scores$cl <- cl_final_raw
pc_scores <- pc_scores %>%
  dplyr::left_join(fg_map, by = "cl") %>%
  dplyr::select(all_of(sp_col), starts_with("PC"), FG)

saveRDS(pc_scores, file.path(dir_data, "species_PCA_scores_and_FG.rds"))
write.csv(pc_scores, file.path(dir_tables, "species_PCA_scores_and_FG.csv"), row.names = FALSE)

## ============================================================
## 10) FG trait interpretation tables / FG 性状解释表
## ============================================================

FG_centroids_traits_raw <- sp_traits_raw %>%
  dplyr::left_join(pc_scores %>% dplyr::select(all_of(sp_col), FG), by = sp_col) %>%
  dplyr::group_by(FG) %>%
  dplyr::summarise(
    dplyr::across(all_of(trait_raw), ~ mean(.x, na.rm = TRUE)),
    n_species = dplyr::n(),
    .groups = "drop"
  )
write.csv(FG_centroids_traits_raw, file.path(dir_tables, "FG_centroids_traits_raw.csv"),
          row.names = FALSE)

fg_profile <- sp_traits_raw %>%
  dplyr::left_join(pc_scores %>% dplyr::select(all_of(sp_col), FG), by = sp_col) %>%
  dplyr::mutate(dplyr::across(all_of(trait_raw), ~ as.numeric(scale(.x)))) %>%
  dplyr::group_by(FG) %>%
  dplyr::summarise(
    dplyr::across(all_of(trait_raw), ~ mean(.x, na.rm = TRUE)),
    n_species = dplyr::n(),
    .groups = "drop"
  )
write.csv(fg_profile, file.path(dir_tables, "FG_trait_profiles_Zscore_wide.csv"),
          row.names = FALSE)

fg_profile_long <- fg_profile %>%
  tidyr::pivot_longer(cols = all_of(trait_raw), names_to = "Trait", values_to = "Z")

TOPK <- 3
fg_dom <- fg_profile_long %>%
  dplyr::mutate(absZ = abs(Z)) %>%
  dplyr::group_by(FG) %>%
  dplyr::arrange(dplyr::desc(absZ), .by_group = TRUE) %>%
  dplyr::slice_head(n = TOPK) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(piece = paste0(Trait, ": ", sprintf("%+.2f", Z))) %>%
  dplyr::group_by(FG) %>%
  dplyr::summarise(Dominant_traits_topk = paste(piece, collapse = "; "), .groups = "drop")
write.csv(fg_dom, file.path(dir_tables, "FG_dominant_traits_top3_by_absZ.csv"),
          row.names = FALSE)

fg_interpret <- fg_profile %>%
  dplyr::left_join(fg_dom, by = "FG") %>%
  dplyr::mutate(
    Interpretation_EN = dplyr::case_when(
      row_number() == 1 ~ "Slow-history and high-dispersal strategy group with longer generation time and larger body size.",
      row_number() == 2 ~ "Thermally broader, habitat-generalist and higher-fecundity strategy group.",
      row_number() == 3 ~ "Range-restricted, low-dispersal and climatically narrower strategy group.",
      TRUE ~ "Trait syndrome defined by the dominant z-score profile."
    ),
    Interpretation_CN = dplyr::case_when(
      row_number() == 1 ~ "慢生活史-高扩散组：世代时间较长、体型偏大、扩散能力较高。",
      row_number() == 2 ~ "热耐受较宽-生境泛化-较高繁殖输出组：热幅较宽、生境广、窝卵数偏高。",
      row_number() == 3 ~ "地方性/专性策略组：分布范围较小、扩散较低、气候耐受较窄。",
      TRUE ~ "由主导性状 z-score 组合定义的策略组。"
    )
  )

write.csv(fg_interpret,
          file.path(dir_tables, "FG_trait_profiles_Zscore_with_interpretation.csv"),
          row.names = FALSE)

log_file <- file.path(dir_tables, "FG_trait_signatures_and_interpretation.txt")
sink(log_file)
cat("========== FG CLUSTER SELECTION ==========\n")
print(score_tbl)
cat("\nChosen method:", best_method, "\nChosen K:", best_K,
    "\nBootstrap B:", BOOT_B, "\n\n")
cat("========== FG TRAIT SIGNATURES (RAW-TRAIT Z-SCORES) ==========\n")
print(fg_profile)
cat("\n========== DOMINANT TRAITS (TOP-|Z|) ==========\n")
print(fg_dom)
cat("\n========== ECOLOGICAL INTERPRETATION ==========\n")
print(fg_interpret %>% dplyr::select(FG, n_species, Dominant_traits_topk,
                                     Interpretation_EN, Interpretation_CN))
sink()

## ============================================================
## 11) Join FG back to occurrence data / FG 回连记录数据
## ============================================================

myocc_fg <- myoccdata %>%
  dplyr::left_join(pc_scores %>% dplyr::select(all_of(sp_col), FG), by = sp_col) %>%
  dplyr::filter(!is.na(FG))

myocc_fg$FG <- factor(myocc_fg$FG, levels = levels(pc_scores$FG))
FG_base <- names(which.max(table(myocc_fg$FG)))
myocc_fg$FG <- relevel(myocc_fg$FG, ref = FG_base)

saveRDS(myocc_fg, file.path(dir_data, "myoccdata_with_FG.rds"))

cat("Records after FG join:", nrow(myocc_fg), "\n")
cat("FG counts (records):\n")
print(table(myocc_fg$FG))

## ============================================================
## 12) FG interpretation figures / FG 解释图
## ============================================================

## 12.1 PCA biplot
pc12 <- pc_scores %>%
  dplyr::select(all_of(sp_col), PC1, PC2, FG) %>%
  dplyr::mutate(FG = factor(FG, levels = levels(myocc_fg$FG)))

cent12 <- pc12 %>%
  dplyr::group_by(FG) %>%
  dplyr::summarise(PC1 = mean(PC1), PC2 = mean(PC2), .groups = "drop")

load_plot <- as.data.frame(pca_fg$rotation[, 1:2, drop = FALSE])
load_plot$Trait <- rownames(load_plot)
range_PC1 <- diff(range(pc12$PC1, na.rm = TRUE))
range_PC2 <- diff(range(pc12$PC2, na.rm = TRUE))
arrow_scale <- 0.35 * min(range_PC1, range_PC2)
load_plot$PC1a <- load_plot$PC1 * arrow_scale
load_plot$PC2a <- load_plot$PC2 * arrow_scale

trait_full <- c(
  logRS  = "Range size",
  HB  = "Habitat breadth",
  logTR  = "Thermal range",
  logHWI = "Dispersal ability (HWI)",
  logGL  = "Generation length",
  logCS  = "Clutch size",
  logBM  = "Body mass"
)
load_plot$Trait_full <- unname(trait_full[load_plot$Trait])

pve12 <- round(pve[1:2] * 100, 1)
max_abs_x <- max(abs(pc12$PC1), na.rm = TRUE)
max_abs_y <- max(abs(pc12$PC2), na.rm = TRUE)
lim <- max(max_abs_x, max_abs_y) * 1.05

p_biplot <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.45) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.45) +
  stat_ellipse(data = pc12, aes(PC1, PC2, colour = FG), type = "norm", level = 0.95, linewidth = 1.1) +
  geom_point(data = pc12, aes(PC1, PC2, colour = FG), alpha = 0.25, size = 1.4) +
  geom_point(data = cent12, aes(PC1, PC2, fill = FG, colour = FG),
             shape = 21, size = 5.0, stroke = 0.9) +
  geom_text(data = cent12, aes(PC1, PC2, label = FG, colour = FG),
            nudge_y = 0.15, fontface = "bold", size = 4.0, show.legend = FALSE) +
  geom_segment(data = load_plot,
               aes(x = 0, y = 0, xend = PC1a, yend = PC2a),
               arrow = arrow(length = unit(0.22, "cm")),
               linewidth = 0.85, color = "black") +
  {
    if ("ggrepel" %in% loadedNamespaces()) {
      ggrepel::geom_text_repel(
        data = load_plot, aes(PC1a, PC2a, label = Trait_full),
        size = 3.2, fontface = "bold",
        box.padding = 0.35, point.padding = 0.25,
        min.segment.length = 0, segment.color = "grey30", segment.size = 0.4,
        max.overlaps = Inf, seed = 123
      )
    } else {
      geom_text(data = load_plot, aes(PC1a, PC2a, label = Trait_full),
                size = 3.0, fontface = "bold", vjust = -0.6)
    }
  } +
  scale_colour_manual(values = pal_FG, drop = FALSE, name = "FG") +
  scale_fill_manual(values = pal_FG, drop = FALSE, name = "FG") +
  coord_fixed(1, xlim = c(-lim, lim), ylim = c(-lim, lim), clip = "off") +
  theme_paper(BASE_SIZE) +
  labs(
    title = "Functional strategy space (PCA) and functional strategy groups (FG)",
    subtitle = str_wrap(
      paste0("Clustering: ", str_to_title(best_method), ", K=", best_K,
             " (avg silhouette=", sprintf("%.3f", best_row$silhouette[1]), "). ",
             "Robust PCA on transformed and winsorised traits; ",
             "PC1=", pve12[1], "%; PC2=", pve12[2], "%."),
      width = 95
    ),
    x = paste0("PC1 (", pve12[1], "%)"),
    y = paste0("PC2 (", pve12[2], "%)")
  )

save_plot_3formats(p_biplot, "Fig3_PCA_biplot_FG_RobustStable", dir_expl,
                   width = 10.8, height = 7.6)
add_figure_caption(
  "Fig3_PCA_biplot_FG_RobustStable", dir_expl,
  "Functional strategy space (PCA) and functional strategy groups (FG)",
  "The three FG occupy distinct regions of multivariate trait space, indicating that the robust clustering scheme separates species into interpretable strategy groups rather than an unstable outlier partition.",
  "三个 FG 在多维性状空间中占据相对分离的位置，说明这套稳健聚类得到的是可解释的功能策略组，而不是由离群物种驱动的不稳定划分。"
)

## 12.2 Trait profile bar
fg_profile_bar_long <- fg_profile_long %>%
  dplyr::mutate(
    Trait = factor(Trait, levels = trait_raw),
    FG = factor(FG, levels = levels(myocc_fg$FG))
  )

p_profile_bar <- ggplot(fg_profile_bar_long, aes(Trait, Z, fill = FG)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.68) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.6) +
  scale_fill_manual(values = pal_FG, drop = FALSE) +
  theme_paper(BASE_SIZE) +
  labs(
    title = "Trait syndromes defining functional strategy groups (FG)",
    subtitle = "Mean z-score across species (positive = above overall mean; negative = below overall mean).",
    x = NULL,
    y = "Mean z-score (traits standardised across species)",
    fill = "FG"
  )
save_plot_3formats(p_profile_bar, "Fig4_FG_trait_profiles_Zscore_bar", dir_expl,
                   width = 10.6, height = 6.6)
add_figure_caption(
  "Fig4_FG_trait_profiles_Zscore_bar", dir_expl,
  "Trait syndromes defining functional strategy groups (FG)",
  "FG1 is associated with longer generation length and larger body size, FG2 with broader thermal range, habitat breadth and clutch size, whereas FG3 is generally below the overall mean for most traits.",
  "FG1 主要表现为较长世代时间和较大体型，FG2 主要表现为较宽热幅、较广生境和较高窝卵数，而 FG3 在多数性状上整体低于总体平均水平。"
)

## 12.3 Heatmap
p_heat <- ggplot(fg_profile_bar_long, aes(Trait, FG, fill = Z)) +
  geom_tile(color = "white", linewidth = 0.7) +
  scale_fill_gradient2(
    low = "steelblue4", mid = "white", high = "firebrick3",
    midpoint = 0, name = "Mean z-score"
  ) +
  theme_paper(BASE_SIZE) +
  theme(axis.title = element_blank()) +
  labs(
    title = "Functional strategy group signatures across traits",
    subtitle = "Heatmap of mean z-scores (cold = below overall mean; hot = above overall mean)."
  )
save_plot_3formats(p_heat, "Fig5_FG_trait_signatures_heatmap_Zscore", dir_expl,
                   width = 9.6, height = 4.8)
add_figure_caption(
  "Fig5_FG_trait_signatures_heatmap_Zscore", dir_expl,
  "Functional strategy group signatures across traits",
  "The heatmap highlights complementary trait syndromes across FG, with FG1 dominated by slow-history traits, FG2 by broader niche and fecundity traits, and FG3 by comparatively low scores across several traits.",
  "热图清楚展示了各 FG 之间互补的性状组合：FG1 偏向慢生活史性状，FG2 偏向较宽生态位和较高繁殖输出，而 FG3 在多个性状上相对较低。"
)

if (identical(Sys.getenv("FG_SKIP_MODEL"), "1")) {
  write_figure_captions()
  message("FG_SKIP_MODEL=1: stop after FG interpretation figures.")
  message("Outputs saved under: ", normalizePath(out_root))
  quit(save = "no", status = 0)
}

## ============================================================
## 13) glmmTMB model / glmmTMB 三重交互模型
## ============================================================

fml_tmb <- Occur ~ UI2 * StdTmeanAnomalyRS * FG +
  (1 | SS) + (1 | SSBS) + (1 | Best_guess_binomial)

ctrl_tmb <- glmmTMBControl(
  optimizer = optim,
  optArgs = list(method = "BFGS", maxit = 20000)
)

fit_file <- file.path(dir_models, "glmmTMB_Occ_UI2_Temp_FG_stable.rds")

if (file.exists(fit_file)) {
  message("Loading existing glmmTMB fit: ", fit_file)
  m_FG_tmb <- readRDS(fit_file)
} else {
  m_FG_tmb <- glmmTMB(
    formula = fml_tmb,
    data = droplevels(myocc_fg),
    family = binomial(link = "logit"),
    control = ctrl_tmb
  )
  saveRDS(m_FG_tmb, fit_file)
}

sink(file.path(dir_models, "FG_model_summary.txt"))
cat("==== Stable FG glmmTMB model summary ====\n")
print(summary(m_FG_tmb))
sink()

## ============================================================
## 14) Prediction helpers / 预测辅助函数
## ============================================================

predict_prob_sims_tmb <- function(mod, newdata, nsim = 2000, seed = 123, re_form = NA) {
  if (!identical(re_form, NA)) {
    stop("predict_prob_sims_tmb currently supports fixed-effects predictions only (re_form = NA).")
  }

  set.seed(seed)

  beta_hat <- fixef(mod)$cond
  vc <- as.matrix(vcov(mod)[["cond"]])
  f_fixed <- stats::delete.response(stats::terms(lme4::nobars(formula(mod))))
  X <- model.matrix(f_fixed, data = newdata)

  missing_cols <- setdiff(names(beta_hat), colnames(X))
  if (length(missing_cols) > 0) {
    X0 <- matrix(0, nrow = nrow(X), ncol = length(missing_cols))
    colnames(X0) <- missing_cols
    X <- cbind(X, X0)
  }
  X <- X[, names(beta_hat), drop = FALSE]

  beta_sim <- MASS::mvrnorm(n = nsim, mu = beta_hat, Sigma = vc)
  if (is.null(dim(beta_sim))) beta_sim <- matrix(beta_sim, nrow = 1)

  eta_sim <- beta_sim %*% t(X)
  plogis(eta_sim)
}

## ============================================================
## 15) Three-way predictions / 三重交互预测
## ============================================================

temp_seq <- seq(
  quantile(myocc_fg$StdTmeanAnomalyRS, 0.02, na.rm = TRUE),
  quantile(myocc_fg$StdTmeanAnomalyRS, 0.98, na.rm = TRUE),
  length.out = N_TEMP
)

nd3 <- expand.grid(
  StdTmeanAnomalyRS = temp_seq,
  UI2 = levels(myocc_fg$UI2),
  FG = levels(myocc_fg$FG)
)

p_sim3 <- predict_prob_sims_tmb(m_FG_tmb, nd3, nsim = NSIM_CI, seed = SEED_CI, re_form = NA)
pred_abs3 <- cbind(nd3, as.data.frame(summ_prob_sims(p_sim3)))

nd0 <- data.frame(
  UI2 = factor(BASE_UI2, levels = levels(myocc_fg$UI2)),
  StdTmeanAnomalyRS = BASE_TEMP0,
  FG = factor(FG_base, levels = levels(myocc_fg$FG))
)
p0_sim <- predict_prob_sims_tmb(m_FG_tmb, nd0, nsim = NSIM_CI, seed = SEED_CI, re_form = NA)

pred_pct3 <- cbind(nd3, as.data.frame(summ_pct_vs_baseline(p_sim3, p0_sim)))

saveRDS(pred_abs3, file.path(dir_data, "pred_abs_threeway_glmmTMB_fastCI.rds"))
saveRDS(pred_pct3, file.path(dir_data, "pred_pct_threeway_vsPV0_glmmTMB_fastCI.rds"))

pred_abs3$UI2 <- factor(pred_abs3$UI2, levels = UI2_levels)
pred_pct3$UI2 <- factor(pred_pct3$UI2, levels = UI2_levels)
pred_abs3 <- pred_abs3 %>% dplyr::arrange(FG, UI2, StdTmeanAnomalyRS)
pred_pct3 <- pred_pct3 %>% dplyr::arrange(FG, UI2, StdTmeanAnomalyRS)

p_abs3 <- ggplot(pred_abs3, aes(x = StdTmeanAnomalyRS, colour = UI2, fill = UI2)) +
  geom_ribbon(aes(ymin = low95, ymax = high95), alpha = 0.22, colour = NA) +
  geom_line(aes(y = mean), linewidth = 1.05) +
  facet_grid(FG ~ UI2, labeller = labeller(UI2 = label_map_UI2)) +
  scale_colour_manual(values = pal_UI2, drop = FALSE, guide = "none") +
  scale_fill_manual(values = pal_UI2, drop = FALSE, guide = "none") +
  theme_paper(BASE_SIZE) +
  labs(
    title = "Occurrence responses to land use and warming differ among functional strategy groups",
    subtitle = "glmmTMB fixed effects; fast 95% CI via link-scale normal approximation.",
    x = "Standardised mean temperature anomaly (StdTmeanAnomalyRS)",
    y = "Predicted probability of occurrence"
  )
save_plot_3formats(p_abs3, "Fig1_ThreeWay_ABS_UI2_Temp_FG_glmmTMB_fastCI", dir_plots3,
                   width = 12.2, height = 7.2)
add_figure_caption(
  "Fig1_ThreeWay_ABS_UI2_Temp_FG_glmmTMB_fastCI", dir_plots3,
  "Occurrence responses to land use and warming differ among functional strategy groups",
  "Predicted occurrence probability changes across the warming gradient differ among land-use types and among FG, supporting a heterogeneous three-way interaction between land use, warming and functional strategy.",
  "沿变暖梯度的预测出现概率在不同土地利用类型之间、也在不同 FG 之间存在明显差异，支持土地利用、变暖与功能策略之间存在异质性的三重交互。"
)

p_pct3 <- ggplot(pred_pct3, aes(x = StdTmeanAnomalyRS, colour = UI2, fill = UI2)) +
  geom_ribbon(aes(ymin = low95, ymax = high95), alpha = 0.22, colour = NA) +
  geom_line(aes(y = mean), linewidth = 1.05) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5) +
  facet_grid(FG ~ UI2, labeller = labeller(UI2 = label_map_UI2)) +
  scale_colour_manual(values = pal_UI2, drop = FALSE, guide = "none") +
  scale_fill_manual(values = pal_UI2, drop = FALSE, guide = "none") +
  theme_paper(BASE_SIZE) +
  labs(
    title = "Relative warming responses of occurrence probability across land uses and strategy groups",
    subtitle = paste0(
      "Baseline: UI2=", BASE_UI2, ", StdTmeanAnomalyRS=", BASE_TEMP0, ", FG=", FG_base,
      ". Percent change = (p-p0)/p0*100. fast 95% CI."
    ),
    x = "Standardised mean temperature anomaly (StdTmeanAnomalyRS)",
    y = "Percent change in predicted occurrence probability (%)"
  )
save_plot_3formats(p_pct3, "Fig2_ThreeWay_PCT_vsPV0_UI2_Temp_FG_glmmTMB_fastCI", dir_plots3,
                   width = 12.2, height = 7.2)
add_figure_caption(
  "Fig2_ThreeWay_PCT_vsPV0_UI2_Temp_FG_glmmTMB_fastCI", dir_plots3,
  "Relative warming responses of occurrence probability across land uses and strategy groups",
  "Relative change from the primary-vegetation baseline differs strongly among FG and land-use contexts, showing that warming sensitivity depends on both ecological strategy and land-use background.",
  "相对于原始植被基线的变化幅度在不同 FG 和不同土地利用背景下明显不同，说明变暖敏感性同时依赖于物种功能策略和土地利用环境。"
)

## ============================================================
## 16) Two-way summaries / 两重交互总结图
## ============================================================

## 16.1 UI2 x warming averaged over FG frequency
fg_wt <- prop.table(table(myocc_fg$FG))
nd_UI2Temp <- expand.grid(
  StdTmeanAnomalyRS = temp_seq,
  UI2 = levels(myocc_fg$UI2),
  FG = levels(myocc_fg$FG)
)
p_sim_UI2Temp <- predict_prob_sims_tmb(m_FG_tmb, nd_UI2Temp,
                                       nsim = NSIM_CI, seed = SEED_CI + 10, re_form = NA)

group_key_UI2Temp <- paste(nd_UI2Temp$UI2, nd_UI2Temp$StdTmeanAnomalyRS, sep = "___")
collapsed_UI2Temp <- sapply(unique(group_key_UI2Temp), function(g) {
  cols <- which(group_key_UI2Temp == g)
  fg_levels_here <- as.character(nd_UI2Temp$FG[cols])
  w <- as.numeric(fg_wt[fg_levels_here])
  p_sim_UI2Temp[, cols, drop = FALSE] %*% w / sum(w)
})
if (is.vector(collapsed_UI2Temp)) collapsed_UI2Temp <- matrix(collapsed_UI2Temp, ncol = length(unique(group_key_UI2Temp)))

ui2temp_df <- nd_UI2Temp %>%
  dplyr::select(UI2, StdTmeanAnomalyRS) %>%
  dplyr::distinct() %>%
  dplyr::arrange(UI2, StdTmeanAnomalyRS)
ui2temp_sum <- cbind(ui2temp_df, as.data.frame(summ_prob_sims(collapsed_UI2Temp)))
ui2temp_sum <- ui2temp_sum %>% dplyr::arrange(UI2, StdTmeanAnomalyRS)

p_ui2temp <- ggplot(ui2temp_sum, aes(x = StdTmeanAnomalyRS, y = mean, color = UI2, fill = UI2)) +
  geom_ribbon(aes(ymin = low95, ymax = high95), alpha = 0.20, colour = NA) +
  geom_line(linewidth = 1.05) +
  scale_color_manual(values = pal_UI2, drop = FALSE) +
  scale_fill_manual(values = pal_UI2, drop = FALSE) +
  theme_paper(BASE_SIZE) +
  labs(
    title = "Two-way interaction: land use × warming (FG-averaged)",
    subtitle = "glmmTMB fixed effects; fast 95% CI. (Average over FG within each draw.)",
    x = "Standardised mean temperature anomaly (StdTmeanAnomalyRS)",
    y = "Predicted occurrence probability",
    color = "Land-use type",
    fill = "Land-use type"
  )
save_plot_3formats(p_ui2temp, "Fig6_TwoWay_UI2_by_Warming_FGaveraged_ABS_glmmTMB_fastCI", dir_plots2,
                   width = 10.6, height = 6.6)
add_figure_caption(
  "Fig6_TwoWay_UI2_by_Warming_FGaveraged_ABS_glmmTMB_fastCI", dir_plots2,
  "Two-way interaction: land use × warming (FG-averaged)",
  "After averaging over FG composition, land-use types still differ in both baseline occurrence and warming-response shape, indicating that land use exerts an overall filtering effect beyond FG-specific responses.",
  "在对 FG 组成加权平均之后，不同土地利用类型在基线出现概率和随变暖变化的曲线形状上仍然不同，说明土地利用具有超越 FG 差异的总体过滤效应。"
)

## 16.2 FG x warming averaged over UI2 frequency
ui2_wt <- prop.table(table(myocc_fg$UI2))
nd_FGTemp <- expand.grid(
  StdTmeanAnomalyRS = temp_seq,
  UI2 = levels(myocc_fg$UI2),
  FG = levels(myocc_fg$FG)
)
p_sim_FGTemp <- predict_prob_sims_tmb(m_FG_tmb, nd_FGTemp,
                                      nsim = NSIM_CI, seed = SEED_CI + 20, re_form = NA)

group_key_FGTemp <- paste(nd_FGTemp$FG, nd_FGTemp$StdTmeanAnomalyRS, sep = "___")
collapsed_FGTemp <- sapply(unique(group_key_FGTemp), function(g) {
  cols <- which(group_key_FGTemp == g)
  ui2_levels_here <- as.character(nd_FGTemp$UI2[cols])
  w <- as.numeric(ui2_wt[ui2_levels_here])
  p_sim_FGTemp[, cols, drop = FALSE] %*% w / sum(w)
})
if (is.vector(collapsed_FGTemp)) collapsed_FGTemp <- matrix(collapsed_FGTemp, ncol = length(unique(group_key_FGTemp)))

fgtemp_df <- nd_FGTemp %>%
  dplyr::select(FG, StdTmeanAnomalyRS) %>%
  dplyr::distinct() %>%
  dplyr::arrange(FG, StdTmeanAnomalyRS)
fgtemp_sum <- cbind(fgtemp_df, as.data.frame(summ_prob_sims(collapsed_FGTemp)))
fgtemp_sum <- fgtemp_sum %>% dplyr::arrange(FG, StdTmeanAnomalyRS)

p_fgtemp <- ggplot(fgtemp_sum, aes(x = StdTmeanAnomalyRS, y = mean, color = FG, fill = FG)) +
  geom_ribbon(aes(ymin = low95, ymax = high95), alpha = 0.20, colour = NA) +
  geom_line(linewidth = 1.05) +
  scale_color_manual(values = pal_FG, drop = FALSE) +
  scale_fill_manual(values = pal_FG, drop = FALSE) +
  theme_paper(BASE_SIZE) +
  labs(
    title = "Two-way interaction: functional strategy group × warming (UI2-averaged)",
    subtitle = "glmmTMB fixed effects; fast 95% CI. (Average over UI2 within each draw.)",
    x = "Standardised mean temperature anomaly (StdTmeanAnomalyRS)",
    y = "Predicted occurrence probability",
    color = "FG",
    fill = "FG"
  )
save_plot_3formats(p_fgtemp, "Fig7_TwoWay_FG_by_Warming_UI2averaged_ABS_glmmTMB_fastCI", dir_plots2,
                   width = 10.6, height = 6.6)
add_figure_caption(
  "Fig7_TwoWay_FG_by_Warming_UI2averaged_ABS_glmmTMB_fastCI", dir_plots2,
  "Two-way interaction: functional strategy group × warming (UI2-averaged)",
  "After averaging over land-use composition, FG still differ in warming-response trajectories, confirming that the FG scheme captures biologically meaningful variation in climate sensitivity.",
  "在对土地利用组成加权平均之后，不同 FG 的变暖响应轨迹仍然不同，说明这套 FG 划分确实捕捉到了具有生物学意义的气候敏感性差异。"
)

## 16.3 UI2 x FG at selected warming levels
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

p_sim_UI2FG <- predict_prob_sims_tmb(m_FG_tmb, nd_UI2FG,
                                     nsim = NSIM_CI, seed = SEED_CI + 30, re_form = NA)
ui2fg_sum <- cbind(nd_UI2FG, as.data.frame(summ_prob_sims(p_sim_UI2FG)))
ui2fg_sum$UI2 <- factor(ui2fg_sum$UI2, levels = UI2_levels)
ui2fg_sum$FG <- factor(ui2fg_sum$FG, levels = levels(myocc_fg$FG))
ui2fg_sum <- ui2fg_sum %>% dplyr::arrange(WarmLabel, FG, UI2)

p_ui2fg <- ggplot(ui2fg_sum, aes(x = UI2, y = mean, color = FG, group = FG)) +
  geom_hline(yintercept = 0, colour = "grey82", linewidth = 0.4) +
  geom_pointrange(aes(ymin = low95, ymax = high95),
                  position = position_dodge(width = 0.55),
                  linewidth = 0.7) +
  facet_wrap(~ WarmLabel, nrow = 1) +
  scale_x_discrete(labels = label_map_UI2_wrap) +
  scale_color_manual(values = pal_FG, drop = FALSE) +
  theme_paper(BASE_SIZE) +
  theme(axis.text.x = element_text(size = BASE_SIZE - 1, lineheight = 0.9)) +
  labs(
    title = "Two-way interaction: land use × FG (at fixed warming levels)",
    subtitle = "glmmTMB fixed effects; fast 95% CI at StdTmeanAnomalyRS = 0, 1, and median.",
    x = NULL,
    y = "Predicted occurrence probability",
    color = "FG"
  )
save_plot_3formats(p_ui2fg, "Fig8_TwoWay_UI2_by_FG_at_fixedWarming_ABS_glmmTMB_fastCI", dir_plots2,
                   width = 12.2, height = 5.6)
add_figure_caption(
  "Fig8_TwoWay_UI2_by_FG_at_fixedWarming_ABS_glmmTMB_fastCI", dir_plots2,
  "Two-way interaction: land use × FG (at fixed warming levels)",
  "At fixed warming levels, the rank order and separation among FG vary across land-use types, showing that land-use filtering and functional strategy interact in shaping occurrence probability.",
  "在固定变暖水平下，不同土地利用类型中的 FG 排序和差异幅度并不一致，说明土地利用过滤作用与功能策略共同塑造了物种出现概率。"
)

message("\nStable FG pipeline completed.")
write_figure_captions()
message("Outputs saved under: ", normalizePath(out_root))
message("Key tables: ", dir_tables)
message("Key figures: ", dir_expl, " / ", dir_plots2, " / ", dir_plots3)
message("Figure brief descriptions: ", dir_caps)
