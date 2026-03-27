###############################################################
## Land-use x warming x multiple traits (joint screen)
## Occurrence (brms) -- Bayesian joint multitrait screening
##
## Purpose
## 1) Provide a fully Bayesian counterpart to the earlier glmmTMB joint-trait screen
## 2) Keep the scientific scope as a two-way joint screen:
##      UI2 x climate, UI2 x trait, climate x trait
##    while avoiding Bayesian stepwise-AIC model selection
## 3) Export posterior diagnostics, coefficient summaries, and the same
##    main plot families used in the earlier multitrait script
##
## Main interpretation
## - This script is a robustness / conditional-effect screen.
## - It is NOT the headline three-way test for RQ3 on its own.
###############################################################

suppressPackageStartupMessages({
  library(brms)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
  library(Hmisc)
  library(ggcorrplot)
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
## 1) Global settings
## ============================================================

out_root <- "brms_occurrence_alltraits_joint2way_medianControls"
dir.create(out_root, showWarnings = FALSE, recursive = TRUE)

dir_corr   <- file.path(out_root, "00_TraitCorrelation")
dir_plots  <- file.path(out_root, "01_Plots_JointModel")
dir_models <- file.path(out_root, "02_Models")
dir_tables <- file.path(out_root, "03_Tables")
dir_caps   <- file.path(out_root, "04_Figure_Captions")
dir_diag   <- file.path(out_root, "98_Diagnostics")
for (d in c(dir_corr, dir_plots, dir_models, dir_tables, dir_caps, dir_diag)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

UI2_levels <- c(
  "Primary vegetation", "Secondary vegetation",
  "Agriculture_Low", "Agriculture_High", "Urban"
)

pal_UI2 <- c("#009E73", "#0072B2", "#E69F00", "#D55E00", "#CC79A7")
names(pal_UI2) <- UI2_levels

TITLE_SIZE <- 14
BASE_SIZE <- 13

BASE_UI2 <- "Primary vegetation"
BASE_TEMP0 <- 0
TEMP_SCENARIOS <- c(0, 1)

trait_vars_cont <- c(
  "RS.rs", "HB.rs", "TR.rs", "HWI.rs", "GL.rs", "CS.rs",
  "Tmin_position", "Tmax_position"
)

trait_labels <- c(
  "RS.rs" = "Range size",
  "HB.rs" = "Habitat breadth",
  "TR.rs" = "Thermal breadth",
  "HWI.rs" = "Dispersal ability (HWI)",
  "GL.rs" = "Generation length",
  "CS.rs" = "Clutch size",
  "Tmin_position" = "Cold-edge thermal position",
  "Tmax_position" = "Warm-edge thermal position"
)

BRMS_CHAINS <- suppressWarnings(as.integer(Sys.getenv("MT_BRMS_CHAINS", "4")))
BRMS_ITER <- suppressWarnings(as.integer(Sys.getenv("MT_BRMS_ITER", "6000")))
BRMS_WARMUP <- suppressWarnings(as.integer(Sys.getenv("MT_BRMS_WARMUP", "3000")))
BRMS_THREADS <- suppressWarnings(as.integer(Sys.getenv("MT_BRMS_THREADS", "2")))
BRMS_CORES <- suppressWarnings(as.integer(Sys.getenv("MT_BRMS_CORES", as.character(min(BRMS_CHAINS, parallel::detectCores())))))
BRMS_ADAPT <- suppressWarnings(as.numeric(Sys.getenv("MT_BRMS_ADAPT_DELTA", "0.99")))
BRMS_TREEDEPTH <- suppressWarnings(as.integer(Sys.getenv("MT_BRMS_MAX_TREEDEPTH", "15")))
BRMS_BACKEND <- Sys.getenv("MT_BRMS_BACKEND", unset = "cmdstanr")
BRMS_GRAINSIZE <- suppressWarnings(as.integer(Sys.getenv("MT_BRMS_GRAINSIZE", "1000")))
BRMS_FORCE_REFIT <- identical(tolower(Sys.getenv("MT_BRMS_FORCE_REFIT", "0")), "1")
BRMS_SEED <- 123
N_DRAWS_PRED <- suppressWarnings(as.integer(Sys.getenv("MT_BRMS_NDRAWS_PRED", "1600")))
SMOKE_N <- suppressWarnings(as.integer(Sys.getenv("MT_BRMS_SMOKE_N", "0")))

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
## 2) Helpers
## ============================================================

save_plot_3formats <- function(p, filename_noext, outdir,
                               width = 12.0, height = 6.6, dpi = 320) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  png_file  <- file.path(outdir, paste0(filename_noext, ".png"))
  pdf_file  <- file.path(outdir, paste0(filename_noext, ".pdf"))
  pptx_file <- file.path(outdir, paste0(filename_noext, ".pptx"))

  if (requireNamespace("ragg", quietly = TRUE)) {
    ggsave(png_file, p, width = width, height = height, dpi = dpi,
           device = ragg::agg_png, bg = "white", limitsize = FALSE)
  } else {
    ggsave(png_file, p, width = width, height = height, dpi = dpi,
           bg = "white", limitsize = FALSE)
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
  cap_csv <- file.path(dir_caps, "multitrait_brms_figure_descriptions_CN_EN.csv")
  cap_txt <- file.path(dir_caps, "multitrait_brms_figure_descriptions_CN_EN.txt")
  utils::write.csv(figure_caption_registry, cap_csv, row.names = FALSE)

  con <- file(cap_txt, open = "wt")
  on.exit(close(con), add = TRUE)
  cat("Multitrait brms figure descriptions (CN + EN)\n\n", file = con)
  for (i in seq_len(nrow(figure_caption_registry))) {
    row <- figure_caption_registry[i, ]
    cat("[", row$filename_noext, "] ", row$title, "\n", sep = "", file = con)
    cat("Path: ", row$figure_dir, "\n", sep = "", file = con)
    cat("EN: ", row$description_en, "\n", sep = "", file = con)
    cat("CN: ", row$description_cn, "\n\n", sep = "", file = con)
  }
}

label_map_UI2 <- c(
  "Primary vegetation" = "Primary vegetation",
  "Secondary vegetation" = "Secondary vegetation",
  "Agriculture_Low" = "Agriculture (low)",
  "Agriculture_High" = "Agriculture (high)",
  "Urban" = "Urban"
)

label_map_UI2_wrap <- c(
  "Primary vegetation" = "Primary\nvegetation",
  "Secondary vegetation" = "Secondary\nvegetation",
  "Agriculture_Low" = "Agriculture\n(low)",
  "Agriculture_High" = "Agriculture\n(high)",
  "Urban" = "Urban"
)

summ_draws <- function(draws_vec) {
  c(
    mean = mean(draws_vec),
    low95 = unname(quantile(draws_vec, 0.025)),
    high95 = unname(quantile(draws_vec, 0.975))
  )
}

percent_change_draws <- function(ep, ep0) {
  ep0v <- as.vector(ep0)
  ep0v[ep0v < 1e-9] <- 1e-9
  sweep(ep, 1, ep0v, FUN = function(p, p0) (p - p0) / p0 * 100)
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

extract_max_treedepth_hits <- function(fit_obj, treedepth) {
  if (inherits(fit_obj, "CmdStanMCMC")) {
    diag_arr <- fit_obj$sampler_diagnostics()
    if (length(dim(diag_arr)) == 3 && "treedepth__" %in% dimnames(diag_arr)[[3]]) {
      return(sum(diag_arr[, , "treedepth__"] >= treedepth, na.rm = TRUE))
    }
  }
  if (inherits(fit_obj, "stanfit")) {
    sampler_list <- rstan::get_sampler_params(fit_obj, inc_warmup = FALSE)
    return(sum(vapply(sampler_list, function(x) sum(x[, "treedepth__"] >= treedepth, na.rm = TRUE), numeric(1))))
  }
  NA_real_
}

pretty_coef_label <- function(x) {
  x <- stringr::str_replace_all(x, c(
    "StdTmeanAnomalyRS" = "Warming",
    "UI2Secondaryvegetation" = "Secondary vegetation",
    "UI2Agriculture_Low" = "Agriculture (low)",
    "UI2Agriculture_High" = "Agriculture (high)",
    "UI2Urban" = "Urban",
    "Tmin_position" = "Cold-edge position",
    "Tmax_position" = "Warm-edge position",
    "RS.rs" = "Range size",
    "HB.rs" = "Habitat breadth",
    "TR.rs" = "Thermal breadth",
    "HWI.rs" = "Dispersal ability",
    "GL.rs" = "Generation length",
    "CS.rs" = "Clutch size",
    ":" = " x "
  ))
  stringr::str_wrap(x, width = 38)
}

## ============================================================
## 3) Data
## ============================================================

myoccdata <- readRDS("myoccdata.rds")
need_cols <- c(
  "Occur", "UI2", "SS", "SSBS", "Best_guess_binomial",
  "StdTmeanAnomalyRS", trait_vars_cont
)
miss_cols <- setdiff(need_cols, names(myoccdata))
if (length(miss_cols) > 0) stop("Missing columns: ", paste(miss_cols, collapse = ", "))

myoccdata <- myoccdata %>%
  dplyr::select(all_of(need_cols)) %>%
  dplyr::mutate(
    Occur = as.numeric(Occur),
    UI2 = factor(UI2, levels = UI2_levels),
    SS = factor(SS),
    SSBS = factor(SSBS),
    Best_guess_binomial = factor(Best_guess_binomial)
  ) %>%
  tidyr::drop_na()

if (is.finite(SMOKE_N) && SMOKE_N > 0 && SMOKE_N < nrow(myoccdata)) {
  set.seed(BRMS_SEED)
  myoccdata <- myoccdata %>%
    dplyr::group_by(UI2) %>%
    dplyr::slice_sample(n = max(1, min(dplyr::n(), ceiling(SMOKE_N / length(UI2_levels))))) %>%
    dplyr::ungroup()
}

trait_med <- lapply(trait_vars_cont, function(v) median(myoccdata[[v]], na.rm = TRUE))
names(trait_med) <- trait_vars_cont
temp_median <- median(myoccdata$StdTmeanAnomalyRS, na.rm = TRUE)

support_df <- myoccdata %>%
  dplyr::count(UI2, name = "n_records") %>%
  dplyr::mutate(label = format(n_records, big.mark = ",", scientific = FALSE))
utils::write.csv(support_df, file.path(dir_tables, "multitrait_brms_UI2_record_support.csv"), row.names = FALSE)

run_manifest <- data.frame(
  metric = c(
    "backend", "chains", "iter", "warmup", "cores", "threads_per_chain",
    "adapt_delta", "max_treedepth", "prediction_draws",
    "smoke_n", "rows_used", "species_used", "seed"
  ),
  value = c(
    BRMS_BACKEND, BRMS_CHAINS, BRMS_ITER, BRMS_WARMUP, BRMS_CORES, BRMS_THREADS,
    BRMS_ADAPT, BRMS_TREEDEPTH, N_DRAWS_PRED,
    SMOKE_N, nrow(myoccdata), dplyr::n_distinct(myoccdata$Best_guess_binomial), BRMS_SEED
  )
)
utils::write.csv(run_manifest, file.path(dir_tables, "multitrait_brms_run_manifest.csv"), row.names = FALSE)

## ============================================================
## 4) Trait correlation
## ============================================================

trait_mat <- myoccdata %>%
  dplyr::select(all_of(trait_vars_cont)) %>%
  as.matrix()
rc <- Hmisc::rcorr(trait_mat, type = "pearson")
cor_mat <- rc$r
p_mat <- rc$P
utils::write.csv(cor_mat, file.path(dir_corr, "trait_correlation_r.csv"))
utils::write.csv(p_mat, file.path(dir_corr, "trait_correlation_p.csv"))

p_corr <- ggcorrplot::ggcorrplot(
  cor_mat,
  method = "square",
  type = "lower",
  lab = TRUE,
  p.mat = p_mat,
  sig.level = 0.05,
  insig = "blank"
) +
  theme_paper(BASE_SIZE) +
  labs(
    title = "Trait correlation structure used in the joint Bayesian multitrait screen",
    subtitle = "Blank cells indicate non-significant correlations at p >= 0.05."
  )
save_plot_3formats(p_corr, "Fig0_Trait_correlation_joint_multitrait_brms", dir_corr,
                   width = 10.0, height = 7.6)
add_figure_caption(
  "Fig0_Trait_correlation_joint_multitrait_brms", dir_corr,
  "Trait correlation structure used in the joint Bayesian multitrait screen",
  "The joint Bayesian multitrait model includes correlated predictors, so the correlation matrix provides important context for interpreting conditional posterior effects.",
  "联合 Bayesian 多性状模型同时包含多个相关预测变量，因此相关矩阵是解释条件后验效应的重要背景信息。"
)

## ============================================================
## 5) Joint brms model
## ============================================================

traits_term <- paste(trait_vars_cont, collapse = " + ")
joint_formula <- bf(
  as.formula(
    paste0(
      "Occur ~ UI2 * StdTmeanAnomalyRS",
      " + UI2 * (", traits_term, ")",
      " + StdTmeanAnomalyRS * (", traits_term, ")",
      " + (1 || SS) + (1 || SSBS) + (1 || Best_guess_binomial)"
    )
  ),
  family = bernoulli(link = "logit"),
  decomp = "QR"
)

pri <- c(
  prior(normal(0, 0.5), class = "b"),
  prior(normal(0, 1.5), class = "Intercept"),
  prior(exponential(2), class = "sd")
)

fit_file <- file.path(dir_models, "brms_occurrence_alltraits_joint2way_cmdstanr.rds")
if (file.exists(fit_file) && !BRMS_FORCE_REFIT) {
  message("Loading existing fit: ", fit_file)
  m_joint <- readRDS(fit_file)
} else {
  brm_args <- list(
    formula = joint_formula,
    data = droplevels(myoccdata),
    backend = BRMS_BACKEND,
    prior = pri,
    chains = BRMS_CHAINS,
    iter = BRMS_ITER,
    warmup = BRMS_WARMUP,
    cores = BRMS_CORES,
    seed = BRMS_SEED,
    control = list(adapt_delta = BRMS_ADAPT, max_treedepth = BRMS_TREEDEPTH),
    save_pars = save_pars(all = FALSE, group = FALSE),
    sparse = TRUE,
    normalize = FALSE,
    init = 0,
    refresh = 100
  )
  if (BRMS_THREADS > 1) {
    brm_args$threads <- threading(BRMS_THREADS, grainsize = BRMS_GRAINSIZE)
  }
  m_joint <- do.call(brm, brm_args)
  saveRDS(m_joint, fit_file)
}

## ============================================================
## 6) Diagnostics and coefficients
## ============================================================

fit_summary <- summary(m_joint)
capture.output(fit_summary, file = file.path(dir_models, "joint_multitrait_brms_model_summary.txt"))

fixef_tab <- as.data.frame(fixef(m_joint, probs = c(0.025, 0.975)))
fixef_tab$term <- rownames(fixef_tab)
rownames(fixef_tab) <- NULL
utils::write.csv(fixef_tab, file.path(dir_tables, "joint_multitrait_brms_fixed_effects_summary.csv"), row.names = FALSE)

diag_summary <- data.frame(
  metric = c("chains", "iter", "warmup", "threads_per_chain", "ndraws_total",
             "max_Rhat_fixed", "min_Bulk_ESS_fixed", "min_Tail_ESS_fixed",
             "divergences", "max_treedepth_hits"),
  value = c(
    BRMS_CHAINS,
    BRMS_ITER,
    BRMS_WARMUP,
    BRMS_THREADS,
    posterior::ndraws(m_joint),
    max(fit_summary$fixed[, "Rhat"], na.rm = TRUE),
    min(fit_summary$fixed[, "Bulk_ESS"], na.rm = TRUE),
    min(fit_summary$fixed[, "Tail_ESS"], na.rm = TRUE),
    extract_divergences(m_joint$fit),
    extract_max_treedepth_hits(m_joint$fit, BRMS_TREEDEPTH)
  )
)
utils::write.csv(diag_summary, file.path(dir_tables, "joint_multitrait_brms_diagnostics_summary.csv"), row.names = FALSE)

coef_focus <- fixef_tab %>%
  dplyr::filter(
    term == "StdTmeanAnomalyRS" |
      stringr::str_detect(term, "StdTmeanAnomalyRS") |
      stringr::str_detect(term, "UI2")
  ) %>%
  dplyr::mutate(label = pretty_coef_label(term)) %>%
  dplyr::arrange(Estimate)

p_coef <- ggplot(coef_focus, aes(x = Estimate, y = reorder(label, Estimate))) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5, colour = "grey55") +
  geom_linerange(aes(xmin = Q2.5, xmax = Q97.5), linewidth = 0.9, colour = "#1F3B70") +
  geom_point(size = 2.6, colour = "#1F3B70") +
  theme_paper(BASE_SIZE) +
  labs(
    title = "Posterior summaries for warming-related and land-use-related coefficients",
    subtitle = "Points are posterior means; horizontal intervals are 95% credible intervals on the logit scale.",
    x = "Posterior estimate (logit scale)",
    y = NULL
  )
save_plot_3formats(p_coef, "Fig1_Joint_multitrait_brms_coefficient_forest", dir_diag,
                   width = 10.6, height = 7.4)
add_figure_caption(
  "Fig1_Joint_multitrait_brms_coefficient_forest", dir_diag,
  "Posterior summaries for warming-related and land-use-related coefficients",
  "The forest plot summarises conditional posterior estimates in the joint multitrait model and helps identify which warming-related coefficients remain important after controlling for the other traits.",
  "森林图总结了联合多性状模型中的条件后验估计，有助于识别在控制其他性状后仍然重要的变暖相关系数。"
)

## ============================================================
## 7) Prediction builders
## ============================================================

posterior_summary_df <- function(ep_mat, newdata) {
  sm <- apply(ep_mat, 2, summ_draws)
  cbind(newdata, as.data.frame(t(sm)))
}

plot_UI2_by_temp_abs_pct <- function(mod, data, ndraws_pred) {
  temp_seq <- seq(
    quantile(data$StdTmeanAnomalyRS, 0.02, na.rm = TRUE),
    quantile(data$StdTmeanAnomalyRS, 0.98, na.rm = TRUE),
    length.out = 260
  )

  nd <- expand.grid(
    UI2 = levels(data$UI2),
    StdTmeanAnomalyRS = temp_seq
  )
  for (v in trait_vars_cont) nd[[v]] <- trait_med[[v]]

  nd0 <- data.frame(UI2 = factor(BASE_UI2, levels = levels(data$UI2)),
                    StdTmeanAnomalyRS = BASE_TEMP0)
  for (v in trait_vars_cont) nd0[[v]] <- trait_med[[v]]

  ep <- posterior_epred(mod, newdata = nd, re_formula = NA, ndraws = ndraws_pred)
  ep0 <- posterior_epred(mod, newdata = nd0, re_formula = NA, ndraws = ndraws_pred)

  df_abs <- posterior_summary_df(ep, nd)
  df_pct <- posterior_summary_df(percent_change_draws(ep, ep0), nd)

  p_abs <- ggplot(df_abs, aes(x = StdTmeanAnomalyRS, y = mean, colour = UI2, fill = UI2)) +
    geom_ribbon(aes(ymin = low95, ymax = high95), alpha = 0.20, colour = NA) +
    geom_line(linewidth = 1.05) +
    scale_colour_manual(values = pal_UI2, drop = FALSE) +
    scale_fill_manual(values = pal_UI2, drop = FALSE) +
    theme_paper(BASE_SIZE) +
    labs(
      title = "Joint multitrait Bayesian screen: land use × warming",
      subtitle = "All continuous traits are fixed at their sample medians.",
      x = "Standardised mean temperature anomaly (StdTmeanAnomalyRS)",
      y = "Predicted probability of occurrence",
      colour = "Land-use type",
      fill = "Land-use type"
    )

  p_pct <- ggplot(df_pct, aes(x = StdTmeanAnomalyRS, y = mean, colour = UI2, fill = UI2)) +
    geom_ribbon(aes(ymin = low95, ymax = high95), alpha = 0.20, colour = NA) +
    geom_line(linewidth = 1.05) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5) +
    scale_colour_manual(values = pal_UI2, drop = FALSE) +
    scale_fill_manual(values = pal_UI2, drop = FALSE) +
    theme_paper(BASE_SIZE) +
    labs(
      title = "Joint multitrait Bayesian screen: relative land use × warming effect",
      subtitle = "Percent change is relative to primary vegetation at warming = 0 with all traits at their medians.",
      x = "Standardised mean temperature anomaly (StdTmeanAnomalyRS)",
      y = "Percent change in predicted occurrence probability (%)",
      colour = "Land-use type",
      fill = "Land-use type"
    )

  list(abs_plot = p_abs, pct_plot = p_pct)
}

plot_UI2_by_trait_fixedTemp_abs_pct <- function(mod, data, trait_var, ndraws_pred) {
  temp_values <- sort(unique(c(TEMP_SCENARIOS, temp_median)))
  temp_labels <- c("0", "1", "median")
  names(temp_labels) <- as.character(temp_values)

  xseq <- seq(
    quantile(data[[trait_var]], 0.02, na.rm = TRUE),
    quantile(data[[trait_var]], 0.98, na.rm = TRUE),
    length.out = 220
  )

  nd0 <- data.frame(UI2 = factor(BASE_UI2, levels = levels(data$UI2)),
                    StdTmeanAnomalyRS = BASE_TEMP0)
  for (v in trait_vars_cont) nd0[[v]] <- trait_med[[v]]
  ep0 <- posterior_epred(mod, newdata = nd0, re_formula = NA, ndraws = ndraws_pred)

  abs_list <- list()
  pct_list <- list()

  for (i in seq_along(temp_values)) {
    tmp <- temp_values[i]
    nd <- expand.grid(UI2 = levels(data$UI2), x = xseq)
    nd$UI2 <- factor(nd$UI2, levels = levels(data$UI2))
    nd[[trait_var]] <- nd$x
    nd$x <- NULL
    nd$StdTmeanAnomalyRS <- tmp
    for (v in trait_vars_cont) {
      if (!v %in% names(nd)) nd[[v]] <- trait_med[[v]]
      if (v != trait_var) nd[[v]] <- trait_med[[v]]
    }

    ep <- posterior_epred(mod, newdata = nd, re_formula = NA, ndraws = ndraws_pred)
    abs_list[[i]] <- posterior_summary_df(ep, nd) %>% dplyr::mutate(temp_label = unname(temp_labels[as.character(tmp)]))
    pct_list[[i]] <- posterior_summary_df(percent_change_draws(ep, ep0), nd) %>% dplyr::mutate(temp_label = unname(temp_labels[as.character(tmp)]))
  }

  df_abs <- dplyr::bind_rows(abs_list)
  df_pct <- dplyr::bind_rows(pct_list)

  p_abs <- ggplot(df_abs, aes(x = .data[[trait_var]], y = mean, colour = UI2, fill = UI2)) +
    geom_ribbon(aes(ymin = low95, ymax = high95), alpha = 0.20, colour = NA) +
    geom_line(linewidth = 1.00) +
    facet_wrap(~ temp_label, nrow = 1) +
    scale_colour_manual(values = pal_UI2, drop = FALSE) +
    scale_fill_manual(values = pal_UI2, drop = FALSE) +
    theme_paper(BASE_SIZE) +
    labs(
      title = paste0("Joint multitrait Bayesian screen: land use × ", trait_labels[[trait_var]]),
      subtitle = "All other continuous traits are fixed at their sample medians.",
      x = trait_labels[[trait_var]],
      y = "Predicted probability of occurrence",
      colour = "Land-use type",
      fill = "Land-use type"
    )

  p_pct <- ggplot(df_pct, aes(x = .data[[trait_var]], y = mean, colour = UI2, fill = UI2)) +
    geom_ribbon(aes(ymin = low95, ymax = high95), alpha = 0.20, colour = NA) +
    geom_line(linewidth = 1.00) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5) +
    facet_wrap(~ temp_label, nrow = 1) +
    scale_colour_manual(values = pal_UI2, drop = FALSE) +
    scale_fill_manual(values = pal_UI2, drop = FALSE) +
    theme_paper(BASE_SIZE) +
    labs(
      title = paste0("Joint multitrait Bayesian screen: relative land use × ", trait_labels[[trait_var]]),
      subtitle = "Percent change is relative to the primary-vegetation baseline at warming = 0 with all traits at their medians.",
      x = trait_labels[[trait_var]],
      y = "Percent change in predicted occurrence probability (%)",
      colour = "Land-use type",
      fill = "Land-use type"
    )

  list(abs_plot = p_abs, pct_plot = p_pct)
}

plot_temp_by_trait_levels_abs_pct <- function(mod, data, trait_var, ndraws_pred) {
  qvals <- quantile(data[[trait_var]], probs = c(0.25, 0.50, 0.75), na.rm = TRUE)
  qdf <- data.frame(
    q_label = factor(c("Q25", "Q50", "Q75"), levels = c("Q25", "Q50", "Q75")),
    q_value = as.numeric(qvals)
  )

  temp_seq <- seq(
    quantile(data$StdTmeanAnomalyRS, 0.02, na.rm = TRUE),
    quantile(data$StdTmeanAnomalyRS, 0.98, na.rm = TRUE),
    length.out = 220
  )

  nd0 <- data.frame(UI2 = factor(BASE_UI2, levels = levels(data$UI2)),
                    StdTmeanAnomalyRS = BASE_TEMP0)
  for (v in trait_vars_cont) nd0[[v]] <- trait_med[[v]]
  ep0 <- posterior_epred(mod, newdata = nd0, re_formula = NA, ndraws = ndraws_pred)

  nd <- expand.grid(
    UI2 = levels(data$UI2),
    StdTmeanAnomalyRS = temp_seq,
    q_label = qdf$q_label
  )
  nd$UI2 <- factor(nd$UI2, levels = levels(data$UI2))
  nd$q_label <- factor(nd$q_label, levels = levels(qdf$q_label))
  nd[[trait_var]] <- qdf$q_value[match(nd$q_label, qdf$q_label)]

  for (v in trait_vars_cont) {
    if (!v %in% names(nd)) nd[[v]] <- trait_med[[v]]
    if (v != trait_var) nd[[v]] <- trait_med[[v]]
  }

  ep <- posterior_epred(mod, newdata = nd, re_formula = NA, ndraws = ndraws_pred)
  df_abs <- posterior_summary_df(ep, nd)
  df_pct <- posterior_summary_df(percent_change_draws(ep, ep0), nd)

  p_abs <- ggplot(df_abs, aes(x = StdTmeanAnomalyRS, y = mean, colour = q_label, fill = q_label)) +
    geom_ribbon(aes(ymin = low95, ymax = high95), alpha = 0.18, colour = NA) +
    geom_line(linewidth = 1.00) +
    facet_wrap(~ UI2, nrow = 1, labeller = labeller(UI2 = label_map_UI2_wrap)) +
    scale_colour_manual(values = c("Q25" = "#1B9E77", "Q50" = "#7570B3", "Q75" = "#D95F02")) +
    scale_fill_manual(values = c("Q25" = "#1B9E77", "Q50" = "#7570B3", "Q75" = "#D95F02")) +
    theme_paper(BASE_SIZE) +
    labs(
      title = paste0("Joint multitrait Bayesian screen: warming × ", trait_labels[[trait_var]]),
      subtitle = "The focal trait is shown at its 25th, 50th, and 75th percentiles; all other continuous traits are fixed at their medians.",
      x = "Standardised mean temperature anomaly (StdTmeanAnomalyRS)",
      y = "Predicted probability of occurrence",
      colour = paste0(trait_labels[[trait_var]], " level"),
      fill = paste0(trait_labels[[trait_var]], " level")
    )

  p_pct <- ggplot(df_pct, aes(x = StdTmeanAnomalyRS, y = mean, colour = q_label, fill = q_label)) +
    geom_ribbon(aes(ymin = low95, ymax = high95), alpha = 0.18, colour = NA) +
    geom_line(linewidth = 1.00) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5) +
    facet_wrap(~ UI2, nrow = 1, labeller = labeller(UI2 = label_map_UI2_wrap)) +
    scale_colour_manual(values = c("Q25" = "#1B9E77", "Q50" = "#7570B3", "Q75" = "#D95F02")) +
    scale_fill_manual(values = c("Q25" = "#1B9E77", "Q50" = "#7570B3", "Q75" = "#D95F02")) +
    theme_paper(BASE_SIZE) +
    labs(
      title = paste0("Joint multitrait Bayesian screen: relative warming × ", trait_labels[[trait_var]]),
      subtitle = "Percent change is relative to primary vegetation at warming = 0 with all traits at their medians.",
      x = "Standardised mean temperature anomaly (StdTmeanAnomalyRS)",
      y = "Percent change in predicted occurrence probability (%)",
      colour = paste0(trait_labels[[trait_var]], " level"),
      fill = paste0(trait_labels[[trait_var]], " level")
    )

  list(abs_plot = p_abs, pct_plot = p_pct)
}

## ============================================================
## 8) Predictions and plotting
## ============================================================

ndraws_pred_use <- min(N_DRAWS_PRED, posterior::ndraws(m_joint))

res_ui2temp <- plot_UI2_by_temp_abs_pct(m_joint, myoccdata, ndraws_pred_use)
save_plot_3formats(res_ui2temp$abs_plot, "Fig2_Joint_multitrait_UI2_by_warming_ABS_brms", dir_plots,
                   width = 10.8, height = 6.6)
save_plot_3formats(res_ui2temp$pct_plot, "Fig3_Joint_multitrait_UI2_by_warming_PCT_brms", dir_plots,
                   width = 10.8, height = 6.6)
add_figure_caption(
  "Fig2_Joint_multitrait_UI2_by_warming_ABS_brms", dir_plots,
  "Joint multitrait Bayesian screen: land use × warming",
  "After controlling the other continuous traits at their medians, land-use types still show distinct posterior warming-response curves.",
  "在将其他连续性状固定为中位数后，不同土地利用类型仍显示出不同的后验变暖响应曲线。"
)
add_figure_caption(
  "Fig3_Joint_multitrait_UI2_by_warming_PCT_brms", dir_plots,
  "Joint multitrait Bayesian screen: relative land use × warming effect",
  "Relative change from the primary-vegetation baseline differs among land-use types even in the joint multitrait Bayesian model.",
  "即使在联合多性状 Bayesian 模型中，相对于原始植被基线的变化幅度在不同土地利用类型之间仍然不同。"
)

for (v in trait_vars_cont) {
  message("Plotting multitrait brms predictions for: ", v)

  r1 <- plot_UI2_by_trait_fixedTemp_abs_pct(m_joint, myoccdata, v, ndraws_pred_use)
  save_plot_3formats(r1$abs_plot, paste0("Fig_UI2x", v, "_TempFixed_ABS_brms"), dir_plots,
                     width = 12.8, height = 5.9)
  save_plot_3formats(r1$pct_plot, paste0("Fig_UI2x", v, "_TempFixed_PCT_brms"), dir_plots,
                     width = 12.8, height = 5.9)

  r2 <- plot_temp_by_trait_levels_abs_pct(m_joint, myoccdata, v, ndraws_pred_use)
  save_plot_3formats(r2$abs_plot, paste0("Fig_TempX", v, "_TraitLevels_ABS_brms"), dir_plots,
                     width = 14.2, height = 5.4)
  save_plot_3formats(r2$pct_plot, paste0("Fig_TempX", v, "_TraitLevels_PCT_brms"), dir_plots,
                     width = 14.2, height = 5.4)

  add_figure_caption(
    paste0("Fig_UI2x", v, "_TempFixed_ABS_brms"), dir_plots,
    paste0("Joint multitrait Bayesian screen: land use × ", trait_labels[[v]]),
    paste0("Conditional posterior predictions show how the land-use effect changes across the observed range of ", trait_labels[[v]], " when the other traits are held at their medians."),
    paste0("条件后验预测显示，当其他性状固定在中位数时，土地利用效应会如何随着 ", trait_labels[[v]], " 的变化而改变。")
  )
  add_figure_caption(
    paste0("Fig_TempX", v, "_TraitLevels_ABS_brms"), dir_plots,
    paste0("Joint multitrait Bayesian screen: warming × ", trait_labels[[v]]),
    paste0("Conditional posterior predictions show whether warming-response curves differ across lower, median, and higher values of ", trait_labels[[v]], "."),
    paste0("条件后验预测展示了 ", trait_labels[[v]], " 的低值、中值和高值对应的变暖响应曲线是否不同。")
  )
}

write_figure_captions()
message("Joint multitrait brms pipeline completed.")
message("Outputs saved under: ", normalizePath(out_root))
