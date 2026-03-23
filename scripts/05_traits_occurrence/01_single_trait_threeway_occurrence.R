##%#############################################################%##
#                                                               #
####  brms (cmdstanr) Occurrence models + 2-way & 3-way plots ####
#                                                               #
##%#############################################################%##
# FULL pipeline (ABS + %Δp) with consistent baseline rule:
#
#   Two response scales produced for EVERY plot:
#     (A) Absolute probability of occurrence (p)
#     (B) Percent change in probability vs baseline (%Δp) = (p - p0)/p0 * 100
#
#   Baseline definition (same for 2-way and 3-way):
#     - Land-use type: Primary vegetation (PV)
#     - Temperature anomaly: StdTmeanAnomalyRS = 0  (≈ no anomaly)
#     - Trait level: lowest quantile group (Q1) representative value
#       (for 2-way: Q1 representative sets the baseline point; x still varies in curves)
#
# UPDATED variable meanings/labels to match the paper figure:
#   - TR: Thermal tolerance breadth (Tmax_sp - Tmin_sp)
#   - Tmin_position: site cold-end position within species tolerance (0–1)
#   - Tmax_position: site warm-end position within species tolerance (0–1)
#
# IMPORTANT:
#   - Tmin_position / Tmax_position are used as-is (NOT re-centered / NOT standardised),
#     consistent with your current modelling setup.
#
# OUTPUT FORMATS unchanged:
#   - PNG + PDF + PPTX (topptx preferred; fallback officer+rvg)
#
# Model diagnostics are placed at the END and COMMENTED OUT by default.
##%#############################################################%##

rm(list = ls())

## ---------- 0) Libraries ----------
suppressPackageStartupMessages({
  library(brms)
  library(cmdstanr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(posterior)
  library(stringr)
})

backend_opt <- "cmdstanr"

## PPTX export: topptx preferred; fallback officer+rvg
have_topptx <- exists("topptx")
if (!have_topptx) {
  if (requireNamespace("officer", quietly = TRUE) && requireNamespace("rvg", quietly = TRUE)) {
    library(officer)
    library(rvg)
  } else {
    message("NOTE: topptx() not found and officer/rvg not installed -> PPTX export will be skipped.")
  }
}

## ---------- 1) Global settings ----------
out_root <- "brms_occurrence_full_pipeline_ABS_and_PCT"
dir.create(out_root, showWarnings = FALSE, recursive = TRUE)

UI2_levels <- c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High","Urban")

pal <- c("#009E73", "#0072B2", "#E69F00", "#D55E00", "#CC79A7")
names(pal) <- UI2_levels

## Plot style toggles (easy to tweak)
CI95_ONLY   <- TRUE   # TRUE: only 95% CI; FALSE: 67% + 95%
TITLE_SIZE  <- 14     # smaller title
BASE_SIZE   <- 13     # overall theme base size

## Quantile grouping for 3-way plots (2/3/4)
Q_GROUPS_3WAY <- 3

## Baseline definition (same across ALL plots)
BASE_UI2   <- "Primary vegetation"
BASE_TEMP0 <- 0       # StdTmeanAnomalyRS = 0 (≈ no anomaly)
BASE_Q     <- 1       # baseline trait group = Q1

## Temperature scenarios for 2-way plots
TEMP_SCENARIOS <- c(0, 1)   # always include 0/1
ADD_TEMP_MEDIAN <- TRUE

## ---------- 2) Labels (UPDATED to match paper definitions) ----------
label_map <- list(
  y_prob  = "Probability of occurrence",
  y_pct   = "Change in probability of occurrence (%)",
  legend_UI2 = "Land-use type",
  temp_mean = "Standardised mean temperature anomaly",
  
  x_LogRS   = "Range size (log-transformed)",
  x_HB      = "Habitat breadth",
  x_HWI     = "Dispersal ability (HWI)",
  x_TR      = "Thermal tolerance breadth (TR)",
  
  ## Paper definition (0–1 position within tolerance)
  x_TminPos = "Thermal position within species thermal tolerance (Tmin_position, 0–1)",
  x_TmaxPos = "Thermal position based on l temperature of coldest month (Tmax_position, 0–1)",
  
  x_GL      = "Generation length",
  x_CS      = "Clutch size"
)

make_title_2way <- function(x_full, temp_full, temp_label, y_mode) {
  ytxt <- if (y_mode == "absolute_prob") "Absolute probability" else "Percent change vs baseline"
  paste0("Two-way interaction: UI2 × ", x_full,
         " (", temp_full, " fixed at ", temp_label, ")",
         " — ", ytxt)
}

make_title_3way <- function(temp_full, mod_full, q_groups, y_mode) {
  ytxt <- if (y_mode == "absolute_prob") "Absolute probability" else "Percent change vs baseline"
  paste0("Three-way interaction: ", temp_full, " × UI2 × ", mod_full,
         " (", q_groups, " quantile groups)",
         " — ", ytxt)
}

## ---------- 3) Export helper ----------
save_and_show_plot <- function(p, png_file, pdf_file = NULL, pptx_file = NULL,
                               width = 9, height = 5.5, dpi = 300) {
  print(p)  # show first
  ggsave(png_file, p, width = width, height = height, dpi = dpi)
  if (!is.null(pdf_file)) ggsave(pdf_file, p, width = width, height = height)
  
  if (!is.null(pptx_file)) {
    if (exists("topptx")) {
      topptx(p, pptx_file)
    } else if (requireNamespace("officer", quietly = TRUE) && requireNamespace("rvg", quietly = TRUE)) {
      doc <- officer::read_pptx()
      doc <- officer::add_slide(doc, layout = "Title and Content", master = "Office Theme")
      doc <- officer::ph_with(doc, rvg::dml(ggobj = p), location = officer::ph_location_fullsize())
      print(doc, target = pptx_file)
    } else {
      message("PPTX export skipped: ", pptx_file)
    }
  }
  invisible(TRUE)
}

## ---------- 4) Core utilities ----------
summ_draws <- function(draws_vec, ci95_only = TRUE) {
  if (isTRUE(ci95_only)) {
    c(mean  = mean(draws_vec),
      low95 = unname(quantile(draws_vec, 0.025)),
      high95= unname(quantile(draws_vec, 0.975)))
  } else {
    c(mean  = mean(draws_vec),
      low67 = unname(quantile(draws_vec, 0.16667)),
      high67= unname(quantile(draws_vec, 0.83333)),
      low95 = unname(quantile(draws_vec, 0.025)),
      high95= unname(quantile(draws_vec, 0.975)))
  }
}

get_required_predictors <- function(formula_obj) {
  v <- all.vars(formula_obj)
  v <- setdiff(v, c("Occur"))
  unique(v)
}

fill_missing_predictors <- function(newdata, data, required_vars) {
  miss <- setdiff(required_vars, names(newdata))
  if (length(miss) == 0) return(newdata)
  
  for (v in miss) {
    if (!v %in% names(data)) next
    if (is.numeric(data[[v]])) {
      newdata[[v]] <- median(data[[v]], na.rm = TRUE)
    } else if (is.factor(data[[v]])) {
      newdata[[v]] <- factor(levels(data[[v]])[1], levels = levels(data[[v]]))
    } else {
      newdata[[v]] <- data[[v]][which(!is.na(data[[v]]))[1]]
    }
  }
  newdata
}

## Robust standardisation with external params table
make_rs_with_params <- function(df, source_name, rs_name, params = list()) {
  if (!source_name %in% names(df)) stop("Source column not found: ", source_name)
  
  x <- df[[source_name]]
  mu <- mean(x, na.rm = TRUE)
  sdv <- stats::sd(x, na.rm = TRUE)
  
  if (!is.finite(mu) || !is.finite(sdv) || sdv == 0) {
    stop("Cannot standardise ", source_name, ": mean/sd not finite or sd=0.")
  }
  
  df[[rs_name]] <- (x - mu) / sdv
  params[[rs_name]] <- list(source = source_name, center = mu, scale = sdv)
  list(df = df, params = params)
}

back_transform_x_params <- function(df, x_var_rs, x_plot_name, params,
                                    mode = c("identity","log_only","exp","square")) {
  mode <- match.arg(mode)
  if (!x_var_rs %in% names(params)) {
    stop("No standardisation parameters found for ", x_var_rs,
         ". Did you create it via make_rs_with_params() ?")
  }
  mu  <- params[[x_var_rs]]$center
  sdv <- params[[x_var_rs]]$scale
  x_bt <- df[[x_var_rs]] * sdv + mu
  
  if (mode == "identity" || mode == "log_only") {
    df[[x_plot_name]] <- x_bt
  } else if (mode == "exp") {
    df[[x_plot_name]] <- exp(x_bt)
  } else if (mode == "square") {
    df[[x_plot_name]] <- x_bt^2
  }
  df
}

ensure_sqrt <- function(df, raw_name, out_name) {
  if (out_name %in% names(df)) return(df)
  if (!raw_name %in% names(df)) stop("Column not found: ", raw_name)
  df[[out_name]] <- sqrt(df[[raw_name]])
  df
}

## Percent change (draw-wise) relative to baseline
percent_change_draws <- function(ep, ep0) {
  ep0v <- as.vector(ep0)
  ep0v[ep0v < 1e-9] <- 1e-9
  sweep(ep, 1, ep0v, FUN = function(p, p0) (p - p0) / p0 * 100)
}

make_baseline_newdata <- function(data, required_vars,
                                  UI2_baseline = BASE_UI2,
                                  temp_var = "StdTmeanAnomalyRS",
                                  temp0 = BASE_TEMP0,
                                  mod_var = NULL,
                                  mod_value = NULL) {
  nd0 <- data.frame(UI2 = factor(UI2_baseline, levels = levels(data$UI2)))
  nd0[[temp_var]] <- temp0
  if (!is.null(mod_var) && !is.null(mod_value)) nd0[[mod_var]] <- mod_value
  nd0 <- fill_missing_predictors(nd0, data, required_vars)
  nd0
}

## Quantile groups: representative value per bin = median within bin
make_quantile_groups <- function(x, q_groups = 3) {
  stopifnot(q_groups %in% c(2,3,4))
  probs <- seq(0, 1, length.out = q_groups + 1)
  brks <- unique(quantile(x, probs = probs, na.rm = TRUE, type = 7))
  
  if (length(brks) < (q_groups + 1)) {
    brks <- unique(pretty(x, n = q_groups))
    brks <- brks[is.finite(brks)]
    if (length(brks) < 3) stop("Not enough unique values to form groups.")
  }
  
  grp <- cut(x, breaks = brks, include.lowest = TRUE, right = TRUE)
  reps <- tapply(x, grp, function(z) median(z, na.rm = TRUE))
  list(reps = reps)
}

## ---------- 5) Data exploration plots ----------
plot_distribution <- function(df, var, xlab, outdir, title = NULL) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  
  if (is.null(title)) title <- paste0("Distribution of ", var)
  
  p <- ggplot(df, aes(x = .data[[var]])) +
    geom_histogram(bins = 60) +
    theme_classic(base_size = BASE_SIZE) +
    labs(title = title, x = xlab, y = "Count")
  
  save_and_show_plot(
    p,
    png_file  = file.path(outdir, paste0("Distribution_", var, ".png")),
    pdf_file  = file.path(outdir, paste0("Distribution_", var, ".pdf")),
    pptx_file = file.path(outdir, paste0("Distribution_", var, ".pptx")),
    width = 8.5, height = 5
  )
  invisible(p)
}

plot_density_by_UI2 <- function(df, var, title, xlab, outdir) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  p <- ggplot(df, aes(x = .data[[var]], color = UI2)) +
    geom_density(linewidth = 1) +
    scale_color_manual(values = pal, name = label_map$legend_UI2) +
    theme_classic(base_size = BASE_SIZE) +
    labs(title = title, x = xlab, y = "Density")
  save_and_show_plot(
    p,
    png_file  = file.path(outdir, paste0("DensityByUI2_", var, ".png")),
    pdf_file  = file.path(outdir, paste0("DensityByUI2_", var, ".pdf")),
    pptx_file = file.path(outdir, paste0("DensityByUI2_", var, ".pptx")),
    width = 9, height = 5.2
  )
  invisible(p)
}

## ---------- 6) 2-way plot function (ABS + %Δp) ----------
plot_two_way_UI2_x_multiTemp_backx_ABS_or_PCT <- function(
    fit, data, std_params,
    tag_prefix, outdir,
    x_var_rs,
    x_plot_name,
    x_full_label,
    x_back_mode = c("identity","log_only","exp","square"),
    temp_var_rs = "StdTmeanAnomalyRS",
    temp_full_label = label_map$temp_mean,
    temp_values = c(0, 1),
    add_temp_median = TRUE,
    n = 220,
    y_mode = c("absolute_prob","percent_change"),
    UI2_baseline = BASE_UI2,
    temp0 = BASE_TEMP0,
    baseline_q = BASE_Q,
    ci95_only = CI95_ONLY
) {
  y_mode <- match.arg(y_mode)
  x_back_mode <- match.arg(x_back_mode)
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  
  fml <- formula(fit)
  required <- get_required_predictors(fml)
  
  ## x sequence (on model scale)
  x_seq_rs <- seq(
    quantile(data[[x_var_rs]], 0.02, na.rm = TRUE),
    quantile(data[[x_var_rs]], 0.98, na.rm = TRUE),
    length.out = n
  )
  
  ## Temperature scenarios
  tv <- temp_values
  temp_med <- median(data[[temp_var_rs]], na.rm = TRUE)
  if (add_temp_median) tv <- sort(unique(c(tv, temp_med)))
  
  temp_labels <- sapply(tv, function(v) {
    if (isTRUE(all.equal(v, 0))) return("0")
    if (isTRUE(all.equal(v, 1))) return("1")
    if (isTRUE(all.equal(v, temp_med))) return("the median")
    format(signif(v, 3), scientific = FALSE)
  }, USE.NAMES = FALSE)
  
  ## Baseline trait value = Q1 representative on MODEL SCALE (x_var_rs)
  qobj <- make_quantile_groups(data[[x_var_rs]], q_groups = 3)  # use 3 groups just to define Q1 robustly
  reps <- qobj$reps
  base_x_value <- as.numeric(reps[baseline_q])
  if (!is.finite(base_x_value)) base_x_value <- quantile(data[[x_var_rs]], 0.25, na.rm = TRUE)
  
  ## Baseline epred (single point)
  nd0 <- make_baseline_newdata(
    data = data, required_vars = required,
    UI2_baseline = UI2_baseline,
    temp_var = temp_var_rs, temp0 = temp0,
    mod_var = x_var_rs, mod_value = base_x_value
  )
  ep0 <- posterior_epred(fit, newdata = nd0, re_formula = NA)  # draws x 1
  
  out_list <- vector("list", length(tv))
  names(out_list) <- paste0("temp_fixed_", temp_labels)
  
  for (i in seq_along(tv)) {
    temp0_i <- tv[i]
    tlabel  <- temp_labels[i]
    
    newdata <- expand.grid(
      UI2 = levels(data$UI2),
      x = x_seq_rs
    )
    names(newdata) <- c("UI2", x_var_rs)
    newdata[[temp_var_rs]] <- temp0_i
    newdata <- fill_missing_predictors(newdata, data, required)
    
    ep <- posterior_epred(fit, newdata = newdata, re_formula = NA)
    
    ep_use <- if (y_mode == "percent_change") percent_change_draws(ep, ep0) else ep
    
    summ <- apply(ep_use, 2, summ_draws, ci95_only = ci95_only)
    summ <- as.data.frame(t(summ))
    df <- cbind(newdata, summ)
    
    ## back-transform x for plotting
    df <- back_transform_x_params(df, x_var_rs = x_var_rs, x_plot_name = x_plot_name,
                                  params = std_params, mode = x_back_mode)
    
    ylab_use <- if (y_mode == "percent_change") label_map$y_pct else label_map$y_prob
    
    p <- ggplot(df, aes(x = .data[[x_plot_name]], group = UI2)) +
      {
        if (isTRUE(ci95_only)) {
          geom_ribbon(aes(ymin = low95, ymax = high95, fill = UI2), alpha = 0.22, colour = NA)
        } else {
          list(
            geom_ribbon(aes(ymin = low95, ymax = high95, fill = UI2), alpha = 0.18, colour = NA),
            geom_ribbon(aes(ymin = low67, ymax = high67, fill = UI2), alpha = 0.35, colour = NA)
          )
        }
      } +
      geom_line(aes(y = mean, color = UI2), linewidth = 1.05) +
      scale_fill_manual(values = pal, name = label_map$legend_UI2) +
      scale_color_manual(values = pal, name = label_map$legend_UI2) +
      theme_classic(base_size = BASE_SIZE) +
      theme(
        legend.position = "right",
        plot.title = element_text(size = TITLE_SIZE, face = "bold")
      ) +
      labs(
        title = make_title_2way(x_full_label, temp_full_label, tlabel, y_mode),
        x = x_full_label,
        y = ylab_use
      )
    
    if (y_mode == "percent_change") {
      p <- p + geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5)
    }
    
    file_tag <- paste0(tag_prefix,
                       "_TempFixed_", str_replace_all(tlabel, "\\s+", "_"),
                       ifelse(y_mode == "percent_change","_PCT","_ABS"))
    
    save_and_show_plot(
      p,
      png_file  = file.path(outdir, paste0(file_tag, ".png")),
      pdf_file  = file.path(outdir, paste0(file_tag, ".pdf")),
      pptx_file = file.path(outdir, paste0(file_tag, ".pptx")),
      width = 9.8, height = 5.6
    )
    
    out_list[[i]] <- list(plot = p, pred = df, temp_fixed = temp0_i, temp_label = tlabel,
                          y_mode = y_mode, baseline_x_value = base_x_value)
  }
  
  out_list
}

## ---------- 7) 2-way plot function for RAW x (Tmin_position / Tmax_position) ----------
plot_two_way_UI2_x_multiTemp_raw_ABS_or_PCT <- function(
    fit, data,
    tag_prefix, outdir,
    x_var, x_full_label,
    temp_var = "StdTmeanAnomalyRS",
    temp_full_label = label_map$temp_mean,
    temp_values = c(0, 1),
    add_temp_median = TRUE,
    n = 220,
    y_mode = c("absolute_prob","percent_change"),
    UI2_baseline = BASE_UI2,
    temp0 = BASE_TEMP0,
    baseline_q = BASE_Q,
    ci95_only = CI95_ONLY
) {
  y_mode <- match.arg(y_mode)
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  
  fml <- formula(fit)
  required <- get_required_predictors(fml)
  
  ## x sequence (raw scale)
  x_seq <- seq(
    quantile(data[[x_var]], 0.02, na.rm = TRUE),
    quantile(data[[x_var]], 0.98, na.rm = TRUE),
    length.out = n
  )
  
  ## Temperature scenarios
  tv <- temp_values
  temp_med <- median(data[[temp_var]], na.rm = TRUE)
  if (add_temp_median) tv <- sort(unique(c(tv, temp_med)))
  
  temp_labels <- sapply(tv, function(v) {
    if (isTRUE(all.equal(v, 0))) return("0")
    if (isTRUE(all.equal(v, 1))) return("1")
    if (isTRUE(all.equal(v, temp_med))) return("the median")
    format(signif(v, 3), scientific = FALSE)
  }, USE.NAMES = FALSE)
  
  ## Baseline x value = Q1 representative (raw scale)
  qobj <- make_quantile_groups(data[[x_var]], q_groups = 3)
  reps <- qobj$reps
  base_x_value <- as.numeric(reps[baseline_q])
  if (!is.finite(base_x_value)) base_x_value <- quantile(data[[x_var]], 0.25, na.rm = TRUE)
  
  ## Baseline epred
  nd0 <- make_baseline_newdata(
    data = data, required_vars = required,
    UI2_baseline = UI2_baseline,
    temp_var = temp_var, temp0 = temp0,
    mod_var = x_var, mod_value = base_x_value
  )
  ep0 <- posterior_epred(fit, newdata = nd0, re_formula = NA)  # draws x 1
  
  out_list <- vector("list", length(tv))
  names(out_list) <- paste0("temp_fixed_", temp_labels)
  
  for (i in seq_along(tv)) {
    temp0_i <- tv[i]
    tlabel  <- temp_labels[i]
    
    newdata <- expand.grid(
      UI2 = levels(data$UI2),
      x = x_seq
    )
    names(newdata) <- c("UI2", x_var)
    newdata[[temp_var]] <- temp0_i
    newdata <- fill_missing_predictors(newdata, data, required)
    
    ep <- posterior_epred(fit, newdata = newdata, re_formula = NA)
    ep_use <- if (y_mode == "percent_change") percent_change_draws(ep, ep0) else ep
    
    summ <- apply(ep_use, 2, summ_draws, ci95_only = ci95_only)
    summ <- as.data.frame(t(summ))
    df <- cbind(newdata, summ)
    
    ylab_use <- if (y_mode == "percent_change") label_map$y_pct else label_map$y_prob
    
    p <- ggplot(df, aes(x = .data[[x_var]], group = UI2)) +
      {
        if (isTRUE(ci95_only)) {
          geom_ribbon(aes(ymin = low95, ymax = high95, fill = UI2), alpha = 0.22, colour = NA)
        } else {
          list(
            geom_ribbon(aes(ymin = low95, ymax = high95, fill = UI2), alpha = 0.18, colour = NA),
            geom_ribbon(aes(ymin = low67, ymax = high67, fill = UI2), alpha = 0.35, colour = NA)
          )
        }
      } +
      geom_line(aes(y = mean, color = UI2), linewidth = 1.05) +
      scale_fill_manual(values = pal, name = label_map$legend_UI2) +
      scale_color_manual(values = pal, name = label_map$legend_UI2) +
      theme_classic(base_size = BASE_SIZE) +
      theme(
        legend.position = "right",
        plot.title = element_text(size = TITLE_SIZE, face = "bold")
      ) +
      labs(
        title = make_title_2way(x_full_label, temp_full_label, tlabel, y_mode),
        x = x_full_label,
        y = ylab_use
      )
    
    if (y_mode == "percent_change") {
      p <- p + geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5)
    }
    
    file_tag <- paste0(tag_prefix,
                       "_TempFixed_", str_replace_all(tlabel, "\\s+", "_"),
                       ifelse(y_mode == "percent_change","_PCT","_ABS"))
    
    save_and_show_plot(
      p,
      png_file  = file.path(outdir, paste0(file_tag, ".png")),
      pdf_file  = file.path(outdir, paste0(file_tag, ".pdf")),
      pptx_file = file.path(outdir, paste0(file_tag, ".pptx")),
      width = 11.5, height = 5.5
    )
    
    out_list[[i]] <- list(plot = p, pred = df, temp_fixed = temp0_i, temp_label = tlabel,
                          y_mode = y_mode, baseline_x_value = base_x_value)
  }
  
  out_list
}

## ---------- 8) 3-way plot function (ABS + %Δp) ----------
## Facet: rows = "Trait (Qk)" (no numeric ranges), cols = UI2
plot_three_way_temp_UI2_modQuantile_ABS_or_PCT <- function(
    fit, data,
    tag, outdir,
    temp_var = "StdTmeanAnomalyRS",
    temp_full_label = label_map$temp_mean,
    mod_var,
    mod_full_label,
    q_groups = 3,
    n = 220,
    y_mode = c("absolute_prob","percent_change"),
    UI2_baseline = BASE_UI2,
    temp0 = BASE_TEMP0,
    baseline_q = BASE_Q,
    ci95_only = CI95_ONLY
) {
  y_mode <- match.arg(y_mode)
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  
  fml <- formula(fit)
  required <- get_required_predictors(fml)
  
  temp_seq <- seq(
    quantile(data[[temp_var]], 0.02, na.rm = TRUE),
    quantile(data[[temp_var]], 0.98, na.rm = TRUE),
    length.out = n
  )
  
  qobj <- make_quantile_groups(data[[mod_var]], q_groups = q_groups)
  reps <- qobj$reps
  
  ## ensure order
  reps <- reps[order(seq_along(reps))]
  q_labels <- paste0(mod_full_label, " (Q", seq_len(length(reps)), ")")
  
  newdata <- expand.grid(
    UI2 = levels(data$UI2),
    Q = seq_len(length(reps)),
    tmp = temp_seq
  )
  names(newdata)[names(newdata) == "tmp"] <- temp_var
  newdata[[mod_var]] <- as.numeric(reps[newdata$Q])
  newdata <- fill_missing_predictors(newdata, data, required)
  
  ep <- posterior_epred(fit, newdata = newdata, re_formula = NA)
  
  base_mod_value <- as.numeric(reps[pmax(1, pmin(length(reps), baseline_q))])
  if (!is.finite(base_mod_value)) base_mod_value <- quantile(data[[mod_var]], 0.25, na.rm = TRUE)
  
  nd0 <- make_baseline_newdata(
    data = data, required_vars = required,
    UI2_baseline = UI2_baseline,
    temp_var = temp_var, temp0 = temp0,
    mod_var = mod_var, mod_value = base_mod_value
  )
  ep0 <- posterior_epred(fit, newdata = nd0, re_formula = NA)
  
  ep_use <- if (y_mode == "percent_change") percent_change_draws(ep, ep0) else ep
  
  summ <- apply(ep_use, 2, summ_draws, ci95_only = ci95_only)
  summ <- as.data.frame(t(summ))
  df <- cbind(newdata, summ)
  df$Qlab <- factor(df$Q, levels = seq_len(length(reps)), labels = q_labels)
  
  ylab_use <- if (y_mode == "percent_change") label_map$y_pct else label_map$y_prob
  
  p <- ggplot(df, aes(x = .data[[temp_var]])) +
    {
      if (isTRUE(ci95_only)) {
        geom_ribbon(aes(ymin = low95, ymax = high95, fill = UI2), alpha = 0.22, colour = NA)
      } else {
        list(
          geom_ribbon(aes(ymin = low95, ymax = high95, fill = UI2), alpha = 0.18, colour = NA),
          geom_ribbon(aes(ymin = low67, ymax = high67, fill = UI2), alpha = 0.35, colour = NA)
        )
      }
    } +
    geom_line(aes(y = mean, color = UI2), linewidth = 1.10) +
    facet_grid(Qlab ~ UI2) +
    scale_color_manual(values = pal, name = label_map$legend_UI2) +
    scale_fill_manual(values = pal, name = label_map$legend_UI2) +
    theme_classic(base_size = BASE_SIZE) +
    theme(
      legend.position = "right",
      plot.title = element_text(size = TITLE_SIZE, face = "bold"),
      strip.background = element_rect(fill = "white", color = "black"),
      strip.text = element_text(face = "bold")
    ) +
    labs(
      title = make_title_3way(temp_full_label, mod_full_label, q_groups, y_mode),
      x = temp_full_label,
      y = ylab_use
    )
  
  if (y_mode == "percent_change") {
    p <- p + geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5)
  }
  
  save_and_show_plot(
    p,
    png_file  = file.path(outdir, paste0(tag, ".png")),
    pdf_file  = file.path(outdir, paste0(tag, ".pdf")),
    pptx_file = file.path(outdir, paste0(tag, ".pptx")),
    width = 12.2, height = 6.8
  )
  
  list(plot = p, pred = df, y_mode = y_mode, q_groups = q_groups,
       baseline_UI2 = UI2_baseline, baseline_temp = temp0, baseline_q = baseline_q)
}

## ---------- 9) Load data ----------
myoccdata <- read.csv("myoccdata.csv")  # <<< change path if needed
myoccdata$Occur <- as.numeric(myoccdata$Occur)
myoccdata$UI2 <- factor(myoccdata$UI2, levels = UI2_levels)
length(unique(myoccdata$Best_guess_binomial))
length(unique(myoccdata$SSBS))

## Random effect IDs must exist
stopifnot(all(c("SS", "SSBS", "Best_guess_binomial") %in% names(myoccdata)))

## ---------- 10) Temperature anomaly standardisation ----------
stopifnot("StdTmeanAnomaly" %in% names(myoccdata))
mu_t <- mean(myoccdata$StdTmeanAnomaly, na.rm = TRUE)
sd_t <- sd(myoccdata$StdTmeanAnomaly, na.rm = TRUE)
if (!is.finite(mu_t) || !is.finite(sd_t) || sd_t == 0) stop("StdTmeanAnomaly cannot be standardised (sd=0 or non-finite).")
myoccdata$StdTmeanAnomalyRS <- (myoccdata$StdTmeanAnomaly - mu_t) / sd_t

## ---------- 11) Build transformed columns (keep your original logic) ----------
stopifnot(all(c("LogRS","HB","TR","LogGL","LogCS") %in% names(myoccdata)))
stopifnot(all(c("Tmin_position","Tmax_position") %in% names(myoccdata)))

## HWI support (optional but required for plots below)
stopifnot("HWI" %in% names(myoccdata))
## if you already have LogHWI in your data, keep it; else create it from HWI.
if (!"LogHWI" %in% names(myoccdata)) {
  myoccdata$LogHWI <- myoccdata$HWI
}

## Keep unchanged: previously you had LogTR <- TR
myoccdata$LogTR <- myoccdata$TR

## Habitat breadth sqrt transform
myoccdata <- ensure_sqrt(myoccdata, raw_name = "HB", out_name = "SqrtHB")

## ---------- 12) Standardise predictors with params table ----------
std_params <- list()
tmp <- make_rs_with_params(myoccdata, "LogRS",   "RS.rs",  std_params); myoccdata <- tmp$df; std_params <- tmp$params
tmp <- make_rs_with_params(myoccdata, "SqrtHB",  "HB.rs",  std_params); myoccdata <- tmp$df; std_params <- tmp$params
tmp <- make_rs_with_params(myoccdata, "LogTR",   "TR.rs",  std_params); myoccdata <- tmp$df; std_params <- tmp$params
tmp <- make_rs_with_params(myoccdata, "LogHWI",  "HWI.rs", std_params); myoccdata <- tmp$df; std_params <- tmp$params
tmp <- make_rs_with_params(myoccdata, "LogGL",   "GL.rs",  std_params); myoccdata <- tmp$df; std_params <- tmp$params
tmp <- make_rs_with_params(myoccdata, "LogCS",   "CS.rs",  std_params); myoccdata <- tmp$df; std_params <- tmp$params

## ---------- 13) Data checks + exploration (view then export) ----------
expl_dir <- file.path(out_root, "00_DataExploration")
dir.create(expl_dir, showWarnings = FALSE, recursive = TRUE)

sink(file.path(expl_dir, "DataSummary.txt"))
cat("N =", nrow(myoccdata), "\n\n")
cat("UI2 counts:\n"); print(table(myoccdata$UI2, useNA = "ifany"))
cat("\nOccur counts:\n"); print(table(myoccdata$Occur, useNA = "ifany"))
cat("\nUnique Best_guess_binomial levels:\n"); print(length(unique(myoccdata$Best_guess_binomial)))
sink()

## Distributions
plot_distribution(myoccdata, "StdTmeanAnomalyRS",
                  xlab  = label_map$temp_mean,
                  outdir= expl_dir,
                  title = "Distribution of standardised mean temperature anomaly")

plot_distribution(myoccdata, "LogRS", xlab = label_map$x_LogRS, outdir = expl_dir,
                  title = "Distribution of log-transformed range size (raw)")

plot_distribution(myoccdata, "HB", xlab = label_map$x_HB, outdir = expl_dir)

plot_distribution(myoccdata, "TR", xlab = label_map$x_TR, outdir = expl_dir,
                  title = "Distribution of thermal tolerance breadth (TR, raw)")

plot_distribution(myoccdata, "HWI", xlab = label_map$x_HWI, outdir = expl_dir,
                  title = "Distribution of dispersal ability (HWI, raw)")

plot_distribution(myoccdata, "Tmin_position", xlab = label_map$x_TminPos, outdir = expl_dir,
                  title = "Distribution of Tmin_position (0–1 cold-end position)")

plot_distribution(myoccdata, "Tmax_position", xlab = label_map$x_TmaxPos, outdir = expl_dir,
                  title = "Distribution of Tmax_position (0–1 warm-end position)")

## Densities by UI2
plot_density_by_UI2(myoccdata, "StdTmeanAnomalyRS",
                    "Density of temperature anomaly by land-use type",
                    label_map$temp_mean, expl_dir)

plot_density_by_UI2(myoccdata, "RS.rs",
                    "Density of standardised log range size by land-use type",
                    "Standardised log range size (RS.rs)", expl_dir)

plot_density_by_UI2(myoccdata, "HB.rs",
                    "Density of standardised habitat breadth by land-use type",
                    "Standardised habitat breadth (HB.rs)", expl_dir)

plot_density_by_UI2(myoccdata, "TR.rs",
                    "Density of standardised TR by land-use type",
                    "Standardised TR (TR.rs)", expl_dir)

plot_density_by_UI2(myoccdata, "HWI.rs",
                    "Density of standardised dispersal ability by land-use type",
                    "Standardised HWI (HWI.rs)", expl_dir)

plot_density_by_UI2(myoccdata, "GL.rs",
                    "Density of standardised generation length by land-use type",
                    "Standardised GL (GL.rs)", expl_dir)

plot_density_by_UI2(myoccdata, "CS.rs",
                    "Density of standardised clutch size by land-use type",
                    "Standardised CS (CS.rs)", expl_dir)

plot_density_by_UI2(myoccdata, "Tmin_position",
                    "Density of Tmin_position (0–1) by land-use type",
                    label_map$x_TminPos, expl_dir)

plot_density_by_UI2(myoccdata, "Tmax_position",
                    "Density of Tmax_position (0–1) by land-use type",
                    label_map$x_TmaxPos, expl_dir)

## ---------- 14) Load fitted models (recommended: do NOT refit here) ----------
## >>> Replace paths to your own RDS files <<<
m_RS_Tmean      <- readRDS("D:/DCC/BIRDLIFE1/brm_OccmodelRSTmean_ub.rds")
m_HB_Tmean      <- readRDS("D:/DCC/BIRDLIFE1/brm_OccmodelHBTmean_ub.rds")
m_TR_Tmean      <- readRDS("D:/DCC/BIRDLIFE1/brm_OccmodelTRTmean_ub.rds")
m_HWI_Tmean     <- readRDS("D:/DCC/BIRDLIFE1/brm_OccmodelHWITmean_ub.rds")
m_GL_Tmean      <- readRDS("D:/DCC/BIRDLIFE1/brm_OccmodelGLTmean_ub.rds")
m_CS_Tmean      <- readRDS("D:/DCC/BIRDLIFE1/brm_OccmodelCSTmean_ub.rds")
m_TminPos_Tmean <- readRDS("D:/DCC/BIRDLIFE1/brm_OccmodelTminPTmean_ub.rds")
m_TmaxPos_Tmean <- readRDS("D:/DCC/BIRDLIFE1/brm_OccmodelTmaxPTmean_ub.rds")

## ---------- 15) Interaction plots (ABS + %Δp) ----------
plot_dir <- file.path(out_root, "01_Plots_ABS_and_PCT")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

plot_dir_3way <- file.path(out_root, "02_ThreeWay_ABS_and_PCT")
dir.create(plot_dir_3way, showWarnings = FALSE, recursive = TRUE)

## ============================================================
## 15.1 Two-way plots (standardised x + back-transform)
##      (for every trait: produce ABS + PCT)
## ============================================================

## Range size (log scale shown)
two_RS_ABS <- plot_two_way_UI2_x_multiTemp_backx_ABS_or_PCT(
  fit = m_RS_Tmean, data = myoccdata, std_params = std_params,
  tag_prefix = "TwoWay_UI2_x_RangeSize_LogScale",
  outdir = plot_dir,
  x_var_rs = "RS.rs",
  x_plot_name = "LogRS_plot",
  x_full_label = label_map$x_LogRS,
  x_back_mode = "log_only",
  temp_var_rs = "StdTmeanAnomalyRS",
  temp_full_label = label_map$temp_mean,
  temp_values = TEMP_SCENARIOS,
  add_temp_median = ADD_TEMP_MEDIAN,
  y_mode = "absolute_prob",
  ci95_only = CI95_ONLY
)
two_RS_PCT <- plot_two_way_UI2_x_multiTemp_backx_ABS_or_PCT(
  fit = m_RS_Tmean, data = myoccdata, std_params = std_params,
  tag_prefix = "TwoWay_UI2_x_RangeSize_LogScale",
  outdir = plot_dir,
  x_var_rs = "RS.rs",
  x_plot_name = "LogRS_plot",
  x_full_label = label_map$x_LogRS,
  x_back_mode = "log_only",
  temp_var_rs = "StdTmeanAnomalyRS",
  temp_full_label = label_map$temp_mean,
  temp_values = TEMP_SCENARIOS,
  add_temp_median = ADD_TEMP_MEDIAN,
  y_mode = "percent_change",
  ci95_only = CI95_ONLY
)

## Habitat breadth (back to HB on original scale via square of sqrtHB)
two_HB_ABS <- plot_two_way_UI2_x_multiTemp_backx_ABS_or_PCT(
  fit = m_HB_Tmean, data = myoccdata, std_params = std_params,
  tag_prefix = "TwoWay_UI2_x_HabitatBreadth",
  outdir = plot_dir,
  x_var_rs = "HB.rs",
  x_plot_name = "HB_plot",
  x_full_label = label_map$x_HB,
  x_back_mode = "square",
  temp_var_rs = "StdTmeanAnomalyRS",
  temp_full_label = label_map$temp_mean,
  temp_values = TEMP_SCENARIOS,
  add_temp_median = ADD_TEMP_MEDIAN,
  y_mode = "absolute_prob",
  ci95_only = CI95_ONLY
)
two_HB_PCT <- plot_two_way_UI2_x_multiTemp_backx_ABS_or_PCT(
  fit = m_HB_Tmean, data = myoccdata, std_params = std_params,
  tag_prefix = "TwoWay_UI2_x_HabitatBreadth",
  outdir = plot_dir,
  x_var_rs = "HB.rs",
  x_plot_name = "HB_plot",
  x_full_label = label_map$x_HB,
  x_back_mode = "square",
  temp_var_rs = "StdTmeanAnomalyRS",
  temp_full_label = label_map$temp_mean,
  temp_values = TEMP_SCENARIOS,
  add_temp_median = ADD_TEMP_MEDIAN,
  y_mode = "percent_change",
  ci95_only = CI95_ONLY
)

## TR (identity back-transform; keep your logic unchanged)
two_TR_ABS <- plot_two_way_UI2_x_multiTemp_backx_ABS_or_PCT(
  fit = m_TR_Tmean, data = myoccdata, std_params = std_params,
  tag_prefix = "TwoWay_UI2_x_TR",
  outdir = plot_dir,
  x_var_rs = "TR.rs",
  x_plot_name = "TR_plot",
  x_full_label = label_map$x_TR,
  x_back_mode = "identity",
  temp_var_rs = "StdTmeanAnomalyRS",
  temp_full_label = label_map$temp_mean,
  temp_values = TEMP_SCENARIOS,
  add_temp_median = ADD_TEMP_MEDIAN,
  y_mode = "absolute_prob",
  ci95_only = CI95_ONLY
)
two_TR_PCT <- plot_two_way_UI2_x_multiTemp_backx_ABS_or_PCT(
  fit = m_TR_Tmean, data = myoccdata, std_params = std_params,
  tag_prefix = "TwoWay_UI2_x_TR",
  outdir = plot_dir,
  x_var_rs = "TR.rs",
  x_plot_name = "TR_plot",
  x_full_label = label_map$x_TR,
  x_back_mode = "identity",
  temp_var_rs = "StdTmeanAnomalyRS",
  temp_full_label = label_map$temp_mean,
  temp_values = TEMP_SCENARIOS,
  add_temp_median = ADD_TEMP_MEDIAN,
  y_mode = "percent_change",
  ci95_only = CI95_ONLY
)

## HWI (dispersal ability) — keep unchanged logic: identity (or exp if you used log(HWI))
two_HWI_ABS <- plot_two_way_UI2_x_multiTemp_backx_ABS_or_PCT(
  fit = m_HWI_Tmean, data = myoccdata, std_params = std_params,
  tag_prefix = "TwoWay_UI2_x_HWI",
  outdir = plot_dir,
  x_var_rs = "HWI.rs",
  x_plot_name = "HWI_plot",
  x_full_label = label_map$x_HWI,
  x_back_mode = "identity",  # if LogHWI is log(HWI), change to "exp"
  temp_var_rs = "StdTmeanAnomalyRS",
  temp_full_label = label_map$temp_mean,
  temp_values = TEMP_SCENARIOS,
  add_temp_median = ADD_TEMP_MEDIAN,
  y_mode = "absolute_prob",
  ci95_only = CI95_ONLY
)
two_HWI_PCT <- plot_two_way_UI2_x_multiTemp_backx_ABS_or_PCT(
  fit = m_HWI_Tmean, data = myoccdata, std_params = std_params,
  tag_prefix = "TwoWay_UI2_x_HWI",
  outdir = plot_dir,
  x_var_rs = "HWI.rs",
  x_plot_name = "HWI_plot",
  x_full_label = label_map$x_HWI,
  x_back_mode = "identity",  # if LogHWI is log(HWI), change to "exp"
  temp_var_rs = "StdTmeanAnomalyRS",
  temp_full_label = label_map$temp_mean,
  temp_values = TEMP_SCENARIOS,
  add_temp_median = ADD_TEMP_MEDIAN,
  y_mode = "percent_change",
  ci95_only = CI95_ONLY
)

## GL
two_GL_ABS <- plot_two_way_UI2_x_multiTemp_backx_ABS_or_PCT(
  fit = m_GL_Tmean, data = myoccdata, std_params = std_params,
  tag_prefix = "TwoWay_UI2_x_GenerationLength",
  outdir = plot_dir,
  x_var_rs = "GL.rs",
  x_plot_name = "GL_plot",
  x_full_label = label_map$x_GL,
  x_back_mode = "exp",
  temp_var_rs = "StdTmeanAnomalyRS",
  temp_full_label = label_map$temp_mean,
  temp_values = TEMP_SCENARIOS,
  add_temp_median = ADD_TEMP_MEDIAN,
  y_mode = "absolute_prob",
  ci95_only = CI95_ONLY
)
two_GL_PCT <- plot_two_way_UI2_x_multiTemp_backx_ABS_or_PCT(
  fit = m_GL_Tmean, data = myoccdata, std_params = std_params,
  tag_prefix = "TwoWay_UI2_x_GenerationLength",
  outdir = plot_dir,
  x_var_rs = "GL.rs",
  x_plot_name = "GL_plot",
  x_full_label = label_map$x_GL,
  x_back_mode = "exp",
  temp_var_rs = "StdTmeanAnomalyRS",
  temp_full_label = label_map$temp_mean,
  temp_values = TEMP_SCENARIOS,
  add_temp_median = ADD_TEMP_MEDIAN,
  y_mode = "percent_change",
  ci95_only = CI95_ONLY
)

## CS
two_CS_ABS <- plot_two_way_UI2_x_multiTemp_backx_ABS_or_PCT(
  fit = m_CS_Tmean, data = myoccdata, std_params = std_params,
  tag_prefix = "TwoWay_UI2_x_ClutchSize",
  outdir = plot_dir,
  x_var_rs = "CS.rs",
  x_plot_name = "CS_plot",
  x_full_label = label_map$x_CS,
  x_back_mode = "exp",
  temp_var_rs = "StdTmeanAnomalyRS",
  temp_full_label = label_map$temp_mean,
  temp_values = TEMP_SCENARIOS,
  add_temp_median = ADD_TEMP_MEDIAN,
  y_mode = "absolute_prob",
  ci95_only = CI95_ONLY
)
two_CS_PCT <- plot_two_way_UI2_x_multiTemp_backx_ABS_or_PCT(
  fit = m_CS_Tmean, data = myoccdata, std_params = std_params,
  tag_prefix = "TwoWay_UI2_x_ClutchSize",
  outdir = plot_dir,
  x_var_rs = "CS.rs",
  x_plot_name = "CS_plot",
  x_full_label = label_map$x_CS,
  x_back_mode = "exp",
  temp_var_rs = "StdTmeanAnomalyRS",
  temp_full_label = label_map$temp_mean,
  temp_values = TEMP_SCENARIOS,
  add_temp_median = ADD_TEMP_MEDIAN,
  y_mode = "percent_change",
  ci95_only = CI95_ONLY
)

## ============================================================
## 15.2 Two-way plots (RAW x: Tmin_position / Tmax_position)
## ============================================================

two_Tmin_ABS <- plot_two_way_UI2_x_multiTemp_raw_ABS_or_PCT(
  fit = m_TminPos_Tmean, data = myoccdata,
  tag_prefix = "TwoWay_UI2_x_Tmin_position",
  outdir = plot_dir,
  x_var = "Tmin_position",
  x_full_label = label_map$x_TminPos,
  temp_var = "StdTmeanAnomalyRS",
  temp_full_label = label_map$temp_mean,
  temp_values = TEMP_SCENARIOS,
  add_temp_median = ADD_TEMP_MEDIAN,
  y_mode = "absolute_prob",
  ci95_only = CI95_ONLY
)
two_Tmin_PCT <- plot_two_way_UI2_x_multiTemp_raw_ABS_or_PCT(
  fit = m_TminPos_Tmean, data = myoccdata,
  tag_prefix = "TwoWay_UI2_x_Tmin_position",
  outdir = plot_dir,
  x_var = "Tmin_position",
  x_full_label = label_map$x_TminPos,
  temp_var = "StdTmeanAnomalyRS",
  temp_full_label = label_map$temp_mean,
  temp_values = TEMP_SCENARIOS,
  add_temp_median = ADD_TEMP_MEDIAN,
  y_mode = "percent_change",
  ci95_only = CI95_ONLY
)

two_Tmax_ABS <- plot_two_way_UI2_x_multiTemp_raw_ABS_or_PCT(
  fit = m_TmaxPos_Tmean, data = myoccdata,
  tag_prefix = "TwoWay_UI2_x_Tmax_position",
  outdir = plot_dir,
  x_var = "Tmax_position",
  x_full_label = label_map$x_TmaxPos,
  temp_var = "StdTmeanAnomalyRS",
  temp_full_label = label_map$temp_mean,
  temp_values = TEMP_SCENARIOS,
  add_temp_median = ADD_TEMP_MEDIAN,
  y_mode = "absolute_prob",
  ci95_only = CI95_ONLY
)
two_Tmax_PCT <- plot_two_way_UI2_x_multiTemp_raw_ABS_or_PCT(
  fit = m_TmaxPos_Tmean, data = myoccdata,
  tag_prefix = "TwoWay_UI2_x_Tmax_position",
  outdir = plot_dir,
  x_var = "Tmax_position",
  x_full_label = label_map$x_TmaxPos,
  temp_var = "StdTmeanAnomalyRS",
  temp_full_label = label_map$temp_mean,
  temp_values = TEMP_SCENARIOS,
  add_temp_median = ADD_TEMP_MEDIAN,
  y_mode = "percent_change",
  ci95_only = CI95_ONLY
)

## ============================================================
## 15.3 Three-way plots (quantile groups; ABS + %Δp)
##      For EVERY trait: RS, HB, TR, HWI, GL, CS, Tmin_pos, Tmax_pos
## ============================================================

## Helper to run ABS+PCT for a 3-way plot
run_threeway_ABS_PCT <- function(fit, data, tag_base, outdir, temp_var, temp_lab,
                                 mod_var, mod_lab, q_groups, ci95_only) {
  abs <- plot_three_way_temp_UI2_modQuantile_ABS_or_PCT(
    fit = fit, data = data,
    tag = paste0(tag_base, "_ABS"),
    outdir = outdir,
    temp_var = temp_var,
    temp_full_label = temp_lab,
    mod_var = mod_var,
    mod_full_label = mod_lab,
    q_groups = q_groups,
    y_mode = "absolute_prob",
    ci95_only = ci95_only
  )
  pct <- plot_three_way_temp_UI2_modQuantile_ABS_or_PCT(
    fit = fit, data = data,
    tag = paste0(tag_base, "_PCT"),
    outdir = outdir,
    temp_var = temp_var,
    temp_full_label = temp_lab,
    mod_var = mod_var,
    mod_full_label = mod_lab,
    q_groups = q_groups,
    y_mode = "percent_change",
    ci95_only = ci95_only
  )
  list(abs = abs, pct = pct)
}

## Range size
three_RS <- run_threeway_ABS_PCT(
  fit = m_RS_Tmean, data = myoccdata,
  tag_base = paste0("ThreeWay_Temp_x_UI2_x_RangeSize_Q", Q_GROUPS_3WAY),
  outdir = plot_dir_3way,
  temp_var = "StdTmeanAnomalyRS", temp_lab = label_map$temp_mean,
  mod_var = "RS.rs", mod_lab = label_map$x_LogRS,
  q_groups = Q_GROUPS_3WAY, ci95_only = CI95_ONLY
)

## Habitat breadth
three_HB <- run_threeway_ABS_PCT(
  fit = m_HB_Tmean, data = myoccdata,
  tag_base = paste0("ThreeWay_Temp_x_UI2_x_HabitatBreadth_Q", Q_GROUPS_3WAY),
  outdir = plot_dir_3way,
  temp_var = "StdTmeanAnomalyRS", temp_lab = label_map$temp_mean,
  mod_var = "HB.rs", mod_lab = label_map$x_HB,
  q_groups = Q_GROUPS_3WAY, ci95_only = CI95_ONLY
)

## TR
three_TR <- run_threeway_ABS_PCT(
  fit = m_TR_Tmean, data = myoccdata,
  tag_base = paste0("ThreeWay_Temp_x_UI2_x_TR_Q", Q_GROUPS_3WAY),
  outdir = plot_dir_3way,
  temp_var = "StdTmeanAnomalyRS", temp_lab = label_map$temp_mean,
  mod_var = "TR.rs", mod_lab = label_map$x_TR,
  q_groups = Q_GROUPS_3WAY, ci95_only = CI95_ONLY
)

## HWI
three_HWI <- run_threeway_ABS_PCT(
  fit = m_HWI_Tmean, data = myoccdata,
  tag_base = paste0("ThreeWay_Temp_x_UI2_x_HWI_Q", Q_GROUPS_3WAY),
  outdir = plot_dir_3way,
  temp_var = "StdTmeanAnomalyRS", temp_lab = label_map$temp_mean,
  mod_var = "HWI.rs", mod_lab = label_map$x_HWI,
  q_groups = Q_GROUPS_3WAY, ci95_only = CI95_ONLY
)

## GL
three_GL <- run_threeway_ABS_PCT(
  fit = m_GL_Tmean, data = myoccdata,
  tag_base = paste0("ThreeWay_Temp_x_UI2_x_GenerationLength_Q", Q_GROUPS_3WAY),
  outdir = plot_dir_3way,
  temp_var = "StdTmeanAnomalyRS", temp_lab = label_map$temp_mean,
  mod_var = "GL.rs", mod_lab = label_map$x_GL,
  q_groups = Q_GROUPS_3WAY, ci95_only = CI95_ONLY
)

## CS
three_CS <- run_threeway_ABS_PCT(
  fit = m_CS_Tmean, data = myoccdata,
  tag_base = paste0("ThreeWay_Temp_x_UI2_x_ClutchSize_Q", Q_GROUPS_3WAY),
  outdir = plot_dir_3way,
  temp_var = "StdTmeanAnomalyRS", temp_lab = label_map$temp_mean,
  mod_var = "CS.rs", mod_lab = label_map$x_CS,
  q_groups = Q_GROUPS_3WAY, ci95_only = CI95_ONLY
)

## Tmin_position (raw 0–1)
three_Tmin <- run_threeway_ABS_PCT(
  fit = m_TminPos_Tmean, data = myoccdata,
  tag_base = paste0("ThreeWay_Temp_x_UI2_x_Tmin_position_Q", Q_GROUPS_3WAY),
  outdir = plot_dir_3way,
  temp_var = "StdTmeanAnomalyRS", temp_lab = label_map$temp_mean,
  mod_var = "Tmin_position", mod_lab = label_map$x_TminPos,
  q_groups = Q_GROUPS_3WAY, ci95_only = CI95_ONLY
)

## Tmax_position (raw 0–1)
three_Tmax <- run_threeway_ABS_PCT(
  fit = m_TmaxPos_Tmean, data = myoccdata,
  tag_base = paste0("ThreeWay_Temp_x_UI2_x_Tmax_position_Q", Q_GROUPS_3WAY),
  outdir = plot_dir_3way,
  temp_var = "StdTmeanAnomalyRS", temp_lab = label_map$temp_mean,
  mod_var = "Tmax_position", mod_lab = label_map$x_TmaxPos,
  q_groups = Q_GROUPS_3WAY, ci95_only = CI95_ONLY
)

message("DONE: Data exploration + ALL traits (RS/HB/TR/HWI/GL/CS/Tmin/Tmax) 2-way and 3-way plots exported (ABS + PCT).")
message("NOTE: Model diagnostics section is at the end and currently commented out.")

## ---------- 16) Model diagnostics (PLACE AT END; COMMENTED OUT) ----------
## Keep this section for later (fast & safe); DO NOT run by default for huge N.

# run_model_diagnostics_fast <- function(
    #     fit, tag, outdir,
#     do_summary  = TRUE,
#     do_cmdstan  = TRUE,
#     do_ppcheck  = TRUE,
#     pp_ndraws   = 60,
#     pp_nsamples = 8000,
#     do_bayesR2  = TRUE,
#     do_loo      = FALSE,
#     loo_subsample_n = 20000,
#     loo_seed    = 1
# ) {
#   dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
#
#   if (isTRUE(do_summary)) {
#     sink(file.path(outdir, paste0(tag, "_summary.txt")))
#     cat("==== brms summary ====\n")
#     print(summary(fit))
#     sink()
#   }
#
#   if (isTRUE(do_cmdstan) && inherits(fit$fit, "CmdStanMCMC")) {
#     sink(file.path(outdir, paste0(tag, "_cmdstan_diagnostic_summary.txt")))
#     cat("==== cmdstanr diagnostic_summary ====\n")
#     print(fit$fit$diagnostic_summary())
#     sink()
#   }
#
#   if (isTRUE(do_ppcheck)) {
#     ppc1 <- pp_check(fit, ndraws = pp_ndraws, nsamples = pp_nsamples)
#     print(ppc1)
#     ggsave(file.path(outdir, paste0(tag, "_ppcheck.png")), ppc1, width = 7.8, height = 5.2, dpi = 300)
#     ggsave(file.path(outdir, paste0(tag, "_ppcheck.pdf")), ppc1, width = 7.8, height = 5.2)
#   }
#
#   if (isTRUE(do_bayesR2)) {
#     sink(file.path(outdir, paste0(tag, "_bayesR2.txt")))
#     cat("==== bayes_R2 ====\n")
#     print(bayes_R2(fit))
#     sink()
#   }
#
#   if (isTRUE(do_loo)) {
#     N <- nrow(fit$data)
#     sink(file.path(outdir, paste0(tag, "_loo.txt")))
#     cat("==== LOO ====\n")
#     cat("N in model data =", N, "\n")
#
#     if (N > loo_subsample_n) {
#       cat("NOTE: Using subsampled LOO with n =", loo_subsample_n, "\n")
#       set.seed(loo_seed)
#       idx <- sample.int(N, loo_subsample_n)
#       ll <- log_lik(fit, newdata = fit$data[idx, , drop = FALSE])
#       loo_res <- loo::loo(ll)
#       print(loo_res)
#     } else {
#       loo_res <- try(loo(fit), silent = TRUE)
#       if (inherits(loo_res, "try-error")) {
#         cat("LOO failed:\n")
#         print(loo_res)
#       } else {
#         print(loo_res)
#       }
#     }
#     sink()
#   }
#
#   invisible(TRUE)
# }
#
# diag_dir <- file.path(out_root, "99_Diagnostics")
# dir.create(diag_dir, showWarnings = FALSE, recursive = TRUE)
# run_model_diagnostics_fast(m_RS_Tmean,      "Occ_RS_Tmean",      diag_dir)
# run_model_diagnostics_fast(m_HB_Tmean,      "Occ_HB_Tmean",      diag_dir)
# run_model_diagnostics_fast(m_TR_Tmean,      "Occ_TR_Tmean",      diag_dir)
# run_model_diagnostics_fast(m_HWI_Tmean,     "Occ_HWI_Tmean",     diag_dir)
# run_model_diagnostics_fast(m_GL_Tmean,      "Occ_GL_Tmean",      diag_dir)
# run_model_diagnostics_fast(m_CS_Tmean,      "Occ_CS_Tmean",      diag_dir)
# run_model_diagnostics_fast(m_TminPos_Tmean, "Occ_TminPos_Tmean", diag_dir)
# run_model_diagnostics_fast(m_TmaxPos_Tmean, "Occ_TmaxPos_Tmean", diag_dir)


##%#############################################################%##
#                                                               #
####  brms Occurrence: 3-way plots for categorical TL / DC   ####
#                                                               #
##%#############################################################%##
# 目标 / Goal:
#   为两个分类变量 TL 与 DC 分别绘制三重交互：
#     UI2 × StdTmeanAnomalyRS × (TL or DC)
#
# 输出 / Output:
#   - ABS: Probability of occurrence (p)
#   - PCT: Percent change vs baseline (%Δp) = (p - p0)/p0*100
#
# 基线 / Baseline (for %Δp):
#   - UI2 = Primary vegetation
#   - StdTmeanAnomalyRS = 0
#   - TL(或DC) = 参照水平（默认选择样本最多的水平）
##%#############################################################%##

rm(list = ls())

## ---------- 0) Libraries / 加载包 ----------
suppressPackageStartupMessages({
  library(brms)
  library(cmdstanr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(posterior)
  library(stringr)
})

## ---------- 1) Global settings / 全局设置 ----------
out_root <- "brms_occurrence_TL_DC_plots"
dir.create(out_root, showWarnings = FALSE, recursive = TRUE)

UI2_levels <- c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High","Urban")

BASE_UI2   <- "Primary vegetation"
BASE_TEMP0 <- 0

## Plot style toggles / 作图开关
CI95_ONLY  <- TRUE      # TRUE: only 95% CI; FALSE: 67%+95%
TITLE_SIZE <- 14
BASE_SIZE  <- 13

pal <- c("#009E73", "#0072B2", "#E69F00", "#D55E00", "#CC79A7")
names(pal) <- UI2_levels

label_map <- list(
  legend_UI2 = "Land-use type",
  y_prob  = "Probability of occurrence",
  y_pct   = "Change in probability of occurrence (%)",
  temp_mean = "Standardised mean temperature anomaly",
  cat_TL = "Trophic level (TL)",
  cat_DC = "Dietary category (DC)"
)

cat_full_label <- function(cat_var) {
  if (cat_var == "TL") return(label_map$cat_TL)
  if (cat_var == "DC") return(label_map$cat_DC)
  return(cat_var)
}

make_title_3way_cat <- function(temp_full, cat_full, y_mode) {
  ytxt <- if (y_mode == "absolute_prob") "Absolute probability" else "Percent change vs baseline"
  paste0("Three-way interaction: ", temp_full, " × UI2 × ", cat_full, " — ", ytxt)
}

## ---------- 2) Export helper / 导出函数 ----------
have_topptx <- exists("topptx")
if (!have_topptx) {
  if (requireNamespace("officer", quietly = TRUE) && requireNamespace("rvg", quietly = TRUE)) {
    library(officer); library(rvg)
  } else {
    message("NOTE: topptx() not found and officer/rvg not installed -> PPTX export will be skipped.")
  }
}

save_and_show_plot <- function(p, png_file, pdf_file = NULL, pptx_file = NULL,
                               width = 12.2, height = 6.8, dpi = 300) {
  print(p)
  ggsave(png_file, p, width = width, height = height, dpi = dpi)
  if (!is.null(pdf_file)) ggsave(pdf_file, p, width = width, height = height)
  
  if (!is.null(pptx_file)) {
    if (exists("topptx")) {
      topptx(p, pptx_file)
    } else if (requireNamespace("officer", quietly = TRUE) && requireNamespace("rvg", quietly = TRUE)) {
      doc <- officer::read_pptx()
      doc <- officer::add_slide(doc, layout = "Title and Content", master = "Office Theme")
      doc <- officer::ph_with(doc, rvg::dml(ggobj = p), location = officer::ph_location_fullsize())
      print(doc, target = pptx_file)
    } else {
      message("PPTX export skipped: ", pptx_file)
    }
  }
  invisible(TRUE)
}

## ---------- 3) Core utilities / 核心工具 ----------
summ_draws <- function(draws_vec, ci95_only = TRUE) {
  if (isTRUE(ci95_only)) {
    c(mean  = mean(draws_vec),
      low95 = unname(quantile(draws_vec, 0.025)),
      high95= unname(quantile(draws_vec, 0.975)))
  } else {
    c(mean  = mean(draws_vec),
      low67 = unname(quantile(draws_vec, 0.16667)),
      high67= unname(quantile(draws_vec, 0.83333)),
      low95 = unname(quantile(draws_vec, 0.025)),
      high95= unname(quantile(draws_vec, 0.975)))
  }
}

get_required_predictors <- function(formula_obj) {
  v <- all.vars(formula_obj)
  v <- setdiff(v, "Occur")
  unique(v)
}

## 用中位数/第一个水平为 newdata 填充缺失变量（稳健）
fill_missing_predictors <- function(newdata, data, required_vars) {
  miss <- setdiff(required_vars, names(newdata))
  if (length(miss) == 0) return(newdata)
  
  for (v in miss) {
    if (!v %in% names(data)) next
    if (is.numeric(data[[v]])) {
      newdata[[v]] <- median(data[[v]], na.rm = TRUE)
    } else if (is.factor(data[[v]])) {
      newdata[[v]] <- factor(levels(data[[v]])[1], levels = levels(data[[v]]))
    } else {
      newdata[[v]] <- data[[v]][which(!is.na(data[[v]]))[1]]
    }
  }
  newdata
}

## %Δp（逐抽样）/ percent change per draw
percent_change_draws <- function(ep, ep0) {
  ep0v <- as.vector(ep0)
  ep0v[ep0v < 1e-9] <- 1e-9
  sweep(ep, 1, ep0v, FUN = function(p, p0) (p - p0) / p0 * 100)
}

## 自动选 baseline 类别水平：默认选样本最多的水平（也最稳定）
choose_baseline_level <- function(df, cat_var) {
  tb <- table(df[[cat_var]])
  names(tb)[which.max(tb)]
}

## ---------- 4) 3-way plot for categorical moderator / 分类变量三重交互绘图函数 ----------
## Facet: rows = CAT level, cols = UI2
plot_three_way_temp_UI2_cat_ABS_or_PCT <- function(
    fit, data,
    tag, outdir,
    cat_var = c("TL","DC"),
    temp_var = "StdTmeanAnomalyRS",
    temp_full_label = label_map$temp_mean,
    n = 220,
    
    ## response scale / 响应尺度
    y_mode = c("absolute_prob","percent_change"),
    
    ## baseline rule / 基线规则
    UI2_baseline = BASE_UI2,
    temp0 = BASE_TEMP0,
    cat_baseline = NULL,     # if NULL -> choose most frequent
    
    ## CI / 置信区间开关
    ci95_only = CI95_ONLY
) {
  y_mode <- match.arg(y_mode)
  cat_var <- match.arg(cat_var)
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  
  ## Ensure factors / 确保是因子
  if (!is.factor(data$UI2)) stop("data$UI2 must be a factor.")
  if (!is.factor(data[[cat_var]])) stop(paste0("data$", cat_var, " must be a factor."))
  
  ## Decide baseline level / 选择基线水平
  if (is.null(cat_baseline)) {
    cat_baseline <- choose_baseline_level(data, cat_var)
  }
  cat_baseline <- as.character(cat_baseline)
  if (!cat_baseline %in% levels(data[[cat_var]])) {
    stop("cat_baseline not in factor levels: ", cat_baseline)
  }
  
  ## Sequences / 温度序列
  temp_seq <- seq(
    quantile(data[[temp_var]], 0.02, na.rm = TRUE),
    quantile(data[[temp_var]], 0.98, na.rm = TRUE),
    length.out = n
  )
  
  ## Required predictors / 模型需要的变量
  required <- get_required_predictors(formula(fit))
  
  ## Build newdata grid / 构造预测网格
  newdata <- expand.grid(
    UI2 = levels(data$UI2),
    CAT = levels(data[[cat_var]]),
    tmp = temp_seq
  )
  names(newdata)[names(newdata) == "tmp"] <- temp_var
  newdata[[cat_var]] <- factor(newdata$CAT, levels = levels(data[[cat_var]]))
  newdata$CAT <- NULL
  
  ## Fill missing predictors / 填补其它协变量（中位数/第一个水平）
  newdata <- fill_missing_predictors(newdata, data, required)
  
  ## Posterior predictions / 后验预测（只用固定效应 re_formula=NA）
  ep <- posterior_epred(fit, newdata = newdata, re_formula = NA)
  
  ## Baseline point: PV + temp=0 + CAT=baseline / 基线点
  nd0 <- data.frame(
    UI2 = factor(UI2_baseline, levels = levels(data$UI2))
  )
  nd0[[temp_var]] <- temp0
  nd0[[cat_var]]  <- factor(cat_baseline, levels = levels(data[[cat_var]]))
  nd0 <- fill_missing_predictors(nd0, data, required)
  
  ep0 <- posterior_epred(fit, newdata = nd0, re_formula = NA)  # draws x 1
  
  ## Choose response scale / 选择 ABS 或 %Δp
  ep_use <- if (y_mode == "percent_change") percent_change_draws(ep, ep0) else ep
  
  ## Summarise / 汇总均值与区间
  summ <- apply(ep_use, 2, summ_draws, ci95_only = ci95_only)
  summ <- as.data.frame(t(summ))
  df <- cbind(newdata, summ)
  
  ## Labels / 标签
  ylab_use <- if (y_mode == "percent_change") label_map$y_pct else label_map$y_prob
  cat_lab  <- cat_full_label(cat_var)
  
  ## Build plot / 作图
  p <- ggplot(df, aes(x = .data[[temp_var]])) +
    {
      if (isTRUE(ci95_only)) {
        geom_ribbon(aes(ymin = low95, ymax = high95, fill = UI2), alpha = 0.22, colour = NA)
      } else {
        list(
          geom_ribbon(aes(ymin = low95, ymax = high95, fill = UI2), alpha = 0.18, colour = NA),
          geom_ribbon(aes(ymin = low67, ymax = high67, fill = UI2), alpha = 0.35, colour = NA)
        )
      }
    } +
    geom_line(aes(y = mean, color = UI2), linewidth = 1.10) +
    facet_grid(reformulate("UI2", response = cat_var)) +  # rows=cat_var, cols=UI2
    scale_color_manual(values = pal, name = label_map$legend_UI2) +
    scale_fill_manual(values = pal, name = label_map$legend_UI2) +
    theme_classic(base_size = BASE_SIZE) +
    theme(
      legend.position = "right",
      plot.title = element_text(size = TITLE_SIZE, face = "bold"),
      strip.background = element_rect(fill = "white", color = "black"),
      strip.text = element_text(face = "bold")
    ) +
    labs(
      title = paste0(
        make_title_3way_cat(temp_full_label, cat_lab, y_mode),
        "\nBaseline for %Δp: UI2=", UI2_baseline, ", Temp=", temp0, ", ", cat_var, "=", cat_baseline
      ),
      x = temp_full_label,
      y = ylab_use
    )
  
  if (y_mode == "percent_change") {
    p <- p + geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5)
  }
  
  ## Save / 导出
  file_tag <- paste0(tag, ifelse(y_mode == "percent_change","_PCT","_ABS"))
  
  save_and_show_plot(
    p,
    png_file  = file.path(outdir, paste0(file_tag, ".png")),
    pdf_file  = file.path(outdir, paste0(file_tag, ".pdf")),
    pptx_file = file.path(outdir, paste0(file_tag, ".pptx")),
    width = 12.2, height = 6.8
  )
  
  invisible(list(
    plot = p,
    pred = df,
    y_mode = y_mode,
    baseline = list(UI2 = UI2_baseline, temp = temp0, cat_var = cat_var, cat_level = cat_baseline)
  ))
}

## ---------- 5) Load data / 读取数据 ----------
## 你可以替换为 readRDS / fread 等
myoccdata <- read.csv("myoccdata.csv")
myoccdata$Occur <- as.numeric(myoccdata$Occur)
myoccdata$UI2   <- factor(myoccdata$UI2, levels = UI2_levels)

## 必要随机效应列检查
stopifnot(all(c("SS","SSBS","Best_guess_binomial") %in% names(myoccdata)))

## TL / DC must be factors / TL/DC 必须是因子
stopifnot(all(c("TL","DC") %in% names(myoccdata)))
myoccdata$TL <- factor(myoccdata$TL)
myoccdata$DC <- factor(myoccdata$DC)

## 选择参照水平（推荐：样本最多，稳定）
## Choose reference (baseline) level (recommended: most frequent)
TL_base <- choose_baseline_level(myoccdata, "TL")
DC_base <- choose_baseline_level(myoccdata, "DC")

## 也可以手动指定，比如：
## TL_base <- "Carnivore"
## DC_base <- "Insectivore"

myoccdata$TL <- relevel(myoccdata$TL, ref = TL_base)
myoccdata$DC <- relevel(myoccdata$DC, ref = DC_base)

## StdTmeanAnomalyRS 必须存在
stopifnot("StdTmeanAnomalyRS" %in% names(myoccdata))

## ---------- 6) Load fitted models / 读取已拟合模型 ----------
## 用你自己的路径替换
m_TL_Tmean <- readRDS("brm_Occmodel_TL_Tmean.rds")
m_DC_Tmean <- readRDS("brm_Occmodel_DC_Tmean.rds")

## ---------- 7) Plot output dirs / 输出目录 ----------
plot_dir_cat <- file.path(out_root, "CAT_ThreeWay_ABS_and_PCT")
dir.create(plot_dir_cat, showWarnings = FALSE, recursive = TRUE)

## ---------- 8) TL plots (ABS + PCT) / TL 三重交互作图 ----------
three_TL_ABS <- plot_three_way_temp_UI2_cat_ABS_or_PCT(
  fit  = m_TL_Tmean,
  data = myoccdata,
  tag  = "ThreeWay_Temp_x_UI2_x_TL",
  outdir = plot_dir_cat,
  cat_var = "TL",
  y_mode  = "absolute_prob",
  cat_baseline = TL_base,
  ci95_only = CI95_ONLY
)

three_TL_PCT <- plot_three_way_temp_UI2_cat_ABS_or_PCT(
  fit  = m_TL_Tmean,
  data = myoccdata,
  tag  = "ThreeWay_Temp_x_UI2_x_TL",
  outdir = plot_dir_cat,
  cat_var = "TL",
  y_mode  = "percent_change",
  cat_baseline = TL_base,
  ci95_only = CI95_ONLY
)

## ---------- 9) DC plots (ABS + PCT) / DC 三重交互作图 ----------
three_DC_ABS <- plot_three_way_temp_UI2_cat_ABS_or_PCT(
  fit  = m_DC_Tmean,
  data = myoccdata,
  tag  = "ThreeWay_Temp_x_UI2_x_DC",
  outdir = plot_dir_cat,
  cat_var = "DC",
  y_mode  = "absolute_prob",
  cat_baseline = DC_base,
  ci95_only = CI95_ONLY
)

three_DC_PCT <- plot_three_way_temp_UI2_cat_ABS_or_PCT(
  fit  = m_DC_Tmean,
  data = myoccdata,
  tag  = "ThreeWay_Temp_x_UI2_x_DC",
  outdir = plot_dir_cat,
  cat_var = "DC",
  y_mode  = "percent_change",
  cat_baseline = DC_base,
  ci95_only = CI95_ONLY
)

message("DONE: TL & DC three-way plots exported (ABS + PCT), style consistent with your pipeline.")



library(dplyr)

## 1) 从 myoccdata 提取物种名
sp_df <- myoccdata %>%
  select(Species = Best_guess_binomial) %>%
  distinct()

## 2) 读取 AVONET BirdLife 性状数据
avonet <- read.csv("AVONET1_BirdLife.csv", stringsAsFactors = FALSE)

## 3) 按物种名匹配并提取所需性状
traits_df <- sp_df %>%
  left_join(
    avonet %>%
      select(
        Species1,
        Habitat,
        Habitat.Density,
        Migration,
        Trophic.Level,
        Trophic.Niche,
        Primary.Lifestyle
      ),
    by = c("Species" = "Species1")
  )

## 4) 查找 Trophic.Niche 为 NA 的物种
na_trophic_species <- traits_df %>%
  filter(is.na(Trophic.Niche)) %>%
  select(Species)

## 5) 快速检查
n_na  <- nrow(na_trophic_species)
prop_na <- mean(is.na(traits_df$Trophic.Niche))

cat("Number of species with NA Trophic.Niche:", n_na, "\n")
cat("Proportion NA:", round(prop_na, 3), "\n")

## 6) 如需保存结果
# write.csv(traits_df, "traits_df_with_AVONET_traits.csv", row.names = FALSE)
# write.csv(na_trophic_species, "species_with_NA_TrophicNiche.csv", row.names = FALSE)
name_map <- tibble::tribble(
  ~Species, ~Species_std,
  "Pycnonotus atriceps", "Rubigula atriceps",
  "Pycnonotus erythropthalmos", "Rubigula erythropthalmos",
  "Pycnonotus cyaniventris", "Rubigula cyaniventris",
  "Hylocharis chrysura", "Chlorostilbon chrysurus",
  "Hylocharis cyanus", "Chlorostilbon cyanus",
  "Lepidopyga coeruleogularis", "Chlorestes coeruleogularis"
)


