suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(officer)
  library(rvg)
})

get_script_dir <- function() {
  cmd <- commandArgs(trailingOnly = FALSE)
  file_arg <- cmd[grepl("^--file=", cmd)]
  if (length(file_arg) > 0) {
    return(normalizePath(dirname(sub("^--file=", "", file_arg[1]))))
  }
  normalizePath(getwd())
}

get_project_root <- function(script_dir = get_script_dir()) {
  normalizePath(file.path(script_dir, "..", ".."))
}

ui2_levels <- c(
  "Primary vegetation",
  "Secondary vegetation",
  "Agriculture_Low",
  "Agriculture_High",
  "Urban"
)

ui2_cols <- c("#009E73", "#0072B2", "#E69F00", "#D55E00", "pink")
names(ui2_cols) <- ui2_levels

theme_paper <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.25, colour = "grey85"),
      panel.border = element_rect(linewidth = 0.4, colour = "black"),
      axis.ticks = element_line(linewidth = 0.3),
      legend.title = element_blank(),
      legend.background = element_blank(),
      strip.background = element_rect(fill = "grey95", colour = "black", linewidth = 0.35),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0),
      plot.subtitle = element_text(hjust = 0),
      plot.caption = element_text(hjust = 0),
      plot.margin = margin(8, 14, 8, 8)
    )
}

save_plot_3formats <- function(p, prefix, out_dir, width_mm = 183, height_mm = 110,
                               dpi_pdf = 300, dpi_png = 600, slide_title = NULL) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  pdf_path <- file.path(out_dir, paste0(prefix, ".pdf"))
  png_path <- file.path(out_dir, paste0(prefix, ".png"))
  ppt_path <- file.path(out_dir, paste0(prefix, ".pptx"))

  ggsave(pdf_path, p, width = width_mm, height = height_mm, units = "mm", dpi = dpi_pdf)
  ggsave(png_path, p, width = width_mm, height = height_mm, units = "mm", dpi = dpi_png)

  doc <- read_pptx()
  doc <- add_slide(doc, layout = "Title and Content", master = "Office Theme")
  doc <- ph_with(
    doc,
    value = if (is.null(slide_title)) prefix else slide_title,
    location = ph_location_type(type = "title")
  )
  doc <- ph_with(doc, dml(ggobj = p), location = ph_location_fullsize())
  print(doc, target = ppt_path)

  invisible(list(pdf = pdf_path, png = png_path, ppt = ppt_path))
}

record_caption <- function(registry, fig_id, file_prefix, title_cn, title_en, key_cn, key_en) {
  registry[[length(registry) + 1]] <- data.frame(
    fig_id = fig_id,
    file_prefix = file_prefix,
    title_cn = title_cn,
    title_en = title_en,
    key_points_cn = key_cn,
    key_points_en = key_en,
    stringsAsFactors = FALSE
  )
  registry
}

write_caption_registry <- function(registry, out_dir, file_stub) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  reg <- do.call(rbind, registry)
  write.csv(reg, file.path(out_dir, paste0(file_stub, ".csv")), row.names = FALSE)

  txt_path <- file.path(out_dir, paste0(file_stub, ".txt"))
  con <- file(txt_path, open = "wt")
  on.exit(close(con), add = TRUE)

  for (i in seq_len(nrow(reg))) {
    writeLines(paste0("[", reg$fig_id[i], "] ", reg$file_prefix[i]), con)
    writeLines(paste0("CN: ", reg$title_cn[i]), con)
    writeLines(paste0("EN: ", reg$title_en[i]), con)
    writeLines(paste0("Key points CN: ", reg$key_points_cn[i]), con)
    writeLines(paste0("Key points EN: ", reg$key_points_en[i]), con)
    writeLines("", con)
  }
}

coef_table_lmer <- function(mod, keep_pattern = "UI2|pair_clim_mean|StdT") {
  sm <- coef(summary(mod))
  out <- data.frame(term = rownames(sm), sm, row.names = NULL, check.names = FALSE)
  out %>% dplyr::filter(grepl(keep_pattern, term))
}

coef_table_glmmtmb <- function(mod, keep_pattern = "UI2|StdT") {
  sm <- summary(mod)$coefficients$cond
  out <- data.frame(term = rownames(sm), sm, row.names = NULL, check.names = FALSE)
  out %>% dplyr::filter(grepl(keep_pattern, term))
}

plot_fixed_effect_forest <- function(df, estimate_col = "Estimate", se_col = "Std. Error",
                                     term_col = "term", title = NULL, subtitle = NULL,
                                     xlab = "Estimate") {
  df2 <- df %>%
    mutate(
      lower = .data[[estimate_col]] - 1.96 * .data[[se_col]],
      upper = .data[[estimate_col]] + 1.96 * .data[[se_col]],
      term = factor(.data[[term_col]], levels = rev(.data[[term_col]]))
    )

  ggplot(df2, aes(x = .data[[estimate_col]], y = term)) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3) +
    geom_errorbar(aes(xmin = lower, xmax = upper), orientation = "y", width = 0.12, linewidth = 0.4) +
    geom_point(size = 2.1, colour = "black") +
    labs(title = title, subtitle = subtitle, x = xlab, y = NULL) +
    theme_paper()
}

build_ui2_support_plot <- function(df, x_var, x_lab, title, subtitle = NULL) {
  ggplot(df, aes(x = .data[[x_var]], fill = UI2)) +
    geom_density(alpha = 0.22, linewidth = 0.4) +
    scale_fill_manual(values = ui2_cols) +
    labs(title = title, subtitle = subtitle, x = x_lab, y = "Density") +
    theme_paper()
}

label_map_ui2 <- c(
  "Primary vegetation" = "Primary vegetation",
  "Secondary vegetation" = "Secondary vegetation",
  "Agriculture_Low" = "Agriculture (low)",
  "Agriculture_High" = "Agriculture (high)",
  "Urban" = "Urban"
)
