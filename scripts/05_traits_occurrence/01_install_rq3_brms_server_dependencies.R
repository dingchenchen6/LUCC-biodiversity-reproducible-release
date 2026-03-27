###############################################################
## Install server dependencies for the RQ3 brms pipeline
###############################################################

cran_repo <- Sys.getenv("RQ3_CRAN_MIRROR", unset = "https://cloud.r-project.org")
options(repos = c(CRAN = cran_repo))

cran_pkgs <- c(
  "remotes",
  "brms",
  "posterior",
  "bayesplot",
  "ggplot2",
  "dplyr",
  "tidyr",
  "stringr",
  "Hmisc",
  "ggcorrplot",
  "officer",
  "rvg",
  "ragg"
)

missing_cran <- cran_pkgs[!vapply(cran_pkgs, requireNamespace, FUN.VALUE = logical(1), quietly = TRUE)]
if (length(missing_cran) > 0) {
  install.packages(missing_cran, dependencies = TRUE)
}

if (!requireNamespace("cmdstanr", quietly = TRUE)) {
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
  }
  remotes::install_github("stan-dev/cmdstanr", upgrade = "never", dependencies = TRUE)
}

suppressPackageStartupMessages(library(cmdstanr))

current_cmdstan <- tryCatch(cmdstanr::cmdstan_path(), error = function(e) "")
if (!nzchar(current_cmdstan) || !dir.exists(current_cmdstan)) {
  cmdstanr::install_cmdstan(
    cores = max(1L, min(8L, parallel::detectCores())),
    quiet = FALSE
  )
}

manifest <- data.frame(
  package = c(cran_pkgs, "cmdstanr", "cmdstan_path"),
  status = c(
    vapply(cran_pkgs, requireNamespace, FUN.VALUE = logical(1), quietly = TRUE),
    requireNamespace("cmdstanr", quietly = TRUE),
    TRUE
  ),
  detail = c(
    vapply(cran_pkgs, function(x) as.character(packageVersion(x)), character(1)),
    if (requireNamespace("cmdstanr", quietly = TRUE)) as.character(packageVersion("cmdstanr")) else NA_character_,
    tryCatch(cmdstanr::cmdstan_path(), error = function(e) NA_character_)
  ),
  stringsAsFactors = FALSE
)

out_dir <- "results/traits_fg/rq3_brms_server_runner"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(manifest, file.path(out_dir, "rq3_brms_server_dependency_manifest.csv"), row.names = FALSE)

cat("RQ3 brms server dependencies are ready.\n")
print(manifest)
