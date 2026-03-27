###############################################################
## RQ3 brms server runner
## Sequential launcher for Bayesian RQ3 modules with log manifest
###############################################################

suppressPackageStartupMessages(library(dplyr))

get_script_path <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) == 0) stop("Cannot determine script path from commandArgs().")
  normalizePath(sub("^--file=", "", file_arg[1]))
}

as_flag <- function(x, default = TRUE) {
  if (!nzchar(x)) return(default)
  tolower(x) %in% c("1", "true", "yes", "y")
}

script_path <- get_script_path()
script_dir <- dirname(script_path)
repo_root <- normalizePath(file.path(script_dir, "..", ".."), winslash = "/", mustWork = TRUE)
setwd(repo_root)
log_dir <- Sys.getenv(
  "RQ3_BRMS_LOG_DIR",
  unset = file.path(repo_root, "results", "traits_fg", "rq3_brms_server_runner")
)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

run_single <- as_flag(Sys.getenv("RQ3_RUN_SINGLE", "1"))
run_multitrait <- as_flag(Sys.getenv("RQ3_RUN_MULTITRAIT", "1"))
run_fg <- as_flag(Sys.getenv("RQ3_RUN_FG", "1"))
stop_on_error <- as_flag(Sys.getenv("RQ3_STOP_ON_ERROR", "1"))

steps <- data.frame(
  name = c("single_trait_threeway", "multitrait_joint_screen", "stable_fg"),
  script = c(
    file.path(script_dir, "01_single_trait_threeway_occurrence.R"),
    file.path(script_dir, "02_multitrait_joint_occurrence_brms.R"),
    file.path(script_dir, "04_fg_stable_brms.R")
  ),
  enabled = c(run_single, run_multitrait, run_fg),
  stringsAsFactors = FALSE
)

manifest <- list()
rscript_bin <- file.path(R.home("bin"), "Rscript")
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

for (i in seq_len(nrow(steps))) {
  step <- steps[i, ]
  if (!isTRUE(step$enabled)) next

  log_file <- file.path(log_dir, paste0(timestamp, "_", step$name, ".log"))
  start_time <- Sys.time()
  exit_status <- system2(
    rscript_bin,
    args = step$script,
    stdout = log_file,
    stderr = log_file
  )
  end_time <- Sys.time()

  manifest[[length(manifest) + 1]] <- data.frame(
    step = step$name,
    script = normalizePath(step$script, winslash = "/", mustWork = FALSE),
    log_file = normalizePath(log_file, winslash = "/", mustWork = FALSE),
    start_time = format(start_time, "%Y-%m-%d %H:%M:%S"),
    end_time = format(end_time, "%Y-%m-%d %H:%M:%S"),
    exit_status = exit_status,
    stringsAsFactors = FALSE
  )

  if (!identical(exit_status, 0L) && isTRUE(stop_on_error)) {
    break
  }
}

manifest_df <- dplyr::bind_rows(manifest)
manifest_file <- file.path(log_dir, paste0(timestamp, "_rq3_brms_server_manifest.csv"))
utils::write.csv(manifest_df, manifest_file, row.names = FALSE)

cat("RQ3 brms server runner finished.\n")
cat("Manifest:", manifest_file, "\n")
