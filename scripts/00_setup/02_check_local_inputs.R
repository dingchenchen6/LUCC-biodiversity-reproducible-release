suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

script_dir <- normalizePath(dirname(sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grepl("^--file=", commandArgs(trailingOnly = FALSE))][1])))
project_root <- normalizePath(file.path(script_dir, "..", ".."))

manifest <- readr::read_csv(file.path(project_root, "manifests", "data_manifest.csv"), show_col_types = FALSE)

check_path <- function(rel_path) {
  abs_path <- file.path(project_root, rel_path)
  data.frame(
    release_path_or_source = rel_path,
    exists = file.exists(abs_path),
    stringsAsFactors = FALSE
  )
}

included <- manifest %>%
  filter(status == "included") %>%
  mutate(check = lapply(release_path_or_source, check_path))

included_checks <- bind_rows(included$check)
print(included_checks)

message("\nSummary:")
message("Included files checked: ", nrow(included_checks))
message("Missing included files: ", sum(!included_checks$exists))

if (any(!included_checks$exists)) {
  stop("Some tracked derived data files are missing from the release.")
}

message("\nExternal datasets are documented but not validated automatically here.")
