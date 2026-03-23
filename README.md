# LUCC Biodiversity Reproducible Release

This repository is a curated, GitHub-ready release of the larger `LUCC` working directory.
It is designed to make the project understandable and reproducible for external users without exposing the entire exploratory workspace.

## Scope

This release reorganises the main LUCC tasks into five layers:

1. Environment and setup
2. Data preparation and spatial matching
3. Main land-use × climate biodiversity models
4. Community diversity and community structure analyses
5. Trait-based occurrence and functional-group analyses

It includes:

- Canonical scripts for the main analysis pipeline
- A curated set of derived `.rds` files small enough to distribute directly
- Key result previews and summary tables
- Clear manifests for tasks, data, and results
- Explicit notes on what is canonical, what is supplementary, and what remains archival

It does **not** attempt to track every exploratory or legacy file from the original workspace.
It also keeps the public Git history lightweight by omitting some large fitted-model binaries when equivalent fixed-effect tables or predictions are already included.

## Repository layout

- [`scripts`](./scripts): analysis scripts organised by stage
- [`data`](./data): directly distributed small derived data + external-data instructions
- [`results`](./results): key result previews, rewritten outputs, and diagnosis documents
- [`docs`](./docs): project overview, run order, and reproducibility notes
- [`manifests`](./manifests): machine-readable inventories for tasks, data, and results

## Quick start

1. Read [`docs/PIPELINE_RUN_ORDER.md`](./docs/PIPELINE_RUN_ORDER.md)
2. Check [`manifests/data_manifest.csv`](./manifests/data_manifest.csv)
3. Run [`scripts/00_setup/02_check_local_inputs.R`](./scripts/00_setup/02_check_local_inputs.R)
4. Start from:
   - `scripts/01_data_prep` for full rebuilds
   - `scripts/04_community` for the cleaned community-analysis release
   - `scripts/05_traits_occurrence/03_fg_stable_glmmtmb.R` for the stable FG main analysis

## Canonical recommendations

### Main land-use × climate biodiversity models

- Use [`scripts/03_main_models/01_run_lu_climate_models_main.R`](./scripts/03_main_models/01_run_lu_climate_models_main.R) as the main non-Bayesian diversity workflow.
- Treat the brms version as a robustness layer rather than the default first-pass model family.

### Community analyses

Use the rewritten community release in [`scripts/04_community`](./scripts/04_community):

- `01_diversity_main_rewrite.R`
- `02_compositional_main_rewrite.R`
- `03_functional_main_rewrite.R`

These scripts correspond to the cleaned outputs in [`results/community`](./results/community).

### Trait-based occurrence analyses

- Stable functional-group main analysis:
  [`scripts/05_traits_occurrence/03_fg_stable_glmmtmb.R`](./scripts/05_traits_occurrence/03_fg_stable_glmmtmb.R)
- Bayesian stable FG extension:
  [`scripts/05_traits_occurrence/04_fg_stable_brms.R`](./scripts/05_traits_occurrence/04_fg_stable_brms.R)

The older single-trait and multi-trait scripts are kept for transparency but should be interpreted together with the diagnosis documents in [`results/traits_fg`](./results/traits_fg).

## Data policy

Small derived data objects used repeatedly across the cleaned analyses are included directly in [`data/derived_public`](./data/derived_public).

Large raw or generated data are **not** all versioned in GitHub. For those, see:

- [`data/raw_external/README.md`](./data/raw_external/README.md)
- [`manifests/data_manifest.csv`](./manifests/data_manifest.csv)

## Reproducibility notes

- Some scripts from the original workspace remain archival and are not fully portable.
- The cleaned community-analysis scripts in `scripts/04_community` were edited to remove machine-specific absolute paths.
- Large Bayesian fit objects are not tracked when they exceed practical GitHub size limits.
- Some large non-Bayesian preview model objects are also omitted from the GitHub release when the same figure can be rebuilt from included coefficient tables and prediction `.rds` files.

For the full status map, see [`docs/LEGACY_STATUS.md`](./docs/LEGACY_STATUS.md).
