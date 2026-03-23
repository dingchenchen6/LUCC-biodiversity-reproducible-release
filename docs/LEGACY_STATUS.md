# Legacy Status

## Canonical in this release

- `scripts/03_main_models/01_run_lu_climate_models_main.R`
- `scripts/04_community/01_diversity_main_rewrite.R`
- `scripts/04_community/02_compositional_main_rewrite.R`
- `scripts/04_community/03_functional_main_rewrite.R`
- `scripts/05_traits_occurrence/03_fg_stable_glmmtmb.R`
- `scripts/05_traits_occurrence/04_fg_stable_brms.R`

## Public release packaging rules

- Small recurrent derived data are included directly.
- Key figure previews and summary tables are included directly.
- Large fitted model binaries may be omitted when the same public-facing figure can be regenerated from tracked prediction tables or fixed-effect summaries.

## Included but not recommended as the default public mainline

- `scripts/03_main_models/02_run_lu_climate_models_brms.R`
  Reason: robustness layer, computationally heavy.

- `scripts/03_main_models/03_run_tropical_extension.R`
  Reason: extension analysis, not the main global pipeline.

- `scripts/05_traits_occurrence/01_single_trait_threeway_occurrence.R`
  Reason: historical single-trait workflow with legacy path dependence and version inconsistency.

- `scripts/05_traits_occurrence/02_multitrait_screening_occurrence.R`
  Reason: useful screening layer, but not the cleanest direct test of the focal three-way RQ.

## Archived in the original workspace, not carried forward as main release scripts

Examples include:

- old duplicate scripts with `_副本`
- exploratory similarity workflow variants
- early or superseded FG clustering scripts
- large intermediate result folders created during method development

These remain in the original `LUCC` workspace but are intentionally not promoted as the GitHub release mainline.
