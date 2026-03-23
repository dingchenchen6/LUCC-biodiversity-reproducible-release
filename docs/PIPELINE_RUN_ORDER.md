# Pipeline Run Order

## Full rebuild

1. `scripts/00_setup/00_install_requirements.R`
2. `scripts/01_data_prep/01_prepare_predicts_data.R`
3. `scripts/01_data_prep/02_prepare_climate_index_maps.R`
4. `scripts/01_data_prep/03_match_predicts_climate_index.R`
5. `scripts/01_data_prep/04_match_predicts_prop_natural_habitat.R`
6. `scripts/03_main_models/01_run_lu_climate_models_main.R`

## Baseline and auxiliary model layers

1. `scripts/02_baselines/01_run_simple_lui_model.R`
2. `scripts/02_baselines/02_landuse_region_models_update.R`
3. `scripts/03_main_models/02_run_lu_climate_models_brms.R`
4. `scripts/03_main_models/03_run_tropical_extension.R`

## Cleaned community-analysis release

1. `scripts/04_community/01_diversity_main_rewrite.R`
2. `scripts/04_community/02_compositional_main_rewrite.R`
3. `scripts/04_community/03_functional_main_rewrite.R`

Note:
- In the public GitHub release, some `results/community/*/source_models/*.rds` files may be omitted to keep the repository lightweight.
- The compositional and functional rewrite scripts can still regenerate the coefficient figure from the included fixed-effect CSV tables.

## Trait occurrence and functional-group release

1. `scripts/05_traits_occurrence/03_fg_stable_glmmtmb.R`
2. `scripts/05_traits_occurrence/04_fg_stable_brms.R`

## Important interpretation note

The scripts in `scripts/05_traits_occurrence/01_*` and `02_*` are kept for transparency and diagnosis, but they should not be treated as the cleanest final release without reading:

- `results/traits_fg/RQ3_trait_analysis_diagnosis_CN.md`
- `results/traits_fg/RQ3_trait_analysis_diagnosis_EN.md`
- `docs/LEGACY_STATUS.md`
