# Trait Occurrence Scripts

## Recommended mainline

- `03_fg_stable_glmmtmb.R`
- `04_fg_stable_brms.R`

These two scripts represent the most stable and interpretable functional-group workflow in the current release.

## Secondary / diagnostic scripts

- `02_multitrait_screening_occurrence.R`

This script is useful for screening and robustness checks, but it should not replace the main three-way scientific test on its own.

## Archival / legacy-path scripts

- `01_single_trait_threeway_occurrence.R`

This script is retained for transparency because it documents an important earlier workflow, but it still contains historical path dependence and version inconsistency noted in the diagnosis documents under `results/traits_fg`.

Please read:

- `results/traits_fg/RQ3_trait_analysis_diagnosis_CN.md`
- `results/traits_fg/RQ3_trait_analysis_diagnosis_EN.md`

before treating the older trait scripts as reproducible public mainline analyses.
