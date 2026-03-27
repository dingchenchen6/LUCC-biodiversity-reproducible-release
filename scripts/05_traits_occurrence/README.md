# Trait Occurrence Scripts

## Recommended Bayesian mainline

- `01_install_rq3_brms_server_dependencies.R`
- `00_run_rq3_brms_server.R`
- `run_rq3_brms_pipeline_nohup.sh`
- `02_multitrait_joint_occurrence_brms.R`
- `04_fg_stable_brms.R`

This is the recommended server-oriented Bayesian workflow for RQ3 in the current release.

## Modules

- `01_single_trait_threeway_occurrence.R`
  This is the original single-trait brms workflow kept for transparency and direct trait-by-trait three-way tests.
- `02_multitrait_joint_occurrence_brms.R`
  This is the clean Bayesian replacement for the earlier glmmTMB multitrait screening script.
- `04_fg_stable_brms.R`
  This is the Bayesian stable-FG mainline with server-friendly defaults, posterior diagnostics, and paper-grade plots.

## Server background run

Recommended background launch:

```bash
cd /path/to/LUCC_reproducible_release_publish
Rscript scripts/05_traits_occurrence/01_install_rq3_brms_server_dependencies.R
bash scripts/05_traits_occurrence/run_rq3_brms_pipeline_nohup.sh
```

Useful environment overrides:

```bash
export RQ3_RUN_SINGLE=1
export RQ3_RUN_MULTITRAIT=1
export RQ3_RUN_FG=1

export MT_BRMS_ITER=6000
export MT_BRMS_WARMUP=3000
export MT_BRMS_THREADS=2
export MT_BRMS_ADAPT_DELTA=0.99
export MT_BRMS_MAX_TREEDEPTH=15

export FG_BRMS_ITER=6000
export FG_BRMS_WARMUP=3000
export FG_BRMS_THREADS=2
export FG_BRMS_ADAPT_DELTA=0.99
export FG_BRMS_MAX_TREEDEPTH=15
```

Logs and manifest files are written under:

- `results/traits_fg/rq3_brms_server_runner`

## Legacy / archival scripts

- `02_multitrait_screening_occurrence.R`

This older glmmTMB + AIC script is retained for transparency but is no longer the recommended public mainline for Bayesian reproducibility.

Please read:

- `results/traits_fg/RQ3_trait_analysis_diagnosis_CN.md`
- `results/traits_fg/RQ3_trait_analysis_diagnosis_EN.md`

before treating the older trait scripts as reproducible public mainline analyses.
