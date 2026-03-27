# RQ3 Bayesian Server Run Guide

This note documents the recommended server-side workflow for the third scientific question (RQ3) under the reproducible release.

## Goal

Run the trait-mediated occurrence analyses as a coherent Bayesian pipeline that is:

- server friendly
- restartable
- logged
- traceable via manifests
- reproducible from GitHub

## Main scripts

- `scripts/05_traits_occurrence/01_single_trait_threeway_occurrence.R`
- `scripts/05_traits_occurrence/02_multitrait_joint_occurrence_brms.R`
- `scripts/05_traits_occurrence/04_fg_stable_brms.R`
- `scripts/05_traits_occurrence/00_run_rq3_brms_server.R`
- `scripts/05_traits_occurrence/run_rq3_brms_pipeline_nohup.sh`

## Recommended background launch

```bash
cd LUCC_reproducible_release_publish
bash scripts/05_traits_occurrence/run_rq3_brms_pipeline_nohup.sh
```

## Recommended robust sampling settings

These defaults are already encoded in the scripts, but can be overridden from the shell:

```bash
export MT_BRMS_CHAINS=4
export MT_BRMS_ITER=6000
export MT_BRMS_WARMUP=3000
export MT_BRMS_THREADS=2
export MT_BRMS_ADAPT_DELTA=0.99
export MT_BRMS_MAX_TREEDEPTH=15

export FG_BRMS_CHAINS=4
export FG_BRMS_ITER=6000
export FG_BRMS_WARMUP=3000
export FG_BRMS_THREADS=2
export FG_BRMS_ADAPT_DELTA=0.99
export FG_BRMS_MAX_TREEDEPTH=15
```

## Selectively running modules

```bash
export RQ3_RUN_SINGLE=1
export RQ3_RUN_MULTITRAIT=1
export RQ3_RUN_FG=1
```

Set any of these to `0` if a module should be skipped.

## Smoke-test mode

For quick checks before a full server run:

```bash
export MT_BRMS_SMOKE_N=5000
export FG_BRMS_SMOKE_N=5000
```

This reduces data volume for a short validation run without changing the full-run code path.

## Outputs

The Bayesian runner writes log files and a manifest to:

- `results/traits_fg/rq3_brms_server_runner`

The multitrait Bayesian outputs are written to:

- `brms_occurrence_alltraits_joint2way_medianControls`

The stable FG Bayesian outputs are written to:

- `brms_occurrence_FG_RobustPCAKmeans_FULL_plusTwoWay`

## Reproducibility notes

- The multitrait Bayesian script is intended as a conditional-effect robustness screen, not as the sole headline RQ3 test.
- The stable FG Bayesian script is the recommended Bayesian mainline for the functional-group synthesis part of RQ3.
- The older glmmTMB multitrait script is kept for transparency but is not the preferred public Bayesian workflow.
