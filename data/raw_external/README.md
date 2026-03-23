# External Raw Data

This folder documents raw or large upstream data that are not fully tracked inside the GitHub release.

## Public sources already referenced in the project

### PREDICTS

- 2016 release:
  <https://data.nhm.ac.uk/dataset/the-2016-release-of-the-predicts-database>
- November 2022 additional release:
  <https://data.nhm.ac.uk/dataset/release-of-data-added-to-the-predicts-database-november-2022>

The project also contains direct `readRDS(url(...))` calls in:

- `scripts/01_data_prep/01_prepare_predicts_data.R`

### CRU climate data

- <https://crudata.uea.ac.uk/cru/data/hrg/>

The original setup notes request the monthly `tmp` and `tmx` products to be placed in a local `0_data` directory before climate-layer preparation.

### Natural habitat layers

- Hoskins et al. layer reference:
  <http://doi.org/10.4225/08/56DCD9249B224>

## Local-only large generated products

Some large outputs are intentionally not tracked in GitHub, for example:

- full `0_data/`
- heavy `brms` fit objects
- large historical result directories from exploratory workflows

See `manifests/data_manifest.csv` for the exact status of each important dataset.
