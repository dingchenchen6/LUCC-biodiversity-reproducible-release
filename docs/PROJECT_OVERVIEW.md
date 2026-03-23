# Project Overview

This release focuses on biodiversity responses to land use, climate anomalies, and functional strategy differences using bird data assembled from PREDICTS-linked workflows and trait tables.

## Main scientific components

1. Data assembly and cleaning
   - prepare bird-focused PREDICTS subsets
   - derive site-level biodiversity metrics

2. Climate and habitat covariates
   - build climate anomaly layers
   - match sites to climate and natural-habitat context

3. Main biodiversity models
   - land-use × climate effects on abundance and richness
   - regional and tropical extensions

4. Community structure
   - compositional dissimilarity
   - functional dissimilarity in ATA space

5. Trait occurrence and functional groups
   - single-trait and multi-trait occurrence modelling
   - stable functional-group derivation and interaction modelling

## Release philosophy

The original workspace contains a large amount of iterative development. This release therefore uses a “curated reproducibility” approach:

- keep canonical code paths
- keep interpretable key results
- keep small recurrent derived data
- document large-data dependencies instead of forcing all large files into GitHub
