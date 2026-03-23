# Diagnostic Checklist for Scientific Question 3 (English)

## Scope

This checklist evaluates the current workspace analyses related to Scientific Question 3 (RQ3: traits/functional strategies as winners and losers under land-use x warming), focusing on three analysis lines:

1. Single-trait three-way interaction models
2. Multi-trait joint two-way interaction models
3. Functional-group (FG) models based on PCA/clustering

The assessment is based on the proposal framing of RQ3 and the current code/results in the workspace, especially:

- `/Users/dingchenchen/lucc/性状交互作用_出现概率.R`
- `/Users/dingchenchen/lucc/性状预测.R`
- `/Users/dingchenchen/lucc/筛选性状交互作用_出现概率.R`
- `/Users/dingchenchen/lucc/性状交互作用_FG_Gower1.R`
- Result folders:
  - `/Users/dingchenchen/lucc/02_ThreeWay_ABS_and_PCT1`
  - `/Users/dingchenchen/lucc/glmmTMB_occurrence_alltraits_AIC_2way_ABS_PCT_medianControls`
  - `/Users/dingchenchen/lucc/glmmTMB_occurrence_FG_GowerClustering_FULL_plusTwoWay`
  - `/Users/dingchenchen/lucc/glmmTMB_occurrence_FG_GowerWinsor_FULL_plusTwoWay`
  - `/Users/dingchenchen/lucc/brms_occurrence_FG3_PCA_rawTraits_FULL_plusTwoWay`

## 1. Overall conclusion

The RQ3 workflow is already rich and ambitious, but the three methods should not be treated as equivalent.

- Single-trait three-way interaction models are the closest match to the original RQ3 question and should be treated as the primary analysis.
- The multi-trait joint model is better interpreted as a conditional-independence or robustness analysis, not as the main RQ3 test.
- FG analysis is valuable for mechanism synthesis, but only stable FG schemes should be used as headline results. The current raw Gower two-group version should not be used as a main result.

## 2. Recommended for main results

### 2.1 Single-trait three-way interaction models

Why they should be primary:

- They directly match the proposal-level RQ3 formulation: `Occurrence ~ UI2 x warming x trait`
- They answer the biologically important question most clearly: which traits mediate vulnerability or resilience under interacting land use and warming
- They are straightforward to map onto the stated hypotheses

Locally verified example:

- `/Users/dingchenchen/lucc/brm_OccmodelRSTmean_ub.rds`

The verified model formula is:

- `Occur ~ UI2 * StdTmeanAnomalyRS * RS.rs + (1 || SS) + (1 || SSBS) + (1 || Best_guess_binomial)`

This RS model uses:

- 444,921 records
- 2,958 species

Strengths:

- It uses the full species pool rather than the reduced complete-case subset
- brms produces posterior summaries that are easy to report as 95% credible intervals
- It aligns best with the core RQ3 narrative

### 2.2 PCA-based three-group FG model

Recommended folder:

- `/Users/dingchenchen/lucc/brms_occurrence_FG3_PCA_rawTraits_FULL_plusTwoWay`

Why this FG scheme is defensible:

- Group sizes are reasonably balanced at the species level: FG1 = 1898, FG2 = 514, FG3 = 433
- Record-level proportions are also acceptable:
  - FG1 = 216,947
  - FG2 = 66,130
  - FG3 = 149,483
- Trait syndromes are interpretable as ecological strategy groups rather than outlier partitions

Key interpretation file:

- `/Users/dingchenchen/lucc/brms_occurrence_FG3_PCA_rawTraits_FULL_plusTwoWay/05_Tables_FG_Interpretation/FG_trait_signatures_and_interpretation.txt`

Current ecological reading:

- FG1: low-dispersal, lower thermal breadth, smaller-ranged, more specialist-like
- FG2: higher dispersal, longer generation length, somewhat larger-bodied
- FG3: wider-ranging, thermally broader, higher reproductive output

This FG structure is much more credible than the raw Gower two-group solution.

## 3. Recommended as supplementary or robustness analyses

### 3.1 Multi-trait joint two-way model

Script:

- `/Users/dingchenchen/lucc/筛选性状交互作用_出现概率.R`

Results:

- `/Users/dingchenchen/lucc/glmmTMB_occurrence_alltraits_AIC_2way_ABS_PCT_medianControls`

The full formula in the script is:

- `Occur ~ UI2 * StdTmeanAnomalyRS + UI2 * (all traits) + StdTmeanAnomalyRS * (all traits) + random effects`

Key limitation:

- This model does not include `UI2 x warming x trait`
- Therefore it does not directly test the core RQ3 trait-mediation hypothesis

What it is good for:

- Checking whether trait associations remain after controlling for other traits
- Assessing which effects are robust in a multivariate setting
- Providing a conditional-independence screen

What it is not good for:

- Serving as the main RQ3 inferential test
- Replacing the single-trait three-way framework

Recommended framing:

- Supplementary robustness analysis
- A secondary multivariate model used to assess whether focal trait signals remain after controlling for other traits

### 3.2 Winsorized Gower two-group FG model

Folder:

- `/Users/dingchenchen/lucc/glmmTMB_occurrence_FG_GowerWinsor_FULL_plusTwoWay`

Why it is useful:

- Extreme values were winsorized
- A minimum cluster size rule was imposed
- The resulting groups are much more stable than the raw Gower two-group version

Current species counts:

- FG1 = 2554
- FG2 = 291

Record-level counts:

- FG1 = 386,964
- FG2 = 45,596

Selection file:

- `/Users/dingchenchen/lucc/glmmTMB_occurrence_FG_GowerWinsor_FULL_plusTwoWay/05_Tables_FG_Interpretation/FG_GowerWinsor_method_K_silhouette_min5pct.csv`

This is a useful sensitivity analysis showing that once outliers are controlled and tiny clusters are disallowed, the FG contrast becomes much more stable.

Still, because it remains a two-group solution, it is less ecologically informative than the PCA-based three-group model.

## 4. Not recommended as a main result

### 4.1 Raw Gower two-group FG model

Folder:

- `/Users/dingchenchen/lucc/glmmTMB_occurrence_FG_GowerClustering_FULL_plusTwoWay`

Selection log:

- `/Users/dingchenchen/lucc/glmmTMB_occurrence_FG_GowerClustering_FULL_plusTwoWay/05_Tables_FG_Interpretation/FG_trait_signatures_and_selection_log.txt`

The main problem is cluster collapse:

- FG1 = 2835 species
- FG2 = 10 species

Record-level imbalance is also extreme:

- FG1 = 427,688
- FG2 = 4,872

The FG2 species include a small set of highly distinctive taxa such as:

- `Egretta garzetta`
- `Falco peregrinus`
- `Himantopus himantopus`
- `Milvus migrans`
- `Pandion haliaetus`
- `Larus argentatus`

This strongly suggests that FG2 is functioning as an outlier cluster rather than a stable ecological strategy group.

The code-level reason is visible in:

- `/Users/dingchenchen/lucc/性状交互作用_FG_Gower1.R`

The clustering is selected by silhouette only, without a minimum cluster-size constraint.

Implications:

- The “best” silhouette solution is not necessarily ecologically meaningful
- A few extreme species can be peeled off into an artificial cluster
- Downstream three-way interaction significance may partly reflect outlier-vs-background separation rather than robust FG ecology

Conclusion:

- The raw Gower two-group FG result should not be used as a headline result
- It can only be mentioned briefly as an exploratory attempt or sensitivity check

## 5. Key code-result inconsistencies

### 5.1 Single-trait plotting script still depends on external paths

In:

- `/Users/dingchenchen/lucc/性状交互作用_出现概率.R`

around the model-loading section, the script still reads:

- `D:/DCC/BIRDLIFE1/...`

This means the current plotting script is not a fully self-contained, workspace-reproducible final version.

### 5.2 Tmin/Tmax position scaling is inconsistent across scripts

In the plotting script:

- `/Users/dingchenchen/lucc/性状交互作用_出现概率.R`

the comments state that `Tmin_position` and `Tmax_position` are used on their raw 0-1 scales.

But in the modeling script:

- `/Users/dingchenchen/lucc/性状预测.R`

the actual fitted models use `Tmin_position.rs` and `Tmax_position.rs`.

This creates a real interpretational risk:

- figure labels and model scales may not match
- any manuscript wording about these trait definitions could become inconsistent

### 5.3 Some trait figures come from mixed code generations

For example:

- `/Users/dingchenchen/lucc/02_ThreeWay_ABS_and_PCT1` contains BM plots
- but the current plotting script no longer includes BM in its active trait loop

This suggests that at least some exported figures were generated from earlier script versions and are not fully synchronized with the current code.

## 6. Cross-method comparability is limited by different sample pools

Different methods are using different subsets of records and species:

- Single-trait RS model: 444,921 records, 2,958 species
- Multi-trait joint model: 432,412 records, 2,845 species
- FG raw-trait model: 432,560 records, 2,845 species

The main losses are driven by missing values in:

- `TR`
- `Tmin_position`
- `Tmax_position`

Therefore:

- effect strengths should not be compared naively across methods
- changes in the species pool are confounded with changes in model structure

Recommended fix:

- either standardize all main analyses to the same 2,845-species complete-case dataset
- or explicitly state that each analytical layer uses a different species pool

## 7. Percent-change plots are not directly comparable across methods

The `%Δp` baselines differ across analytical frameworks:

- Single-trait models: PV + temp = 0 + trait Q1
- Multi-trait model: PV + temp = 0 + all traits at median
- FG model: PV + temp = 0 + FG base group

So:

- `%Δp` can be interpreted within each analysis
- `%Δp` should not be compared across methods
- if cross-method synthesis is needed, absolute probabilities, contrasts, or marginal slopes are safer

## 8. Uncertainty intervals are not equivalent across model families

For the brms single-trait models:

- uncertainty is based on posterior draws and can be treated as proper posterior credible intervals

For the glmmTMB multi-trait and FG models:

- uncertainty bands are currently based on fixed-effect covariance draws or link-scale normal approximation

Relevant code locations:

- multi-trait model: `/Users/dingchenchen/lucc/筛选性状交互作用_出现概率.R`
- FG glmmTMB model: `/Users/dingchenchen/lucc/性状交互作用_FG_Gower1.R`

This means:

- the glmmTMB intervals are fine for exploratory visualization
- but they should not be described as fully comparable to brms posterior intervals

## 9. Final ranking of methods

### Recommended for main text

1. Single-trait three-way brms models
2. PCA-based three-group FG model

### Recommended for supplements / robustness

1. Multi-trait joint two-way glmmTMB model
2. Winsorized Gower two-group FG model

### Not recommended as main evidence

1. Raw Gower two-group FG model
2. Older figures whose generating code version cannot be cleanly traced

## 10. Most useful next steps

If the goal is a manuscript-ready or presentation-ready RQ3 section, the most valuable next steps are:

1. Unify the single-trait scripts, model objects, and exported figures into one reproducible version.
2. Decide whether `Tmin_position/Tmax_position` are raw or standardized, and keep that decision consistent everywhere.
3. Reframe the multi-trait joint model explicitly as a supplementary robustness analysis.
4. Use PCA-based FG3 as the main FG result, winsorized Gower FG2 as sensitivity, and drop raw Gower FG2 from the main narrative.
5. If strict cross-method comparison is desired, rerun all key models on the same complete-case species pool.

## One-sentence summary

The most defensible RQ3 storyline is:

"Use single-trait three-way interaction models to test trait-mediated vulnerability directly, use a stable FG scheme to synthesize mechanisms, treat the multivariate joint model as a robustness check, and avoid using the raw Gower two-group result as main evidence."
