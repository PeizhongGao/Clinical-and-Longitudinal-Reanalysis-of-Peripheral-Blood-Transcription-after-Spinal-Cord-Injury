# Analysis report

## Comparison and units

This exploratory revision audit evaluates a newly rebuilt candidate-union LASSO, not the historical fixed nine-gene signature. Each classification comparison contains 38 SCI and 10 controls (48 independent participants). The four observed source/contrast configurations each use 20 repeats of five outer and three inner folds. Training-only screening, scaling and tuning precede held-out prediction. CV repeats and folds are dependent computational repeats.

## Complete numeric summary

| Source | Contrast | Mean repeat AUC ± SD | Repeat range | Mean Brier | Intercept-only fits |
| --- | --- | --- | --- | --- | --- |
| Raw counts | SCI–HC | 0.979 ± 0.015 | 0.926–1.000 | 0.053 | 0/100 |
| Raw counts | SCI–TC | 0.864 ± 0.068 | 0.747–0.958 | 0.115 | 17/100 |
| Archived processed | SCI–HC | 0.998 ± 0.003 | 0.989–1.000 | 0.023 | 0/100 |
| Archived processed | SCI–TC | 0.868 ± 0.072 | 0.668–0.958 | 0.121 | 9/100 |

The raw-count analysis is primary. Archived processed-source results remain conditional on unreconstructed upstream whole-cohort processing. AUC split SD is not a population confidence interval. Conditional participant-bootstrap intervals for fixed averaged OOF predictions are provided separately; those ensemble-prediction metrics are not substituted for the primary mean-repeat AUC.

## Full-pipeline permutation evidence

| Contrast | First-repeat OOF AUC | Permutations | Exceedances | Permutation P | Holm2 |
| --- | --- | --- | --- | --- | --- |
| SCI_HC | 0.974 | 999 | 0 | 0.001 | 0.002 |
| SCI_TC | 0.921 | 999 | 0 | 0.001 | 0.002 |

The observed and null comparators use the same first nested-CV schedule. Each of 999 permutations reruns the entire newly specified supervised pipeline with labels permuted within sex. The null is marginal class association conditional on sex; age, injury severity, batch and treatment are not controlled. The fixed historical heuristic candidate recipe remains data-informed in origin. The P values do not correct unrecorded researcher choices or supply independent validation.

## Clinical resampling evidence

| BTM | SCI–TC estimate (original 95% CI) | Residual-permutation P | Wild-bootstrap P | Separate Holm4 |
| --- | --- | --- | --- | --- |
| M37.1 | -0.163 (-0.651, 0.326) | 0.4374 | 0.4489 | 1.000 / 1.000 |
| M163 | 0.019 (-0.287, 0.325) | 0.8783 | 0.8875 | 1.000 / 1.000 |
| M7.0 | 0.144 (-0.496, 0.783) | 0.5982 | 0.6006 | 1.000 / 1.000 |
| M223 | 0.242 (-0.367, 0.852) | 0.3555 | 0.3526 | 1.000 / 1.000 |

All eight Holm4 values equal 1. The original estimates/HC3 standard errors were reproduced before 9,999 residual permutations and 9,999 wild bootstrap draws. These two sensitivities do not provide evidence to strengthen the clinical SCI–TC claim. They do not establish equivalence or absence of association. The maximum full-model leverage is 0.701; only seven trauma controls contribute. Nonsignificant diagnostics do not prove residual exchangeability, normality or homoskedasticity.

## Claim Candidates

- Claim: The newly specified raw-count nested workflow shows internal class discrimination.
  - Source evidence: nested_performance_summary.csv; full_pipeline_permutation_results.csv.
  - Allowed wording: Exploratory, internally evaluated discrimination under the fixed heuristic recipe and stated sex-conditional null.
  - Forbidden stronger wording: Independently validated SCI diagnostic biomarker, recovery predictor, or clinically independent SCI-specific effect.
  - Uncertainty: Small control groups, historical heuristic provenance, missing batch metadata, split dependence and no independent evaluation cohort.
  - Next check: Review the supplementary methods and preserve all source/split/null records.
  - Decision: keep as supplementary methodological evidence.
- Claim: Selection is not an invariant reconstruction of the old nine-gene model.
  - Source evidence: nested_outer_fit_audit.csv; nested_gene_selection_frequencies.csv; nested_selection_overlap.csv.
  - Allowed wording: The selected gene counts and identities vary across training partitions; SCI–TC has 17 intercept-only raw-count fits among 100 outer fits.
  - Forbidden stronger wording: Recovery of a fixed, stable original five-plus-four signature or formal false-selection control.
  - Uncertainty: The 100 fits and pairwise Jaccard values are dependent.
  - Next check: Keep the historical signature withdrawn.
  - Decision: revise historical stability claims.
- Claim: Alternative fixed-score resampling does not strengthen clinical BTM differences.
  - Source evidence: clinical_resampling_results.csv and diagnostics.
  - Allowed wording: Both additional resampling sensitivities leave the four adjusted clinical estimates uncertain.
  - Forbidden stronger wording: Equal groups, proof of no biological effect, or rescue of the historical selected-axis permutation.
  - Uncertainty: Fitted nuisance models, approximate exchangeability, small TC group and high leverage.
  - Next check: Retain the existing cautious clinical interpretation.
  - Decision: keep boundary.
