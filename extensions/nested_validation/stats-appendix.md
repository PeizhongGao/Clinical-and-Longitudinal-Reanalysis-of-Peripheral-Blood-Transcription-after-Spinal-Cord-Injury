# Statistical appendix

The immutable PROTOCOL.md and audit/protocol_input_lock.json define questions, inputs, seeds, preprocessing, models, tuning, null hypotheses, family definitions and boundaries before new outcome fits.

Independent units: participants. Each classifier comparison has 48; the clinical models have 43. No longitudinal samples enter these resampling analyses. The 20 CV repeats, 100 outer fits, 999 permutations and 9,999 clinical draws are not participant sample sizes.

Primary classification metric: mean of 20 repeat-level OOF AUCs, with split SD and range. Additional metrics: average precision, Brier score and log loss, with an outer-training intercept baseline. Five thousand stratified participant resamples of fixed mean OOF predictions produce explicitly conditional approximate intervals, excluding training uncertainty. Do not attach these intervals to the mean-repeat AUC.

Classifier permutation statistic: the AUC from the predeclared first five-fold/three-fold schedule, recomputed identically for each sex-stratified label permutation. P = (exceedances + 1)/(999 + 1). Holm2 is a new family comprising the two raw-count comparisons. Monte Carlo binomial intervals describe null-tail estimation uncertainty, not predictive-performance uncertainty.

Clinical statistic: absolute HC3-studentized SCI-minus-TC coefficient in the existing nine-column M2 design. Approximate Freedman–Lane permutations operate on nuisance residuals within sex after removal of the group column; common row maps preserve cross-module dependence. Wild bootstrap draws use common Rademacher weights on HC2-adjusted nuisance residuals. The full model and HC3 statistic are refitted for every draw. Each method has its own four-BTM Holm4 family. Existing BH16 results, pointwise confidence intervals and 93 original BH families are unchanged.

Failure handling: no resample is silently dropped. Empty selected models remain intercept-only. Expected glmnet small-class warnings are counted; unexpected warnings or nonzero convergence codes stop the run. Training-state invariance tests and independent Python/Numpy metric/HC3 reconstruction are recorded in audit/analysis_verification.json. Known limits include high leverage, small control groups, unavailable batch/treatment information, historical candidate-rule provenance and no external participant validation. Diagnostic P values cannot establish model assumptions.

Numerical verification found rank artifacts from approximately 1e-16 differences in intercept-only probabilities. The implementation now returns the training class prior exactly for an empty model, and participant predictions use 17-digit exports. All observed repeats and classifier permutations were regenerated under the unchanged locked protocol. The correction and prior local computation are documented in audit/numerical_implementation_correction.md.
