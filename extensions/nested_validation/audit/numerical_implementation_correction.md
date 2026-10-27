# Intercept-only numerical correction

During independent reconstruction of AUC from exported participant predictions, five of the first 80 repeat-level AUCs disagreed with their in-memory R metrics. Inspection found approximately 1e-16 prediction differences within models classified as having no selected coefficients. Such floating-point differences create arbitrary ranks, and ordinary CSV precision can collapse those ranks. They do not represent meaningful participant discrimination.

The implementation was corrected to enforce the locked protocol's intercept-only behavior: any lambda path column whose predictor coefficients are all below the existing 1e-10 selection tolerance returns the exact outer or inner training class prior for every held-out participant. The same rule is applied inside tuning and outer prediction. Participant probabilities are exported with 17 significant digits; Python verification reads them with round-trip precision. No AUC-specific rounding, favorable-repeat selection, new candidate rule, lambda grid or significance stopping rule was introduced.

All 80 observed nested repeats and all 1,998 complete-pipeline label permutations were rerun from the same inputs, seeds and locked protocol after this correction. The first computation was preserved in the local internal_pre_numerical_fix directory for audit and is not used as the final scientific result. Clinical fixed-BTM resampling is unaffected because it does not use the classifier. Final completion, empty-model equality and metric reconstruction are checked in analysis_verification.json.

This is an implementation repair discovered during verification, not a retrospectively preregistered decision. The original protocol and input lock remain byte-for-byte unchanged.
