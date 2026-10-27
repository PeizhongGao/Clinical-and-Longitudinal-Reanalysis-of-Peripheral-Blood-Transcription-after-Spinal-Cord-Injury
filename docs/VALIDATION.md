# Reproduction validation

Validation date: 26 October 2026. The standard workflow was executed from a copied
package in a temporary directory **outside the original SCI project**. Input, code
and output paths were resolved within that copied package/run directory. The original
project remained untouched. All 19 recorded execution stages completed.

All **21 validation groups passed**, covering **357,077 numerical values**,
with absolute tolerance 1e-8 plus relative tolerance 1e-8 times the reference magnitude.
The original two 15-gene membership tables were also checked by gene identity/module.
These are numerical reproducibility checks, not independent scientific validation.

## Recomputed and checked

- Discovery and full-count longitudinal 15+15 score projection and the 500-iteration
  revision-stage processed-matrix selection sensitivity.
- Full-count TMM normalization, 16 BTM score series and retained BTM gene memberships.
- ODC405 participant clinical extraction/mapping and 27 fixed-sample score models.
- All 88 BTM clinical models and 240 longitudinal contrasts, with their existing BH
  families, including the complete designated primary16 family.
- Original-score longitudinal and same-platform/reuse sensitivity results.
- 17,498 gene models and 4,032 CAMERA GO sets under each of M0/M1/M2.
- Main Figure 1 numerical source data and generation of the supporting figure files.

The separate `bridge_input_reconstruction.json` verifies that source expression,
clinical membership, model designs and GO indices reproduce the locked packet.
The standard workflow uses fixed input definitions/checkpoints where documented;
it is not a FASTQ-to-results pipeline.

## Explicit boundaries

- The original count-based 500-iteration signature selection is an optional mode.
  Its startup and early iterations were exercised, then the optional test was stopped
  after a bounded execution window. It was **not** counted as a completed validation.
  The delivered standard workflow uses the original frozen selection checkpoint.
- Historical clustering and historical preranked-GO permutations are included as
  optional `extended` steps, with the original seeds and iteration counts; they were
  not fully rerun for this packaging validation. Standard-mode Figure S4 uses saved
  clustering results. Extended dependencies are documented separately.
- No clean-machine dependency installation or fresh renv restore was certified.
  The recorded installed environment was used. Linux/macOS execution was not tested.
- Source-study expression preprocessing and supplied CIBERSORTx fractions were reused,
  not re-estimated. No sequence alignment or new functional/predictive validation was performed.
- All source versions, cohorts, effect interpretation and multiplicity boundaries
  remain those of the existing revision; reproducibility does not raise evidence level.

## Evidence

See `validation/relocated_standard_validation.json`, `validation/relocated_run_status.json`,
`validation/executed_source_match.csv`, the per-stage logs and `validation/SUMMARY.json`.
All delivered active scripts were checked against the script hashes recorded by the run.
The earlier in-project standard test is retained separately for audit, not counted as
independent biological evidence. Numerical reference values were never replaced.

## Upload bundle 1.0.0 verification

The 1.0.0 upload bundle was rerun on 26 October 2026 with outputs in a new temporary
folder outside the original project. All 19 scientific stages and
21 validation groups passed (357,077 compared numerical values).
All 10 named main/supplementary figures were exported. The runner now marks a run
complete only after export succeeds, and six source/parent output targets were
confirmed to be rejected before staging.

The executed code/input fingerprint matches the delivered bundle. All scientific
analysis scripts remain identical to those in the earlier relocated test; only the
execution wrapper and packaging tools were improved. Earlier test files retain
their actual version identifiers. Evidence is in validation/upload_v1_*.json,
validation/upload_v1_*.csv and validation/upload_v1_logs/.

All optional-analysis, fresh-installation and cross-platform limits stated above
continue to apply. Code completeness is separate from full optional test coverage.

## Extension validation (28 October 2026)

See `validation/revision_extension_20261028.json` for the relocated-package
verification, nested pilot and seven-figure rerender. These checks are separate
from the historical standard-pipeline validation. The complete 1,998 classifier
permutations are supplied from the executed analysis and checked from their saved
outputs; they were not all refitted during this synchronization.
