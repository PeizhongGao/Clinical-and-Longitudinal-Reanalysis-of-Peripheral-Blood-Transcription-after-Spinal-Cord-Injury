# Archived source-language records

The README, user instructions and current validation documentation are in English.
Author-facing Chinese upload guides and internal submission-management notes are
kept outside this public repository. Some historical source records retain their
original language to preserve their content, timestamps and checksums.

These records do not need to be removed to use the English reproduction workflow.
Deleting a locked protocol can prevent a run or invalidate its source checks.
The table below provides English navigation; it is not a replacement protocol,
a new preregistration or a claim that every archived paragraph has been translated.

All paths below are relative to the repository root.

| Source record | Role and English interpretation |
|---|---|
| `inputs/project/revision_scope_lock_20261018/bridge43/analysis_protocol_zh.md` | Locked October 18 revision-stage protocol for the same 43 participants, 17,498 genes and 4,032 GO sets under M0/M1/M2. The runner records its hash. Its structured companion `analysis_spec.json` specifies the models and parameters largely in English. |
| `inputs/project/revision_scope_lock_20261018/bridge43/analysis_spec.json` | Structured bridge specification; an original historical source-path field contains Chinese text. Preserve the source identifier and portability lock. |
| `inputs/project/revision_integrated_20261012/fixed_sample/analysis-protocol.md` | October 12 fixed-sample score models: the existing 48/43/41-person sets, frozen 15+15 scores, OLS/HC3 inference, all 27 score models and sensitivity diagnostics. |
| `inputs/project/revision_R2_20261008/fixed_sample_20261012/analysis-protocol.md` | Corresponding staged copy of the fixed-sample protocol, used in the executable input/provenance chain. |
| `inputs/project/revision_external_adjusted_20261016/analysis-protocol.md` | October 16 revision-stage external-module plan: four BTM modules, two definitions, 88 clinical coefficients, 240 longitudinal contrasts and the designated 16-test primary family. External definitions are not independent participants. |
| `inputs/project/revision_R1_20261007/tables/author_provenance_confirmation.json` | Original confirmation that contemporaneous candidate-pool/threshold-formation records were not retained. Existing code and lists do not establish prospective definition or recover missing historical justification. |
| `reference_results/supplementary_data/Supplementary_Data_7_Multiplicity_and_disposition/Claim_evidence_boundaries.json` | Frozen interpretation boundaries: model sensitivity, the fixed GO panel, imprecise clinical BTM estimates, limited longitudinal associations, participant reuse and withdrawn predictive/postselection claims. |
| `historical/auxiliary_scripts/revision_R1_20261007/scripts/lock_external_module_protocol.py` | Historical source that generated a Chinese protocol. It is retained for provenance and is not a current execution entry point. |

## Current scope and access information

The complete negative and uncertain results remain in the package. Preserving an
archived record does not reinstate withdrawn LASSO/prediction or postselection
permutation claims. The 15+15 scores remain exploratory; BTM uses external gene
definitions in the recorded participant sets. See README.md and VALIDATION.md for
the current inferential and reproduction boundaries.

Historical statements that treated the published ODC405 dataset as restricted are
superseded by DATA_PROVENANCE.md and the retained public source/license metadata.
Protocol status statements refer to when the records were written; the actual
completed runs are documented in validation/ and VALIDATION.md.
