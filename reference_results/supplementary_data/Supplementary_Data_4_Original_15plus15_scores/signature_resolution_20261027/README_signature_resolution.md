# Frozen discovery score provenance and reconciliation

This 27 October 2026 audit reconciles archived inputs and saved outputs. It introduces no new expression analysis, score selection, participant comparison or significance test. The frozen 30-gene list and all reported numerical results are retained.

## What can be verified

The byte-preserved 4_archived_source.Rmd contains 59 myeloid/innate-associated and 47 lymphocyte-associated candidates. The locally available GSE151371_stability_selection_all_candidates.csv has 103 rows (59 and 44). Candidate_to_frozen_signature_audit.csv reports every candidate, marks the three unrepresented entries as unavailable and never assigns them zero frequencies. All 30 saved genes are present in the appropriate pool and in the archived frequency table.

The workflow draws 500 stratified 80% subsamples without replacement. Each iteration requires concordant SCI–TC and SCI–HC directions, ranks eligible genes by nominal SCI–TC P values and selects at most 15 per component. Retention uses frequency >=0.40, directional consistency >=0.80 and full-data direction concordance. The archived fallback applies below five genes and the final maximum is 20 per component. Applying these retention rules to the archived frequency table reproduces the saved 15+15 membership and order; neither size branch changes it. Therefore 15+15 describes the saved outcome, not a forced final size. These are repeated subsampling-based ranking rules, not formal false-selection control.

The two component scores are means of control-referenced gene z scores; their difference has coefficients +1 and -1. Equal numbers of genes are not required by that formula and are not presented as the historical reason for choosing the iteration cap. Existing threshold and top-N sensitivity files remain local, conditional checks; they do not identify optimal settings or prove convergence.

## Historical boundary

The source, JSON and frequency table verify their present archived contents. They do not independently establish the first candidate-formation date, original gene-by-gene literature sources, first inspection of GSE151371 results, or the historical rationale for 15. File hashes, self-entered dates, random seeds and script line order are not evidence of pre-data specification. Later annotation is explicitly contemporary interpretation. The original code is preserved, including legacy terminology. The separately labelled review copy changes comments only and was not executed. Legacy JSON keys such as validation_dataset are preserved source fields, not current claims of independent validation.

## Relationship to BTM and other analyses

The published BTM definitions, frozen discovery score and nested classifier have distinct roles. BTM members were not selected from the frozen 30-gene list; this is external gene-definition evidence in the existing cohorts, not independent participant confirmation. The nested LASSO in Data 9 evaluates a different training procedure and does not validate the fixed 30-gene score. GSE293559 represents 14/15 genes in each component, missing TLR5 and ITK, and includes documented participant reuse.

BTM_proxy_exclusion_dependency_audit.csv compares the original MCP-counter marker sets under BTM-only and BTM-plus-original30 exclusions. No extra marker is excluded by the original30 in either cell marker set. Thus the marker membership used in the existing BTM proxy analysis is unchanged by removing this redundant original30 term. The current coefficients, uncertainty intervals and multiplicity families are preserved; this membership audit is not a new robustness experiment.

## Files to inspect

- Candidate_to_frozen_signature_audit.csv: all 106 candidates and archived record coverage.
- Frozen_30_gene_selection_audit.csv: all frozen genes and recorded frequencies/direction consistency.
- Archived_threshold_reconciliation.csv: the complete 3 x 3 threshold grid for both components, computed from archived summaries.
- Archived_rule_line_map.json: exact lines of the preserved source.
- BTM_proxy_exclusion_dependency_audit.csv: all original marker entries and exclusion flags.
- Signature_resolution_audit.json and Source_evidence_manifest.csv: checks and content hashes.
- Source files and the comment-only review copy: inspectable historical material, not new analysis entry points.

Source-specific licenses remain applicable; this addition grants no new blanket license. This clarification is supplied with Supplementary Data 4 and has not been pushed to the previously cited GitHub snapshot.
