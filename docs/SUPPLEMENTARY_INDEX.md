# Supplementary Materials Index

Clinical and Longitudinal Reanalysis of Peripheral Blood Transcription after Spinal Cord Injury

This index identifies nine supplementary tables, eleven figures and nine data archives. Tables S1–S7, Figures S1–S9 and their original numerical results are retained. Tables S8–S9, Figures S10–S11 and Data 9 contain the newly executed nested/resampling audit. Data 7 includes CURRENT_REVISION_NOTE.md, which distinguishes historical withdrawal from the new extension.

## Public source and reproducibility record

The public code and supporting-data repository is https://github.com/PeizhongGao/Clinical-and-Longitudinal-Reanalysis-of-Peripheral-Blood-Transcription-after-Spinal-Cord-Injury. The verified analysis snapshot is 6415b0bdc3dff91d3949e940c363ed22f48f387c (package metadata version 1.0.0). Source data are from GEO GSE151371, GSE226238 and GSE293559; GSE177034 is metadata-only. Public ODC-SCI 405 data (DOI 10.34945/F5QC7J, CC BY 4.0) and the derived clinical crosswalk are available in the repository with source attribution. Data 8 gives exact paths and distinguishes the inventory from the actual public source code.

## Supplementary tables

| Item | File | Contents and boundary |

| --- | --- | --- |

| Supplementary Table S1 | Table_S1_Cohort_characteristics.csv | Cohort-level characteristics, denominators, control definitions, availability and exclusions. Aggregate fields only; no participant-level clinical linkage. |

| Supplementary Table S2 | Table_S2_Longitudinal_visit_flow.csv | Planned and observed longitudinal SCI observations by visit; three missing 12-month observations, with reasons unavailable. |

| Supplementary Table S3 | Table_S3_Fixed_immune_GO_panel.csv | All eight protocol-fixed immune GO terms under M0, M1 and M2, with complete-set BH correction across 4,032 sets per model. |

| Supplementary Table S4 | Table_S4_BTM_primary16.csv | Four published BTM clinical coefficients and twelve dual-disjoint-proxy longitudinal contrasts in the designated BH16 extension family. |

| Supplementary Table S5 | Table_S5_GSE293559_reuse_sensitivities.csv | All-score estimates for full, reused-SCI-excluded, same-platform and same-platform/reused-SCI-excluded scopes. |

| Supplementary Table S6 | Table_S6_Source_study_comparison.csv | Source-grounded separation of findings already reported in the original studies from the additional estimands and sensitivity analyses in this reanalysis. |

| Supplementary Table S7 | Table_S7_Multiplicity_family_registry.csv | The 93 retained named BH families, denominators, P/q fields and stated error-control boundaries. |

| Supplementary Table S8 | Table_S8_Nested_LASSO_audit.csv | All four observed source/contrast configurations, mean repeat AUC and split variation, baseline metrics, empty-model counts and the two raw-count complete-pipeline permutation tests. Conditional fixed-OOF intervals and every repeat are in Data 9. |

| Supplementary Table S9 | Table_S9_Clinical_BTM_resampling.csv | Four fixed clinical BTM effects with original HC3 intervals/BH16 and all residual-permutation/wild-bootstrap P and separate Holm4 values; 9,999 draws per method. |

## Supplementary figures

| Item | File | Caption and interpretation |

| --- | --- | --- |

| Supplementary Figure S1 | Figure_S1_Discovery_heatmap.png | Frozen discovery-derived gene heatmap. All 30 genes across 58 participants (38 SCI/10 TC/10 HC). Gene-wise z scores are display transformations. Selection used these group labels; apparent separation is descriptive. Original manuscript Figure 3. |

| Supplementary Figure S2 | Figure_S2_Original_score_longitudinal.png | Original score longitudinal observations and model uncertainty. Full-count 15+15 projection in 10 SCI participants/37 observations, with nine controls providing reference parameters. Participant connections and later-minus-acute model CIs are retained. Acute/3/6/12 months contribute 10/10/10/seven observations. Original manuscript Figure 5. |

| Supplementary Figure S3 | Figure_S3_Original_score_fixed_models.png | Original score sensitivity within fixed participant sets. SCI-minus-TC component coefficients with HC3 pointwise 95% CIs. Sets contain 48 (38/10), 43 (36/seven), or 41 (36/five) participants. Lists and reference parameters remain fixed. Postselection diagnostic evidence. Original manuscript Figure 8. |

| Supplementary Figure S4 | Figure_S4_Clustering_diagnostics.png | Conditional gene clustering diagnostics. Diagnostics cluster 2,509 previously selected gene profiles across 58 participants, rather than participant immune states. Silhouette, gap and bootstrap summaries do not consistently support k = 3. Resampling conditions on the selected gene universe. Original manuscript Figure 4. |

| Supplementary Figure S5 | Figure_S5_Clinical_support.png | Clinical covariate support in the discovery comparison. Observed age, injury severity and draw-time support across fixed SCI/TC subsets. The 43-person ISS intersection contains two SCI and three TC participants. Original fixed-sample covariate-support figure. |

| Supplementary Figure S6 | Figure_S6_BTM_composition_sensitivity.png | BTM longitudinal composition sensitivity across all visits. Published BTMs under time-only, official dual-proxy and disjoint dual-proxy models. Ten SCI participants/37 observations, REML random intercepts and pointwise 95% CIs. Full definitions and multiplicity families are in Supplementary Data 3. |

| Supplementary Figure S7 | Figure_S7_Additional_cohort_reuse.png | Additional cohort and partial reuse sensitivities. Original 36-SCI/nine-control and reused-SCI-excluded 26/nine comparisons retain all nine reused controls. Individual scores and uncertainty remain conditional on the supplied score transformation. Hedges g uses participant-bootstrap CIs; Wilcoxon P values test ranks. No trauma comparison is present. |

| Supplementary Figure S8 | Figure_S8_Fixed_immune_GO.png | Fixed immune GO panel across clinical models. All eight revision-fixed immune GO terms are shown under M0, M1 and M2 for the same 36 SCI/seven TC participants. Each q value is adjusted across the full 4,032-set family; no displayed term passed BH q < 0.05. Competitive-test direction is not pathway activity or a biological effect size. |

| Supplementary Figure S9 | Figure_S9_BTM_primary16_all.png | Complete clinical and longitudinal BTM family. All four clinically adjusted SCI-minus-TC coefficients and all 12 later-minus-acute contrasts in the designated BH16 family are shown with pointwise 95% CIs. Color identifies family-adjusted q < 0.05. Ten SCI participants contributed 37 observations; longitudinal estimates use both disjoint expression-derived proxies. Complete estimates and definitions are in Supplementary Table S4 and Data 3. |

| Supplementary Figure S10 | Figure_S10_Nested_LASSO_and_full_pipeline_permutation.png | All 20 repeat AUCs per source/contrast (38 SCI/10 controls) and all 999 null AUCs for each raw-count comparison. Black ticks show repeat means, not independent-repeat confidence intervals. Permutation comparator is the first nested repeat. Holm2 P = 0.002/0.002. Sex-conditional internal association only. |

| Supplementary Figure S11 | Figure_S11_Clinical_BTM_resampling.png | Four original M2 clinical estimates with HC3 pointwise 95% CIs and residual-permutation/wild-bootstrap P values. n = 43 (36 SCI/7 TC); each method has 9,999 draws. Both separate Holm4 families equal 1. Does not establish equivalence. |

## Supplementary data archives

### Supplementary Data 1

Supplementary_Data_1_Cohort_curation.zip

Cohort search, screening, curation rules and public GEO sample/reuse mapping. The clinical crosswalk is separately public in the repository; see Data 8.

Revision_search_log.csv

Dataset_screening.csv

Metadata_curation_rules.md

Public_GEO_sample_and_reuse_manifest.csv

### Supplementary Data 2

Supplementary_Data_2_Gene_GO_results.zip

Fixed43 gene and GO results, annotation sources and model summary. Historical_ members retain different original participant sets and model roles; their filenames distinguish them from the fixed43 results.

Genes_M0.csv

Genes_M1.csv

Genes_M2.csv

GO_CAMERA_M0.csv

GO_CAMERA_M1.csv

GO_CAMERA_M2.csv

Fixed_eight_GO.csv

GO_membership.csv

GO_tested_gene_universe.csv

Historical_full48_preranked_GO.csv

annotation_package_metadata.csv

annotation_sources.csv

bridge43_model_summary.csv

Historical_SCI_vs_TC_processed_limma.csv

Historical_SCI_vs_HC_processed_limma.csv

Historical_TC_vs_HC_processed_limma.csv

Historical_ODC405_DE_model_summary.csv

Historical_ODC405_DE_SCI_vs_HC_age_sex_prior_CNS.csv

Historical_ODC405_DE_SCI_vs_HC_unadjusted.csv

Historical_ODC405_DE_SCI_vs_TC_age_sex.csv

Historical_ODC405_DE_SCI_vs_TC_clinical_complete_case.csv

Historical_ODC405_DE_SCI_vs_TC_clinical_plus_CIBERSORTx.csv

Historical_ODC405_DE_SCI_vs_TC_clinical_plus_draw_time.csv

Historical_ODC405_DE_SCI_vs_TC_unadjusted.csv

Historical_ODC405_DE_TC_vs_HC_age_sex.csv

Historical_ODC405_DE_TC_vs_HC_unadjusted.csv

### Supplementary Data 3

Supplementary_Data_3_BTM_results.zip

Complete clinical and longitudinal BTM results, diagnostics, memberships, coverage and mappings. The primary16 family and other explicitly named families remain distinct.

BTM_primary16.csv

BTM_discovery_all88.csv

BTM_longitudinal_all240.csv

BTM_longitudinal_diagnostics.csv

BTM_members_annotated.csv

BTM_used_genes.csv

BTM_coverage.csv

BTM_discovery_symbol_mapping.csv

BTM_NCBI_symbol_mapping.csv

Marker_overlap_coverage.csv

### Supplementary Data 4

Supplementary_Data_4_Original_15plus15_scores.zip

Original 15+15 score analyses, selection sensitivities, later annotation audits, reference parameters, component correlations, alternative balance scales and descriptive three-group comparisons. All remain exploratory and conditional on discovery selection. The signature_resolution_20261027 folder adds a reconciliation of all 106 candidates, 103 archived frequency records and the frozen 30-gene output, with source hashes and a BTM marker-dependency audit. It changes no score or reported model estimate.

#### Additional provenance files in signature_resolution_20261027

signature_resolution_20261027/4_archived_source.Rmd

signature_resolution_20261027/4_comment_clarified_review_copy.Rmd

signature_resolution_20261027/Archived_rule_line_map.json

signature_resolution_20261027/Archived_threshold_reconciliation.csv

signature_resolution_20261027/BTM_members_annotated.csv

signature_resolution_20261027/BTM_proxy_exclusion_dependency_audit.csv

signature_resolution_20261027/Candidate_to_frozen_signature_audit.csv

signature_resolution_20261027/Frozen_30_gene_selection_audit.csv

signature_resolution_20261027/GSE151371_stability_selection_all_candidates.csv

signature_resolution_20261027/MCPcounter_genes.txt

signature_resolution_20261027/README_signature_resolution.md

signature_resolution_20261027/Signature_resolution_audit.json

signature_resolution_20261027/Source_evidence_manifest.csv

signature_resolution_20261027/candidate_gene_pools.csv

signature_resolution_20261027/locked_signature_v1.json

Candidate_gene_pools.csv

Frozen_15plus15_gene_list.csv

Selection_parameter_sensitivity.csv

Fixed_sample_models27.csv

Longitudinal_contrasts36.csv

Longitudinal_sensitivity42.csv

Clinical_composition_models15.csv

Paired_sensitivities9.csv

Longitudinal_model_diagnostics.csv

component_correlations.csv

balance_scale_correlations.csv

balance_sensitivity_scores.csv

discovery_three_group_descriptive.csv

locked_gene_biological_annotation.csv

candidate_annotation_followup.csv

processed_matrix_topN_B_sensitivity.csv

R1_11_topN_B_frequency_summary.csv

standardization_parameters.csv

NCBI_standardization_parameters.csv

### Supplementary Data 5

Supplementary_Data_5_Clinical_technical_sensitivities.zip

Clinical and technical sensitivity analyses, covariate availability/support and BTM TC-deletion diagnostics.

Covariate_availability_summary.csv

Clinical_and_platform_sensitivity11.csv

Technical_PC_associations25.csv

Aggregate_PC_covariate_associations.csv

Clinical_category_counts.csv

Covariate_common_support.csv

Technical_availability.csv

BTM_TC_deletion_diagnostics.csv

Original_score_TC_leave_one_out_models.csv

Original_score_TC_leave_one_out_summary.csv

Original_score_model_assumption_diagnostics.csv

Original_score_all_participant_influence.csv

Original_score_design_column_VIF.csv

Original_score_covariate_center_constants.csv

### Supplementary Data 6

Supplementary_Data_6_Reuse_and_additional_cohort.zip

Additional-cohort reuse and same-platform sensitivity outputs; all nine controls remain reused.

GSE293559_reuse_sensitivity6.csv

Reuse_sensitivity6.csv

Same_platform_sensitivity6.csv

Technical_PC_associations.csv

Technical_PCA_variance.csv

### Supplementary Data 7

Supplementary_Data_7_Multiplicity_and_disposition.zip

Multiplicity registry, complete family membership, withdrawal boundaries and detailed conditional gene-clustering diagnostics.

Statistical_family_registry.csv

Statistical_family_membership.csv

Claim_evidence_boundaries.json

Withdrawn_analysis_disposition.csv

Cluster_metrics.csv

Cluster_decision.csv

Cluster_bootstrap_summary.csv

cluster_gap.csv

cluster_initializations.csv

cluster_bootstrap_ARI.csv

cluster_bootstrap.csv

### Supplementary Data 8

Supplementary_Data_8_Reproducibility_inventory.zip

Public repository and verified commit record, source-to-portable mapping, software environments, seeds and verification scope. The included repository manifest describes the verified analysis commit, not this supplementary ZIP.

Reproducibility_inventory.json

renv.lock

requirements.txt

Random_seeds.json

Original_analysis_sessionInfo.txt

Verified_pipeline_sessionInfo.txt

Source_to_portable_code_mapping.csv

Verified_commit_manifest.csv

Verified_commit_SHA256SUMS

Historical_protocol_fixed43.md

Historical_protocol_BTM_extension.md

Historical_protocol_fixed_scores.md

Historical_initial_BTM_protocol_lock.json

### Supplementary Data 9

Supplementary_Data_9_Nested_validation_and_resampling.zip

New nested-classification and fixed-BTM resampling audit. Includes the locked protocol, source code, input snapshots, all observed and null outputs, diagnostic checks and provenance. Current interpretation remains conditional on the stated recipe and exchangeability/error-model assumptions.

MANIFEST.csv

PROTOCOL.md

README.md

R_sessionInfo.txt

analysis-report.md

audit/R_dependency_versions.csv

audit/analysis_verification.json

audit/clinical_resampling_verification.json

audit/implementation_tests.json

audit/method_reference_verification.json

audit/numerical_implementation_correction.md

audit/protocol_input_lock.json

audit/python_environment.json

figure-catalog.md

inputs/btm_frozen_models.csv

inputs/btm_model_inputs.csv

inputs/historical_source.Rmd

inputs/participants58.csv

inputs/processed_expression.csv

inputs/raw_counts.csv

requirements-analysis.txt

scripts/make_figures.R

scripts/nested_core.R

scripts/run_clinical_resampling.R

scripts/run_nested.R

scripts/summarize_results.R

scripts/verify_results.py

stats-appendix.md

tables/clinical_design_matrix.csv

tables/clinical_resampling_assignments.rds

tables/clinical_resampling_diagnostics.csv

tables/clinical_resampling_inputs.csv

tables/clinical_resampling_null.csv

tables/clinical_resampling_results.csv

tables/full_pipeline_permutation_null.csv

tables/full_pipeline_permutation_results.csv

tables/nested_conditional_intervals.csv

tables/nested_gene_selection_frequencies.csv

tables/nested_inner_screen_audit.csv

tables/nested_mean_oof_predictions.csv

tables/nested_oof_predictions.csv

tables/nested_outer_fit_audit.csv

tables/nested_performance_summary.csv

tables/nested_repeat_metrics.csv

tables/nested_selection_overlap.csv

tables/nested_split_assignments.csv

tables/new_resampling_family_registry.csv

tables/permutation_label_assignments.csv

tables/permutation_outer_fit_audit.csv

## Interpretation and access boundaries

The original numerical outputs, including the 19 October fixed SCI–TC freeze, remain unchanged. Data 9 contains newly executed revision-stage nested and resampling results, with a protocol locked before these new fits.

The eight displayed immune GO terms are a subset of the full 4,032-set family. Their nonsignificance does not mean that all GO tests were nonsignificant. Pointwise confidence intervals do not establish equivalence.

Published BTMs check external gene definitions within existing cohorts. Discovery-derived 15+15 scores and the partially reused additional cohort do not provide independent participant validation.

The old nine-gene model and supervised postselection P values remain withdrawn. Data 9 provides a newly specified nested candidate-union LASSO, full-pipeline sex-stratified permutations and fixed-BTM residual/wild resampling; these do not restore the old model or provide independent validation.

The repository is publicly accessible. No separately archived DOI-bearing release or Git tag has been created, and original-code licensing remains an author decision. Source-specific licenses remain in force. The current repository tree includes Data 9 under extensions/nested_validation and the Data 4 provenance addition under reference_results/supplementary_data/Supplementary_Data_4_Original_15plus15_scores/signature_resolution_20261027. They were not part of the earlier analysis snapshot.

Data 8 documents standard-pipeline numerical checks and the optional historical reruns, clean-machine installation and cross-platform checks that were not all performed.

## Expanded repository paths

Data 1-8 remain under `reference_results/supplementary_data/`; Data 9 is expanded
under `extensions/nested_validation/`. Supplementary figures and tables are under
`reference_results/manuscript/supplement/`. The seven additional redraws are in
`extensions/figure_redraw/` and have not replaced formal supplementary numbering.
