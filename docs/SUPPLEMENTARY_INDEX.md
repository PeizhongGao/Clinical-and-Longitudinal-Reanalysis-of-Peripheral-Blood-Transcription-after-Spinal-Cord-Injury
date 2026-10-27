# Supplementary Materials Index

Clinical and Longitudinal Reanalysis of Peripheral Blood Transcription after Spinal Cord Injury

This index identifies seven supplementary tables, nine figures and eight data archives. Archive members are listed below. Numerical source files and figure assets are retained unchanged; the archive coverage and reproducibility documentation were updated on 27 October 2026.

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

Original 15+15 score analyses, selection sensitivities, later annotation audits, reference parameters, component correlations, alternative balance scales and descriptive three-group comparisons. All remain exploratory and conditional on discovery selection.

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

## Interpretation and access boundaries

All numerical outputs are existing results. The main-text SCI–TC numbers come from the unchanged 19 October 2026 freeze. No models were rerun for this packaging update.

The eight displayed immune GO terms are a subset of the full 4,032-set family. Their nonsignificance does not mean that all GO tests were nonsignificant. Pointwise confidence intervals do not establish equivalence.

Published BTMs check external gene definitions within existing cohorts. Discovery-derived 15+15 scores and the partially reused additional cohort do not provide independent participant validation.

Withdrawn LASSO and supervised postselection-permutation outputs are not reinstated as current results. No nested LASSO or repaired label-permutation test is claimed.

The repository is publicly accessible. No separately archived DOI-bearing release or Git tag has been created, and original-code licensing remains an author decision. Source-specific licenses remain in force.

Data 8 documents standard-pipeline numerical checks and the optional historical reruns, clean-machine installation and cross-platform checks that were not all performed.

