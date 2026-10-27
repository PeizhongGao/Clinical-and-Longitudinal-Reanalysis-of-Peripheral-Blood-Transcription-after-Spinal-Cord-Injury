# ODC-SCI 405 clinical and cell-composition sensitivity protocol

Locked on 12 October 2026 before fitting the new models. This is a revision-stage analysis after the original study, candidate selection and prior results were known; it is not prospective prespecification.

## Purpose and fixed source

Use the author-supplied ODC-SCI 405 Full Data Package and its dictionary to add individual-level covariates and CIBERSORTx-derived cell-abundance estimates for the 58 GSE151371 participants. Keep the GEO-supplied ComBat-normalized expression matrix and the original 17,498-gene tested universe used in the R1 reproduction. Do not use the ODC long table as a replacement expression matrix or change the frozen 15+15 signature. Link by the exact group-plus-number sample key (HC01–HC10, TC01–TC10, SCI01–SCI38), verify against the GEO title and existing GSM manifest, and report the crosswalk.

## Covariates and missing data

Report sample-level availability by group for age, sex, race, prior CNS pathology, concurrent TBI, Injury Severity Score (ISS), blood-draw time, AIS and neurological level. Preserve source-coded `Unknown` values as an explicit category. Do not impute missing numeric values. AIS/NLI apply only to SCI participants and are not treated as common covariates in SCI-versus-control models. Infection, treatment, batch, RIN and processing date are unavailable in this package; the GEO record describes ComBat processing, but the batch labels are not supplied for independent audit. Race will be summarized and assessed for overlap; sparse/non-overlapping categories will not enter the primary models.

## Fixed differential-expression models

Use limma empirical-Bayes models on the original 17,498-gene matrix. First reproduce all three unadjusted pairwise fits (`expression ~ group`) as a source check. For the primary SCI-versus-TC comparison, fit, in order: (1) `expression ~ group + age + sex`; (2) complete cases with `expression ~ group + age + sex + ISS + concurrent_TBI + prior_CNS_pathology`; (3) the same model plus blood-draw time, on its complete-case subset; and (4) the same complete-case clinical model plus standardized CIBERSORTx neutrophil and aggregate T-cell estimates. The CIBERSORTx terms are expression-derived proxy sensitivities, not measured counts or independent covariates. Do not interpret attenuation as a cell-intrinsic effect.

For SCI-versus-HC, fit `expression ~ group + age + sex + prior_CNS_pathology` using available observations; for TC-versus-HC, fit `expression ~ group + age + sex`, because injury severity, TBI and post-injury draw time are not defined comparably for uninjured controls. Center continuous covariates within each fitted subset. Check design rank, residual degrees of freedom and group sample counts. Apply BH correction across the 17,498 genes separately within each model/contrast. Keep existing thresholds and pathway analyses unchanged; these fits are covariate sensitivities, not replacement discovery analyses.

## Frozen score and cell-abundance sensitivities

Reuse the existing TC-referenced, frozen innate/adaptive/balance scores. Fit participant-level score models for SCI versus TC under the same unadjusted, clinical, draw-time and clinical-plus-CIBERSORTx specifications. Use HC3 robust standard errors for these score regressions and BH correction across the three scores within each specification. Label these selected-discovery scores as post-selection sensitivities, not independent tests. Summarize the 22 CIBERSORTx relative-abundance estimates by group, emphasizing neutrophils and the sum of T-cell subsets; verify the reported components sum to approximately one and retain all participants. Correct the 22 cell-type SCI-versus-TC tests as one BH family and the two focused neutrophil/aggregate-T comparisons as a separate BH family. Although the dictionary labels the components as percentages, preserve the observed supplied scale when their sum is approximately one. These are algorithmic estimates from bulk expression, not measured blood counts or validated cell fractions.

Run PCA on the original 17,498-gene discovery matrix and report exploratory associations of the first five PCs with group and available mapped clinical/cell variables. This does not test unavailable batch, RIN or processing variables. Provide a PCA plot and a score-model effect plot; use these as supplementary diagnostics rather than changing the main biological claims.

## Interpretation boundary

No clinical values are fabricated, no new genes are selected, and no original analysis is overwritten. Clinical adjustment and estimated-cell sensitivities can strengthen but cannot eliminate residual clinical/technical confounding. If model rank, cell counts or mapping checks fail, retain the blocker and do not promote the affected analysis.
