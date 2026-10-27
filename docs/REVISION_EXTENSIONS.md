# Revision extensions added 28 October 2026

Package metadata version 1.1.0 extends the earlier analysis snapshot
`6415b0bdc3dff91d3949e940c363ed22f48f387c` and documentation commit `bb91f70`.
This version number does not imply a Git tag, GitHub Release or archive DOI.
The original-code license remains pending author decision; third-party rights remain applicable.

## Included materials

| Component | Public location | Evidence role |
|---|---|---|
| Data 4 provenance addition | `reference_results/supplementary_data/Supplementary_Data_4_Original_15plus15_scores/signature_resolution_20261027/` | Reconciles 106 candidates, 103 archived selection records and the frozen 30-gene output; does not establish historical prespecification |
| Data 9 complete analysis extension | `extensions/nested_validation/` | New nested LASSO, complete-pipeline permutations and fixed-BTM resampling, including original protocol, inputs, code, all outputs and environment records |
| Tables S8-S9 and Figures S10-S11 | `reference_results/manuscript/supplement/` | Saved summary tables and figures for the new analyses |
| Seven redrawn figures | `extensions/figure_redraw/` | R drawing code, source CSVs, vector PDF/SVG and 600-dpi TIFF; descriptive or conditional evidence as stated in each legend |

The redrawn Figure 1-7 filenames identify the redraw set. They are not a new formal
manuscript numbering scheme. The supervisor reading copy numbers its BTM main figure
and these seven displays consecutively as Figures 1-8. Formal manuscript and
supplementary numbering must be resolved in the submission package separately.

## Reproduce safely

The existing `reproduce.py` standard workflow remains unchanged. Use the new runner
from the repository root; each command requires a new output directory:

```sh
python reproduce_revision.py --task verify --output outputs/revision-check
python reproduce_revision.py --task pilot --output outputs/revision-pilot
python reproduce_revision.py --task figures --output outputs/redrawn-figures
python reproduce_revision.py --task full --output outputs/revision-full
```

Supply `--rscript /path/to/Rscript` when needed. The runner copies work to the output
directory and does not write into the published reference inputs or results.
`verify` independently recalculates saved statistics and checks split/permutation
records; it does not refit classifiers. `pilot` additionally executes the original
implementation checks and one nested repeat. `full` refits all observed repetitions,
1,998 full-pipeline null runs and clinical resampling distributions; it is computationally
expensive and uses the archived six-worker permutation setting. `figures` renders
all seven displays from the included source CSVs without refitting biological models.

R and Python dependencies are recorded in the root environment files and
`extensions/nested_validation/audit/R_dependency_versions.csv`. Original analysis
code and frozen inputs are retained byte-for-byte. The figure drawing copy changes
only the output/source-directory resolution; its original source is preserved in
`extensions/figure_redraw/provenance/`.

## Interpretation and historical records

The old 5+4 model and its postselection inference remain withdrawn. The new nested
classifier is a separate internal audit, not validation of the frozen 15+15 score.
Published BTM definitions are not an independent participant cohort. New results
do not alter the original BH16 or the 93 named historical testing families.

Archived Data 7/Data 8 records and original notebooks describe their original
versions. Statements there that analyses had not yet been executed or uploaded are
historical; the current scope is described here and in the root README. The Data 9
README and inventory have been adapted for repository use, with original versions
preserved under `provenance/revision_extension_20261028/`.

Participant-linked clinical inputs are derived from the already-public ODC-SCI 405
dataset (DOI 10.34945/F5QC7J, CC BY 4.0) and the existing public crosswalk; dataset
attribution is retained. No manuscript drafts or reviewer correspondence are included
in this update.
