# Code completeness for the upload bundle

This audit concerns the retained manuscript and supplementary analyses, not every
historical exploration or manuscript-formatting utility in the original workspace.

## Current finding

No execution source is missing from the declared reproduction workflow. The main
Figure 1, Figures S1-S9 and retained numerical tables have analysis/display code or
an explicitly documented fixed reference source. See `metadata/CODE_COMPLETENESS.csv`
for every pipeline entry and `metadata/FIGURE_TABLE_CODE_MAP.csv` for the manuscript mapping.
`check_upload.py` independently reads the entry points from the runner and checks
their presence, package hashes and Python syntax.

## Included code and actual test coverage

| Component | Code included | Validation status |
|---|---|---|
| Sample metadata and overlap curation | Yes | Standard workflow passed |
| Original 15+15 projection and revision-stage sensitivity | Yes | Standard workflow passed |
| Full-count longitudinal TMM and scores | Yes | Standard workflow passed |
| BTM mapping, scoring and all 88 clinical / 240 longitudinal results | Yes | Standard workflow passed |
| ODC clinical mapping and 27 fixed-sample models | Yes | Standard workflow passed |
| Original-score, platform and participant-reuse sensitivities | Yes | Standard workflow passed |
| 43-participant gene/GO analyses and locked-input reconstruction checks | Yes | Standard workflow passed |
| Main and supporting figure rendering and result export | Yes | Standard workflow passed |
| Original count-based 500-iteration signature selection | Yes | Optional; startup/early iterations exercised, not fully rerun |
| Historical clustering and whole-cohort preranked GO permutations | Yes | Optional extended profile; not fully rerun in the packaging test |
| Dependency restoration | Yes, with R lockfile and Python pins | Current environment checked; clean-machine restoration not certified |

The completed standard test used a package relocated outside the original project:
19 analysis stages and 21 validation groups passed; 357,077 numerical values were
compared with frozen references. This is reproducibility evidence, not external
biological validation. The original selected gene identities were checked as well.

## Intentional scope exclusions, not lost source files

- FASTQ alignment, source-study ComBat preprocessing and CIBERSORTx estimation:
  this secondary-analysis package begins at supplied public count/processed matrices
  and cell fractions. It does not claim to provide or rerun those upstream pipelines.
- Withdrawn nine-gene LASSO/prediction and postselection label-permutation branches:
  excluded from current execution; no nested LASSO was performed or fabricated.
- Old duplicate exploration notebooks, Word/report generation and reviewer-letter
  formatting utilities are not required to reproduce the retained scientific analyses.
- Fixed definitions, selection checkpoints and historical diagnostic results are
  explicitly identified as reused inputs; their presence is not misreported as a new fit.

## Remaining author/deposit fields

The original-code license remains undecided at the author's request. Repository
ownership, actual publication URL and archive DOI can only be filled at deposition.
They are publication metadata decisions, not missing code. This package has not
been uploaded, and no permanent identifier has been invented.
