# Clinical and Longitudinal Reanalysis of Peripheral Blood Transcription after Spinal Cord Injury

**Manuscript:** Neurology International, neurolint-4585750.  
**Package version:** 1.0.0, prepared 26 October 2026.  
**Materials:** analysis code, supporting data and reproducibility records for the manuscript.
Original-code licensing is awaiting the authors' decision; no permanent archive DOI is
recorded. See `LICENSE_STATUS.md` for the current licensing status.

This package supplies actual analysis source, public-study data inputs, metadata
curation, frozen gene definitions, complete retained statistical families, reference
results, and an executable reproduction workflow. It is not merely a script inventory.
Read `docs/RELEASE_CHECKLIST.md` for deposition guidance and
`docs/CODE_COMPLETENESS.md` for the code audit.
Run `python check_upload.py` to verify the package after a clone or download.

## What the workflow reproduces

- Public-source sample curation and overlap mapping for the retained cohorts.
- The original frozen 15+15 score projection and revision-stage selection sensitivity.
- GSE226238 TMM normalization from the complete available count matrix.
- Four published BTM definitions and the candidate-exclusion sensitivity, with symbol
  mapping, coverage, score calculation and composition-proxy construction.
- ODC405 clinical extraction/crosswalk, fixed-sample models and original-score
  clinical, visit, platform and reuse sensitivities.
- All 88 BTM clinical models and 240 longitudinal contrasts, including the designated
  16-test family. Non-significant results are retained.
- The fixed 43-participant, 17,498-gene and 4,032-GO-set analyses under M0/M1/M2.
- Main Figure 1 and supporting figures, with numeric validation against frozen results.

The default `standard` profile reuses the **original signature-selection checkpoint**
and **historical clustering/preranked-GO checkpoints**. It recomputes current matrix-based
analysis and the revision-stage top-N/iteration sensitivity. The original count-based
selection can also be run with `--original-selection`. The costly historical clustering
and preranked GO procedures are included in `extended`; they are explicitly separated
from the current inferential core. See `docs/VALIDATION.md` for what was actually tested.

## Quick start

Use Python 3.12 and the recorded R 4.6.1 environment. Install Python dependencies:

```shell
python -m pip install -r environment/requirements.txt
```

Restore R package versions from the provided 102-package lockfile:

```shell
Rscript --vanilla environment/restore.R .
```

The restored R library is `environment/r_library/`; the runner automatically passes it
to the child R processes. If using an existing environment, exact versions of limma
3.68.4 and org.Hs.eg.db 3.23.1 are checked before numerical analysis.

From this folder, run:

```shell
python reproduce.py --output outputs/run01
```

If `Rscript` is not on PATH, supply its installed executable using `--rscript`.
To rerun the original count-based selection as well:

```shell
python reproduce.py --output outputs/run02 --original-selection
```

For historical clustering and preranked-GO permutations:

```shell
python -m pip install -r environment/requirements-extended.txt
python reproduce.py --output outputs/full --profile extended --original-selection
```

Always use a new empty output directory. `--resume` is restricted to the same profile,
selection setting, package inputs and code; changed inputs require a fresh run.
Outputs are generated in `OUTPUT/work`, logs in `OUTPUT/logs`, validation in
`OUTPUT/validation.json`, and conveniently named manuscript figures/tables in
`OUTPUT/results`. The runner never writes to `inputs/` or `reference_results/`.

## Folder structure

| Folder/file | Purpose |
|---|---|
| `reproduce.py` | Single command-line execution entry point |
| `code/` | Orchestrator support and portable copies of the actual analysis scripts |
| `inputs/project/` | Source matrices, curated metadata, fixed definitions and named checkpoints |
| `reference_results/` | Unchanged numerical reference outputs and all final supporting results |
| `environment/`, `renv.lock` | Installation instructions and recorded dependencies |
| `metadata/` | File inventory, provenance mapping, seeds and figure/table crosswalk |
| `validation/` | Test results and stated remaining reproducibility limits |
| `licenses/` | Third-party notices and license texts |
| `historical/code/4.Rmd` | Original signature notebook; supplementary plotting blocks retained for provenance |
| `provenance/` | Prior source/consistency audit evidence |
| `docs/ARCHIVED_RECORDS.md` | English guide to preserved source-language protocol and evidence records |
| `upload_archives/` | Optional archive for a repository deposit; excluded from Git to avoid duplication |

The older date-stamped directory names inside `inputs/project` and `code/project` are
retained to preserve the established script-to-file relationships. They are internal
relative names, not requirements to have the author's original project directory.

## Scientific scope

This is an exploratory reanalysis of public count/processed-expression matrices,
not a FASTQ alignment pipeline. Discovery processing and source-study cell-fraction
estimates are reused at the documented input level. Published BTM definitions provide
an external gene-definition check, not independent participants. The 15+15 scores are
postselection exploratory. GSE293559 reuses 10 SCI participants and all 9 controls from
GSE226238; exclusions do not create fully independent controls. The longitudinal
comparison lacks matched trauma controls over time and does not establish SCI-specific
recovery, causality or a cell-intrinsic mechanism.

The withdrawn 5+4 LASSO/prediction and postselection label-permutation claims are not
restored. No nested LASSO was performed. Complete retained families and negative
findings remain available, including the full BTM16 and fixed immune-GO panels.
The unchanged original supplemental workbooks and their `original_*` imports can
contain superseded nine-gene/permutation sheets. Those records are source provenance,
not current statistical evidence; the runner does not refit the withdrawn models.
Historical protocol-generation and superseded QA utilities are explicitly filed under
`historical/auxiliary_scripts`, rather than presented as current entry points.

## Reuse, citation and deposition

See `docs/DATA_PROVENANCE.md`, `licenses/THIRD_PARTY_NOTICES.md`, and `CITATION.cff`.
ODC405 is a reused public CC BY 4.0 dataset; its participant-linked public clinical
fields and the derived crosswalk are included with attribution, not described as
prohibited-to-redistribute data. Dataset-specific rights continue to apply.

The code and supporting data are publicly available at https://github.com/PeizhongGao/Clinical-and-Longitudinal-Reanalysis-of-Peripheral-Blood-Transcription-after-Spinal-Cord-Injury.
Anonymous access and the downloaded archive were verified on 27 October 2026.
The analysis snapshot cited in the manuscript is commit `6415b0bdc3dff91d3949e940c363ed22f48f387c`.
Version 1.0.0 is the package metadata version; no Git tag, GitHub Release or
separate archive DOI has been assigned. The original-code license remains
awaiting author decision. See `docs/AVAILABILITY_STATUS.md` for the verified scope.

The 27 October documentation update also supplies the previously omitted
supplementary annotation, score sensitivity and clustering tables without
changing scientific calculations. `docs/SUPPLEMENTARY_INDEX.md` lists the
current archive contents; `reference_results/supplementary_data` contains
their expanded files. Historical validation records retain their original dates.

GitHub import preparation and byte-preserving file handling are documented in
`docs/GITHUB_IMPORT.md`. The article title above follows the current revised manuscript.
