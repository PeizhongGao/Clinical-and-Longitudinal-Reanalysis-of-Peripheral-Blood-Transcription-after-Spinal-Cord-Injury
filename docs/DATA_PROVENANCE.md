# Data provenance and access routes

This is a secondary analysis of reused public studies. The package contains no newly
generated sequencing reads. Source matrices and public clinical fields are retained
at their published analysis-input level; missing information is not invented.

| Source | Packaged input | Access/provenance |
|---|---|---|
| GSE151371 | Gene counts, published ComBat-processed expression, public GEO annotations | GEO accession GSE151371; processed-expression and count pipelines have distinct roles |
| GSE226238 | NCBI count matrix and matching GRCh38 annotation; author-processed workbook/matrix; public sample and clinical mapping | GEO accession GSE226238; 10 SCI participants/37 visits plus 9 control libraries |
| GSE293559 | Published normalized expression and curated library metadata | GEO accession GSE293559; 36 SCI/9 controls; 10 SCI and all 9 controls are reused from GSE226238 |
| GSE177034 | Public sample annotation only | Historical non-SCI metadata audit; excluded from retained signature inference; its count matrix is deliberately not needed in this package |
| ODC-SCI 405 | Original public long table, data dictionary, source ZIP and derived clinical mapping | DOI 10.34945/F5QC7J, version 1.0; CC BY 4.0 |
| BTM | `BTM_for_GSEA_20131008.gmt`, mapping/coverage/used-member tables | shuzhao-li/BTM at commit e973bf374956d40367487492d825de97f96ae095; original 2013 definitions, not BTM Plus |
| MCPcounter | Published marker table and original R implementation | ebecht/MCPcounter source at commit b6eac73e91c246fcff0bb1a5c68a816cd588fc48; upstream GPL-3 |
| GO | Frozen GO BP terms, gene memberships, annotation metadata and 43-participant input indices | Fixed existing 4,032-set membership; source annotation versions recorded in renv.lock and input metadata |

## ODC405 attribution and transformations

Kyritsis, N., Torres-Espin, A., Schupp, P. G., and colleagues (2021).
*Acute blood RNA sequencing from SCI human patients*, version 1.0.
Open Data Commons for Spinal Cord Injury. DOI: **10.34945/F5QC7J**.
The full creator list and license metadata are provided in
`metadata/ODC405_datacite_source.json`; this shortened display citation does not
replace that complete record.

Source record: https://odc-sci.org/data/405  
License: https://creativecommons.org/licenses/by/4.0/  
Data-use rules: https://odc-sci.org/about/odc-data-use  
Terms: https://odc-sci.org/about/odc-terms

The source CSV is a long expression/clinical table (1,015,058 rows, 61 fields),
not 1,015,058 participants. It covers 58 public study-specific participant IDs:
38 SCI, 10 trauma controls and 10 healthy controls. The source dictionary is included.
The analysis extracts participant-level public covariates, checks within-ID consistency,
matches the study-specific group/number to the expression sample and GEO GSM, and
preserves the crosswalk and missingness. Clinical complete-case sets contain 43 or 41
participants under the recorded definitions. Unknown categories are not silently
turned into clinically observed absence.

Public source-study CIBERSORTx fractions are reused as published inputs. This package
does not rerun CIBERSORTx or redistribute a proprietary CIBERSORTx implementation.
The analyses derive summaries and composition-adjusted models from those supplied
fields and from explicitly identified expression-based marker proxies.

The public ODC405 record is CC BY 4.0. Its registration/login requirement is not used
as a basis to claim that redistribution of this published dataset is prohibited.
The original dataset license and attribution continue to apply to reused data;
we do not apply the authors' prospective software license to those data.

## Sample identity and exclusions

`dataset + sample` identifies a library within an accession. `subject` identifies the
participant for longitudinal models; `gsm` identifies the public specimen record.
An accession reusing a specimen is not a new independent participant. The curated
manifest, its executable curation code and detailed rules are provided.

The GSE293559 full comparison is 36 SCI/9 controls. Removing reused SCI leaves
26 SCI/9 controls; all nine controls are still reused. Same-platform restrictions
and missing 12-month visits are documented in the manifests. No new exclusion rule
was chosen during this packaging operation.

## Input tiers and numerical authority

`inputs/project` deliberately retains source files and named checkpoints so readers
can inspect both upstream and downstream steps. The standard runner regenerates
major intermediate states and statistical outputs in a new workspace. It explicitly
reuses fixed gene/GO definitions, original supplemental-score inputs and historical
checkpoints where stated; no checkpoint should be mistaken for a newly observed dataset.

`reference_results` is the unchanged numerical comparison target. Do not overwrite it
with a new run. The October 19 numerical freeze remains the manuscript authority.
`metadata/source_to_release_manifest.csv` records original and public-copy hashes.
`manifest.csv` and `SHA256SUMS` describe the final packaged files.

## Packaging amendments

Only the public copies were relocated. Absolute author-machine paths in provenance
records are replaced by `PROJECT_ROOT` or relative identifiers. The executable entry
point resolves input/output locations from the package and supplied output argument.
The public bridge manifest is explicitly a **post-results portability lock**, not a
new prospective protocol; the original protocol hash is retained. Historical notes
that classified ODC405 as restricted are superseded by this source-specific assessment.

No changes were made to statistical methods, participant inclusion, gene lists,
GO membership, multiplicity families or the original frozen results.
