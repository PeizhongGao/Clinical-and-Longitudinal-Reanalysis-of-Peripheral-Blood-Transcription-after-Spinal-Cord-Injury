# Third-party attribution and rights

## BTM definitions

- Source repository: https://github.com/shuzhao-li/BTM
- Pinned revision: `e973bf374956d40367487492d825de97f96ae095`
- File: `BTM/datasets/BTM_for_GSEA_20131008.gmt`
- Associated article: Li et al., *Nature Immunology* 15, 195–204 (2014), DOI
  `10.1038/ni.2789`.
- The upstream MIT license is retained verbatim in `third_party/BTM_LICENSE.txt`.
  Its copyright notice is the upstream notice, not a declaration of project authorship.

## MCPcounter

- Source: https://github.com/ebecht/MCPcounter
- Pinned R-source revision: `b6eac73e91c246fcff0bb1a5c68a816cd588fc48`
- File: `Source/R/MCPcounter.R`
- Upstream authors in DESCRIPTION: Etienne Becht, Aurelien de Reynies, Florent Petitprez.
- Version 1.2.0 DESCRIPTION specifies `License: GPL-3`.
- See `third_party/MCP_DESCRIPTION.txt` and `third_party/GPL-3.0.txt`.
- The marker file and implementation are third-party components. No claim is made
  that an undecided license for the authors' scripts replaces the upstream GPL.

## ODC-SCI 405

See `docs/DATA_PROVENANCE.md` and the complete DataCite record. Dataset DOI
`10.34945/F5QC7J`; version 1.0; CC BY 4.0. Source files and public participant-linked
covariate extraction are supplied with attribution and transformation descriptions.

## Other public sources and software dependencies

GEO matrices and annotations retain their source accessions and third-party attribution.
The original supplemental workbooks retain their original contents. This deposit does
not assign a blanket new license to them. GO/organism annotation versions are fixed in
the supplied input metadata. Installed R/Python libraries are not bundled; they remain
subject to their upstream licenses, as recorded in package metadata and the environment
lock. The optional gseapy dependency is BSD-3-Clause and is obtained through its normal
package distribution rather than a copied machine-specific binary.
