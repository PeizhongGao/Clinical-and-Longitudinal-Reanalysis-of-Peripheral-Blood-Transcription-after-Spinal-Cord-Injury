# Licensing status of the public package

The authors have explicitly deferred choosing a license for their original analysis
code. This public package therefore does **not** grant an MIT, GPL or other new open-source
license for author-owned code. Resolve this item before describing the release as
open source or selecting a repository-wide license field. Code visibility and an
open-source reuse license are different matters.

Do not apply one blanket license to every file in this mixed code/data deposit.

- ODC-SCI accession 405: the public source dataset is CC BY 4.0. Preserve its
  attribution, dataset DOI, license reference and the documented extraction changes.
- BTM: the pinned upstream repository license is retained in
  `licenses/third_party/BTM_LICENSE.txt`. Its original copyright notice is unchanged.
- MCPcounter: the included source is identified by its pinned commit. The upstream
  package DESCRIPTION specifies GPL-3; the DESCRIPTION and GPL-3 text are included.
  The project does not relicense that third-party component.
- GEO source matrices, annotations and original-study tables: source accessions and
  provenance are supplied. No new ownership or blanket CC license is asserted over
  these third-party materials.
- R/Python dependencies are not redistributed as installed binaries. Their respective
  upstream licenses apply when installed.

Public GitHub access is verified; see `docs/AVAILABILITY_STATUS.md`. Version 1.0.0
is package metadata. No separate archive DOI, Git tag or GitHub Release has been
assigned, and no original-code open-source license is claimed.
