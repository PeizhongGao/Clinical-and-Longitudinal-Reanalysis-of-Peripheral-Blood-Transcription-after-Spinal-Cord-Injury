# GitHub import preparation — 27 October 2026

This repository folder was populated from the validated 1.0.0 upload bundle.
Author-facing Chinese guides are archived outside this repository; current user
documentation is in English. Source-language protocols are retained as described in
ARCHIVED_RECORDS.md. The README and citation title match the current revised manuscript.
This is a packaging/documentation update; all numerical analysis sources, inputs
and frozen reference results remain identical to the validated upload bundle.

The repository-specific .gitattributes sets `* -text` to prevent automatic end-of-line
conversion. The source inputs and their hashes are part of the scientific record;
changing LF/CRLF during a Git checkout would otherwise invalidate the byte checks.
Normal text diffs remain available. Do not replace this with automatic normalization
without explicitly rebuilding the manifests and rechecking the input locks.

Run `python check_upload.py` after a clone or download to check the files against
manifest.csv and SHA256SUMS. The existing validation records retain the actual
versions and dates on which numerical tests were performed.

The first public push and unauthenticated access verification were completed on
27 October 2026 at analysis commit `6415b0bdc3dff91d3949e940c363ed22f48f387c`.
Current access information is in AVAILABILITY_STATUS.md. The original-code license
remains undecided, and no archive DOI, Git tag or GitHub Release is assigned.
Historic NOT_UPLOADED records retain their recorded preparation state.
