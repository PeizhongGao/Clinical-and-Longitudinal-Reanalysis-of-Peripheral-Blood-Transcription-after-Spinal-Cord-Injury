# Public package and remaining release choices

The code and supporting-data repository is public. Anonymous access and downloaded
file hashes were verified on 27 October 2026; see `AVAILABILITY_STATUS.md`.

## Included and checked

- Actual scientific source and public source matrices, metadata and checkpoints.
- Curation, overlap and exclusion rules, fixed gene lists and random seeds.
- Complete retained results, including negative and uncertain estimates.
- Executable workflow, environment records, dependency lock and source/file hashes.
- Recorded standard-pipeline numerical validation and clearly stated untested scope.
- Current supplementary archive contents indexed in `SUPPLEMENTARY_INDEX.md`.

## Author choices still open

- Select an original-code reuse license; this was explicitly deferred.
- Optionally create a Git tag/Release and separate archival deposit for a DOI.

Version 1.0.0 is package metadata, not a claim that a Git tag or DOI-bearing release
exists. Do not apply a blanket software license to third-party data or components.
Use the verified repository and analysis commit in correspondence.

## Updating the repository

Use this directory as the repository root and respect `.gitignore` and the
byte-preserving `.gitattributes`. Rebuild manifests after an intentional package
change, then run `python check_upload.py`. Keep manuscripts, reviewer correspondence,
local caches and working outputs outside this public repository.
