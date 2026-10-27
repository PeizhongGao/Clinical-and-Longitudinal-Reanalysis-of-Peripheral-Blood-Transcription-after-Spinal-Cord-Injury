# Local release checklist

## Prepared in this bundle

- Actual source scripts, source matrices and explicitly labeled intermediate inputs.
- Metadata curation, participant overlap/exclusion rules, gene lists and seeds.
- Main and supplementary numerical results, including negative and uncertain findings.
- A runnable entry point that stages an isolated workspace and stops on failed steps.
- R environment lock, Python dependency lists and restore instructions.
- Original/public file mapping, file checksums, data provenance and third-party notices.
- Numerical test records and stated limits; no invented permanent identifiers.

## Author decisions deliberately left open

- Select a license for author-owned original code. The author explicitly deferred it.
- Select repository ownership and confirm the final public release name/version.
- After actual publication, verify anonymous access and the archived version/DOI.
- Insert the real URL/version/DOI in manuscript and reviewer response R2-M19.

Preparing a folder does not fulfil the permanent-public-deposition part of R2-M19.
Do not state that this is already an open-source public archive while license and
publication decisions are pending.

## Upload guidance

Use this directory as the Git repository root, respecting `.gitignore`. A ZIP in
`upload_archives` is provided for a general-purpose archive. It contains the same
release files, excluding installed environments, run workspaces and recursively nested
archives. Do not upload the entire parent SCI project, original reviewer correspondence,
old Word drafts, caches or unrelated withdrawn-analysis outputs.
