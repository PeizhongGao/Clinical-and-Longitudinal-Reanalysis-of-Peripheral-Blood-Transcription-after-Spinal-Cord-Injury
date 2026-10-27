# Seven redrawn scientific figures

This directory includes R code, twenty direct source CSVs, seven vector PDFs,
seven SVGs and seven 600-dpi, RGB, LZW-compressed TIFFs. Width is capped at 170 mm.
Arial, cohort colors, panel typography and statistical boundaries are consistent
across the set. English captions are in `Figure_legends.md`.

From the repository root, render into a new output directory:

```sh
python reproduce_revision.py --task figures --output outputs/redrawn-figures
```

The portable `make_figures.R` uses its first argument as a directory containing
`source_data/`, and writes PDF/TIFF/SVG and QA previews there. Use the root runner
to protect the published reference images. It requires data.table, ggplot2,
patchwork, jsonlite, scales, svglite and ragg in the recorded R environment.

Only plotting is rerun; saved estimates, uncertainty intervals and group counts
are not refitted. The original project-specific extraction and verification
scripts are under `provenance/` for traceability, not portable entry points.
The CSV snapshots allow rendering without the author's original project tree.

These are redrawn versions of the original Figure 1-7 set, not final submission
figure numbers. The current formal manuscript uses a separate BTM main figure;
see `docs/REVISION_EXTENSIONS.md` for the numbering boundary.
