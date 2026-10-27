# Supplementary Data 9 — nested and resampling audit

This archive contains newly executed revision-stage analyses. It supersedes earlier statements that no nested or alternative resampling computation had been performed; it does not reinstate the historical nine-gene model or old P values. The original 15+15 and BTM results remain unchanged.

Start with PROTOCOL.md and analysis-report.md. Tables contain all observed predictions, split/candidate/selection audits, permutation labels and statistics, fixed-score null results and diagnostics. Scripts implement every new calculation; inputs are byte-preserved snapshots, with SHA-256 values in audit/protocol_input_lock.json. The historical_source.Rmd input is provenance only, not a current analysis entry point.

Run from this extracted directory with R 4.6.1 and the recorded packages:

    Rscript --vanilla scripts/run_nested.R . pilot
    Rscript --vanilla scripts/run_nested.R . observed
    Rscript --vanilla scripts/run_nested.R . permutations
    Rscript --vanilla scripts/run_clinical_resampling.R .
    Rscript --vanilla scripts/summarize_results.R .
    Rscript --vanilla scripts/make_figures.R .
    python scripts/verify_results.py

The permutation command is resumable and uses six local workers. Independent Python checks require numpy and pandas. Existing result tables are included; no execution is needed to inspect them. Classifier permutations assume sex-conditional label exchangeability; fixed-score residual permutations are approximate and use separate Holm4 families. Neither provides external clinical validation.

Source data: GSE151371 raw counts and archived normalized/ComBat matrix; public ODC-SCI dataset 405 (DOI 10.34945/F5QC7J, CC BY 4.0) and the derived public sample crosswalk. Clinical BTM inputs derive from the existing reanalysis with published BTM definitions. Source-specific attribution/licensing remains applicable; the authors have not selected a reuse license for original project code. No new blanket license is granted.

The prior public analysis commit remains 6415b0bdc3dff91d3949e940c363ed22f48f387c. This new extension is supplied in Data 9 and is not represented as already included in that commit or pushed to GitHub. Data 7 records the historical withdrawal; Data 8 covers the earlier reproducibility snapshot.
