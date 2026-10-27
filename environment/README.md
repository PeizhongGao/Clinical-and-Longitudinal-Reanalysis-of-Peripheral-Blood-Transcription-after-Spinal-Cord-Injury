# Recorded computational environment

The tested host uses Windows 11, R 4.6.1 and Python 3.12.14. `renv.lock` records
102 direct/transitive non-base R packages and their source metadata. In particular,
the fixed annotation and numerical pipeline require org.Hs.eg.db 3.23.1, limma
3.68.4, edgeR 4.10.1, lmerTest 3.2-1 and emmeans 2.0.3.

`Rscript --vanilla environment/restore.R .` installs into `environment/r_library`.
It may need internet access and compiler/system libraries when binary packages are
unavailable. The entry point automatically uses this library if it exists.
Do not silently substitute the latest annotation database: it can change symbol
mapping and GO membership. The original fixed definitions are also included.

Standard Python dependencies are pinned in `requirements.txt`. An isolated virtual
environment is recommended. The optional extended GO script uses gseapy 1.3.1;
`requirements-extended.txt` requests that package and its dependencies. It no longer
loads a Windows-specific extension by the author's personal filesystem path.

The recorded environment and the relocated-package tests are not a claim that a
fresh installation has been verified on every OS. No Linux/macOS execution or
full clean-library download/restore has been certified. Arial/Cairo/font libraries
affect visual rendering; numerical validation is independent of pixel-identical
font rendering. PNG/PDF/SVG/TIFF reference figures are provided.

R may emit Windows locale startup warnings; these were observed on the tested host.
Actual stage failures stop the runner and are recorded in the logs.
