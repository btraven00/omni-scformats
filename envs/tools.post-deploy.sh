#!/usr/bin/env bash
# Seurat 5.4.0 (no R 4.5 conda build) and the GitHub-only converters, from a
# CRAN snapshot dated to crosscell's benchmark run (results/*.json: 2026-02-12).
set -euo pipefail
Rscript --vanilla -e '
  options(repos = c(CRAN = "https://packagemanager.posit.co/cran/2026-02-12"), Ncpus = 8)
  remotes::install_version("Seurat", "5.4.0", upgrade = "never")
  # crosscell installed both from main; pinned to main as of 2026-02-12.
  remotes::install_github("settylab/convert2anndata@550bfb07fc", upgrade = "never")
  remotes::install_github("xleizi/easySCF/r@0a5a667720", upgrade = "never")
  for (p in c("Seurat", "zellkonverter", "anndataR", "convert2anndata", "easySCFr"))
    cat(p, as.character(packageVersion(p)), "\n")
'
