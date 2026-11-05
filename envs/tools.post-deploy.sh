#!/usr/bin/env bash
# The two GitHub-only converters, at main HEAD as of 2026-11-05 (unchanged
# since the paper's run). Missing R deps come from current CRAN.
set -euo pipefail
Rscript --vanilla -e '
  options(repos = c(CRAN = "https://cloud.r-project.org"), Ncpus = 8)
  remotes::install_github("settylab/convert2anndata@550bfb07fc", upgrade = "never")
  remotes::install_github("xleizi/easySCF/r@0a5a667720", upgrade = "never")
  for (p in c("Seurat", "zellkonverter", "anndataR", "convert2anndata", "easySCFr"))
    cat(p, as.character(packageVersion(p)), "\n")
'
