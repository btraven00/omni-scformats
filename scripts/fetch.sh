#!/usr/bin/env bash
# Fetch crosscell's benchmark data (Zenodo 10.5281/zenodo.18610876) with hapiq
# and unpack it to $CROSSCELL_DATA/data/generated/, where modules/data.py looks.
#   scripts/fetch.sh                                   # all 5 archives, 6.4 GB
#   scripts/fetch.sh h5ad_small seurat_v5              # a slice
# Archives: seurat_v4 seurat_v5 h5ad_small h5ad_spatial h5ad_cellxgene
set -euo pipefail
dest="${CROSSCELL_DATA:-$(cd "$(dirname "$0")/.." && pwd)/data}"
hapiq="${HAPIQ:-$HOME/phd/hapiq/bin/hapiq}"
doi=10.5281/zenodo.18610876
[[ $# -gt 0 ]] || set -- seurat_v4 seurat_v5 h5ad_small h5ad_spatial h5ad_cellxgene
mkdir -p "$dest/zenodo"
# --filename-pattern is a single Go glob (no braces), so one call per file.
for f in manifest.json "${@/%/_datasets.tar.gz}"; do
  "$hapiq" download zenodo "$doi" --out "$dest/zenodo" -y --skip-existing --filename-pattern "$f"
done
for f in "${@/%/_datasets.tar.gz}"; do
  tar -xzf "$dest/zenodo/$f" -C "$dest"   # members are data/generated/<file>
done
# sha256 of each file is checked against manifest.json when it is staged.
ls "$dest/data/generated" | wc -l
