#!/usr/bin/env bash
# Fetch crosscell's benchmark data with hapiq and unpack it to
# $CROSSCELL_DATA/data/generated/, where modules/data.py looks.
#   scripts/fetch.sh                                   # all 5 archives, 6.4 GB
#   scripts/fetch.sh h5ad_small seurat_v5              # a slice
# Archives: seurat_v4 seurat_v5 h5ad_small h5ad_spatial h5ad_cellxgene
#
# Pinned to the version record (18610877), not the concept DOI (18610876, which
# follows the latest version), and to the md5s it published on 2026-02-11:
# hapiq checks Zenodo's md5 on its own, --hash makes a new upload fail here.
set -euo pipefail
dest="${CROSSCELL_DATA:-$(cd "$(dirname "$0")/.." && pwd)/data}"
hapiq="${HAPIQ:-$HOME/phd/hapiq/bin/hapiq}"
record=10.5281/zenodo.18610877
declare -A md5=(
  [manifest.json]=23c39874fc69a98c58786762a0a62f72
  [seurat_v4_datasets.tar.gz]=e82e7876130ff68abe45a92bf42fd30a
  [seurat_v5_datasets.tar.gz]=86fa98b5f42eeb16b1ca06bd685ada0b
  [h5ad_small_datasets.tar.gz]=e7663b581b4991a7ae35d36134191efd
  [h5ad_spatial_datasets.tar.gz]=1457638aab62870e8e80ef2867a41dfb
  [h5ad_cellxgene_datasets.tar.gz]=d24bae3c6ade178f3e443441a7dfdc69
)
[[ $# -gt 0 ]] || set -- seurat_v4 seurat_v5 h5ad_small h5ad_spatial h5ad_cellxgene
mkdir -p "$dest/zenodo"
# --filename-pattern is a single Go glob (no braces), so one call per file,
# which is also what --hash needs.
for f in manifest.json "${@/%/_datasets.tar.gz}"; do
  "$hapiq" download zenodo "$record" --out "$dest/zenodo" -y --skip-existing \
    --filename-pattern "$f" --hash "md5:${md5[$f]}"
done
for f in "${@/%/_datasets.tar.gz}"; do
  tar -xzf "$dest/zenodo/$f" -C "$dest"   # members are data/generated/<file>
done
# sha256 of each file is checked against manifest.json when it is staged.
ls "$dest/data/generated" | wc -l
