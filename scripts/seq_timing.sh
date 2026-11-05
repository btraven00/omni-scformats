#!/usr/bin/env bash
# Re-time every conversion one at a time, then put those timings in out/.
#   scripts/seq_timing.sh           (after a full run in out/; envs and data reused)
# The concurrent run in out/ keeps its outputs and checks; only the perf records
# (perf.json + .runs/) are replaced by the sequential ones, and SUMMARY rebuilt.
# Conversion outputs don't depend on timing, so the checks stay valid.
set -euo pipefail
cd "$(dirname "$0")/.."
# Data + one conversion stage per pass (--until keeps a stage and its ancestors;
# an obfilter would keep the gather stages and fail on their pruned inputs).
for stage in TO_H5AD TO_SEURAT; do
  scripts/run.sh --until "$stage" --out-dir out-seq --cores 2 -- \
    --conda-prefix "$PWD/out/.snakemake/conda" --rerun-triggers mtime
done
for d in out-seq/*/TO_*; do
  dst="out/${d#out-seq/}"
  [[ -d "$dst" ]] || { echo "no counterpart for $d in out/" >&2; exit 1; }
  cp "$d"/*.perf.json "$dst"/
  rm -rf "$dst/.runs" && cp -r "$d/.runs" "$dst/.runs"
done
rm -rf out/SUMMARY.*
scripts/run.sh --cores 8 -- --rerun-triggers mtime
