#!/usr/bin/env bash
# Steps 4-5 of extended/PLAN.md: gates, negative control, then plan-extended.yaml with the
# conversions one at a time (--cores 2, as seq_timing.sh) and the rest at --cores 16.
#   READERS=<readers env> scripts/run_extended.sh
set -euo pipefail
cd "$(dirname "$0")/.."
gen="$CROSSCELL_DATA/data/generated"
python3 - "$gen/extended_inputs.json" <<'PY'
import json, sys
m = json.load(open(sys.argv[1]))
bad = [k for k, v in m.items() if k.endswith(".rds") and not v.get("gate")]
print("gates:", "all pass" if not bad else f"FAIL {bad}")
sys.exit(1 if bad else 0)
PY
PATH="$READERS/bin:$PATH" python scripts/negctl.py
export PLAN=plan-extended.yaml
common=(--out-dir out-ext)
for stage in TO_H5AD TO_SEURAT; do
  scripts/run.sh --until "$stage" "${common[@]}" --cores 2 -- --conda-prefix "$PWD/out/.snakemake/conda" --rerun-triggers mtime
done
scripts/run.sh "${common[@]}" --cores 16 -- --conda-prefix "$PWD/out/.snakemake/conda" --rerun-triggers mtime
python3 scripts/verdict.py out-ext/SUMMARY.tables.default/results.tsv || true
echo "=== extended done $(date)"
