#!/usr/bin/env bash
# Step 3 of extended/PLAN.md on a big host: CELLxGENE inputs 3 at a time, the
# ladder in its own lane (one subsample pass), then the manifest.
#   READERS=<readers conda env> scripts/build_extended.sh
set -euo pipefail
cd "$(dirname "$0")/.."
export PATH="$READERS/bin:$PATH" OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2
logs="${CROSSCELL_DATA:-data}/logs/build"; mkdir -p "$logs"
b() { python scripts/build_extended.py one "$1" > "$logs/$1.log" 2>&1 && echo "done $1" || echo "FAILED $1"; }
export -f b; export logs
python scripts/build_extended.py list | grep '^ext_' | xargs -P 3 -I{} bash -c 'b {}' &
python scripts/build_extended.py list | grep '^ladder_' | while read -r n; do b "$n"; done &
wait
python scripts/build_extended.py manifest > /dev/null
echo "=== build done $(date)"
