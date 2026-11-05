#!/usr/bin/env bash
# ob run, plus the conda post-deploy scripts ob does not carry over.
#   scripts/run.sh [ob run flags...]      e.g. --filter slice.obfilter.yaml --cores 8
# ob copies envs/<x>.yml to out/.envs/<id>.yaml and nothing else, so snakemake
# never finds envs/<x>.post-deploy.sh. Seurat 5.4.0 (no R 4.5 conda build) and
# the GitHub-only converters live in that script, so copy it next to the copied
# yaml before snakemake creates the env.
set -euo pipefail
cd "$(dirname "$0")/.."
# The venv's ob binary, not `uv run`: an activated venv leaks into every rule
# (reticulate and `python3` then pick the venv's python, not the conda env's).
OB="${OB:-$HOME/lab/omnibenchmark/.venv/bin/ob}"
unset VIRTUAL_ENV
$OB run benchmark.yaml --dry --dirty "$@" >/dev/null
out=out; prev=""; for a in "$@"; do [[ "$prev" == "--out-dir" ]] && out="$a"; prev="$a"; done
cp envs/tools.post-deploy.sh "$out/.envs/tools.post-deploy.sh"
exec $OB run benchmark.yaml --dirty -k "$@"
