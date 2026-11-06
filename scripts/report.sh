#!/usr/bin/env bash
# Render report/report.html from report/results.tsv (+ tables.md) with Quarto.
#   scripts/report.sh                 # after copying a SUMMARY's outputs into report/
# Figures are redrawn from results.tsv, so edits to modules/report.py need no rerun.
set -euo pipefail
cd "$(dirname "$0")/../report"
exec uv run --no-project -q --python 3.12 --with jupyter,matplotlib,scipy,"argcomplete<3.6" quarto render report.qmd
