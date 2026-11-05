#!/usr/bin/env python3
"""Global gather: every perf.json and check.json -> results.tsv (one row per
test case x tool) and pass_rates.md (the crosscell README table, plus a
validated column). Pass = conversion exited 0 on every run, as the original."""
import csv, json, sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from ob import parse

one, many, out = parse()
# The gather rule passes no --output (omnibenchmark HEAD); fall back to the plan's filenames.
res = Path(out.get("results_tsv", Path(one["output_dir"]) / "results.tsv"))
md = Path(out.get("pass_rates_md", Path(one["output_dir"]) / "pass_rates.md"))

rows = {}
for p in many.get("perf", []):
    r = json.loads(Path(p).read_text())
    rows[(r["test_id"], r["tool"])] = r
for p in many.get("check", []):
    c = json.loads(Path(p).read_text())
    rows[(c["test_id"], c["tool"])].update(check=c["status"], failed=",".join(c.get("failed", [])),
                                           scx_diff_agrees=c.get("scx_diff_agrees", ""))

cols = ["test_id", "tool", "direction", "status", "check", "failed", "scx_diff_agrees", "conversion_time_seconds",
        "peak_memory_mb", "output_file_size_mb", "all_times", "all_mems", "phases", "output", "error"]
with res.open("w", newline="") as fh:
    w = csv.DictWriter(fh, cols, delimiter="\t", extrasaction="ignore")
    w.writeheader()
    for k in sorted(rows):
        err = (rows[k].get("error") or "").strip().splitlines()
        w.writerow({**rows[k], "error": err[-1] if err else "", "phases": json.dumps(rows[k].get("phases", {}))})

tools = ["crosscell", "zellkonverter", "anndataR", "convert2anndata", "easySCF", "scx"]
dirs = ["rds_to_h5ad", "h5ad_to_rds"]
n, ok, val = Counter(), Counter(), Counter()
for r in rows.values():
    k = (r["tool"], r["direction"])
    n[k] += 1
    ok[k] += r["status"] == "success"
    val[k] += r.get("check") == "valid"
cell = lambda c, t, ds: (lambda a, b: f"{a}/{b} ({100 * a // b}%)" if b else "N/A")(sum(c[t, d] for d in ds), sum(n[t, d] for d in ds))
lines = ["| Tool | rds_to_h5ad | h5ad_to_rds | Total | Faithful (all properties) |", "|---|---|---|---|---|"]
for t in tools:
    if any(n[t, d] for d in dirs):
        lines.append(f"| {t} | {cell(ok, t, dirs[:1])} | {cell(ok, t, dirs[1:])} | {cell(ok, t, dirs)} | {cell(val, t, dirs)} |")
md.write_text("\n".join(lines) + "\n")
print("\n".join(lines))
