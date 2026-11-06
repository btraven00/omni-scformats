#!/usr/bin/env python3
"""Global gather: perf.json, check.json and the paper-fidelity JSON of every
node -> results.tsv (one row per test case x tool) and, via report.py,
tables.md + tools/phases figures. Converted = exited 0 on every run, as the
original; faithful = passed every property of our check."""
import csv, json, sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import report
from ob import parse

one, many, out = parse()
# The gather rule passes no --output (omnibenchmark HEAD): write next to output_dir.
odir = Path(one["output_dir"])

rows = {}
for p in many.get("perf", []):
    r = json.loads(Path(p).read_text())
    rows[(r["test_id"], r["tool"])] = r
for p in many.get("check", []):
    c = json.loads(Path(p).read_text())
    rows[(c["test_id"], c["tool"])].update(check=c["status"], failed=",".join(c.get("failed", [])),
                                           properties=c.get("properties", []),
                                           failed_detail=json.dumps({x["property"]: x.get("detail", "")
                                                                     for x in c.get("properties", [])
                                                                     if x["status"] == "fail"}),
                                           scx_diff_agrees=c.get("scx_diff_agrees", ""))
paper = [r for p in many.get("paper_fidelity", []) for r in json.loads(Path(p).read_text())]

cols = ["test_id", "tool", "direction", "status", "check", "failed", "failed_detail", "scx_diff_agrees", "conversion_time_seconds",
        "peak_memory_mb", "output_file_size_mb", "all_times", "all_mems", "phases", "output", "error"]
with (odir / "results.tsv").open("w", newline="") as fh:
    w = csv.DictWriter(fh, cols, delimiter="\t", extrasaction="ignore")
    w.writeheader()
    for k in sorted(rows):
        err = (rows[k].get("error") or "").strip().splitlines()
        w.writerow({**rows[k], "error": err[-1] if err else "", "phases": json.dumps(rows[k].get("phases", {}))})

report.write(rows, paper, odir)
print((odir / "tables.md").read_text())
