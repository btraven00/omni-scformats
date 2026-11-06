#!/usr/bin/env python3
"""Mechanical verdicts for extended/PLAN.md §4 (C1 speed, C2 memory, C3 coverage).

    python3 scripts/verdict.py <SUMMARY dir of plan-extended.yaml>/results.tsv
"""
import csv, sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "modules"))
import report as R

rows = {}
for r in csv.DictReader(open(sys.argv[1]), delimiter="\t"):
    r["conversion_time_seconds"] = float(r["conversion_time_seconds"] or 0)
    rows[(r["test_id"], r["tool"])] = r
meta = R.ext_meta()

# C1: per bin x direction x other tool, scx faster on >= 80 % of datasets both converted
print("C1 speed: scx faster than each tool on >= 80 % of datasets both converted")
c1 = True
for b in ("S", "M", "L"):
    for d in R.DIRS:
        ok = defaultdict(dict)
        for (c, t), r in rows.items():
            if meta.get(c, {}).get("bin") == b and r["direction"] == d and r["status"] == "success":
                ok[t][c] = r["conversion_time_seconds"]
        for t in sorted(set(ok) - {"scx"}):
            both = [c for c in ok[t] if c in ok["scx"]]
            if not both:
                print(f"  {b} {d} vs {t}: no dataset both converted")
                continue
            wins = sum(ok["scx"][c] < ok[t][c] for c in both)
            p = wins / len(both) >= 0.8
            c1 &= p
            print(f"  {b} {d} vs {t}: {wins}/{len(both)} {'PASS' if p else 'FAIL'}")
print(f"C1 {'PASS' if c1 else 'FAIL'}\n")

# C2: ladder, scx peak-RSS log-log slope, 95 % CI upper bound < 1, per direction
print("C2 memory: scx ladder slope of log peak RSS on log cells, 95 % CI upper < 1")
s = R.ladder_slopes(rows, meta)
for d in R.DIRS:
    v = s.get((d, "scx", "rss"))
    print(f"  {d}: " + (f"slope {v[0]:.2f} [{v[1]:.2f}, {v[2]:.2f}] n={v[3]} {'PASS' if v[2] < 1 else 'FAIL'}"
                        if v else "fewer than 3 converted rungs: NOT EVALUABLE"))

# C3: coverage, reported as is
print("\nC3 coverage (reported, no criterion): scx converted / cases")
for b in ("S", "M", "L", "ladder"):
    for d in R.DIRS:
        rs = [r for (c, t), r in rows.items() if t == "scx" and r["direction"] == d and meta.get(c, {}).get("bin") == b]
        if rs:
            print(f"  {b} {d}: {sum(r['status'] == 'success' for r in rs)}/{len(rs)}")
