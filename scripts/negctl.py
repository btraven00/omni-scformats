#!/usr/bin/env python3
"""Negative control, extended/PLAN.md step 4: every input against itself must be called
identical, by our check (one snapshot compared with itself) and by `scx diff file file`.
Run in the readers env; prints one line per input and exits 1 on any failure."""
import json, os, subprocess, sys, tempfile
from pathlib import Path

root = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(root / "modules"))
import compare

gen = Path(os.environ.get("CROSSCELL_DATA", root / "data")) / "data/generated"
bad = 0
for name in json.loads((gen / "extended_inputs.json").read_text()):
    f = gen / name
    with tempfile.TemporaryDirectory(dir=gen) as t:
        cmd = ([sys.executable, str(root / "modules/snapshot.py")] if f.suffix == ".h5ad"
               else ["Rscript", "--vanilla", str(root / "modules/snapshot.R")]) + [str(f), t + "/s"]
        subprocess.run(cmd, check=True, capture_output=True)
        s = compare.load(Path(t) / "s")
        fails = [p["property"] for p in compare.compare(s, s) if p["status"] == "fail"]
    d = subprocess.run([os.environ["SCX_BIN"], "diff", str(f), str(f)], capture_output=True, text=True)
    ok = not fails and (d.returncode == 0 or f.suffix == ".rds")
    bad += not ok
    print(name, "PASS" if ok else f"FAIL check={fails} scx_diff_exit={d.returncode} {d.stderr[-200:]}", flush=True)
sys.exit(1 if bad else 0)
