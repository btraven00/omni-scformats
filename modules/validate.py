#!/usr/bin/env python3
"""Fidelity check of one conversion: input and output are each read by their
format's own reader (anndata via snapshot.py; readRDS + SeuratObject, or rhdf5
+ BPCells for .h5seurat, via snapshot.R) into a neutral snapshot, and
compare.py checks them property by property (names, obs/var columns with type
class, values and NA, every matrix: exact if integer-valued, else rtol 1e-6).

Cross-check: `scx diff --json input output` is recorded next to it, with
whether its overall verdict agrees. scx is never the judge: a disagreement is
a bug in one of the two, to be looked at.

status: valid (no property failed), invalid, not_run (empty output: the
conversion failed upstream), error (a reference reader failed on the input).

The original (crosscell validate_output.py) checked n_cells x n_genes only.
"""
import json, os, shutil, subprocess, sys, tempfile
from pathlib import Path

here = Path(__file__).parent
sys.path.insert(0, str(here))
import compare
from ob import parse


def snap(f, d):
    """Snapshot f into d with its reference reader; -> error text or None."""
    cmd = ([sys.executable, str(here / "snapshot.py")] if f.suffix == ".h5ad"
           else ["Rscript", "--vanilla", str(here / "snapshot.R")]) + [str(f), str(d)]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
    return None if r.returncode == 0 else (r.stdout + r.stderr)[-2000:]


def scx_diff(a, b):
    scx = os.environ.get("SCX_BIN") or shutil.which("scx")
    if not scx:
        return {"exit": None, "error": "scx not found (SCX_BIN)"}
    r = subprocess.run([scx, "diff", "--json", str(a), str(b)], capture_output=True, text=True, timeout=3600)
    try:
        return {"exit": r.returncode, "report": json.loads(r.stdout)}
    except json.JSONDecodeError:
        return {"exit": r.returncode, "error": (r.stdout + r.stderr)[-1000:]}


one, _, out = parse()
declared = Path(one.get("out_h5ad") or one["out_seurat"])
f, src = declared.resolve(), Path(one.get("rds") or one["h5ad"]).resolve()
rec = {"test_id": src.name.rsplit(".", 1)[0], "tool": declared.stem, "output": f.name}

if f.stat().st_size == 0:
    rec["status"] = "not_run"
else:
    # Snapshots hold every matrix in float64: next to the output, removed after.
    with tempfile.TemporaryDirectory(dir=Path(out["check"]).parent) as t:
        err = snap(src, Path(t) / "in")
        if err:
            rec.update(status="error", error="reference reader failed on the input: " + err)
        elif (err := snap(f, Path(t) / "out")):
            rec.update(status="invalid", error="reference reader failed on the output: " + err)
        else:
            props = compare.compare(compare.load(Path(t) / "in"), compare.load(Path(t) / "out"))
            fails = [p["property"] for p in props if p["status"] == "fail"]
            rec.update(status="invalid" if fails else "valid", failed=fails, properties=props)
    rec["scx_diff"] = scx_diff(src, f)
    if rec["status"] in ("valid", "invalid") and rec["scx_diff"].get("exit") in (0, 1):
        rec["scx_diff_agrees"] = (rec["scx_diff"]["exit"] == 0) == (rec["status"] == "valid")

Path(out["check"]).write_text(json.dumps(rec, indent=1) + "\n")
print(rec["test_id"], rec["tool"], rec["status"], rec.get("failed", ""))
