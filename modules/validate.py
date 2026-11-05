#!/usr/bin/env python3
"""Check one converted file's n_cells x n_genes, as crosscell's
validate_output.py: anndata.read_h5ad for .h5ad, readRDS for .rds. Deviation:
the original tried `crosscell inspect` first for RDS; here R's own reader is
the only judge, so crosscell never grades its own output.
An empty file is a failed conversion upstream: status not_run."""
import json, subprocess, sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from ob import parse

one, _, out = parse()
f = Path(one.get("out_h5ad") or one["out_seurat"])
exp = json.loads(Path(one.get("rds_expected") or one["h5ad_expected"]).read_text())
rec = {"test_id": exp["test_id"], "tool": f.stem, "expected": [exp["n_cells"], exp["n_genes"]]}

if f.stat().st_size == 0:
    rec["status"] = "not_run"
else:
    if f.suffix == ".h5ad":
        cmd = [sys.executable, "-c", "import anndata,sys; a=anndata.read_h5ad(sys.argv[1]); print(a.n_obs, a.n_vars)", str(f)]
    else:
        cmd = ["Rscript", "-e", "suppressMessages(library(SeuratObject)); o <- readRDS(commandArgs(TRUE)[1]); "
               "if (!inherits(o, 'Seurat')) stop('not a Seurat object'); cat(ncol(o), nrow(o))", str(f)]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
    if r.returncode != 0:
        rec.update(status="invalid", error=r.stderr[-2000:])
    else:
        got = [int(x) for x in r.stdout.split()[-2:]]
        rec.update(got=got, status="valid" if got == rec["expected"] else "invalid")
Path(out["check"]).write_text(json.dumps(rec, indent=1) + "\n")
print(rec["test_id"], rec["tool"], rec["status"])
