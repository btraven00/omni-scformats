#!/usr/bin/env python3
"""Check one converted file's n_cells x n_genes against its input's, both read
by the format's own reader: anndata for .h5ad, readRDS + SeuratObject for
.rds, BPCells::open_matrix_hdf5 for scx's BPCells-backed .h5seurat.

As crosscell's validate_output.py, except: (1) the expectation is the input
file read the same way, not the test matrix (wrong for 4 of its 13 h5ad
cases); (2) no `crosscell inspect` first, so no converter grades output.
An empty output is a failed conversion upstream: status not_run."""
import json, subprocess, sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from ob import parse

READ = {
    ".h5ad": [sys.executable, "-c", "import anndata,sys; a=anndata.read_h5ad(sys.argv[1]); print(a.n_obs, a.n_vars)"],
    ".h5seurat": ["Rscript", "-e", "f <- commandArgs(TRUE)[1]; l <- rhdf5::h5ls(f); "
                  "g <- with(l, paste0(group, '/', name))[grepl('^/assays/[^/]+$', l$group) & l$name == 'counts'][1]; "
                  "m <- BPCells::open_matrix_hdf5(f, sub('^/', '', g)); cat(ncol(m), nrow(m))"],
    ".rds": ["Rscript", "-e", "suppressMessages(library(SeuratObject)); o <- readRDS(commandArgs(TRUE)[1]); "
             "if (!inherits(o, 'Seurat')) stop('not a Seurat object'); cat(ncol(o), nrow(o))"],
}


def dims(f):
    """[n_cells, n_genes], or the reader's error text."""
    cmd = READ.get(f.suffix, READ[".rds"]) + [str(f)]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
    return [int(x) for x in r.stdout.split()[-2:]] if r.returncode == 0 else r.stderr[-2000:]


one, _, out = parse()
declared = Path(one.get("out_h5ad") or one["out_seurat"])
f, src = declared.resolve(), Path(one.get("rds") or one["h5ad"])
rec = {"test_id": src.name.rsplit(".", 1)[0], "tool": declared.stem, "output": f.name}
if f.stat().st_size == 0:
    rec["status"] = "not_run"
else:
    want, got = dims(src), dims(f)
    if isinstance(want, str):
        rec.update(status="error", error="reference reader failed on the input: " + want)
    elif isinstance(got, str):
        rec.update(expected=want, status="invalid", error=got)
    else:
        rec.update(expected=want, got=got, status="valid" if got == want else "invalid")
Path(out["check"]).write_text(json.dumps(rec, indent=1) + "\n")
print(rec["test_id"], rec["tool"], rec["status"])
