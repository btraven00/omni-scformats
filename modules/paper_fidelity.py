#!/usr/bin/env python3
"""The CrossCell paper's matrix fidelity, as published
(benchmark/notebooks/compute_fidelity.ipynb): for one RDS test case, each
tool's h5ad X against crosscell's h5ad X as the reference, by position, with
their 7 metrics; crosscell's own row is its h5ad -> rds -> h5ad round trip.
get_X_matrix, compute_metrics and safe_read_h5ad are copied verbatim.

Reported next to our fidelity check, not instead of it: the reference is one
of the tools, not the input, and cells/genes are matched by position, so a
converter that renames them (or is the reference) is not penalised.

Gather over out_h5ad, grouped by RDS case: --out_h5ad <one file per tool>.
"""
import json, os, subprocess, sys, tempfile, warnings
from pathlib import Path

import anndata as ad
import numpy as np
import scipy.sparse as sp
from scipy.stats import pearsonr, spearmanr

sys.path.insert(0, str(Path(__file__).parent))
from ob import parse

warnings.filterwarnings("ignore")


def get_X_matrix(adata):
    """Extract expression matrix from AnnData, handling X=None (anndataR puts data in layers)."""
    if adata.X is not None:
        return adata.X
    if adata.layers:
        for key in ['counts', 'data', 'X']:
            if key in adata.layers:
                return adata.layers[key]
        first_key = list(adata.layers.keys())[0]
        return adata.layers[first_key]
    return None


def compute_metrics(X1, X2, max_n=500000):
    """Compute 7 fidelity metrics between two expression matrices."""
    if X1.shape != X2.shape:
        return None, f'shape mismatch {X1.shape} vs {X2.shape}'
    total = X1.shape[0] * X1.shape[1]
    if total <= max_n:
        f1 = np.asarray(X1.todense() if sp.issparse(X1) else X1).flatten()
        f2 = np.asarray(X2.todense() if sp.issparse(X2) else X2).flatten()
    else:
        np.random.seed(42)
        ri = np.random.randint(0, X1.shape[0], max_n)
        ci = np.random.randint(0, X1.shape[1], max_n)
        if sp.issparse(X1):
            f1 = np.array([X1[r, c] for r, c in zip(ri, ci)]).flatten()
        else:
            f1 = X1[ri, ci].flatten()
        if sp.issparse(X2):
            f2 = np.array([X2[r, c] for r, c in zip(ri, ci)]).flatten()
        else:
            f2 = X2[ri, ci].flatten()
    if np.std(f1) == 0 or np.std(f2) == 0:
        pr = 1.0 if np.allclose(f1, f2) else 0.0
        sr = pr
    else:
        pr, _ = pearsonr(f1, f2)
        sr, _ = spearmanr(f1, f2)
    diff = np.abs(f1 - f2)
    mse = float(np.mean(diff**2))
    exact = np.sum(np.isclose(f1, f2, rtol=1e-7, atol=1e-10))
    return {
        'pearson_r': float(pr), 'spearman_rho': float(sr),
        'mse': mse, 'rmse': float(np.sqrt(mse)),
        'max_diff': float(np.max(diff)), 'mean_diff': float(np.mean(diff)),
        'exact_match_pct': float(exact / len(f1) * 100),
    }, None


def safe_read_h5ad(h5ad_path, tool_name=''):
    """Read H5AD file, handling non-standard formats from different tools."""
    import h5py
    try:
        with h5py.File(h5ad_path, 'r') as f:
            keys = list(f.keys())
    except Exception as e:
        return None, f'not valid HDF5: {e}'
    if 'X' not in keys and 'assay' in keys:
        try:
            from easySCFpy import loadH5
            adata = loadH5(h5ad_path)
            try:
                if adata.raw is not None:
                    adata = ad.AnnData(X=adata.raw.X)
            except Exception:
                pass
            return adata, None
        except Exception as e:
            return None, f'easySCF loadH5 error: {e}'
    if 'X' not in keys and 'obs' not in keys:
        return None, f'non-standard H5 format, keys={keys}'
    try:
        adata = ad.read_h5ad(h5ad_path)
        return adata, None
    except Exception as e:
        return None, str(e)


def metrics(ref_path, path, tool):
    a_ref, err = safe_read_h5ad(ref_path, 'CrossCell')
    if err:
        return {"status": "ref_read_error", "error": err[:300]}
    a, err = safe_read_h5ad(path, tool)
    if err:
        return {"status": "read_error", "error": err[:300]}
    X_ref, X = get_X_matrix(a_ref), get_X_matrix(a)
    if X_ref is None or X is None:
        return {"status": "no_X"}
    m, err = compute_metrics(X_ref, X)
    return {"status": "success", **m} if m else {"status": "metrics_error", "error": err}


def roundtrip(cc_h5ad):
    """crosscell h5ad -> rds -> h5ad, compared with the first h5ad (their section 5)."""
    cc = os.environ.get("CROSSCELL_BIN", "crosscell")
    with tempfile.TemporaryDirectory() as t:
        rds, back = f"{t}/rt.rds", f"{t}/rt.h5ad"
        for i, o, f in ((cc_h5ad, rds, "seurat"), (rds, back, "anndata")):
            r = subprocess.run([cc, "convert", "-i", i, "-o", o, "-f", f], capture_output=True, text=True, timeout=600)
            if r.returncode != 0:
                return {"status": "roundtrip_failed", "error": r.stderr[-300:]}
        return {**metrics(cc_h5ad, back, "CrossCell"), "test_type": "roundtrip"}


one, many, out = parse()
files = {Path(p).stem: p for p in many.get("out_h5ad", [])}
rows = []
cc = files.get("crosscell")
case = Path(cc or next(iter(files.values()))).resolve().parent.parent.name.split(".")[1]  # RDS.<case>.<hash>
for tool, p in sorted(files.items()):
    if Path(p).stat().st_size == 0:
        res = {"status": "failed"}  # conversion failed upstream
    elif cc is None or Path(cc).stat().st_size == 0:
        res = {"status": "no_reference"}
    elif tool == "crosscell":
        res = roundtrip(cc)
    else:
        res = {**metrics(cc, p, tool), "test_type": "cross_tool"}
    rows.append({"test_id": case, "tool": tool, **res})

dst = out.get("paper_fidelity") or str(Path(one["output_dir"]) / f"{one.get('name', case)}.paper_fidelity.json")
Path(dst).write_text(json.dumps(rows, indent=1) + "\n")
print(case, [(r["tool"], r["status"]) for r in rows])
