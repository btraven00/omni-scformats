#!/usr/bin/env python3
"""Snapshot an .h5ad with anndata, the format's own reader, into the neutral
layout snapshot.R also writes, so compare.py can check R- and Python-read
files property by property:

  obs_names.txt  var_names.txt                 one name per line
  obs/<k>.* var/<k>.*  per column k (0-based, column order):
      .name  .class (int|float|bool|string|categorical)  .na (uint8, 1 = NA)
      .f64 (float64 LE; int/float/bool)  .txt (string)
      .codes (int32 LE, 0-based, -1 = NA) .levels.txt .ordered (categorical)
  mats/<k>.*  per matrix: .name .shape ("n_obs n_vars") .dtype .var_names.txt
      .indptr (float64) .indices (float64) .data (float64): CSR, cells x genes

    snapshot.py <file.h5ad> <out_dir>
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp


def lines(p, xs):
    Path(p).write_text("".join(f"{x}\n" for x in xs))


def column(d, k, name, s):
    p = d / str(k)
    lines(f"{p}.name", [name])
    na = s.isna().to_numpy()
    na.astype(np.uint8).tofile(f"{p}.na")
    dt = s.dtype
    if isinstance(dt, pd.CategoricalDtype):
        cls = "categorical"
        s.cat.codes.to_numpy().astype("<i4").tofile(f"{p}.codes")
        lines(f"{p}.levels.txt", [str(c) for c in s.cat.categories])
        lines(f"{p}.ordered", [int(dt.ordered)])
    elif pd.api.types.is_bool_dtype(dt):
        cls = "bool"
    elif pd.api.types.is_integer_dtype(dt):
        cls = "int"
    elif pd.api.types.is_float_dtype(dt):
        cls = "float"
    else:
        cls = "string"
        lines(f"{p}.txt", ["" if m else str(v) for v, m in zip(s, na)])
    if cls in ("int", "float", "bool"):
        v = s.astype("float64").to_numpy(dtype="<f8", copy=True)
        v[na] = 0
        v.tofile(f"{p}.f64")
    lines(f"{p}.class", [cls])


def matrix(d, k, name, m, var_names):
    p = d / str(k)
    dtype = str(m.dtype)
    m = sp.csr_matrix(m)
    m.sort_indices()
    lines(f"{p}.name", [name])
    lines(f"{p}.shape", [f"{m.shape[0]} {m.shape[1]}"])
    lines(f"{p}.dtype", [dtype])
    lines(f"{p}.var_names.txt", var_names)
    for part in ("indptr", "indices", "data"):
        getattr(m, part).astype("<f8").tofile(f"{p}.{part}")


def main(f, out):
    import anndata
    write(anndata.read_h5ad(f), out)


def write(a, out):
    out = Path(out)
    for sub in ("obs", "var", "mats"):
        (out / sub).mkdir(parents=True, exist_ok=True)
    lines(out / "obs_names.txt", a.obs_names)
    lines(out / "var_names.txt", a.var_names)
    for slot, df in (("obs", a.obs), ("var", a.var)):
        for k, c in enumerate(df.columns):
            column(out / slot, k, str(c), df[c])
    mats = [("X", a.X, a.var_names)] if a.X is not None else []
    # anndata >= 0.13 lists X as layers[None]; it is already "X".
    mats += [(f"layers/{n}", a.layers[n], a.var_names) for n in a.layers if n is not None]
    if a.raw is not None:
        mats.append(("raw/X", a.raw.X, a.raw.var_names))
    for k, (n, m, v) in enumerate(mats):
        matrix(out / "mats", k, n, m, v)


if __name__ == "__main__":
    main(*sys.argv[1:3])
