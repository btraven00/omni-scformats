#!/usr/bin/env python3
"""Compare two snapshots (snapshot.py / snapshot.R layout), input vs output,
property by property. Each property is pass, note or fail.

Fail: cell or gene names changed or reordered; an input obs/var column
missing; a column's type class changed (int/float/bool/string/categorical,
`ordered` included); values or NA mask changed; an input matrix with no
equal matrix in the output. Integer-valued matrices must match exactly, float
matrices within rtol 1e-6 (a f64 -> f32 downcast passes, and is noted).
Note: categorical level order, columns only in the output, float32 storage.
A passing matrix records where it landed (e.g. RNA/counts = X): slot names
differ by format, so the mapping is data for the report, not a verdict.

    compare.py <input_snapshot> <output_snapshot>   -> JSON list on stdout
"""
import json, sys
from pathlib import Path

import numpy as np
import scipy.sparse as sp

RTOL = 1e-6


def names(p):
    t = Path(p).read_text()
    return t.split("\n")[:-1] if t else []


def load(d):
    d = Path(d)
    snap = {"obs_names": names(d / "obs_names.txt"), "var_names": names(d / "var_names.txt")}
    for slot in ("obs", "var"):
        cols = {}
        for f in sorted((d / slot).glob("*.name"), key=lambda f: int(f.stem)):
            p = str(f)[: -len(".name")]
            cls = Path(p + ".class").read_text().strip()
            na = np.fromfile(p + ".na", dtype=np.uint8).astype(bool)
            if cls == "categorical":
                codes = np.fromfile(p + ".codes", dtype="<i4")
                lv = names(p + ".levels.txt")
                vals = np.array([lv[c] if c >= 0 else "" for c in codes], dtype=object)
                c = {"levels": lv, "ordered": Path(p + ".ordered").read_text().strip() == "1"}
            elif cls == "string":
                vals, c = np.array(names(p + ".txt"), dtype=object), {}
            else:
                vals, c = np.fromfile(p + ".f64", dtype="<f8"), {}
            cols[f.read_text().rstrip("\n")] = {"class": cls, "na": na, "values": vals, **c}
        snap[slot] = cols
    mats = {}
    for f in sorted((d / "mats").glob("*.name"), key=lambda f: int(f.stem)):
        p = str(f)[: -len(".name")]
        n_obs, n_vars = map(int, Path(p + ".shape").read_text().split())
        parts = [np.fromfile(f"{p}.{x}", dtype="<f8") for x in ("data", "indices", "indptr")]
        m = sp.csr_matrix((parts[0], parts[1].astype(np.int64), parts[2].astype(np.int64)), shape=(n_obs, n_vars))
        mats[f.read_text().rstrip("\n")] = {"m": m, "var_names": names(p + ".var_names.txt"),
                                            "dtype": Path(p + ".dtype").read_text().strip()}
    snap["mats"] = mats
    return snap


def r(prop, status, detail=""):
    return {"property": prop, "status": status, "detail": detail}


def same_names(prop, a, b):
    """-> (result, permutation of b onto a, or None to align by position)."""
    if a == b:
        return r(prop, "pass"), None
    if sorted(a) == sorted(b):
        pos = {x: i for i, x in enumerate(b)}
        return r(prop, "fail", "same names, different order"), np.array([pos[x] for x in a]) if len(set(a)) == len(a) else None
    diff = [(x, y) for x, y in zip(a, b) if x != y]
    return r(prop, "fail", f"{len(a)} vs {len(b)} names; {len(diff)} differ by position, e.g. {diff[:3]}"), None


def columns(slot, a, b, perm):
    out = []
    for name, ca in a.items():
        prop = f"{slot}/{name}"
        cb = b.get(name)
        if cb is None:
            out.append(r(prop, "fail", "column dropped"))
            continue
        if len(cb["values"]) != len(ca["values"]):
            out.append(r(prop, "fail", f"length {len(ca['values'])} -> {len(cb['values'])}"))
            continue
        take = (lambda x: x[perm]) if perm is not None else (lambda x: x)
        vb, nb = take(cb["values"]), take(cb["na"])
        if ca["class"] != cb["class"] or ca.get("ordered") != cb.get("ordered"):
            out.append(r(prop, "fail", f"type {ca['class']}{' ordered' if ca.get('ordered') else ''} -> "
                                       f"{cb['class']}{' ordered' if cb.get('ordered') else ''}"))
            continue
        if not np.array_equal(ca["na"], nb):
            out.append(r(prop, "fail", f"NA mask differs in {int((ca['na'] != nb).sum())} rows"))
            continue
        ok = ~ca["na"]
        if ca["class"] == "float":
            bad = ~np.isclose(vb[ok], ca["values"][ok], rtol=RTOL, atol=0, equal_nan=True)
        else:
            bad = vb[ok] != ca["values"][ok]
        if bad.any():
            i = np.flatnonzero(ok)[np.flatnonzero(bad)[:3]]
            out.append(r(prop, "fail", f"{int(bad.sum())} values differ, e.g. "
                                       f"{[(str(ca['values'][j]), str(vb[j])) for j in i]}"))
        elif ca["class"] == "categorical" and ca["levels"] != cb["levels"]:
            out.append(r(prop, "note", "level order differs" if sorted(ca["levels"]) == sorted(cb["levels"])
                         else f"levels {len(ca['levels'])} -> {len(cb['levels'])} (unused levels changed)"))
        else:
            out.append(r(prop, "pass"))
    extra = [n for n in b if n not in a]
    if extra:
        out.append(r(f"{slot}/+", "note", f"only in output: {extra}"))
    return out


def aligned(mb, rows, var_a):
    """mb's matrix with rows permuted onto a, columns onto var_a; None if impossible."""
    m = mb["m"] if rows is None else mb["m"][rows]
    if mb["var_names"] == var_a:
        return m
    if sorted(mb["var_names"]) == sorted(var_a) and len(set(var_a)) == len(var_a):
        pos = {x: i for i, x in enumerate(mb["var_names"])}
        return m[:, [pos[x] for x in var_a]]
    return m if m.shape[1] == len(var_a) else None  # names changed: by position


def equal(a, b):
    """-> (equal?, detail). Exact when a is integer-valued, else rtol."""
    if a.shape != b.shape:
        return False, f"shape {a.shape} vs {b.shape}"
    a, b = a.tocsr(), b.tocsr()
    a.sort_indices(); b.sort_indices()
    integral = np.array_equal(a.data, np.round(a.data))
    d = (a - b).tocsr()
    d.eliminate_zeros()
    if integral:
        return d.nnz == 0, f"{d.nnz} entries differ (exact, integer-valued)"
    over = (abs(d) - abs(a).maximum(abs(b)) * RTOL).tocsr()
    n = int((over.data > 0).sum())
    return n == 0, f"{n} entries differ beyond rtol {RTOL}"


def matrices(a, b, rows):
    out = []
    for name, ma in a["mats"].items():
        hits, misses = [], []
        for nb, mb in b["mats"].items():
            m = aligned(mb, rows, ma["var_names"])
            if m is None or m.shape != ma["m"].shape:
                continue
            ok, detail = equal(ma["m"], m)
            (hits if ok else misses).append(nb if ok else f"{nb}: {detail}")
        out.append(r(f"matrix/{name}", "pass", f"= {hits}") if hits else
                   r(f"matrix/{name}", "fail", "; ".join(misses) or "no output matrix of the same shape"))
    for nb, mb in b["mats"].items():
        if mb["dtype"] in ("float32",):
            out.append(r(f"matrix/{nb}.dtype", "note", "stored as float32"))
    return out


def compare(a, b):
    res, rows = same_names("obs_names", a["obs_names"], b["obs_names"])
    out = [res]
    res, vperm = same_names("var_names", a["var_names"], b["var_names"])
    out.append(res)
    out += columns("obs", a["obs"], b["obs"], rows)
    out += columns("var", a["var"], b["var"], vperm)
    out += matrices(a, b, rows)
    return out


if __name__ == "__main__":
    print(json.dumps(compare(load(sys.argv[1]), load(sys.argv[2])), indent=1))
