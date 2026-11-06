#!/usr/bin/env python3
"""The draw of extended/PLAN.md §3, run once on the frozen CELLxGENE listing.

    python3 -I scripts/select_datasets.py extended/cxg_datasets.json.gz > config/extended.json
"""
import gzip, json, sys

import numpy as np

SEED, K, PER_ASSAY = 20261106, 5, 2
BINS = {"S": (5_000, 20_000), "M": (50_000, 200_000), "L": (500_000, 2_000_000)}
LADDER = [10_000, 30_000, 100_000, 300_000, 1_000_000]
NNZ_MAX = 2**31 - 1

listing = json.load(gzip.open(sys.argv[1]))


def h5ad(d):
    return next((a for a in d.get("assets", []) if a.get("filetype") == "H5AD"), None)


def assay(d):
    return d["assay"][0]["label"] if d.get("assay") else ""


def ordered(lo, hi):
    c = sorted((d for d in listing if lo <= d["cell_count"] < hi and h5ad(d) and not d.get("tombstone")),
               key=lambda d: d["dataset_id"])
    return [c[i] for i in np.random.default_rng(SEED).permutation(len(c))]


def row(d, bin_):
    nnz = round(d["cell_count"] * d["mean_genes_per_cell"])
    a = h5ad(d)
    return {"bin": bin_, "dataset_id": d["dataset_id"], "collection_id": d["collection_id"],
            "title": d["title"], "assay": assay(d), "organism": [o["label"] for o in d["organism"]],
            "n_cells": d["cell_count"], "nnz_est": nnz, "rds": nnz < NNZ_MAX,
            "url": a["url"], "filesize": a.get("filesize")}


seen, picks, order = set(), [], {}
for b, (lo, hi) in BINS.items():
    order[b] = ordered(lo, hi)
    per_assay = {}
    for d in order[b]:
        if sum(p["bin"] == b for p in picks) == K:
            break
        if d["collection_id"] in seen or per_assay.get(assay(d), 0) == PER_ASSAY:
            continue
        seen.add(d["collection_id"])
        per_assay[assay(d)] = per_assay.get(assay(d), 0) + 1
        picks.append(row(d, b))

src = next((p for p in picks if p["bin"] == "L" and p["n_cells"] >= 1_000_000), None)
if src is None:
    src = next(row(d, "L") for d in order["L"]
               if d["cell_count"] >= 1_000_000 and d["collection_id"] not in seen)
ladder = [{"rung": n, "n_cells": n, "nnz_est": round(n * src["nnz_est"] / src["n_cells"]),
           "rds": n * src["nnz_est"] / src["n_cells"] < NNZ_MAX} for n in LADDER]
json.dump({"seed": SEED, "datasets": picks, "ladder": {"source": src["dataset_id"], "rungs": ladder},
           "reserve_order": {b: [d["dataset_id"] for d in o[:50]] for b, o in order.items()}},
          sys.stdout, indent=1)
print()
