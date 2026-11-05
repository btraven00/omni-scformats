"""snapshot.py + compare.py on synthetic corruptions (the h5ad read itself is
exercised by the slice run). Run with the readers env:
    <env>/bin/python modules/test_compare.py
"""
import sys, tempfile
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).parent))
import compare, snapshot


def base():
    rng = np.random.default_rng(0)
    X = sp.random(20, 8, density=0.3, format="csr", random_state=0)
    X.data = np.round(X.data * 10) + 1  # counts
    obs = pd.DataFrame({
        "cond": pd.Categorical(rng.choice(["ctrl", "kd"], 20)),
        "n": rng.integers(0, 100, 20),
        "score": rng.random(20),
        "flag": rng.random(20) > 0.5,
        "donor": pd.array(rng.choice(["a", "b", None], 20), dtype=object),
    }, index=[f"bc{i}" for i in range(20)])
    return ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=[f"g{j}" for j in range(8)]))


def verdicts(a, b):
    with tempfile.TemporaryDirectory() as t:
        for x, d in ((a, "a"), (b, "b")):
            snapshot.write(x, f"{t}/{d}")  # in memory: write_h5ad would re-categorise strings
        return {p["property"]: p["status"] for p in compare.compare(compare.load(f"{t}/a"), compare.load(f"{t}/b"))}


a = base()
v = verdicts(a, base())
assert set(v.values()) == {"pass"}, v

b = base(); b.obs_names = [f"Cell_{i + 1}" for i in range(20)]
v = verdicts(a, b)
assert v["obs_names"] == "fail" and v["obs/cond"] == "pass" and v["matrix/X"] == "pass", v  # aligned by position

b = base(); b.obs.loc[b.obs_names[:5], "cond"] = "kd"; b.obs.loc[b.obs_names[5:10], "cond"] = "ctrl"
assert verdicts(a, b)["obs/cond"] == "fail"

b = base(); del b.obs["score"]
assert verdicts(a, b)["obs/score"] == "fail"

b = base(); b.obs["n"] = b.obs["n"].astype(float)
assert verdicts(a, b)["obs/n"] == "fail"  # int -> float: type class changed

b = base(); b.obs["cond"] = b.obs["cond"].astype(str)
assert verdicts(a, b)["obs/cond"] == "fail"  # categorical -> string

b = base(); b.obs["cond"] = b.obs["cond"].cat.reorder_categories(["kd", "ctrl"])
assert verdicts(a, b)["obs/cond"] == "note"  # level order only

b = base(); b.X = b.X.astype(np.float32)
v = verdicts(a, b)
assert v["matrix/X"] == "pass" and v["matrix/X.dtype"] == "note", v

b = base(); b.X = b.X.copy(); b.X.data[0] += 1
assert verdicts(a, b)["matrix/X"] == "fail"

b = base()[::-1].copy()
v = verdicts(a, b)
assert v["obs_names"] == "fail" and v["obs/n"] == "pass" and v["matrix/X"] == "pass", v  # reordered, realigned

b = base(); b.obs["extra"] = 1
assert verdicts(a, b)["obs/+"] == "note"
print("ok")
