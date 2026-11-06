#!/usr/bin/env python3
"""Build the extended plan's inputs (extended/PLAN.md §3, step 3), in the readers env.

    build_extended.py list              # one input name per line: ext_S1 ... ext_L5, ladder_10000 ...
    build_extended.py one <name>        # download if needed, stage h5ad, build rds, gate
    build_extended.py manifest          # sha256 + gate verdicts -> extended_inputs.json

Files land in $CROSSCELL_DATA/data/generated/ next to the paper's, so modules/data.py
stages both; raw downloads in $CROSSCELL_DATA/cxg/<dataset_id>.h5ad (hapiq url).
"""
import hashlib, json, os, subprocess, sys, tempfile
from pathlib import Path

root = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(root / "modules"))
import compare

cfg = json.loads((root / "config/extended.json").read_text())
data = Path(os.environ.get("CROSSCELL_DATA", root / "data"))
gen, cxg = data / "data/generated", data / "cxg"
HAPIQ = os.environ.get("HAPIQ", "hapiq")


def names():
    out, n = {}, {}
    for d in cfg["datasets"]:
        n[d["bin"]] = n.get(d["bin"], 0) + 1
        out[f"ext_{d['bin']}{n[d['bin']]}"] = d
    for r in cfg["ladder"]["rungs"]:
        out[f"ladder_{r['rung']}"] = r
    return out


def fetch(dataset_id):
    d = next(x for x in cfg["datasets"] if x["dataset_id"] == dataset_id)
    f = cxg / Path(d["url"]).name
    for _ in range(2):  # PLAN §3.4: one retry
        if f.is_file() and f.stat().st_size == d["filesize"]:
            return f
        subprocess.run([HAPIQ, "download", "url", d["url"], "--out", str(cxg), "-y"], check=False)
    sys.exit(f"download failed: {d['url']} (plumbing, PLAN §3.4: take the next in order)")


def snap(f, d):
    cmd = ([sys.executable, str(root / "modules/snapshot.py")] if f.suffix == ".h5ad"
           else ["Rscript", "--vanilla", str(root / "modules/snapshot.R")]) + [str(f), str(d)]
    subprocess.run(cmd, check=True, capture_output=True, text=True)


def gate(h5ad, rds):
    """Reading of PLAN §4's gate: compare(h5ad as reference, rds) so the columns
    CreateSeuratObject adds (orig.ident, nCount_RNA, nFeature_RNA) are extra, not failures."""
    with tempfile.TemporaryDirectory(dir=gen) as t:
        snap(h5ad, Path(t) / "a"); snap(rds, Path(t) / "b")
        props = compare.compare(compare.load(Path(t) / "a"), compare.load(Path(t) / "b"))
    fails = [f"{p['property']}: {p.get('detail', '')}" for p in props if p["status"] == "fail"]
    (gen / f"{rds.stem}.gate.json").write_text(json.dumps({"pass": not fails, "fails": fails, "properties": props}, indent=1))
    print(rds.stem, "gate", "PASS" if not fails else f"FAIL {fails}", flush=True)


def one(name):
    item = names()[name]
    h5ad, rds = gen / f"{name}.h5ad", gen / f"{name}.rds"
    if name.startswith("ext_"):
        src = fetch(item["dataset_id"])
        if not h5ad.exists():
            os.link(src, h5ad)
    elif not h5ad.exists():  # ladder rung: nested subsample of the source
        import anndata as ad, numpy as np
        a = ad.read_h5ad(fetch(cfg["ladder"]["source"]))
        perm = np.random.default_rng(cfg["seed"]).permutation(a.n_obs)
        for r in sorted(cfg["ladder"]["rungs"], key=lambda r: r["rung"]):
            f = gen / f"ladder_{r['rung']}.h5ad"
            if not f.exists():
                a[np.sort(perm[: r["rung"]])].copy().write_h5ad(f)
                print("wrote", f, flush=True)
    if item.get("rds", True) and not rds.exists():
        subprocess.run(["Rscript", "--vanilla", str(root / "scripts/build_rds.R"), str(h5ad), str(rds) + ".tmp"],
                       check=True)
        os.replace(str(rds) + ".tmp", rds)
    if rds.exists() and not (gen / f"{name}.gate.json").exists():
        gate(h5ad, rds)


def sha(f):
    h = hashlib.sha256()
    with open(f, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 24), b""):
            h.update(chunk)
    return h.hexdigest()


cmd = sys.argv[1]
if cmd == "list":
    # ladder rungs share one subsample pass: list one of them first
    print("\n".join(names()))
elif cmd == "one":
    one(sys.argv[2])
elif cmd == "manifest":
    out = {}
    for n, item in names().items():
        for ext in ("h5ad", "rds"):
            f = gen / f"{n}.{ext}"
            if f.exists():
                out[f.name] = {"sha256": sha(f)}
        g = gen / f"{n}.gate.json"
        if g.exists():
            out[f"{n}.rds"]["gate"] = json.loads(g.read_text())["pass"]
    (gen / "extended_inputs.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out, indent=1))
