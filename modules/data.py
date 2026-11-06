#!/usr/bin/env python3
"""Stage one test case: hardlink $CROSSCELL_DATA/data/generated/<file> (laid
out by scripts/fetch.sh) and check its sha256 against the Zenodo manifest.
"""
import hashlib, json, os, shutil, sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from ob import parse

one, _, out = parse()
root = Path(os.environ.get("CROSSCELL_DATA", Path(__file__).resolve().parent.parent / "data"))
src = root / "data/generated" / one["file"]
if not src.is_file():
    sys.exit(f"dataset not found: {src}\nrun scripts/fetch.sh, or set CROSSCELL_DATA")

manifest = json.loads((root / "zenodo/manifest.json").read_text())
want = {d["file"]: d["sha256"] for arch in manifest["datasets"].values() for d in arch}
ext = root / "data/generated/extended_inputs.json"  # plan-extended.yaml: scripts/build_extended.py
if ext.is_file():
    want.update({k: v["sha256"] for k, v in json.loads(ext.read_text()).items()})
h = hashlib.sha256()
with open(src, "rb") as f:
    for chunk in iter(lambda: f.read(1 << 20), b""):
        h.update(chunk)
if h.hexdigest() != want.get(one["file"]):
    sys.exit(f"sha256 mismatch for {src}: {h.hexdigest()} != {want.get(one['file'])}")

dst = Path(next(iter(out.values())))
dst.parent.mkdir(parents=True, exist_ok=True)
dst.unlink(missing_ok=True)
try:
    os.link(src, dst)  # same fs: no copy charged to the run
except OSError:
    shutil.copy2(src, dst)
