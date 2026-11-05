# omni-scformats

The [CrossCell](https://github.com/biodancerwangzhi/crosscell) conversion benchmark, rerun as an
omnibenchmark plan (api 0.8). Same 55 test cases from the same Zenodo record, same five tools
with their original commands, plus scx.

```
RDS (42) --> TO_H5AD   (crosscell zellkonverter anndataR convert2anndata easySCF scx) --> CHECK_H5AD   --+
H5AD (13) -> TO_SEURAT (crosscell zellkonverter anndataR scx)                            --> CHECK_SEURAT -+-> SUMMARY
```

## Run

You need omnibenchmark at HEAD (the 0.7.0 release can't read api 0.8 plans) and three binaries
on hand:

```bash
scripts/fetch.sh                        # 6.4 GB from Zenodo, via hapiq
export CROSSCELL_DATA=$PWD/data
export CROSSCELL_BIN=$PWD/data/crosscell   # the paper's prebuilt 0.1.0, from crosscell/benchmark/
micromamba create -y -p data/scx-env -c conda-forge -c https://prefix.dev/edge scx=0.5.0.post147
export SCX_BIN=$PWD/data/scx-env/bin/scx
export DENET=$(command -v denet)

scripts/run.sh --cores 32                                 # full run, on omni or hainan
scripts/run.sh --cores 8 --filter slice.obfilter.yaml     # small slice for a laptop
```

The slice only needs `scripts/fetch.sh h5ad_small seurat_v5`. Results end up in
`out/SUMMARY.tables.*/`: `results.tsv` (one row per case and tool) and `pass_rates.md`.

`benchmark.yaml` is generated from the paper's test matrix; edit `scripts/gen_plan.py`, not the YAML.

## How it differs from the paper

- **Current tool versions**, from conda, not the paper's Docker image. Their numbers are a
  baseline, not a target (they used anndataR 1.0, for example). One pin: anndata 0.12.9 in the
  tools env, because convert2anndata calls `AnnData(dtype=...)`, which 0.13 removed. The
  checks read with current anndata.
- **Phases.** R steps are split into load / read / convert / write, so you can see how much of
  an h5ad → rds run is just `saveRDS`. Memory is GNU time's peak plus a denet trace.
- **scx writes `.h5seurat`, not `.rds`.** `.rds` isn't an interchange format; scx writes a file
  R can read (BPCells-backed counts plus metadata) and stops there.
- **Fidelity, not just shape.** The paper only checked n_cells × n_genes, against a test matrix
  that's wrong for 4 of its 13 h5ad cases. Here input and output are read with each format's
  own reader (anndata, readRDS, BPCells) and compared on names, every obs/var column (type,
  values, missing values) and every matrix. `scx diff` runs alongside as a second opinion.
  (Their `results/validation.json`, where that shape check should land, is an empty list; and
  the 3 CELLxGENE rows of the test matrix list 33,694 genes each, the 10x GRCh38 count, not
  the files' 17,222 / 32,330 / 58,232.)
- **IDs** lose their dots (`celegans.embryo` → `celegans_embryo`); api 0.8 requires it.

Not ported: the paper's round-trip, robustness and scalability notebooks.

## Gotchas

- Use `scripts/run.sh`, not `ob run`: ob doesn't copy `envs/tools.post-deploy.sh` (which
  installs convert2anndata and easySCF from GitHub), and running ob through `uv run` leaks its
  venv's Python into the tools.
- zellkonverter's basilisk compiles Python 3.14 from source unless you give it one:
  `export BASILISK_CUSTOM_PYTHON_3_14_0=<conda python=3.14.0>/bin/python`.
- Modules run from the committed tree, so commit before running.
