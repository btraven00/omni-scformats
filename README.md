# omni-scformats: the CrossCell benchmark, reproduced

The CrossCell paper benchmark ([crosscell](https://github.com/biodancerwangzhi/crosscell)
`benchmark/`, runner `scripts/run_benchmark.py` recovered from git `a910b12`), ported to an
omnibenchmark plan at `api_version: "0.8.0"`. Same 55 test cases (42 Seurat RDS v4/v5,
13 H5AD) from the same Zenodo record, the same 5 tools with their exact commands, the same
measurement (3 runs, median wall, GNU-time ru_maxrss, 600 s timeout), the same
n_cells x n_genes check. Plus scx arms.

```
RDS (42) --> TO_H5AD   (crosscell zellkonverter anndataR convert2anndata easySCF scx) --> CHECK_H5AD   --+
H5AD (13) -> TO_SEURAT (crosscell zellkonverter anndataR scx_picklerick scx_bpcells)    --> CHECK_SEURAT -+-> SUMMARY (gather)
```

`benchmark.yaml` is generated: `python3 scripts/gen_plan.py > benchmark.yaml` from
`config/benchmark_testcases.json` (crosscell@af3ea57, verbatim).

## Run

Needs omnibenchmark with api 0.8 (`~/lab/omnibenchmark` HEAD; the released `ob` 0.7.0 rejects
the plan) and three side-loaded binaries.

```bash
scripts/fetch.sh                                 # hapiq: Zenodo 10.5281/zenodo.18610876, 6.4 GB
export CROSSCELL_DATA=$PWD/data                  # where fetch.sh unpacked; modules run from .modules/, so set it
export CROSSCELL_BIN=$PWD/data/crosscell         # cp ../crosscell/benchmark/crosscell; chmod +x  (the paper's 0.1.0 binary)
export SCX_BIN=~/phd/scx/target/release/scx      # cargo build --release -p scx-cli
export DENET=$(command -v denet)
ob run benchmark.yaml --dirty -k --cores 16      # full run: omni / hainan, not the laptop
ob run benchmark.yaml --dirty -k --cores 8 --filter slice.obfilter.yaml   # laptop slice: 3 v5 RDS + 2 small H5AD
```

The slice needs only `scripts/fetch.sh h5ad_small seurat_v5`.
Results: `out/SUMMARY.tables.*/results.tsv` (one row per case x tool, phases as JSON) and
`pass_rates.md` (the README table of the paper, plus a validated column).

## What differs from the paper, and why

- **Phases.** Every R step is wrapped in obkit anchors (`load/read/convert/write`), joined to
  a denet trace (50 ms, children included). Total wall and ru_maxrss are the paper's numbers;
  phases let time-to-object be read without `saveRDS`, which is most of an h5ad -> rds run.
- **scx arms.** scx does not write `.rds` (scx `docs/rds-write.md`). Its h5ad -> Seurat arms
  build the object in R and `saveRDS` there: `scx_picklerick` (`read_h5ad(as = "Seurat")`),
  `scx_bpcells` (`scx convert -> .bpcells`, `BPCells::open_matrix_dir`, `CreateSeuratObject`;
  X + dimnames only, no meta.data). `scx` itself does rds -> h5ad.
- **Validation** uses R's `readRDS` only; the original tried `crosscell inspect` first.
- **Envs.** conda instead of the paper's Docker image. Seurat has no R 4.5 conda build, so
  Seurat 5.4.0 and the GitHub-only tools (convert2anndata `550bfb0`, easySCFr `0a5a667`,
  main as of the paper's run date) install in `envs/tools.post-deploy.sh` from a dated CRAN
  snapshot. scx arms get their own env: picklerick needs a newer Bioconductor generation than
  zellkonverter 1.20.1.
- **Ids.** api 0.8 ids are `[A-Za-z0-9_]`: `celegans.embryo` -> `celegans_embryo`,
  `thp1.eccite` -> `thp1_eccite`.

Not ported: the paper's fidelity, robustness, round-trip and scalability notebooks
(`benchmark/notebooks/compute_*`), only the conversion benchmark.
