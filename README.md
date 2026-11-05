# omni-scformats: the CrossCell benchmark, reproduced

The CrossCell paper benchmark ([crosscell](https://github.com/biodancerwangzhi/crosscell)
`benchmark/`, runner `scripts/run_benchmark.py` recovered from git `a910b12`), ported to an
omnibenchmark plan at `api_version: "0.8.0"`. Same 55 test cases (42 Seurat RDS v4/v5,
13 H5AD) from the same Zenodo record, the same 5 tools with their exact commands, the same
measurement (3 runs, median wall, GNU-time ru_maxrss, 600 s timeout) on current tool versions, the same
n_cells x n_genes check. Plus scx arms.

```
RDS (42) --> TO_H5AD   (crosscell zellkonverter anndataR convert2anndata easySCF scx) --> CHECK_H5AD   --+
H5AD (13) -> TO_SEURAT (crosscell zellkonverter anndataR scx)                            --> CHECK_SEURAT -+-> SUMMARY (gather)
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
micromamba create -y -p data/scx-env -c conda-forge -c https://prefix.dev/edge scx=0.5.0.post147
export SCX_BIN=$PWD/data/scx-env/bin/scx         # pinned scx (7247a06); its HDF5 pin clashes with BPCells', so not in envs/
export DENET=$(command -v denet)
scripts/run.sh --cores 16                                      # full run: omni / hainan, not the laptop
scripts/run.sh --cores 8 --filter slice.obfilter.yaml          # laptop slice: 3 v5 RDS + 2 small H5AD
```

The slice needs only `scripts/fetch.sh h5ad_small seurat_v5`.
Results: `out/SUMMARY.tables.*/results.tsv` (one row per case x tool, phases as JSON) and
`pass_rates.md` (the README table of the paper, plus a validated column).

## What differs from the paper, and why

- **Phases.** Every R step is wrapped in obkit anchors (`load/read/convert/write`), joined to
  a denet trace (50 ms, children included). Total wall and ru_maxrss are the paper's numbers;
  phases let time-to-object be read without `saveRDS`, which is most of an h5ad -> rds run.
- **scx arm**, the binary alone, both directions. h5ad -> Seurat writes scx's default
  `.h5seurat` (BPCells-packed counts + meta.data), a file R reads, not an `.rds`: `.rds` is not
  an interchange format (scx `docs/rds-write.md`), and opening the file is the consumer's job,
  not the converter's. The original tools keep writing `.rds`. TO_SEURAT's declared
  `<tool>.seurat` is a symlink to what the tool wrote.
- **Validation is value-level**, not n_cells x n_genes. `modules/snapshot.{py,R}` read input
  and output with the format's own reader (anndata; `readRDS` + SeuratObject; rhdf5 + BPCells
  for `.h5seurat`) into one neutral layout; `modules/compare.py` checks names, every obs/var
  column (type class, values, NA mask) and every matrix (exact when integer-valued, rtol 1e-6
  otherwise). Fail: names changed or reordered, a column dropped, a type class changed, a
  matrix not found in the output. Notes: level order, extra columns, float32 storage; a passing
  matrix records where it landed (`RNA/counts = X`). `scx diff --json` runs alongside as a
  cross-check; disagreement is flagged, not resolved. Test: `modules/test_compare.py`.
  The original's test-matrix dims were also wrong for 4 of 13 h5ad cases.
- **Versions.** Current releases, not the paper's: conda envs, unpinned, instead of its Docker
  image (e.g. anndataR 1.2.x vs the paper's 1.0.0, Seurat 5.5 vs 5.4). The paper's pass rates
  are a baseline to compare against, not a target. convert2anndata and easySCFr are GitHub-only
  and install in `envs/tools.post-deploy.sh`, which ob does not carry into `out/.envs/`:
  `scripts/run.sh` puts it there. Checks and the scx arm run in `envs/readers.yml`
  (BPCells does not solve alongside the tools env).
- **Ids.** api 0.8 ids are `[A-Za-z0-9_]`: `celegans.embryo` -> `celegans_embryo`,
  `thp1.eccite` -> `thp1_eccite`.

Not ported: the paper's fidelity, robustness, round-trip and scalability notebooks
(`benchmark/notebooks/compute_*`), only the conversion benchmark.
