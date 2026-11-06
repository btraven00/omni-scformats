# CrossCell benchmark, reproduced

Per tool and direction: cases converted / cases, cases passing every fidelity property,
then wall time and peak RSS (GNU time ru_maxrss) as the median over converted cases of
each case's median of 3 runs; range = min–max over those cases. "excl. start-up" subtracts
R + package loading (4–10 s per run), which a user converting inside an R session
doesn't pay. The 42 rds cases are 11 datasets in up to 4 variants (v4/v5 × raw/processed),
so per-dataset n is 11, not 42.

## Seurat .rds → h5ad

| tool | converted | faithful | wall | excl. start-up | range | peak RSS |
|---|--:|--:|--:|--:|--:|--:|
| scx | 36/42 | 14/42 | 2.5 s | 2.5 s | 0.6–7.8 | 246 MB |
| crosscell | 42/42 | 0/42 | 5.2 s | 5.2 s | 0.6–22.0 | 1.6 GB |
| anndataR | 41/42 | 34/42 | 12.8 s | 9.1 s | 6.0–29.3 | 2.0 GB |
| zellkonverter | 41/42 | 8/42 | 19.3 s | 10.8 s | 12.2–36.2 | 2.6 GB |
| convert2anndata | 19/42 | 4/42 | 23.1 s | 14.8 s | 14.7–44.5 | 2.9 GB |
| easySCF | 41/42 | 0/42 | 18.2 s | 14.6 s | 6.3–51.1 | 2.1 GB |

## h5ad → Seurat-readable

| tool | converted | faithful | wall | excl. start-up | range | peak RSS |
|---|--:|--:|--:|--:|--:|--:|
| scx | 8/13 | 3/13 | 1.1 s | 1.1 s | 0.6–6.6 | 152 MB |
| scConvert | 12/13 | 6/13 | 0.7 s | 0.7 s | 0.6–1.1 | 28 MB |
| crosscell | 13/13 | 2/13 | 16.7 s | 16.7 s | 0.6–346.5 | 342 MB |
| anndataR | 3/13 | 1/13 | 32.7 s | 28.9 s | 18.5–175.3 | 2.9 GB |
| zellkonverter | 13/13 | 2/13 | 21.1 s | 12.8 s | 10.7–289.5 | 1.2 GB |

## By phase

Median over converted cases of each phase's seconds, and of its denet-sampled peak.
R tools: obkit anchors in each step. scx: its own log timestamps. crosscell: its
printed step durations laid end to end (long steps are rounded to whole seconds).

| tool | direction | start-up (R + packages / binary) | read | convert | matrix stream (scx: read, convert, write in one pass) | write (saveRDS / h5ad; scx: finalize, flush, hashes) |
|---|---|--:|--:|--:|--:|--:|
| scx | rds_to_h5ad |  | 1.0 s / 47 MB |  | 1.1 s / 256 MB | 0.2 s / 43 MB |
| scx | h5ad_to_rds |  | 0.2 s / 118 MB |  | 0.3 s / 194 MB | 0.4 s / 75 MB |
| crosscell | rds_to_h5ad |  | 3.0 s / 1.3 GB | 0.3 s / 1.4 GB |  | 1.0 s / 1.1 GB |
| crosscell | h5ad_to_rds |  | 0.2 s / 298 MB | 0.1 s / 666 MB |  | 15.0 s / 311 MB |
| anndataR | rds_to_h5ad | 3.7 s / 419 MB | 2.6 s / 880 MB | 4.5 s / 2.0 GB |  | 1.9 s / 2.0 GB |
| anndataR | h5ad_to_rds | 3.8 s / 418 MB | 4.3 s / 1.2 GB | 6.2 s / 2.8 GB |  | 17.8 s / 2.9 GB |
| zellkonverter | rds_to_h5ad | 8.6 s / 745 MB | 3.7 s / 1.2 GB | 1.8 s / 1.6 GB |  | 4.8 s / 2.5 GB |
| zellkonverter | h5ad_to_rds | 8.4 s / 932 MB | 3.7 s / 1.1 GB | 1.9 s / 1.2 GB |  | 6.4 s / 1.2 GB |
| convert2anndata | rds_to_h5ad | 8.4 s / 724 MB | 3.0 s / 1.0 GB | 10.1 s / 2.9 GB |  | 1.1 s / 2.9 GB |
| easySCF | rds_to_h5ad | 3.7 s / 413 MB | 2.7 s / 875 MB | 0.6 s / 1.2 GB |  | 10.7 s / 2.1 GB |

## The paper's own checks, as published

Pass rates as in its README: converted = every run exited 0 (no output check).

| Tool | rds_to_h5ad | h5ad_to_rds |
|---|--:|--:|
| scx | 36/42 | 8/13 |
| scConvert | N/A | 12/13 |
| crosscell | 42/42 | 13/13 |
| anndataR | 41/42 | 3/13 |
| zellkonverter | 41/42 | 13/13 |
| convert2anndata | 19/42 | N/A |
| easySCF | 41/42 | N/A |

Matrix fidelity (compute_fidelity.ipynb): each tool's h5ad X against crosscell's h5ad X,
matched by position, means over cases; crosscell's row is its own h5ad → rds → h5ad
round trip. Cells and genes are not matched by name and the reference is a tool, not
the input, so this cannot see renamed cells or a wrong reference.

| tool | cases | Pearson r | Spearman ρ | MSE | exact match | not computed |
|---|--:|--:|--:|--:|--:|---|
| scx | 36/42 | 1.000000 | 1.000000 | 1.17e-12 | 100.0% | failed |
| crosscell | 42/42 | 1.000000 | 1.000000 | 0.00e+00 | 100.0% |  |
| anndataR | 41/42 | 1.000000 | 1.000000 | 0.00e+00 | 100.0% | failed |
| zellkonverter | 41/42 | 1.000000 | 1.000000 | 0.00e+00 | 100.0% | failed |
| convert2anndata | 19/42 | 1.000000 | 1.000000 | 1.11e-12 | 100.0% | failed |
| easySCF | 36/42 | 1.000000 | 1.000000 | 0.00e+00 | 100.0% | failed, read_error |

Type retention (compute_type_fidelity.ipynb): the 50-cell object with one column per
type, converted to h5ad; the class read back per column (values not checked).

| tool | retention | int_col | dbl_col | bool_col | str_col | cat_col |
|---|--:|---|---|---|---|---|
| scx | 100% | int | float | bool | str | category |
| crosscell | 100% | int | float | bool | str | category |
| anndataR | 100% | int | float | bool | str | category |
| zellkonverter | 80% | int | float | bool | str→category ✖ | category |
| convert2anndata | 80% | int | float | bool | str→category ✖ | category |
| easySCF | 0% | READ_ERROR | READ_ERROR | READ_ERROR | READ_ERROR | READ_ERROR |
