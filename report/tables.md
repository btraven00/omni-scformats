# CrossCell benchmark, reproduced

Per tool and direction: cases converted / cases, cases passing every fidelity property,
then wall time and peak RSS (GNU time ru_maxrss) as the median over converted cases of
each case's median of 3 runs; range = min–max over those cases.

## Seurat .rds → h5ad

| tool | converted | faithful | wall | range | peak RSS |
|---|--:|--:|--:|--:|--:|
| scx | 36/42 | 13/42 | 4.0 s | 0.7–11.4 | 246 MB |
| crosscell | 42/42 | 0/42 | 5.5 s | 1.1–23.0 | 1.6 GB |
| anndataR | 41/42 | 35/42 | 20.3 s | 8.2–49.2 | 2.0 GB |
| zellkonverter | 41/42 | 9/42 | 27.6 s | 15.4–61.0 | 2.6 GB |
| convert2anndata | 19/42 | 5/42 | 24.9 s | 16.2–48.8 | 2.9 GB |
| easySCF | 41/42 | 0/42 | 23.8 s | 9.5–59.7 | 2.1 GB |

## h5ad → Seurat-readable

| tool | converted | faithful | wall | range | peak RSS |
|---|--:|--:|--:|--:|--:|
| scx | 8/13 | 3/13 | 1.8 s | 1.0–8.8 | 147 MB |
| crosscell | 13/13 | 2/13 | 17.0 s | 1.2–362.5 | 342 MB |
| anndataR | 3/13 | 1/13 | 35.4 s | 20.5–181.0 | 2.9 GB |
| zellkonverter | 13/13 | 2/13 | 32.8 s | 12.1–300.7 | 1.2 GB |

## By phase

Median over converted cases of each phase's seconds, and of its denet-sampled peak.
R tools: obkit anchors in each step. scx: its own log timestamps. crosscell: its
printed step durations laid end to end (long steps are rounded to whole seconds).

| tool | direction | start-up (R + packages / binary) | read | convert | matrix stream (scx: read, convert, write in one pass) | write (saveRDS / h5ad; scx: finalize, flush, hashes) |
|---|---|--:|--:|--:|--:|--:|
| scx | rds_to_h5ad |  | 1.2 s / 51 MB |  | 1.7 s / 267 MB | 0.3 s / 43 MB |
| scx | h5ad_to_rds |  | 0.3 s / 107 MB |  | 0.6 s / 321 MB | 0.6 s / 76 MB |
| crosscell | rds_to_h5ad |  | 3.0 s / 1.3 GB | 0.3 s / 1.7 GB |  | 1.0 s / 1.2 GB |
| crosscell | h5ad_to_rds |  | 0.2 s / 352 MB | 0.1 s / 850 MB |  | 16.0 s / 311 MB |
| anndataR | rds_to_h5ad | 5.1 s / 419 MB | 4.2 s / 881 MB | 6.1 s / 2.0 GB |  | 3.1 s / 2.0 GB |
| anndataR | h5ad_to_rds | 4.0 s / 426 MB | 5.0 s / 1.2 GB | 7.6 s / 2.8 GB |  | 18.1 s / 2.9 GB |
| zellkonverter | rds_to_h5ad | 9.9 s / 746 MB | 4.6 s / 1.2 GB | 3.0 s / 1.6 GB |  | 6.0 s / 2.5 GB |
| zellkonverter | h5ad_to_rds | 9.9 s / 932 MB | 5.4 s / 1.1 GB | 2.6 s / 1.2 GB |  | 6.4 s / 1.2 GB |
| convert2anndata | rds_to_h5ad | 8.8 s / 723 MB | 3.2 s / 1.0 GB | 11.3 s / 2.9 GB |  | 1.4 s / 2.9 GB |
| easySCF | rds_to_h5ad | 4.2 s / 414 MB | 3.8 s / 875 MB | 0.9 s / 1.3 GB |  | 11.6 s / 2.1 GB |

## The paper's own checks, as published

Pass rates as in its README: converted = every run exited 0 (no output check).

| Tool | rds_to_h5ad | h5ad_to_rds |
|---|--:|--:|
| scx | 36/42 | 8/13 |
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
