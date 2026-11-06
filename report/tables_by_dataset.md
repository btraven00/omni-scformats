# By dataset

As above, but one row per dataset: the 4 variants of an rds dataset (v4/v5 × raw/processed)
count once, with the median over its converted variants. Faithful = every variant faithful.
Converted = at least one variant converted. n = 11 rds and 13 h5ad datasets.

## Seurat .rds → h5ad

| tool | converted | faithful | wall | excl. start-up | range | peak RSS |
|---|--:|--:|--:|--:|--:|--:|
| scx | 9/11 | 0/11 | 3.1 s | 3.1 s | 0.7–7.2 | 246 MB |
| crosscell | 11/11 | 0/11 | 6.8 s | 6.8 s | 0.7–16.8 | 2.0 GB |
| anndataR | 11/11 | 7/11 | 15.9 s | 12.2 s | 7.0–27.0 | 2.6 GB |
| zellkonverter | 11/11 | 0/11 | 21.1 s | 12.5 s | 13.0–33.9 | 3.4 GB |
| convert2anndata | 10/11 | 0/11 | 25.6 s | 17.1 s | 16.4–36.8 | 4.2 GB |
| easySCF | 11/11 | 0/11 | 14.7 s | 11.0 s | 6.3–42.8 | 1.7 GB |

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
| scx | rds_to_h5ad |  | 1.8 s / 46 MB |  | 1.1 s / 256 MB | 0.2 s / 23 MB |
| scx | h5ad_to_rds |  | 0.2 s / 118 MB |  | 0.3 s / 194 MB | 0.4 s / 75 MB |
| crosscell | rds_to_h5ad |  | 5.0 s / 1.8 GB | 0.2 s / 1.4 GB |  | 1.0 s / 1.2 GB |
| crosscell | h5ad_to_rds |  | 0.2 s / 298 MB | 0.1 s / 666 MB |  | 15.0 s / 311 MB |
| anndataR | rds_to_h5ad | 3.7 s / 419 MB | 4.5 s / 1.3 GB | 4.9 s / 2.6 GB |  | 2.2 s / 2.6 GB |
| anndataR | h5ad_to_rds | 3.8 s / 418 MB | 4.3 s / 1.2 GB | 6.2 s / 2.8 GB |  | 17.8 s / 2.9 GB |
| zellkonverter | rds_to_h5ad | 8.6 s / 744 MB | 5.1 s / 1.6 GB | 2.4 s / 2.2 GB |  | 5.0 s / 3.4 GB |
| zellkonverter | h5ad_to_rds | 8.4 s / 932 MB | 3.7 s / 1.1 GB | 1.9 s / 1.2 GB |  | 6.4 s / 1.2 GB |
| convert2anndata | rds_to_h5ad | 8.4 s / 724 MB | 3.6 s / 1.2 GB | 11.3 s / 4.2 GB |  | 1.4 s / 3.8 GB |
| easySCF | rds_to_h5ad | 3.7 s / 413 MB | 4.5 s / 1.3 GB | 0.7 s / 1.4 GB |  | 8.1 s / 1.7 GB |
