# Extended benchmark: datasets drawn by rule, plus a size ladder

`plan.yaml` reproduces the CrossCell paper on its own 55 cases. Those are small (most under
100k cells), and the 42 rds cases are 11 datasets in up to 4 variants. `plan-extended.yaml`
asks the same question on datasets nobody picked by hand.

## 1. Question

How do the converters compare in time, memory and fidelity on datasets drawn blindly, by a rule
fixed before the draw, across three size bins (~10k, ~100k, ~1M cells), and how do time and memory
grow with size when only size changes?

The benchmarker develops one of the tools (scx), so the dataset choice cannot be left to judgment.
CELLxGENE Discover is the frame: thousands of public datasets with one schema, a cell count and
mean genes per cell in its listing, and a direct h5ad URL for each.

## 2. Data facts

- [paper] CELLxGENE curation API `GET /curation/v1/datasets` lists every public dataset with
  `cell_count`, `mean_genes_per_cell`, `assay`, `organism`, `collection_id`, and an H5AD asset.
  Snapshot used: `extended/cxg_datasets.json.gz`, fetched 2026-11-06 with hapiq (`00c0a23`),
  sha256 of the uncompressed JSON `1333f4f3d7c97285ccb0778634b45e376fb2eef0f7c93c8f058ceb6c99db5221`,
  2,269 datasets.
- [paper] CELLxGENE offers H5AD only (no Seurat .rds download), so every rds input is built here.
- [verify] Schema-5 h5ad: `X` normalized, `raw/X` raw counts when the submitter gave them.
- [verify] nnz ≈ `cell_count × mean_genes_per_cell`; an R `dgCMatrix` holds at most 2^31 − 1 nonzeros.

## 3. Construction

**Bins** (by `cell_count`): S = [5,000, 20,000), M = [50,000, 200,000), L = [500,000, 2,000,000].

**Draw** (`scripts/select_datasets.py`, run once on the snapshot):
1. Candidates: every listed dataset with an H5AD asset, in a bin, not tombstoned.
2. Order: sort a bin's candidates by `dataset_id`, permute with `numpy.random.default_rng(20261106)`.
3. Accept greedily in that order until the bin has 5: skip a dataset whose `collection_id` is
   already accepted (any bin; one dataset per study) or whose primary assay (first `assay` label)
   already has 2 accepted in that bin.
4. Rescue, registered now: a dataset that fails to download (after one retry) or that anndata cannot
   read is a plumbing failure: log it, drop it, take the next in order. Nothing else replaces a draw:
   not size, not a tool failing on it, not its results.

**Inputs per dataset** (one recipe, no converter under test involved):
- h5ad: the CELLxGENE file as downloaded (h5ad → Seurat arm).
- rds (rds → h5ad arm): SeuratObject v5, assay `RNA`, `counts` = `raw/X` (else `X`), `data` = `X`,
  meta.data = obs (categoricals as factors), feature metadata = var. Built with rhdf5 + SeuratObject,
  `saveRDS` default compression. Not built when nnz ≥ 2^31 − 1: logged as "not representable as one
  dgCMatrix", and the rds arm has no case for that dataset.

**Ladder**: source = the first accepted L dataset with `cell_count` ≥ 1,000,000 (if none, the first
such candidate in L's order that passes the collection rule). Rungs 10k, 30k, 100k, 300k, 1M cells:
one permutation of the cells (`default_rng(20261106)`), rung n = its first n cells in file order,
so rungs are nested. All genes kept. h5ad and rds per rung by the recipe above.

**Measurement**: as `plan.yaml` (GNU time + denet, obkit phases, 3 runs, median; first failing run
stops a case), timeout 3,600 s, conversions run one at a time (`--cores 2`). Fidelity check as
`plan.yaml`. A timeout or out-of-memory is a result for that size, not a reason to drop it.

## 4. Pre-registration

**Unit**: the dataset (the ladder is one dataset and is reported apart, never pooled with the bins).

**Gate (plumbing)**: for every built rds, our reader's snapshot of the rds must match the source
h5ad on cell names, gene names, obs columns and the counts matrix (`modules/compare.py`). A failure
is a builder bug: fix, log, rebuild, before any conversion is timed.

**Negative control**: `scx diff` and our check must call a file identical to itself identical
(input vs input, every dataset). A failure means the checker is wrong, not the tools.

**What is reported, per bin and direction**: converted n, faithful n, median wall and peak RSS over
converted datasets; paired per-dataset time ratio of each tool to scx. Ladder: time and peak RSS
against nnz per tool, log-log slope with 95 % CI.

**Claims the application note may make, each mechanical** (`scripts/verdict.py`, run on results):
- C1 speed: scx is faster than each other tool on ≥ 80 % of the datasets both converted, in every
  bin and direction. Otherwise the note states the bins and tools where it is not.
- C2 memory: on the ladder, scx's log–log slope of peak RSS on nnz has its 95 % CI upper bound < 1
  (sub-linear). Otherwise the note says memory grows linearly (or worse).
- C3 coverage: scx's converted count is reported per bin as is; no claim of "converts everything".
- Fidelity is reported per tool and property kind, with no win criterion.

**Decision table**:
| outcome | action |
|---|---|
| gate fails | fix the builder, rebuild, rerun the gate; no timings before it passes |
| C1 fails in a bin | report that bin as a loss; no re-draw, no extra datasets |
| C2 fails | report linear memory growth |
| a tool fails at L size | that is the L result for that tool |
| fewer than 5 accepted in a bin | report the smaller n; do not widen the bin |

**Expected outcome** (written before the draw): scx fastest in S and M on ≥ 80 % of datasets,
both directions; at L the R tools time out or exceed 100 GB; the 1M rds rung is not representable
for most sources; scx's memory slope on the ladder below 1; scx fails some datasets whose
primary assay is spatial, and changes bool obs columns to integers (known).

## 5. Null and power

No hypothesis test on C1: it is a count over ≤ 5 datasets per bin, reported as such. The ladder
slope CI comes from the OLS fit on 5 rungs × 3 runs.

## 6. Work plan

1. Freeze this file and `scripts/select_datasets.py` (commit, then write the hash below).
   Check: the commit exists and the snapshot's sha256 matches §2.
2. Run the draw. Check: ≤ 5 per bin, one per collection, ≤ 2 per assay per bin; write
   `config/extended.json`.
3. Download (hapiq, md5 from the listing where given) and build inputs on hainan. Check: the gate.
4. Negative control. Check: identical called identical.
5. Generate `plan-extended.yaml`; run conversions one at a time, then checks and summary.
6. `scripts/verdict.py` prints C1–C3; record "Step 6 done" with the tables verbatim.

FROZEN at commit: (filled in by the next commit)
