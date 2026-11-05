#!/usr/bin/env Rscript
# The CrossCell paper's type-fidelity input (benchmark/notebooks/
# compute_type_fidelity.ipynb, verbatim apart from the output path): a
# 50-cell Seurat object with one metadata column per type. Staged as one more
# RDS test case, `type_test`, so every tool converts it like the others.
a <- commandArgs(TRUE)
out <- sub("^[^=]*=", "", a[which(a == "--output") + 1])

suppressPackageStartupMessages({
    library(Seurat)
    library(Matrix)
})

set.seed(42)
n_cells <- 50
n_genes <- 100

counts <- Matrix(rpois(n_genes * n_cells, lambda=2), nrow=n_genes, ncol=n_cells, sparse=TRUE)
rownames(counts) <- paste0("Gene", seq_len(n_genes))
colnames(counts) <- paste0("Cell", seq_len(n_cells))

obj <- CreateSeuratObject(counts=counts, project="type_test")

obj$int_col <- as.integer(sample(1:10, n_cells, replace=TRUE))
obj$dbl_col <- rnorm(n_cells, mean=5.0, sd=1.0)
obj$bool_col <- sample(c(TRUE, FALSE), n_cells, replace=TRUE)
obj$str_col <- sample(c("alpha", "beta", "gamma"), n_cells, replace=TRUE)
obj$cat_col <- factor(sample(c("TypeA", "TypeB", "TypeC"), n_cells, replace=TRUE))

saveRDS(obj, out)
