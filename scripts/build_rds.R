#!/usr/bin/env Rscript
# The extended plan's rds recipe (extended/PLAN.md §3): one CELLxGENE h5ad ->
# SeuratObject v5 .rds, read with rhdf5 only (no converter under test).
#   Rscript scripts/build_rds.R in.h5ad out.rds
# counts = raw/X (else X), data = X, meta.data = obs, feature metadata = var.
suppressPackageStartupMessages({ library(rhdf5); library(Matrix); library(SeuratObject) })
a <- commandArgs(TRUE); f <- a[1]; out <- a[2]
l <- h5ls(f)
path <- paste0(sub("/$", "", l$group), "/", l$name)
has <- function(p) p %in% path
attrs <- function(p) h5readAttributes(f, p)

# anndata dataframe group: column-order attr, _index, categoricals as categories + codes
read_df <- function(g) {
  at <- attrs(g)
  idx <- as.character(h5read(f, paste0(g, "/", at[["_index"]])))
  cols <- as.character(at[["column-order"]])
  df <- lapply(cols, function(c) {
    q <- paste0(g, "/", c)
    if (l$otype[path == q] == "H5I_GROUP") {
      lv <- as.character(h5read(f, paste0(q, "/categories")))
      codes <- as.integer(h5read(f, paste0(q, "/codes")))
      factor(ifelse(codes < 0, NA, lv[codes + 1L]), levels = lv)
    } else {
      x <- h5read(f, q)
      if (l$dclass[path == q] == "ENUM" && all(as.character(x) %in% c("FALSE", "TRUE", NA)))
        as.logical(as.character(x)) else as.vector(x)
    }
  })
  names(df) <- cols
  df <- as.data.frame(df, check.names = FALSE, stringsAsFactors = FALSE)
  rownames(df) <- idx
  df
}

# genes x cells dgCMatrix from an anndata cells x genes matrix (csr = gene-major csc)
read_mat <- function(g, genes, cells) {
  enc <- attrs(g)[["encoding-type"]]
  if (is.null(enc) || enc == "array") {  # dense, stored cells x genes
    m <- as(t(h5read(f, g)), "CsparseMatrix")
  } else {
    sh <- as.integer(attrs(g)[["shape"]])
    i <- as.integer(h5read(f, paste0(g, "/indices"))); p <- as.integer(h5read(f, paste0(g, "/indptr")))
    x <- as.numeric(h5read(f, paste0(g, "/data"))); d <- if (enc == "csr_matrix") rev(sh) else sh
    # anndata allows unsorted indices within a row; dgCMatrix does not: sparseMatrix sorts
    # (extended/PLAN.md step 3 plumbing fix 3)
    m <- tryCatch(new("dgCMatrix", i = i, p = p, x = x, Dim = d),
                  error = function(e) sparseMatrix(i = i, p = p, x = x, dims = d, index1 = FALSE, repr = "C"))
    rm(i, p, x)
    if (enc == "csc_matrix") m <- t(m)
  }
  dimnames(m) <- list(genes, cells)
  m
}

obs <- read_df("/obs"); var <- read_df("/var")
cells <- rownames(obs); genes <- rownames(var)
data <- read_mat("/X", genes, cells)
counts <- data
if (has("/raw/X")) {
  rg <- rownames(read_df("/raw/var"))
  counts <- read_mat("/raw/X", rg, cells)
  if (!identical(rg, genes)) {  # schema 5 says raw.var == var; logged if not
    message("raw/var differs from var: counts subset to var's ", length(genes), " of ", length(rg), " genes")
    counts <- counts[genes, ]
  }
}
assay <- CreateAssay5Object(counts = counts, data = data)
assay <- AddMetaData(assay, var)
obj <- CreateSeuratObject(counts = assay, assay = "RNA", meta.data = obs)
# CreateSeuratObject recomputes nCount_/nFeature_<assay>, overwriting a source obs column of
# that name: put the source's values back (extended/PLAN.md, step 3 plumbing fix 1)
obj@meta.data[names(obs)] <- obs
saveRDS(obj, out)
cat(out, ncol(obj), "cells", nrow(obj), "genes\n")
