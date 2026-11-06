# Snapshot a Seurat .rds (readRDS + SeuratObject) or a BPCells-backed
# .h5seurat (rhdf5 + BPCells) into the layout described in snapshot.py.
# Base R I/O only: the readers env has no jsonlite.
#   Rscript snapshot.R <file.rds|file.h5seurat> <out_dir>
a <- commandArgs(TRUE); f <- a[1]; out <- a[2]
for (s in c("obs", "var", "mats")) dir.create(file.path(out, s), recursive = TRUE, showWarnings = FALSE)

column <- function(dir, k, name, x) {
  p <- file.path(dir, k)
  writeLines(name, paste0(p, ".name"))
  na <- is.na(x)
  writeBin(as.raw(na), paste0(p, ".na"))
  cls <- if (is.factor(x)) "categorical" else if (is.logical(x)) "bool" else
    if (is.integer(x)) "int" else if (is.double(x)) "float" else "string"
  if (cls == "categorical") {
    writeBin(ifelse(na, -1L, as.integer(x) - 1L), paste0(p, ".codes"), size = 4, endian = "little")
    writeLines(levels(x), paste0(p, ".levels.txt"))
    writeLines(as.character(as.integer(is.ordered(x))), paste0(p, ".ordered"))
  } else if (cls == "string") {
    writeLines(ifelse(na, "", as.character(x)), paste0(p, ".txt"))
  } else {
    writeBin(as.double(ifelse(na, 0, x)), paste0(p, ".f64"), size = 8, endian = "little")
  }
  writeLines(cls, paste0(p, ".class"))
}

columns <- function(dir, df) for (k in seq_along(df)) column(dir, k - 1L, names(df)[k], df[[k]])

# m: genes x cells. Its CSC arrays are the CSR arrays of cells x genes.
matrix_ <- function(k, name, m, dtype) {
  m <- as(as(m, "CsparseMatrix"), "dgCMatrix")
  p <- file.path(out, "mats", k)
  writeLines(name, paste0(p, ".name"))
  writeLines(paste(ncol(m), nrow(m)), paste0(p, ".shape"))
  writeLines(dtype, paste0(p, ".dtype"))
  writeLines(rownames(m), paste0(p, ".var_names.txt"))
  for (part in list(list("indptr", m@p), list("indices", m@i), list("data", m@x)))
    writeBin(as.double(part[[2]]), paste0(p, ".", part[[1]]), size = 8, endian = "little")
}

if (grepl("\\.h5seurat$", f)) {
  # Two layouts in this benchmark. scx: root cell.names, BPCells-packed
  # assays/<a>/<layer>, assays/<a>/meta.features. scConvert (Seurat v5 style):
  # cell names in meta.data/<_index>, dgCMatrix groups (data/indices/indptr +
  # dims) under assays/<a>/layers/<layer>, feature metadata in a root `var`.
  suppressPackageStartupMessages({ library(rhdf5); library(BPCells); library(Matrix) })
  l <- h5ls(f)
  path <- paste0(sub("/$", "", l$group), "/", l$name)
  has <- function(p) p %in% path
  attrs <- function(p) tryCatch(h5readAttributes(f, p), error = function(e) list())
  assays <- l$name[l$group == "/assays" & l$otype == "H5I_GROUP"]
  act <- attrs("/")[["active.assay"]]
  act <- if (is.null(act)) assays[1] else act
  idx <- attrs("/meta.data")[["_index"]]
  cells <- as.character(h5read(f, if (has("/cell.names")) "cell.names" else paste0("meta.data/", idx)))
  writeLines(cells, file.path(out, "obs_names.txt"))
  writeLines(as.character(h5read(f, paste0("assays/", act, "/features"))), file.path(out, "var_names.txt"))
  # A data.frame group. Column order from `colnames` (SeuratDisk) or
  # `column-order` (anndata); factors are levels + 1-based values (SeuratDisk)
  # or categories + 0-based codes, -1 = NA (anndata).
  read_df <- function(g) {
    if (!has(g)) return(data.frame())
    a <- attrs(g)
    cols <- as.character(c(a[["colnames"]], a[["column-order"]]))
    if (!length(cols)) cols <- setdiff(l$name[l$group == g], c(a[["_index"]], "_index", "__categories"))
    df <- lapply(cols, function(c) {
      q <- paste0(g, "/", c)
      if (l$otype[path == q] == "H5I_GROUP") {
        kids <- l$name[l$group == q]
        if ("levels" %in% kids) {
          lv <- as.character(h5read(f, paste0(q, "/levels")))
          factor(lv[as.integer(h5read(f, paste0(q, "/values")))], levels = lv)
        } else {
          lv <- as.character(h5read(f, paste0(q, "/categories")))
          codes <- as.integer(h5read(f, paste0(q, "/codes")))
          factor(ifelse(codes < 0, NA, lv[codes + 1L]), levels = lv)
        }
      } else if (has(paste0(g, "/__categories/", c))) {
        # anndata 0.7 categoricals, as copied from old h5ad inputs: 0-based codes
        # plus levels in <df>/__categories/<col>
        lv <- as.character(h5read(f, paste0(g, "/__categories/", c)))
        codes <- as.integer(h5read(f, q))
        factor(ifelse(codes < 0, NA, lv[codes + 1L]), levels = lv)
      } else {
        x <- h5read(f, q)
        # rhdf5 returns HDF5 booleans (h5py's FALSE/TRUE enum) as text
        if (l$dclass[path == q] == "ENUM" && all(as.character(x) %in% c("FALSE", "TRUE", NA)))
          as.logical(as.character(x)) else as.vector(x)
      }
    })
    names(df) <- cols
    as.data.frame(df, check.names = FALSE, stringsAsFactors = FALSE)
  }
  columns(file.path(out, "obs"), read_df("/meta.data"))
  vg <- paste0("/assays/", act, "/meta.features")
  columns(file.path(out, "var"), read_df(if (has(vg)) vg else "/var"))
  k <- 0L
  for (a in assays) {
    feats <- as.character(h5read(f, paste0("assays/", a, "/features")))
    base <- paste0("/assays/", a)
    groups <- c(l$name[l$group == base & l$otype == "H5I_GROUP"])
    if ("layers" %in% groups) groups <- c(setdiff(groups, "layers"),
                                          paste0("layers/", l$name[l$group == paste0(base, "/layers")]))
    dense <- l$name[l$group == paste0(base, "/layers") & l$otype == "H5I_DATASET"]
    layers <- c(groups, if (length(dense)) paste0("layers/", dense))
    for (lay in setdiff(unique(layers), "meta.features")) {
      g <- paste0(base, "/", lay)
      kids <- l$name[l$group == g]
      # A layer that can't be read as genes x cells with these names is skipped
      # with a message: the comparison then reports that matrix as not found,
      # and every other property is still checked.
      m <- tryCatch({
        m <- if (l$otype[path == g] == "H5I_DATASET") {
          x <- h5read(f, g)                      # dense
          if (!identical(dim(x), c(length(feats), length(cells)))) x <- t(x)
          as(Matrix(x, sparse = TRUE), "dgCMatrix")
        } else if (all(c("data", "indices", "indptr") %in% kids)) {
          new("dgCMatrix", i = as.integer(h5read(f, paste0(g, "/indices"))),
              p = as.integer(h5read(f, paste0(g, "/indptr"))),
              x = as.double(h5read(f, paste0(g, "/data"))), Dim = as.integer(attrs(g)[["dims"]]))
        } else as(open_matrix_hdf5(f, sub("^/", "", g)), "dgCMatrix")
        if (!identical(dim(m), c(length(feats), length(cells))))
          stop("dims ", paste(dim(m), collapse = " x "), ", but ", length(feats), " features x ", length(cells), " cells")
        dimnames(m) <- list(feats, cells)
        m
      }, error = function(e) { message("skipped layer ", g, ": ", conditionMessage(e)); NULL })
      if (is.null(m)) next
      matrix_(k, paste0(a, "/", sub("^layers/", "", lay)), m, "h5seurat"); k <- k + 1L
    }
  }
} else {
  suppressPackageStartupMessages(library(SeuratObject))
  o <- readRDS(f)
  # Update old objects; an object Update rejects (e.g. a malformed DimReduc) is
  # still read as is: the slots compared here don't depend on it.
  md <- o@meta.data
  o <- tryCatch(suppressMessages(suppressWarnings(UpdateSeuratObject(o))), error = function(e) o)
  # UpdateSeuratObject recomputes nCount_/nFeature_ columns: the reader reports what
  # the file stores, so the stored meta.data goes back
  if (identical(rownames(md), Cells(o))) o@meta.data <- md
  writeLines(Cells(o), file.path(out, "obs_names.txt"))
  writeLines(rownames(o), file.path(out, "var_names.txt"))
  columns(file.path(out, "obs"), o[[]])
  columns(file.path(out, "var"), o[[DefaultAssay(o)]][[]])
  k <- 0L
  for (a in Assays(o)) for (lay in Layers(o[[a]])) {
    if (grepl("^scale", lay)) next  # dense, HVG subset: not a cells x genes slot
    m <- LayerData(o, assay = a, layer = lay)
    if (!length(m)) next
    matrix_(k, paste0(a, "/", lay), m, class(m)[1]); k <- k + 1L
  }
}
