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
  suppressPackageStartupMessages({ library(rhdf5); library(BPCells) })
  l <- h5ls(f)
  assays <- l$name[l$group == "/assays" & l$otype == "H5I_GROUP"]
  act <- tryCatch(h5readAttributes(f, "/")[["active.assay"]], error = function(e) NULL)
  act <- if (is.null(act)) assays[1] else act
  writeLines(as.character(h5read(f, "cell.names")), file.path(out, "obs_names.txt"))
  writeLines(as.character(h5read(f, paste0("assays/", act, "/features"))), file.path(out, "var_names.txt"))
  # A data.frame group: one dataset per column, factors as a group of 1-based
  # `values` + `levels` (SeuratDisk layout). Column order from `colnames`.
  read_df <- function(g) {
    if (!any(paste0(l$group, "/", l$name) == g | paste0(l$group, l$name) == g)) return(data.frame())
    cols <- as.character(h5readAttributes(f, g)[["colnames"]])
    df <- lapply(cols, function(c) {
      if (l$otype[l$group == g & l$name == c] == "H5I_GROUP") {
        lv <- as.character(h5read(f, paste0(g, "/", c, "/levels")))
        factor(lv[as.integer(h5read(f, paste0(g, "/", c, "/values")))], levels = lv)
      } else as.vector(h5read(f, paste0(g, "/", c)))
    })
    names(df) <- cols
    as.data.frame(df, check.names = FALSE, stringsAsFactors = FALSE)
  }
  columns(file.path(out, "obs"), read_df("/meta.data"))
  columns(file.path(out, "var"), read_df(paste0("/assays/", act, "/meta.features")))
  k <- 0L
  for (a in assays) for (lay in l$name[l$group == paste0("/assays/", a) & l$otype == "H5I_GROUP"]) {
    if (lay == "meta.features") next
    g <- paste0("assays/", a, "/", lay)
    m <- open_matrix_hdf5(f, g)
    feats <- as.character(h5read(f, paste0("assays/", a, "/features")))
    cells <- as.character(h5read(f, "cell.names"))
    m <- as(m, "dgCMatrix"); dimnames(m) <- list(feats, cells)
    matrix_(k, paste0(a, "/", lay), m, "h5seurat:bpcells"); k <- k + 1L
  }
} else {
  suppressPackageStartupMessages(library(SeuratObject))
  o <- suppressMessages(suppressWarnings(UpdateSeuratObject(readRDS(f))))
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
