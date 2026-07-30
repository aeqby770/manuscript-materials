ADAPTER_CONFIG <- list(
  threads = parallel::detectCores(),
  seed = 12345L,
  output_name = "scorpius_pseudotime.csv"
)
set.seed(ADAPTER_CONFIG$seed)
Sys.setenv(OMP_NUM_THREADS = ADAPTER_CONFIG$threads)
data.table::setDTthreads(ADAPTER_CONFIG$threads)

suppressPackageStartupMessages({
  library(Matrix)

# path bootstrap (relative; set LLR_ASSETS or LLR_PIPELINE to override)
.llr_cmd <- commandArgs(trailingOnly = FALSE)
.llr_file <- sub("^--file=", "", grep("^--file=", .llr_cmd, value = TRUE)[1])
if (!length(.llr_file) || is.na(.llr_file) || !nzchar(.llr_file)) {
  stop("Run this script via Rscript so --file= is set")
}
.llr_here <- dirname(normalizePath(.llr_file, winslash = "/", mustWork = TRUE))
source(file.path(.llr_here, "../config/llr_bootstrap.R"), local = FALSE)

  library(data.table)
  library(SCORPIUS)
  library(matrixStats)
})
llr_source_config()
llr_source_rel("adapters/_mem_guard.R")
llr_source_rel("adapters/_benchmark_io.R")
args <- commandArgs(trailingOnly = TRUE)
ds <- args[1]
outdir <- args[2]
if (is.na(ds)) stop("Usage: Rscript run_scorpius.R <dataset> <outdir>")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
rd <- EnsureExport(ds, need = c("expr.mtx", "cells.txt"))
# Official SCORPIUS: rows = cells (samples), cols = genes — see vignettes/anndata.Rmd & R/dimensionality_reduction.R
x <- Matrix::readMM(file.path(rd, "expr.mtx"))
x <- Matrix::t(x)
if (!inherits(x, "CsparseMatrix")) x <- methods::as(x, "CsparseMatrix")
cells <- readLines(file.path(rd, "cells.txt"), warn = FALSE, encoding = "UTF-8")
if (nrow(x) != length(cells)) {
  stop("SCORPIUS: expr rows (cells) ", nrow(x), " != cells.txt ", length(cells))
}
if (anyNA(cells) || any(!nzchar(cells)) || anyDuplicated(cells)) {
  stop("SCORPIUS: cells.txt must contain one complete, unique cell_id sequence")
}
rownames(x) <- cells
ndim <- min(3L, nrow(x) - 1L)
if (ndim < 1L) stop("SCORPIUS: need >= 2 cells")
# reduce_dimensionality() builds an O(n_cells^2) cell-cell distance matrix.
mem_require_under_gb(mem_square_gb(nrow(x), copies = 2), what = paste0("SCORPIUS/", ds))
space <- reduce_dimensionality(x, dist = "spearman", ndim = ndim)
res <- infer_trajectory(space)
pt <- as.numeric(res$time)
if (length(pt) != nrow(x)) {
  stop("SCORPIUS: pseudotime length (", length(pt), ") != cells (", nrow(x), ")")
}
if (any(!is.finite(pt)) || diff(range(pt)) <= 0) {
  stop("SCORPIUS: pseudotime must be finite for every cell and have non-zero range")
}
out_file <- file.path(outdir, ADAPTER_CONFIG$output_name)
fwrite(
  data.table(cell_id = cells, pseudotime = pt),
  out_file
)
written <- fread(
  out_file,
  encoding = "UTF-8",
  colClasses = c(cell_id = "character", pseudotime = "numeric")
)
if (!identical(names(written), c("cell_id", "pseudotime")) ||
    !identical(as.character(written$cell_id), cells) ||
    anyDuplicated(written$cell_id) ||
    any(!is.finite(written$pseudotime)) ||
    diff(range(written$pseudotime)) <= 0) {
  stop("SCORPIUS: written CSV failed the full-cell pseudotime contract")
}
cat("OK:", length(pt), "cells\n")
