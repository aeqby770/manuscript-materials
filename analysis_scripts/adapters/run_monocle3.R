suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(monocle3)
})
set.seed(12345)
Sys.setenv(OMP_NUM_THREADS = parallel::detectCores())
data.table::setDTthreads(parallel::detectCores())
SCRIPT_PATH <- normalizePath(
  sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1L]),
  winslash = "/", mustWork = TRUE
)
PIPELINE_DIR <- normalizePath(file.path(dirname(SCRIPT_PATH), ".."), winslash = "/", mustWork = TRUE)
source(file.path(PIPELINE_DIR, "config", "00_config.R"))
source(file.path(dirname(SCRIPT_PATH), "_benchmark_io.R"))

args <- commandArgs(trailingOnly = TRUE)
ds_name <- args[1]
outdir <- args[2]
if (is.na(ds_name)) stop("Usage: Rscript run_monocle3.R <dataset> <outdir>")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

rd <- EnsureExport(ds_name, need = c("expr.mtx", "meta.csv", "genes.txt", "cells.txt"))
expr <- as(readMM(file.path(rd, "expr.mtx")), "CsparseMatrix")
meta <- as.data.frame(fread(file.path(rd, "meta.csv")))
genes <- readLines(file.path(rd, "genes.txt"))
cells <- readLines(file.path(rd, "cells.txt"))
rownames(expr) <- genes
colnames(expr) <- cells
idx <- match(cells, meta$cell_id)
if (any(is.na(idx))) idx <- match(cells, as.character(meta[[1L]]))
meta <- meta[idx, , drop = FALSE]
rownames(meta) <- cells

# No subsampling: monocle3 runs on the full cell set. The expression matrix is
# kept sparse (CsparseMatrix) so memory is bounded by the UMAP / principal graph
# rather than a dense O(n_cells^2) copy.

gm <- data.frame(gene_short_name = genes, row.names = genes)
cds <- new_cell_data_set(expr, cell_metadata = meta, gene_metadata = gm)
cds <- preprocess_cds(cds, num_dim = 30)
cds <- reduce_dimension(cds, reduction_method = "UMAP")
cds <- cluster_cells(cds, reduction_method = "UMAP", k = 50, cluster_method = "louvain")
cds <- learn_graph(cds, learn_graph_control = list(prune_graph = FALSE))
root_cells <- PickMonocle3RootCells(cds, meta)
# Restrict root cells to the current cell_data_set.
root_cells <- intersect(root_cells, colnames(cds))
if (length(root_cells) < 1L) {
  # Use the first cell of the largest cluster when no root cell remains.
  cl_cds <- NULL
  root_cl <- "first_cell"
  if (!is.null(cds@clusters) && "UMAP" %in% names(cds@clusters)) {
    cl_cds <- as.character(cds@clusters$UMAP$clusters)
  }
  if (is.null(cl_cds) || length(cl_cds) == 0L) {
    root_cells <- colnames(cds)[1L]
  } else {
    tab <- sort(table(cl_cds), decreasing = TRUE)
    root_cl <- names(tab)[1L]
    root_cells <- colnames(cds)[cl_cds == root_cl]
  }
  cat("monocle3: root_cells not found, using fallback cluster", root_cl, "\n")
}

# order_cells() leaves every partition without a root at Inf. Keep the selected
# roots and add one expression-only anchor for each uncovered partition.
partition <- as.character(partitions(cds, reduction_method = "UMAP")[cells])
detected <- as.numeric(Matrix::colSums(expr > 0))
root_partitions <- unique(partition[match(root_cells, cells)])
missing_partitions <- setdiff(unique(partition), root_partitions)
partition_anchors <- data.table(
  cell = cells,
  partition = partition,
  detected = detected,
  position = seq_along(cells)
)[partition %in% missing_partitions][
  order(partition, -detected, position)
][!duplicated(partition), cell]
root_cells <- unique(c(root_cells, partition_anchors))
cat(
  "monocle3 partitions:", length(unique(partition)),
  "supplemental expression roots:", length(partition_anchors), "\n"
)
cds <- order_cells(cds, root_cells = root_cells)
pt <- cds@principal_graph_aux$UMAP$pseudotime[cells]
stopifnot(
  identical(names(pt), cells),
  !anyDuplicated(cells),
  all(is.finite(pt)),
  diff(range(pt)) > 0
)
out_file <- file.path(outdir, "monocle3_pseudotime.csv")
fwrite(data.table(cell_id = cells, pseudotime = as.numeric(pt)), out_file)
written <- fread(out_file, colClasses = c(cell_id = "character"))
stopifnot(
  identical(names(written), c("cell_id", "pseudotime")),
  identical(written$cell_id, cells),
  !anyDuplicated(written$cell_id),
  all(is.finite(written$pseudotime)),
  diff(range(written$pseudotime)) > 0
)
# Scale-extension evaluation consumes the verified full-cell pseudotime CSV.
switch(
  as.character(ds_name %in% BENCHMARK_SCALE_EXTENSION_DATASETS),
  "TRUE" = NULL,
  "FALSE" = saveRDS(cds, file.path(outdir, "monocle3_cds.rds"))
)
cat("OK:", length(pt), "cells,", sum(is.finite(pt)), "finite\n")
