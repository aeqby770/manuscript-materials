# Official TSCAN API: exprmclust + TSCANorder (genes x cells, log-normalized).
# Reference: TSCAN package vignette; same logic as run_tscan_official_batch.R.
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
  library(TSCAN)
})
llr_source_config()
llr_source_rel("adapters/_mem_guard.R")
llr_source_rel("adapters/_benchmark_io.R")
set.seed(12345)
Sys.setenv(OMP_NUM_THREADS = parallel::detectCores())
data.table::setDTthreads(parallel::detectCores())
args <- commandArgs(trailingOnly = TRUE)
ds <- args[1]
outdir <- args[2]
TSCAN_HVG_LIMIT <- 2000L
if (is.na(ds)) stop("Usage: Rscript run_tscan.R <dataset> <outdir>")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
rd <- EnsureExport(ds, need = c("expr_lognorm.mtx", "genes.txt", "cells.txt"))
expr_s <- readMM(file.path(rd, "expr_lognorm.mtx"))
genes <- fread(file.path(rd, "genes.txt"), header = FALSE)[["V1"]]
cells <- fread(file.path(rd, "cells.txt"), header = FALSE)[["V1"]]
genes <- as.character(genes)
cells <- as.character(cells)
if (nrow(expr_s) != length(genes) || ncol(expr_s) != length(cells)) {
  stop("TSCAN: expression dimensions ", nrow(expr_s), " x ", ncol(expr_s),
       " do not match genes/cells ", length(genes), " x ", length(cells))
}
if (anyNA(cells) || any(!nzchar(cells)) || anyDuplicated(cells)) {
  stop("TSCAN: cells.txt must contain one complete, unique cell_id sequence")
}
rownames(expr_s) <- genes
colnames(expr_s) <- cells

# Drop zero-variance genes, then retain a fixed high-variance representation before
# TSCAN's PCA-based exprmclust step. This keeps every cell while avoiding a dense
# all-gene matrix on large panels.
gene_var <- Matrix::rowMeans(expr_s^2) - Matrix::rowMeans(expr_s)^2
keep <- which(is.finite(gene_var) & gene_var > 1e-10)
keep <- keep[order(gene_var[keep], decreasing = TRUE, method = "radix")]
keep <- head(keep, min(TSCAN_HVG_LIMIT, length(keep)))
# exprmclust() densifies the retained genes and runs its own PCA.
mem_require_under_gb(mem_dense_gb(ncol(expr_s), length(keep), copies = 2), what = paste0("TSCAN/", ds))
expr <- as.matrix(expr_s[keep, , drop = FALSE])
rm(expr_s); gc()
if (nrow(expr) < 10L || ncol(expr) < 10L) {
  stop("TSCAN: insufficient genes/cells after filtering (", nrow(expr), " x ", ncol(expr), ")")
}

mclust_res <- exprmclust(expr, clusternum = NULL, modelNames = "EII")
saveRDS(mclust_res, file.path(outdir, "tscan_mclust.rds"))
order_res <- TSCANorder(mclust_res, orderonly = FALSE, flip = TRUE, listbranch = FALSE)
if (!nrow(order_res)) stop("TSCAN: official backbone ordering is empty")

# TSCANorder() ranks one endpoint-to-endpoint backbone. Extend the same TSCAN
# PCA-center/MST geometry to every branch by projecting each cell onto its
# nearest incident MST edge, then measuring weighted tree distance from the
# root selected by the official flipped backbone.
cluster_names <- names(mclust_res$clusterid)
pca_names <- rownames(mclust_res$pcareduceres)
if (is.null(cluster_names) || is.null(pca_names) || anyDuplicated(cluster_names) ||
    anyDuplicated(pca_names) || anyNA(match(cells, cluster_names)) ||
    anyNA(match(cells, pca_names))) {
  stop("TSCAN: model cell IDs do not exactly cover canonical cells")
}

cell_cluster <- as.integer(mclust_res$clusterid[match(cells, cluster_names)])
cell_pca <- mclust_res$pcareduceres[match(cells, pca_names), , drop = FALSE]
centers <- mclust_res$clucenter
mst <- mclust_res$MSTtree
adjacency <- as.matrix(igraph::as_adjacency_matrix(mst, sparse = FALSE))
allowed_edges <- adjacency[cell_cluster, , drop = FALSE] > 0
center_distance_sq <-
  rowSums(cell_pca^2) +
  rep(rowSums(centers^2), each = nrow(cell_pca)) -
  2 * (cell_pca %*% t(centers))
center_distance_sq[!allowed_edges] <- Inf
edge_neighbor <- max.col(-center_distance_sq, ties.method = "first")

root_cluster <- as.integer(order_res$State[[1]])
root_distance <- as.numeric(igraph::distances(
  mst,
  v = root_cluster,
  to = igraph::V(mst),
  weights = igraph::E(mst)$weight
))
edge_parent <- ifelse(
  root_distance[cell_cluster] <= root_distance[edge_neighbor],
  cell_cluster,
  edge_neighbor
)
edge_child <- ifelse(edge_parent == cell_cluster, edge_neighbor, cell_cluster)
parent_center <- centers[edge_parent, , drop = FALSE]
edge_vector <- centers[edge_child, , drop = FALSE] - parent_center
edge_length <- sqrt(rowSums(edge_vector^2))
edge_fraction <- rowSums((cell_pca - parent_center) * edge_vector) / edge_length^2
pt <- root_distance[edge_parent] + edge_fraction * edge_length
pt <- (pt - min(pt)) / diff(range(pt))
if (any(!is.finite(pt)) || diff(range(pt)) <= 0) {
  stop("TSCAN: pseudotime must be finite for every cell and have non-zero range")
}

out_file <- file.path(outdir, "tscan_pseudotime.csv")
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
  stop("TSCAN: written CSV failed the full-cell pseudotime contract")
}
cat("OK:", length(pt), "cells,", nrow(centers), "clusters (official TSCAN full tree)\n")
