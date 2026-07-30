suppressPackageStartupMessages({library(Matrix);library(data.table);library(SLICER);library(lle)})

# path bootstrap (relative; set LLR_ASSETS or LLR_PIPELINE to override)
.llr_cmd <- commandArgs(trailingOnly = FALSE)
.llr_file <- sub("^--file=", "", grep("^--file=", .llr_cmd, value = TRUE)[1])
if (!length(.llr_file) || is.na(.llr_file) || !nzchar(.llr_file)) {
  stop("Run this script via Rscript so --file= is set")
}
.llr_here <- dirname(normalizePath(.llr_file, winslash = "/", mustWork = TRUE))
source(file.path(.llr_here, "../config/llr_bootstrap.R"), local = FALSE)

llr_source_config()
llr_source_rel("adapters/_mem_guard.R")
args<-commandArgs(trailingOnly=TRUE);ds<-args[1];outdir<-args[2]
if(is.na(ds))stop("Usage: Rscript run_slicer.R <dataset> <outdir>")
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
rd<-file.path(CONFIG$dir$r_export,ds)
if(!file.exists(file.path(rd,"expr_lognorm.mtx")))
  sys::exec_wait(CONFIG$exe$rscript,
    c(file.path(llr_pipeline_dir(), "adapters/00_export_R.R"),ds))

# Official SLICER workflow (README.md):
# genes = select_genes(traj)
# k = select_k(traj[,genes], kmin=5)
# traj_lle = lle(traj[,genes], m=2, k)$Y
# traj_graph = conn_knn_graph(traj_lle, 5)
# ends = find_extreme_cells(traj_graph, traj_lle)
# start = 1
# cells_ordered = cell_order(traj_graph, start)

expr_s <- readMM(file.path(rd, "expr_lognorm.mtx"))
genes_all <- readLines(file.path(rd, "genes.txt"), warn = FALSE)
cells <- readLines(file.path(rd, "cells.txt"), warn = FALSE)
rownames(expr_s) <- genes_all
colnames(expr_s) <- cells

# No subsampling: SLICER runs on the full cell set (official workflow).
# Drop only zero-variance genes (no top-N cap) — SLICER's own select_genes()
# performs the gene selection, so it receives the full informative gene set.
gene_var <- Matrix::rowMeans(expr_s^2) - Matrix::rowMeans(expr_s)^2
keep <- which(is.finite(gene_var) & gene_var > 1e-10)
# SLICER densifies the kept genes and builds an O(n_cells^2) kNN graph / LLE.
mem_require_under_gb(max(mem_dense_gb(ncol(expr_s), length(keep), copies = 2),
                        mem_square_gb(ncol(expr_s))), what = paste0("SLICER/", ds))
expr <- as.matrix(expr_s[keep, , drop = FALSE])
rm(expr_s); gc()

# Transpose to cells x genes (SLICER expects samples as rows)
expr_t <- t(expr)

# Step 1: select_genes (official API)
genes_idx <- SLICER::select_genes(expr_t)
if (length(genes_idx) < 10) {
  # select_genes() is SLICER's gene-selection step. Keep its input contract rather
  # than replacing it with a generic variance filter.
  stop("SLICER: select_genes() returned too few genes (", length(genes_idx),
       "); the official gene-selection step returned too few genes.")
}
expr_sel <- expr_t[, genes_idx, drop = FALSE]

# Step 2: select_k (official API; full search range, SLICER default kmax=500)
kmax_search <- min(500L, nrow(expr_sel) - 2L)
k <- SLICER::select_k(expr_sel, kmin = 5, kmax = kmax_search, by = 5)

# Step 3: lle (official API)
traj_lle <- lle::lle(expr_sel, m = 2, k = k)$Y

# Step 4: conn_knn_graph (official API)
traj_graph <- SLICER::conn_knn_graph(traj_lle, k = 5)

# Step 5: find_extreme_cells (official API)
ends <- SLICER::find_extreme_cells(traj_graph, traj_lle)
start <- ends[1]

# Step 6: cell_order (official API) - returns ordered cell indices.
cells_ordered <- SLICER::cell_order(traj_graph, start)

# Convert the returned cell order to pseudotime in [0, 1].
pt <- rep(NA_real_, nrow(expr_t))
pt[cells_ordered] <- seq(0, 1, length.out = length(cells_ordered))

fwrite(data.table(cell_id=cells, pseudotime=pt), file.path(outdir,"slicer_pseudotime.csv"))
cat("OK:", length(pt), "cells, k=", k, ", selected_genes=", length(genes_idx), "\n")
