suppressPackageStartupMessages({library(slingshot); library(data.table); library(Matrix)})
set.seed(12345)
Sys.setenv(OMP_NUM_THREADS = parallel::detectCores())
data.table::setDTthreads(parallel::detectCores())
SCRIPT_PATH <- normalizePath(
  sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1L]),
  winslash = "/", mustWork = TRUE
)
PIPELINE_DIR <- normalizePath(file.path(dirname(SCRIPT_PATH), ".."), winslash = "/", mustWork = TRUE)
WORK_ROOT <- normalizePath(file.path(PIPELINE_DIR,".."),winslash="/",mustWork=TRUE)
Sys.setenv(LLR_WORK_ROOT=WORK_ROOT,LLR_PIPELINE=PIPELINE_DIR,LLR_ASSETS=WORK_ROOT)
source(file.path(PIPELINE_DIR, "config", "00_config.R"))
source(file.path(dirname(SCRIPT_PATH), "_benchmark_io.R"))

args <- commandArgs(trailingOnly = TRUE)
ds <- args[1]
outdir <- args[2]
if (is.na(ds)) stop("Usage: Rscript run_slingshot.R <dataset> <outdir>")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Slingshot runs in a multi-dimensional PCA space computed transparently here
# (the slingshot vignette runs on PCA), instead of an externally precomputed 2D
# embedding. No subsampling, no cluster merging/capping.
rd <- EnsureExport(ds, need = c("expr_lognorm.mtx", "meta.csv", "cells.txt", "genes.txt"))
expr <- as(Matrix::readMM(file.path(rd, "expr_lognorm.mtx")), "CsparseMatrix")  # genes x cells
meta <- fread(file.path(rd, "meta.csv"))
genes <- readLines(file.path(rd, "genes.txt"), warn = FALSE)
cells <- readLines(file.path(rd, "cells.txt"), warn = FALSE)
rownames(expr) <- genes
colnames(expr) <- cells

# Top variable genes -> PCA reducedDim (sparse row variance, no dense full matrix).
gene_var <- Matrix::rowMeans(expr^2) - Matrix::rowMeans(expr)^2
keep <- which(is.finite(gene_var) & gene_var > 1e-10)
top_n <- min(2000L, length(keep))
keep <- keep[order(gene_var[keep], decreasing = TRUE)[seq_len(top_n)]]
expr_hvg <- t(as.matrix(expr[keep, , drop = FALSE]))  # cells x genes
n_pc <- min(30L, nrow(expr_hvg) - 1L, ncol(expr_hvg) - 1L)
pca <- prcomp(expr_hvg, rank. = n_pc, center = TRUE, scale. = FALSE)$x
rownames(pca) <- cells

cl <- PickClusterLabels(meta, cells)
if (is.null(cl)) {
  cl <- ComputeFallbackClusters(expr, cells)
  cat("slingshot: using fallback PCA+kmeans clusters\n")
}

# Official slingshot workflow on the full cluster set: getLineages() then
# getCurves() with vendor defaults (shrink=TRUE, maxit=15, reassign/reweight=TRUE).
# Prefer vendor default dist.method="slingshot" (Mahalanobis). If singular
# covariance blocks lineage discovery, fall back to "simple" (Euclidean) and
# record the reason in slingshot_dist_method.txt.
lin <- NULL
dist_used <- "slingshot"
fallback_reason <- NA_character_
try_sl <- try(
  getLineages(pca, clusterLabels = factor(cl), reducedDim = "PCA", dist.method = "slingshot"),
  silent = TRUE
)
if (inherits(try_sl, "try-error")) {
  fallback_reason <- conditionMessage(attr(try_sl, "condition"))
  cat("slingshot: dist.method=slingshot failed:", fallback_reason, "\n")
  cat("slingshot: retrying with dist.method=simple (Euclidean centroids)\n")
  dist_used <- "simple"
  lin <- getLineages(pca, clusterLabels = factor(cl), reducedDim = "PCA",
                     dist.method = "simple")
} else {
  lin <- try_sl
}
sce <- getCurves(lin)
pt <- slingPseudotime(sce, na = FALSE)
cw <- slingCurveWeights(sce, as.probs = TRUE)
valid <- is.finite(pt) & is.finite(cw)
weighted_pt <- pt * cw
weighted_pt[!valid] <- 0
valid_weights <- cw
valid_weights[!valid] <- 0
weighted_sum <- rowSums(weighted_pt)
weight_sum <- rowSums(valid_weights)
global_pt <- weighted_sum / weight_sum
fwrite(data.table(cell_id = cells, pseudotime = global_pt), file.path(outdir, "slingshot_pseudotime.csv"))
fwrite(as.data.table(pt), file.path(outdir, "slingshot_lineage_pseudotime.csv"))
fwrite(as.data.table(cw), file.path(outdir, "slingshot_weights.csv"))
meta_lines <- c(
  paste0("dist.method=", dist_used),
  paste0("n_lineages=", ncol(pt)),
  paste0("n_clusters=", length(unique(cl))),
  "pseudotime_collapse=curve_weighted_mean"
)
if (!is.na(fallback_reason)) {
  meta_lines <- c(meta_lines, paste0("fallback_reason=", fallback_reason))
}
writeLines(meta_lines, file.path(outdir, "slingshot_dist_method.txt"))
cat("OK:", dist_used, ":", ncol(pt), "lineages,", nrow(pt), "cells\n")
