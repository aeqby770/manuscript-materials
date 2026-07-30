suppressPackageStartupMessages({library(Matrix); library(data.table); library(FateID)})
SCRIPT_PATH <- normalizePath(
  sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1L]),
  winslash = "/", mustWork = TRUE
)
PIPELINE_DIR <- normalizePath(file.path(dirname(SCRIPT_PATH), ".."), winslash = "/", mustWork = TRUE)
source(file.path(PIPELINE_DIR, "config", "00_config.R"))
source(file.path(dirname(SCRIPT_PATH), "_benchmark_io.R"))
source(file.path(dirname(SCRIPT_PATH), "_mem_guard.R"))

set.seed(12345)
Sys.setenv(OMP_NUM_THREADS = parallel::detectCores())
data.table::setDTthreads(parallel::detectCores())

HVG_TOP <- 2000L
KM_HVG <- 500L

args <- commandArgs(trailingOnly = TRUE)
ds <- args[1]
outdir <- args[2]
if (is.na(ds)) stop("Usage: Rscript run_fateid.R <dataset> <outdir>")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

rd <- file.path(CONFIG$dir$r_export, ds)
if (!file.exists(file.path(rd, "expr_lognorm.mtx"))) {
  sys::exec_wait(CONFIG$exe$rscript,
    c(file.path(dirname(SCRIPT_PATH), "00_export_R.R"), ds))
}

expr_s <- readMM(file.path(rd, "expr_lognorm.mtx"))
genes <- readLines(file.path(rd, "genes.txt"), warn = FALSE)
cells <- readLines(file.path(rd, "cells.txt"), warn = FALSE)
rownames(expr_s) <- make.unique(genes)
colnames(expr_s) <- cells[seq_len(ncol(expr_s))]

meta <- fread(file.path(rd, "meta.csv"))
cl <- PickClusterLabels(meta, cells)

gene_var <- Matrix::rowMeans(expr_s^2) - Matrix::rowMeans(expr_s)^2
keep_var <- which(is.finite(gene_var) & gene_var > 1e-10)
if (is.null(cl)) {
  top_km <- min(KM_HVG, length(keep_var))
  keep_km <- keep_var[order(gene_var[keep_var], decreasing = TRUE)[seq_len(top_km)]]
  expr_km <- t(as.matrix(expr_s[keep_km, , drop = FALSE]))
  n_cl <- max(3L, min(10L, as.integer(ncol(expr_s) / 100L)))
  km <- stats::kmeans(expr_km, centers = n_cl, nstart = 20, iter.max = 100)
  cl <- as.character(km$cluster)
  names(cl) <- cells
  rm(expr_km)
  gc()
  cat("FateID: fallback kmeans on ", top_km, " HVG, k=", n_cl, "\n", sep = "")
}

top_n <- min(HVG_TOP, length(keep_var))
keep <- keep_var[order(gene_var[keep_var], decreasing = TRUE)[seq_len(top_n)]]
mem_require_under_gb(
  mem_dense_gb(ncol(expr_s), length(keep), copies = 3),
  what = paste0("FateID/", ds)
)
if (nrow(expr_s) > HVG_TOP) {
  cat("FateID: HVG ", length(keep), "/", nrow(expr_s), " genes\n", sep = "")
}

expr <- as.matrix(expr_s[keep, , drop = FALSE])
rm(expr_s)
gc()
colnames(expr) <- cells

if (is.null(names(cl))) names(cl) <- cells
y_chr <- as.character(cl[cells])
y_chr[y_chr == ""] <- NA_character_
lev <- sort(unique(y_chr[!is.na(y_chr)]))
if (length(lev) < 2L) stop("FateID: need >=2 cluster labels, got ", length(lev))
y_map <- stats::setNames(seq_along(lev), lev)
y <- as.integer(y_map[y_chr])
names(y) <- cells

expr_t <- t(expr)
pc <- stats::prcomp(expr_t, rank. = 2)$x[, 1]
cluster_pc <- tapply(pc, y, mean)
tar <- as.integer(names(cluster_pc)[order(cluster_pc)[c(1, length(cluster_pc))]])
if (any(is.na(tar))) stop("FateID: invalid target clusters (check cluster label mapping)")

fb <- FateID::fateBias(
  expr,
  y,
  tar,
  z = NULL,
  minnr = 5,
  minnrh = 10,
  adapt = TRUE,
  confidence = 0.75,
  nbfactor = 5,
  use.dist = FALSE,
  seed = 12345,
  nbtree = NULL,
  verbose = FALSE
)

fateid_row_index <- as.integer(rownames(fb$probs)) + 1L
probs <- as.data.table(fb$probs)
probs[, cell_id := cells[fateid_row_index]]
probs[, fateid_cluster := y_chr[fateid_row_index]]
setcolorder(probs, c("cell_id", "fateid_cluster", setdiff(names(probs), c("cell_id", "fateid_cluster"))))
fwrite(probs, file.path(outdir, "fateid_bias.csv"))
saveRDS(fb, file.path(outdir, "fateid_bias.rds"))
cat("OK:", nrow(probs), "cells, targets=", paste(tar, collapse = ","), "\n")
