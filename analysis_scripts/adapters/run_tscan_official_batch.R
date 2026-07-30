
# path bootstrap (relative; set LLR_ASSETS or LLR_PIPELINE to override)
.llr_cmd <- commandArgs(trailingOnly = FALSE)
.llr_file <- sub("^--file=", "", grep("^--file=", .llr_cmd, value = TRUE)[1])
if (!length(.llr_file) || is.na(.llr_file) || !nzchar(.llr_file)) {
  stop("Run this script via Rscript so --file= is set")
}
.llr_here <- dirname(normalizePath(.llr_file, winslash = "/", mustWork = TRUE))
source(file.path(.llr_here, "../config/llr_bootstrap.R"), local = FALSE)

set.seed(12345)
Sys.setenv(OMP_NUM_THREADS = parallel::detectCores())
data.table::setDTthreads(parallel::detectCores())

llr_source_config()
TSCAN_CONFIG <- list(
  base_output = file.path(CONFIG$dir$algorithm_out, "TSCAN_official"),
  base_export = CONFIG$dir$r_export,
  implementation_class = "OFFICIAL_TSCAN_API_EXPRMCLUST_EII",
  risk_class = "MEDIUM"
)

library(Matrix)
library(data.table)
library(TSCAN)

split_ds <- function(x) {
  x <- trimws(x)
  x <- x[nzchar(x)]
  unique(x)
}

args <- commandArgs(trailingOnly = TRUE)
target_datasets <- NULL
for (a in args) {
  if (grepl("^--datasets-file=", a)) {
    f <- sub("^--datasets-file=", "", a)
    if (!file.exists(f)) stop("datasets file not found: ", f)
    target_datasets <- split_ds(readLines(f, warn = FALSE))
  }
  if (grepl("^--datasets=", a)) {
    raw <- sub("^--datasets=", "", a)
    sep <- if (grepl(";", raw, fixed = TRUE)) ";" else ","
    target_datasets <- split_ds(unlist(strsplit(raw, sep, fixed = TRUE)))
  }
}
if (is.null(target_datasets) && nzchar(Sys.getenv("TSCAN_DATASETS"))) {
  target_datasets <- split_ds(strsplit(Sys.getenv("TSCAN_DATASETS"), ",", fixed = TRUE)[[1]])
}
if (is.null(target_datasets)) {
  target_datasets <- names(datasets)
}

run_tscan <- function(ds) {
  cat(sprintf("\n[TSCAN_official] %s ... ", ds))
  if (!ds %in% names(datasets)) {
    cat("SKIP (unknown dataset)\n")
    stop("TSCAN_official: unknown dataset ", ds)
  }
  rd <- file.path(TSCAN_CONFIG$base_export, ds)
  od <- file.path(TSCAN_CONFIG$base_output, ds)
  marker <- file.path(od, "tscan_official_pseudotime.csv")
  if (file.exists(marker) && file.info(marker)$size > 100) {
    cat("CACHED\n")
    return(invisible(marker))
  }
  dir.create(od, recursive = TRUE, showWarnings = FALSE)

  if (!file.exists(file.path(rd, "expr_lognorm.mtx"))) {
    status <- system2(CONFIG$exe$rscript,
      args = c(file.path(llr_pipeline_dir(), "adapters/00_export_R.R"), ds),
      stdout = file.path(CONFIG$dir$log, sprintf("tscan_official_export_%s.log", ds)),
      stderr = file.path(CONFIG$dir$log, sprintf("tscan_official_export_%s.err", ds)))
    if (status != 0) { cat("FAIL (export)\n"); return(invisible(NULL)) }
  }

  expr_s <- readMM(file.path(rd, "expr_lognorm.mtx"))
  genes <- fread(file.path(rd, "genes.txt"), header = FALSE)[["V1"]]
  cells <- fread(file.path(rd, "cells.txt"), header = FALSE)[["V1"]]
  meta <- fread(file.path(rd, "meta.csv"))
  rownames(expr_s) <- genes
  colnames(expr_s) <- cells

  # Drop only zero-variance genes (no top-N cap): exprmclust() does its own PCA.
  gene_var <- Matrix::rowMeans(expr_s^2) - Matrix::rowMeans(expr_s)^2
  keep <- which(is.finite(gene_var) & gene_var > 1e-10)
  expr <- as.matrix(expr_s[keep, , drop = FALSE])
  rm(expr_s); gc()

  mclust_res <- exprmclust(expr, clusternum = NULL, modelNames = "EII")
  saveRDS(mclust_res, file.path(od, "tscan_official_mclust.rds"))

  order_res <- TSCANorder(mclust_res, orderonly = FALSE, flip = TRUE, listbranch = FALSE)
  pt <- order_res$Pseudotime
  names(pt) <- order_res$sample_name

  pt_df <- data.frame(cell = names(pt), pseudotime = pt, stringsAsFactors = FALSE)
  fwrite(pt_df, marker)

  gt_vec <- meta[["true_pseudotime"]]
  names(gt_vec) <- as.character(seq_len(nrow(meta)))
  common_cells <- intersect(names(pt), names(gt_vec))
  n_valid <- length(common_cells)
  if (n_valid > 2) {
    pt_vals <- pt[common_cells]
    gt_vals <- gt_vec[common_cells]
    spear <- abs(cor(pt_vals, gt_vals, method = "spearman"))
    pear <- cor(pt_vals, gt_vals, method = "pearson")
    rmse_val <- sqrt(mean((pt_vals - gt_vals)^2))
  } else {
    spear <- pear <- rmse_val <- NA_real_
  }

  summary <- data.table(
    algorithm = "TSCAN_official",
    dataset = ds,
    n_cells = ncol(expr),
    n_genes = nrow(expr),
    n_valid = n_valid,
    n_total = nrow(meta),
    n_clusters = mclust_res$clusternum,
    Spearman = spear,
    Pearson = pear,
    RMSE = rmse_val,
    implementation_class = TSCAN_CONFIG$implementation_class,
    misrepresentation_risk = TSCAN_CONFIG$risk_class,
    metric_layer = "pseudotime"
  )
  fwrite(summary, file.path(od, "tscan_official_summary.csv"))
  cat("OK\n")
  invisible(marker)
}

results <- data.table(algorithm = character(), dataset = character(), status = character())
for (ds in target_datasets) {
  out <- run_tscan(ds)
  # Status must reflect the real pseudotime output, not be hard-coded "OK".
  ok <- is.character(out) && length(out) == 1L && file.exists(out) && file.info(out)$size > 100
  results <- rbind(results, data.table(
    algorithm = "TSCAN_official", dataset = ds, status = if (ok) "OK" else "FAIL"
  ), fill = TRUE)
}
cat("\n===== TSCAN_official Summary =====\n")
print(results[, .N, by = .(algorithm, status)])
