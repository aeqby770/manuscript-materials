# Task E: blind-root vs user-root (GT earliest cells) Spearman on 6 real GT datasets
# path bootstrap (relative; set LLR_ASSETS to override)
.llr_cmd <- commandArgs(trailingOnly = FALSE)
.llr_file <- sub("^--file=", "", grep("^--file=", .llr_cmd, value = TRUE)[1])
if (!length(.llr_file) || is.na(.llr_file) || !nzchar(.llr_file)) {
  stop("Run via Rscript so --file= is set")
}
.llr_here <- dirname(normalizePath(.llr_file, winslash = "/", mustWork = TRUE))
source(file.path(.llr_here, "../config/llr_bootstrap.R"), local = FALSE)
.llr_assets <- llr_assets_dir()
.llr_config_env <- new.env(parent = globalenv())
sys.source(
  file.path(.llr_assets, "pipeline", "config", "00_config.R"),
  envir = .llr_config_env
)
PROJECT_CONFIG <- .llr_config_env$CONFIG
rm(.llr_config_env)

.libPaths(c(
  Sys.getenv("LLR064_RLIB", unset = "F:/biotrade/0.64LLR/rlib"),
  .libPaths()
))
set.seed(12345)
Sys.setenv(OMP_NUM_THREADS = as.character(parallel::detectCores()))
data.table::setDTthreads(parallel::detectCores())

library(data.table)
library(Matrix)
library(LineaLinkR)
stopifnot(identical(
  as.character(packageVersion("LineaLinkR")),
  PROJECT_CONFIG$llr$package_version
))

CONFIG <- list(
  paths = list(
    workspace = PROJECT_CONFIG$dir$linealinkr,
    ref = PROJECT_CONFIG$dir$reference,
    out = file.path(PROJECT_CONFIG$dir$supplement_tables, "user_root_contrast")
  ),
  contract_datasets = c(
    "lung", "dentategyrus", "monocle3_worm_l2",
    "gastrulation_erythroid", "zebrafish_neural_dev", "larry_invitro"
  ),
  seed = 12345,
  knn_k = 50L
)
CONFIG$datasets <- CONFIG$contract_datasets

dataset_arg <- grep("^--datasets=", commandArgs(trailingOnly = TRUE), value = TRUE)
selective_run <- length(dataset_arg) == 1L
if (selective_run) {
  CONFIG$datasets <- strsplit(
    sub("^--datasets=", "", dataset_arg),
    ",",
    fixed = TRUE
  )[[1L]]
}
stopifnot(
  length(CONFIG$datasets) > 0L,
  !anyDuplicated(CONFIG$datasets),
  all(CONFIG$datasets %in% CONFIG$contract_datasets)
)

dir.create(CONFIG$paths$out, recursive = TRUE, showWarnings = FALSE)

.extract_ordinal_time <- function(vals) {
  s <- as.character(vals)
  num <- as.numeric(sub(".*?(-?[0-9]+(?:\\.[0-9]+)?).*", "\\1", s, perl = TRUE))
  num[!grepl("[0-9]", s)] <- NA_real_
  num
}

RebuildNnFromSparse <- function(ds, n_cells, input_manifest) {
  dataset_dir <- file.path(CONFIG$paths$workspace, ds)
  dist_path <- file.path(dataset_dir, "knn_dist.mtx")
  conn_path <- file.path(dataset_dir, "knn_conn.mtx")
  knn_manifest <- readRDS(file.path(dataset_dir, "knn_precompute_manifest.rds"))
  manifest_matches_rds <- identical(
    knn_manifest[names(input_manifest)],
    input_manifest
  )
  if (!manifest_matches_rds) {
    stop("kNN decoded-content manifest does not match current RDS input for ", ds)
  }
  if (!identical(as.integer(knn_manifest$knn_neighbors), as.integer(CONFIG$knn_k))) {
    stop("kNN manifest has unexpected k for ", ds)
  }
  knn_dist <- as(Matrix::readMM(dist_path), "CsparseMatrix")
  knn_conn <- as(Matrix::readMM(conn_path), "CsparseMatrix")
  knn_contract_ok <- identical(dim(knn_dist), c(n_cells, n_cells)) &&
    identical(dim(knn_conn), c(n_cells, n_cells)) &&
    all(is.finite(knn_dist@x)) && all(knn_dist@x >= 0) &&
    all(is.finite(knn_conn@x)) && all(knn_conn@x >= 0) &&
    all(Matrix::rowSums(knn_conn) > 0)
  if (!knn_contract_ok) {
    stop("kNN decoded matrix contract does not match current RDS for ", ds)
  }
  k <- as.integer(min(CONFIG$knn_k, max(2L, n_cells - 1L)))
  sm <- Matrix::summary(knn_dist)
  sm <- sm[order(sm$i, sm$x), , drop = FALSE]
  nn_idx <- matrix(0L, nrow = n_cells, ncol = k)
  nn_dists <- matrix(Inf, nrow = n_cells, ncol = k)
  split_idx <- split(seq_len(nrow(sm)), sm$i)
  for (i_chr in names(split_idx)) {
    i <- as.integer(i_chr)
    rows <- split_idx[[i_chr]]
    take <- seq_len(min(k, length(rows)))
    nn_idx[i, take] <- as.integer(sm$j[rows[take]])
    nn_dists[i, take] <- as.numeric(sm$x[rows[take]])
  }
  list(
    nn_idx = nn_idx,
    nn_dists = nn_dists,
    k = k,
    source_matrix = normalizePath(dist_path, winslash = "/", mustWork = TRUE),
    manifest_dataset = as.character(knn_manifest$dataset)
  )
}

LoadGt <- function(ds, cells) {
  ref <- data.table::fread(file.path(CONFIG$paths$ref, paste0("ref_", ds, ".csv")), encoding = "UTF-8")
  cid <- if ("cell_id" %in% names(ref)) "cell_id" else names(ref)[1]
  gt_col <- NULL
  numeric_vals <- NULL
  for (col in c(
    "true_pseudotime", "timepoint", "day", "age(days)", "stage",
    "time_point", "time_info", "plate"
  )) {
    if (!col %in% names(ref)) next
    vals <- ref[[col]]
    vals_chr <- as.character(vals)
    numeric_tokens <- is.na(vals_chr) | grepl("^[-+]?(?:[0-9]+\\.?[0-9]*|\\.[0-9]+)$", vals_chr)
    nv <- if (all(numeric_tokens)) as.numeric(vals_chr) else .extract_ordinal_time(vals_chr)
    if (sum(!is.na(nv)) >= 10 && length(unique(nv[!is.na(nv)])) >= 2) {
      gt_col <- col
      numeric_vals <- nv
      break
    }
  }
  if (is.null(gt_col)) stop("no GT col for ", ds)
  gt <- data.table(
    cell = as.character(ref[[cid]]),
    gt = as.numeric(numeric_vals),
    gt_source_column = gt_col
  )
  gt <- gt[cell %in% cells & is.finite(gt)]
  if (nrow(gt) < 10) stop("too few GT cells for ", ds)
  gt
}

SpearmanSafe <- function(a, b) {
  ok <- is.finite(a) & is.finite(b)
  if (sum(ok) < 10) return(NA_real_)
  as.numeric(cor(a[ok], b[ok], method = "spearman"))
}

rows <- list()
for (ds in CONFIG$datasets) {
  message("=== ", ds)
  rds_path <- file.path(CONFIG$paths$workspace, ds, "linealinkr_complete.rds")
  rds <- readRDS(rds_path)
  extract_cells <- readLines(
    file.path(CONFIG$paths$workspace, ds, "extract", paste0(ds, "_cells.txt")),
    encoding = "UTF-8",
    warn = FALSE
  )
  stopifnot(
    identical(as.character(rds$dataset), ds),
    identical(as.character(rds$package_version), PROJECT_CONFIG$llr$package_version),
    identical(as.character(rds$benchmark_protocol), PROJECT_CONFIG$llr$benchmark_protocol),
    identical(as.character(rds$result_revision), PROJECT_CONFIG$llr$result_revision),
    identical(rds$input_manifest$content_schema, "decoded_h5ad_content_v1"),
    identical(
      as.character(rds$input_manifest$canonical_path),
      as.character(PROJECT_CONFIG$dataset_registry[[ds]])
    ),
    identical(as.character(rds$input_manifest$cells), as.character(rds$cells)),
    identical(as.character(rds$cells), as.character(extract_cells)),
    length(rds$pseudotime) == length(rds$cells),
    all(is.finite(rds$pseudotime)),
    nrow(rds$lineage_weights) == length(rds$cells),
    ncol(rds$lineage_weights) == as.integer(rds$n_lineages),
    all(is.finite(rds$lineage_weights))
  )
  cells <- as.character(rds$cells)
  n_cells <- length(cells)
  blind_pt <- as.numeric(rds$pseudotime)
  names(blind_pt) <- cells
  gt <- LoadGt(ds, cells)

  t0 <- min(gt$gt, na.rm = TRUE)
  early <- gt[gt <= t0 + 1e-8]
  early_cells <- intersect(early$cell, cells)
  # Stable user-root: earliest GT cells, pick median index in the LLR cell order
  early_idx <- sort(match(early_cells, cells))
  user_root_idx <- as.integer(early_idx[ceiling(length(early_idx) / 2)])
  user_root_cell <- cells[user_root_idx]

  nn <- RebuildNnFromSparse(ds, n_cells, rds$input_manifest)
  user_pt <- LineaLinkR::ComputePseudotimeFuzzy(
    nn$nn_idx, nn$nn_dists, start_cell = as.integer(user_root_idx), rank_normalize = TRUE
  )
  names(user_pt) <- cells

  common <- intersect(gt$cell, cells)
  gt_v <- gt[match(common, cell)]$gt
  blind_v <- blind_pt[common]
  user_v <- user_pt[common]
  rho_blind <- SpearmanSafe(blind_v, gt_v)
  rho_user <- SpearmanSafe(user_v, gt_v)

  rows[[length(rows) + 1L]] <- data.table(
    dataset = ds,
    n_cells = n_cells,
    n_gt = length(common),
    gt_source_column = unique(gt$gt_source_column),
    gt_min = t0,
    n_early = length(early_cells),
    user_root_cell = user_root_cell,
    blind_start_cell = as.integer(rds$start_cell),
    package_version = as.character(utils::packageVersion("LineaLinkR")),
    benchmark_protocol = as.character(rds$benchmark_protocol),
    result_revision = as.character(rds$result_revision),
    rds_source = normalizePath(rds_path, winslash = "/", mustWork = TRUE),
    knn_source = nn$source_matrix,
    knn_neighbors = nn$k,
    knn_manifest_dataset = nn$manifest_dataset,
    Spearman_blind = rho_blind,
    Spearman_user_root = rho_user,
    Spearman_user_root_abs = abs(rho_user),
    delta_user_minus_blind = rho_user - rho_blind
  )
  fwrite(
    data.table(cell = cells, pseudotime_blind = blind_pt, pseudotime_user_root = as.numeric(user_pt)),
    file.path(CONFIG$paths$out, paste0(ds, "_blind_vs_user_pseudotime.csv"))
  )
  message(sprintf("  blind=%.3f user=%.3f delta=%.3f", rho_blind, rho_user, rho_user - rho_blind))
  rm(rds, nn)
  gc(FALSE)
}

out <- rbindlist(rows)
summary_path <- file.path(CONFIG$paths$out, "user_root_vs_blind_spearman.csv")
if (selective_run) {
  previous <- fread(summary_path, encoding = "UTF-8")
  retained <- previous[!dataset %in% CONFIG$datasets]
  if (nrow(retained) > 0L) {
    required_retained_columns <- c(
      "dataset", "n_cells", "blind_start_cell", "package_version",
      "benchmark_protocol", "result_revision", "rds_source", "knn_source"
    )
    stopifnot(
      all(required_retained_columns %in% names(retained)),
      !anyDuplicated(retained$dataset),
      all(retained$dataset %in% CONFIG$contract_datasets),
      all(retained$package_version == PROJECT_CONFIG$llr$package_version),
      all(retained$benchmark_protocol == PROJECT_CONFIG$llr$benchmark_protocol),
      all(retained$result_revision == PROJECT_CONFIG$llr$result_revision),
      all(retained$rds_source == file.path(
        CONFIG$paths$workspace,
        retained$dataset,
        "linealinkr_complete.rds"
      )),
      all(retained$knn_source == file.path(
        CONFIG$paths$workspace,
        retained$dataset,
        "knn_dist.mtx"
      ))
    )
    out <- rbindlist(list(retained, out), use.names = TRUE)
  }
}
out[, dataset_order := match(dataset, CONFIG$contract_datasets)]
setorder(out, dataset_order)
out[, dataset_order := NULL]
stopifnot(
  nrow(out) == length(CONFIG$contract_datasets),
  !anyDuplicated(out$dataset),
  identical(out$dataset, CONFIG$contract_datasets)
)
fwrite(out, summary_path)
print(out)
message("mean blind=", mean(out$Spearman_blind), " mean user=", mean(out$Spearman_user_root))
