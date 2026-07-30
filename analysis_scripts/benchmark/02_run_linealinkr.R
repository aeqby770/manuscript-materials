#!/usr/bin/env Rscript
# =============================================================================
# 02_run_linealinkr.R
# LineaLinkR benchmark runner, separate from the competitor entry point.
# =============================================================================
# Usage:
#   Rscript 02_run_linealinkr.R                         # full panel
#   Rscript 02_run_linealinkr.R --datasets pancreas     # one dataset
#   Rscript 02_run_linealinkr.R --entropy-only          # entropy and forks only
#   Rscript 02_run_linealinkr.R --technical-smoke --datasets GSE243917_RA --max-cells=15000
# =============================================================================

suppressPackageStartupMessages({
  library(Matrix)

# path bootstrap (relative; set LLR_ASSETS or LLR_PIPELINE to override)
.llr_cmd <- commandArgs(trailingOnly = FALSE)
.llr_file <- sub("^--file=", "", grep("^--file=", .llr_cmd, value = TRUE)[1])
if (!length(.llr_file) || is.na(.llr_file) || !nzchar(.llr_file)) {
  stop("Run this script via Rscript so --file= is set")
}
.llr_here <- dirname(normalizePath(.llr_file, winslash = "/", mustWork = TRUE))
# Resolve the isolated 0.6.4 checkout for direct entry points.
.llr_work_root <- normalizePath(file.path(.llr_here, "../.."), winslash = "/", mustWork = TRUE)
Sys.setenv(
  LLR064_WORK_ROOT = .llr_work_root,
  LLR_PIPELINE = file.path(.llr_work_root, "pipeline"),
  LLR_ASSETS = .llr_work_root,
  LLR064_RLIB = file.path(.llr_work_root, "rlib")
)
source(file.path(.llr_here, "../config/llr_bootstrap.R"), local = FALSE)

  library(data.table)
})

llr_source_config()

llr_package_library <- Sys.getenv(
  "LLR064_RLIB",
  unset = file.path(CONFIG$dir$root, "rlib")
)
.libPaths(c(llr_package_library, .libPaths()))

`%||%` <- function(x, y) if (is.null(x)) y else x

IsValidClusterCol <- function(meta, cn) {
  if (!cn %in% names(meta)) return(FALSE)
  vals <- as.character(meta[[cn]])
  vals[is.na(vals) | trimws(vals) == ""] <- NA_character_
  vals <- vals[!is.na(vals)]
  nu <- uniqueN(vals)
  if (nu < 3L) return(FALSE)
  n <- length(vals)
  id_like <- grepl("patient|donor|sample_id|sample_name|status|diagnosis|orig\\.ident", cn, ignore.case = TRUE)
  if (id_like && !grepl("cluster", cn, ignore.case = TRUE)) return(FALSE)
  if (nu > max(80L, ceiling(0.2 * n)) && !grepl("cluster", cn, ignore.case = TRUE)) return(FALSE)
  TRUE
}

FindBiologicalSampleCol <- function(meta) {
  candidates <- grep(
    "patient|donor|sample_id|sample_name|orig\\.ident|replicate",
    names(meta), ignore.case = TRUE, value = TRUE
  )
  candidates <- candidates[!grepl(
    "status|condition|diagnosis|cluster|cell|barcode|state|type",
    candidates, ignore.case = TRUE
  )]
  valid <- vapply(candidates, function(cn) {
    values <- as.character(meta[[cn]])
    complete <- !is.na(values) & nzchar(trimws(values))
    counts <- table(values[complete])
    n_levels <- length(counts)
    n_levels >= 2L && n_levels <= min(100L, max(2L, floor(sum(complete) / 10L))) &&
      all(counts >= 10L)
  }, logical(1L))
  if (any(valid)) candidates[[which(valid)[[1L]]]] else NA_character_
}

# Load LineaLinkR.
library(LineaLinkR)

# Command-line arguments.
args <- commandArgs(trailingOnly = TRUE)
filter_datasets <- NULL
entropy_only    <- FALSE
no_temporal     <- FALSE
no_entropy      <- FALSE
max_cells       <- NULL
technical_smoke <- FALSE

split_ds <- function(x) {
  x <- trimws(x)
  x <- x[nzchar(x)]
  unique(x)
}
for (a in args) {
  if (grepl("^--datasets-file=", a)) {
    f <- sub("^--datasets-file=", "", a)
    if (!file.exists(f)) stop("datasets file not found: ", f)
    filter_datasets <- split_ds(readLines(f, warn = FALSE))
  }
  if (grepl("^--datasets=", a)) {
    raw <- sub("^--datasets=", "", a)
    sep <- if (grepl(";", raw, fixed = TRUE)) ";" else ","
    filter_datasets <- split_ds(unlist(strsplit(raw, sep, fixed = TRUE)))
  }
  if (grepl("^--max-cells=", a)) {
    max_cells <- as.integer(sub("^--max-cells=", "", a))
  }
  if (a == "--entropy-only") entropy_only <- TRUE
  if (a == "--no-temporal")  no_temporal  <- TRUE
  if (a == "--no-entropy")   no_entropy   <- TRUE
  if (a == "--technical-smoke") technical_smoke <- TRUE
}
if (!is.null(max_cells) && !technical_smoke) {
  stop("--max-cells is restricted to --technical-smoke; formal results require all cells")
}
if (technical_smoke && is.null(max_cells)) max_cells <- 2000L

# Initialization.
OUT <- if (technical_smoke) {
  file.path(CONFIG$dir$cache, "smoke")
} else {
  CONFIG$dir$linealinkr
}
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

LlrDsDir <- function(ds_name) {
  d <- file.path(OUT, ds_name)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}
LlrExtractDir <- function(ds_name) {
  d <- file.path(OUT, ds_name, "extract")
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

set.seed(CONFIG$seed)
Sys.setenv(OMP_NUM_THREADS = CONFIG$n_threads)
data.table::setDTthreads(CONFIG$n_threads)

cat(sprintf("OMP threads: %s | seed: %d\n",
    Sys.getenv("OMP_NUM_THREADS"), CONFIG$seed))

# =============================================================================
# Helper functions (trusted RDS contract: protocol, revision, and core config).
# =============================================================================
LlrExpectedMeta <- function() {
  cfg <- CONFIG$llr %||% list()
  list(
    benchmark_protocol = as.character(cfg$benchmark_protocol),
    result_revision = as.character(cfg$result_revision),
    root_strategy = cfg$root_strategy %||% "attractor_energy",
    root_potency_mode = cfg$root_potency_mode %||% "detected_genes_high",
    root_endpoint_alignment_min_gain = as.numeric(
      cfg$root_endpoint_alignment_min_gain %||% 0.20
    ),
    root_endpoint_require_terminal_leaf = isTRUE(
      cfg$root_endpoint_require_terminal_leaf
    ),
    root_ordering_min_gain = as.numeric(
      cfg$root_ordering_min_gain %||% 0.10
    ),
    root_representation_alignment_min_gain = as.numeric(
      cfg$root_representation_alignment_min_gain %||% 0.00
    ),
    run_no_temporal = isTRUE(no_temporal),
    run_no_entropy = isTRUE(no_entropy),
    run_entropy_only = isTRUE(entropy_only),
    root_direction_min_potency_cor = as.numeric(
      cfg$root_direction_min_potency_cor %||% 0.10
    ),
    root_ordering_fallback_threshold = as.numeric(
      cfg$root_ordering_fallback_threshold %||% 0.45
    ),
    root_fallback_k_max = as.integer(cfg$root_fallback_k_max %||% 25L),
    root_fallback_k_min = as.integer(cfg$root_fallback_k_min %||% 8L),
    root_fallback_cells_per_cluster = as.integer(
      cfg$root_fallback_cells_per_cluster %||% 50L
    ),
    root_fallback_nstart = as.integer(cfg$root_fallback_nstart %||% 10L),
    root_fallback_kmeans_algorithm = cfg$root_fallback_kmeans_algorithm %||%
      "Hartigan-Wong",
    root_fallback_selection_score_mode =
      cfg$root_fallback_selection_score_mode %||% "confidence_regularized",
    root_fallback_weights = as.numeric(
      cfg$root_fallback_weights %||% c(0.30, 0.25, 0.20, 0.15, 0.10)
    ),
    lineage_embedding = cfg$lineage_embedding %||% "pca2",
    lineage_topology_k = as.integer(cfg$lineage_topology_k %||% 15L),
    lineage_topology_nstart = as.integer(cfg$lineage_topology_nstart %||% 10L),
    root_top_fraction = as.numeric(cfg$root_top_fraction %||% 0.05),
    root_n_attractors = as.integer(cfg$root_n_attractors %||% 0L),
    root_max_iter = as.integer(cfg$root_max_iter %||% 300L),
    root_tol = as.numeric(cfg$root_tol %||% 1e-5),
    root_weights = as.numeric(cfg$root_weights %||% c(0.25, 0.20, 0.30, 0.15, 0.10)),
    start_cell_method = cfg$start_cell_method %||% "auto",
    prefer_root_label = isTRUE(cfg$prefer_root_label),
    rank_normalize = isTRUE(cfg$rank_normalize),
    force_kmeans_refine = isTRUE(cfg$force_kmeans_refine),
    knn_neighbors = as.integer(cfg$knn_neighbors %||% 30L),
    lr_scale_factor = as.numeric(cfg$lr_scale_factor %||% 10000),
    augment_mst_bridges = isTRUE(cfg$augment_mst_bridges),
    connect_components = isTRUE(cfg$connect_components),
    package_version = as.character(cfg$package_version),
    root_ensemble_k_grid = as.integer(cfg$root_ensemble_k_grid %||% c(30L, 50L, 70L)),
    root_ensemble_pc_grid = as.integer(cfg$root_ensemble_pc_grid %||% c(10L, 20L, 30L)),
    root_ensemble_resamples = as.integer(cfg$root_ensemble_resamples %||% 20L),
    root_ensemble_min_resamples = as.integer(cfg$root_ensemble_min_resamples %||% 3L),
    root_ensemble_adaptive_resampling = isTRUE(cfg$root_ensemble_adaptive_resampling),
    root_ensemble_subsample_fraction = as.numeric(cfg$root_ensemble_subsample_fraction %||% 0.8),
    root_ensemble_direction_margin = as.numeric(cfg$root_ensemble_direction_margin %||% 0.10),
    n_bins = as.integer(cfg$n_bins %||% 12L),
    min_n_eff = as.numeric(cfg$min_n_eff %||% 5.0),
    fork_window_bins = as.integer(cfg$fork_window_bins %||% 2L),
    fork_activity_floor = as.numeric(cfg$fork_activity_floor %||% 1e-4),
    fork_sample_boot = as.integer(cfg$fork_sample_boot %||% 999L)
  )
}

LlrMetaSame <- function(found, expected) {
  if (is.null(found)) return(FALSE)
  if (is.numeric(expected)) {
    return(isTRUE(all.equal(as.numeric(found), expected, tolerance = 1e-12, check.attributes = FALSE)))
  }
  identical(found, expected)
}

LlrCanonicalPath <- function(ds_name) {
  DatasetPath(ds_name)
}

LlrDecodedInputManifest <- function(ds_name, h5, expr, meta, genes, cells, matrix_source) {
  row_anchor_index <- unique(as.integer(round(seq(
    1L, nrow(expr), length.out = min(17L, nrow(expr))
  ))))
  col_anchor_index <- unique(as.integer(round(seq(
    1L, ncol(expr), length.out = min(17L, ncol(expr))
  ))))
  matrix_values <- as.numeric(expr@x)
  finite_values <- matrix_values[is.finite(matrix_values)]
  meta_anchor_index <- col_anchor_index
  meta_anchor_values <- as.data.frame(lapply(
    meta[meta_anchor_index],
    as.character
  ), stringsAsFactors = FALSE)
  list(
    content_schema = "decoded_h5ad_content_v1",
    dataset = ds_name,
    canonical_path = normalizePath(h5, winslash = "/", mustWork = TRUE),
    matrix_source = as.character(matrix_source),
    n_cells = as.integer(ncol(expr)),
    n_genes = as.integer(nrow(expr)),
    cells = as.character(cells),
    genes = as.character(genes),
    matrix_nnz = as.numeric(Matrix::nnzero(expr)),
    matrix_stored_values = as.numeric(length(matrix_values)),
    matrix_finite_values = as.numeric(length(finite_values)),
    matrix_min = as.numeric(min(c(0, finite_values))),
    matrix_max = as.numeric(max(c(0, finite_values))),
    matrix_sum = as.numeric(sum(finite_values)),
    matrix_sum_squares = as.numeric(sum(finite_values ^ 2)),
    row_anchor_index = row_anchor_index,
    row_anchor_genes = as.character(genes[row_anchor_index]),
    row_anchor_sums = as.numeric(Matrix::rowSums(expr[row_anchor_index, , drop = FALSE])),
    col_anchor_index = col_anchor_index,
    col_anchor_cells = as.character(cells[col_anchor_index]),
    col_anchor_sums = as.numeric(Matrix::colSums(expr[, col_anchor_index, drop = FALSE])),
    metadata_columns = names(meta),
    metadata_classes = vapply(
      meta,
      function(x) paste(class(x), collapse = "/"),
      character(1L)
    ),
    metadata_anchor_cells = as.character(cells[meta_anchor_index]),
    metadata_anchor_values = meta_anchor_values
  )
}

LlrInputManifestTrusted <- function(found, expected) {
  if (is.null(found) || !is.list(found)) return(FALSE)
  all(vapply(names(expected), function(nm) LlrMetaSame(found[[nm]], expected[[nm]]), logical(1L)))
}

LlrRdsComplete <- function(x) {
  n <- as.integer(x$n_cells %||% length(x$pseudotime))
  isTRUE(length(n) == 1L && is.finite(n) && n > 0L) &&
    is.character(x$dataset) && length(x$dataset) == 1L && nzchar(x$dataset) &&
    is.character(x$cells) && length(x$cells) == n && !anyNA(x$cells) && uniqueN(x$cells) == n &&
    is.numeric(x$pseudotime) && length(x$pseudotime) == n && any(is.finite(x$pseudotime)) &&
    is.matrix(x$lineage_weights) && nrow(x$lineage_weights) == n && ncol(x$lineage_weights) >= 1L &&
    all(is.finite(x$lineage_weights)) &&
    (is.null(x$lr_scores) || is.data.frame(x$lr_scores)) &&
    (is.null(x$temporal_lr) || (is.list(x$temporal_lr) && is.data.frame(x$temporal_lr$scores) && is.list(x$temporal_lr$intermediate))) &&
    identical(as.character(x$package_version), "0.6.4") &&
    is.numeric(x$root_support) && length(x$root_support) == 1L &&
    is.character(x$direction_status) && length(x$direction_status) == 1L &&
    "root_alternative" %in% names(x) &&
    "pseudotime_alternative" %in% names(x) &&
    !is.null(x$fork_contrasts) && !is.null(x$fork_evidence) &&
    is.integer(x$terminal_anchor_cluster) && length(x$terminal_anchor_cluster) == 1L &&
    (is.na(x$terminal_anchor_cluster) || x$terminal_anchor_cluster >= 1L) &&
    is.logical(x$terminal_leaf_eligible) && length(x$terminal_leaf_eligible) == 1L &&
    is.integer(x$terminal_leaf_clusters) && all(x$terminal_leaf_clusters >= 1L) &&
    is.integer(x$selector_topology_k) && length(x$selector_topology_k) == 1L &&
    x$selector_topology_k >= 2L &&
    is.data.frame(x$root_selection) && nrow(x$root_selection) == 4L &&
    all(c(
      "terminal_anchor_is_high_end",
      "terminal_anchor_is_leaf",
      "terminal_leaf_eligible"
    ) %in% names(x$root_selection)) &&
    identical(
      as.logical(x$root_selection$terminal_leaf_eligible),
      as.logical(x$root_selection$terminal_anchor_is_high_end) &
        as.logical(x$root_selection$terminal_anchor_is_leaf)
    ) &&
    is.character(x$root_cluster_source) && length(x$root_cluster_source) == 1L && nzchar(x$root_cluster_source) &&
    is.integer(x$root_clusters) && length(x$root_clusters) == n && all(x$root_clusters >= 1L) &&
    is.integer(x$root_cluster_empty_start_warnings) && length(x$root_cluster_empty_start_warnings) == 1L && x$root_cluster_empty_start_warnings >= 0L &&
    is.character(x$lineage_cluster_source) && length(x$lineage_cluster_source) == 1L && nzchar(x$lineage_cluster_source) &&
    is.integer(x$lineage_clusters) && length(x$lineage_clusters) == n && all(x$lineage_clusters >= 1L) &&
    is.integer(x$lineage_topology_k_requested) && length(x$lineage_topology_k_requested) == 1L && x$lineage_topology_k_requested >= 2L &&
    is.integer(x$lineage_topology_k_occupied) && length(x$lineage_topology_k_occupied) == 1L && x$lineage_topology_k_occupied >= 1L &&
    !is.null(x$config_snapshot$algorithm) &&
    !is.null(x$input_manifest)
}

LlrRdsIsTrusted <- function(x, expected_manifest, expected_cells) {
  if (is.null(x) || !is.list(x)) return(FALSE)
  exp <- LlrExpectedMeta()
  saved <- x$config_snapshot %||% list()
  force_auto_clusters <- expected_manifest$dataset %chin% BENCHMARK_FORCE_AUTO_CLUSTER_DATASETS
  expected_cluster_policy <- c(
    "metadata_or_expression_auto_clusters",
    "expression_auto_clusters"
  )[1L + as.integer(force_auto_clusters)]
  cluster_input_policy_ok <- !force_auto_clusters || (
    identical(
      as.character(x$cluster_input_policy %||% saved$cluster_input_policy),
      expected_cluster_policy
    ) &&
      identical(as.character(x$celltype_column), "auto_clusters")
  )
  saved$benchmark_protocol <- x$benchmark_protocol %||% saved$benchmark_protocol
  saved$result_revision <- x$result_revision %||% saved$result_revision
  saved$root_potency_mode <- x$root_potency_mode %||% saved$root_potency_mode
  saved$lineage_embedding <- x$lineage_embedding %||% saved$lineage_embedding
  cluster_input_policy_ok &&
    all(vapply(names(exp), function(nm) LlrMetaSame(saved[[nm]], exp[[nm]]), logical(1L))) &&
    LlrRdsComplete(x) &&
    identical(as.character(x$cells), as.character(expected_cells)) &&
    LlrInputManifestTrusted(x$input_manifest, expected_manifest)
}

LlrTrustCache <- function(cache_file, expected_manifest, expected_cells) {
  if (!file.exists(cache_file)) return(FALSE)
  x <- readRDS(cache_file)
  if (!LlrRdsIsTrusted(x, expected_manifest, expected_cells)) {
    cat(sprintf("  STALE CACHE (trust/manifest): %s\n", cache_file))
    return(FALSE)
  }
  TRUE
}

# Read a canonical H5AD and convert it to a LineaLinkR dgCMatrix.
ReadH5ADforLLR <- function(ds_name) {
  h5 <- LlrCanonicalPath(ds_name)
  if (!file.exists(h5)) stop(sprintf("Canonical file not found: %s", h5))

  # Use Python to extract counts + metadata
  py <- sprintf(r"(
import scanpy as sc, numpy as np, pandas as pd, os, gzip, json

inpath = r'%s'
outdir = r'%s'
dsname = '%s'

adata = sc.read(inpath)

# Expression matrix: use counts layer if available, else X; keep sparse.
if 'counts' in adata.layers:
    mat = adata.layers['counts']
    matrix_source = 'layers/counts'
else:
    mat = adata.X
    matrix_source = 'X'
from scipy import sparse
if sparse.issparse(mat):
    mat = mat.T.tocsr()
else:
    mat = sparse.csr_matrix(np.asarray(mat).T)

# Metadata
meta = adata.obs.copy()
for col in ['clusters', 'cluster']:
    if col in meta.columns and (('cell_type' not in meta.columns) or meta['cell_type'].isna().all()):
        meta['cell_type'] = meta[col].astype(str)
for col in ['cell_type', 'clusters']:
    if col in meta.columns and (('cluster' not in meta.columns) or meta['cluster'].isna().all()):
        meta['cluster'] = meta[col].astype(str)
if 'pseudotime' not in meta.columns:
    meta['pseudotime'] = np.nan

# Save as MatrixMarket + metadata CSV
from scipy.io import mmwrite
mmwrite(os.path.join(outdir, f'{dsname}_matrix.mtx'), mat)

meta.to_csv(
    os.path.join(outdir, f'{dsname}_meta.csv'),
    index=True,
    index_label='barcode_id'
)
np.savetxt(os.path.join(outdir, f'{dsname}_genes.txt'), adata.var_names, fmt='%%s')
np.savetxt(os.path.join(outdir, f'{dsname}_cells.txt'), adata.obs_names, fmt='%%s')
with open(os.path.join(outdir, f'{dsname}_matrix_source.txt'), 'w', encoding='utf-8', newline='\n') as stream:
    stream.write(matrix_source + '\n')

print('OK')
)", h5, LlrExtractDir(ds_name), ds_name)

  ext <- LlrExtractDir(ds_name)
  mtx_path <- file.path(ext, sprintf("%s_matrix.mtx", ds_name))
  meta_path <- file.path(ext, sprintf("%s_meta.csv", ds_name))
  genes_path <- file.path(ext, sprintf("%s_genes.txt", ds_name))
  cells_path <- file.path(ext, sprintf("%s_cells.txt", ds_name))
  matrix_source_path <- file.path(ext, sprintf("%s_matrix_source.txt", ds_name))
  manifest_path <- file.path(ext, sprintf("%s_input_content_manifest.rds", ds_name))
  unlink(c(mtx_path, meta_path, genes_path, cells_path, matrix_source_path, manifest_path), force = TRUE)
  extract_status <- sys::exec_wait(CONFIG$exe$python, args = c("-X", "utf8", "-c", py),
    std_out = file.path(CONFIG$dir$log, sprintf("llr_extract_%s.log", ds_name)),
    std_err = file.path(CONFIG$dir$log, sprintf("llr_extract_%s.err", ds_name)))
  if (extract_status != 0L) stop(sprintf("Python h5ad extraction failed for %s", ds_name))

  if (!file.exists(mtx_path)) stop("Matrix export failed")

  expr <- as(readMM(mtx_path), "CsparseMatrix")  # genes x cells, CSC required by Rust bindings
  meta <- fread(meta_path, header = TRUE)
  genes <- readLines(genes_path, encoding = "UTF-8", warn = FALSE)
  cells <- readLines(cells_path, encoding = "UTF-8", warn = FALSE)
  matrix_source <- trimws(readLines(matrix_source_path, encoding = "UTF-8", warn = FALSE)[1L])

  if (ncol(expr) != length(cells) || nrow(meta) != length(cells)) {
    stop(sprintf(
      "Extracted cell alignment mismatch for %s: matrix=%d metadata=%d cells=%d",
      ds_name, ncol(expr), nrow(meta), length(cells)
    ))
  }
  meta_cell_col <- if ("barcode_id" %in% names(meta)) "barcode_id" else names(meta)[1]
  if (!identical(as.character(meta[[meta_cell_col]]), as.character(cells))) {
    stop(sprintf("Extracted metadata cell order mismatch for %s", ds_name))
  }
  if (nrow(expr) != length(genes)) {
    stop(sprintf(
      "Extracted gene alignment mismatch for %s: matrix=%d genes=%d",
      ds_name, nrow(expr), length(genes)
    ))
  }

  colnames(expr) <- cells
  rownames(expr) <- genes
  input_manifest <- LlrDecodedInputManifest(
    ds_name, h5, expr, meta, genes, cells, matrix_source
  )
  saveRDS(input_manifest, manifest_path)

  list(expr = expr, meta = meta, genes = genes, cells = cells, input_manifest = input_manifest)
}

# Build the LR database (gene symbols to row indices).
BuildLRDB <- function(expr, species = "human") {
  if (!species %in% c("human", "mouse")) {
    return(list(ligand_rows = list(), receptor_rows = list()))
  }
  cc_db <- if (species == "mouse") CellChat::CellChatDB.mouse else CellChat::CellChatDB.human
  lr_pairs <- cc_db$interaction
  ligand_rows <- list()
  receptor_rows <- list()
  gene_names <- rownames(expr)
  for (i in seq_len(nrow(lr_pairs))) {
    lig <- strsplit(lr_pairs$ligand[i], "_")[[1]]
    rec <- strsplit(lr_pairs$receptor[i], "_")[[1]]
    lig_idx <- match(lig, gene_names)
    rec_idx <- match(rec, gene_names)
    if (any(is.na(lig_idx)) || any(is.na(rec_idx))) next
    ligand_rows[[length(ligand_rows) + 1]] <- lig_idx
    receptor_rows[[length(receptor_rows) + 1]] <- rec_idx
  }
  list(ligand_rows = ligand_rows, receptor_rows = receptor_rows)
}

# =============================================================================
# Run the complete LineaLinkR workflow for one dataset.
# =============================================================================
RunLineaLinkROnDataset <- function(ds_name) {
  cat(sprintf("\n===== LineaLinkR on %s =====\n", ds_name))
  result_dir <- LlrDsDir(ds_name)
  log_file <- file.path(CONFIG$dir$log, sprintf("llr_%s.log", ds_name))

  # Decode the canonical input before checking reusable results.
  cat("  Step 1/7: Reading data ... ")
  dat <- ReadH5ADforLLR(ds_name)
  expr <- dat$expr
  meta <- dat$meta
  cells <- colnames(expr)
  if (!is.null(max_cells) && ncol(expr) > max_cells) {
    cell_idx <- sample.int(ncol(expr), max_cells)
    expr <- expr[, cell_idx, drop = FALSE]
    meta <- meta[cell_idx]
    cells <- colnames(expr)
  }
  cat(sprintf("done [%d genes x %d cells]\n", nrow(expr), ncol(expr)))

  # Reuse requires the decoded input contract and exact result-cell order.
  cache_file <- file.path(result_dir, "linealinkr_complete.rds")
  if (LlrTrustCache(cache_file, dat$input_manifest, cells)) {
    cat(sprintf("  CACHED (trusted): %s\n", cache_file))
    return(readRDS(cache_file))
  }
  stale_archive_dir <- file.path(
    CONFIG$dir$archive,
    "legacy_linealinkr_cache",
    ds_name
  )
  dir.create(stale_archive_dir, recursive = TRUE, showWarnings = FALSE)
  run_attempt <- 1L
  if (file.exists(cache_file)) {
    stale_file <- file.path(stale_archive_dir, "linealinkr_complete.rds")
    if (file.exists(stale_file)) {
      stop(sprintf("Fixed legacy archive already exists: %s", stale_file))
    }
    if (!file.rename(cache_file, stale_file)) {
      stop(sprintf("Could not archive stale cache: %s", cache_file))
    }
  }

  ext <- LlrExtractDir(ds_name)
  knn_matrix_path <- file.path(ext, sprintf("%s_matrix.mtx", ds_name))
  if (!is.null(max_cells) && ncol(expr) < ncol(dat$expr)) {
    knn_matrix_path <- file.path(result_dir, sprintf("%s_matrix_downsampled.mtx", ds_name))
    writeMM(expr, knn_matrix_path)
  }

  # Step 2: prepare cell types.
  cat("  Step 2/7: Cell types ... ")
  force_auto_clusters <- ds_name %chin% BENCHMARK_FORCE_AUTO_CLUSTER_DATASETS
  cluster_input_policy <- c(
    "metadata_or_expression_auto_clusters",
    "expression_auto_clusters"
  )[1L + as.integer(force_auto_clusters)]
  # Prefer explicit cluster/cell_type columns; then any "*_clusters" pattern;
  # never fall back to names(meta)[1] (that is the row-name column "V1").
  ct_candidates <- c("clusters", "cluster", "seurat_clusters", "clusters_fine", "clusters_coarse",
                     "paul15_clusters", "louvain", "leiden",
                     "cell_type", "celltype", "ct", "ct_subtype", "state_info", "cell_state",
                     "celltype_label", "cell.label", "lineage", "subtype")
  ct_candidates <- ct_candidates[!force_auto_clusters]
  ct_col <- NA_character_
  for (cn in ct_candidates) {
    if (IsValidClusterCol(meta, cn)) { ct_col <- cn; break }
  }
  if (is.na(ct_col)) {
    cluster_cols <- grep("cluster", names(meta), ignore.case = TRUE, value = TRUE)
    cluster_cols <- cluster_cols[!force_auto_clusters]
    for (cn in cluster_cols) {
      if (IsValidClusterCol(meta, cn)) { ct_col <- cn; break }
    }
  }
  if (is.na(ct_col)) {
    cl_csv <- file.path(result_dir, sprintf("%s_auto_clusters.csv", ds_name))
    cluster_manifest_path <- file.path(result_dir, "auto_cluster_manifest.rds")
    expected_cluster_manifest <- c(
      dat$input_manifest,
      list(
        cluster_revision = "linealinkr064_scanpy_leiden_pca30_neighbors30",
        cluster_input_policy = cluster_input_policy,
        seed = CONFIG$seed
      )
    )
    cluster_cache_ok <- FALSE
    cat(if (cluster_cache_ok) {
      "  no cell-type column, using cached unsupervised clusters ... "
    } else {
      "  no cell-type column, running leiden clustering ... "
    })
    cluster_py <- sprintf(r"(
import scanpy as sc, numpy as np, pandas as pd
from scipy.io import mmread
from sklearn.cluster import KMeans
a = sc.AnnData(mmread(r'%s').tocsr().T)
sc.pp.normalize_total(a)
sc.pp.log1p(a)
sc.pp.highly_variable_genes(a, n_top_genes=2000)
a = a[:, a.var.highly_variable].copy()
sc.pp.pca(a, n_comps=30)
sc.pp.neighbors(a, n_neighbors=30)
# Adaptive resolution: try multiple values, pick the first with >=3 clusters
best_labels = None
for res in [0.1, 0.2, 0.3, 0.5, 0.8, 1.0, 1.5]:
    sc.tl.leiden(a, resolution=res, key_added='auto_clusters', random_state=%d, flavor='igraph', n_iterations=2, directed=False)
    n = a.obs['auto_clusters'].nunique()
    print(f'resolution={res}: {n} clusters')
    if n >= 3:
        best_labels = a.obs['auto_clusters'].astype(str).to_numpy()
        break
if best_labels is None:
    # Leiden failed: use kmeans
    km = KMeans(n_clusters=8, random_state=%d, n_init=10).fit(a.obsm['X_pca'])
    best_labels = np.array([f'C{{c}}' for c in km.labels_])
    print('Used kmeans k=8 (leiden failed)')
pd.DataFrame({'auto_clusters': best_labels}, index=a.obs_names.astype(str)).to_csv(r'%s')
print('CLUSTERS_OK')
)", knn_matrix_path, CONFIG$seed, CONFIG$seed, file.path(result_dir, sprintf("%s_auto_clusters.csv", ds_name)))
    if (!cluster_cache_ok) {
      cl_status <- sys::exec_wait(CONFIG$exe$python, args = c("-X", "utf8", "-c", cluster_py),
        std_out = log_file,
        std_err = file.path(CONFIG$dir$log, sprintf("llr_%s_cluster.err", ds_name)))
      if (cl_status != 0L) {
        cat(sprintf("FAIL (leiden clustering failed for %s)\n", ds_name))
        return(invisible(NULL))
      }
      saveRDS(expected_cluster_manifest, cluster_manifest_path)
    }
    if (file.exists(cl_csv)) {
      cl_dt <- fread(cl_csv)
      # CSV row names may be numeric indices, not barcodes.
      # If row count matches, assign by position (both follow scanpy obs order).
      if (nrow(cl_dt) == length(cells)) {
        meta$auto_clusters <- cl_dt$auto_clusters
      } else {
        meta$auto_clusters <- cl_dt$auto_clusters[match(cells, cl_dt[[1]])]
      }
      ct_col <- "auto_clusters"
      cat(sprintf("done [%d clusters]\n", length(unique(meta[[ct_col]]))))
    } else {
      cat(sprintf("FAIL (no cluster CSV for %s)\n", ds_name))
      return(invisible(NULL))
    }
  }
  ct_labels <- as.character(meta[[ct_col]])
  ct_labels[is.na(ct_labels) | ct_labels == ""] <- "unassigned"
  celltypes <- as.integer(as.factor(ct_labels))
  n_types <- max(celltypes, na.rm = TRUE)
  n_cells <- ncol(expr)
  if(is.na(n_types) || is.infinite(n_types) || n_types < 1) {
    stop(sprintf("No valid cell-type labels for %s", ds_name))
  }
  cat(sprintf("%d types (from '%s')\n", n_types, ct_col))

  # For lineage extraction, use finer clusters if available
  fine_col <- if (force_auto_clusters) ct_col
              else if ("clusters_fine" %in% names(meta)) "clusters_fine"
              else if ("clusters_coarse" %in% names(meta)) "clusters_coarse"
              else ct_col
  lineage_labels <- as.character(meta[[fine_col]])
  lineage_labels[is.na(lineage_labels) | lineage_labels == ""] <- "unassigned"
  lineage_cl <- as.integer(as.factor(lineage_labels))
  lineage_cluster_source <- sprintf("metadata:%s", fine_col)
  n_lineage_cl <- max(lineage_cl, na.rm = TRUE)
  # If too few clusters (<8) for MST to capture branching, auto-refine via
  # kmeans on the UMAP embedding (computed later in Step 3).  Flag for deferred
  # refinement — actual kmeans runs after UMAP is available.
  need_kmeans_refine <- (n_lineage_cl < 8)
  cat(sprintf("  Lineage clusters: %d (from '%s')%s\n",
    n_lineage_cl, fine_col,
    if(need_kmeans_refine) " [will refine via kmeans after UMAP]" else ""))

  llr_cfg <- if (!is.null(CONFIG$llr)) CONFIG$llr else list(
    start_cell_method = "mst_extremum", rank_normalize = FALSE,
    knn_neighbors = 30L, augment_mst_bridges = TRUE
  )
  knn_k <- as.integer(min(llr_cfg$knn_neighbors %||% 30L, max(2L, ncol(expr) - 1L)))

  # Step 3: compute pseudotime.
  cat("  Step 3/7: Pseudotime (fuzzy Dijkstra) ... ")
  # Scanpy-equivalent normalization, Seurat-HVG, PCA, and neighbour backends
  # without importing the full Scanpy stack. UMAP is generated only for figures.
  knn_script <- file.path(
    CONFIG$dir$manuscript,
    "taskM_knn_precompute_one.py"
  )
  knn_manifest_path <- file.path(result_dir, "knn_precompute_manifest.rds")
  expected_knn_manifest <- c(
    dat$input_manifest,
    list(
      preprocessor_contract = "taskM_pca30_hvg2000_knn_umap",
      seed = CONFIG$seed,
      knn_neighbors = knn_k,
      write_umap = TRUE
    )
  )
  knn_required <- file.path(
    result_dir,
    c("knn_dist.mtx", "knn_conn.mtx", "pca.csv", "umap.csv")
  )
  unlink(knn_required, force = TRUE)
  knn_status <- sys::exec_wait(
    CONFIG$exe$python,
    args = c(
      "-X", "utf8", knn_script, ds_name, result_dir,
      knn_matrix_path, as.character(knn_k), "with_umap"
    ),
    std_out = log_file,
    std_err = file.path(CONFIG$dir$log, sprintf("llr_%s.err", ds_name)))
  if (knn_status != 0L) {
    stop(sprintf("Python kNN computation failed for %s", ds_name))
  }

  # Decode every output before the preprocessor contract is recorded.
  knn_dist_path <- file.path(result_dir, "knn_dist.mtx")
  knn_conn_path <- file.path(result_dir, "knn_conn.mtx")

  if (!all(file.exists(knn_required))) {
    stop(sprintf("kNN computation did not write required MatrixMarket outputs for %s", ds_name))
  }

  knn_dist <- as(readMM(knn_dist_path), "CsparseMatrix")
  knn_conn <- as(readMM(knn_conn_path), "CsparseMatrix")
  pca_embedding <- as.matrix(fread(file.path(result_dir, "pca.csv"), header = FALSE))
  umap <- as.matrix(fread(file.path(result_dir, "umap.csv"), header = FALSE))
  n_cells <- ncol(expr)
  k <- knn_k
  knn_contract_ok <- identical(dim(knn_dist), c(n_cells, n_cells)) &&
    identical(dim(knn_conn), c(n_cells, n_cells)) &&
    nrow(pca_embedding) == n_cells && ncol(pca_embedding) >= 2L &&
    nrow(umap) == n_cells && ncol(umap) == 2L &&
    all(is.finite(knn_dist@x)) && all(knn_dist@x >= 0) &&
    all(is.finite(knn_conn@x)) && all(knn_conn@x >= 0) &&
    all(Matrix::rowSums(knn_conn) > 0) &&
    all(is.finite(pca_embedding)) && all(is.finite(umap))
  if (!knn_contract_ok) {
    stop(sprintf("Decoded kNN/PCA/UMAP contract failed for %s", ds_name))
  }
  saveRDS(expected_knn_manifest, knn_manifest_path)

  # Reconstruct nearest-neighbour index and distance matrices.
  # Vectorized: summary() gives (i, j, x) triplets for all non-zero entries.
  sm <- summary(knn_dist)
  sm_i <- sm$i; sm_j <- sm$j; sm_x <- sm$x
  ord <- order(sm_i, sm_x, decreasing = c(FALSE, FALSE))
  sm_i <- sm_i[ord]; sm_j <- sm_j[ord]; sm_x <- sm_x[ord]
  run_starts <- c(0L, which(diff(sm_i) != 0L), length(sm_i))
  nn_idx <- matrix(0L, nrow = n_cells, ncol = k)
  nn_dists <- matrix(0.0, nrow = n_cells, ncol = k)
  for (g in seq_len(length(run_starts) - 1L)) {
    i <- sm_i[run_starts[g] + 1L]
    s <- run_starts[g] + 1L
    e <- run_starts[g + 1L]
    take <- seq_len(min(k, e - s + 1L))
    nn_idx[i, take] <- sm_j[s + take - 1L]
    nn_dists[i, take] <- sm_x[s + take - 1L]
  }

  # ── Step 3b: UMAP (needed before blind root + lineage weights) ─────
  is_sim_dataset <- grepl("^(sim_|simulated_)", ds_name)
  blind_emb <- pca_embedding
  lineage_embedding_mode <- llr_cfg$lineage_embedding %||% "pca2"
  geom_emb <- switch(
    lineage_embedding_mode,
    pca2 = blind_emb[, seq_len(min(2L, ncol(blind_emb))), drop = FALSE],
    pca = blind_emb,
    umap = umap,
    stop(sprintf("Unsupported lineage_embedding: %s", lineage_embedding_mode))
  )
  selector_topology_k <- min(
    as.integer(llr_cfg$lineage_topology_k %||% 15L),
    n_cells
  )
  lineage_topology <- LineaLinkR::ClusterLineageTopology(
    geom_emb,
    k = selector_topology_k,
    nstart = as.integer(llr_cfg$lineage_topology_nstart %||% 10L),
    iter_max = 500L,
    seed = CONFIG$seed
  )

  # Make the blind root search see a connected graph; otherwise the score can
  # collapse to -Inf and the root picker degenerates to the first method.
  if (isTRUE(llr_cfg$connect_components)) {
    cc <- LineaLinkR::ConnectKnnComponents(nn_idx, nn_dists, geom_emb)
    nn_idx <- cc$nn_idx
    nn_dists <- cc$nn_dists
    if (!is.null(cc$n_components_before) && cc$n_components_before > 1L) {
      cat(sprintf("  [connected %d comps] ", cc$n_components_before))
    }
  }

  # Optional k-means refinement (off in the canonical submission protocol).
  force_km <- isTRUE(llr_cfg$force_kmeans_refine)
  root_cluster_empty_start_warnings <- 0L
  if (need_kmeans_refine && force_km) {
    k_target <- min(25L, as.integer(n_cells / 50))
    k_target <- max(k_target, 8L)
    root_topology <- LineaLinkR::ClusterLineageTopology(
      geom_emb,
      k = k_target,
      nstart = 10,
      iter_max = 500L,
      seed = CONFIG$seed
    )
    lineage_cl <- root_topology$labels
    n_lineage_cl <- uniqueN(lineage_cl)
    root_cluster_empty_start_warnings <-
      root_topology$empty_start_warnings
    lineage_cluster_source <- sprintf(
      "kmeans:%s:requested_k=%d:occupied_k=%d",
      lineage_embedding_mode,
      k_target,
      n_lineage_cl
    )
    cat(sprintf("[kmeans k=%d/%d occupied] ", n_lineage_cl, k_target))
    if (root_cluster_empty_start_warnings > 0L) {
      cat(sprintf(
        "[audited root empty starts=%d] ",
        root_cluster_empty_start_warnings
      ))
    }
  }

  # Root protocol: never use GT / Palantir pseudotime columns
  start_cfg <- llr_cfg$start_cell_method %||% "auto"
  rank_norm <- isTRUE(llr_cfg$rank_normalize)
  root_anchor_used <- FALSE
  root_strategy <- llr_cfg$root_strategy %||% "attractor_energy"
  root_potency_mode <- llr_cfg$root_potency_mode %||% "detected_genes_high"
  root_diagnostics <- NULL
  root_potency_vector <- NULL
  direction_diagnostics <- NULL
  root_selection <- NULL
  root_ensemble_perturbations <- NULL
  used_confidence_regularized_fallback <- FALSE
  root_adaptive_stopped <- FALSE
  root_resamples_used <- NA_integer_
  root_support <- NA_real_
  direction_margin <- NA_real_
  direction_status <- "unresolved"
  root_evidence_mode <- "blind"
  root_input_source <- "blind_expression_pca_knn_no_gt"
  selected_root_potency_mode <- root_potency_mode
  root_alternative <- NULL
  alternative_start_cell <- NA_integer_
  alternative_start_cluster <- NA_integer_
  pseudotime_alternative <- NULL
  precomputed_pseudotime <- NULL
  endpoint_anchor_cell <- NA_integer_
  terminal_anchor_cluster <- NA_integer_
  terminal_leaf_eligible <- NA
  terminal_leaf_clusters <- integer()
  cat("  Step 3c/7: Blind root ... ")
  raw_root_labels <- as.character(meta[[fine_col]])
  root_mask <- raw_root_labels %in% c("root", "Root", "STEM", "stem", "early", "Early")
  if (identical(root_strategy, "root_ensemble")) {
    detected_genes <- as.numeric(Matrix::colSums(expr > 0))
    attractor_diagnostics <- LineaLinkR::InferStartCellAttractor(
      blind_emb,
      lineage_cl,
      nn_idx,
      nn_dists,
      potency = rep(0, n_cells),
      top_fraction = llr_cfg$root_top_fraction %||% 0.05,
      n_attractors = llr_cfg$root_n_attractors %||% 0L,
      max_iter = llr_cfg$root_max_iter %||% 300L,
      tol = llr_cfg$root_tol %||% 1e-5,
      weights = llr_cfg$root_weights %||% c(0.25, 0.20, 0.30, 0.15, 0.10),
      seed = CONFIG$seed
    )
    root_k_grid <- unique(pmin(
      as.integer(llr_cfg$root_ensemble_k_grid %||% c(30L, 50L, 70L)),
      n_cells - 1L
    ))
    root_pc_grid <- unique(pmin(
      as.integer(llr_cfg$root_ensemble_pc_grid %||% c(10L, 20L, 30L)),
      ncol(blind_emb)
    ))
    root_ensemble <- LineaLinkR::InferRootEnsemble(
      blind_emb,
      lineage_cl,
      nn_idx,
      nn_dists,
      fate_entropy = attractor_diagnostics$scores$fate_entropy,
      commitment = attractor_diagnostics$scores$commitment,
      potency = NULL,
      k_grid = root_k_grid,
      pc_grid = root_pc_grid,
      n_resamples = llr_cfg$root_ensemble_resamples %||% 20L,
      min_resamples = llr_cfg$root_ensemble_min_resamples %||% 3L,
      adaptive_resampling = isTRUE(
        llr_cfg$root_ensemble_adaptive_resampling
      ),
      subsample_fraction = llr_cfg$root_ensemble_subsample_fraction %||% 0.8,
      min_direction_margin = llr_cfg$root_ensemble_direction_margin %||% 0.10,
      seed = CONFIG$seed
    )
    start_cell <- root_ensemble$start_cell
    root_support <- root_ensemble$root_support
    direction_margin <- root_ensemble$direction_margin
    direction_status <- root_ensemble$direction_status
    root_evidence_mode <- root_ensemble$evidence_mode
    root_selection <- root_ensemble$candidates
    root_ensemble_perturbations <- root_ensemble$perturbations
    root_adaptive_stopped <- root_ensemble$adaptive_stopped
    root_resamples_used <- root_ensemble$resamples_used
    alternative_region <- root_selection$candidate_region[[2L]]
    alternative_region_cells <- which(lineage_cl == alternative_region)
    alternative_region_center <- colMeans(
      blind_emb[alternative_region_cells, , drop = FALSE]
    )
    alternative_region_distance <- rowSums((
      blind_emb[alternative_region_cells, , drop = FALSE] -
        matrix(
          alternative_region_center,
          nrow = length(alternative_region_cells),
          ncol = ncol(blind_emb),
          byrow = TRUE
        )
    )^2)
    alternative_start_cell <- alternative_region_cells[
      which.min(alternative_region_distance)
    ]
    alternative_start_cluster <- as.integer(lineage_cl[alternative_start_cell])
    root_alternative <- root_selection[2L, , drop = FALSE]
    root_alternative$start_cell <- as.integer(alternative_start_cell)
    root_diagnostics <- list(
      detected_genes = detected_genes,
      attractor = attractor_diagnostics,
      root_region = root_ensemble$root_region,
      k_grid = root_k_grid,
      pc_grid = root_pc_grid,
      n_resamples = llr_cfg$root_ensemble_resamples %||% 20L,
      min_resamples = llr_cfg$root_ensemble_min_resamples %||% 3L,
      adaptive_resampling = isTRUE(
        llr_cfg$root_ensemble_adaptive_resampling
      ),
      adaptive_stopped = root_ensemble$adaptive_stopped,
      resamples_used = root_ensemble$resamples_used,
      direction_score_components = root_ensemble$direction_score_components,
      confidence_components = root_ensemble$confidence_components,
      subsample_fraction = llr_cfg$root_ensemble_subsample_fraction %||% 0.8
    )
    selected_root_potency_mode <- "unanchored_zero"
    root_anchor_used <- TRUE
    start_method <- sprintf(
      "root_ensemble:%s:support=%.3f:margin=%.3f",
      direction_status, root_support, direction_margin
    )
    cat(sprintf("%s -> cell %d\n", start_method, start_cell))
  } else if (identical(root_strategy, "attractor_energy")) {
    detected_genes <- as.numeric(Matrix::colSums(expr > 0))
    if (identical(root_potency_mode, "detected_genes_ordering")) {
      high_potency <- log1p(detected_genes)
      low_potency <- -high_potency
      high_diagnostics <- LineaLinkR::InferStartCellAttractor(
        blind_emb,
        lineage_cl,
        nn_idx,
        nn_dists,
        potency = high_potency,
        top_fraction = llr_cfg$root_top_fraction %||% 0.05,
        n_attractors = llr_cfg$root_n_attractors %||% 0L,
        max_iter = llr_cfg$root_max_iter %||% 300L,
        tol = llr_cfg$root_tol %||% 1e-5,
        weights = llr_cfg$root_weights %||% c(0.25, 0.20, 0.30, 0.15, 0.10),
        seed = CONFIG$seed
      )
      low_diagnostics <- LineaLinkR::InferStartCellAttractor(
        blind_emb,
        lineage_cl,
        nn_idx,
        nn_dists,
        potency = low_potency,
        top_fraction = llr_cfg$root_top_fraction %||% 0.05,
        n_attractors = llr_cfg$root_n_attractors %||% 0L,
        max_iter = llr_cfg$root_max_iter %||% 300L,
        tol = llr_cfg$root_tol %||% 1e-5,
        weights = llr_cfg$root_weights %||% c(0.25, 0.20, 0.30, 0.15, 0.10),
        seed = CONFIG$seed
      )
      candidate_anchor_cells <- c(
        high_diagnostics$start_cell,
        high_diagnostics$cellwise_start_cell,
        low_diagnostics$start_cell,
        low_diagnostics$cellwise_start_cell
      )
      candidate_pt_raw <- LineaLinkR:::ComputePseudotimeFuzzyMulti(
        nn_idx,
        nn_dists,
        candidate_anchor_cells,
        rank_normalize = rank_norm
      )
      candidate_pt <- cbind(
        LineaLinkR::OrientPseudotimeByPotency(
          candidate_pt_raw[, 1L],
          high_potency,
          llr_cfg$root_direction_min_potency_cor %||% 0.10
        ),
        LineaLinkR::OrientPseudotimeByPotency(
          candidate_pt_raw[, 2L],
          high_potency,
          llr_cfg$root_direction_min_potency_cor %||% 0.10
        ),
        LineaLinkR::OrientPseudotimeByPotency(
          candidate_pt_raw[, 3L],
          high_potency,
          llr_cfg$root_direction_min_potency_cor %||% 0.10
        ),
        LineaLinkR::OrientPseudotimeByPotency(
          candidate_pt_raw[, 4L],
          high_potency,
          llr_cfg$root_direction_min_potency_cor %||% 0.10
        )
      )
      high_component_score <- high_diagnostics$components$median_score[
        high_diagnostics$components$selected
      ][1L]
      low_component_score <- low_diagnostics$components$median_score[
        low_diagnostics$components$selected
      ][1L]
      root_selection <- data.table(
        candidate = c(
          "detected_genes_high_state",
          "detected_genes_high_cellwise",
          "detected_genes_low_state",
          "detected_genes_low_cellwise"
        ),
        potency_mode = rep(
          c("detected_genes_high", "detected_genes_low"),
          each = 2L
        ),
        root_representation = rep(
          c("state_aggregate", "cellwise_component"),
          2L
        ),
        anchor_cell = candidate_anchor_cells,
        blind_score = c(
          LineaLinkR::BlindEmbeddingScore(candidate_pt[, 1L], blind_emb),
          LineaLinkR::BlindEmbeddingScore(candidate_pt[, 2L], blind_emb),
          LineaLinkR::BlindEmbeddingScore(candidate_pt[, 3L], blind_emb),
          LineaLinkR::BlindEmbeddingScore(candidate_pt[, 4L], blind_emb)
        ),
        ordering_score = c(
          LineaLinkR::BlindOrderingScore(candidate_pt[, 1L], blind_emb),
          LineaLinkR::BlindOrderingScore(candidate_pt[, 2L], blind_emb),
          LineaLinkR::BlindOrderingScore(candidate_pt[, 3L], blind_emb),
          LineaLinkR::BlindOrderingScore(candidate_pt[, 4L], blind_emb)
        ),
        potency_spearman = c(
          stats::cor(candidate_pt[, 1L], high_potency, method = "spearman"),
          stats::cor(candidate_pt[, 2L], high_potency, method = "spearman"),
          stats::cor(candidate_pt[, 3L], high_potency, method = "spearman"),
          stats::cor(candidate_pt[, 4L], high_potency, method = "spearman")
        ),
        finite_fraction = colMeans(is.finite(candidate_pt)),
        component_score = c(
          as.numeric(high_component_score),
          as.numeric(high_diagnostics$cellwise_component_score),
          as.numeric(low_component_score),
          as.numeric(low_diagnostics$cellwise_component_score)
        ),
        tie_rank = 1:4
      )
      candidate_root_cells <- max.col(-t(candidate_pt), ties.method = "first")
      candidate_root_clusters <- as.integer(
        lineage_topology$labels[candidate_root_cells]
      )
      candidate_anchor_clusters <- as.integer(
        lineage_topology$labels[candidate_anchor_cells]
      )
      candidate_lineages <- lapply(
        candidate_root_clusters,
        function(root_cluster) LineaLinkR::ExtractLineages(
          geom_emb,
          lineage_topology$labels,
          as.integer(root_cluster)
        )
      )
      candidate_leaf_clusters <- lapply(
        candidate_lineages,
        function(paths) sort(unique(as.integer(vapply(
          paths,
          function(path) as.integer(path[[length(path)]]),
          integer(1L)
        ))))
      )
      candidate_anchor_pseudotime <- candidate_pt[cbind(
        candidate_anchor_cells,
        seq_len(ncol(candidate_pt))
      )]
      candidate_max_pseudotime <- vapply(
        seq_len(ncol(candidate_pt)),
        function(index) max(candidate_pt[, index], na.rm = TRUE),
        numeric(1L)
      )
      candidate_anchor_is_high_end <- candidate_anchor_pseudotime >=
        candidate_max_pseudotime - sqrt(.Machine$double.eps)
      candidate_anchor_is_leaf <- vapply(
        seq_along(candidate_anchor_clusters),
        function(index) as.integer(candidate_anchor_clusters[[index]]) %in%
          as.integer(candidate_leaf_clusters[[index]]),
        logical(1L)
      )
      terminal_leaf_gate <- candidate_anchor_is_high_end &
        candidate_anchor_is_leaf
      root_selection[, `:=`(
        terminal_root_cluster = as.integer(candidate_root_clusters),
        terminal_anchor_cluster = as.integer(candidate_anchor_clusters),
        terminal_leaf_clusters = vapply(
          candidate_leaf_clusters,
          function(clusters) paste(sort(clusters), collapse = ";"),
          character(1L)
        ),
        terminal_anchor_is_high_end = candidate_anchor_is_high_end,
        terminal_anchor_is_leaf = candidate_anchor_is_leaf,
        terminal_leaf_eligible =
          !isTRUE(llr_cfg$root_endpoint_require_terminal_leaf) |
          terminal_leaf_gate,
        selector_topology_k = as.integer(selector_topology_k)
      )]
      endpoint_selection <-
        LineaLinkR::SelectRootEndpointByPotencyAlignment(
        root_selection,
        default_mode = "detected_genes_high",
        min_alignment_gain =
          llr_cfg$root_endpoint_alignment_min_gain %||% 0.20,
        min_ordering_gain = llr_cfg$root_ordering_min_gain %||% 0.10,
        min_representation_alignment_gain =
          llr_cfg$root_representation_alignment_min_gain %||% 0.00
      )
      selected_root_potency_mode <- endpoint_selection$selected_candidate
      selected_row <- root_selection[
        candidate == selected_root_potency_mode
      ][1L]
      selected_index <- match(
        selected_root_potency_mode,
        root_selection$candidate
      )
      precomputed_pseudotime <- candidate_pt[, selected_index]
      endpoint_anchor_cell <- as.integer(selected_row$anchor_cell)
      terminal_anchor_cluster <- as.integer(
        selected_row$terminal_anchor_cluster
      )
      terminal_leaf_eligible <- as.logical(
        selected_row$terminal_leaf_eligible
      )
      terminal_leaf_clusters <- as.integer(strsplit(
        selected_row$terminal_leaf_clusters,
        ";",
        fixed = TRUE
      )[[1L]])
      oriented_root_cells <- which(
        precomputed_pseudotime == min(precomputed_pseudotime, na.rm = TRUE)
      )
      start_cell <- oriented_root_cells[[order(
        -high_potency[oriented_root_cells],
        oriented_root_cells
      )[1L]]]
      root_diagnostics <- list(
        detected_genes_high = high_diagnostics,
        detected_genes_low = low_diagnostics
      )[[endpoint_selection$selected_potency_mode]]
      root_region_n <- max(
        2L,
        ceiling(n_cells * (llr_cfg$root_top_fraction %||% 0.05))
      )
      root_diagnostics$terminal_anchor_cell <- endpoint_anchor_cell
      root_diagnostics$start_cell <- as.integer(start_cell)
      root_diagnostics$method <- ifelse(
        identical(
          endpoint_selection$selected_potency_mode,
          "detected_genes_low"
        ),
        "attractor_terminal_distance_oriented",
        ifelse(
          identical(
            endpoint_selection$selected_root_representation,
            "cellwise_component"
          ),
          "attractor_energy_cellwise_fallback",
          "attractor_energy"
        )
      )
      root_diagnostics$root_region <- order(
        precomputed_pseudotime,
        -high_potency,
        seq_len(n_cells)
      )[seq_len(root_region_n)]
      root_diagnostics$root_cluster <- as.integer(
        lineage_topology$labels[[start_cell]]
      )
      root_diagnostics$terminal_anchor_cluster <- terminal_anchor_cluster
      root_diagnostics$terminal_leaf_eligible <- terminal_leaf_eligible
      root_diagnostics$terminal_leaf_clusters <- terminal_leaf_clusters
      root_diagnostics$selector_topology_k <- as.integer(selector_topology_k)
      root_diagnostics$potency_mode <- root_potency_mode
      root_diagnostics$selected_potency_mode <- selected_root_potency_mode
      root_diagnostics$root_selection <- root_selection
      root_diagnostics$endpoint_selection <- endpoint_selection
      root_potency_vector <- high_potency
      root_anchor_used <- TRUE
      direction_status <- "expression_potency_anchored"
      root_evidence_mode <- "anchored_bidirectional_potency_ordering"
      root_input_source <- paste0(
        "expression_detected_genes_pca_knn_no_gt:",
        endpoint_selection$endpoint_selection_status
      )
      start_method <- sprintf(
        paste0(
          "attractor_energy:%s:%s:alignment_gain=%.3f:",
          "%s:ordering_gain=%.3f"
        ),
        selected_root_potency_mode,
        endpoint_selection$endpoint_selection_status,
        endpoint_selection$endpoint_alignment_gain,
        endpoint_selection$representation_selection_status,
        endpoint_selection$representation_ordering_gain
      )
      root_ordering_confidence <- max(root_selection$ordering_score)
      use_confidence_regularized_fallback <- root_ordering_confidence <
        (llr_cfg$root_ordering_fallback_threshold %||% 0.45)
      if (use_confidence_regularized_fallback) {
        fallback_k <- max(
          as.integer(llr_cfg$root_fallback_k_min %||% 8L),
          min(
            as.integer(llr_cfg$root_fallback_k_max %||% 25L),
            as.integer(
              n_cells /
                as.integer(llr_cfg$root_fallback_cells_per_cluster %||% 50L)
            )
          )
        )
        set.seed(CONFIG$seed)
        fallback_kmeans <- stats::kmeans(
          blind_emb,
          centers = fallback_k,
          nstart = as.integer(llr_cfg$root_fallback_nstart %||% 10L),
          iter.max = 500L,
          algorithm = llr_cfg$root_fallback_kmeans_algorithm %||%
            "Hartigan-Wong"
        )
        fallback_clusters <- as.integer(factor(fallback_kmeans$cluster))
        fallback_attractor <- LineaLinkR::InferStartCellAttractor(
          blind_emb,
          fallback_clusters,
          nn_idx,
          nn_dists,
          potency = rep(0, n_cells),
          top_fraction = llr_cfg$root_top_fraction %||% 0.05,
          n_attractors = llr_cfg$root_n_attractors %||% 0L,
          max_iter = llr_cfg$root_max_iter %||% 300L,
          tol = llr_cfg$root_tol %||% 1e-5,
          weights = llr_cfg$root_fallback_weights %||%
            c(0.30, 0.25, 0.20, 0.15, 0.10),
          seed = CONFIG$seed
        )
        fallback_k_grid <- unique(pmin(
          as.integer(llr_cfg$root_ensemble_k_grid %||% c(30L, 50L, 70L)),
          n_cells - 1L
        ))
        fallback_pc_grid <- unique(pmin(
          as.integer(llr_cfg$root_ensemble_pc_grid %||% c(10L, 20L, 30L)),
          ncol(blind_emb)
        ))
        fallback_ensemble <- LineaLinkR::InferRootEnsemble(
          blind_emb,
          fallback_clusters,
          nn_idx,
          nn_dists,
          fate_entropy = fallback_attractor$scores$fate_entropy,
          commitment = fallback_attractor$scores$commitment,
          potency = NULL,
          k_grid = fallback_k_grid,
          pc_grid = fallback_pc_grid,
          n_resamples = llr_cfg$root_ensemble_resamples %||% 20L,
          min_resamples = llr_cfg$root_ensemble_min_resamples %||% 3L,
          adaptive_resampling = isTRUE(
            llr_cfg$root_ensemble_adaptive_resampling
          ),
          selection_score_mode =
            llr_cfg$root_fallback_selection_score_mode %||%
              "confidence_regularized",
          subsample_fraction =
            llr_cfg$root_ensemble_subsample_fraction %||% 0.8,
          min_direction_margin =
            llr_cfg$root_ensemble_direction_margin %||% 0.10,
          seed = CONFIG$seed
        )
        lineage_cl <- fallback_clusters
        n_lineage_cl <- uniqueN(lineage_cl)
        lineage_cluster_source <- sprintf(
          paste0(
            "kmeans:%s:algorithm=%s:requested_k=%d:",
            "occupied_k=%d:purpose=low_confidence_root_fallback"
          ),
          lineage_embedding_mode,
          llr_cfg$root_fallback_kmeans_algorithm %||% "Hartigan-Wong",
          fallback_k,
          n_lineage_cl
        )
        start_cell <- fallback_ensemble$start_cell
        precomputed_pseudotime <- NULL
        endpoint_anchor_cell <- NA_integer_
        terminal_anchor_cluster <- NA_integer_
        terminal_leaf_eligible <- NA
        terminal_leaf_clusters <- integer()
        selected_root_potency_mode <-
          "confidence_regularized_root_ensemble"
        used_confidence_regularized_fallback <- TRUE
        root_support <- fallback_ensemble$root_support
        direction_margin <- fallback_ensemble$direction_margin
        direction_status <- fallback_ensemble$direction_status
        root_evidence_mode <-
          "ambiguous_confidence_regularized_root_ensemble"
        root_input_source <- paste0(
          "expression_pca_knn_stability_no_gt:ordering_confidence=",
          sprintf("%.3f", root_ordering_confidence)
        )
        root_ensemble_perturbations <- fallback_ensemble$perturbations
        root_adaptive_stopped <- fallback_ensemble$adaptive_stopped
        root_resamples_used <- fallback_ensemble$resamples_used
        root_potency_vector <- NULL
        root_anchor_used <- TRUE
        root_diagnostics$confidence_regularized_ensemble <- list(
          trigger_ordering_confidence = root_ordering_confidence,
          trigger_threshold =
            llr_cfg$root_ordering_fallback_threshold %||% 0.45,
          root_region = fallback_ensemble$root_region,
          root_support = fallback_ensemble$root_support,
          direction_margin = fallback_ensemble$direction_margin,
          direction_status = fallback_ensemble$direction_status,
          selection_score_mode = fallback_ensemble$selection_score_mode,
          candidates = fallback_ensemble$candidates,
          kmeans_algorithm = llr_cfg$root_fallback_kmeans_algorithm %||%
            "Hartigan-Wong",
          requested_k = fallback_k,
          occupied_k = n_lineage_cl,
          k_grid = fallback_k_grid,
          pc_grid = fallback_pc_grid
        )
        start_method <- sprintf(
          paste0(
            "confidence_regularized_root_ensemble:",
            "ordering=%.3f<threshold=%.3f:support=%.3f:margin=%.3f"
          ),
          root_ordering_confidence,
          llr_cfg$root_ordering_fallback_threshold %||% 0.45,
          root_support,
          direction_margin
        )
      }
      cat(sprintf("%s -> cell %d\n", start_method, start_cell))
    } else {
    root_modes <- if (identical(root_potency_mode, "detected_genes_auto")) {
      c("detected_genes_high", "detected_genes_low", "none")
    } else {
      root_potency_mode
    }
    root_candidates <- lapply(root_modes, function(mode) {
      potency <- switch(
        mode,
        detected_genes_high = log1p(detected_genes),
        detected_genes_low = -log1p(detected_genes),
        none = rep(0, n_cells),
        stop(sprintf("Unsupported root_potency_mode: %s", mode))
      )
      diag <- LineaLinkR::InferStartCellAttractor(
        blind_emb,
        lineage_cl,
        nn_idx,
        nn_dists,
        potency = potency,
        top_fraction = llr_cfg$root_top_fraction %||% 0.05,
        n_attractors = llr_cfg$root_n_attractors %||% 0L,
        max_iter = llr_cfg$root_max_iter %||% 300L,
        tol = llr_cfg$root_tol %||% 1e-5,
        weights = llr_cfg$root_weights %||% c(0.25, 0.20, 0.30, 0.15, 0.10),
        seed = CONFIG$seed
      )
      candidate_pt <- LineaLinkR::ComputePseudotimeFuzzy(
        nn_idx, nn_dists, diag$start_cell, rank_normalize = rank_norm
      )
      finite_fraction <- mean(is.finite(candidate_pt))
      blind_score <- if (finite_fraction < 0.99) -Inf else LineaLinkR::BlindEmbeddingScore(candidate_pt, blind_emb)
      component_score <- diag$components$median_score[diag$components$selected][1]
      list(
        mode = mode,
        diagnostics = diag,
        blind_score = blind_score,
        finite_fraction = finite_fraction,
        component_score = as.numeric(component_score)
      )
    })
    root_selection <- data.table(
      potency_mode = vapply(root_candidates, `[[`, character(1L), "mode"),
      start_cell = vapply(root_candidates, function(z) z$diagnostics$start_cell, integer(1L)),
      blind_score = vapply(root_candidates, `[[`, numeric(1L), "blind_score"),
      finite_fraction = vapply(root_candidates, `[[`, numeric(1L), "finite_fraction"),
      component_score = vapply(root_candidates, `[[`, numeric(1L), "component_score")
    )
    root_selection[!is.finite(blind_score), blind_score := -Inf]
    if (!any(is.finite(root_selection$blind_score))) {
      stop(sprintf("No finite blind root candidate for %s", ds_name))
    }
    root_selection[, tie_rank := seq_len(.N)]
    setorder(root_selection, -blind_score, -component_score, tie_rank)
    selected_root_potency_mode <- root_selection$potency_mode[[1L]]
    root_diagnostics <- root_candidates[[match(selected_root_potency_mode, root_modes)]]$diagnostics
    root_diagnostics$potency_mode <- root_potency_mode
    root_diagnostics$selected_potency_mode <- selected_root_potency_mode
    root_diagnostics$root_selection <- root_selection
    root_potency_vector <- switch(
      selected_root_potency_mode,
      detected_genes_high = log1p(detected_genes),
      detected_genes_low = -log1p(detected_genes),
      none = rep(0, n_cells),
      NULL
    )
    start_cell <- root_diagnostics$start_cell
    root_anchor_used <- TRUE
    direction_status <- "expression_potency_anchored"
    root_evidence_mode <- "anchored"
    root_input_source <- "expression_detected_genes_pca_knn_no_gt"
    start_method <- sprintf(
      "attractor_energy:%s:%d_attractors:%d_region",
      selected_root_potency_mode,
      length(root_diagnostics$attractor_clusters),
      length(root_diagnostics$root_region)
    )
    cat(sprintf("%s -> cell %d\n", start_method, start_cell))
    }
  } else if (identical(root_strategy, "principal_min")) {
    start_method <- "principal_min:cluster_pc1_min"
    start_cell <- LineaLinkR::InferStartCell(
      blind_emb,
      lineage_cl,
      method = "cluster_pc1_min",
      seed = CONFIG$seed
    )
    root_anchor_used <- TRUE
    cat(sprintf("%s -> cell %d\n", start_method, start_cell))
  } else if (isTRUE(llr_cfg$prefer_root_label) && any(root_mask)) {
    root_idx <- which(root_mask)
    root_center <- colMeans(blind_emb[root_idx, , drop = FALSE])
    d2_root <- rowSums((blind_emb[root_idx, , drop = FALSE] - matrix(root_center, nrow = length(root_idx), ncol = ncol(blind_emb), byrow = TRUE))^2)
    start_cell <- root_idx[which.min(d2_root)]
    root_anchor_used <- TRUE
    root_evidence_mode <- "independent_user_root"
    root_input_source <- sprintf("user_supplied_metadata_root:%s", fine_col)
    direction_status <- "externally_anchored"
    start_method <- sprintf("root_label:%s:center", raw_root_labels[start_cell])
    cat(sprintf("%s -> cell %d\n", start_method, start_cell))
  } else if (identical(start_cfg, "auto")) {
    auto <- LineaLinkR::InferStartCellAuto(
      blind_emb, lineage_cl, nn_idx, nn_dists,
      rank_normalize = rank_norm,
      seed = CONFIG$seed
    )
    start_cell <- auto$start_cell
    rank_norm <- isTRUE(auto$rank_normalize)
    start_method <- sprintf(
      "auto:%s:%s(%.3f)",
      auto$method,
      if (rank_norm) "rank" else "linear",
      auto$blind_score
    )
    cat(sprintf("%s -> cell %d\n", start_method, start_cell))
  } else {
    start_method <- start_cfg
    start_cell <- LineaLinkR::InferStartCell(
      blind_emb,
      lineage_cl,
      method = start_method,
      seed = CONFIG$seed
    )
    cat(sprintf("%s -> cell %d\n", start_method, start_cell))
  }

  start_cluster <- as.integer(lineage_cl[start_cell])
  base_nn_idx <- nn_idx
  base_nn_dists <- nn_dists
  if (isTRUE(llr_cfg$augment_mst_bridges)) {
    br <- LineaLinkR::AugmentKnnWithMstBridges(nn_idx, nn_dists, geom_emb, lineage_cl, start_cluster)
    nn_idx <- br$nn_idx
    nn_dists <- br$nn_dists
    if (!is.null(precomputed_pseudotime)) {
      endpoint_raw_pseudotime <- LineaLinkR::ComputePseudotimeFuzzy(
        nn_idx,
        nn_dists,
        endpoint_anchor_cell,
        rank_normalize = rank_norm
      )
      precomputed_pseudotime <- LineaLinkR::OrientPseudotimeByPotency(
        endpoint_raw_pseudotime,
        high_potency,
        llr_cfg$root_direction_min_potency_cor %||% 0.10
      )
      start_cell <- which.min(precomputed_pseudotime)
      start_cluster <- as.integer(lineage_cl[start_cell])
    }
    if (!is.null(root_alternative)) {
      alternative_br <- LineaLinkR::AugmentKnnWithMstBridges(
        base_nn_idx,
        base_nn_dists,
        geom_emb,
        lineage_cl,
        alternative_start_cluster
      )
      pseudotime_alternative <- LineaLinkR::ComputePseudotimeFuzzy(
        alternative_br$nn_idx,
        alternative_br$nn_dists,
        alternative_start_cell,
        rank_normalize = rank_norm
      )
    }
    cat("  [MST kNN bridges] ")
  } else if (!is.null(root_alternative)) {
    pseudotime_alternative <- LineaLinkR::ComputePseudotimeFuzzy(
      base_nn_idx,
      base_nn_dists,
      alternative_start_cell,
      rank_normalize = rank_norm
    )
  }

  cat("  Step 3d/7: Pseudotime ... ")
  pseudotime <- if (is.null(precomputed_pseudotime)) {
    LineaLinkR::ComputePseudotimeFuzzy(
      nn_idx,
      nn_dists,
      start_cell,
      rank_normalize = rank_norm
    )
  } else {
    precomputed_pseudotime
  }
  if (length(pseudotime) != n_cells || !any(is.finite(pseudotime))) {
    stop(sprintf("Invalid pseudotime output for %s", ds_name))
  }
  if (!isTRUE(root_anchor_used)) {
    pseudotime <- LineaLinkR::OrientPseudotimeBlind(pseudotime, blind_emb)
  }
  direction_diagnostics <- LineaLinkR::AuditPseudotimeDirection(
    pseudotime,
    potency = root_potency_vector,
    embedding = blind_emb,
    n_pc = llr_cfg$root_direction_n_pc %||% 3L,
    endpoint_fraction = llr_cfg$root_direction_endpoint_fraction %||% 0.05,
    min_potency_delta = llr_cfg$root_direction_min_potency_delta %||% 0.10,
    min_blind_delta = llr_cfg$root_direction_min_blind_delta %||% 0.05
  )
  direction_diagnostics$root_strategy <- root_strategy
  direction_diagnostics$root_potency_mode <- root_potency_mode
  direction_diagnostics$selected_root_potency_mode <- selected_root_potency_mode
  direction_diagnostics$ensemble_root_support <- root_support
  direction_diagnostics$ensemble_direction_margin <- direction_margin
  direction_diagnostics$ensemble_direction_status <- direction_status
  if (!identical(root_strategy, "root_ensemble") &&
      !isTRUE(used_confidence_regularized_fallback) &&
      !identical(direction_status, "externally_anchored")) {
    direction_status <- direction_diagnostics$direction_status
    direction_margin <- abs(direction_diagnostics$blind_score_margin)
  }
  if (!is.null(root_diagnostics)) {
    root_diagnostics$direction_diagnostics <- direction_diagnostics
  }
  cat("done\n")

  # ── Step 4: Lineage weights ───────────────────────────────────────
  cat("  Step 4/7: Lineage weights ... ")

  root_clusters <- as.integer(lineage_cl)
  root_cluster_source <- lineage_cluster_source
  if (lineage_topology$empty_start_warnings > 0L) {
    cat(sprintf(
      "[audited empty starts=%d] ",
      lineage_topology$empty_start_warnings
    ))
  }
  lineage_cl <- lineage_topology$labels
  start_cluster <- as.integer(lineage_cl[start_cell])
  lineage_cluster_source <- sprintf(
    "kmeans:%s:algorithm=%s:requested_k=%d:occupied_k=%d:purpose=lineage_topology",
    lineage_embedding_mode,
    lineage_topology$algorithm,
    lineage_topology$requested_k,
    lineage_topology$occupied_k
  )
  lineage_weights <- LineaLinkR::ComputeLineageWeightsTrunkAware(geom_emb, lineage_cl, start_cluster)
  if (!is.matrix(lineage_weights) || nrow(lineage_weights) != n_cells || ncol(lineage_weights) < 1L) {
    stop(sprintf("Invalid lineage weights for %s", ds_name))
  }
  n_lineages <- ncol(lineage_weights)
  lineage_paths <- LineaLinkR::ExtractLineages(
    geom_emb, lineage_cl, start_cluster
  )
  sibling_map <- if (length(lineage_paths) >= 2L) {
    LineaLinkR::BuildSiblingContrastMap(lineage_paths)
  } else {
    data.frame(
      focal_lineage = integer(),
      sibling_lineage = integer(),
      fork_cluster = integer(),
      shared_prefix_length = integer()
    )
  }
  cat(sprintf("done [%d lineages]\n", n_lineages))

  # ── Step 4b: Per-lineage pseudotime ─────────────────────────────────
  # Temporal LR uses one global ordering; weights retain lineage membership.
  pseudotime_per_lineage <- NULL
  if (n_lineages >= 2) {
    cat("  Step 4b: Per-lineage pseudotime ... shared global pseudotime\n")
  }

  # Step 5: static LR scoring.
  cat("  Step 5/7: LR scores (static) ... ")
  sp <- CONFIG$dataset_species[[ds_name]]
  if (is.null(sp)) sp <- "human"
  lr_scores <- NULL; lr_db <- NULL; expr_lr <- NULL
  if (!sp %in% c("human", "mouse")) {
    cat(sprintf("SKIP (non-human/mouse species: %s)\n", sp))
  } else {
    expr_lr <- expr %*% Matrix::Diagonal(x = as.numeric(llr_cfg$lr_scale_factor) / Matrix::colSums(expr))
    expr_lr@x <- log1p(expr_lr@x)
    lr_db <- BuildLRDB(expr_lr, species = sp)
    if (length(lr_db$ligand_rows) < 1L) {
      cat(sprintf("SKIP (no LR pairs matched for species=%s)\n", sp))
      lr_db <- NULL
    } else {
      lr_scores <- LineaLinkR::ComputeLRScores(expr_lr, celltypes, lr_db, n_perms = 100L)
      cat(sprintf("done [%d pairs]\n", nrow(lr_scores)))
    }
  }

  # ── Step 6: Temporal LR scoring ───────────────────────────────────
  temporal_lr <- NULL
  if(!no_temporal && !is.null(lr_db)) {
    cat("  Step 6/7: Temporal LR (binned) ... ")
    temporal_pseudotime <- if (!is.null(pseudotime_per_lineage)) pseudotime_per_lineage else pseudotime
    temporal_lr <- LineaLinkR::ComputeTemporalLRWithIntermediate(
      expr_lr, celltypes, lr_db, temporal_pseudotime, lineage_weights,
      n_bins = llr_cfg$n_bins %||% 12L,
      min_n_eff = llr_cfg$min_n_eff %||% 5.0
    )
    cat(sprintf("done [%d rows]\n", nrow(temporal_lr$scores)))
  } else if (is.null(lr_db)) { cat("  Step 6/7: SKIP (no LR DB for this species)\n")
  } else { cat("  Step 6/7: ABLATION (no temporal LR)\n") }

  # ── Step 7: Entropy + bifurcation + role scores ──────────────────
  entropy_raw <- h_tilde <- bif_df <- NULL
  fork_contrasts <- data.frame(status = "not_evaluable")
  fork_evidence <- data.frame(status = "not_evaluable")
  fork_sample_contrasts <- NULL
  fork_sample_summary <- NULL
  fork_sample_column <- NA_character_
  fork_inference_unit <- "cell_resampling_only"
  if(!no_entropy && !no_temporal && !is.null(temporal_lr)) {
    cat("  Step 7/7: Lineage entropy & bifurcation ... ")
    entropy_raw <- LineaLinkR::ComputeLineageEntropyRaw(temporal_lr$intermediate)
    h_tilde <- LineaLinkR::SmoothMA3(entropy_raw)
    bif_df <- LineaLinkR::FindBifurcationBins(h_tilde)
    cat(sprintf("b*: %s\n",
      paste(sprintf("L%d=%s", bif_df$lineage,
        ifelse(is.na(bif_df$b_star), "NA", bif_df$b_star)),
        collapse = ", ")))
  } else if (is.null(temporal_lr)) { cat("  Step 7/7: SKIP (no temporal LR)\n")
  } else { cat("  Step 7/7: ABLATION (no entropy)\n") }

  if (!is.null(temporal_lr) && !is.null(bif_df) &&
      nrow(temporal_lr$scores) > 0L && nrow(sibling_map) > 0L) {
    fork_contrasts <- LineaLinkR::ComputeForkContrasts(
      temporal_lr$scores,
      bif_df,
      sibling_map,
      window_bins = llr_cfg$fork_window_bins %||% 2L,
      activity_floor = llr_cfg$fork_activity_floor %||% 1e-4,
      min_n_eff = llr_cfg$min_n_eff %||% 5.0
    )
    fork_sample_column <- FindBiologicalSampleCol(meta)
    if (!is.na(fork_sample_column)) {
      sample_labels <- as.character(meta[[fork_sample_column]])
      pooled_bounds <- attr(temporal_lr$scores, "bin_bounds")
      pooled_hill_scale <- attr(temporal_lr$scores, "hill_scale")
      sample_contrast_tables <- lapply(unique(sample_labels), function(block) {
        keep <- which(sample_labels == block)
        score_table <- LineaLinkR::ComputeTemporalLRFixedBounds(
          expr_lr[, keep, drop = FALSE],
          celltypes[keep],
          lr_db,
          pseudotime[keep],
          lineage_weights[keep, , drop = FALSE],
          bin_bounds = pooled_bounds,
          hill_scale = pooled_hill_scale,
          n_bins = llr_cfg$n_bins %||% 12L,
          min_n_eff = llr_cfg$min_n_eff %||% 5.0,
          return_complete_grid = TRUE
        )
        score_table$block_label <- block
        contrast_table <- LineaLinkR::ComputeForkContrasts(
          score_table,
          bif_df,
          sibling_map,
          sample_col = "block_label",
          window_bins = llr_cfg$fork_window_bins %||% 2L,
          activity_floor = llr_cfg$fork_activity_floor %||% 1e-4,
          min_n_eff = llr_cfg$min_n_eff %||% 5.0
        )
        rm(score_table)
        invisible(gc(verbose = FALSE))
        contrast_table
      })
      fork_sample_contrasts <- rbindlist(sample_contrast_tables, fill = TRUE)
      if (any(fork_sample_contrasts$evaluable %in% TRUE &
              is.finite(fork_sample_contrasts$sibling_contrast))) {
        fork_sample_summary <- LineaLinkR::SummarizeForkContrastsBySample(
          fork_sample_contrasts,
          n_boot = llr_cfg$fork_sample_boot %||% 999L,
          seed = CONFIG$seed
        )
        fork_inference_unit <- "biological_sample_bootstrap"
      } else {
        fork_inference_unit <- "biological_samples_no_evaluable_contrast"
      }
    }
    fork_evidence <- LineaLinkR::ComposeForkEvidence(
      fork_contrasts,
      sample_support = fork_sample_summary
    )
  }

  # Map pair indices to gene names.
  pair_gene_map <- NULL
  if (!is.null(lr_scores) && !is.null(lr_db) && nrow(lr_scores) > 0) {
    cc_db <- if (sp == "mouse") CellChat::CellChatDB.mouse else CellChat::CellChatDB.human
    lr_inter <- cc_db$interaction
    gene_names <- rownames(expr)
    # Pre-count valid LR pairs for pre-allocation
    valid_mask <- vapply(seq_len(nrow(lr_inter)), function(i) {
      lig <- strsplit(lr_inter$ligand[i], "_")[[1]]
      rec <- strsplit(lr_inter$receptor[i], "_")[[1]]
      !any(is.na(match(lig, gene_names))) && !any(is.na(match(rec, gene_names)))
    }, logical(1L))
    n_valid <- sum(valid_mask)
    pid_vec  <- seq_len(n_valid)
    lig_vec  <- character(n_valid)
    rec_vec  <- character(n_valid)
    pid <- 1L
    for (i in seq_len(nrow(lr_inter))) {
      if (!valid_mask[i]) next
      lig_vec[pid] <- lr_inter$ligand[i]
      rec_vec[pid] <- lr_inter$receptor[i]
      pid <- pid + 1L
    }
    pair_gene_map <- data.table(pair_idx = pid_vec, ligand_genes = lig_vec,
                                 receptor_genes = rec_vec)
  }

  # Save the result.
  exp_meta <- LlrExpectedMeta()
  result <- list(
    dataset          = ds_name,
    timestamp        = Sys.time(),
    result_timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    run_attempt      = run_attempt,
    package_version  = as.character(packageVersion("LineaLinkR")),
    cells            = cells,
    start_cell       = start_cell,
    start_method     = start_method,
    benchmark_protocol = exp_meta$benchmark_protocol,
    result_revision  = exp_meta$result_revision,
    root_potency_mode = exp_meta$root_potency_mode,
    selected_root_potency_mode = selected_root_potency_mode,
    endpoint_anchor_cell = endpoint_anchor_cell,
    terminal_anchor_cluster = terminal_anchor_cluster,
    terminal_leaf_eligible = terminal_leaf_eligible,
    terminal_leaf_clusters = terminal_leaf_clusters,
    selector_topology_k = as.integer(selector_topology_k),
    root_support     = as.numeric(root_support),
    direction_margin = as.numeric(direction_margin),
    direction_status = as.character(direction_status),
    root_adaptive_stopped = root_adaptive_stopped,
    root_resamples_used = root_resamples_used,
    root_candidates = root_selection,
    root_alternative = root_alternative,
    pseudotime_alternative = pseudotime_alternative,
    root_evidence_mode = root_evidence_mode,
    root_input_source = root_input_source,
    lineage_embedding = exp_meta$lineage_embedding,
    n_cells          = n_cells,
    n_genes          = nrow(expr),
    n_types          = n_types,
    n_lineages       = n_lineages,
    celltype_column   = ct_col,
    cluster_input_policy = cluster_input_policy,
    root_cluster_source = root_cluster_source,
    root_clusters    = root_clusters,
    root_cluster_empty_start_warnings =
      as.integer(root_cluster_empty_start_warnings),
    lineage_cluster_source = lineage_cluster_source,
    lineage_clusters = as.integer(lineage_cl),
    lineage_topology_k_requested = lineage_topology$requested_k,
    lineage_topology_k_occupied = lineage_topology$occupied_k,
    lineage_topology_empty_start_warnings =
      lineage_topology$empty_start_warnings,
    pseudotime       = pseudotime,
    pseudotime_per_lineage = pseudotime_per_lineage,
    lineage_weights  = lineage_weights,
    lineage_paths    = lineage_paths,
    sibling_map      = sibling_map,
    lr_scores        = lr_scores,
    temporal_lr      = temporal_lr,
    entropy_raw      = entropy_raw,
    h_tilde          = h_tilde,
    bif_df           = bif_df,
    fork_contrasts   = fork_contrasts,
    fork_evidence    = fork_evidence,
    fork_sample_contrasts = fork_sample_contrasts,
    fork_sample_summary = fork_sample_summary,
    fork_sample_column = fork_sample_column,
    fork_inference_unit = fork_inference_unit,
    legacy_role_features = NULL,
    legacy_role_scores = NULL,
    root_diagnostics = root_diagnostics,
    direction_diagnostics = direction_diagnostics,
    root_selection   = root_selection,
    root_ensemble_perturbations = root_ensemble_perturbations,
    input_manifest   = dat$input_manifest,
    pair_gene_map    = pair_gene_map,
    config_snapshot  = c(exp_meta, list(
      seed      = CONFIG$seed,
      n_threads = CONFIG$n_threads,
      n_bins    = llr_cfg$n_bins %||% 12L,
      cluster_input_policy = cluster_input_policy,
      selected_root_potency_mode = selected_root_potency_mode,
      algorithm = sprintf("LineaLinkR v%s %s rev=%s potency=%s geom=%s topology_k=%d",
        packageVersion("LineaLinkR"), exp_meta$benchmark_protocol,
        exp_meta$result_revision, exp_meta$root_potency_mode,
        exp_meta$lineage_embedding, exp_meta$lineage_topology_k)
    ))
  )

  saveRDS(result, cache_file)
  cat(sprintf("\n  >>> Saved: %s\n", cache_file))
  result
}

# =============================================================================
# Main workflow.
# =============================================================================
RunLineaLinkRBenchmark <- function() {
  cat("\n===== LineaLinkR Benchmark =====\n")
  cat(sprintf("Packages loaded: LineaLinkR v%s\n",
    as.character(packageVersion("LineaLinkR"))))

  ds_names <- names(datasets)
  if (!is.null(filter_datasets)) {
    ds_names <- filter_datasets
  }

  cat(sprintf("Datasets: %d\n", length(ds_names)))
  for (ds in ds_names) cat(sprintf("  - %s\n", ds))

  all_results <- list()

  for (ds in ds_names) {
    cat(sprintf("\n%s\n", paste(rep("=", 60), collapse = "")))
    res <- RunLineaLinkROnDataset(ds)
    if (is.null(res)) next
    all_results[[ds]] <- res
  }

  # Rebuild the summary from every current-contract RDS. A focused run must not
  # truncate the workspace index to only the selected datasets.
  summary_file <- file.path(OUT, "linealinkr_summary.csv")
  summary_paths <- file.path(OUT, names(datasets), "linealinkr_complete.rds")
  summary_paths <- summary_paths[file.exists(summary_paths)]
  summary_results <- lapply(summary_paths, readRDS)
  expected_summary_meta <- LlrExpectedMeta()
  summary_current <- vapply(summary_results, function(r) {
    saved <- r$config_snapshot %||% list()
    saved$benchmark_protocol <- r$benchmark_protocol %||% saved$benchmark_protocol
    saved$result_revision <- r$result_revision %||% saved$result_revision
    LlrRdsComplete(r) &&
      identical(
        as.character(saved$benchmark_protocol),
        expected_summary_meta$benchmark_protocol
      ) &&
      identical(
        as.character(saved$result_revision),
        expected_summary_meta$result_revision
      )
  }, logical(1L))
  summary_results <- summary_results[summary_current]
  summary <- rbindlist(lapply(summary_results, function(r) {
    data.table(
      dataset     = r$dataset,
      n_cells     = r$n_cells,
      n_genes     = r$n_genes,
      n_types     = r$n_types,
      n_lineages  = r$n_lineages,
      has_pseudotime = !is.null(r$pseudotime),
      has_lr_scores  = !is.null(r$lr_scores),
      has_temporal   = !is.null(r$temporal_lr),
      has_bifurcation= !is.null(r$bif_df)
    )
  }), fill = TRUE)
  fwrite(summary, summary_file)

  cat(sprintf("\n===== Done =====\n"))
  cat(sprintf("Summary: %s\n", summary_file))
  cat(sprintf("Individual results: %s/<dataset>/\n", OUT))
  print(summary)
}

RunLineaLinkRBenchmark()
