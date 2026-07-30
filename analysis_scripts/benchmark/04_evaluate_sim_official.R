suppressPackageStartupMessages({
  library(data.table)

# path bootstrap (relative; set LLR_ASSETS or LLR_PIPELINE to override)
.llr_cmd <- commandArgs(trailingOnly = FALSE)
.llr_file <- sub("^--file=", "", grep("^--file=", .llr_cmd, value = TRUE)[1])
if (!length(.llr_file) || is.na(.llr_file) || !nzchar(.llr_file)) {
  stop("Run this script via Rscript so --file= is set")
}
.llr_here <- dirname(normalizePath(.llr_file, winslash = "/", mustWork = TRUE))
source(file.path(.llr_here, "../config/llr_bootstrap.R"), local = FALSE)

})
llr_source_config()
`%||%` <- function(x, y) if (is.null(x)) y else x
set.seed(CONFIG$seed)
Sys.setenv(OMP_NUM_THREADS = CONFIG$n_threads)
data.table::setDTthreads(CONFIG$n_threads)

args <- commandArgs(trailingOnly = TRUE)
datasets_file <- file.path(llr_pipeline_dir(), "sim_datasets.txt")
for (a in args) {
  if (grepl("^--datasets-file=", a)) datasets_file <- sub("^--datasets-file=", "", a)
}

DATASETS <- trimws(readLines(datasets_file, warn = FALSE))
DATASETS <- DATASETS[nzchar(DATASETS)]
LLR <- CONFIG$dir$linealinkr
REFD <- CONFIG$dir$reference
OUT <- CONFIG$dir$table1_main
PYTHON <- CONFIG$exe$python
PROTOCOL_TAG <- gsub("[^A-Za-z0-9_]+", "_", CONFIG$llr$benchmark_protocol %||% "unknown")
dir.create(REFD, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

LlrExpectedMeta <- function() {
  cfg <- CONFIG$llr %||% list()
  list(
    package_version = cfg$package_version %||% "0.6.4",
    benchmark_protocol = cfg$benchmark_protocol,
    result_revision = cfg$result_revision %||% "unknown",
    root_strategy = cfg$root_strategy %||% "attractor_energy",
    root_potency_mode = cfg$root_potency_mode %||% "detected_genes_high",
    root_endpoint_require_terminal_leaf = isTRUE(cfg$root_endpoint_require_terminal_leaf),
    lineage_embedding = cfg$lineage_embedding %||% "pca2",
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
    connect_components = isTRUE(cfg$connect_components)
  )
}

LlrMetaSame <- function(found, expected) {
  if (is.null(found)) return(FALSE)
  if (is.numeric(expected)) {
    return(isTRUE(all.equal(as.numeric(found), expected, tolerance = 1e-12, check.attributes = FALSE)))
  }
  identical(found, expected)
}

CanonicalPath <- function(ds) {
  DatasetPath(ds)
}

LlrRdsComplete <- function(x) {
  n <- as.integer(x$n_cells %||% length(x$pseudotime))
  isTRUE(length(n) == 1L && is.finite(n) && n > 0L) &&
    is.character(x$dataset) && length(x$dataset) == 1L && nzchar(x$dataset) &&
    is.character(x$cells) && length(x$cells) == n && !anyNA(x$cells) && uniqueN(x$cells) == n &&
    is.numeric(x$pseudotime) && length(x$pseudotime) == n && all(is.finite(x$pseudotime)) &&
    min(x$pseudotime) >= -1e-12 && max(x$pseudotime) <= 1 + 1e-12 &&
    diff(range(x$pseudotime)) > 0 &&
    is.matrix(x$lineage_weights) && nrow(x$lineage_weights) == n && ncol(x$lineage_weights) >= 1L &&
    all(is.finite(x$lineage_weights)) && all(x$lineage_weights >= 0) &&
    all(rowSums(x$lineage_weights) > 0) &&
    (is.null(x$lr_scores) || is.data.frame(x$lr_scores)) &&
    (is.null(x$temporal_lr) || (is.list(x$temporal_lr) && is.data.frame(x$temporal_lr$scores) && is.list(x$temporal_lr$intermediate))) &&
    identical(as.character(x$package_version), "0.6.4") &&
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
    is.character(x$lineage_cluster_source) && length(x$lineage_cluster_source) == 1L && nzchar(x$lineage_cluster_source) &&
    !is.null(x$config_snapshot$algorithm) &&
    !is.null(x$input_manifest)
}

LlrRdsTrusted <- function(x) {
  exp <- LlrExpectedMeta()
  saved <- x$config_snapshot %||% list()
  saved$package_version <- x$package_version %||% saved$package_version
  saved$benchmark_protocol <- x$benchmark_protocol %||% saved$benchmark_protocol
  saved$result_revision <- x$result_revision %||% saved$result_revision
  saved$root_potency_mode <- x$root_potency_mode %||% saved$root_potency_mode
  saved$lineage_embedding <- x$lineage_embedding %||% saved$lineage_embedding
  ref <- ExportRef(x$dataset)
  ref_cell_col <- RefCellIdCol(ref)
  ref_cells <- as.character(ref[[ref_cell_col]])
  input_manifest <- x$input_manifest %||% list()
  canonical_path <- normalizePath(
    CanonicalPath(x$dataset), winslash = "/", mustWork = TRUE
  )
  manifest_contract <- identical(input_manifest$content_schema, "decoded_h5ad_content_v1") &&
    identical(input_manifest$dataset, x$dataset) &&
    identical(input_manifest$canonical_path, canonical_path) &&
    identical(as.integer(input_manifest$n_cells), length(ref_cells)) &&
    identical(as.character(input_manifest$cells), ref_cells) &&
    is.character(input_manifest$genes) && length(input_manifest$genes) >= 1L &&
    !anyNA(input_manifest$genes) && uniqueN(input_manifest$genes) == length(input_manifest$genes) &&
    is.numeric(input_manifest$matrix_nnz) && input_manifest$matrix_nnz >= 0 &&
    is.numeric(input_manifest$matrix_finite_values) &&
    input_manifest$matrix_finite_values == input_manifest$matrix_stored_values &&
    all(is.finite(c(
      input_manifest$matrix_min,
      input_manifest$matrix_max,
      input_manifest$matrix_sum,
      input_manifest$matrix_sum_squares,
      input_manifest$row_anchor_sums,
      input_manifest$col_anchor_sums
    )))
  all(vapply(names(exp), function(nm) LlrMetaSame(saved[[nm]], exp[[nm]]), logical(1L))) &&
    LlrRdsComplete(x) &&
    manifest_contract &&
    identical(as.character(x$cells), ref_cells)
}

RefCellIdCol <- function(ref) {
  candidates <- c("cell_id", "CellID", "cell", "barcode")
  hit <- intersect(candidates, names(ref))
  if (length(hit)) hit[[1L]] else names(ref)[[1L]]
}

REF_CACHE <- new.env(parent = emptyenv())
ExportRef <- function(ds) {
  if (exists(ds, envir = REF_CACHE, inherits = FALSE)) {
    return(get(ds, envir = REF_CACHE, inherits = FALSE))
  }
  csv_out <- file.path(REFD, sprintf("ref_%s.csv", ds))
  unlink(csv_out, force = TRUE)
  py <- sprintf("import anndata as ad\na=ad.read_h5ad(r'%s', backed='r')\nref=a.obs.copy()\nif 'cell_id' in ref.columns:\n    ref=ref.rename(columns={'cell_id':'obs_cell_id'})\nref.index.name='cell_id'\nref.to_csv(r'%s')\n", CanonicalPath(ds), csv_out)
  status <- sys::exec_wait(PYTHON, args = c("-X", "utf8", "-c", py))
  if (status != 0L || !file.exists(csv_out)) {
    stop(sprintf("reference export failed for %s", ds))
  }
  ref <- fread(csv_out)
  assign(ds, ref, envir = REF_CACHE)
  ref
}

ARI <- function(pred, truth) {
  tbl <- table(pred, truth)
  a <- sum(choose(tbl, 2))
  b <- sum(choose(rowSums(tbl), 2)) - a
  c2 <- sum(choose(colSums(tbl), 2)) - a
  numerator <- a - (a + b) * (a + c2) / choose(sum(tbl), 2)
  denominator <- (a + b + a + c2) / 2 - (a + b) * (a + c2) / choose(sum(tbl), 2)
  if (denominator == 0) 0 else numerator / denominator
}

BranchF1 <- function(pred, truth) {
  truth_int <- as.integer(as.factor(truth))
  pred_labels <- sort(unique(pred))
  truth_labels <- sort(unique(truth_int))
  f1s <- vapply(truth_labels, function(tl) {
    max(vapply(pred_labels, function(pl) {
      tp <- sum(pred == pl & truth_int == tl)
      fp <- sum(pred == pl & truth_int != tl)
      fn <- sum(pred != pl & truth_int == tl)
      pre <- tp / (tp + fp)
      rec <- tp / (tp + fn)
      ifelse(pre + rec == 0, 0, 2 * pre * rec / (pre + rec))
    }, numeric(1L)))
  }, numeric(1L))
  mean(f1s)
}

pseudotime_all <- data.table()
lineage_all <- data.table()
branch_all <- data.table()
manifest_all <- data.table()
root_audit_all <- data.table()

for (ds in DATASETS) {
  cat(sprintf("[%s] ", ds))
  ref <- ExportRef(ds)
  cid <- RefCellIdCol(ref)
  if (!"true_pseudotime" %in% names(ref)) stop(sprintf("%s lacks true_pseudotime", ds))
  if (!"true_branch" %in% names(ref)) stop(sprintf("%s lacks true_branch", ds))
  ref_pt <- data.table(cell_id = as.character(ref[[cid]]), true_pseudotime = as.numeric(ref$true_pseudotime))
  ref_branch <- data.table(cell_id = as.character(ref[[cid]]), true_branch = as.factor(ref$true_branch))

  rds_path <- file.path(LLR, ds, "linealinkr_complete.rds")
  llr <- readRDS(rds_path)
  if (!LlrRdsTrusted(llr)) stop(sprintf("untrusted LineaLinkR RDS for %s", ds))
  if (is.null(llr$cells) || !length(llr$cells)) stop(sprintf("%s lacks saved cell ids", ds))
  if (is.null(llr$lineage_weights)) stop(sprintf("%s lacks lineage_weights", ds))

  pt_dt <- data.table(cell_id = as.character(llr$cells), pseudotime = as.numeric(llr$pseudotime))
  merged <- merge(ref_pt, pt_dt, by = "cell_id", all = FALSE)
  if (nrow(merged) != nrow(ref_pt)) stop(sprintf("exact cell_id alignment failed for %s", ds))
  valid <- is.finite(merged$pseudotime) & is.finite(merged$true_pseudotime)
  if (sum(valid) < 50L) stop(sprintf("too few finite pseudotime cells for %s", ds))
  rng_pt <- diff(range(merged$pseudotime[valid]))
  rng_gt <- diff(range(merged$true_pseudotime[valid]))
  if (rng_pt == 0 || rng_gt == 0) stop(sprintf("zero pseudotime range for %s", ds))
  pt01 <- (merged$pseudotime[valid] - min(merged$pseudotime[valid])) / rng_pt
  gt01 <- (merged$true_pseudotime[valid] - min(merged$true_pseudotime[valid])) / rng_gt
  spearman_native <- cor(pt01, gt01, method = "spearman")
  pearson_native <- cor(pt01, gt01, method = "pearson")
  rmse_native <- sqrt(mean((pt01 - gt01)^2))

  pred_lineage <- max.col(llr$lineage_weights, ties.method = "first")
  pred_dt <- data.table(cell_id = as.character(llr$cells), pred_lineage = pred_lineage)
  branch_merged <- merge(ref_branch, pred_dt, by = "cell_id", all = FALSE)
  if (nrow(branch_merged) != nrow(ref_branch)) stop(sprintf("exact branch alignment failed for %s", ds))
  lineage_metric_scope <- if (grepl("^kmeans:", llr$lineage_cluster_source %||% "")) {
    "unsupervised_kmeans_lineage_audit"
  } else {
    "annotation_conditioned_not_blind"
  }

  pseudotime_all <- rbind(pseudotime_all, data.table(
    dataset = ds,
    algorithm = "LineaLinkR",
    ground_truth = "true_pseudotime",
    gt_tier = "independent_simulation_truth",
    benchmark_protocol = CONFIG$llr$benchmark_protocol,
    result_revision = CONFIG$llr$result_revision,
    n_cells = sum(valid),
    n_aligned = nrow(merged),
    Spearman = spearman_native,
    Pearson = pearson_native,
    RMSE = rmse_native
  ))
  lineage_all <- rbind(lineage_all, data.table(
    dataset = ds,
    algorithm = "LineaLinkR",
    reference = "true_branch",
    metric_scope = lineage_metric_scope,
    lineage_cluster_source = llr$lineage_cluster_source,
    ARI = ARI(branch_merged$pred_lineage, branch_merged$true_branch)
  ))
  branch_all <- rbind(branch_all, data.table(
    dataset = ds,
    algorithm = "LineaLinkR",
    reference = "true_branch",
    metric_scope = lineage_metric_scope,
    lineage_cluster_source = llr$lineage_cluster_source,
    Branch_F1 = BranchF1(branch_merged$pred_lineage, branch_merged$true_branch)
  ))
  manifest_all <- rbind(manifest_all, data.table(
    dataset = ds,
    algorithm = "LineaLinkR",
    rds_path = rds_path,
    trusted = TRUE,
    benchmark_protocol = llr$benchmark_protocol,
    result_revision = llr$result_revision,
    root_potency_mode = llr$root_potency_mode,
    selected_root_potency_mode = llr$selected_root_potency_mode %||% NA_character_,
    lineage_embedding = llr$lineage_embedding,
    celltype_column = llr$celltype_column %||% NA_character_,
    lineage_cluster_source = llr$lineage_cluster_source,
    n_cells = llr$n_cells,
    n_genes = llr$n_genes,
    n_lineages = llr$n_lineages,
    package_config = llr$config_snapshot$algorithm %||% NA_character_
  ), fill = TRUE)
  if (!is.null(llr$root_selection)) {
    ra <- as.data.table(llr$root_selection)
    ra[, dataset := ds]
    setcolorder(ra, c("dataset", setdiff(names(ra), "dataset")))
    root_audit_all <- rbind(root_audit_all, ra, fill = TRUE)
  }
  cat(sprintf("S=%.4f Branch_F1=%.4f root=%s lin=%d\n",
              pseudotime_all[dataset == ds, Spearman],
              branch_all[dataset == ds, Branch_F1],
              llr$selected_root_potency_mode %||% llr$root_potency_mode,
              llr$n_lineages))
}

fwrite(pseudotime_all, file.path(OUT, sprintf("table1_pseudotime_%s_sim16_official.csv", PROTOCOL_TAG)))
fwrite(lineage_all, file.path(OUT, sprintf("table1_lineage_%s_sim16_official.csv", PROTOCOL_TAG)))
fwrite(branch_all, file.path(OUT, sprintf("table1_branch_f1_%s_sim16_official.csv", PROTOCOL_TAG)))
fwrite(manifest_all, file.path(OUT, sprintf("linealinkr_manifest_%s_sim16_official.csv", PROTOCOL_TAG)))
fwrite(root_audit_all, file.path(OUT, sprintf("linealinkr_root_selection_%s_sim16_official.csv", PROTOCOL_TAG)))

cat("\n===== SIM16 OFFICIAL SUMMARY =====\n")
cat(sprintf("Mean Spearman: %.4f\n", mean(pseudotime_all$Spearman)))
cat(sprintf("Median Spearman: %.4f\n", median(pseudotime_all$Spearman)))
cat(sprintf("Mean Pearson: %.4f\n", mean(pseudotime_all$Pearson)))
cat(sprintf("Mean RMSE: %.4f\n", mean(pseudotime_all$RMSE)))
cat(sprintf("Negative Spearman: %d/%d\n", sum(pseudotime_all$Spearman < 0), nrow(pseudotime_all)))
cat(sprintf("Spearman < 0.5: %d/%d\n", sum(pseudotime_all$Spearman < 0.5), nrow(pseudotime_all)))
cat(sprintf("Mean Branch F1 (lineage audit): %.4f\n", mean(branch_all$Branch_F1)))
cat(sprintf("Mean ARI (lineage audit): %.4f\n", mean(lineage_all$ARI)))
