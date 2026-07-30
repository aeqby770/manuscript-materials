suppressPackageStartupMessages({library(data.table);library(Matrix)})
SCRIPT_PATH <- normalizePath(
  sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1L]),
  winslash = "/", mustWork = TRUE
)
PIPELINE_DIR <- normalizePath(file.path(dirname(SCRIPT_PATH), ".."), winslash = "/", mustWork = TRUE)
source(file.path(PIPELINE_DIR, "config", "00_config.R"))
`%||%` <- function(x, y) if (is.null(x)) y else x
source(file.path(PIPELINE_DIR, "qa", "evaluator_contract_064.R"))
set.seed(CONFIG$seed); Sys.setenv(OMP_NUM_THREADS=CONFIG$n_threads); data.table::setDTthreads(CONFIG$n_threads)

COMPETITOR_RESULT_ROOT <- normalizePath(
  CONFIG$dir$competitor_baseline, winslash = "/", mustWork = TRUE
)
OUT <- normalizePath(
  file.path(COMPETITOR_RESULT_ROOT, "by_algorithm"),
  winslash = "/", mustWork = TRUE
)
OUT_MISC <- normalizePath(
  file.path(COMPETITOR_RESULT_ROOT, "misc"),
  winslash = "/", mustWork = TRUE
)
EVAL <- CONFIG$dir$table1_main
SUPP <- CONFIG$dir$supplement_tables
REFD <- CONFIG$dir$reference
LLR <- CONFIG$dir$linealinkr
COMPETITOR_R_EXPORT_ROOT <- normalizePath(
  CONFIG$dir$r_export, winslash = "/", mustWork = TRUE
)
dir.create(EVAL, recursive = TRUE, showWarnings = FALSE)
dir.create(SUPP, recursive = TRUE, showWarnings = FALSE)
AssertEvaluationRootSeparation(
  OUT, c(EVAL, SUPP, REFD, LLR)
)
AssertEvaluationRootSeparation(
  OUT_MISC, c(EVAL, SUPP, REFD, LLR)
)
stopifnot(
  identical(CONFIG$llr$package_version, "0.6.4"),
  is.character(CONFIG$llr$benchmark_protocol),
  length(CONFIG$llr$benchmark_protocol) == 1L,
  nzchar(CONFIG$llr$benchmark_protocol),
  is.character(CONFIG$llr$result_revision),
  length(CONFIG$llr$result_revision) == 1L,
  nzchar(CONFIG$llr$result_revision)
)
BENCHMARK_PROTOCOL <- CONFIG$llr$benchmark_protocol
LLR_RESULT_REVISION <- CONFIG$llr$result_revision
LLR_ROOT_POTENCY_MODE <- CONFIG$llr$root_potency_mode
LLR_LINEAGE_EMBEDDING <- CONFIG$llr$lineage_embedding
PROTOCOL_TAG <- gsub("[^A-Za-z0-9_]+", "_", BENCHMARK_PROTOCOL)
RUN_ARGS <- commandArgs(trailingOnly = TRUE)
ALGO_ARG <- RUN_ARGS[grepl("^--algos=", RUN_ARGS)]
INCREMENTAL_MODE <- length(ALGO_ARG) == 1L
RECOMPUTE_ALGOS <- if (INCREMENTAL_MODE) {
  strsplit(sub("^--algos=", "", ALGO_ARG), ",", fixed = TRUE)[[1L]]
} else {
  character(0)
}
BASE_TABLE_PATH <- file.path(EVAL, sprintf("table1_pseudotime_%s.csv", PROTOCOL_TAG))
BASE_M1 <- if (INCREMENTAL_MODE) fread(BASE_TABLE_PATH) else data.table()

LlrExpectedMeta <- function() {
  cfg <- CONFIG$llr
  list(
    package_version = cfg$package_version,
    benchmark_protocol = cfg$benchmark_protocol,
    result_revision = cfg$result_revision,
    root_strategy = cfg$root_strategy %||% "attractor_energy",
    root_potency_mode = cfg$root_potency_mode %||% "detected_genes_high",
    root_endpoint_require_terminal_leaf = isTRUE(
      cfg$root_endpoint_require_terminal_leaf
    ),
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

LLR_TRUST_CACHE <- new.env(parent = emptyenv())
LlrRdsTrusted <- function(path) {
  cache_key <- normalizePath(path, winslash = "/", mustWork = FALSE)
  if (exists(cache_key, envir = LLR_TRUST_CACHE, inherits = FALSE)) {
    return(get(cache_key, envir = LLR_TRUST_CACHE, inherits = FALSE))
  }
  if (!file.exists(path)) return(FALSE)
  x <- readRDS(path)
  exp <- LlrExpectedMeta()
  saved <- x$config_snapshot %||% list()
  saved$package_version <- x$package_version %||% saved$package_version
  saved$benchmark_protocol <- x$benchmark_protocol %||% saved$benchmark_protocol
  saved$result_revision <- x$result_revision %||% saved$result_revision
  saved$root_potency_mode <- x$root_potency_mode %||% saved$root_potency_mode
  saved$lineage_embedding <- x$lineage_embedding %||% saved$lineage_embedding
  force_auto_clusters <- x$dataset %chin% BENCHMARK_FORCE_AUTO_CLUSTER_DATASETS
  expected_cluster_policy <- c(
    "metadata_or_expression_auto_clusters",
    "expression_auto_clusters"
  )[1L + as.integer(force_auto_clusters)]
  cluster_input_contract <- !force_auto_clusters || (
    identical(
      as.character(x$cluster_input_policy %||% saved$cluster_input_policy),
      expected_cluster_policy
    ) &&
      identical(as.character(x$celltype_column), "auto_clusters")
  )
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
  ValidateLlr064Identity(
    x,
    exp$package_version,
    exp$benchmark_protocol,
    exp$result_revision,
    path = cache_key
  )
  trusted <- cluster_input_contract &&
    all(vapply(names(exp), function(nm) LlrMetaSame(saved[[nm]], exp[[nm]]), logical(1L))) &&
    LlrRdsComplete(x) &&
    manifest_contract &&
    identical(as.character(x$cells), ref_cells)
  assign(cache_key, trusted, envir = LLR_TRUST_CACHE)
  trusted
}
PYTHON <- CONFIG$exe$python %||% "E:/Python/biotrade/Scripts/python.exe"

# Official/adapters only; proxy inline implementations excluded (STREAM, Wishbone, PHLOWER baseline).
PSEUDOTIME_ALGOS_INCLUDED <- BENCHMARK_PSEUDOTIME_ALGORITHMS
PSEUDOTIME_ALGOS_EXCLUDED <- character(0)
# Algorithms whose pseudotime axis has no user- or biology-defined root;
# signed rho is retained separately and the headline score is orientation-aligned.
NO_ROOT_ALGOS <- c("slingshot", "SCORPIUS", "pCreode", "TSCAN")
ROOT_INPUT_MODE <- c(
  LineaLinkR = "blind_data_driven",
  slingshot = "root_free_native",
  DPT = "cluster_label_heuristic",
  SCORPIUS = "root_free_native",
  MultistageOT = "time_or_embedding_anchors",
  scFates = "internal_graph_root",
  pCreode = "arbitrary_graph_node",
  monocle3 = "cluster_label_or_timepoint_heuristic",
  TSCAN = "root_free_longest_path",
  Palantir = "cluster_label_heuristic",
  StaVia = "cluster_label_heuristic",
  scTour = "internal_latent_time"
)
MIN_PSEUDOTIME_CELL_COVERAGE <- 1.00
TI_PLANNED <- CONFIG$algorithm_groups$ti_trajectory
ORDERING_DATASETS <- BENCHMARK_ORDERING_DATASETS
CORE_DATASETS <- BENCHMARK_COMPLETE_CORE_DATASETS
BRANCH_DATASETS <- BENCHMARK_BRANCH_DATASETS
BRANCH_TERMINAL_MIN_PSEUDOTIME <- BENCHMARK_BRANCH_TERMINAL_MIN_PSEUDOTIME
stopifnot(
  identical(PSEUDOTIME_ALGOS_INCLUDED, c("LineaLinkR", TI_PLANNED)),
  !anyDuplicated(PSEUDOTIME_ALGOS_INCLUDED),
  all(ORDERING_DATASETS %in% names(datasets)),
  all(CORE_DATASETS %in% ORDERING_DATASETS),
  all(BRANCH_DATASETS %in% ORDERING_DATASETS),
  identical(names(BRANCH_TERMINAL_MIN_PSEUDOTIME), BRANCH_DATASETS)
)

# Metrics consume the current fair-matrix decision, then independently decode
# each trajectory below. A raw CSV cannot bypass the current run registry.
FAIR_MATRIX_PATH <- file.path(
  CONFIG$dir$paper_tables, "_registry", "fair_comparison_matrix.csv"
)
FAIR_MATRIX_REQUIRED_COLUMNS <- c(
  "algorithm", "dataset", "dataset_role", "output_status",
  "n_reference", "n_output", "n_aligned", "n_missing", "n_unexpected",
  "n_duplicate", "n_nonfinite", "nonzero_range", "current_run_status",
  "benchmark_protocol", "result_revision"
)
stopifnot(file.exists(FAIR_MATRIX_PATH))
FAIR_MATRIX <- fread(FAIR_MATRIX_PATH, showProgress = FALSE)
stopifnot(
  all(FAIR_MATRIX_REQUIRED_COLUMNS %in% names(FAIR_MATRIX)),
  nrow(FAIR_MATRIX) == length(PSEUDOTIME_ALGOS_INCLUDED) * length(ORDERING_DATASETS),
  uniqueN(FAIR_MATRIX, by = c("algorithm", "dataset")) == nrow(FAIR_MATRIX),
  setequal(FAIR_MATRIX$algorithm, PSEUDOTIME_ALGOS_INCLUDED),
  setequal(FAIR_MATRIX$dataset, ORDERING_DATASETS),
  all(FAIR_MATRIX$benchmark_protocol == BENCHMARK_PROTOCOL),
  all(FAIR_MATRIX$result_revision == LLR_RESULT_REVISION)
)
FAIR_CORE <- FAIR_MATRIX[dataset_role == "complete_core"]
FAIR_CORE_COMPETITORS <- FAIR_CORE[algorithm != "LineaLinkR"]
stopifnot(
  nrow(FAIR_CORE) == length(PSEUDOTIME_ALGOS_INCLUDED) * length(CORE_DATASETS),
  all(FAIR_CORE$output_status == "complete_output"),
  all(FAIR_CORE$n_output == FAIR_CORE$n_reference),
  all(FAIR_CORE$n_aligned == FAIR_CORE$n_reference),
  all(FAIR_CORE$n_missing == 0L),
  all(FAIR_CORE$n_unexpected == 0L),
  all(FAIR_CORE$n_duplicate == 0L),
  all(FAIR_CORE$n_nonfinite == 0L),
  all(FAIR_CORE$nonzero_range),
  all(FAIR_CORE_COMPETITORS$current_run_status == "OK")
)
setkey(FAIR_MATRIX, algorithm, dataset)
FairOutputEligible <- function(algo, ds) {
  row <- FAIR_MATRIX[list(algo, ds)]
  nrow(row) == 1L && identical(row$output_status[[1L]], "complete_output")
}

RefCellIdCol <- function(ref) {
  if ("cell_id" %in% names(ref)) return("cell_id")
  names(ref)[1L]
}

# ── Helper: export reference data from the canonical h5ad ─────────
REF_CACHE <- new.env(parent = emptyenv())
ExportRef <- function(ds_name) {
  if (exists(ds_name, envir = REF_CACHE, inherits = FALSE)) {
    return(get(ds_name, envir = REF_CACHE, inherits = FALSE))
  }
  csv_out <- file.path(REFD, sprintf("ref_%s.csv", ds_name))
  canon <- DatasetPath(ds_name)
  unlink(csv_out, force = TRUE)
  py <- sprintf(r"(
import anndata as ad
a = ad.read_h5ad(r'%s', backed='r')
ref = a.obs.copy()
if 'cell_id' in ref.columns:
    ref = ref.rename(columns={'cell_id': 'obs_cell_id'})
ref.index.name = 'cell_id'
ref.to_csv(r'%s')
print('OK')
)", canon, csv_out)
  status <- sys::exec_wait(PYTHON, args = c("-X", "utf8", "-c", py))
  if (status != 0L || !file.exists(csv_out)) {
    stop(sprintf("Reference export failed for %s", ds_name))
  }
  ref <- fread(csv_out)
  assign(ds_name, ref, envir = REF_CACHE)
  ref
}

# ── Helpers: robust ordinal-time parsing (M1 fix) ─────────────────
# Extract embedded numeric time from semi-structured labels:
#   "day 10"->10, "d14_PBS"->14, "3dpf"->3, "E8.5"->8.5, "E7.25"->7.25
.extract_ordinal_time <- function(vals) {
  s <- as.character(vals)
  num <- suppressWarnings(as.numeric(
    sub(".*?(-?[0-9]+(?:\\.[0-9]+)?).*", "\\1", s, perl = TRUE)
  ))
  num[!grepl("[0-9]", s)] <- NA_real_
  num
}

# Natural-sorted ordinal codes; NEVER bare lexical as.factor (mis-orders
# "day 2" vs "day 10", "E9.5" vs "E12.5"). Preserves full input length.
.natural_ordinal_codes <- function(vals) {
  s <- as.character(vals)
  keep <- !is.na(s) & s != ""
  lev <- unique(s[keep])
  if (requireNamespace("gtools", quietly = TRUE)) {
    lev <- gtools::mixedsort(lev)
  } else {
    key <- suppressWarnings(as.numeric(
      sub(".*?(-?[0-9]+(?:\\.[0-9]+)?).*", "\\1", lev, perl = TRUE)
    ))
    lev <- lev[order(key, lev, na.last = TRUE)]
  }
  out <- rep(NA_real_, length(s))
  out[keep] <- match(s[keep], lev)
  out
}

# ── Helper: find ground truth pseudotime column ───────────────────
FindGroundTruth <- function(ref, ds_name) {
  if (!ds_name %in% ORDERING_DATASETS) return(NULL)
  # Prefer independent time labels; Palantir is only an external reference GT.
  registered_gt_cols <- c(
    lung = "time_point",
    larry_invitro = "time_info",
    human_fetal_liver_organoid = "timepoint",
    human_fetal_pancreas = "pcw",
    zebrafish_tracerseq = "time_hours"
  )
  gt_specs <- list(
    list(col = "true_pseudotime", tier = "independent"),
    list(col = "timepoint", tier = "independent"),
    list(col = "day", tier = "independent"),
    list(col = "age(days)", tier = "independent"),
    list(col = "stage", tier = "independent")
  )
  registered_col <- unname(registered_gt_cols[ds_name])
  if (length(registered_col) == 1L && !is.na(registered_col)) {
    gt_specs <- c(list(list(col = registered_col, tier = "independent")), gt_specs)
  }
  cid <- RefCellIdCol(ref)
  for (spec in gt_specs) {
    col <- spec$col
    if (!col %in% names(ref)) next
    vals <- ref[[col]]
    numeric_vals <- suppressWarnings(as.numeric(vals))
    if (sum(!is.na(numeric_vals)) < 10) {
      # M1 fix: parse embedded ordinal time BEFORE any factor coding.
      extracted <- .extract_ordinal_time(vals)
      if (sum(!is.na(extracted)) >= 10 &&
          length(unique(extracted[!is.na(extracted)])) >= 2) {
        numeric_vals <- extracted
      } else {
        numeric_vals <- .natural_ordinal_codes(vals)
      }
    }
    valid <- !is.na(numeric_vals)
    if (sum(valid) > 10 && length(unique(numeric_vals[valid])) >= 2) {
      return(list(
        column = col,
        values = numeric_vals,
        cell_id = ref[[cid]],
        gt_tier = spec$tier
      ))
    }
  }
  return(NULL)
}

# ── Helper: find cluster column ───────────────────────────────────
FindClusters <- function(ref) {
  ct_cols <- c("true_branch","true_branch_numeric","branch_truth","lineage_truth")
  for(col in ct_cols) {
    if(col %in% names(ref)) {
      vals <- ref[[col]]
      valid <- !is.na(vals) & vals != ""
      if(sum(valid) > 0 && length(unique(vals[valid])) >= 2)
        return(list(column=col, values=as.factor(vals), cell_id=ref[[RefCellIdCol(ref)]]))
    }
  }
  return(NULL)
}

PseudoPath <- function(algo, ds) {
  switch(algo,
    LineaLinkR = file.path(LLR, ds, "linealinkr_complete.rds"),
    slingshot = file.path(OUT, "slingshot", ds, "slingshot_pseudotime.csv"),
    DPT = file.path(OUT_MISC, sprintf("dpt_%s_pseudotime.csv", ds)),
    SCORPIUS = file.path(OUT, "SCORPIUS", ds, "scorpius_pseudotime.csv"),
    MultistageOT = file.path(
      OUT_MISC,
      sprintf("multistageot_%s_pseudotime.csv", ds)
    ),
    scFates = file.path(OUT_MISC, sprintf("scfates_%s_pseudotime.csv", ds)),
    FateID = file.path(OUT, "FateID", ds, "fateid_bias.csv"),
    pCreode = file.path(OUT_MISC, sprintf("pcreode_%s_pseudotime.csv", ds)),
    monocle3 = file.path(OUT, "monocle3", ds, "monocle3_pseudotime.csv"),
    TSCAN = file.path(OUT, "TSCAN", ds, "tscan_pseudotime.csv"),
    Palantir = file.path(OUT_MISC, sprintf("palantir_%s_result.h5ad", ds)),
    STREAM = file.path(OUT_MISC, sprintf("stream_%s_pseudotime.csv", ds)),
    Wishbone = file.path(OUT_MISC, sprintf("wishbone_%s_pseudotime.csv", ds)),
    PHLOWER = file.path(OUT_MISC, sprintf("phlower_%s_pseudotime.csv", ds)),
    StaVia = file.path(OUT_MISC, sprintf("stavia_%s_pseudotime.csv", ds)),
    scTour = file.path(OUT_MISC, sprintf("sctour_%s_pseudotime.csv", ds)),
    NULL
  )
}

ReadPseudoFromH5ad <- function(path, obs_cols = c("pal_pt", "pseudotime", "dpt_pseudotime")) {
  col_expr <- paste(sprintf("'%s'", obs_cols), collapse = ",")
  py <- sprintf(r"(
import anndata as ad
import pandas as pd
import sys
a = ad.read_h5ad(r'%s')
idx = a.obs_names.astype(str)
column = next((name for name in [%s] if name in a.obs.columns), None)
if column is None:
    raise KeyError('No registered pseudotime column')
out = pd.DataFrame({'cell_id': idx, 'pseudotime': pd.to_numeric(a.obs[column], errors='coerce')})
sys.stdout.write(out.to_csv(index=False, lineterminator='\n'))
)", path, col_expr)
  decoded <- sys::exec_internal(
    PYTHON,
    args = c("-X", "utf8", "-c", py),
    error = FALSE
  )
  if (!identical(decoded$status, 0L)) {
    stop(sprintf("Python pseudotime extraction failed for %s (exit code %d)", path, decoded$status))
  }
  csv_text <- rawToChar(decoded$stdout)
  if (!nzchar(csv_text)) return(NULL)
  dt <- fread(text = csv_text)
  if (!all(c("cell_id", "pseudotime") %in% names(dt))) return(NULL)
  data.table(cell_id = as.character(dt$cell_id), pseudotime = as.numeric(dt$pseudotime), cell_alignment_mode = "explicit_h5ad_obs_names")
}

NormalizePseudoDT <- function(dt) {
  if (is.null(dt) || !nrow(dt)) return(NULL)
  if ("cell" %in% names(dt) && !"cell_id" %in% names(dt)) dt[, cell_id := as.character(cell)]
  pt_col <- intersect(
    c("pseudotime", "pt", "pal_pt", "dpt_pseudotime", "Lineage1", "Lineage2"),
    names(dt)
  )[1]
  if (is.na(pt_col) || !"cell_id" %in% names(dt)) return(NULL)
  out <- data.table(cell_id = as.character(dt$cell_id), pseudotime = as.numeric(dt[[pt_col]]), cell_alignment_mode = "explicit_cell_id")
  out[is.finite(pseudotime)]
}

ReadSlingshotPseudoDT <- function(ds) {
  p <- file.path(OUT, "slingshot", ds, "slingshot_pseudotime.csv")
  if (!file.exists(p)) return(NULL)
  raw <- fread(p)
  strict <- NormalizePseudoDT(raw)
  if (!is.null(strict)) return(strict)
  cells <- readLines(file.path(COMPETITOR_R_EXPORT_ROOT, ds, "cells.txt"), warn = FALSE)
  weights <- as.matrix(fread(file.path(OUT, "slingshot", ds, "slingshot_weights.csv")))
  pt_cols <- grep("^Lineage[0-9]+$", names(raw), value = TRUE)
  pt <- as.matrix(raw[, ..pt_cols])
  valid <- is.finite(pt) & is.finite(weights)
  weighted_pt <- pt * weights
  weighted_pt[!valid] <- 0
  valid_weights <- weights
  valid_weights[!valid] <- 0
  weighted_sum <- rowSums(weighted_pt)
  weight_sum <- rowSums(valid_weights)
  global_pt <- weighted_sum / weight_sum
  data.table(cell_id = as.character(cells), pseudotime = as.numeric(global_pt), cell_alignment_mode = "cells_txt_row_order_contract")[is.finite(pseudotime)]
}

ReadMonocle3PseudoDT <- function(ds) {
  csv <- file.path(OUT, "monocle3", ds, "monocle3_pseudotime.csv")
  if (file.exists(csv) && file.info(csv)$size > 500) {
    return(NormalizePseudoDT(fread(csv)))
  }
  rds <- file.path(OUT, "monocle3", ds, "monocle3_cds.rds")
  if (!file.exists(rds)) return(NULL)
  cds <- readRDS(rds)
  pt <- NULL
  if (requireNamespace("monocle3", quietly = TRUE)) {
    pt <- monocle3::pseudotime(cds)
  }
  if (is.null(pt) && !is.null(cds@principal_graph_aux$UMAP$pseudotime)) {
    pt <- cds@principal_graph_aux$UMAP$pseudotime
  }
  if (is.null(pt) || !length(pt)) return(NULL)
  data.table(
    cell_id = as.character(names(pt)),
    pseudotime = as.numeric(pt),
    cell_alignment_mode = "explicit_monocle3_names"
  )[is.finite(pseudotime)]
}

# ── Helper: read pseudotime as cell_id + pseudotime ────────────────
ReadPseudoDTUncached <- function(algo, ds) {
  if (algo == "FateID") return(NULL)
  if (algo == "slingshot") return(ReadSlingshotPseudoDT(ds))
  if (algo == "monocle3") {
    pt <- ReadMonocle3PseudoDT(ds)
    if (!is.null(pt)) return(pt)
  }
  p <- PseudoPath(algo, ds)
  if (is.null(p) || !length(p) || is.na(p) || !file.exists(p)) return(NULL)
  if (algo == "LineaLinkR" && !LlrRdsTrusted(p)) return(NULL)
  if (grepl("\\.rds$", p, ignore.case = TRUE)) {
    rds <- readRDS(p)
    ids <- rds$cells
    if (is.null(ids) || !length(ids)) return(NULL)
    return(data.table(cell_id = as.character(ids), pseudotime = as.numeric(rds$pseudotime), cell_alignment_mode = "explicit_rds_cells"))
  }
  if (grepl("\\.h5ad$", p, ignore.case = TRUE)) {
    return(ReadPseudoFromH5ad(p))
  }
  raw <- fread(p)
  out <- NormalizePseudoDT(raw)
  if (!is.null(out)) return(out)
  # Fallback: CSV has no cell_id column; try r_export cells.txt for barcodes.
  cells_file <- file.path(COMPETITOR_R_EXPORT_ROOT, ds, "cells.txt")
  if (file.exists(cells_file) && nrow(raw) > 0L) {
    cells_vec <- readLines(cells_file)
    if (length(cells_vec) == nrow(raw)) {
      pt_col <- intersect(c("pseudotime","pt","pal_pt","dpt_pseudotime","Lineage1","Lineage2"), names(raw))[1]
      if (!is.na(pt_col)) {
       out <- data.table(cell_id = cells_vec, pseudotime = as.numeric(raw[[pt_col]]), cell_alignment_mode = "cells_txt_row_order_contract")
        return(out[is.finite(pseudotime)])
      }
    }
  }
  NULL
}

PSEUDO_READ_CACHE <- new.env(parent = emptyenv())
ReadPseudoDT <- function(algo, ds) {
  cache_key <- paste(algo, ds, sep = "::")
  if (exists(cache_key, envir = PSEUDO_READ_CACHE, inherits = FALSE)) {
    return(get(cache_key, envir = PSEUDO_READ_CACHE, inherits = FALSE)$value)
  }
  out <- ReadPseudoDTUncached(algo, ds)
  assign(cache_key, list(value = out), envir = PSEUDO_READ_CACHE)
  out
}

HasPseudoOutput <- function(algo, ds) {
  if (algo == "FateID") return(FALSE)
  p <- PseudoPath(algo, ds)
  !is.null(p) && length(p) == 1L && !is.na(p) && file.exists(p)
}

AuditPseudoCoverage <- function(algo, ds, ref) {
  has_artifact <- HasPseudoOutput(algo, ds)
  ref_ids <- unique(as.character(ref[[RefCellIdCol(ref)]]))
  pt_dt <- ReadPseudoDT(algo, ds)
  output_ids <- if (is.null(pt_dt)) character(0) else unique(as.character(pt_dt$cell_id))
  n_output <- length(output_ids)
  duplicate_ids <- !is.null(pt_dt) && anyDuplicated(pt_dt$cell_id) > 0L
  n_aligned <- sum(output_ids %in% ref_ids)
  n_missing <- sum(!ref_ids %in% output_ids)
  n_unexpected <- sum(!output_ids %in% ref_ids)
  coverage <- n_aligned / length(ref_ids)
  status <- fcase(
    !has_artifact, "missing_output",
    is.null(pt_dt) || nrow(pt_dt) == 0L, "unreadable_output",
    duplicate_ids, "duplicate_cell_ids",
    n_unexpected > 0L, "unexpected_cell_ids",
    n_missing > 0L || coverage < MIN_PSEUDOTIME_CELL_COVERAGE, "partial_output",
    default = "complete_output"
  )
  data.table(
    algorithm = algo,
    dataset = ds,
    has_artifact = has_artifact,
    n_reference = length(ref_ids),
    n_output_finite = n_output,
    n_aligned_finite = n_aligned,
    n_missing_cells = n_missing,
    n_unexpected_cells = n_unexpected,
    cell_coverage = coverage,
    coverage_status = status
  )
}

AlignPseudoToRef <- function(algo, ds, ref, gt) {
  pt_dt <- ReadPseudoDT(algo, ds)
  if (is.null(pt_dt) || nrow(pt_dt) == 0L) return(NULL)
  if (anyDuplicated(pt_dt$cell_id) > 0L) return(NULL)
  ref_dt <- data.table(cell_id = as.character(gt$cell_id), gt = as.numeric(gt$values))
  if (!setequal(as.character(pt_dt$cell_id), ref_dt$cell_id)) return(NULL)

  merge_by_id <- function(pt) merge(ref_dt, pt, by = "cell_id", all = FALSE)

  merged <- merge_by_id(pt_dt)
  coverage <- nrow(merged) / nrow(ref_dt)
  if (
    nrow(merged) < 50L ||
    coverage < MIN_PSEUDOTIME_CELL_COVERAGE
  ) return(NULL)
  list(
    pt = merged$pseudotime,
    gt = merged$gt,
    n_aligned = nrow(merged),
    n_reference = nrow(ref_dt),
    cell_coverage = coverage,
    cell_alignment_mode = paste(sort(unique(merged$cell_alignment_mode)), collapse = ";")
  )
}

# ── Helper: ARI computation ───────────────────────────────────────
ARI <- function(pred, truth) {
  tbl <- table(pred, truth)
  a <- sum(choose(tbl, 2)); b <- sum(choose(rowSums(tbl), 2)) - a
  c2 <- sum(choose(colSums(tbl), 2)) - a; d <- choose(sum(tbl), 2) - a - b - c2
  numerator <- a - (a+b)*(a+c2)/choose(sum(tbl), 2)
  denominator <- (a+b+a+c2)/2 - (a+b)*(a+c2)/choose(sum(tbl), 2)
  if(denominator == 0) 0 else numerator / denominator
}

NMI <- function(pred, truth) {
  tbl <- table(pred, truth); n <- sum(tbl)
  p_i <- rowSums(tbl)/n; p_j <- colSums(tbl)/n
  mi <- sum(tbl/n * log((tbl/n) / (p_i %*% t(p_j)) + 1e-10), na.rm=TRUE)
  h_i <- -sum(p_i * log(p_i + 1e-10)); h_j <- -sum(p_j * log(p_j + 1e-10))
  if(h_i == 0 || h_j == 0) 0 else 2 * mi / (h_i + h_j)
}

Jaccard <- function(pred, truth) {
  pred <- as.integer(as.factor(pred))
  truth <- as.integer(as.factor(truth))
  pred_labels <- sort(unique(pred))
  truth_labels <- sort(unique(truth))
  jac <- sapply(truth_labels, function(tl) {
    max(sapply(pred_labels, function(pl) {
      tp <- sum(pred == pl & truth == tl)
      fp <- sum(pred == pl & truth != tl)
      fn <- sum(pred != pl & truth == tl)
      den <- tp + fp + fn
      ifelse(den == 0, 0, tp / den)
    }))
  })
  mean(jac)
}

# ══════════════════════════════════════════════════════════════════
# MAIN: Loop through all datasets
# ══════════════════════════════════════════════════════════════════

pseudotime_algos <- if (INCREMENTAL_MODE) {
  intersect(PSEUDOTIME_ALGOS_INCLUDED, RECOMPUTE_ALGOS)
} else {
  PSEUDOTIME_ALGOS_INCLUDED
}

# The active registry has already removed duplicate cell collections.  Keep
# the flagging hook empty so old duplicate aliases cannot enter summary tables.
DUPLICATE_SETS <- list()
m1_all <- data.table()
m2_all <- data.table()

for(ds in ORDERING_DATASETS) {
  cat(sprintf("\n[%s] ", ds))
  ref <- ExportRef(ds)
  if(is.null(ref)) { cat("NO REF\n"); next }
  
  gt <- FindGroundTruth(ref, ds)
  branch_ref <- ref
  if (ds %in% BRANCH_DATASETS) {
    branch_threshold <- unname(BRANCH_TERMINAL_MIN_PSEUDOTIME[[ds]])
    branch_ref <- ref[
      is.finite(as.numeric(true_pseudotime)) &
        as.numeric(true_pseudotime) >= branch_threshold
    ]
  }
  cl <- if (ds %in% BRANCH_DATASETS) FindClusters(branch_ref) else NULL
  
  cat(sprintf("time=%s cluster=%s n=%d\n",
    if(is.null(gt)) "NO" else gt$column,
    if(is.null(cl)) "NO" else cl$column, nrow(ref)))
  
  if(is.null(gt) && is.null(cl)) next
  
  # ── Metric 1: Pseudotime accuracy (cell_id join, no oracle lineage pick) ──
  if(!is.null(gt)) {
    for(algo in pseudotime_algos) {
      if (!FairOutputEligible(algo, ds)) next
      aligned <- AlignPseudoToRef(algo, ds, ref, gt)
      if(is.null(aligned)) next
      pt <- aligned$pt
      gt_vals_local <- aligned$gt
      valid <- !is.na(pt) & !is.na(gt_vals_local) & is.finite(pt) & is.finite(gt_vals_local)
      if(sum(valid) < 50) next
      rng_pt <- diff(range(pt[valid])); rng_gt <- diff(range(gt_vals_local[valid]))
      if(rng_pt == 0 || rng_gt == 0) next
      pt01 <- (pt[valid] - min(pt[valid])) / rng_pt
      gt01 <- (gt_vals_local[valid] - min(gt_vals_local[valid])) / rng_gt
      sp_signed <- cor(pt01, gt01, method="spearman")
      # M2: orientation alignment for algorithms without a deterministic root
      # (arbitrary pseudotime direction). Rooted methods keep raw signed rho.
      flip <- algo %in% NO_ROOT_ALGOS && is.finite(sp_signed) && sp_signed < 0
      pt01o <- if (flip) 1 - pt01 else pt01
      sp <- cor(pt01o, gt01, method="spearman")
      pe <- cor(pt01o, gt01, method="pearson")
      ke <- cor(pt01o, gt01, method="kendall")
      rmse <- sqrt(mean((pt01o - gt01)^2))
      mae <- mean(abs(pt01o - gt01))
       m1_all <- rbind(m1_all,
         data.table(dataset=ds, algorithm=algo, ground_truth=gt$column,
                    gt_tier=gt$gt_tier %||% "unknown",
                    benchmark_protocol=BENCHMARK_PROTOCOL,
                    root_input_mode=unname(ROOT_INPUT_MODE[algo]),
                    direction_policy=ifelse(algo %in% NO_ROOT_ALGOS, "gt_aligned_if_negative", "signed"),
                    n_cells=sum(valid), n_aligned=aligned$n_aligned,
                   n_reference=aligned$n_reference,
                   cell_coverage=aligned$cell_coverage,
                   cell_alignment_mode=aligned$cell_alignment_mode,
                   Spearman=sp, Pearson=pe, Kendall=ke, RMSE=rmse, MAE=mae,
                   Spearman_signed=sp_signed, orientation_flipped=flip,
                   root_defined=!(algo %in% NO_ROOT_ALGOS)))
    }
  }
  
  # ── Metric 2: Lineage ARI ──
  if(!is.null(cl)) {
    cl_dt <- data.table(cell_id = as.character(cl$cell_id), truth = cl$values)
    # slingshot — cell_id merge via r_export cells.txt (same standard as LLR)
    sl_lw <- file.path(OUT, "slingshot", ds, "slingshot_weights.csv")
    if(file.exists(sl_lw)) {
      cw <- as.matrix(fread(sl_lw))
      cells_file <- file.path(COMPETITOR_R_EXPORT_ROOT, ds, "cells.txt")
      if(file.exists(cells_file) && nrow(cw) > 0L) {
        cells_vec <- readLines(cells_file, warn = FALSE)
        if(length(cells_vec) == nrow(cw)) {
          sl_dt <- data.table(cell_id = cells_vec, pred = apply(cw, 1, which.max))
          aligned <- merge(cl_dt, sl_dt, by = "cell_id", all = FALSE)
          if(nrow(aligned) >= 50L) {
            m2_all <- rbind(m2_all,
              data.table(dataset=ds, algorithm="slingshot", reference=cl$column,
                         lineage_cluster_source=NA_character_,
                         n_aligned=nrow(aligned),
                         ARI=ARI(aligned$pred, aligned$truth),
                         NMI=NMI(aligned$pred, aligned$truth),
                         Jaccard=Jaccard(aligned$pred, aligned$truth)))
          }
        }
      }
    }
    # LineaLinkR
    llr_rds <- file.path(LLR, ds, "linealinkr_complete.rds")
    if(LlrRdsTrusted(llr_rds)) {
      llr <- readRDS(llr_rds)
      if(!is.null(llr$lineage_weights)) {
        pred_dt <- data.table(
          cell_id = as.character(llr$cells),
          pred = apply(llr$lineage_weights, 1, which.max)
        )
        aligned <- merge(cl_dt, pred_dt, by = "cell_id", all = FALSE)
        if(nrow(aligned) >= 50L) {
        m2_all <- rbind(m2_all,
          data.table(dataset=ds, algorithm="LineaLinkR", reference=cl$column,
                     lineage_cluster_source=llr$lineage_cluster_source %||% NA_character_,
                     n_aligned=nrow(aligned),
                     ARI=ARI(aligned$pred, aligned$truth),
                     NMI=NMI(aligned$pred, aligned$truth),
                     Jaccard=Jaccard(aligned$pred, aligned$truth)))
        }
      }
    }
  }
}

if (INCREMENTAL_MODE) {
  cat("\n[stage] merge_incremental_pseudotime\n"); flush.console()
  retained <- BASE_M1[!algorithm %in% c(RECOMPUTE_ALGOS, "FateID")]
  m1_all <- rbindlist(list(retained, m1_all), use.names = TRUE, fill = TRUE)
}
m1_all[algorithm == "Palantir" & (is.na(cell_alignment_mode) | !nzchar(cell_alignment_mode)), cell_alignment_mode := "explicit_h5ad_obs_names"]
m1_all[(is.na(cell_alignment_mode) | !nzchar(cell_alignment_mode)), cell_alignment_mode := fcase(
  algorithm == "LineaLinkR", "explicit_rds_cells",
  algorithm %in% c("slingshot", "MultistageOT", "scFates"), "cells_txt_row_order_contract",
  algorithm %in% c("DPT", "SCORPIUS", "monocle3", "TSCAN", "StaVia", "scTour"), "explicit_cell_id",
  default = "unresolved"
)]
m1_all[, root_input_mode := unname(ROOT_INPUT_MODE[algorithm])]
m1_all[, direction_policy := ifelse(algorithm %in% NO_ROOT_ALGOS, "gt_aligned_if_negative", "signed")]
m1_all[, root_defined := !(algorithm %in% NO_ROOT_ALGOS)]

# ══════════════════════════════════════════════════════════════════
# Metric 3: Branch F1 (per-branch classification)
# ══════════════════════════════════════════════════════════════════
m3_all <- data.table()

BranchF1ManyToOne <- function(pred, truth) {
  # Label IDs from trajectory tools are arbitrary; score each true branch by
  # its best one-vs-rest predicted lineage instead of comparing raw IDs.
  truth_int <- as.integer(as.factor(truth))
  pred_labels <- sort(unique(pred))
  truth_labels <- sort(unique(truth_int))
  f1s <- sapply(truth_labels, function(tl) {
    max(sapply(pred_labels, function(pl) {
      tp <- sum(pred == pl & truth_int == tl)
      fp <- sum(pred == pl & truth_int != tl)
      fn <- sum(pred != pl & truth_int == tl)
      pre <- tp / (tp + fp)
      rec <- tp / (tp + fn)
      ifelse(pre + rec == 0, 0, 2 * pre * rec / (pre + rec))
    }))
  })
  mean(f1s)
}

BranchF1OneToOne <- function(pred, truth) {
  tbl <- table(as.factor(truth), as.factor(pred))
  precision <- sweep(tbl, 2L, colSums(tbl), "/")
  recall <- sweep(tbl, 1L, rowSums(tbl), "/")
  f1 <- 2 * precision * recall / (precision + recall)
  f1[!is.finite(f1)] <- 0
  n_match <- max(nrow(f1), ncol(f1))
  score <- matrix(0, nrow = n_match, ncol = n_match)
  score[seq_len(nrow(f1)), seq_len(ncol(f1))] <- f1
  assignment <- clue::solve_LSAP(score, maximum = TRUE)
  mean(score[cbind(seq_len(n_match), as.integer(assignment))])
}

for(ds in BRANCH_DATASETS) {
  cat(sprintf("[branch:%s]\n", ds)); flush.console()
  ref <- ExportRef(ds)
  if(is.null(ref)) next
  branch_threshold <- unname(BRANCH_TERMINAL_MIN_PSEUDOTIME[[ds]])
  branch_ref <- ref[
    is.finite(as.numeric(true_pseudotime)) &
      as.numeric(true_pseudotime) >= branch_threshold
  ]
  cl <- FindClusters(branch_ref)
  if(is.null(cl)) next
  cl_dt <- data.table(cell_id = as.character(cl$cell_id), truth = cl$values)
  
  # slingshot — cell_id merge via r_export cells.txt
  sl_lw <- file.path(OUT, "slingshot", ds, "slingshot_weights.csv")
  if(file.exists(sl_lw)) {
    cw <- as.matrix(fread(sl_lw))
    cells_file <- file.path(COMPETITOR_R_EXPORT_ROOT, ds, "cells.txt")
    if(file.exists(cells_file) && nrow(cw) > 0L) {
      cells_vec <- readLines(cells_file, warn = FALSE)
      if(length(cells_vec) == nrow(cw)) {
        sl_dt <- data.table(cell_id = cells_vec, pred = apply(cw, 1, which.max))
        aligned <- merge(cl_dt, sl_dt, by = "cell_id", all = FALSE)
        if(nrow(aligned) >= 50L) {
          m3_all <- rbind(m3_all,
            data.table(
              dataset=ds, algorithm="slingshot", lineage_cluster_source=NA_character_,
              n_aligned=nrow(aligned),
              Branch_F1=BranchF1OneToOne(aligned$pred, aligned$truth),
              Branch_F1_many_to_one_legacy=BranchF1ManyToOne(aligned$pred, aligned$truth)
            ))
        }
      }
    }
  }
  
  llr_rds <- file.path(LLR, ds, "linealinkr_complete.rds")
  if (LlrRdsTrusted(llr_rds)) {
    llr <- readRDS(llr_rds)
    if(!is.null(llr$lineage_weights)) {
      pred_dt <- data.table(
        cell_id = as.character(llr$cells),
        pred = apply(llr$lineage_weights, 1, which.max)
      )
      aligned <- merge(cl_dt, pred_dt, by = "cell_id", all = FALSE)
      if(nrow(aligned) >= 50L) {
       m3_all <- rbind(m3_all,
         data.table(dataset=ds, algorithm="LineaLinkR",
                    lineage_cluster_source=llr$lineage_cluster_source %||% NA_character_,
                    n_aligned=nrow(aligned),
                    Branch_F1=BranchF1OneToOne(aligned$pred, aligned$truth),
                    Branch_F1_many_to_one_legacy=BranchF1ManyToOne(aligned$pred, aligned$truth)))
      }
    }
  }
}
cat("[stage] branch_metrics_complete\n"); flush.console()
if (nrow(m3_all) > 0L) {
  fwrite(m3_all, file.path(EVAL, sprintf("table1_branch_f1_%s.csv", PROTOCOL_TAG)))
}

# ══════════════════════════════════════════════════════════════════
# Metric 4: Algorithm coverage (how many datasets each algorithm runs on)
# ══════════════════════════════════════════════════════════════════
m4_all <- data.table()
cat("[stage] coverage_start\n"); flush.console()
all_ds <- ORDERING_DATASETS
coverage_rows <- list()
coverage_index <- 0L
for(algo in c(PSEUDOTIME_ALGOS_INCLUDED, PSEUDOTIME_ALGOS_EXCLUDED)) {
  cat(sprintf("[coverage:%s]\n", algo)); flush.console()
  for(ds in all_ds) {
    coverage_index <- coverage_index + 1L
    ref <- ExportRef(ds)
    coverage_rows[[coverage_index]] <- if (is.null(ref)) {
      data.table(
        algorithm = algo,
        dataset = ds,
        has_artifact = HasPseudoOutput(algo, ds),
        n_reference = NA_integer_,
        n_output_finite = NA_integer_,
        n_aligned_finite = NA_integer_,
        n_missing_cells = NA_integer_,
        n_unexpected_cells = NA_integer_,
        cell_coverage = NA_real_,
        coverage_status = "reference_unavailable"
      )
    } else {
      AuditPseudoCoverage(algo, ds, ref)
    }
  }
}
coverage_all <- rbindlist(coverage_rows)
m4_all <- coverage_all[algorithm %in% PSEUDOTIME_ALGOS_INCLUDED, .(
  datasets_with_artifact = sum(has_artifact),
  datasets_complete_output = sum(coverage_status == "complete_output"),
  datasets_partial_output = sum(coverage_status == "partial_output"),
  datasets_unreadable_output = sum(
    coverage_status %in% c(
      "unreadable_output",
      "duplicate_cell_ids",
      "unexpected_cell_ids",
      "reference_unavailable"
    )
  ),
  median_cell_coverage = median(cell_coverage, na.rm = TRUE),
  minimum_cell_coverage = min(cell_coverage, na.rm = TRUE)
), by = algorithm][order(-datasets_complete_output, -datasets_with_artifact)]
cat("[stage] coverage_complete\n"); flush.console()

# ══════════════════════════════════════════════════════════════════
# Gap matrix: planned TI algorithms × datasets
# ══════════════════════════════════════════════════════════════════
cat("[stage] gap_matrix_start\n"); flush.console()
gt_cache <- setNames(vector("list", length(all_ds)), all_ds)
for (ds in all_ds) {
  ref <- ExportRef(ds)
  gt_cache[[ds]] <- if (!is.null(ref)) FindGroundTruth(ref, ds) else NULL
}
gap_rows <- list()
idx <- 0L
for (algo in c(PSEUDOTIME_ALGOS_INCLUDED, PSEUDOTIME_ALGOS_EXCLUDED)) {
  for (ds in all_ds) {
    idx <- idx + 1L
    coverage_row <- coverage_all[algorithm == algo & dataset == ds][1L]
    has_out <- coverage_row$has_artifact
    has_gt <- !is.null(gt_cache[[ds]])
    ev <- m1_all[algorithm == algo & dataset == ds]
    gap_rows[[idx]] <- data.table(
      algorithm = algo,
      dataset = ds,
      in_included_panel = algo %in% PSEUDOTIME_ALGOS_INCLUDED,
      planned_ti = algo %in% TI_PLANNED,
      has_output = has_out,
      coverage_status = coverage_row$coverage_status,
      n_reference = coverage_row$n_reference,
      n_output_finite = coverage_row$n_output_finite,
      n_aligned_finite = coverage_row$n_aligned_finite,
      cell_coverage = coverage_row$cell_coverage,
      has_numeric_gt = has_gt,
      spearman_evaluated = nrow(ev) > 0L && is.finite(ev$Spearman[1]),
      spearman = if (nrow(ev)) ev$Spearman[1] else NA_real_,
       gap_type = fifelse(
        algo %in% PSEUDOTIME_ALGOS_EXCLUDED, "excluded_non_pseudotime_or_proxy",
        fifelse(!has_out, "missing_output",
          fifelse(coverage_row$coverage_status != "complete_output",
            coverage_row$coverage_status,
            fifelse(!has_gt, "output_no_gt",
              fifelse(nrow(ev) == 0L, "evaluation_failed", "evaluated")))))
    )
  }
}
gap_all <- rbindlist(gap_rows)
for (algo in setdiff(TI_PLANNED, c(PSEUDOTIME_ALGOS_INCLUDED, PSEUDOTIME_ALGOS_EXCLUDED))) {
  for (ds in all_ds) {
    gap_all <- rbind(gap_all, data.table(
      algorithm = algo, dataset = ds, in_included_panel = FALSE, planned_ti = TRUE,
      has_output = FALSE, coverage_status = "algorithm_not_run",
      n_reference = NA_integer_, n_output_finite = 0L,
      n_aligned_finite = 0L, cell_coverage = NA_real_,
      has_numeric_gt = FALSE, spearman_evaluated = FALSE,
      spearman = NA_real_, gap_type = "algorithm_not_run"
    ))
  }
}

gap_summary <- rbindlist(list(
  gap_all[in_included_panel == TRUE, .(
    n_datasets = .N,
    n_with_output = sum(has_output),
    n_with_gt = sum(has_numeric_gt),
    n_evaluated = sum(spearman_evaluated),
    n_complete_output = sum(coverage_status == "complete_output"),
    n_partial_output = sum(coverage_status == "partial_output"),
    n_unreadable_output = sum(
      coverage_status %in% c("unreadable_output", "duplicate_cell_ids", "unexpected_cell_ids")
    ),
    n_missing_output = sum(gap_type == "missing_output"),
    n_output_no_gt = sum(gap_type == "output_no_gt"),
    n_evaluation_failed = sum(gap_type == "evaluation_failed")
  ), by = algorithm],
  gap_all[planned_ti == TRUE & in_included_panel == FALSE & algorithm %in% PSEUDOTIME_ALGOS_EXCLUDED, .(
    n_datasets = .N, n_with_output = sum(has_output), n_with_gt = sum(has_numeric_gt),
    n_evaluated = sum(spearman_evaluated),
    n_complete_output = sum(coverage_status == "complete_output"),
    n_partial_output = sum(coverage_status == "partial_output"),
    n_unreadable_output = sum(
      coverage_status %in% c("unreadable_output", "duplicate_cell_ids", "unexpected_cell_ids")
    ),
    n_missing_output = sum(gap_type == "missing_output"),
    n_output_no_gt = sum(gap_type == "output_no_gt"),
    n_evaluation_failed = sum(gap_type == "evaluation_failed")
  ), by = algorithm],
  gap_all[gap_type == "algorithm_not_run", .(
    n_datasets = .N, n_with_output = 0L, n_with_gt = 0L, n_evaluated = 0L,
    n_complete_output = 0L, n_partial_output = 0L,
    n_unreadable_output = 0L, n_missing_output = .N,
    n_output_no_gt = 0L, n_evaluation_failed = 0L
  ), by = algorithm]
), fill = TRUE)
gap_summary <- unique(gap_summary, by = "algorithm")

planned_not_run <- setdiff(TI_PLANNED, c(PSEUDOTIME_ALGOS_INCLUDED, PSEUDOTIME_ALGOS_EXCLUDED))
cat("[stage] gap_matrix_complete\n"); flush.console()

# ══════════════════════════════════════════════════════════════════
# Results
# ══════════════════════════════════════════════════════════════════
cat("[stage] result_tables_start\n"); flush.console()

cat("\n\n========================================================\n")
cat("  MULTI-DATASET BENCHMARK RESULTS\n")
cat("========================================================\n\n")

if(nrow(m1_all) > 0) {
  cat("Metric 1 — Pseudotime Spearman ρ\n")
  cat("(ordered by algorithm, all datasets):\n\n")
  m1_cast <- dcast(m1_all, algorithm ~ dataset, value.var="Spearman")
  print(m1_cast, digits=3)
}

if(nrow(m2_all) > 0) {
  cat("\nMetric 2 — Lineage ARI & NMI\n\n")
  m2a <- dcast(m2_all, algorithm ~ dataset, value.var="ARI")
  m2n <- dcast(m2_all, algorithm ~ dataset, value.var="NMI")
  cat("ARI:\n"); print(m2a, digits=3)
  cat("\nNMI:\n"); print(m2n, digits=3)
}

  # Mark duplicate datasets (keep first in each set, flag rest)
  dup_flag <- function(dt, dup_sets) {
    dt[, is_duplicate := FALSE]
    for(dset in dup_sets) {
      keep <- dset[1]
      drop <- dset[-1]
      dt[dataset %in% drop, is_duplicate := TRUE]
    }
    dt
  }
  m1_all <- dup_flag(m1_all, DUPLICATE_SETS)
  m2_all <- dup_flag(m2_all, DUPLICATE_SETS)

  stopifnot(
    "The fixed complete-core matrix has no evaluable pseudotime rows" = all(c(
      "dataset", "algorithm", "Spearman", "n_reference", "n_aligned",
      "cell_coverage", "is_duplicate"
    ) %in% names(m1_all))
  )
  m1_dedup <- m1_all[is_duplicate == FALSE]
  core_complete <- m1_dedup[
    dataset %in% CORE_DATASETS & algorithm %in% PSEUDOTIME_ALGOS_INCLUDED
  ]
  core_key_counts <- core_complete[, .N, by = .(dataset, algorithm)]
  stopifnot(
    "The fixed complete-core matrix is not one complete row per configured dataset and algorithm" =
      nrow(core_complete) == length(CORE_DATASETS) * length(PSEUDOTIME_ALGOS_INCLUDED) &&
      nrow(core_key_counts) == length(CORE_DATASETS) * length(PSEUDOTIME_ALGOS_INCLUDED) &&
      all(core_key_counts$N == 1L) &&
      setequal(core_complete$dataset, CORE_DATASETS) &&
      setequal(core_complete$algorithm, PSEUDOTIME_ALGOS_INCLUDED) &&
      all(is.finite(core_complete$Spearman)) &&
      all(core_complete$n_reference == core_complete$n_aligned) &&
      all(core_complete$cell_coverage == MIN_PSEUDOTIME_CELL_COVERAGE)
  )
  core_complete[, dataset_rank := frank(-Spearman, ties.method = "average"), by = dataset]
  core_summary <- core_complete[, .(
    n_common = .N,
    mean_spearman = mean(Spearman),
    median_spearman = median(Spearman),
    mean_rank = mean(dataset_rank)
  ), by = .(algorithm, root_input_mode, direction_policy)][order(-mean_spearman)]
  friedman_result <- friedman.test(Spearman ~ algorithm | dataset, data = core_complete)
  friedman_summary <- data.table(
    n_datasets = length(CORE_DATASETS),
    n_algorithms = length(PSEUDOTIME_ALGOS_INCLUDED),
    statistic = unname(friedman_result$statistic),
    df = unname(friedman_result$parameter),
    p_value = friedman_result$p.value
  )
  fwrite(core_complete, file.path(EVAL, sprintf("table1_common_core_values_%s.csv", PROTOCOL_TAG)))
  fwrite(core_summary, file.path(EVAL, sprintf("table1_common_core_summary_%s.csv", PROTOCOL_TAG)))
  fwrite(friedman_summary, file.path(EVAL, sprintf("table1_common_core_friedman_%s.csv", PROTOCOL_TAG)))
  
  # Deduplicated summary (for paper: excludes duplicate datasets)
  if(nrow(m1_all) > 0) {
    extension_values <- m1_dedup[
      dataset %in% ORDERING_DATASETS &
        algorithm %in% PSEUDOTIME_ALGOS_INCLUDED &
        is.finite(Spearman)
    ]
    llr_pair <- extension_values[
      algorithm == "LineaLinkR", .(dataset, linealinkr_spearman = Spearman)
    ]
    pairwise <- merge(
      extension_values[
        algorithm != "LineaLinkR",
        .(dataset, algorithm, competitor_spearman = Spearman)
      ],
      llr_pair,
      by = "dataset",
      all = FALSE
    )
    pairwise[, delta_linealinkr_minus_competitor := linealinkr_spearman - competitor_spearman]
    pairwise_summary <- pairwise[, .(
      n_planned_extension_datasets = length(ORDERING_DATASETS),
      n_common = .N,
      linealinkr_mean = mean(linealinkr_spearman),
      competitor_mean = mean(competitor_spearman),
      mean_delta = mean(delta_linealinkr_minus_competitor),
      median_delta = median(delta_linealinkr_minus_competitor),
      linealinkr_wins = sum(delta_linealinkr_minus_competitor > 0),
      ties = sum(delta_linealinkr_minus_competitor == 0),
      linealinkr_losses = sum(delta_linealinkr_minus_competitor < 0)
    ), by = algorithm]
    pairwise_summary[, algorithm_order := match(algorithm, PSEUDOTIME_ALGOS_INCLUDED)]
    setorder(pairwise_summary, algorithm_order)
    pairwise_summary[, algorithm_order := NULL]
    fwrite(pairwise, file.path(EVAL, sprintf("table1_pairwise_values_%s.csv", PROTOCOL_TAG)))
    fwrite(pairwise_summary, file.path(EVAL, sprintf("table1_pairwise_summary_%s.csv", PROTOCOL_TAG)))
  }
  
  m1_all[, panel := "included_official"]
  fwrite(m1_all, file.path(EVAL, sprintf("table1_pseudotime_%s.csv", PROTOCOL_TAG)))
  fwrite(
    m1_all[, .(
      dataset, algorithm, Spearman, Spearman_signed, orientation_flipped,
      root_input_mode, direction_policy, cell_alignment_mode,
      n_reference, n_aligned, cell_coverage
    )],
    file.path(SUPP, "tableS_pseudotime_signed_orientation.csv")
  )

  # The manuscript summary is restricted to one identical dataset set for every
  # headline method. Unequal marginal means remain excluded from paper tables.
  summary_all <- core_complete[, .(
    mean_spearman = mean(Spearman),
    n_spearman = .N
  ), by = algorithm]
  summary_scope_long <- core_complete[, .(
    mean_spearman_scope = mean(Spearman),
    n_spearman_scope = .N
  ), by = .(algorithm, scope = fifelse(grepl("^sim", dataset), "sim", "real"))]
  summary_scope <- dcast(
    summary_scope_long,
    algorithm ~ scope,
    value.var = c("mean_spearman_scope", "n_spearman_scope")
  )
  setnames(
    summary_scope,
    c("mean_spearman_scope_real", "mean_spearman_scope_sim", "n_spearman_scope_real", "n_spearman_scope_sim"),
    c("mean_spearman_real", "mean_spearman_sim", "n_spearman_real", "n_spearman_sim"),
    skip_absent = TRUE
  )
  summary_lineage <- m2_all[is_duplicate == FALSE, .(
    mean_lineage_ari = mean(ARI),
    n_lineage_ari = .N
  ), by = algorithm]
  summary_branch <- m3_all[, .(
    mean_branch_f1 = mean(Branch_F1),
    n_branch_f1 = .N,
    mean_branch_f1_many_to_one_legacy = mean(Branch_F1_many_to_one_legacy)
  ), by = algorithm]
  table1_summary <- merge(summary_all, summary_scope, by = "algorithm", all = TRUE)
  table1_summary <- merge(table1_summary, summary_lineage, by = "algorithm", all = TRUE)
  table1_summary <- merge(table1_summary, summary_branch, by = "algorithm", all = TRUE)
  setorder(table1_summary, -mean_spearman)
  fwrite(table1_summary, file.path(CONFIG$dir$assets, "manuscript", "table1_benchmark_summary.csv"))
  fwrite(m2_all, file.path(EVAL, sprintf("table1_lineage_%s.csv", PROTOCOL_TAG)))
  if (nrow(m3_all) > 0) {
    cat("\nMetric 3 — Branch F1 (mean per-branch)\n\n")
    m3_cast <- dcast(m3_all, algorithm ~ dataset, value.var = "Branch_F1")
    print(m3_cast, digits = 3)
    fwrite(m3_all, file.path(EVAL, sprintf("table1_branch_f1_%s.csv", PROTOCOL_TAG)))
  }
  if (nrow(m4_all) > 0) {
    cat("\nMetric 4 — Algorithm Coverage (datasets with pseudotime output)\n\n")
    print(m4_all)
    fwrite(m4_all, file.path(EVAL, "tableS1_algorithm_coverage.csv"))
  }
  fwrite(gap_all, file.path(SUPP, "benchmark_gap_matrix.csv"))
  fwrite(gap_summary, file.path(SUPP, "benchmark_gap_summary.csv"))
  fwrite(coverage_all, file.path(SUPP, "benchmark_output_coverage.csv"))
  reg_dir <- file.path(CONFIG$dir$paper_tables, "_registry")
  dir.create(reg_dir, recursive = TRUE, showWarnings = FALSE)
  fwrite(gap_summary, file.path(reg_dir, "benchmark_gap_summary.csv"))
  fwrite(coverage_all, file.path(reg_dir, "benchmark_output_coverage.csv"))

  discovered <- rbindlist(lapply(PSEUDOTIME_ALGOS_INCLUDED, function(a) {
    rbindlist(lapply(all_ds, function(d) {
      p <- PseudoPath(a, d)
      data.table(
        algorithm = a,
        dataset = d,
        has_artifact = HasPseudoOutput(a, d),
        source_path = if (length(p) && !is.na(p) && file.exists(p)) p else NA_character_
      )
    }))
  }))
  fwrite(discovered, file.path(reg_dir, "benchmark_discovered_outputs.csv"))

  cat("\n── Benchmark gap (included TI panel) ──\n")
  print(gap_summary[algorithm %in% PSEUDOTIME_ALGOS_INCLUDED][order(algorithm)], digits = 3)
cat(sprintf("\nSaved: %s\n", EVAL))
cat(sprintf("Gap tables: %s\n", SUPP))
