suppressPackageStartupMessages({library(Matrix);library(data.table);library(CellChat)})

# path bootstrap (relative; set LLR_ASSETS or LLR_PIPELINE to override)
.llr_cmd <- commandArgs(trailingOnly = FALSE)
.llr_file <- sub("^--file=", "", grep("^--file=", .llr_cmd, value = TRUE)[1])
if (!length(.llr_file) || is.na(.llr_file) || !nzchar(.llr_file)) {
  stop("Run this script via Rscript so --file= is set")
}
.llr_here <- dirname(normalizePath(.llr_file, winslash = "/", mustWork = TRUE))
source(file.path(.llr_here, "../../config/llr_bootstrap.R"), local = FALSE)

llr_source_config()
.libPaths(c(
  Sys.getenv("LLR064_RLIB", unset = file.path(CONFIG$dir$root, "rlib")),
  .libPaths()
))
`%||%` <- function(x, y) if (is.null(x)) y else x
llr_cfg <- CONFIG$llr
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

NormalizeExtractMeta <- function(meta, cells) {
  if (nrow(meta) > 0L && any(as.character(meta[[2]]) == "paul15_clusters", na.rm = TRUE)) {
    hdr <- as.character(meta[1])
    colnames(meta) <- if (trimws(hdr[1]) == "") c("barcode_id", hdr[-1]) else hdr
    meta <- meta[-1]
  }
  idx_col <- if ("barcode_id" %in% names(meta)) "barcode_id" else names(meta)[1]
  if (length(cells) > 0L && idx_col %in% names(meta)) meta <- meta[match(cells, as.character(meta[[idx_col]]))]
  if (!"clusters" %in% names(meta) && "paul15_clusters" %in% names(meta)) meta$clusters <- meta$paul15_clusters
  meta
}

library(LineaLinkR)
stopifnot(identical(
  as.character(packageVersion("LineaLinkR")),
  CONFIG$llr$package_version
))

set.seed(CONFIG$seed)
Sys.setenv(OMP_NUM_THREADS = CONFIG$n_threads)
data.table::setDTthreads(CONFIG$n_threads)

args <- commandArgs(trailingOnly = TRUE)
ds_name <- args[1]
if (is.na(ds_name)) stop("Usage: Rscript run_linealinkr_ablation.R <dataset>")

ABL_DIR <- file.path(CONFIG$dir$cache, "linealinkr_ablation")
dir.create(ABL_DIR, recursive = TRUE, showWarnings = FALSE)

sp <- CONFIG$dataset_species[[ds_name]]
if (is.null(sp)) sp <- "human"

baseline_rds <- file.path(CONFIG$dir$linealinkr, ds_name, "linealinkr_complete.rds")
if (!file.exists(baseline_rds)) stop("Baseline not found. Run 02_run_linealinkr.R first.")
full <- readRDS(baseline_rds)
stopifnot(
  identical(as.character(full$dataset), ds_name),
  identical(as.character(full$package_version), CONFIG$llr$package_version),
  identical(as.character(full$benchmark_protocol), CONFIG$llr$benchmark_protocol),
  identical(as.character(full$result_revision), CONFIG$llr$result_revision),
  identical(full$input_manifest$content_schema, "decoded_h5ad_content_v1"),
  identical(
    normalizePath(full$input_manifest$canonical_path, winslash = "/", mustWork = TRUE),
    normalizePath(DatasetPath(ds_name), winslash = "/", mustWork = TRUE)
  )
)

extract <- file.path(CONFIG$dir$linealinkr, ds_name, "extract")
matrix_mtx <- file.path(extract, sprintf("%s_matrix.mtx", ds_name))
meta_csv <- file.path(extract, sprintf("%s_meta.csv", ds_name))
genes_txt <- file.path(extract, sprintf("%s_genes.txt", ds_name))
cells_txt <- file.path(extract, sprintf("%s_cells.txt", ds_name))

expr <- as(readMM(matrix_mtx), "CsparseMatrix")
meta <- fread(meta_csv, header = TRUE)
genes <- readLines(genes_txt)
cells <- readLines(cells_txt)
stopifnot(identical(as.character(full$cells), as.character(cells)))
rownames(expr) <- genes
colnames(expr) <- cells
expr_lr <- expr %*% Matrix::Diagonal(x = as.numeric(llr_cfg$lr_scale_factor %||% 10000) / Matrix::colSums(expr))
expr_lr@x <- log1p(expr_lr@x)

meta <- NormalizeExtractMeta(meta, cells)
ct_col <- as.character(full$celltype_column)
if (identical(ct_col, "auto_clusters")) {
  auto_clusters <- fread(file.path(
    dirname(baseline_rds),
    sprintf("%s_auto_clusters.csv", ds_name)
  ))
  stopifnot(nrow(auto_clusters) == length(cells))
  meta$auto_clusters <- as.character(auto_clusters$auto_clusters)
}
stopifnot(
  length(ct_col) == 1L,
  nzchar(ct_col),
  ct_col %in% names(meta)
)
ct_values <- as.character(meta[[ct_col]])
ct_values[is.na(ct_values) | trimws(ct_values) == ""] <- "unassigned"
celltypes <- as.integer(as.factor(ct_values))
n_types <- max(celltypes, na.rm = TRUE)
n_cells <- ncol(expr)

fine_col <- if (identical(ct_col, "auto_clusters")) {
  ct_col
} else if ("clusters_fine" %in% names(meta)) {
  "clusters_fine"
} else if ("clusters_coarse" %in% names(meta)) {
  "clusters_coarse"
} else {
  ct_col
}
fine_values <- as.character(meta[[fine_col]])
fine_values[is.na(fine_values) | trimws(fine_values) == ""] <- "unassigned"
lineage_cl_base <- as.integer(as.factor(fine_values))

# start_cell is selected inside ComputeUpstream after UMAP.
start_cell <- 1L

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
    ligand_rows[[length(ligand_rows) + 1L]] <- lig_idx
    receptor_rows[[length(receptor_rows) + 1L]] <- rec_idx
  }
  list(ligand_rows = ligand_rows, receptor_rows = receptor_rows)
}

BuildPairGeneMap <- function(expr, species = "human") {
  if (!species %in% c("human", "mouse")) {
    return(data.table(pair_idx = integer(), ligand_genes = character(), receptor_genes = character()))
  }
  cc_db <- if (species == "mouse") CellChat::CellChatDB.mouse else CellChat::CellChatDB.human
  lr_pairs <- cc_db$interaction
  gene_names <- rownames(expr)
  pair_genes <- data.table(pair_idx = integer(), ligand_genes = character(), receptor_genes = character())
  pid <- 0L
  for (i in seq_len(nrow(lr_pairs))) {
    lig <- strsplit(lr_pairs$ligand[i], "_")[[1]]
    rec <- strsplit(lr_pairs$receptor[i], "_")[[1]]
    if (any(is.na(match(lig, gene_names))) || any(is.na(match(rec, gene_names)))) next
    pid <- pid + 1L
    pair_genes <- rbind(pair_genes, data.table(pair_idx = pid, ligand_genes = lr_pairs$ligand[i], receptor_genes = lr_pairs$receptor[i]))
  }
  pair_genes
}

RunScanpyKNN <- function(matrix_path, out_dir, n_neighbors) {
  py <- sprintf(r"(
import scanpy as sc, numpy as np
from scipy.io import mmread, mmwrite
from scipy.sparse import csr_matrix
np.random.seed(%d)
sc.settings.seed = %d
a = sc.AnnData(mmread(r'%s').tocsr().T)
sc.pp.normalize_total(a)
sc.pp.log1p(a)
sc.pp.highly_variable_genes(a, n_top_genes=2000)
a = a[:, a.var.highly_variable].copy()
sc.pp.pca(a, n_comps=30)
sc.pp.neighbors(a, n_neighbors=%d)
pca30 = a.obsm['X_pca'][:, :min(30, a.obsm['X_pca'].shape[1])]
mmwrite(r'%s', csr_matrix(a.obsp['distances']))
mmwrite(r'%s', csr_matrix(a.obsp['connectivities']))
np.savetxt(r'%s', pca30, delimiter=',')
np.savetxt(r'%s', pca30[:, :2], delimiter=',')
print('KNN_OK')
)", CONFIG$seed, CONFIG$seed, matrix_path, n_neighbors,
    file.path(out_dir, "knn_dist.mtx"),
    file.path(out_dir, "knn_conn.mtx"),
    file.path(out_dir, "pca.csv"),
    file.path(out_dir, "umap.csv"))
  tf <- file.path(out_dir, sprintf("_scanpy_knn_k%d.py", n_neighbors))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  writeLines(py, tf, useBytes = TRUE)
  out <- system2(CONFIG$exe$python, args = tf, stdout = TRUE, stderr = TRUE)
  unlink(tf)
  if (!any(grepl("KNN_OK", out))) stop(paste(out, collapse = "\n"))
}

ComputeUpstream <- function(
  out_dir,
  n_neighbors = 30L,
  reuse_baseline_graph = FALSE,
  root_clusters = NULL,
  root_cluster_source = NULL,
  ablation_id = "graph_k"
) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  baseline_k <- identical(
    as.integer(n_neighbors),
    as.integer(llr_cfg$knn_neighbors %||% 50L)
  )
  if (baseline_k && !reuse_baseline_graph && is.null(root_clusters)) {
    return(list(
      dataset = ds_name,
      package_version = CONFIG$llr$package_version,
      benchmark_protocol = CONFIG$llr$benchmark_protocol,
      result_revision = CONFIG$llr$result_revision,
      execution_status = "baseline_reuse",
      ablation_id = "baseline_reuse",
      n_cells = n_cells,
      n_genes = nrow(expr),
      n_types = n_types,
      pseudotime = full$pseudotime,
      lineage_weights = full$lineage_weights,
      pseudotime_per_lineage = full$pseudotime_per_lineage,
      n_lineages = full$n_lineages,
      start_cell = full$start_cell,
      root_cluster_source = full$root_cluster_source %||% NA_character_,
      root_clusters = full$root_clusters,
      lineage_cluster_source = full$lineage_cluster_source %||% NA_character_,
      lineage_clusters = full$lineage_clusters,
      root_strategy = llr_cfg$root_strategy %||% "attractor_energy",
      root_potency_mode = llr_cfg$root_potency_mode %||%
        "detected_genes_ordering",
      selected_root_potency_mode =
        full$selected_root_potency_mode %||% NA_character_,
      root_selection = full$root_selection %||% NULL,
      endpoint_selection =
        full$root_diagnostics$endpoint_selection %||% NULL,
      root_support = full$root_support %||% NA_real_,
      direction_margin = full$direction_margin %||% NA_real_,
      direction_status = full$direction_status %||% NA_character_
    ))
  }
  graph_dir <- if (reuse_baseline_graph) dirname(baseline_rds) else out_dir
  if (!reuse_baseline_graph) RunScanpyKNN(matrix_mtx, graph_dir, n_neighbors)
  knn_dist <- as(readMM(file.path(graph_dir, "knn_dist.mtx")), "CsparseMatrix")
  knn_conn <- as(readMM(file.path(graph_dir, "knn_conn.mtx")), "CsparseMatrix")
  # Reconstruct nn_idx/nn_dists because this helper is not part of the exported API.
  sm <- summary(knn_dist)
  sm_i <- sm$i; sm_j <- sm$j; sm_x <- sm$x
  ord <- order(sm_i, sm_x, decreasing = c(FALSE, FALSE))
  sm_i <- sm_i[ord]; sm_j <- sm_j[ord]; sm_x <- sm_x[ord]
  run_starts <- c(0L, which(diff(sm_i) != 0L), length(sm_i))
  nn_idx <- matrix(0L, nrow = n_cells, ncol = n_neighbors)
  nn_dists <- matrix(0.0, nrow = n_cells, ncol = n_neighbors)
  for (g in seq_len(length(run_starts) - 1L)) {
    i <- sm_i[run_starts[g] + 1L]
    s <- run_starts[g] + 1L
    e <- run_starts[g + 1L]
    take <- seq_len(min(n_neighbors, e - s + 1L))
    nn_idx[i, take] <- sm_j[s + take - 1L]
    nn_dists[i, take] <- sm_x[s + take - 1L]
  }
  neigh <- list(nn_idx = nn_idx, nn_dists = nn_dists)
  blind_emb <- as.matrix(fread(file.path(graph_dir, "pca.csv"), header = FALSE))
  geom_emb <- blind_emb[, seq_len(min(2L, ncol(blind_emb))), drop = FALSE]
  if (isTRUE(llr_cfg$connect_components)) {
    cc <- LineaLinkR::ConnectKnnComponents(neigh$nn_idx, neigh$nn_dists, geom_emb)
    neigh$nn_idx <- cc$nn_idx
    neigh$nn_dists <- cc$nn_dists
  }
  root_cl <- as.integer(root_clusters %||% full$root_clusters)
  lineage_cl <- as.integer(full$lineage_clusters)
  stopifnot(length(root_cl) == n_cells, length(lineage_cl) == n_cells)
  root_cluster_source <- root_cluster_source %||% paste0(
    full$root_cluster_source %||% NA_character_,
    ":fixed_for_graph_k"
  )
  lineage_cluster_source <- paste0(
    full$lineage_cluster_source %||% NA_character_,
    ":fixed_for_ablation"
  )
  root_potency_mode <-
    llr_cfg$root_potency_mode %||% "detected_genes_ordering"
  detected_genes <- as.numeric(Matrix::colSums(expr > 0))
  high_potency <- log1p(detected_genes)
  low_potency <- -high_potency
  high_diagnostics <- LineaLinkR::InferStartCellAttractor(
    blind_emb, root_cl, neigh$nn_idx, neigh$nn_dists,
    potency = high_potency,
    top_fraction = llr_cfg$root_top_fraction %||% 0.05,
    n_attractors = llr_cfg$root_n_attractors %||% 0L,
    max_iter = llr_cfg$root_max_iter %||% 300L,
    tol = llr_cfg$root_tol %||% 1e-5,
    weights = llr_cfg$root_weights %||% c(0.25, 0.20, 0.30, 0.15, 0.10),
    seed = CONFIG$seed
  )
  low_diagnostics <- LineaLinkR::InferStartCellAttractor(
    blind_emb, root_cl, neigh$nn_idx, neigh$nn_dists,
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
    neigh$nn_idx,
    neigh$nn_dists,
    candidate_anchor_cells,
    rank_normalize = isTRUE(llr_cfg$rank_normalize)
  )
  candidate_pt <- cbind(
    LineaLinkR::OrientPseudotimeByPotency(
      candidate_pt_raw[, 1L], high_potency,
      llr_cfg$root_direction_min_potency_cor %||% 0.10
    ),
    LineaLinkR::OrientPseudotimeByPotency(
      candidate_pt_raw[, 2L], high_potency,
      llr_cfg$root_direction_min_potency_cor %||% 0.10
    ),
    LineaLinkR::OrientPseudotimeByPotency(
      candidate_pt_raw[, 3L], high_potency,
      llr_cfg$root_direction_min_potency_cor %||% 0.10
    ),
    LineaLinkR::OrientPseudotimeByPotency(
      candidate_pt_raw[, 4L], high_potency,
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
    potency_spearman = c(
      stats::cor(candidate_pt[, 1L], high_potency, method = "spearman"),
      stats::cor(candidate_pt[, 2L], high_potency, method = "spearman"),
      stats::cor(candidate_pt[, 3L], high_potency, method = "spearman"),
      stats::cor(candidate_pt[, 4L], high_potency, method = "spearman")
    ),
    finite_fraction = colMeans(is.finite(candidate_pt)),
    ordering_score = c(
      LineaLinkR::BlindOrderingScore(candidate_pt[, 1L], blind_emb),
      LineaLinkR::BlindOrderingScore(candidate_pt[, 2L], blind_emb),
      LineaLinkR::BlindOrderingScore(candidate_pt[, 3L], blind_emb),
      LineaLinkR::BlindOrderingScore(candidate_pt[, 4L], blind_emb)
    ),
    component_score = c(
      as.numeric(high_component_score),
      as.numeric(high_diagnostics$cellwise_component_score),
      as.numeric(low_component_score),
      as.numeric(low_diagnostics$cellwise_component_score)
    ),
    tie_rank = 1:4
  )
  endpoint_selection <- LineaLinkR::SelectRootEndpointByPotencyAlignment(
    root_selection,
    default_mode = "detected_genes_high",
    min_alignment_gain =
      llr_cfg$root_endpoint_alignment_min_gain %||% 0.20,
    min_ordering_gain = llr_cfg$root_ordering_min_gain %||% 0.10,
    min_representation_alignment_gain =
      llr_cfg$root_representation_alignment_min_gain %||% 0.00
  )
  selected_index <- match(
    endpoint_selection$selected_candidate,
    root_selection$candidate
  )
  pseudotime <- candidate_pt[, selected_index]
  oriented_root_cells <- which(
    pseudotime == min(pseudotime, na.rm = TRUE)
  )
  start_cell <- oriented_root_cells[[order(
    -high_potency[oriented_root_cells],
    oriented_root_cells
  )[1L]]]
  root_start_cluster <- as.integer(root_cl[start_cell])
  if (isTRUE(llr_cfg$augment_mst_bridges)) {
    br <- LineaLinkR::AugmentKnnWithMstBridges(neigh$nn_idx, neigh$nn_dists, geom_emb, root_cl, root_start_cluster)
    neigh$nn_idx <- br$nn_idx
    neigh$nn_dists <- br$nn_dists
  }
  lineage_start_cluster <- as.integer(lineage_cl[start_cell])
  lineage_weights <- LineaLinkR::ComputeLineageWeightsTrunkAware(
    geom_emb,
    lineage_cl,
    lineage_start_cluster
  )
  n_lineages <- ncol(lineage_weights)
  # This ablation uses the current global-pseudotime API.
  pseudotime_per_lineage <- NULL
  direction_diagnostics <- LineaLinkR::AuditPseudotimeDirection(
    pseudotime,
    potency = high_potency,
    embedding = blind_emb
  )
  list(
    dataset = ds_name,
    package_version = CONFIG$llr$package_version,
    benchmark_protocol = CONFIG$llr$benchmark_protocol,
    result_revision = CONFIG$llr$result_revision,
    execution_status = "pipeline_rerun_done",
    ablation_id = ablation_id,
    n_cells = n_cells,
    n_genes = nrow(expr),
    n_types = n_types,
    pseudotime = pseudotime,
    lineage_weights = lineage_weights,
    pseudotime_per_lineage = pseudotime_per_lineage,
    n_lineages = n_lineages,
    start_cell = start_cell,
    root_cluster_source = root_cluster_source,
    root_clusters = root_cl,
    lineage_cluster_source = lineage_cluster_source,
    lineage_clusters = lineage_cl,
    root_strategy = llr_cfg$root_strategy %||% "attractor_energy",
    root_potency_mode = root_potency_mode,
    selected_root_potency_mode = endpoint_selection$selected_candidate,
    root_selection = root_selection,
    endpoint_selection = endpoint_selection,
    direction_diagnostics = direction_diagnostics,
    root_support = NA_real_,
    direction_margin = abs(
      direction_diagnostics$blind_score_margin %||% NA_real_
    ),
    direction_status = direction_diagnostics$direction_status
  )
}

SummarizeObject <- function(obj) {
  lr_n <- if (is.null(obj$lr_scores)) 0L else nrow(obj$lr_scores)
  lr_sig <- if (is.null(obj$lr_scores)) 0L else sum(obj$lr_scores$pvalue < 0.05, na.rm = TRUE)
  temporal_n <- if (is.null(obj$temporal_lr) || is.null(obj$temporal_lr$scores)) 0L else nrow(obj$temporal_lr$scores)
  temporal_pairs <- if (is.null(obj$temporal_lr) || is.null(obj$temporal_lr$scores) || temporal_n == 0L) 0L else uniqueN(obj$temporal_lr$scores$pair_idx)
  bif_n <- if (is.null(obj$bif_df)) 0L else sum(obj$bif_df$reason == "crossing", na.rm = TRUE)
  data.table(
    n_lr_pairs = lr_n,
    n_lr_sig = lr_sig,
    n_temporal_rows = temporal_n,
    n_temporal_pairs = temporal_pairs,
    n_bifurcation_points = bif_n,
    n_lineages = obj$n_lineages,
    has_perlineage_pt = !is.null(obj$pseudotime_per_lineage)
  )
}

RunModuleAblation <- function(ablation_id) {
  module_dir <- file.path(ABL_DIR, ds_name, ablation_id)
  upstream <- ComputeUpstream(module_dir, as.integer(llr_cfg$knn_neighbors %||% 50L))
  if (ablation_id == "no_lineage_specific_scoring") {
    upstream$lineage_weights <- matrix(
      1.0,
      nrow = n_cells,
      ncol = 1L,
      dimnames = list(cells, "global")
    )
    upstream$n_lineages <- 1L
    upstream$pseudotime_per_lineage <- NULL
    upstream$lineage_cluster_source <-
      "uniform_single_lineage:no_lineage_specific_scoring"
  }
  lr_scores <- NULL
  temporal_lr <- NULL
  entropy_raw <- NULL
  h_tilde <- NULL
  bif_df <- NULL
  pair_gene_map <- NULL
  if (ablation_id != "no_lr_module" && sp %in% c("human", "mouse")) {
    lr_db <- BuildLRDB(expr_lr, sp)
    pair_gene_map <- BuildPairGeneMap(expr_lr, sp)
    lr_scores <- LineaLinkR::ComputeLRScores(expr_lr, celltypes, lr_db, n_perms = 100L)
    if (ablation_id != "no_temporal_binning") {
      pt_input <- if (!is.null(upstream$pseudotime_per_lineage)) upstream$pseudotime_per_lineage else upstream$pseudotime
      temporal_lr <- LineaLinkR::ComputeTemporalLRWithIntermediate(expr_lr, celltypes, lr_db, pt_input, upstream$lineage_weights, n_bins = 12L)
      if (ablation_id != "no_bifurcation_module") {
        entropy_raw <- LineaLinkR::ComputeLineageEntropyRaw(temporal_lr$intermediate)
        h_tilde <- LineaLinkR::SmoothMA3(entropy_raw)
        bif_df <- LineaLinkR::FindBifurcationBins(h_tilde)
      }
    }
  }
  result <- list(
    dataset = ds_name,
    timestamp = Sys.time(),
    ablation_id = ablation_id,
    execution_status = "pipeline_rerun_done",
    n_cells = n_cells,
    n_genes = nrow(expr),
    n_types = n_types,
    n_lineages = upstream$n_lineages,
    start_cell = upstream$start_cell,
    lineage_cluster_source = upstream$lineage_cluster_source,
    package_version = "0.6.4",
    benchmark_protocol = CONFIG$llr$benchmark_protocol,
    result_revision = CONFIG$llr$result_revision,
    root_strategy = upstream$root_strategy,
    root_support = upstream$root_support,
    direction_margin = upstream$direction_margin,
    direction_status = upstream$direction_status,
    pseudotime = upstream$pseudotime,
    pseudotime_per_lineage = upstream$pseudotime_per_lineage,
    lineage_weights = upstream$lineage_weights,
    lr_scores = lr_scores,
    temporal_lr = temporal_lr,
    entropy_raw = entropy_raw,
    h_tilde = h_tilde,
    bif_df = bif_df,
    pair_gene_map = pair_gene_map,
    config_snapshot = list(
      seed = CONFIG$seed,
      n_threads = CONFIG$n_threads,
      n_bins = 12L,
      benchmark_protocol = CONFIG$llr$benchmark_protocol,
      result_revision = CONFIG$llr$result_revision,
      algorithm = list(
        name = "LineaLinkR ablation module rerun",
        ablation_id = ablation_id,
        root_strategy = upstream$root_strategy,
        graph_k = as.integer(llr_cfg$knn_neighbors %||% 50L)
      )
    )
  )
  saveRDS(result, file.path(module_dir, "linealinkr_complete.rds"))
  fwrite(data.table(cell = cells, pseudotime = as.numeric(result$pseudotime)), file.path(module_dir, sprintf("%s_pseudotime.csv", ablation_id)))
  fwrite(SummarizeObject(result)[, `:=`(dataset = ds_name, ablation = ablation_id, execution_status = "pipeline_rerun_done")], file.path(module_dir, "module_summary.csv"))
  if (!is.null(result$temporal_lr) && !is.null(result$temporal_lr$scores)) fwrite(as.data.table(result$temporal_lr$scores), file.path(module_dir, sprintf("%s_temporal_lr.csv", ablation_id)))
  if (!is.null(result$bif_df)) fwrite(as.data.table(result$bif_df), file.path(module_dir, sprintf("%s_bifurcation.csv", ablation_id)))
  cat(sprintf("[MODULE] %s/%s: lr=%d temporal=%d bif=%d\n", ds_name, ablation_id, SummarizeObject(result)$n_lr_pairs, SummarizeObject(result)$n_temporal_rows, SummarizeObject(result)$n_bifurcation_points))
  result
}

ExtractGT <- function() {
  gt_csv <- file.path(ABL_DIR, sprintf("gt_%s.csv", ds_name))
  canon <- DatasetPath(ds_name)
  py <- sprintf(r"(
import anndata as ad, pandas as pd, numpy as np
a = ad.read_h5ad(r'%s')
gt_col = 'NOGT'
gt_values = np.full(a.n_obs, np.nan)
gt_valid = np.zeros(a.n_obs, dtype=bool)
for c in ['true_pseudotime','timepoint','day','stage']:
    if c in a.obs.columns:
        v = pd.to_numeric(a.obs[c], errors='coerce')
        vf = v.notna() & np.isfinite(v)
        if vf.sum() > 10 and v[vf].nunique() >= 2:
            gt_col = c
            gt_values = v.to_numpy()
            gt_valid = vf.to_numpy()
            break
branch_col = 'NOBRANCH'
branch_values = pd.Series([''] * a.n_obs, index=a.obs_names, dtype='object')
branch_valid = np.zeros(a.n_obs, dtype=bool)
for c in ['true_branch','true_branch_numeric','branch_truth','lineage_truth']:
    if c in a.obs.columns:
        bv = a.obs[c].astype('string')
        bf = bv.notna() & (bv != '')
        if bf.sum() > 0 and bv[bf].nunique() >= 2:
            branch_col = c
            branch_values = bv.fillna('')
            branch_valid = bf.to_numpy()
            break
pd.DataFrame({'cell':a.obs_names.astype(str),'gt':gt_values,'valid':gt_valid,'branch':branch_values.astype(str).to_numpy(),'branch_valid':branch_valid}).to_csv(r'%s',index=False)
print('GT=' + gt_col)
print('BRANCH=' + branch_col)
)", canon, gt_csv)
  tf <- file.path(ABL_DIR, ds_name, "_extract_gt.py")
  dir.create(dirname(tf), recursive = TRUE, showWarnings = FALSE)
  writeLines(py, tf, useBytes = TRUE)
  out <- system2(CONFIG$exe$python, args = tf, stdout = TRUE, stderr = TRUE)
  unlink(tf)
  gt_line <- sub("^GT=", "", grep("^GT=", out, value = TRUE)[1])
  branch_line <- sub("^BRANCH=", "", grep("^BRANCH=", out, value = TRUE)[1])
  gt_col <- if (is.na(gt_line) || gt_line == "NOGT") NA_character_ else gt_line
  branch_col <- if (is.na(branch_line) || branch_line == "NOBRANCH") NA_character_ else branch_line
  if (file.exists(gt_csv)) return(list(column = gt_col, branch_column = branch_col, table = fread(gt_csv)))
  list(column = NA_character_, branch_column = NA_character_, table = data.table(cell = character(), gt = numeric(), valid = logical(), branch = character(), branch_valid = logical()))
}

ComputePTMetrics <- function(pt_vec, gt_obj) {
  gt <- gt_obj$table
  if (!nrow(gt)) return(list(sp = NA_real_, pr = NA_real_, rmse = NA_real_, mae = NA_real_, nv = 0L))
  idx <- match(cells, gt$cell)
  gt_valid <- as.logical(gt$valid[idx])
  gt_valid[is.na(gt_valid)] <- FALSE
  gt_values <- as.numeric(gt$gt[idx])
  v <- gt_valid & is.finite(pt_vec) & !is.na(pt_vec)
  nv <- sum(v)
  if (nv < 10L) return(list(sp = NA_real_, pr = NA_real_, rmse = NA_real_, mae = NA_real_, nv = nv))
  list(
    sp = suppressWarnings(cor(gt_values[v], pt_vec[v], method = "spearman")),
    pr = suppressWarnings(cor(gt_values[v], pt_vec[v], method = "pearson")),
    rmse = sqrt(mean((gt_values[v] - pt_vec[v])^2)),
    mae = mean(abs(gt_values[v] - pt_vec[v])),
    nv = nv
  )
}

ComputeARI <- function(pred, truth) {
  tbl <- table(pred, truth)
  a <- sum(choose(tbl, 2))
  b <- sum(choose(rowSums(tbl), 2)) - a
  c2 <- sum(choose(colSums(tbl), 2)) - a
  expected <- (a + b) * (a + c2) / choose(sum(tbl), 2)
  denominator <- (a + b + a + c2) / 2 - expected
  if (denominator == 0) 0 else (a - expected) / denominator
}

ComputeBranchF1 <- function(pred, truth) {
  tbl <- as.matrix(table(pred, truth))
  tp <- tbl
  fp <- rowSums(tbl) - tbl
  fn <- matrix(colSums(tbl), nrow = nrow(tbl), ncol = ncol(tbl), byrow = TRUE) - tbl
  precision <- tp / (tp + fp)
  recall <- tp / (tp + fn)
  f1 <- 2 * precision * recall / (precision + recall)
  f1[!is.finite(f1)] <- 0
  mean(f1[cbind(max.col(t(f1), ties.method = "first"), seq_len(ncol(f1)))])
}

ComputeLineageMetrics <- function(lineage_weights, gt_obj) {
  gt <- gt_obj$table
  if (is.null(lineage_weights) || !nrow(gt)) return(list(ari = NA_real_, branch_f1 = NA_real_, nv = 0L))
  idx <- match(cells, gt$cell)
  truth <- as.character(gt$branch[idx])
  branch_valid <- as.logical(gt$branch_valid[idx])
  branch_valid[is.na(branch_valid)] <- FALSE
  lw <- as.matrix(lineage_weights)
  valid <- branch_valid & !is.na(truth) & truth != "" & is.finite(rowSums(lw))
  nv <- sum(valid)
  if (nv < 10L) return(list(ari = NA_real_, branch_f1 = NA_real_, nv = nv))
  pred <- max.col(lw, ties.method = "first")
  list(ari = ComputeARI(pred[valid], truth[valid]), branch_f1 = ComputeBranchF1(pred[valid], truth[valid]), nv = nv)
}

PTDiffers <- function(pt_vec) {
  if (is.null(pt_vec) || length(pt_vec) != length(full$pseudotime)) return(TRUE)
  !isTRUE(all.equal(as.numeric(pt_vec), as.numeric(full$pseudotime), tolerance = 1e-8, check.attributes = FALSE))
}

RunGraphKAblation <- function(compute_fn) {
  upstream <- compute_fn()
  list(
    pt = upstream$pseudotime,
    m = ComputePTMetrics(upstream$pseudotime, gt_obj),
    ml = ComputeLineageMetrics(upstream$lineage_weights, gt_obj),
    status = "pipeline_rerun_done",
    root_strategy = upstream$root_strategy,
    root_potency_mode = upstream$root_potency_mode,
    selected_root_potency_mode = upstream$selected_root_potency_mode,
    root_support = upstream$root_support,
    direction_margin = upstream$direction_margin,
    direction_status = upstream$direction_status
  )
}

ComputeGraphK10 <- function() {
  graph_dir <- file.path(ABL_DIR, ds_name, "graph_k_10")
  upstream <- ComputeUpstream(
    graph_dir,
    10L,
    ablation_id = "graph_k_10"
  )
  saveRDS(upstream, file.path(graph_dir, "graph_k_10_upstream.rds"))
  fwrite(data.table(cell = cells, pseudotime = as.numeric(upstream$pseudotime)), file.path(graph_dir, "graph_k_10_pseudotime.csv"))
  upstream
}

ComputeGraphK30 <- function() {
  graph_dir <- file.path(ABL_DIR, ds_name, "graph_k_30")
  upstream <- ComputeUpstream(
    graph_dir,
    30L,
    ablation_id = "graph_k_30"
  )
  saveRDS(upstream, file.path(graph_dir, "graph_k_30_upstream.rds"))
  fwrite(data.table(cell = cells, pseudotime = as.numeric(upstream$pseudotime)), file.path(graph_dir, "graph_k_30_pseudotime.csv"))
  upstream
}

DropDatasetRows <- function(dt, ds_name) {
  if (!is.data.table(dt) || nrow(dt) == 0L || !"dataset" %in% names(dt)) return(data.table())
  dt[dataset != ds_name]
}

RecomputeBins <- function(lr_db) {
  bin_grid <- c(5L, 10L, 12L, 20L)
  pt_input <- if (!is.null(full$pseudotime_per_lineage)) full$pseudotime_per_lineage else full$pseudotime
  rows <- rbindlist(lapply(bin_grid, function(nb) {
    tlr <- LineaLinkR::ComputeTemporalLRWithIntermediate(expr_lr, celltypes, lr_db, pt_input, full$lineage_weights, n_bins = nb)
    data.table(n_bins = nb, total_rows = nrow(tlr$scores), total_unique_pairs = uniqueN(tlr$scores$pair_idx), bins_used = uniqueN(tlr$scores$bin), dataset = ds_name, recompute_status = "real_ComputeTemporalLRWithIntermediate")
  }))
  out_path <- file.path(ABL_DIR, "ablation_bins_comparison.csv")
  old <- if (file.exists(out_path)) fread(out_path) else data.table()
  old <- DropDatasetRows(old, ds_name)
  fwrite(rbindlist(list(old, rows), fill = TRUE), out_path)
  rows
}

module_ids <- c(
  "no_lr_module",
  "no_temporal_binning",
  "no_bifurcation_module",
  "no_lineage_specific_scoring"
)
module_objects <- lapply(module_ids, RunModuleAblation)
names(module_objects) <- module_ids

no_kmeans_dir <- file.path(ABL_DIR, ds_name, "no_kmeans_refinement")
no_kmeans_upstream <- ComputeUpstream(
  no_kmeans_dir,
  as.integer(llr_cfg$knn_neighbors %||% 50L),
  reuse_baseline_graph = TRUE,
  root_clusters = lineage_cl_base,
  root_cluster_source = sprintf(
    "pre_refinement:%s:no_kmeans_refinement",
    fine_col
  ),
  ablation_id = "no_kmeans_refinement"
)
no_kmeans_upstream$config_snapshot <- list(
  seed = CONFIG$seed,
  graph_k = as.integer(llr_cfg$knn_neighbors %||% 50L),
  reused_baseline_knn_pca = TRUE,
  root_cluster_source = no_kmeans_upstream$root_cluster_source,
  lineage_cluster_source = no_kmeans_upstream$lineage_cluster_source
)
saveRDS(
  no_kmeans_upstream,
  file.path(no_kmeans_dir, "no_kmeans_refinement_upstream.rds")
)
fwrite(
  data.table(
    cell = cells,
    pseudotime = as.numeric(no_kmeans_upstream$pseudotime)
  ),
  file.path(no_kmeans_dir, "no_kmeans_refinement_pseudotime.csv")
)

gt_obj <- ExtractGT()
m_full <- ComputePTMetrics(full$pseudotime, gt_obj)
module_metrics <- lapply(module_objects, function(obj) ComputePTMetrics(obj$pseudotime, gt_obj))
m_full_lineage <- ComputeLineageMetrics(full$lineage_weights, gt_obj)
module_lineage_metrics <- lapply(module_objects, function(obj) ComputeLineageMetrics(obj$lineage_weights, gt_obj))
m_no_kmeans <- ComputePTMetrics(no_kmeans_upstream$pseudotime, gt_obj)
m_no_kmeans_lineage <- ComputeLineageMetrics(
  no_kmeans_upstream$lineage_weights,
  gt_obj
)

k10_res <- RunGraphKAblation(ComputeGraphK10)
k10_pt <- k10_res$pt
m_k10 <- k10_res$m
m_k10_lineage <- k10_res$ml
k30_res <- RunGraphKAblation(ComputeGraphK30)
k30_pt <- k30_res$pt
m_k30 <- k30_res$m
m_k30_lineage <- k30_res$ml

lr_db_full <- if (sp %in% c("human", "mouse")) BuildLRDB(expr_lr, sp) else list(ligand_rows = list(), receptor_rows = list())
bins_dt <- if (sp %in% c("human", "mouse")) RecomputeBins(lr_db_full) else data.table(n_bins = c(5L, 10L, 12L, 20L), total_rows = rep(0L, 4), total_unique_pairs = rep(0L, 4), bins_used = rep(0L, 4), dataset = ds_name, recompute_status = "skipped_non_human_mouse")
bin_temporal <- setNames(bins_dt$total_rows, paste0("bins_", bins_dt$n_bins))

summary_full <- SummarizeObject(full)
summary_modules <- rbindlist(lapply(names(module_objects), function(id) {
  SummarizeObject(module_objects[[id]])[, `:=`(
    dataset = ds_name,
    benchmark_protocol = CONFIG$llr$benchmark_protocol,
    result_revision = CONFIG$llr$result_revision,
    ablation = id,
    metric_layer = "robustness_ablation_module",
    execution_status = "pipeline_rerun_done"
  )]
}), fill = TRUE)

module_metrics_path <- file.path(ABL_DIR, "ablation_module_metrics.csv")
old_module_metrics <- if (file.exists(module_metrics_path)) fread(module_metrics_path) else data.table()
old_module_metrics <- DropDatasetRows(old_module_metrics, ds_name)
full_module_row <- summary_full[, `:=`(
  dataset = ds_name,
  benchmark_protocol = CONFIG$llr$benchmark_protocol,
  result_revision = CONFIG$llr$result_revision,
  ablation = "full",
  metric_layer = "robustness_ablation_module",
  execution_status = "baseline_observed"
)]
fwrite(
  rbindlist(list(old_module_metrics, full_module_row, summary_modules), fill = TRUE),
  module_metrics_path
)

metric_for <- function(id, slot) module_metrics[[id]][[slot]]
lineage_metric_for <- function(id, slot) module_lineage_metrics[[id]][[slot]]
summary_for <- function(id, col) SummarizeObject(module_objects[[id]])[[col]]
path_for <- function(id) file.path(ABL_DIR, ds_name, id, "linealinkr_complete.rds")

ablation_dt <- data.table(
  dataset = ds_name,
  benchmark_protocol = CONFIG$llr$benchmark_protocol,
  result_revision = CONFIG$llr$result_revision,
  ablation_id = c(
    "full", module_ids, "no_kmeans_refinement",
    "bins_5", "bins_10", "bins_12", "bins_20",
    "graph_k_10", "graph_k_30"
  ),
  Spearman = c(
    m_full$sp,
    metric_for("no_lr_module", "sp"),
    metric_for("no_temporal_binning", "sp"),
    metric_for("no_bifurcation_module", "sp"),
    metric_for("no_lineage_specific_scoring", "sp"),
    m_no_kmeans$sp,
    rep(m_full$sp, 4L),
    m_k10$sp, m_k30$sp
  ),
  Pearson = c(
    m_full$pr,
    metric_for("no_lr_module", "pr"),
    metric_for("no_temporal_binning", "pr"),
    metric_for("no_bifurcation_module", "pr"),
    metric_for("no_lineage_specific_scoring", "pr"),
    m_no_kmeans$pr,
    rep(m_full$pr, 4L),
    m_k10$pr, m_k30$pr
  ),
  RMSE = c(
    m_full$rmse,
    metric_for("no_lr_module", "rmse"),
    metric_for("no_temporal_binning", "rmse"),
    metric_for("no_bifurcation_module", "rmse"),
    metric_for("no_lineage_specific_scoring", "rmse"),
    m_no_kmeans$rmse,
    rep(m_full$rmse, 4L),
    m_k10$rmse, m_k30$rmse
  ),
  MAE = c(
    m_full$mae,
    metric_for("no_lr_module", "mae"),
    metric_for("no_temporal_binning", "mae"),
    metric_for("no_bifurcation_module", "mae"),
    metric_for("no_lineage_specific_scoring", "mae"),
    m_no_kmeans$mae,
    rep(m_full$mae, 4L),
    m_k10$mae, m_k30$mae
  ),
  n_valid = c(
    m_full$nv,
    metric_for("no_lr_module", "nv"),
    metric_for("no_temporal_binning", "nv"),
    metric_for("no_bifurcation_module", "nv"),
    metric_for("no_lineage_specific_scoring", "nv"),
    m_no_kmeans$nv,
    rep(m_full$nv, 4L),
    m_k10$nv, m_k30$nv
  ),
  ARI = c(
    m_full_lineage$ari,
    lineage_metric_for("no_lr_module", "ari"),
    lineage_metric_for("no_temporal_binning", "ari"),
    lineage_metric_for("no_bifurcation_module", "ari"),
    lineage_metric_for("no_lineage_specific_scoring", "ari"),
    m_no_kmeans_lineage$ari,
    rep(m_full_lineage$ari, 4L),
    m_k10_lineage$ari, m_k30_lineage$ari
  ),
  Branch_F1 = c(
    m_full_lineage$branch_f1,
    lineage_metric_for("no_lr_module", "branch_f1"),
    lineage_metric_for("no_temporal_binning", "branch_f1"),
    lineage_metric_for("no_bifurcation_module", "branch_f1"),
    lineage_metric_for("no_lineage_specific_scoring", "branch_f1"),
    m_no_kmeans_lineage$branch_f1,
    rep(m_full_lineage$branch_f1, 4L),
    m_k10_lineage$branch_f1, m_k30_lineage$branch_f1
  ),
  n_branch_valid = c(
    m_full_lineage$nv,
    lineage_metric_for("no_lr_module", "nv"),
    lineage_metric_for("no_temporal_binning", "nv"),
    lineage_metric_for("no_bifurcation_module", "nv"),
    lineage_metric_for("no_lineage_specific_scoring", "nv"),
    m_no_kmeans_lineage$nv,
    rep(m_full_lineage$nv, 4L),
    m_k10_lineage$nv, m_k30_lineage$nv
  ),
  n_total = n_cells,
  pseudotime_differs_from_full = c(
    FALSE,
    PTDiffers(module_objects$no_lr_module$pseudotime),
    PTDiffers(module_objects$no_temporal_binning$pseudotime),
    PTDiffers(module_objects$no_bifurcation_module$pseudotime),
    PTDiffers(module_objects$no_lineage_specific_scoring$pseudotime),
    PTDiffers(no_kmeans_upstream$pseudotime),
    rep(FALSE, 4L),
    PTDiffers(k10_pt), PTDiffers(k30_pt)
  ),
  graph_k_rerun_done = c(rep(FALSE, 10L), TRUE, TRUE),
  root_strategy = c(
    rep(llr_cfg$root_strategy %||% "attractor_energy", 5L),
    no_kmeans_upstream$root_strategy,
    rep(llr_cfg$root_strategy %||% "attractor_energy", 4L),
    k10_res$root_strategy,
    k30_res$root_strategy
  ),
  root_potency_mode = c(
    rep(llr_cfg$root_potency_mode %||% "detected_genes_ordering", 5L),
    no_kmeans_upstream$root_potency_mode,
    rep(llr_cfg$root_potency_mode %||% "detected_genes_ordering", 4L),
    k10_res$root_potency_mode,
    k30_res$root_potency_mode
  ),
  selected_root_potency_mode = c(
    rep(full$selected_root_potency_mode %||% NA_character_, 5L),
    no_kmeans_upstream$selected_root_potency_mode,
    rep(full$selected_root_potency_mode %||% NA_character_, 4L),
    k10_res$selected_root_potency_mode,
    k30_res$selected_root_potency_mode
  ),
  root_support = c(
    rep(full$root_support %||% NA_real_, 5L),
    no_kmeans_upstream$root_support,
    rep(full$root_support %||% NA_real_, 4L),
    k10_res$root_support, k30_res$root_support
  ),
  direction_margin = c(
    rep(full$direction_margin %||% NA_real_, 5L),
    no_kmeans_upstream$direction_margin,
    rep(full$direction_margin %||% NA_real_, 4L),
    k10_res$direction_margin, k30_res$direction_margin
  ),
  direction_status = c(
    rep(full$direction_status %||% NA_character_, 5L),
    no_kmeans_upstream$direction_status,
    rep(full$direction_status %||% NA_character_, 4L),
    k10_res$direction_status, k30_res$direction_status
  ),
  gt_column = gt_obj$column,
  branch_truth_column = gt_obj$branch_column,
  metric_layer = "robustness_ablation",
  module_n_lr_pairs = c(
    summary_full$n_lr_pairs,
    summary_for("no_lr_module", "n_lr_pairs"),
    summary_for("no_temporal_binning", "n_lr_pairs"),
    summary_for("no_bifurcation_module", "n_lr_pairs"),
    summary_for("no_lineage_specific_scoring", "n_lr_pairs"),
    rep(NA_real_, 7L)
  ),
  module_n_temporal_rows = c(
    summary_full$n_temporal_rows,
    summary_for("no_lr_module", "n_temporal_rows"),
    summary_for("no_temporal_binning", "n_temporal_rows"),
    summary_for("no_bifurcation_module", "n_temporal_rows"),
    summary_for("no_lineage_specific_scoring", "n_temporal_rows"),
    NA_real_,
    unname(bin_temporal["bins_5"]),
    unname(bin_temporal["bins_10"]),
    unname(bin_temporal["bins_12"]),
    unname(bin_temporal["bins_20"]),
    NA_real_, NA_real_
  ),
  module_n_bif_points = c(
    summary_full$n_bifurcation_points,
    summary_for("no_lr_module", "n_bifurcation_points"),
    summary_for("no_temporal_binning", "n_bifurcation_points"),
    summary_for("no_bifurcation_module", "n_bifurcation_points"),
    summary_for("no_lineage_specific_scoring", "n_bifurcation_points"),
    rep(NA_real_, 7L)
  ),
  execution_status = c(
    "baseline_observed",
    rep("pipeline_rerun_done", 5L),
    rep("real_temporal_lr_recompute", 4L),
    k10_res$status, k30_res$status
  ),
  module_output_path = c(
    baseline_rds,
    path_for("no_lr_module"),
    path_for("no_temporal_binning"),
    path_for("no_bifurcation_module"),
    path_for("no_lineage_specific_scoring"),
    file.path(no_kmeans_dir, "no_kmeans_refinement_upstream.rds"),
    rep(NA_character_, 4L),
    file.path(ABL_DIR, ds_name, "graph_k_10", "graph_k_10_pseudotime.csv"),
    file.path(ABL_DIR, ds_name, "graph_k_30", "graph_k_30_pseudotime.csv")
  )
)

ordering_metrics_path <- file.path(ABL_DIR, "ablation_ordering_metrics.csv")
old_ordering_metrics <- if (file.exists(ordering_metrics_path)) fread(ordering_metrics_path) else data.table()
old_ordering_metrics <- DropDatasetRows(old_ordering_metrics, ds_name)
fwrite(rbindlist(list(old_ordering_metrics, ablation_dt), fill = TRUE), ordering_metrics_path)

out_path <- file.path(ABL_DIR, sprintf("ablation_%s.csv", ds_name))
fwrite(ablation_dt, out_path)
cat(sprintf("[ABLATION] %s: %d rows, full Spearman=%.4f, k10 Spearman=%.4f\n", ds_name, nrow(ablation_dt), m_full$sp, m_k10$sp))
cat("DONE\n")
