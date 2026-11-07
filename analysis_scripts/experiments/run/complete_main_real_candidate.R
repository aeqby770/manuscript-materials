# Sourced by run_linealinkr.R inside one dataset call. Only unchanged root,
# geometry and reference paths/weights are reused. Static LR is dependency checked.
COMPLETE <- list(
  input = file.path(CONFIG$dir$cache, anchor_workspace, ds_name, "linealinkr_complete.rds"),
  output = file.path(result_dir, "linealinkr_complete.rds"),
  source_root = Sys.getenv("LLR_SOURCE_ROOT", unset = "F:/biotrade/LineaLinkR_0.6.5"),
  threads = CONFIG$n_threads
)
Sys.setenv(OMP_NUM_THREADS = COMPLETE$threads)
data.table::setDTthreads(COMPLETE$threads)
started <- proc.time()[["elapsed"]]
x <- readRDS(COMPLETE$input)
x$input_extract_dir <- sub(COMPLETE$source_root, .llr_work_root, x$input_extract_dir, fixed = TRUE)
x$knn_precompute_dir <- sub(COMPLETE$source_root, .llr_work_root, x$knn_precompute_dir, fixed = TRUE)
cluster_input_policy <- c("metadata_or_expression_auto_clusters", "expression_auto_clusters")[
  1L + as.integer(ds_name %chin% BENCHMARK_FORCE_AUTO_CLUSTER_DATASETS)]
stopifnot(identical(x$dataset, ds_name), main_candidate,
  identical(x$cluster_input_policy, cluster_input_policy),
  !program_lineage, !resident_sender, !persistent_rna_states,
  x$config_snapshot$seed == CONFIG$seed)
stopifnot(!is.null(x$lineage_reference_weights) ||
  is.null(x$config_snapshot$lineage_weight_mode) ||
  identical(x$config_snapshot$lineage_weight_mode,"reference_geometry"))
llr_cfg <- CONFIG$llr
extract <- x$input_extract_dir
genes <- readLines(file.path(extract, paste0(ds_name, "_genes.txt")), encoding = "UTF-8")
cells <- readLines(file.path(extract, paste0(ds_name, "_cells.txt")), encoding = "UTF-8")
meta <- fread(file.path(extract, paste0(ds_name, "_meta.csv")))
count_entries <- fread(file.path(extract, paste0(ds_name, "_matrix.mtx")),
  header = FALSE, col.names = c("i", "j", "value"))
count_dimensions <- as.integer(unlist(count_entries[1L, .(i, j)], use.names = FALSE))
stopifnot(nrow(count_entries) - 1L == count_entries$value[[1L]])
expr <- sparseMatrix(i = count_entries$i[-1L], j = count_entries$j[-1L],
  x = as.numeric(count_entries$value[-1L]), dims = count_dimensions)
rm(count_entries)
invisible(gc(verbose = FALSE))
rownames(expr) <- genes
colnames(expr) <- cells
stopifnot(identical(cells, x$cells), ncol(expr) == x$n_cells, nrow(expr) == x$n_genes,
  identical(x$package_version, as.character(packageVersion("LineaLinkR"))))
geometry <- switch(as.character(isTRUE(x$geometry_centers_lr) || isTRUE(x$geometry_excludes_lr)),
  "TRUE" = file.path(x$knn_precompute_dir, "structural_geometry"),
  "FALSE" = x$knn_precompute_dir)
geom_emb <- as.matrix(fread(file.path(geometry, "pca.csv"), header = FALSE))
state_clock_conn <- as(as(readMM(file.path(geometry, "knn_conn.mtx")), "generalMatrix"), "CsparseMatrix")
stopifnot(nrow(geom_emb) == length(cells), nrow(state_clock_conn) == length(cells))
n_cells <- length(cells)
start_cell <- x$start_cell
selected_root_potency_mode <- x$selected_root_potency_mode
root_evidence_mode <- x$root_evidence_mode
root_input_source <- x$root_input_source
lineage_cl <- x$lineage_clusters
lineage_paths <- x$lineage_paths
lineage_reference_weights <- x$lineage_reference_weights %||% x$lineage_weights
lineage_assignment_weights <- x$lineage_assignment_weights %||% x$lineage_priming_weights
lineage_weights <- lineage_reference_weights
lineage_priming_weights <- LineaLinkR::CombineLineageMemberships(
  lineage_reference_weights, lineage_assignment_weights)
switch(component_ablation, local_branch_reference =, local_point_reference = {
  stopifnot(identical(dim(lineage_priming_weights),dim(lineage_weights)),
    all(is.finite(lineage_priming_weights)),
    max(abs(rowSums(lineage_priming_weights)-1))<1e-5)
  lineage_weights <- lineage_assignment_weights
  x$lineage_weights <- lineage_weights
  x$lineage_membership <- max.col(lineage_weights,ties.method="first")
}, supported_branch_reference = {
  lineage_weights <- LineaLinkR::CombineLineageMemberships(
    lineage_reference_weights, lineage_assignment_weights)
  x$lineage_weights <- lineage_weights
  x$lineage_membership <- max.col(lineage_weights,ties.method="first")
})
switch(component_ablation, observed_branch_reference =, uniform_branch_reference = {
  native_path_table <- data.table(node=unlist(lineage_paths,use.names=FALSE),
    lineage=rep(seq_along(lineage_paths),lengths(lineage_paths)))
  native_path_nodes <- sort(unique(native_path_table$node))
  native_path_incidence <- sparseMatrix(i=match(native_path_table$node,native_path_nodes),
    j=native_path_table$lineage,x=1,dims=c(length(native_path_nodes),length(lineage_paths)))
  native_path_mask <- as.matrix(native_path_incidence[
    match(lineage_cl,native_path_nodes),,drop=FALSE])
  native_path_reference <- switch(component_ablation,
    uniform_branch_reference=matrix(1,nrow(lineage_reference_weights),ncol(lineage_reference_weights)),
    lineage_reference_weights)
  native_path_weights <- native_path_reference*native_path_mask
  stopifnot(all(is.finite(native_path_weights)),all(rowSums(native_path_weights)>0))
  lineage_weights <- native_path_weights/rowSums(native_path_weights)
  x$lineage_weights <- lineage_weights
  x$lineage_membership <- max.col(lineage_weights,ties.method='first')
})
x$config_snapshot$component_ablation <- component_ablation
x$config_snapshot$activity_transform <- llr_cfg$event_activity_transform
x$config_snapshot$context_error_weighting <- llr_cfg$context_error_weighting
x$config_snapshot$late_bin_merging_enabled <- enable_late_bin_merging
x$config_snapshot$lineage_weight_mode <- CONFIG$llr$lineage_weight_mode
x$config_snapshot$lineage_curve_weight_mode <- CONFIG$llr$lineage_curve_weight_mode
x$config_snapshot$lineage_persistence_weight_mode <- ifelse(use_assignment_layout,
  "reference_geometry",CONFIG$llr$lineage_curve_weight_mode)
n_lineages <- x$n_lineages
sibling_map <- x$sibling_map
switch(as.character(enable_knn_topology), "TRUE" = {
  switch(as.character(nzchar(lineage_state_column)), "TRUE" = {
    lineage_cl <- as.integer(factor(as.character(meta[[lineage_state_column]])))
  })
  start_cluster <- as.integer(lineage_cl[start_cell])
  context_tree <- LineaLinkR::BuildKnnClusterTree(state_clock_conn, lineage_cl, start_cluster)
  lineage_paths <- LineaLinkR::ExtractLineagesFromKnnGraph(
    state_clock_conn, lineage_cl, start_cluster)$paths
  lineage_weights <- LineaLinkR::ComputeLineageWeightsFromKnnGraph(
    state_clock_conn, geom_emb, lineage_cl, start_cluster, tree=context_tree)
  stopifnot(all(rowSums(lineage_weights)>0))
  lineage_weights <- lineage_weights/rowSums(lineage_weights)
  lineage_reference_weights <- lineage_assignment_weights <- lineage_weights
  lineage_priming_weights <- lineage_weights
  n_lineages <- ncol(lineage_weights)
  sibling_map <- switch(as.character(n_lineages>1L),
    "TRUE"=llr_event_methods$BuildSiblingContrastMap(lineage_paths),
    "FALSE"=data.frame(focal_lineage=integer(),sibling_lineage=integer(),
      fork_cluster=integer(),shared_prefix_length=integer()))
  x$knn_graph_topology <- TRUE
  x$config_snapshot$knn_graph_topology <- TRUE
  x$lineage_clusters <- as.integer(lineage_cl)
  x$lineage_paths <- lineage_paths
  x$lineage_weights <- lineage_weights
  x$lineage_priming_weights <- lineage_priming_weights
  x$lineage_membership <- max.col(lineage_weights,ties.method="first")
  x$n_lineages <- n_lineages
  x$sibling_map <- sibling_map
  x$lineage_cluster_source <- paste0("observed_knn_graph:", 
    ifelse(nzchar(lineage_state_column),lineage_state_column,x$lineage_cluster_source))
})
x$lineage_reference_weights <- lineage_reference_weights
x$lineage_assignment_weights <- lineage_assignment_weights
x$lineage_priming_weights <- lineage_priming_weights
x$lineage_weights <- lineage_weights
x$config_snapshot$lineage_priming_mode <- CONFIG$llr$lineage_priming_mode
x$config_snapshot$lineage_assignment_mode <- CONFIG$llr$lineage_assignment_mode
pseudotime <- switch(as.character(is.null(x$state_clock)),
  "TRUE" = x$pseudotime,
  "FALSE" = x$state_clock$geometry_anchor_order %||% x$state_clock$anchor_order)
pseudotime_per_lineage <- NULL
lineage_cluster_source <- x$lineage_cluster_source
root_diagnostics <- x$root_diagnostics
direction_diagnostics <- x$direction_diagnostics
celltypes <- x$celltypes
transition_fragment_labels <- integer()
program_profiles <- NULL
if (rescore_events) {
  pseudotime <- x$pseudotime
  pseudotime_per_lineage <- x$pseudotime_per_lineage
  state_clock <- x$state_clock
  direction_status <- x$direction_status
  branch_clock_calibration <- x$branch_clock_calibration
} else {
  source(file.path(.llr_work_root,"experiments/run/compute_state_clock.R"),
    local=TRUE,encoding="UTF-8")
}

# Count normalization, LR subunit definitions and missing-resource handling use
# the same current entry point; static scores are refreshed when the database changes.
sp <- CONFIG$dataset_species[[ds_name]] %||% "human"
input_scale_requested <- unname(ifelse(
  ds_name %in% names(CONFIG$dataset_lr_input_scale),
  CONFIG$dataset_lr_input_scale[ds_name], x$lr_input_scale))
expr_lr <- lr_db <- NULL
switch(as.character(sp %in% c("human", "mouse", "zebrafish", "celegans")), "TRUE" = {
  expr_lr <- PrepareLRExpression(expr, llr_cfg$lr_scale_factor, input_scale_requested)
  lr_db <- BuildLRDB(expr_lr, species = sp)
  switch(as.character(length(lr_db$ligand_rows) == 0L), "TRUE" = lr_db <- NULL)
})
# A changed database invalidates both static scores and their pair indices.
reuse_static_lr <- identical(as.data.frame(lr_db$pair_gene_map),
  as.data.frame(x$pair_gene_map)) && identical(input_scale_requested, x$lr_input_scale)
lr_scores <- x$lr_scores
switch(as.character(reuse_static_lr), "FALSE" = {
  cat("  Step 5/7: LR scores (static; current database) ... ")
  lr_scores <- LineaLinkR::ComputeLRScores(expr_lr, celltypes, lr_db, n_perms = 100L)
  x$lr_scores <- lr_scores
  x$pair_gene_map <- lr_db$pair_gene_map
  cat(sprintf("done [%d molecular pairs; %d context scores]\n",
    nrow(x$pair_gene_map), nrow(lr_scores)))
})
x$static_lr_reused <- reuse_static_lr
x$lr_input_scale <- input_scale_requested
x$config_snapshot$lr_database_version <- ifelse(
  sp %in% c("human", "mouse", "zebrafish"),
  as.character(packageVersion("CellChat")), NA_character_)
event_expression_mean_mode <- x$event_expression_mean_mode
lr_detection_intensity_rho <- x$lr_detection_intensity_rho
# Release obsolete large dynamic tables before recomputing every affected stage.
x$temporal_lr <- x$fork_contrasts <- x$fork_evidence <- NULL
x$fork_sample_contrasts <- x$fork_sample_summary <- NULL
invisible(gc(verbose = FALSE))
source(file.path(.llr_work_root, "experiments/run/compute_communication.R"),
  local = TRUE, encoding = "UTF-8")

old_revision <- x$result_revision
x$pseudotime <- pseudotime
x$pseudotime_origin_cell <- which.min(pseudotime)
x$pseudotime_per_lineage <- pseudotime_per_lineage
x$branch_clock_calibration <- branch_clock_calibration
x$state_clock <- state_clock
x$lineage_cluster_source <- lineage_cluster_source
x$direction_status <- direction_status
x$direction_diagnostics <- direction_diagnostics
x$root_diagnostics <- root_diagnostics
x$lineage_curve_weights <- lineage_curve_weights
x$lineage_shape_weights <- lineage_shape_weights
x$shape_temporal_lr <- shape_temporal_lr
x$reference_temporal_lr <- reference_temporal_lr
x$temporal_lr <- temporal_lr
x$sample_temporal_lr <- sample_temporal_lr
x$calendar_temporal_lr <- calendar_temporal_lr
x$calendar_fork_dynamics <- calendar_fork_dynamics
x$calendar_inference_status <- calendar_inference_status
x$entropy_raw <- entropy_raw
x$h_tilde <- h_tilde
x$bif_df <- bif_df
x$fork_bins <- fork_bins
x$fork_contrasts <- fork_contrasts
x$fork_evidence <- fork_evidence
x$ligand_nominations <- attr(fork_evidence,"ligand_nominations")
x$fork_sample_contrasts <- fork_sample_contrasts
x$fork_sample_summary <- fork_sample_summary
x$fork_sample_column <- fork_sample_column
x$fork_inference_unit <- fork_inference_unit
x$event_expression_mean_mode <- event_expression_mean_mode
x$lr_detection_intensity_rho <- lr_detection_intensity_rho
x$result_revision <- x$config_snapshot$result_revision <- CONFIG$llr$result_revision
x$config_snapshot$no_fork_control <- CONFIG$llr$no_fork_control
x$config_snapshot$prefork_weighting <- CONFIG$llr$prefork_weighting
x$config_snapshot$sample_amplitude_rule <- CONFIG$llr$sample_amplitude_rule
x$config_snapshot$persistence_operator <- CONFIG$llr$persistence_operator
x$config_snapshot$persistence_sibling_scope <- CONFIG$llr$persistence_sibling_scope
x$config_snapshot$algorithm <- sub(old_revision, x$result_revision,
  x$config_snapshot$algorithm, fixed = TRUE)
x$upstream_result <- COMPLETE$input
x$upstream_result_timestamp <- x$result_timestamp
x$updated_stages <- paste0(
  ifelse(rescore_events,"fixed upstream geometry and clock; event amplitudes, persistence, dynamics and sample support",
    "anchored expression clock; temporal LR; fork contrasts, dynamics, sample support and ranking"),
  ifelse(reuse_static_lr, "", "; static LR and current database pair mapping"))
x$reused_stages <- paste0(
  ifelse(component_ablation %in% c("local_branch_reference","supported_branch_reference"),
    "unchanged full-cell preprocessing, selected roots, global geometry and reference paths; current configured lineage memberships replace geometric weights",
    "unchanged full-cell preprocessing, selected roots, global geometry, reference paths/weights"),
  ifelse(reuse_static_lr, ", static LR", ""))
x$timestamp <- Sys.time()
x$result_timestamp <- format(x$timestamp, "%Y-%m-%dT%H:%M:%S%z")
x$updated_stages_elapsed_seconds <- proc.time()[["elapsed"]] - started
stopifnot(LlrRdsComplete(x), length(x$pseudotime) == x$n_cells,
  all(is.finite(x$pseudotime)), all(is.finite(x$lineage_weights)),
  max(abs(rowSums(x$lineage_weights) - 1)) < 1e-5)
cache_tmp <- file.path(result_dir, sprintf(".linealinkr_complete.rds.tmp-%s", Sys.getpid()))
saveRDS(x, cache_tmp)
saved_check <- readRDS(cache_tmp)
stopifnot(LlrRdsComplete(saved_check), identical(saved_check$pseudotime, pseudotime))
stopifnot(file.rename(cache_tmp, COMPLETE$output))
cat(sprintf("\n  >>> Saved: %s [affected stages %.1fs]\n", COMPLETE$output, x$updated_stages_elapsed_seconds))
