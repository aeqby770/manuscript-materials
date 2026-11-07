#!/usr/bin/env Rscript
# =============================================================================
# run_linealinkr.R
# LineaLinkR benchmark runner, separate from the competitor entry point.
# =============================================================================
# Usage:
#   Rscript run_linealinkr.R --main-065 --late-bin-persistence
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
# Resolve the isolated 0.6.5 revision checkout for direct entry points.
.llr_work_root <- normalizePath(file.path(.llr_here, "../.."), winslash = "/", mustWork = TRUE)
main_requested <- !any(grepl("^--(variant|run-tag)=", commandArgs(trailingOnly = TRUE)))
late_bin_persistence <- main_requested || "--late-bin-persistence" %in% commandArgs(trailingOnly = TRUE)
continuous_fork_priming <- "--continuous-fork-priming" %in% commandArgs(trailingOnly = TRUE)
persistent_rna_states <- "--persistent-rna-states" %in% commandArgs(trailingOnly = TRUE)
program_lineage <- "--program-lineage" %in% commandArgs(trailingOnly = TRUE) || persistent_rna_states
state_time_column <- Sys.getenv("LLR_STATE_TIME_COLUMN", "")
state_condition_column <- Sys.getenv("LLR_STATE_CONDITION_COLUMN", "")
resident_sender <- "--resident-sender" %in% commandArgs(trailingOnly = TRUE)
resident_stage_column <- Sys.getenv("LLR_RESIDENT_CONTEXT_COLUMN", "")
resident_time_column <- Sys.getenv("LLR_RESIDENT_TIME_COLUMN", "")
resident_condition_column <- Sys.getenv("LLR_RESIDENT_CONDITION_COLUMN", "")
anchor_workspace <- Sys.getenv("LLR_ANCHOR_WORKSPACE", unset = "none") # Historical direction-anchor input; outputs use feasible_clock.
anchor_workspace <- switch(anchor_workspace, none = "", anchor_workspace)
precompute_workspace <- Sys.getenv("LLR_PRECOMPUTE_WORKSPACE", unset = "")
stopifnot(!nzchar(precompute_workspace))
Sys.setenv(
  LLR_WORK_ROOT = .llr_work_root,
  LLR_PIPELINE = file.path(.llr_work_root, "experiments"),
  LLR_ASSETS = .llr_work_root
)
Sys.setenv(LLR_RLIB = Sys.getenv("LLR_PACKAGE_LIBRARY",unset=file.path(.llr_work_root,"library")))
source(file.path(.llr_here, "../config/llr_bootstrap.R"), local = FALSE)

  library(data.table)
})

llr_source_config()
if (late_bin_persistence)
  CONFIG$llr$result_revision <- paste0(CONFIG$llr$result_revision,
                                       "_late_bin_persistence")
if (continuous_fork_priming)
  CONFIG$llr$result_revision <- paste0(CONFIG$llr$result_revision,
                                       "_continuous_fork_priming")
switch(as.character(program_lineage), "TRUE" = {
  CONFIG$llr$result_revision <- paste0(CONFIG$llr$result_revision, "_program_lineage")
})
switch(as.character(persistent_rna_states), "TRUE" = {
  CONFIG$llr$result_revision <- paste0(CONFIG$llr$result_revision, "_persistent_rna_states")
})
switch(as.character(resident_sender), "TRUE" = {
  CONFIG$llr$result_revision <- paste0(CONFIG$llr$result_revision, "_resident_sender")
})
external_dataset <- Sys.getenv("LLR_EXTERNAL_DATASET", unset = "")
switch(as.character(nzchar(external_dataset)), "TRUE" = {
  CONFIG$dataset_narrative_registry[[external_dataset]] <- normalizePath(
    Sys.getenv("LLR_EXTERNAL_H5AD"), winslash = "/", mustWork = TRUE)
  CONFIG$dataset_species[[external_dataset]] <- Sys.getenv("LLR_EXTERNAL_SPECIES")
})
source(file.path(.llr_here, "../config/lr_database.R"), encoding = "UTF-8")

llr_package_library <- Sys.getenv(
  "LLR_RLIB",
  unset = file.path(CONFIG$dir$root, "library")
)
.libPaths(c(llr_package_library, .libPaths()))

`%||%` <- function(x, y) if (is.null(x)) y else x

IsValidClusterCol <- function(meta, cn, min_groups = 2L) {
  if (!cn %in% names(meta)) return(FALSE)
  vals <- as.character(meta[[cn]])
  vals[is.na(vals) | trimws(vals) == ""] <- NA_character_
  vals <- vals[!is.na(vals)]
  nu <- uniqueN(vals)
  # Two states (e.g. a single sender and a single receiver population) are a
  # legitimate LR input; do not force callers to invent a third category.
  if (nu < min_groups) return(FALSE)
  n <- length(vals)
  id_like <- grepl("patient|donor|sample_id|sample_name|status|diagnosis|orig\\.ident", cn, ignore.case = TRUE)
  if (id_like && !grepl("cluster", cn, ignore.case = TRUE)) return(FALSE)
  if (nu > max(80L, ceiling(0.2 * n)) && !grepl("cluster", cn, ignore.case = TRUE)) return(FALSE)
  TRUE
}

# Load LineaLinkR.
library(LineaLinkR, lib.loc = llr_package_library)
llr_event_methods <- asNamespace("LineaLinkR")
profile_file <- Sys.getenv("LLR_PROFILE_OUTPUT", "")
switch(as.character(nzchar(profile_file)), "TRUE" = {
  utils::Rprof(profile_file, memory.profiling = TRUE)
})

# Explicit variants define cache paths and scorer parameters together.
args <- commandArgs(trailingOnly = TRUE)
arg_at <- regexpr("=", args, fixed = TRUE)
arg_keys <- substr(args, 1L, ifelse(arg_at < 0L, nchar(args), arg_at - 1L))
arg_values <- substr(args, ifelse(arg_at < 0L, nchar(args) + 1L, arg_at + 1L), nchar(args))
opts <- setNames(arg_values, arg_keys)
allowed <- c("--variant", "--run-tag", "--datasets", "--datasets-file", "--panel",
             "--main-065", "--late-bin-persistence", "--continuous-fork-priming", "--program-lineage", "--resident-sender",
             "--persistent-rna-states", "--shared-parent-signed",
             "--entropy-only", "--no-temporal", "--no-entropy", "--max-cells",
             "--technical-smoke", "--validate-only", "--rescore-events", "--help", "--root-cells-file",
             "--lineage-topology-k", "--root-topology-k", "--knn-topology",
             "--ordering-low-anchor", "--terminal-sibling-refine",
             "--internal-transition-cluster", "--internal-transition-pocket",
             "--timepoint-branchpoint-root", "--timepoint-dominant-early-leaf",
             "--timepoint-source-flow-guard", "--local-branch-residual",
             "--root-geometric-gate", "--separate-priming-weights",
             "--donor-split-geometry", "--donor-center-lr-geometry",
             "--dropout-aware-event-mean",
             "--sparse-overlap-pc1", "--preserve-selected-root")
stopifnot(!anyDuplicated(arg_keys), all(arg_keys %in% allowed))
stopifnot(!persistent_rna_states || nzchar(state_time_column))
main_candidate <- main_requested || "--main-065" %in% arg_keys
rescore_events <- "--rescore-events" %in% arg_keys
stopifnot(!rescore_events || (main_candidate && nzchar(anchor_workspace)))
shared_parent_signed <- main_candidate
stopifnot(!shared_parent_signed || main_candidate)
event_contrast_method <- llr_event_methods$ComputeForkContrastsEvent
stopifnot(!main_candidate ||
  !any(c("--variant", "--run-tag") %in% arg_keys))
if ("--help" %in% arg_keys) {
  cat("--main-065 --late-bin-persistence --datasets=pancreas,dentategyrus | --panel=ordering\n",
      "--validate-only checks selection and the installed scorer without running data.\n")
  quit(save = "no", status = 0L)
}
variants <- fread(file.path(.llr_work_root, "experiments/config/variants.csv"))
variant_id <- if (main_candidate) "L1D0contTargetedInterval" else "L1D0cont"
if ("--run-tag" %in% arg_keys) {
  aliases <- c(baseline = "0.6.5", A = "A", B = "B", AB = "AB",
               L0D0 = "L0D0fix", L0D1 = "L0D1fix", L1D0 = "L1D0fix", L1D1 = "L1D1fix")
  variant_id <- unname(aliases[opts[["--run-tag"]]])
}
if ("--variant" %in% arg_keys) variant_id <- opts[["--variant"]]
variant <- variants[variant == variant_id]
stopifnot(nrow(variant) == 1L, variant$runnable %in% TRUE)
run_tag <- variant_id
variant_id <- switch(run_tag,
  L1D0contBranchBalanceInterval = "L1D0contTargetedInterval",
  L1D0contMaturityMatchedInterval = "L1D0contTargetedInterval",
  L1D0contDonorSplitGeometryInterval = "L1D0contTargetedInterval",
  L1D0contDonorBranchBalanceInterval = "L1D0contDonorStructuralInterval",
  L1D0contDonorBalanceNoSpurInterval = "L1D0contDonorStructuralInterval",
  run_tag)
donor_split_geometry_mode <- run_tag %in% c(
  "L1D0contDonorSplitGeometryInterval", "L1D0contMaturityMatchedInterval") ||
  "--donor-split-geometry" %in% arg_keys
donor_center_lr_mode <- main_candidate ||
  "--donor-center-lr-geometry" %in% arg_keys
stopifnot(!donor_center_lr_mode || !donor_split_geometry_mode)
filter_datasets <- character()
if ("--datasets" %in% arg_keys) {
  filter_datasets <- trimws(strsplit(opts[["--datasets"]], ",", fixed = TRUE)[[1L]])
}
if ("--datasets-file" %in% arg_keys) {
  filter_datasets <- trimws(readLines(opts[["--datasets-file"]], encoding = "UTF-8"))
}
if ("--panel" %in% arg_keys) {
  stopifnot(opts[["--panel"]] == "ordering")
  filter_datasets <- BENCHMARK_ORDERING_DATASETS
}
filter_datasets <- unique(filter_datasets[nzchar(filter_datasets)])
stopifnot(all(filter_datasets %in% c(names(CONFIG$dataset_registry),
                                    names(CONFIG$dataset_narrative_registry))))
if (!length(filter_datasets)) stop("Specify --datasets=... or --panel=ordering")
if ("--root-cells-file" %in% arg_keys) {
  stopifnot(length(filter_datasets) == 1L,
            main_candidate || run_tag %in% c(
              "L1D0contRoot", "L1D0contBranchBalanceInterval",
              "L1D0contMaturityMatchedInterval"))
  CONFIG$llr$user_root_cells_file <- normalizePath(opts[["--root-cells-file"]],
    winslash = "/", mustWork = TRUE)
  CONFIG$llr$user_root_cells <- sort(readLines(CONFIG$llr$user_root_cells_file,
    encoding = "UTF-8"))
  CONFIG$llr$user_root_evidence_source <- Sys.getenv("LLR_ROOT_PRIOR_SOURCE", "")
  stopifnot(length(CONFIG$llr$user_root_cells) > 0L,
    all(nzchar(CONFIG$llr$user_root_cells)), !anyDuplicated(CONFIG$llr$user_root_cells))
}
if (variant_id == "L1D0contRoot") stopifnot("--root-cells-file" %in% arg_keys)
entropy_only <- "--entropy-only" %in% arg_keys
no_temporal <- "--no-temporal" %in% arg_keys
no_entropy <- "--no-entropy" %in% arg_keys
technical_smoke <- "--technical-smoke" %in% arg_keys
max_cells <- NULL
if ("--max-cells" %in% arg_keys) max_cells <- as.integer(opts[["--max-cells"]])
if (!is.null(max_cells) && !technical_smoke) stop("--max-cells requires --technical-smoke")
if (technical_smoke && is.null(max_cells)) max_cells <- 2000L
if ("--lineage-topology-k" %in% arg_keys)
  CONFIG$llr$lineage_topology_k <- as.integer(opts[["--lineage-topology-k"]])
if ("--root-topology-k" %in% arg_keys)
  CONFIG$llr$root_topology_k <- as.integer(opts[["--root-topology-k"]])
enable_merged_windows <- variant_id %in% c("B", "AB")
enable_knn_topology <- variant_id %in% c("A", "AB") || "--knn-topology" %in% arg_keys
enable_ordering_low_anchor <- "--ordering-low-anchor" %in% arg_keys
enable_terminal_sibling_refine <- "--terminal-sibling-refine" %in% arg_keys
enable_internal_transition_pocket <- "--internal-transition-pocket" %in% arg_keys
enable_timepoint_branchpoint_root <- main_candidate ||
  "--timepoint-branchpoint-root" %in% arg_keys
enable_timepoint_dominant_early_leaf <- main_candidate ||
  "--timepoint-dominant-early-leaf" %in% arg_keys
enable_timepoint_source_flow_guard <- "--timepoint-source-flow-guard" %in% arg_keys
enable_local_branch_residual <- main_candidate || continuous_fork_priming ||
  "--local-branch-residual" %in% arg_keys
enable_root_geometric_gate <- main_candidate || continuous_fork_priming ||
  "--root-geometric-gate" %in% arg_keys
enable_separate_priming_weights <- main_candidate || continuous_fork_priming ||
  "--separate-priming-weights" %in% arg_keys
enable_dropout_aware_event_mean <- main_candidate ||
  "--dropout-aware-event-mean" %in% arg_keys
enable_sparse_overlap_pc1 <- main_candidate ||
  "--sparse-overlap-pc1" %in% arg_keys
enable_preserve_selected_root <- main_candidate ||
  "--preserve-selected-root" %in% arg_keys
enable_structural_event_rank <- main_candidate
enable_connected_long_trunk_orientation <- main_candidate
enable_component_axis_repair <- main_candidate
component_ablation <- Sys.getenv("LLR_COMPONENT_ABLATION", "full")
stopifnot(component_ablation %in% c("full", "no_component_phase_alignment",
  "no_crossing_repair", "no_terminal_continuation", "no_curve_assignment", "no_donor_geometry",
  "no_coherent_clock_signal", "no_sample_temporal_conditioning",
  "no_pulse_evidence", "no_donor_local_conditioning", "descendant_prefork",
  "positive_only_persistence", "strict_common_late_window", "legacy_event_baseline", "pairwise_late_reference", "no_late_bin_merging", "log_activity_contrasts",
  "hill_activity_contrasts", "no_context_error_weighting", "legacy_log_error_weighting", "relative_pulse_only", "local_branch_reference", "supported_branch_reference", "observed_branch_reference", "uniform_branch_reference", "local_curve_reference", "local_point_reference", "assigned_shape_reference", "assigned_layout_reference", "geometric_time_reference"))
use_assignment_layout <- main_candidate && !component_ablation %in% c(
  "local_branch_reference", "supported_branch_reference", "observed_branch_reference",
  "uniform_branch_reference", "local_curve_reference", "local_point_reference",
  "assigned_shape_reference", "geometric_time_reference")
enable_shared_parent <- main_candidate &&
  !component_ablation %in% c("descendant_prefork", "legacy_event_baseline", "legacy_log_error_weighting")
signed_persistence <- main_candidate &&
  !component_ablation %in% c("positive_only_persistence", "legacy_event_baseline", "legacy_log_error_weighting")
enable_late_bin_merging <- main_candidate &&
  !component_ablation %in% c("legacy_event_baseline", "no_late_bin_merging")
enable_shared_induction_clock <- main_candidate &&
  component_ablation != "no_component_phase_alignment"
enable_centroid_phase <- main_candidate &&
  component_ablation != "no_component_phase_alignment"
CONFIG$llr$component_ablation <- component_ablation
CONFIG$llr$lineage_priming_mode <- "geometric_local_envelope"
CONFIG$llr$lineage_assignment_mode <- ifelse(continuous_fork_priming,
  "continuous_local_branch_residual", "local_branch_assignment")
CONFIG$llr$lineage_weight_mode <- switch(component_ablation,
  local_branch_reference = "local_branch_assignment",
  supported_branch_reference = "supported_local_assignment",
  observed_branch_reference = "observed_path_geometry",
  uniform_branch_reference = "observed_path_uniform",
  local_curve_reference = "reference_geometry",
  local_point_reference = "local_branch_assignment", "reference_geometry")
CONFIG$llr$lineage_curve_weight_mode <- switch(as.character(use_assignment_layout),
  "TRUE" = "local_branch_assignment",
  switch(component_ablation,
    local_curve_reference = "local_branch_assignment",
    assigned_shape_reference = "reference_geometry",
    local_point_reference = "reference_geometry", CONFIG$llr$lineage_weight_mode))
switch(component_ablation, local_branch_reference =, supported_branch_reference =, observed_branch_reference =, uniform_branch_reference =, local_curve_reference =, local_point_reference =, assigned_shape_reference =, assigned_layout_reference = {
  CONFIG$llr$result_revision <- paste0(CONFIG$llr$result_revision, "_", component_ablation)
})
switch(as.character(main_candidate), "TRUE" = {
  CONFIG$llr$result_revision <- paste0(CONFIG$llr$result_revision,
    "_current_package_sample_missing_20261007_supported_priming")
})
CONFIG$llr$root_endpoint_terminal_segment <- main_candidate
CONFIG$llr$sample_amplitude_rule <- "equal_weight_within_sample_paired_event"
CONFIG$llr$no_fork_control <-
  "same_joint_temporal_background_without_local_sibling_difference"
CONFIG$llr$prefork_weighting <- ifelse(enable_shared_parent, "shared", "descendant")
CONFIG$llr$persistence_operator <- ifelse(signed_persistence,
  "same_signed_local_and_late_change", "positive_persistence")
CONFIG$llr$persistence_sibling_scope <- switch(component_ablation,
  strict_common_late_window="native",legacy_event_baseline="all_descendants",
  legacy_log_error_weighting="all_descendants",
  pairwise_late_reference="pairwise",
  ifelse(main_candidate,"native","all_descendants"))
CONFIG$llr$late_bin_merge_factor <- 2L
CONFIG$llr$event_activity_transform <- ifelse(main_candidate &&
  !component_ablation %in% c("hill_activity_contrasts", "legacy_event_baseline"), "log_exposure", "hill")
CONFIG$llr$context_error_weighting <- main_candidate &&
  !component_ablation %in% c("no_context_error_weighting", "log_activity_contrasts", "legacy_event_baseline")
CONFIG$llr$bridge_clock_max_gap <- 0.5
CONFIG$llr$long_trunk_orientation_min_ratio <- 2
CONFIG$llr$component_axis_alignment_min <- 0.5
CONFIG$llr$shared_induction_component_r <- 0.25
CONFIG$llr$shared_induction_component_fraction <- 0.75
CONFIG$llr$shared_induction_min_genes <- 20L
stopifnot(!enable_timepoint_dominant_early_leaf || enable_timepoint_branchpoint_root)
stopifnot(!enable_timepoint_source_flow_guard || enable_timepoint_branchpoint_root)
stopifnot(!enable_root_geometric_gate || enable_local_branch_residual)
stopifnot(!enable_separate_priming_weights || enable_root_geometric_gate)
internal_transition_pocket_rule <- list(
  boundary_depletion = 0.1,
  min_time_separation_iqr = 1.0,
  min_fragment_n = 20L,
  min_fragment_fraction = 0.05
)
internal_transition_cluster <- if ("--internal-transition-cluster" %in% arg_keys)
  as.integer(opts[["--internal-transition-cluster"]]) else NULL
stopifnot(!enable_terminal_sibling_refine || enable_knn_topology)
stopifnot(is.null(internal_transition_cluster) || enable_knn_topology)
stopifnot(!enable_internal_transition_pocket || enable_knn_topology)
stopifnot(!enable_internal_transition_pocket || is.null(internal_transition_cluster))
enable_event_mode <- variant$event_scope != "none"
event_scope_val <- variant$event_scope
event_boundary_val <- variant$boundary
event_anchor_val <- variant$anchor_rule
event_scoring_val <- variant$scoring_mode
event_activity_val <- variant$activity_mode
event_expression_mean_mode <- NA_character_
lr_detection_intensity_rho <- NA_real_
CONFIG$llr$min_n_eff <- variant$min_n_eff
if (variant_id == "L1D0contSparsePearsonDiffusion")
  CONFIG$llr$sparse_diffusion_dims <- 10L
if (variant_id %in% c("L1D0contAxisCoherent",
                      "L1D0contRootBridgeLocalRank",
                      "L1D0contTypeCoherentAnchor",
                      "L1D0contLowFlowLeaf",
                      "L1D0contBranchSource",
                      "L1D0contSampleRankRoot",
                      "L1D0contSampleRankLowRoot",
                      "L1D0contSampleDominanceRoot",
                      "L1D0contTargeted",
                       "L1D0contTargetedInterval",
                       "L1D0contTargetedDonorAdaptiveInterval",
                       "L1D0contDonorStructuralInterval",
                       "L1D0contDualAttractorInterval",
                      "L1D0contDepthResidualInterval",
                      "L1D0contDepthResidualLeafInterval",
                      "L1D0contDepthResidualNoDirectInterval",
                      "L1D0contDepthResidualOrderedLeafInterval")) {
  CONFIG$llr$root_axis_min_alignment <- 0.80
  CONFIG$llr$root_axis_min_gain <- 0.50
}
CONFIG$llr$root_max_iter <- fcoalesce(as.integer(variant$root_max_iter),
                                    as.integer(CONFIG$llr$root_max_iter))
if (variant$time_method %in% c("endpoint_anchor", "terminal_distance",
                               "bridge_terminal_distance",
                               "bridge_terminal_direct_root",
                               "raw_multibridge_terminal_direct",
                               "tangent_multibridge_terminal_direct")) {
  CONFIG$llr$root_endpoint_alignment_min_gain <- 0.10
}
CONFIG$llr$result_revision <- paste0(CONFIG$llr$result_revision,
  "_root_time_from_start",
  switch(variant_id,
    L1D0contBridgeFloor = "_bridge_floor_cost",
    L1D0contCombined = "_bridge_terminal_distance",
    L1D0contCombinedDirectRoot = "_bridge_terminal_direct_root",
    L1D0contMetricConsistent = "_raw_multibridge_terminal_direct_root",
    L1D0contTangentBridge = "_tangent_multibridge_terminal_direct_root",
    L1D0contTangentLocalRank = "_tangent_multibridge_component_local_rank",
    L1D0contStarLocalRank = "_fixed_main_bridge_component_local_rank",
    L1D0contStarEarlyLocalRank = "_early_fixed_main_bridge_component_local_rank",
    L1D0contStarEarlyPotencyRank = "_early_fixed_main_bridge_component_potency_rank",
    L1D0contStarCorePotencyRank = "_core_fixed_main_bridge_component_potency_rank",
    L1D0contAxisCoherent = "_axis_coherent_root_score_tangent_multibridge",
    L1D0contRootBridgeLocalRank = "_axis_root_tangent_bridge_local_rank",
    L1D0contTypeCoherentAnchor = "_type_coherent_anchor_terminal_topology_tangent_bridge_local_rank",
    L1D0contLowFlowLeaf = "_lowstate_max_source_flow_leaf_root",
    L1D0contBranchSource = "_branch_local_source_tangent_bridge_local_rank",
    L1D0contSampleRankRoot = "_within_sample_detected_rank_root_tangent_bridge_local_rank",
    L1D0contSampleRankLowRoot = "_within_sample_detected_rank_low_root_tangent_bridge_local_rank",
    L1D0contSampleDominanceRoot = "_sample_dominance_detected_rank_root_tangent_bridge_local_rank",
    L1D0contTargeted = "_targeted_root_and_branch_source_tangent_bridge_local_rank",
    L1D0contTargetedInterval = "_targeted_source_leaf_segment_interval",
    L1D0contTargetedDonorAdaptiveInterval = "_targeted_donor_adaptive_source_leaf_segment_interval",
    L1D0contDonorStructuralInterval = "_donor_adaptive_graph_spur_segment_interval",
    L1D0contDualAttractorInterval = "_dual_attractor_source_leaf_segment_interval",
    L1D0contDepthResidualInterval = "_depth_residual_root_segment_interval",
    L1D0contDepthResidualLeafInterval = "_depth_residual_source_leaf_segment_interval",
    L1D0contDepthResidualNoDirectInterval = "_depth_residual_no_direct_segment_interval",
    L1D0contDepthResidualOrderedLeafInterval = "_depth_residual_ordered_leaf_segment_interval",
    L1D0contSparsePC10 = "_sparse_pc10_bridge_terminal_distance",
    L1D0contSparsePearson = "_sparse_pearson_pc10_bridge_terminal_distance",
    L1D0contSparsePearsonDiffusion = "_sparse_pearson_pc10_diffusion_time",
    L1D0contEndpointAnchored = "_terminal_endpoint_anchor",
    L1D0contTerminalDistance = "_terminal_distance_oriented",
    L1D0contPrecise = "_fuzzy_f64_accumulation",
    L1D0contFateSolve = "_ilu0_gmres_candidate_availability",
    L1D0contRetention = "_ilu0_gmres_terminal_retention",
    L1D0contRetentionNoLRGeometry = "_terminal_retention_lr_independent_geometry",
    L1D0contRetentionSourceLeaf = "_retention_source_dominant_leaf_rooted_time",
    L1D0contRetentionSourceLeafInterval = "_retention_source_leaf_segment_interval",
    L1D0contRetentionLu = "_sparse_lu_terminal_retention", ""))
CONFIG$llr$result_revision <- paste0(CONFIG$llr$result_revision,
  switch(as.character(run_tag %in% c("L1D0contBranchBalanceInterval",
                                      "L1D0contMaturityMatchedInterval",
                                      "L1D0contDonorSplitGeometryInterval",
                                      "L1D0contDonorBranchBalanceInterval",
                                      "L1D0contDonorBalanceNoSpurInterval")),
    "TRUE" = "_branch_balance_fraction", "FALSE" = ""))
CONFIG$llr$result_revision <- paste0(CONFIG$llr$result_revision,
  switch(as.character(run_tag == "L1D0contMaturityMatchedInterval"),
    "TRUE" = "_time_matched_sibling_terminal_maturity_factor_source_branch_compatibility", "FALSE" = ""))
CONFIG$llr$result_revision <- paste0(CONFIG$llr$result_revision,
  switch(as.character(run_tag == "L1D0contDonorBalanceNoSpurInterval"),
    "TRUE" = "_retain_backward_tips", "FALSE" = ""))
CONFIG$llr$result_revision <- paste0(CONFIG$llr$result_revision,
  switch(as.character(donor_split_geometry_mode),
    "TRUE" = "_donor_split_geometry", "FALSE" = ""))
CONFIG$llr$result_revision <- paste0(CONFIG$llr$result_revision,
  switch(as.character(donor_center_lr_mode),
    "TRUE" = "_donor_center_lr_geometry", "FALSE" = ""))
if (enable_ordering_low_anchor)
  CONFIG$llr$result_revision <- paste0(CONFIG$llr$result_revision,
    "_ordering_low_anchor")
if (enable_terminal_sibling_refine)
  CONFIG$llr$result_revision <- paste0(CONFIG$llr$result_revision,
    "_terminal_sibling_refine")
if (enable_internal_transition_pocket)
  CONFIG$llr$result_revision <- paste0(CONFIG$llr$result_revision,
    "_internal_transition_pocket_window")
if (!is.null(internal_transition_cluster))
  CONFIG$llr$result_revision <- paste0(CONFIG$llr$result_revision,
    "_internal_transition_cluster", internal_transition_cluster)
CONFIG$llr$result_revision <- paste0(CONFIG$llr$result_revision,
  switch(as.character(enable_local_branch_residual),
    "TRUE" = "_local_branch_residual_scaled", "FALSE" = ""))
CONFIG$llr$result_revision <- paste0(CONFIG$llr$result_revision,
  switch(as.character(enable_root_geometric_gate),
    "TRUE" = "_root_geometric_gate", "FALSE" = ""))
CONFIG$llr$result_revision <- paste0(CONFIG$llr$result_revision,
  switch(as.character(enable_separate_priming_weights),
    "TRUE" = "_separate_priming_weights", "FALSE" = ""))
CONFIG$llr$result_revision <- paste0(CONFIG$llr$result_revision,
  switch(as.character(enable_sparse_overlap_pc1),
    "TRUE" = "_sparse_overlap_pc1_edge_cost", "FALSE" = ""))
CONFIG$llr$result_revision <- paste0(CONFIG$llr$result_revision,
  switch(as.character(enable_preserve_selected_root),
    "TRUE" = "_preserve_selected_root", "FALSE" = ""))
CONFIG$llr$result_revision <- paste0(CONFIG$llr$result_revision,
  switch(as.character(enable_structural_event_rank),
    "TRUE" = "_structural_event_rank", "FALSE" = ""))
CONFIG$llr$result_revision <- paste0(CONFIG$llr$result_revision,
  switch(as.character(enable_connected_long_trunk_orientation),
    "TRUE" = "_connected_long_trunk_orientation", "FALSE" = ""))
CONFIG$llr$result_revision <- paste0(CONFIG$llr$result_revision,
  switch(as.character(enable_component_axis_repair),
    "TRUE" = "_component_axis_repair", "FALSE" = ""))
CONFIG$llr$result_revision <- paste0(CONFIG$llr$result_revision,
  switch(as.character(enable_shared_induction_clock),
    "TRUE" = "_shared_induction_clock", "FALSE" = ""))
CONFIG$llr$result_revision <- paste0(CONFIG$llr$result_revision,
  switch(as.character(enable_centroid_phase),
    "TRUE" = "_centroid_phase_pca2", "FALSE" = ""))
root_candidate_options <- switch(variant_id,
  L1D0contPrecise = list(on_unreachable = "return", fate_solver = "sparse_lu"),
  L1D0contTargeted = list(on_unreachable = "return"),
  L1D0contTargetedInterval = list(on_unreachable = "return"),
  L1D0contTargetedDonorAdaptiveInterval = list(on_unreachable = "return"),
  L1D0contDonorStructuralInterval = list(on_unreachable = "return"),
  L1D0contDualAttractorInterval = list(on_unreachable = "return"),
  L1D0contDepthResidualInterval = list(on_unreachable = "return"),
  L1D0contDepthResidualLeafInterval = list(on_unreachable = "return"),
  L1D0contDepthResidualNoDirectInterval = list(on_unreachable = "return"),
  L1D0contDepthResidualOrderedLeafInterval = list(on_unreachable = "return"),
  L1D0contFateSolve = list(on_unreachable = "return"),
  L1D0contRetention = list(on_unreachable = "return", attractor_selection = "retention"),
  L1D0contRetentionNoLRGeometry = list(on_unreachable = "return",
    attractor_selection = "retention"),
  L1D0contRetentionSourceLeaf = list(on_unreachable = "return",
    attractor_selection = "retention"),
  L1D0contRetentionSourceLeafInterval = list(on_unreachable = "return",
    attractor_selection = "retention"),
  L1D0contRetentionLu = list(on_unreachable = "return",
    attractor_selection = "retention", fate_solver = "sparse_lu"),
  list())
if (enable_event_mode) {
  stopifnot(all(c("activity_mode", "scoring_mode", "bif_df") %in%
                  names(formals(LineaLinkR::ComputeForkContrastsEvent))))
}
CONFIG$llr$result_revision <- switch(as.character(main_candidate),
  "TRUE" = paste0("linealinkr065_within_sample_fixed_sisters",
    ifelse(use_assignment_layout, "_assigned_layout_reference", ""),
    ifelse(component_ablation %in% c("full", "assigned_layout_reference"), "",
      paste0("_", component_ablation))),
  "FALSE" = CONFIG$llr$result_revision)
lr_abundance_mode <- Sys.getenv("LLR_LR_ABUNDANCE", "log1p")
switch(as.character(lr_abundance_mode %in% c("linear", "sqrt", "positive_reference")), "TRUE" = {
  linear_datasets <- unique(c(names(CONFIG$dataset_registry),
    names(CONFIG$dataset_narrative_registry), filter_datasets))
  provided_datasets <- names(CONFIG$dataset_lr_input_scale)[
    CONFIG$dataset_lr_input_scale == "provided_expression"]
  CONFIG$dataset_lr_input_scale[linear_datasets] <- ifelse(
    linear_datasets %in% provided_datasets, "provided_expression",
    switch(lr_abundance_mode, linear = "counts_linear", sqrt = "counts_sqrt",
      positive_reference = "counts_positive_reference"))
  CONFIG$llr$result_revision <- paste0(CONFIG$llr$result_revision, "_", lr_abundance_mode, "_abundance")
})
lineage_state_column <- Sys.getenv("LLR_LINEAGE_STATE_COLUMN", "")
switch(as.character(enable_knn_topology), "TRUE" = {
  CONFIG$llr$result_revision <- paste0(CONFIG$llr$result_revision, "_observed_graph_context",
    ifelse(nzchar(lineage_state_column), paste0("_", lineage_state_column), ""))
})
output_scope <- switch(as.character(main_candidate),
  "TRUE" = switch(as.character(component_ablation == "full"),
    "TRUE" = "final065",
    "FALSE" = file.path("final065_ablation", component_ablation)),
  "FALSE" = file.path("research", variant_id))
OUT <- file.path(CONFIG$dir$cache, Sys.getenv("LLR_OUTPUT_WORKSPACE", output_scope))
if ("--validate-only" %in% arg_keys) {
  invisible(vapply(filter_datasets,DatasetPath,character(1L)))
  cat("variant=", run_tag, "; workspace=", OUT,
      "; datasets=", paste(filter_datasets, collapse = ","),
      "; library=", find.package("LineaLinkR"),
      "; root_max_iter=", CONFIG$llr$root_max_iter, "\n", sep = "")
  quit(save = "no", status = 0L)
}
stopifnot(!any(file.exists(file.path(OUT, filter_datasets, "linealinkr_complete.rds"))))

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

LlrDsDir <- function(ds_name) {
  d <- file.path(OUT, ds_name)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}
LlrExtractDir <- function(ds_name) {
  d <- file.path(switch(as.character(nzchar(precompute_workspace)),
    "TRUE" = precompute_workspace, "FALSE" = OUT), ds_name, "extract")
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

set.seed(CONFIG$seed)
Sys.setenv(OMP_NUM_THREADS = CONFIG$n_threads,
  OPENBLAS_NUM_THREADS = CONFIG$n_threads, MKL_NUM_THREADS = CONFIG$n_threads)
data.table::setDTthreads(CONFIG$n_threads)

cat(sprintf("OMP threads: %s | seed: %d\n",
    Sys.getenv("OMP_NUM_THREADS"), CONFIG$seed))

# =============================================================================
# Helper functions (trusted RDS contract: protocol, revision, and core config).
# =============================================================================
LlrExpectedMeta <- function(knn_topology = enable_knn_topology) {
  cfg <- CONFIG$llr %||% list()
  expected <- list(
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
    root_endpoint_terminal_segment = isTRUE(
      cfg$root_endpoint_terminal_segment
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
    fork_sample_boot = as.integer(cfg$fork_sample_boot %||% 999L),
    # Algorithm-revision variants: cache trust keys on the actual switches.
    run_tag = as.character(run_tag),
    event_anchor_rule = event_anchor_val,
    event_scoring_mode = event_scoring_val,
    event_activity_mode = event_activity_val,
    merge_fork_windows = isTRUE(enable_merged_windows),
    knn_graph_topology = isTRUE(knn_topology)
  )
  if (length(cfg$user_root_cells)) expected$user_root_cells <- cfg$user_root_cells
  switch(as.character(resident_sender), "TRUE" = {
    expected$resident_context_column <- resident_stage_column
    expected$resident_time_column <- resident_time_column
    expected$resident_condition_column <- resident_condition_column
  })
  switch(as.character(persistent_rna_states), "TRUE" = {
    expected$state_time_column <- state_time_column
    expected$state_condition_column <- state_condition_column
  })
  if ("--root-topology-k" %in% arg_keys) {
    expected$root_topology_k <- as.integer(cfg$root_topology_k)
    expected$root_topology_role <- "selector_and_branch_source"
  }
  expected
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
  selected_blind_candidates <-
    identical(x$config_snapshot$root_strategy, "attractor_energy") &&
    identical(x$root_potency_mode, "detected_genes_ordering") &&
    !(x$root_evidence_mode %in% c("independent_user_root", "observed_root_prior"))
  root_selection_ok <- if (selected_blind_candidates) {
    is.data.frame(x$root_selection) && nrow(x$root_selection) == 4L &&
      all(c("terminal_anchor_is_high_end", "terminal_anchor_is_leaf",
        "terminal_leaf_eligible") %in% names(x$root_selection)) &&
      identical(as.logical(x$root_selection$terminal_leaf_eligible),
        (as.logical(x$root_selection$terminal_anchor_is_high_end) &
          as.logical(x$root_selection$terminal_anchor_is_leaf)) %in% TRUE)
  } else {
    is.null(x$root_selection) || is.data.frame(x$root_selection)
  }
  isTRUE(length(n) == 1L && is.finite(n) && n > 0L) &&
    is.character(x$dataset) && length(x$dataset) == 1L && nzchar(x$dataset) &&
    is.character(x$cells) && length(x$cells) == n && !anyNA(x$cells) && uniqueN(x$cells) == n &&
    is.numeric(x$pseudotime) && length(x$pseudotime) == n && any(is.finite(x$pseudotime)) &&
    is.matrix(x$lineage_weights) && nrow(x$lineage_weights) == n && ncol(x$lineage_weights) >= 1L &&
    all(is.finite(x$lineage_weights)) &&
    (is.null(x$lr_scores) || is.data.frame(x$lr_scores)) &&
    (is.null(x$temporal_lr) || (is.list(x$temporal_lr) && is.data.frame(x$temporal_lr$scores) && is.list(x$temporal_lr$intermediate))) &&
    identical(as.character(x$package_version), "0.6.5") &&
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
    root_selection_ok &&
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

LlrRdsIsTrusted <- function(x, expected_manifest, expected_cells,
                            knn_topology = enable_knn_topology) {
  if (is.null(x) || !is.list(x)) return(FALSE)
  exp <- LlrExpectedMeta(knn_topology)
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
  expected_lr_scale <- unname(ifelse(
    expected_manifest$dataset %in% names(CONFIG$dataset_lr_input_scale),
    CONFIG$dataset_lr_input_scale[expected_manifest$dataset], "counts"))
  cluster_input_policy_ok &&
    identical(x$lr_input_scale %||% "counts", expected_lr_scale) &&
    all(vapply(names(exp), function(nm) LlrMetaSame(saved[[nm]], exp[[nm]]), logical(1L))) &&
    LlrRdsComplete(x) &&
    identical(as.character(x$cells), as.character(expected_cells)) &&
    LlrInputManifestTrusted(x$input_manifest, expected_manifest)
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
exclusions = pd.read_csv(r'%s', encoding='utf-8')
excluded_genes = exclusions.loc[exclusions['dataset'].eq(dsname), 'gene']
adata = adata[:, ~adata.var_names.isin(excluded_genes)].copy()

# Expression matrix: use the deposited count layer; keep sparse.
if 'counts' in adata.layers:
    mat = adata.layers['counts']
    matrix_source = 'layers/counts'
elif 'raw' in adata.layers:
    mat = adata.layers['raw']
    matrix_source = 'layers/raw'
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
if 'Sample ID' in meta.columns and 'sample_id' not in meta.columns:
    meta['sample_id'] = meta['Sample ID'].astype(str)
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
with open(os.path.join(outdir, f'{dsname}_matrix.mtx'), 'wb') as stream:
    mmwrite(stream, mat)

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
)", h5, LlrExtractDir(ds_name), ds_name, file.path(CONFIG$dir$pipeline, "config/feature_exclusions.csv"))

  ext <- LlrExtractDir(ds_name)
  mtx_path <- file.path(ext, sprintf("%s_matrix.mtx", ds_name))
  meta_path <- file.path(ext, sprintf("%s_meta.csv", ds_name))
  genes_path <- file.path(ext, sprintf("%s_genes.txt", ds_name))
  cells_path <- file.path(ext, sprintf("%s_cells.txt", ds_name))
  matrix_source_path <- file.path(ext, sprintf("%s_matrix_source.txt", ds_name))
  manifest_path <- file.path(ext, sprintf("%s_input_content_manifest.rds", ds_name))
  switch(as.character(nzchar(precompute_workspace)), "FALSE" = {
    unlink(c(mtx_path, meta_path, genes_path, cells_path, matrix_source_path, manifest_path), force = TRUE)
    extract_status <- sys::exec_wait(CONFIG$exe$python, args = c("-X", "utf8", "-c", py),
      std_out = file.path(CONFIG$dir$log, sprintf("llr_extract_%s.log", ds_name)),
      std_err = file.path(CONFIG$dir$log, sprintf("llr_extract_%s.err", ds_name)))
    if (extract_status != 0L) stop(sprintf("Python h5ad extraction failed for %s", ds_name))
  })

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
  switch(as.character(nzchar(precompute_workspace)),
    "TRUE" = stopifnot(LlrInputManifestTrusted(readRDS(manifest_path), input_manifest)),
    "FALSE" = saveRDS(input_manifest, manifest_path))

  list(expr = expr, meta = meta, genes = genes, cells = cells, input_manifest = input_manifest)
}

# =============================================================================
# Run the complete LineaLinkR workflow for one dataset.
# =============================================================================
# Blind cross-branch clock calibration. Eligible sibling subtrees share an
# early-segment IQR. The default affine map fixes the common fork median;
# translating subtree medians would move that shared anchor. Subtrees whose
# early median is disjoint from the pool (> gate*ref) are left untouched.
# No ground truth is used.
LlrBuildBranchForks <- function(paths) {
  forks <- list()
  allcl <- unique(unlist(paths))
  for (X in allcl) {
    users <- which(vapply(paths, function(p) X %in% p, TRUE))
    if (length(users) < 2) next
    pos <- vapply(users, function(b) match(X, paths[[b]]), 0)
    nxt <- mapply(function(b, i) if (i < length(paths[[b]])) paths[[b]][i + 1]
      else NA_integer_, users, pos)
    kids <- unique(nxt[!is.na(nxt)])
    if (length(kids) < 2) next
    post <- lapply(kids, function(k) unique(unlist(lapply(users[nxt == k],
      function(b) {
        p <- paths[[b]]
        p[(pos[which(users == b)] + 1):length(p)]
      }))))
    forks[[length(forks) + 1]] <- list(F = X, depth = min(pos),
      kids = kids, post = post)
  }
  forks[order(vapply(forks, function(fk) fk$depth, 0))]
}

LlrCalibrateBranchClock <- function(w, cl, paths, gate = 0.6, loc = 0.0,
    min_n = 20L, clamp = c(0.2, 8)) {
  forks <- LlrBuildBranchForks(paths)
  n_cal <- 0L
  for (fk in forks) {
    tF <- median(w[cl == fk$F], na.rm = TRUE)
    E <- lapply(fk$kids, function(k) which(cl == k))
    if (any(vapply(E, length, 0) < min_n))
      E <- lapply(fk$post,
        function(pc) which(cl %in% pc[seq_len(min(2L, length(pc)))]))
    if (any(vapply(E, length, 0) < min_n)) next
    idxs <- lapply(fk$post, function(pc) which(cl %in% pc))
    sc <- vapply(E, function(ix) IQR(w[ix], na.rm = TRUE), 0)
    sc[sc <= 0] <- NA_real_
    if (any(!is.finite(sc))) next
    ref <- median(sc)
    em <- vapply(E, function(ix) median(w[ix], na.rm = TRUE), 0)
    pem <- median(w[unlist(E)], na.rm = TRUE)
    for (s in seq_along(idxs)) {
      if (abs(em[s] - pem) > gate * ref) next
      a <- min(max(ref / sc[s], clamp[1]), clamp[2])
      idx <- idxs[[s]]
      w[idx] <- tF + a * (w[idx] - tF) + loc * (pem - em[s])
      n_cal <- n_cal + 1L
    }
  }
  attr(w, "n_subtrees_calibrated") <- n_cal
  w
}

RunLineaLinkROnDataset <- function(ds_name) {
  cat(sprintf("\n===== LineaLinkR on %s =====\n", ds_name))
  result_dir <- LlrDsDir(ds_name)
  log_file <- file.path(CONFIG$dir$log, sprintf("llr_%s.log", ds_name))

  switch(as.character(nzchar(anchor_workspace)), "TRUE" = {
    source(file.path(.llr_work_root, "experiments/run/complete_main_real_candidate.R"),
      local = TRUE, encoding = "UTF-8")
    return(x)
  })
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

  geometry_excludes_lr <- variant_id == "L1D0contRetentionNoLRGeometry"
  geometry_centers_lr <- FALSE
  geometry_donor_signal <- NULL
  if (component_ablation != "no_donor_geometry" &&
      (variant_id %in% c("L1D0contRetentionNoLRGeometry",
                         "L1D0contTargetedDonorAdaptiveInterval",
                         "L1D0contDonorStructuralInterval") ||
      donor_split_geometry_mode || donor_center_lr_mode)) {
    geometry_db <- BuildLRDB(expr, species = CONFIG$dataset_species[[ds_name]] %||% "human")
    lr_rows <- unique(unlist(c(geometry_db$ligand_rows,
      geometry_db$receptor_rows), use.names = FALSE))
    if (variant_id %in% c("L1D0contTargetedDonorAdaptiveInterval",
                          "L1D0contDonorStructuralInterval") ||
        donor_split_geometry_mode || donor_center_lr_mode) {
      sample_col <- FindBiologicalSampleCol(meta, ds_name)
      if (!is.na(sample_col) && length(lr_rows) > 0L) {
        sample_labels <- as.character(meta[[sample_col]])
        sample_keep <- which(!is.na(sample_labels) & nzchar(trimws(sample_labels)))
        donor_code <- as.integer(factor(sample_labels[sample_keep]))
        donor_n <- tabulate(donor_code)
        normalized <- PrepareLRExpression(expr[, sample_keep, drop = FALSE])
        donor_membership <- sparseMatrix(
          i = seq_along(donor_code), j = donor_code, x = 1,
          dims = c(length(donor_code), length(donor_n))
        )
        overall_mean <- Matrix::rowMeans(normalized)
        donor_mean <- sweep(as.matrix(normalized %*% donor_membership),
                            2L, donor_n, "/")
        between <- rowSums(sweep(
          sweep(donor_mean, 1L, overall_mean, "-") ^ 2,
          2L, donor_n / length(donor_code), "*"
        ))
        total <- Matrix::rowMeans(normalized ^ 2) - overall_mean ^ 2
        donor_r2 <- between / total
        donor_r2[!is.finite(donor_r2)] <- 0
        hvg_rows <- head(order(total, decreasing = TRUE), 2000L)
        nonlr_hvg_rows <- setdiff(hvg_rows, lr_rows)
        background_upper <- unname(quantile(
          donor_r2[nonlr_hvg_rows], 0.99
        ))
        strong_lr_donor_signal <- median(donor_r2[lr_rows]) > background_upper
        geometry_excludes_lr <- strong_lr_donor_signal && !donor_center_lr_mode
        geometry_centers_lr <- strong_lr_donor_signal && donor_center_lr_mode &&
          grepl("donor|patient", sample_col, ignore.case = TRUE)
        geometry_donor_signal <- list(
          sample_column = sample_col,
          lr_median_r2 = median(donor_r2[lr_rows]),
          nonlr_hvg_q99_r2 = background_upper,
          excludes_lr = geometry_excludes_lr,
          centers_lr = geometry_centers_lr
        )
        cat(sprintf("  [LR donor R2 %.4f; non-LR HVG q99 %.4f; exclude LR %s; center LR %s]\n",
          geometry_donor_signal$lr_median_r2,
          geometry_donor_signal$nonlr_hvg_q99_r2,
          geometry_excludes_lr, geometry_centers_lr))
        rm(normalized, donor_membership, donor_mean)
      }
    }
  }
  use_knn_topology <- enable_knn_topology ||
    (variant_id == "L1D0contDonorStructuralInterval" && geometry_excludes_lr) ||
    (donor_split_geometry_mode && geometry_excludes_lr) ||
    geometry_centers_lr

  # Every execution starts from an empty output directory.
  cache_file <- file.path(result_dir, "linealinkr_complete.rds")
  run_attempt <- 1L

  ext <- LlrExtractDir(ds_name)
  knn_matrix_path <- file.path(ext, sprintf("%s_matrix.mtx", ds_name))
  if (!is.null(max_cells) && ncol(expr) < ncol(dat$expr)) {
    knn_matrix_path <- file.path(result_dir, sprintf("%s_matrix_downsampled.mtx", ds_name))
    writeMM(expr, knn_matrix_path)
    writeLines(rownames(expr), file.path(result_dir, sprintf("%s_genes.txt", ds_name)), useBytes = TRUE)
  }
  if (variant_id %in% c("L1D0contRetentionNoLRGeometry",
                        "L1D0contTargetedDonorAdaptiveInterval",
                        "L1D0contDonorStructuralInterval")) {
    if (geometry_excludes_lr) {
      geometry_rows <- setdiff(seq_len(nrow(expr)), lr_rows)
      knn_matrix_path <- file.path(result_dir, sprintf("%s_geometry_matrix.mtx", ds_name))
      writeMM(expr[geometry_rows, , drop = FALSE], knn_matrix_path)
      writeLines(rownames(expr)[geometry_rows],
        file.path(result_dir, sprintf("%s_genes.txt", ds_name)), useBytes = TRUE)
      cat(sprintf("  [geometry excludes %d LR genes]\n", length(lr_rows)))
    }
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
  requested_ct_col <- unname(CONFIG$dataset_celltype_columns[ds_name])
  ct_candidates <- unique(c(requested_ct_col[!is.na(requested_ct_col)], ct_candidates))
  ct_candidates <- ct_candidates[!force_auto_clusters]
  ct_col <- NA_character_
  for (cn in ct_candidates) {
    if (IsValidClusterCol(meta, cn, min_groups = ifelse(cn %chin% requested_ct_col, 1L, 2L))) { ct_col <- cn; break }
  }
  if (is.na(ct_col)) {
    cluster_cols <- grep("cluster", names(meta), ignore.case = TRUE, value = TRUE)
    cluster_cols <- cluster_cols[!force_auto_clusters]
    for (cn in cluster_cols) {
      if (IsValidClusterCol(meta, cn, min_groups = ifelse(cn %chin% requested_ct_col, 1L, 2L))) { ct_col <- cn; break }
    }
  }
  if (is.na(ct_col)) {
    cl_csv <- file.path(result_dir, sprintf("%s_auto_clusters.csv", ds_name))
    cluster_manifest_path <- file.path(result_dir, "auto_cluster_manifest.rds")
    expected_cluster_manifest <- c(
      dat$input_manifest,
      list(
        cluster_revision = "linealinkr065_scanpy_leiden_pca30_neighbors30_igraph_int64",
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
# igraph requests unsigned 32-bit draws; NumPy defaults to int32 on Windows.
sc._utils.RNGIgraph.randint = lambda rng, low, high: int(rng._rng.randint(low, high, dtype=np.int64))
with open(r'%s', 'rb') as stream:
    a = sc.AnnData(mmread(stream).tocsr().T)
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
  llr_cfg <- if (!is.null(CONFIG$llr)) CONFIG$llr else list(
    start_cell_method = "mst_extremum", rank_normalize = FALSE,
    knn_neighbors = 30L, augment_mst_bridges = TRUE,
    root_partition_max_clusters = 8L
  )
  # Auto-refine the lineage partition via kmeans on the UMAP embedding
  # (computed later in Step 3) when the input clustering is too coarse for MST
  # branching.  The cut-off is an empirical operating point, not a derived
  # quantity.  Flag for deferred refinement - actual kmeans runs after UMAP.
  need_kmeans_refine <-
    (n_lineage_cl <= as.integer(llr_cfg$root_partition_max_clusters))
  cat(sprintf("  Lineage clusters: %d (from '%s')%s\n",
    n_lineage_cl, fine_col,
    if(need_kmeans_refine) " [will refine via kmeans after UMAP]" else ""))

  knn_k <- as.integer(min(llr_cfg$knn_neighbors %||% 30L, max(2L, ncol(expr) - 1L)))

  # Step 3: compute pseudotime.
  cat("  Step 3/7: Pseudotime (fuzzy Dijkstra) ... ")
  # Scanpy-equivalent normalization, Seurat-HVG, PCA, and neighbour backends
  # without importing the full Scanpy stack. UMAP is generated only for figures.
  knn_script <- file.path(
    CONFIG$dir$pipeline, "run",
    "prepare_knn.py"
  )
  knn_dir <- file.path(switch(as.character(nzchar(precompute_workspace)),
    "TRUE" = precompute_workspace, "FALSE" = OUT), ds_name)
  knn_manifest_path <- file.path(knn_dir, "knn_precompute_manifest.rds")
  expected_knn_manifest <- c(
    dat$input_manifest,
    list(
      preprocessor_contract = switch(variant_id,
        L1D0contRetentionNoLRGeometry = "taskM_pca30_hvg2000_knn_umap_lr_independent_geometry",
        L1D0contTargetedDonorAdaptiveInterval = c(
          "taskM_pca30_hvg2000_knn_umap",
          "taskM_pca30_hvg2000_knn_umap_lr_independent_geometry"
        )[1L + as.integer(geometry_excludes_lr)],
        L1D0contDonorStructuralInterval = c(
          "taskM_pca30_hvg2000_knn_umap",
          "taskM_pca30_hvg2000_knn_umap_lr_independent_geometry"
        )[1L + as.integer(geometry_excludes_lr)],
        L1D0contSparsePC10 = "taskM_pca10_hvg2000_knn_umap",
        L1D0contSparsePearson = "taskM_pearson_pca10_knn_umap",
        L1D0contSparsePearsonDiffusion = "taskM_pearson_pca10_knn_umap",
        "taskM_pca30_hvg2000_knn_umap"),
      seed = CONFIG$seed,
      knn_neighbors = knn_k,
      write_umap = TRUE
    )
  )
  knn_required <- file.path(
    knn_dir,
    c("knn_dist.mtx", "knn_conn.mtx", "pca.csv", "umap.csv")
  )
  switch(as.character(nzchar(precompute_workspace)),
    "TRUE" = stopifnot(LlrInputManifestTrusted(readRDS(knn_manifest_path), expected_knn_manifest)),
    "FALSE" = {
  unlink(knn_required, force = TRUE)
  knn_status <- sys::exec_wait(
    CONFIG$exe$python,
    args = c(
      "-X", "utf8", knn_script, ds_name, result_dir,
      knn_matrix_path, as.character(knn_k), "with_umap",
      if (variant_id %in% c("L1D0contSparsePC10",
                            "L1D0contSparsePearson",
                            "L1D0contSparsePearsonDiffusion")) "10",
      if (variant_id %in% c("L1D0contSparsePearson",
                            "L1D0contSparsePearsonDiffusion")) "pearson"
    ),
    std_out = log_file,
    std_err = file.path(CONFIG$dir$log, sprintf("llr_%s.err", ds_name)))
  if (knn_status != 0L) {
    stop(sprintf("Python kNN computation failed for %s", ds_name))
  }
  })

  # Decode every output before the preprocessor contract is recorded.
  knn_dist_path <- file.path(knn_dir, "knn_dist.mtx")
  knn_conn_path <- file.path(knn_dir, "knn_conn.mtx")

  if (!all(file.exists(knn_required))) {
    stop(sprintf("kNN computation did not write required MatrixMarket outputs for %s", ds_name))
  }

  knn_dist <- as(readMM(knn_dist_path), "CsparseMatrix")
  knn_conn <- as(readMM(knn_conn_path), "CsparseMatrix")
  pca_embedding <- as.matrix(fread(file.path(knn_dir, "pca.csv"), header = FALSE))
  umap <- as.matrix(fread(file.path(knn_dir, "umap.csv"), header = FALSE))
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
  switch(as.character(nzchar(precompute_workspace)),
    "FALSE" = saveRDS(expected_knn_manifest, knn_manifest_path))

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
    if (!nzchar(Sys.getenv("LLR_NO_FIX_KNN_SELF"))) {
      cand <- s + seq_len(min(k + 1L, e - s + 1L)) - 1L
      keep <- cand[sm_j[cand] != i]
      take <- seq_len(min(k, length(keep)))
      nn_idx[i, take] <- sm_j[keep[take]]
      nn_dists[i, take] <- sm_x[keep[take]]
    } else {
      take <- seq_len(min(k, e - s + 1L))
      nn_idx[i, take] <- sm_j[s + take - 1L]
      nn_dists[i, take] <- sm_x[s + take - 1L]
    }
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
    as.integer(llr_cfg$root_topology_k %||% llr_cfg$lineage_topology_k %||% 15L),
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
  bridge_component_labels <- NULL
  native_component_count <- NA_integer_
  if (variant$time_method %in% c("bridge_floor_cost",
                                 "bridge_terminal_distance",
                                 "bridge_terminal_direct_root",
                                 "raw_multibridge_terminal_direct",
                                 "tangent_multibridge_terminal_direct")) {
    bridge_component_labels <- LineaLinkR:::KnnComponentLabels(nn_idx)
  }
  if (isTRUE(llr_cfg$connect_components)) {
    bridge_candidates <- NULL
    if (variant_id %in% c("L1D0contStarEarlyLocalRank",
                          "L1D0contStarEarlyPotencyRank",
                          "L1D0contStarCorePotencyRank")) {
      bridge_depth <- data.table(
        component = bridge_component_labels,
        detected = as.numeric(Matrix::colSums(expr > 0))
      )
      bridge_depth[, eligible :=
        data.table::frank(-detected, ties.method = "first") <=
          max(1L, ceiling(.N * (llr_cfg$root_top_fraction %||% 0.05))),
        by = component]
      if (variant_id == "L1D0contStarCorePotencyRank") {
        bridge_edges <- which(nn_idx > 0L, arr.ind = TRUE)
        bridge_graph <- igraph::graph_from_data_frame(
          data.frame(from = as.character(bridge_edges[, 1L]),
                     to = as.character(nn_idx[bridge_edges])),
          directed = TRUE,
          vertices = data.frame(name = as.character(seq_len(n_cells))))
        bridge_strong <- igraph::components(bridge_graph, mode = "strong")
        bridge_depth[, core_size := bridge_strong$csize[
          bridge_strong$membership[as.character(seq_len(n_cells))]]]
        bridge_depth[, eligible := eligible & core_size == max(core_size),
          by = component]
      }
      bridge_candidates <- bridge_depth$eligible
    }
    cc <- LineaLinkR::ConnectKnnComponents(nn_idx, nn_dists, geom_emb,
      fixed_main = variant_id %in% c("L1D0contStarLocalRank",
        "L1D0contStarEarlyLocalRank", "L1D0contStarEarlyPotencyRank",
        "L1D0contStarCorePotencyRank"),
      bridge_candidates = bridge_candidates)
    nn_idx <- cc$nn_idx
    nn_dists <- cc$nn_dists
    native_component_count <- as.integer(cc$n_components_before)
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
  cat("  Step 3c/7: Root selection ... ")
  raw_root_labels <- as.character(meta[[fine_col]])
  root_mask <- raw_root_labels %in% c("root", "Root", "STEM", "stem", "early", "Early")
  if (length(llr_cfg$user_root_cells)) {
    stopifnot(all(llr_cfg$user_root_cells %in% cells))
    root_mask <- cells %in% llr_cfg$user_root_cells
  }
  if (length(llr_cfg$user_root_cells) || (isTRUE(llr_cfg$prefer_root_label) && any(root_mask))) {
    root_idx <- which(root_mask)
    root_center <- colMeans(blind_emb[root_idx, , drop = FALSE])
    d2_root <- rowSums((blind_emb[root_idx, , drop = FALSE] -
      matrix(root_center, nrow = length(root_idx), ncol = ncol(blind_emb), byrow = TRUE))^2)
    start_cell <- root_idx[which.min(d2_root)]
    root_anchor_used <- TRUE
    root_evidence_mode <- ifelse(nzchar(llr_cfg$user_root_evidence_source %||% ""),
      "observed_root_prior", "independent_user_root")
    root_input_source <- if (length(llr_cfg$user_root_cells)) {
      paste0(ifelse(nzchar(llr_cfg$user_root_evidence_source %||% ""),
        paste0(llr_cfg$user_root_evidence_source, ":"), ""), llr_cfg$user_root_cells_file)
    } else {
      sprintf("user_supplied_metadata_root:%s", fine_col)
    }
    direction_status <- "externally_anchored"
    start_method <- sprintf("supplied_root_group:center:n=%d", length(root_idx))
    cat(sprintf("%s -> cell %d\n", start_method, start_cell))
  } else if (identical(root_strategy, "root_ensemble")) {
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
      switch(as.character(variant_id %in% c(
        "L1D0contDepthResidualInterval",
        "L1D0contDepthResidualLeafInterval",
        "L1D0contDepthResidualNoDirectInterval",
        "L1D0contDepthResidualOrderedLeafInterval")),
        "TRUE" = {
          high_potency <- as.numeric(stats::lm.fit(
            cbind(1, log1p(as.numeric(Matrix::colSums(expr)))),
            high_potency
          )$residuals)
        })
      sample_col <- FindBiologicalSampleCol(meta, ds_name)
      sample_rank_active <- variant_id %in% c(
        "L1D0contSampleRankRoot", "L1D0contSampleRankLowRoot") &&
        !is.na(sample_col)
      if (variant_id %in% c("L1D0contSampleRankRoot",
                            "L1D0contSampleRankLowRoot",
                            "L1D0contSampleDominanceRoot",
                            "L1D0contTargeted",
                            "L1D0contTargetedInterval",
                            "L1D0contTargetedDonorAdaptiveInterval",
                            "L1D0contDonorStructuralInterval",
                            "L1D0contDualAttractorInterval",
                            "L1D0contDepthResidualInterval",
                            "L1D0contDepthResidualLeafInterval",
                            "L1D0contDepthResidualNoDirectInterval",
                            "L1D0contDepthResidualOrderedLeafInterval") &&
          !is.na(sample_col)) {
        sample_potency <- data.table(
          sample = as.character(meta[[sample_col]]),
          potency = high_potency
        )
        sample_potency[, potency_rank :=
          data.table::frank(potency, ties.method = "average") / .N,
          by = sample]
      }
      if (sample_rank_active)
        high_potency <- sample_potency$potency_rank
      low_potency <- -high_potency
      high_diagnostics <- do.call(LineaLinkR::InferStartCellAttractor, c(list(
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
      ), root_candidate_options))
      if (variant_id %in% c("L1D0contSampleDominanceRoot",
                            "L1D0contTargeted",
                            "L1D0contTargetedInterval",
                            "L1D0contTargetedDonorAdaptiveInterval",
                            "L1D0contDonorStructuralInterval",
                            "L1D0contDualAttractorInterval",
                            "L1D0contDepthResidualInterval",
                            "L1D0contDepthResidualLeafInterval",
                            "L1D0contDepthResidualNoDirectInterval",
                            "L1D0contDepthResidualOrderedLeafInterval") &&
          !is.na(sample_col)) {
        high_root_samples <- as.character(meta[[sample_col]])[
          c(high_diagnostics$start_cell,
            high_diagnostics$cellwise_start_cell)]
        sample_depth <- sample_potency[,
          .(median_potency = median(potency)), by = sample]
        sample_rank_active <- !anyNA(high_root_samples) &&
          uniqueN(high_root_samples) == 1L &&
          nrow(sample_depth) >= 3L &&
          sample_depth[sample == high_root_samples[1L], median_potency] >
          stats::median(sample_depth$median_potency) +
            stats::IQR(sample_depth$median_potency)
        if (sample_rank_active) {
          high_potency <- sample_potency$potency_rank
          low_potency <- -high_potency
          high_diagnostics <- do.call(
            LineaLinkR::InferStartCellAttractor, c(list(
              blind_emb, lineage_cl, nn_idx, nn_dists,
              potency = high_potency,
              top_fraction = llr_cfg$root_top_fraction %||% 0.05,
              n_attractors = llr_cfg$root_n_attractors %||% 0L,
              max_iter = llr_cfg$root_max_iter %||% 300L,
              tol = llr_cfg$root_tol %||% 1e-5,
              weights = llr_cfg$root_weights %||%
                c(0.25, 0.20, 0.30, 0.15, 0.10),
              seed = CONFIG$seed
            ), root_candidate_options))
        }
      }
      low_diagnostics <- do.call(LineaLinkR::InferStartCellAttractor, c(list(
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
      ), root_candidate_options))
      switch(as.character(main_candidate && all(is.na(c(
        high_diagnostics$start_cell, high_diagnostics$cellwise_start_cell,
        low_diagnostics$start_cell, low_diagnostics$cellwise_start_cell)))), "TRUE" = {
        coverage_options <- c(root_candidate_options,
          list(complete_terminal_coverage = TRUE))
        high_diagnostics <- do.call(LineaLinkR::InferStartCellAttractor, c(list(
          blind_emb, lineage_cl, nn_idx, nn_dists, potency = high_potency,
          top_fraction = llr_cfg$root_top_fraction %||% 0.05,
          n_attractors = llr_cfg$root_n_attractors %||% 0L,
          max_iter = llr_cfg$root_max_iter %||% 300L,
          tol = llr_cfg$root_tol %||% 1e-5,
          weights = llr_cfg$root_weights %||% c(0.25, 0.20, 0.30, 0.15, 0.10),
          seed = CONFIG$seed), coverage_options))
        low_diagnostics <- do.call(LineaLinkR::InferStartCellAttractor, c(list(
          blind_emb, lineage_cl, nn_idx, nn_dists, potency = low_potency,
          top_fraction = llr_cfg$root_top_fraction %||% 0.05,
          n_attractors = llr_cfg$root_n_attractors %||% 0L,
          max_iter = llr_cfg$root_max_iter %||% 300L,
          tol = llr_cfg$root_tol %||% 1e-5,
          weights = llr_cfg$root_weights %||% c(0.25, 0.20, 0.30, 0.15, 0.10),
          seed = CONFIG$seed), coverage_options))
        cat(sprintf("[terminal coverage: high +%d, low +%d] ",
          length(high_diagnostics$terminal_coverage_added),
          length(low_diagnostics$terminal_coverage_added)))
      })
      candidate_anchor_cells <- c(
        high_diagnostics$start_cell,
        high_diagnostics$cellwise_start_cell,
        low_diagnostics$start_cell,
        low_diagnostics$cellwise_start_cell
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
        available = !is.na(candidate_anchor_cells),
        status = rep(c(high_diagnostics$status %||% "available",
                       low_diagnostics$status %||% "available"), each = 2L),
        n_unreachable = rep(c(length(high_diagnostics$unreachable_cells),
                              length(low_diagnostics$unreachable_cells)), each = 2L),
        component_score = c(
          as.numeric(high_component_score),
          as.numeric(high_diagnostics$cellwise_component_score),
          as.numeric(low_component_score),
          as.numeric(low_diagnostics$cellwise_component_score)
        ),
        tie_rank = 1:4,
        terminal_anchor_is_high_end = NA,
        terminal_anchor_is_leaf = NA,
        terminal_leaf_eligible = FALSE,
        terminal_descendant_leaves = NA_integer_,
        terminal_segment_eligible = FALSE,
        terminal_root_cluster = NA_integer_,
        terminal_anchor_cluster = NA_integer_,
        terminal_leaf_clusters = NA_character_,
        selector_topology_k = as.integer(selector_topology_k)
      )
      valid_candidate_index <- root_selection[available == TRUE, tie_rank]
      if (!length(valid_candidate_index)) {
        stop("No available root candidate: neither terminal set covers the directed graph")
      }
      candidate_pt_raw <- LineaLinkR:::ComputePseudotimeFuzzyMulti(
        nn_idx, nn_dists, candidate_anchor_cells[valid_candidate_index],
        rank_normalize = rank_norm
      )
      candidate_time_table <- data.table(
        tie_rank = rep(valid_candidate_index, each = n_cells),
        pseudotime = as.vector(candidate_pt_raw)
      )
      candidate_time_table[, pseudotime := LineaLinkR::OrientPseudotimeByPotency(
        pseudotime, high_potency,
        llr_cfg$root_direction_min_potency_cor %||% 0.10
      ), by = tie_rank]
      candidate_valid_pt <- matrix(candidate_time_table$pseudotime, nrow = n_cells)
      candidate_pt <- matrix(NA_real_, nrow = n_cells, ncol = nrow(root_selection))
      candidate_pt[, valid_candidate_index] <- candidate_valid_pt
      candidate_metrics <- candidate_time_table[, .(
        blind_score = LineaLinkR::BlindEmbeddingScore(pseudotime, blind_emb),
        ordering_score = LineaLinkR::BlindOrderingScore(pseudotime, blind_emb),
        potency_spearman = stats::cor(pseudotime, high_potency, method = "spearman"),
        finite_fraction = mean(is.finite(pseudotime)),
        root_cell = which.min(pseudotime)
      ), by = tie_rank]
      root_selection <- merge(root_selection, candidate_metrics,
        by = "tie_rank", all.x = TRUE, sort = FALSE)
      setorder(root_selection, tie_rank)
      candidate_root_cells <- candidate_metrics$root_cell
      candidate_root_clusters <- as.integer(
        lineage_topology$labels[candidate_root_cells]
      )
      candidate_anchor_clusters <- as.integer(
        lineage_topology$labels[candidate_anchor_cells[valid_candidate_index]]
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
      candidate_anchor_pseudotime <- candidate_valid_pt[cbind(
        candidate_anchor_cells[valid_candidate_index],
        seq_len(ncol(candidate_valid_pt))
      )]
      candidate_max_pseudotime <- vapply(
        seq_len(ncol(candidate_valid_pt)),
        function(index) max(candidate_valid_pt[, index], na.rm = TRUE),
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
      candidate_paths <- unlist(candidate_lineages, recursive = FALSE,
        use.names = FALSE)
      candidate_path_nodes <- data.table(
        candidate_index = rep(rep(seq_along(candidate_lineages),
          lengths(candidate_lineages)), lengths(candidate_paths)),
        path = rep(seq_along(candidate_paths), lengths(candidate_paths)),
        cluster = as.integer(unlist(candidate_paths, use.names = FALSE))
      )
      candidate_path_nodes[, terminal_leaf := cluster[.N], by = path]
      candidate_descendant_leaves <- tabulate(unique(candidate_path_nodes[
        cluster == candidate_anchor_clusters[candidate_index],
        .(candidate_index, terminal_leaf)])$candidate_index,
        nbins = length(candidate_lineages))
      root_selection[valid_candidate_index, `:=`(
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
        terminal_descendant_leaves = candidate_descendant_leaves,
        terminal_segment_eligible = candidate_anchor_is_high_end &
          candidate_descendant_leaves == 1L,
        selector_topology_k = as.integer(selector_topology_k)
      )]
      endpoint_candidates <- copy(root_selection)
      endpoint_candidates[, terminal_leaf_eligible := terminal_leaf_eligible |
        (isTRUE(llr_cfg$root_endpoint_terminal_segment) &
          terminal_segment_eligible)]
      endpoint_selection <-
        LineaLinkR::SelectRootEndpointByPotencyAlignment(
        endpoint_candidates,
        default_mode = "detected_genes_high",
        min_alignment_gain =
          llr_cfg$root_endpoint_alignment_min_gain %||% 0.20,
        min_ordering_gain = llr_cfg$root_ordering_min_gain %||% 0.10,
        min_representation_alignment_gain =
          llr_cfg$root_representation_alignment_min_gain %||% 0.00
      )
      direct_low_root <- FALSE
      if (variant$time_method %in% c("bridge_terminal_direct_root",
                                     "raw_multibridge_terminal_direct",
                                     "tangent_multibridge_terminal_direct")) {
        high_state <- root_selection[candidate == "detected_genes_high_state"]
        high_cellwise <- root_selection[
          candidate == "detected_genes_high_cellwise"]
        low_state <- root_selection[candidate == "detected_genes_low_state"]
        hs_lineage_idx <- match(high_state$tie_rank, valid_candidate_index)
        hc_lineage_idx <- match(high_cellwise$tie_rank, valid_candidate_index)
        same_type_coherent <- variant_id %in% c(
          "L1D0contTypeCoherentAnchor",
          "L1D0contLowFlowLeaf",
          "L1D0contSampleDominanceRoot",
          "L1D0contTargeted",
          "L1D0contTargetedInterval",
          "L1D0contTargetedDonorAdaptiveInterval",
          "L1D0contDonorStructuralInterval",
          "L1D0contDualAttractorInterval",
          "L1D0contDepthResidualInterval",
          "L1D0contDepthResidualLeafInterval",
          "L1D0contDepthResidualNoDirectInterval",
          "L1D0contDepthResidualOrderedLeafInterval") &&
          ct_col != "auto_clusters" && n_types > 2L &&
          isTRUE(high_state$available) &&
          isTRUE(high_cellwise$available) &&
          length(candidate_lineages[[hs_lineage_idx]]) ==
          length(candidate_lineages[[hc_lineage_idx]]) &&
          ct_labels[high_state$anchor_cell] != "unassigned" &&
          identical(ct_labels[high_state$anchor_cell],
                    ct_labels[high_cellwise$anchor_cell])
        same_type_anchor <- same_type_coherent &&
          high_state$terminal_root_cluster %in%
            candidate_leaf_clusters[[hc_lineage_idx]]
        quiescent_leaf_hub <- same_type_coherent &&
          !same_type_anchor &&
          !isTRUE(high_state$terminal_leaf_eligible) &&
          isTRUE(high_state$terminal_descendant_leaves > 1L) &&
          isTRUE(high_cellwise$terminal_leaf_eligible) &&
          high_cellwise$ordering_score > high_state$ordering_score &&
          high_cellwise$component_score > high_state$component_score && {
            umi_counts <- as.numeric(Matrix::colSums(expr))
            genes_per_umi <- detected_genes / pmax(umi_counts, 1)
            hub_cl <- high_state$terminal_anchor_cluster
            leaf_cls <- candidate_leaf_clusters[[hs_lineage_idx]]
            stats::median(genes_per_umi[lineage_topology$labels == hub_cl]) <
              stats::median(genes_per_umi[lineage_topology$labels %in% leaf_cls])
          }
        switch(as.character((same_type_anchor || quiescent_leaf_hub) &&
          !nzchar(Sys.getenv('LLR_DISABLE_SAME_TYPE_ANCHOR'))), "TRUE" = {
          endpoint_selection$selected_candidate <-
            "detected_genes_high_cellwise"
          endpoint_selection$selected_potency_mode <- "detected_genes_high"
          endpoint_selection$selected_root_representation <-
            "cellwise_component"
          endpoint_selection$selected_terminal_leaf_eligible <-
            high_cellwise$terminal_leaf_eligible
          endpoint_selection$representation_selection_status <- ifelse(
            quiescent_leaf_hub,
            "same_celltype_quiescent_leaf_hub",
            "same_celltype_terminal_refinement")
          endpoint_selection$representation_ordering_gain <-
            high_cellwise$ordering_score - high_state$ordering_score
        })
        direct_low_root <-
          !variant_id %in% c("L1D0contDepthResidualNoDirectInterval",
                            "L1D0contDepthResidualOrderedLeafInterval") &&
          isTRUE(low_state$available) &&
          !isTRUE(low_state$terminal_leaf_eligible) &&
          high_state$blind_score < 0 && low_state$blind_score > 0 &&
          low_state$blind_score - high_state$blind_score >= 0.20
        if (direct_low_root) {
          endpoint_selection$selected_candidate <- "detected_genes_low_state"
          endpoint_selection$selected_potency_mode <- "detected_genes_low"
          endpoint_selection$selected_root_representation <- "state_aggregate"
          endpoint_selection$selected_terminal_leaf_eligible <- FALSE
          endpoint_selection$endpoint_selection_status <- "blind_direct_low_root"
          endpoint_selection$endpoint_topology_gate_status <- "nonterminal_direct_root"
          endpoint_selection$representation_selection_status <- "state_retained"
          endpoint_selection$direct_root_blind_gain <-
            low_state$blind_score - high_state$blind_score
        }
        axis_switch <- variant_id %in% c("L1D0contAxisCoherent",
          "L1D0contRootBridgeLocalRank",
          "L1D0contTypeCoherentAnchor",
          "L1D0contLowFlowLeaf",
          "L1D0contBranchSource",
          "L1D0contSampleRankRoot",
          "L1D0contSampleRankLowRoot",
          "L1D0contSampleDominanceRoot",
          "L1D0contTargeted",
          "L1D0contTargetedInterval",
          "L1D0contTargetedDonorAdaptiveInterval",
          "L1D0contDonorStructuralInterval",
          "L1D0contDualAttractorInterval",
          "L1D0contDepthResidualInterval",
          "L1D0contDepthResidualLeafInterval",
          "L1D0contDepthResidualNoDirectInterval",
          "L1D0contDepthResidualOrderedLeafInterval") &&
          isTRUE(high_state$available) && isTRUE(low_state$available) &&
          isTRUE(low_state$terminal_leaf_eligible) && {
            high_axis <- abs(stats::cor(candidate_pt[, high_state$tie_rank],
              blind_emb[, 1L], method = "spearman"))
            low_axis <- abs(stats::cor(candidate_pt[, low_state$tie_rank],
              blind_emb[, 1L], method = "spearman"))
            low_axis >= llr_cfg$root_axis_min_alignment &&
              low_axis - high_axis >= llr_cfg$root_axis_min_gain &&
              high_diagnostics$scores$root_score[low_state$root_cell] >=
                high_diagnostics$scores$root_score[high_state$root_cell]
          }
        switch(as.character(axis_switch), "TRUE" = {
          endpoint_selection$selected_candidate <- "detected_genes_low_state"
          endpoint_selection$selected_potency_mode <- "detected_genes_low"
          endpoint_selection$selected_root_representation <- "state_aggregate"
          endpoint_selection$selected_terminal_leaf_eligible <- TRUE
          endpoint_selection$endpoint_selection_status <- "axis_coherent_root"
          endpoint_selection$representation_selection_status <- "state_retained"
          endpoint_selection$axis_alignment_gain <- low_axis - high_axis
          endpoint_selection$axis_root_score_gain <-
            high_diagnostics$scores$root_score[low_state$root_cell] -
            high_diagnostics$scores$root_score[high_state$root_cell]
        })
        if (sample_rank_active &&
            variant_id == "L1D0contSampleRankLowRoot" &&
            isTRUE(low_state$available)) {
          endpoint_selection$selected_candidate <- "detected_genes_low_state"
          endpoint_selection$selected_potency_mode <- "detected_genes_low"
          endpoint_selection$selected_root_representation <- "state_aggregate"
          endpoint_selection$selected_terminal_leaf_eligible <-
            low_state$terminal_leaf_eligible
          endpoint_selection$endpoint_selection_status <-
            "within_sample_low_root"
          endpoint_selection$representation_selection_status <- "state_retained"
          direct_low_root <- FALSE
        }
        if (sample_rank_active &&
            variant_id != "L1D0contSampleRankLowRoot" &&
            isTRUE(high_state$available)) {
          endpoint_selection$selected_candidate <- "detected_genes_high_state"
          endpoint_selection$selected_potency_mode <- "detected_genes_high"
          endpoint_selection$selected_root_representation <- "state_aggregate"
          endpoint_selection$selected_terminal_leaf_eligible <-
            high_state$terminal_leaf_eligible
          endpoint_selection$endpoint_selection_status <-
            "within_sample_high_root"
          endpoint_selection$representation_selection_status <- "state_retained"
          direct_low_root <- FALSE
        }
        if (enable_ordering_low_anchor &&
            endpoint_selection$selected_candidate ==
              "detected_genes_low_state" &&
            low_state$ordering_score - high_state$ordering_score >=
              llr_cfg$root_ordering_min_gain) {
          direct_low_root <- TRUE
          endpoint_selection$endpoint_selection_status <-
            "ordering_supported_direct_low_root"
          endpoint_selection$selected_terminal_leaf_eligible <- FALSE
          endpoint_selection$direct_root_ordering_gain <-
            low_state$ordering_score - high_state$ordering_score
        }
      }
      selected_root_potency_mode <- endpoint_selection$selected_candidate
      selected_row <- root_selection[
        candidate == selected_root_potency_mode
      ][1L]
      terminal_segment_root <- isTRUE(llr_cfg$root_endpoint_terminal_segment) &&
        isTRUE(selected_row$terminal_segment_eligible) &&
        !isTRUE(selected_row$terminal_leaf_eligible) &&
        identical(endpoint_selection$endpoint_selection_status,
          "strong_alignment_switch")
      switch(as.character(terminal_segment_root), "TRUE" = {
        endpoint_selection$endpoint_selection_status <-
          "terminal_segment_alignment_switch"
        endpoint_selection$endpoint_topology_gate_status <-
          "terminal_segment_switch_allowed"
      })
      selected_index <- match(
        selected_root_potency_mode,
        root_selection$candidate
      )
      precomputed_pseudotime <- candidate_pt[, selected_index]
      endpoint_anchor_cell <- as.integer(selected_row$anchor_cell)
      terminal_anchor_cluster <- as.integer(
        selected_row$terminal_anchor_cluster
      )
      terminal_leaf_eligible <- ifelse(direct_low_root, FALSE,
        as.logical(selected_row$terminal_leaf_eligible) | terminal_segment_root)
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
      switch(as.character(terminal_segment_root), "TRUE" = {
        start_cell <- as.integer(selected_row$root_cell)
      })
      quiescent_leaf_hub_active <- exists("quiescent_leaf_hub", inherits = FALSE) &&
        isTRUE(quiescent_leaf_hub) &&
        identical(endpoint_selection$representation_selection_status,
                  "same_celltype_quiescent_leaf_hub")
      switch(as.character(quiescent_leaf_hub_active), "TRUE" = {
        hc_paths <- candidate_lineages[[hc_lineage_idx]]
        hs_leaves <- candidate_leaf_clusters[[hs_lineage_idx]]
        stub_paths <- Filter(
          function(p) length(p) == 2L && p[2L] %in% hs_leaves,
          hc_paths)
        start_cell <- switch(as.character(length(stub_paths) >= 1L),
          "TRUE" = {
            leaf_cl <- stub_paths[[1L]][2L]
            ix <- which(lineage_topology$labels == leaf_cl)
            ctr <- colMeans(blind_emb[ix, , drop = FALSE])
            as.integer(ix[which.min(rowSums((blind_emb[ix, , drop = FALSE] -
              matrix(ctr, length(ix), ncol(blind_emb), byrow = TRUE))^2))])
          },
          "FALSE" = as.integer(selected_row$root_cell))
        pt_hub <- as.numeric(candidate_pt_raw[, hs_lineage_idx])
        pt_rev_hub <- (max(pt_hub, na.rm = TRUE) - pt_hub) /
          diff(range(pt_hub, na.rm = TRUE))
        leaf_medoids <- vapply(hs_leaves, function(cid) {
          ix <- which(lineage_topology$labels == cid)
          ctr <- colMeans(blind_emb[ix, , drop = FALSE])
          as.integer(ix[which.min(rowSums((blind_emb[ix, , drop = FALSE] -
            matrix(ctr, length(ix), ncol(blind_emb), byrow = TRUE))^2))])
        }, integer(1L))
        d_leaves <- do.call(pmin, as.data.frame(
          LineaLinkR:::ComputePseudotimeFuzzyMulti(
            nn_idx, nn_dists, leaf_medoids, rank_normalize = rank_norm)))
        pt_te <- d_leaves / pmax(d_leaves + pt_hub, .Machine$double.eps)
        pt_te <- (pt_te - min(pt_te, na.rm = TRUE)) /
          diff(range(pt_te, na.rm = TRUE))
        precomputed_pseudotime <- 0.5 * pt_rev_hub + 0.5 * pt_te
        precomputed_pseudotime <- (precomputed_pseudotime -
          min(precomputed_pseudotime, na.rm = TRUE)) /
          diff(range(precomputed_pseudotime, na.rm = TRUE))
        endpoint_anchor_cell <- as.integer(high_state$anchor_cell)
      })
      if (!quiescent_leaf_hub_active &&
          (terminal_segment_root ||
          !(variant$time_method %in% c("terminal_distance",
                                       "bridge_terminal_distance",
                                       "bridge_terminal_direct_root",
                                       "raw_multibridge_terminal_direct",
                                       "tangent_multibridge_terminal_direct") &&
            identical(endpoint_selection$selected_potency_mode,
                      "detected_genes_low")))) {
        precomputed_pseudotime <- LineaLinkR::ComputePseudotimeFuzzy(
          nn_idx,
          nn_dists,
          start_cell,
          rank_normalize = rank_norm
        )
      }
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
        ifelse(terminal_segment_root, "terminal_segment_selected_point_root",
          ifelse(direct_low_root, "attractor_low_state_direct_root",
                 "attractor_terminal_distance_oriented")),
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
      root_diagnostics$terminal_segment_root <- terminal_segment_root
      root_diagnostics$terminal_descendant_leaves <-
        as.integer(selected_row$terminal_descendant_leaves)
      root_diagnostics$terminal_leaf_clusters <- terminal_leaf_clusters
      root_diagnostics$selector_topology_k <- as.integer(selector_topology_k)
      root_diagnostics$potency_mode <- root_potency_mode
      root_diagnostics$selected_potency_mode <- selected_root_potency_mode
      root_diagnostics$root_selection <- root_selection
      root_diagnostics$candidate_unreachable_cells <- list(
        detected_genes_high = high_diagnostics$unreachable_cells,
        detected_genes_low = low_diagnostics$unreachable_cells
      )
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
      root_ordering_confidence <- max(root_selection$ordering_score, na.rm = TRUE)
      use_confidence_regularized_fallback <- root_ordering_confidence <
        (llr_cfg$root_ordering_fallback_threshold %||% 0.45) &&
        !direct_low_root &&
        !(variant$time_method %in% c("endpoint_anchor", "terminal_distance",
                                     "bridge_terminal_distance",
                                     "bridge_terminal_direct_root",
                                     "raw_multibridge_terminal_direct",
                                     "tangent_multibridge_terminal_direct") &&
          identical(endpoint_selection$selected_potency_mode,
                    "detected_genes_low") &&
          isTRUE(endpoint_selection$selected_terminal_leaf_eligible) &&
          endpoint_selection$endpoint_alignment_gain >=
            llr_cfg$root_endpoint_alignment_min_gain)
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
        fallback_attractor <- do.call(LineaLinkR::InferStartCellAttractor, c(list(
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
        ), root_candidate_options))
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

  switch(as.character(
    variant_id %in% c("L1D0contRetentionSourceLeaf",
                      "L1D0contRetentionSourceLeafInterval",
                      "L1D0contDepthResidualLeafInterval",
                      "L1D0contDepthResidualNoDirectInterval",
                      "L1D0contDepthResidualOrderedLeafInterval") &&
      isTRUE(llr_cfg$connect_components) &&
      cc$n_components_before == 1L &&
      selected_root_potency_mode %in% c("detected_genes_high_state",
                                       "detected_genes_high_cellwise")
  ), "TRUE" = {
    source_leaf_clusters <- lineage_topology$labels
    source_leaf_paths <- LineaLinkR::ExtractLineages(
      geom_emb, source_leaf_clusters,
      as.integer(source_leaf_clusters[start_cell]))
    source_leaf_tips <- vapply(source_leaf_paths, tail, integer(1L), 1L)
    source_leaf_scores <- data.table(
      cluster = source_leaf_clusters,
      potency = root_diagnostics$scores$potency,
      source_flow = root_diagnostics$scores$source_flow
    )[, .(potency = median(potency),
          source_flow = median(source_flow)), by = cluster]
    source_leaf_root <- source_leaf_scores[
      cluster == source_leaf_clusters[start_cell]]
    source_leaf_tip <- source_leaf_scores[
      cluster %in% source_leaf_tips][order(-potency, -source_flow)][1L]
    switch(as.character(
      source_leaf_tip$cluster != source_leaf_root$cluster &&
        source_leaf_tip$potency > source_leaf_root$potency &&
        source_leaf_tip$source_flow > source_leaf_root$source_flow
    ), "TRUE" = {
      source_leaf_cells <- which(
        source_leaf_clusters == source_leaf_tip$cluster)
      source_leaf_center <- colMeans(
        pca_embedding[source_leaf_cells, , drop = FALSE])
      source_leaf_distance <- rowSums((
        pca_embedding[source_leaf_cells, , drop = FALSE] -
          matrix(source_leaf_center, nrow = length(source_leaf_cells),
                 ncol = ncol(pca_embedding), byrow = TRUE)
      )^2)
      source_leaf_cell <- source_leaf_cells[which.min(source_leaf_distance)]
      source_leaf_ordering_gain <- switch(as.character(
        variant_id == "L1D0contDepthResidualOrderedLeafInterval"),
        "TRUE" = {
          source_leaf_pt <- LineaLinkR::ComputePseudotimeFuzzy(
            nn_idx, nn_dists, source_leaf_cell, rank_normalize = rank_norm)
          LineaLinkR::BlindOrderingScore(source_leaf_pt, blind_emb) -
            LineaLinkR::BlindOrderingScore(precomputed_pseudotime, blind_emb)
        },
        "FALSE" = NA_real_)
      switch(as.character(
        variant_id != "L1D0contDepthResidualOrderedLeafInterval" ||
          source_leaf_ordering_gain > 0), "TRUE" = {
        original_start_cell <- start_cell
        start_cell <- source_leaf_cell
        precomputed_pseudotime <- switch(as.character(
          variant_id == "L1D0contDepthResidualOrderedLeafInterval"),
          "TRUE" = source_leaf_pt, "FALSE" = NULL)
        root_diagnostics$source_leaf_root <- list(
          previous_root = original_start_cell,
          leaf_cluster = source_leaf_tip$cluster,
          potency_gain = source_leaf_tip$potency - source_leaf_root$potency,
          source_flow_gain = source_leaf_tip$source_flow -
            source_leaf_root$source_flow,
          ordering_gain = source_leaf_ordering_gain
        )
        start_method <- paste0(start_method, ":source_dominant_leaf")
        cat(sprintf("  [source leaf root: %d -> %d] ",
          original_start_cell, start_cell))
      })
    })
  })
  root_override <- NULL
  if (variant_id %in% c("L1D0contLowFlowLeaf", "L1D0contTargeted",
                        "L1D0contTargetedInterval",
                        "L1D0contTargetedDonorAdaptiveInterval",
                        "L1D0contDonorStructuralInterval",
                        "L1D0contDualAttractorInterval",
                        "L1D0contDepthResidualInterval",
                        "L1D0contDepthResidualLeafInterval",
                        "L1D0contDepthResidualNoDirectInterval",
                        "L1D0contDepthResidualOrderedLeafInterval") &&
      !(variant_id %in% c("L1D0contTargetedDonorAdaptiveInterval",
                          "L1D0contDonorStructuralInterval") &&
        geometry_excludes_lr) &&
      selected_root_potency_mode == "detected_genes_low_state" &&
      identical(endpoint_selection$endpoint_selection_status,
                "strong_alignment_switch") &&
      isTRUE(llr_cfg$connect_components) &&
      cc$n_components_before == 1L) {
    flow_clusters <- lineage_topology$labels
    flow_paths <- LineaLinkR::ExtractLineages(
      geom_emb, flow_clusters, as.integer(flow_clusters[start_cell]))
    flow_tips <- vapply(flow_paths, tail, integer(1L), 1L)
    flow_scores <- data.table(cluster = flow_clusters,
      source_flow = root_diagnostics$scores$source_flow)[,
        .(source_flow = median(source_flow)), by = cluster]
    flow_leaf <- flow_scores[cluster %in% flow_tips][order(-source_flow)][1L]
    if (flow_leaf$cluster != flow_clusters[start_cell]) {
      flow_cells <- which(flow_clusters == flow_leaf$cluster)
      flow_center <- colMeans(pca_embedding[flow_cells, , drop = FALSE])
      flow_distance <- rowSums((pca_embedding[flow_cells, , drop = FALSE] -
        matrix(flow_center, nrow = length(flow_cells),
          ncol = ncol(pca_embedding), byrow = TRUE))^2)
      flow_cell <- flow_cells[which.min(flow_distance)]
      flow_pt <- LineaLinkR::ComputePseudotimeFuzzy(
        nn_idx, nn_dists, flow_cell, rank_normalize = rank_norm
      )
      flow_ordering_gain <- LineaLinkR::BlindOrderingScore(flow_pt, blind_emb) -
        LineaLinkR::BlindOrderingScore(precomputed_pseudotime, blind_emb)
      flow_time_consistent <- TRUE
      flow_time_columns <- intersect(c("time_days", "age(days)",
        "time_hours", "timepoint"), names(meta))
      switch(as.character(enable_timepoint_source_flow_guard &&
        length(flow_time_columns) > 0L), "TRUE" = {
        flow_observed_time <- as.numeric(gsub("[^0-9.]", "",
          as.character(meta[[flow_time_columns[[1L]]]])))
        flow_first_time <- min(flow_observed_time, na.rm = TRUE)
        switch(as.character(length(unique(flow_observed_time[
          is.finite(flow_observed_time)])) > 1L &&
          flow_observed_time[start_cell] == flow_first_time &&
          flow_observed_time[flow_cell] == flow_first_time), "TRUE" = {
          flow_time_consistent <- stats::cor(flow_pt,
            flow_observed_time, method = "spearman", use = "complete.obs") >=
            stats::cor(precomputed_pseudotime, flow_observed_time,
              method = "spearman", use = "complete.obs")
        })
      })
      switch(as.character(!enable_preserve_selected_root &&
        flow_ordering_gain > 0 &&
        flow_time_consistent), "TRUE" = {
        old_start_cell <- start_cell
        start_cell <- flow_cell
        precomputed_pseudotime <- flow_pt
        root_diagnostics$flow_leaf_root <- list(
          previous_root = old_start_cell,
          leaf_cluster = flow_leaf$cluster,
          source_flow = flow_leaf$source_flow,
          ordering_score_gain = flow_ordering_gain
        )
        start_method <- "lowstate_max_source_flow_leaf_root"
        cat(sprintf("  [source-flow leaf root %d -> %d] ",
          old_start_cell, start_cell))
      })
      switch(as.character(variant_id %in% c(
        "L1D0contTargetedInterval",
        "L1D0contTargetedDonorAdaptiveInterval",
        "L1D0contDonorStructuralInterval",
        "L1D0contDualAttractorInterval",
        "L1D0contDepthResidualInterval",
        "L1D0contDepthResidualLeafInterval",
        "L1D0contDepthResidualNoDirectInterval",
        "L1D0contDepthResidualOrderedLeafInterval") &&
        !nzchar(Sys.getenv('LLR_DISABLE_OCH_ROOT')) &&
        flow_ordering_gain <= 0 &&
        root_selection[candidate == "detected_genes_high_state",
                       ordering_score] >
        root_selection[candidate == "detected_genes_low_state",
                       ordering_score]), "TRUE" = {
        high_row <- root_selection[candidate == "detected_genes_high_state"]
        root_override <- list(
          start_cell = as.integer(high_row$root_cell),
          pseudotime = candidate_pt[, high_row$tie_rank],
          diagnostics = high_diagnostics,
          method = "ordering_consistent_high_root",
          ordering_gain = high_row$ordering_score -
            root_selection[candidate == "detected_genes_low_state",
                           ordering_score],
          source_leaf = FALSE
        )
      })
    }
  }

  switch(as.character(!enable_preserve_selected_root &&
    variant_id %in% c(
    "L1D0contTargetedInterval",
    "L1D0contTargetedDonorAdaptiveInterval",
    "L1D0contDonorStructuralInterval",
    "L1D0contDualAttractorInterval",
    "L1D0contDepthResidualInterval",
    "L1D0contDepthResidualLeafInterval",
    "L1D0contDepthResidualNoDirectInterval",
    "L1D0contDepthResidualOrderedLeafInterval") &&
    is.null(root_override) &&
    selected_root_potency_mode == "detected_genes_low_state" &&
    LineaLinkR::BlindOrderingScore(precomputed_pseudotime, blind_emb) <
      (llr_cfg$root_ordering_fallback_threshold %||% 0.45)), "TRUE" = {
    retention_diagnostics <- LineaLinkR::InferStartCellAttractor(
      blind_emb, lineage_cl, nn_idx, nn_dists,
      potency = high_potency,
      top_fraction = llr_cfg$root_top_fraction %||% 0.05,
      n_attractors = llr_cfg$root_n_attractors %||% 0L,
      max_iter = llr_cfg$root_max_iter %||% 300L,
      tol = llr_cfg$root_tol %||% 1e-5,
      weights = llr_cfg$root_weights %||% c(0.25, 0.20, 0.30, 0.15, 0.10),
      seed = CONFIG$seed,
      on_unreachable = "return",
      attractor_selection = "retention"
    )
    switch(as.character(retention_diagnostics$status == "available"),
      "TRUE" = {
      retention_cell <- retention_diagnostics$start_cell
      retention_pt <- LineaLinkR::ComputePseudotimeFuzzy(
        nn_idx, nn_dists, retention_cell, rank_normalize = rank_norm)
      retention_ordering <- LineaLinkR::BlindOrderingScore(
        retention_pt, blind_emb)
      retention_leaf <- FALSE
      retention_clusters <- lineage_topology$labels
      retention_paths <- LineaLinkR::ExtractLineages(
        geom_emb, retention_clusters,
        as.integer(retention_clusters[retention_cell]))
      retention_tips <- vapply(retention_paths, tail, integer(1L), 1L)
      retention_scores <- data.table(
        cluster = retention_clusters,
        potency = retention_diagnostics$scores$potency,
        source_flow = retention_diagnostics$scores$source_flow
      )[, .(potency = median(potency),
            source_flow = median(source_flow)), by = cluster]
      retention_root <- retention_scores[
        cluster == retention_clusters[retention_cell]]
      retention_tip <- retention_scores[
        cluster %in% retention_tips][order(-potency, -source_flow)][1L]
      switch(as.character(
        retention_tip$cluster != retention_root$cluster &&
          retention_tip$potency > retention_root$potency &&
          retention_tip$source_flow > retention_root$source_flow
      ), "TRUE" = {
        leaf_cells <- which(retention_clusters == retention_tip$cluster)
        leaf_center <- colMeans(pca_embedding[leaf_cells, , drop = FALSE])
        leaf_distance <- rowSums((pca_embedding[leaf_cells, , drop = FALSE] -
          matrix(leaf_center, nrow = length(leaf_cells),
                 ncol = ncol(pca_embedding), byrow = TRUE))^2)
        leaf_cell <- leaf_cells[which.min(leaf_distance)]
        leaf_pt <- LineaLinkR::ComputePseudotimeFuzzy(
          nn_idx, nn_dists, leaf_cell, rank_normalize = rank_norm)
        leaf_ordering <- LineaLinkR::BlindOrderingScore(leaf_pt, blind_emb)
        switch(as.character(leaf_ordering > retention_ordering), "TRUE" = {
          retention_cell <- leaf_cell
          retention_pt <- leaf_pt
          retention_ordering <- leaf_ordering
          retention_leaf <- TRUE
        })
      })
      current_ordering <- LineaLinkR::BlindOrderingScore(
        precomputed_pseudotime, blind_emb)
      retention_time_consistent <- TRUE
      retention_time_columns <- intersect(c("time_days", "age(days)",
        "time_hours", "timepoint"), names(meta))
      switch(as.character(enable_timepoint_source_flow_guard &&
        length(retention_time_columns) > 0L), "TRUE" = {
        retention_observed_time <- as.numeric(gsub("[^0-9.]", "",
          as.character(meta[[retention_time_columns[[1L]]]])))
        retention_first_time <- min(retention_observed_time, na.rm = TRUE)
        switch(as.character(length(unique(retention_observed_time[
          is.finite(retention_observed_time)])) > 1L &&
          retention_observed_time[start_cell] == retention_first_time &&
          retention_observed_time[retention_cell] == retention_first_time),
          "TRUE" = {
          retention_time_consistent <- stats::cor(retention_pt,
            retention_observed_time, method = "spearman",
            use = "complete.obs") >=
            stats::cor(precomputed_pseudotime, retention_observed_time,
              method = "spearman", use = "complete.obs")
        })
      })
      switch(as.character(retention_ordering > current_ordering &&
        retention_time_consistent), "TRUE" = {
        root_override <- list(
          start_cell = as.integer(retention_cell),
          pseudotime = retention_pt,
          diagnostics = retention_diagnostics,
          method = ifelse(retention_leaf,
            "retention_source_leaf_ordering_rescue",
            "retention_ordering_rescue"),
          ordering_gain = retention_ordering - current_ordering,
          source_leaf = retention_leaf
        )
      })
    })
  })
  switch(as.character(variant_id == "L1D0contDualAttractorInterval" &&
    is.null(root_override) &&
    selected_root_potency_mode == "detected_genes_high_state"), "TRUE" = {
    retention_diagnostics <- LineaLinkR::InferStartCellAttractor(
      blind_emb, lineage_cl, nn_idx, nn_dists,
      potency = high_potency,
      top_fraction = llr_cfg$root_top_fraction %||% 0.05,
      n_attractors = llr_cfg$root_n_attractors %||% 0L,
      max_iter = llr_cfg$root_max_iter %||% 300L,
      tol = llr_cfg$root_tol %||% 1e-5,
      weights = llr_cfg$root_weights %||% c(0.25, 0.20, 0.30, 0.15, 0.10),
      seed = CONFIG$seed,
      on_unreachable = "return",
      attractor_selection = "retention"
    )
    switch(as.character(retention_diagnostics$status == "available" &&
      retention_diagnostics$start_cell != start_cell), "TRUE" = {
      retention_pt <- LineaLinkR::ComputePseudotimeFuzzy(
        nn_idx, nn_dists, retention_diagnostics$start_cell,
        rank_normalize = rank_norm)
      retention_gain <- LineaLinkR::BlindOrderingScore(
        retention_pt, blind_emb) -
        LineaLinkR::BlindOrderingScore(precomputed_pseudotime, blind_emb)
      switch(as.character(retention_gain > 0), "TRUE" = {
        root_override <- list(
          start_cell = as.integer(retention_diagnostics$start_cell),
          pseudotime = retention_pt,
          diagnostics = retention_diagnostics,
          method = "retention_ordering_root",
          ordering_gain = retention_gain,
          source_leaf = FALSE
        )
      })
    })
  })
  switch(as.character(!is.null(root_override)), "TRUE" = {
    previous_root <- start_cell
    start_cell <- root_override$start_cell
    precomputed_pseudotime <- root_override$pseudotime
    selected_root_potency_mode <- "detected_genes_high_state"
    endpoint_anchor_cell <- start_cell
    endpoint_selection$selected_candidate <- selected_root_potency_mode
    endpoint_selection$selected_potency_mode <- "detected_genes_high"
    endpoint_selection$selected_root_representation <- "state_aggregate"
    endpoint_selection$endpoint_selection_status <- root_override$method
    endpoint_selection$endpoint_topology_gate_status <-
      root_override$method
    endpoint_selection$representation_selection_status <-
      root_override$method
    endpoint_selection$representation_ordering_gain <-
      root_override$ordering_gain
    endpoint_selection$selected_terminal_leaf_eligible <-
      root_override$source_leaf
    root_diagnostics <- root_override$diagnostics
    root_diagnostics$potency_mode <- root_potency_mode
    root_diagnostics$selected_potency_mode <- selected_root_potency_mode
    root_diagnostics$root_selection <- root_selection
    root_diagnostics$endpoint_selection <- endpoint_selection
    root_diagnostics$start_cell <- start_cell
    root_diagnostics$terminal_anchor_cell <- start_cell
    root_diagnostics$root_region <- order(
      precomputed_pseudotime, -high_potency, seq_len(n_cells))[
        seq_len(max(2L, ceiling(n_cells *
          (llr_cfg$root_top_fraction %||% 0.05))))]
    root_diagnostics$root_cluster <- as.integer(
      lineage_topology$labels[start_cell])
    root_diagnostics$root_replacement <- list(
      previous_root = previous_root,
      ordering_gain = root_override$ordering_gain,
      source_leaf = root_override$source_leaf)
    root_potency_vector <- high_potency
    terminal_anchor_cluster <- root_diagnostics$root_cluster
    terminal_leaf_eligible <- root_override$source_leaf
    terminal_leaf_clusters <- as.integer(vapply(
      LineaLinkR::ExtractLineages(
        geom_emb, lineage_topology$labels,
        as.integer(terminal_anchor_cluster)),
      tail, integer(1L), 1L))
    root_diagnostics$terminal_anchor_cluster <- terminal_anchor_cluster
    root_diagnostics$terminal_leaf_eligible <- terminal_leaf_eligible
    root_diagnostics$terminal_leaf_clusters <- terminal_leaf_clusters
    root_input_source <- paste0(
      "expression_detected_genes_pca_knn_no_gt:", root_override$method)
    start_method <- root_override$method
    cat(sprintf("  [%s %d -> %d] ",
      start_method, previous_root, start_cell))
  })

  switch(as.character(enable_timepoint_branchpoint_root &&
    (identical(root_evidence_mode, "anchored_bidirectional_potency_ordering") ||
      (enable_timepoint_dominant_early_leaf &&
        identical(root_evidence_mode,
          "ambiguous_confidence_regularized_root_ensemble"))) &&
    !isTRUE(llr_cfg$augment_mst_bridges) &&
    length(intersect(c("time_days", "age(days)", "time_hours", "time_info"), names(meta))) > 0L &&
    !is.null(root_diagnostics$scores$potency)), "TRUE" = {
    source_time_col <- intersect(c("time_days", "age(days)", "time_hours", "time_info"),
      names(meta))[[1L]]
    source_time <- as.numeric(meta[[source_time_col]])
    source_first_time <- min(source_time, na.rm = TRUE)
    source_early_fraction <- mean(source_time == source_first_time,
      na.rm = TRUE)
    source_clusters <- lineage_topology$labels
    source_tree <- LineaLinkR::BuildKnnClusterTree(
      knn_conn, source_clusters, as.integer(source_clusters[start_cell]))
    source_degree <- tabulate(
      c(source_tree$edges$from, source_tree$edges$to),
      nbins = max(source_clusters))
    source_states <- data.table(
      cluster = source_clusters,
      time = source_time,
      potency = root_diagnostics$scores$potency,
      fate_entropy = root_diagnostics$scores$fate_entropy,
      attractor_distance = root_diagnostics$scores$nearest_attractor_distance
    )[, .(earliest_cells = sum(time == source_first_time, na.rm = TRUE),
          earliest_fraction = mean(time == source_first_time, na.rm = TRUE),
          potency = median(potency),
          fate_entropy = median(fate_entropy),
          attractor_distance = median(attractor_distance)), by = cluster]
    source_hubs <- source_states[source_degree[cluster] >= 3L &
      earliest_fraction > source_early_fraction][
      order(potency, cluster)]
    switch(as.character(nrow(source_hubs) > 0L &&
      identical(root_evidence_mode,
        "anchored_bidirectional_potency_ordering")), "TRUE" = {
      current_cluster <- as.integer(source_clusters[start_cell])
      source_neighbors <- c(
        source_tree$edges$to[source_tree$edges$from == current_cluster],
        source_tree$edges$from[source_tree$edges$to == current_cluster])
      source_majority <- source_hubs[
        cluster %in% source_neighbors &
          earliest_cells > 0.5 * sum(source_time == source_first_time,
                                    na.rm = TRUE)]
      source_hub <- switch(as.character(source_degree[current_cluster] == 1L &&
        nrow(source_majority) > 0L),
        "TRUE" = source_majority[order(-earliest_cells, potency, cluster)][1L],
        "FALSE" = source_hubs[1L])
      source_cluster <- source_hub$cluster[[1L]]
      source_graph <- igraph::graph_from_data_frame(
        source_tree$edges[, c("from", "to")], directed = FALSE,
        vertices = data.frame(name = as.character(seq_len(max(source_clusters)))))
      source_hops <- as.numeric(igraph::distances(
        source_graph, v = as.character(current_cluster),
        to = as.character(source_cluster)))
      fate_gain <- source_hub$fate_entropy[[1L]] -
        source_states[cluster == current_cluster, fate_entropy]
      distance_gain <- source_hub$attractor_distance[[1L]] -
        source_states[cluster == current_cluster, attractor_distance]
      source_cells <- which(source_clusters == source_cluster &
        source_time == source_first_time)
      source_center <- colMeans(blind_emb[source_cells, , drop = FALSE])
      source_distance <- rowSums((blind_emb[source_cells, , drop = FALSE] -
        matrix(source_center, nrow = length(source_cells),
          ncol = ncol(blind_emb), byrow = TRUE))^2)
      source_cell <- source_cells[which.min(source_distance)]
      source_pt <- LineaLinkR::ComputePseudotimeFuzzy(
        nn_idx, nn_dists, source_cell, rank_normalize = rank_norm)
      current_pt <- LineaLinkR::ComputePseudotimeFuzzy(
        nn_idx, nn_dists, start_cell, rank_normalize = rank_norm)
      source_ordering <- LineaLinkR::BlindOrderingScore(source_pt, blind_emb)
      current_ordering <- LineaLinkR::BlindOrderingScore(current_pt, blind_emb)
      switch(as.character((source_hops > 2L ||
        (source_hops >= 1L &&
          source_hub$earliest_fraction[[1L]] > 0.5 &&
          source_states[cluster == current_cluster, earliest_fraction] <=
            source_early_fraction) ||
        (source_degree[current_cluster] == 1L && source_hops == 1L &&
          source_hub$earliest_cells[[1L]] >
            0.5 * sum(source_time == source_first_time, na.rm = TRUE))) &&
        source_cell != start_cell &&
        source_time[source_cell] == source_first_time &&
        (source_ordering > current_ordering ||
          (fate_gain > stats::IQR(source_states$fate_entropy) &&
            distance_gain > stats::IQR(source_states$attractor_distance)))),
        "TRUE" = {
        previous_root <- start_cell
        start_cell <- source_cell
        precomputed_pseudotime <- source_pt
        endpoint_anchor_cell <- as.integer(source_cell)
        terminal_anchor_cluster <- as.integer(source_cluster)
        terminal_leaf_eligible <- FALSE
        terminal_leaf_clusters <- integer()
        root_diagnostics$start_cell <- as.integer(source_cell)
        root_diagnostics$root_cluster <- as.integer(source_cluster)
        root_diagnostics$root_region <- as.integer(source_cells)
        root_diagnostics$branchpoint_source <- list(
          previous_root = previous_root,
          root_cluster = as.integer(source_cluster),
          tree_hops = source_hops,
          time_column = source_time_col,
          source_time = source_time[source_cell],
          earliest_fraction = source_hub$earliest_fraction[[1L]],
          ordering_gain = source_ordering - current_ordering,
          fate_entropy_gain = fate_gain,
          attractor_distance_gain = distance_gain)
        root_evidence_mode <- "observed_timepoint_branchpoint_source"
        root_input_source <- paste0("expression_pca_knn_plus_", source_time_col)
        start_method <- paste0(start_method, ":timepoint_branchpoint_source")
        cat(sprintf("  [timepoint branchpoint root: %d -> %d] ",
          previous_root, start_cell))
      })
    })
    switch(as.character(nrow(source_hubs) == 0L &&
      source_time[start_cell] > source_first_time), "TRUE" = {
      source_leaves <- source_states[source_degree[cluster] == 1L &
        earliest_fraction > source_early_fraction &
        earliest_fraction > 0.5][
          order(-earliest_cells, potency, cluster)]
      switch(as.character(nrow(source_leaves) > 0L), "TRUE" = {
        source_cluster <- source_leaves$cluster[[1L]]
        source_neighbor <- c(
          source_tree$edges$to[source_tree$edges$from == source_cluster],
          source_tree$edges$from[source_tree$edges$to == source_cluster])
        switch(as.character(source_states[
          cluster == source_neighbor, earliest_fraction] >
          source_early_fraction ||
          (enable_timepoint_dominant_early_leaf &&
            source_leaves$earliest_cells[[1L]] >
              sum(source_time == source_first_time, na.rm = TRUE) / 2)),
          "TRUE" = {
          source_cells <- which(source_clusters == source_cluster &
            source_time == source_first_time)
          source_center <- colMeans(blind_emb[source_cells, , drop = FALSE])
          source_distance <- rowSums((blind_emb[source_cells, , drop = FALSE] -
            matrix(source_center, nrow = length(source_cells),
              ncol = ncol(blind_emb), byrow = TRUE))^2)
          source_cell <- source_cells[which.min(source_distance)]
          source_graph <- igraph::graph_from_data_frame(
            source_tree$edges[, c("from", "to")], directed = FALSE,
            vertices = data.frame(name = as.character(
              seq_len(max(source_clusters)))))
          source_hops <- as.numeric(igraph::distances(
            source_graph, v = as.character(source_clusters[start_cell]),
            to = as.character(source_cluster)))
          switch(as.character(source_hops > 2L), "TRUE" = {
            source_pt <- LineaLinkR::ComputePseudotimeFuzzy(
              nn_idx, nn_dists, source_cell, rank_normalize = rank_norm)
            previous_root <- start_cell
            start_cell <- source_cell
            precomputed_pseudotime <- source_pt
            endpoint_anchor_cell <- as.integer(source_cell)
            terminal_anchor_cluster <- as.integer(source_cluster)
            terminal_leaf_eligible <- FALSE
            terminal_leaf_clusters <- integer()
            root_diagnostics$start_cell <- as.integer(source_cell)
            root_diagnostics$root_cluster <- as.integer(source_cluster)
            root_diagnostics$root_region <- as.integer(source_cells)
            root_diagnostics$early_leaf_source <- list(
              previous_root = previous_root,
              root_cluster = as.integer(source_cluster),
              next_cluster = as.integer(source_neighbor),
              tree_hops = source_hops,
              time_column = source_time_col,
              source_time = source_time[source_cell],
              earliest_fraction = source_leaves$earliest_fraction[[1L]])
            root_evidence_mode <- "observed_timepoint_early_leaf_source"
            root_input_source <- paste0("expression_pca_knn_plus_",
              source_time_col)
            start_method <- paste0(start_method,
              ":timepoint_early_leaf_source")
            cat(sprintf("  [timepoint early leaf root: %d -> %d] ",
              previous_root, start_cell))
          })
        })
      })
    })
  })

  start_cluster <- as.integer(lineage_cl[start_cell])
  if (enable_sparse_overlap_pc1) {
    detected <- expr > 0
    median_shared <- median(Matrix::colSums(
      detected * detected[, nn_idx[, 1L], drop = FALSE]
    ))
    if (median_shared < ncol(pca_embedding)) {
      valid_edge <- nn_idx > 0L
      pc1_gap <- abs(pca_embedding[row(nn_idx)[valid_edge], 1L] -
                       pca_embedding[nn_idx[valid_edge], 1L])
      nn_dists[valid_edge] <- nn_dists[valid_edge] *
        (1 + pc1_gap / median(pc1_gap))
      precomputed_pseudotime <- NULL
      cat(sprintf("  [shared genes %.0f < %d PCs: PC1 edge cost] ",
                  median_shared, ncol(pca_embedding)))
    }
  }
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
      precomputed_pseudotime <- LineaLinkR::ComputePseudotimeFuzzy(
        nn_idx,
        nn_dists,
        start_cell,
        rank_normalize = rank_norm
      )
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

  trunk_orientation_flip <- FALSE
  trunk_time_origin_cell <- NA_integer_
  switch(as.character(enable_connected_long_trunk_orientation &&
    isTRUE(native_component_count == 1L) &&
    identical(root_evidence_mode, "anchored_bidirectional_potency_ordering")),
    "TRUE" = {
      tree_clusters <- as.integer(lineage_topology$labels)
      selected_cluster <- as.integer(tree_clusters[start_cell])
      tree_paths <- LineaLinkR::ExtractLineages(
        geom_emb, tree_clusters, selected_cluster)
      path_nodes <- data.table(
        node = unlist(tree_paths, use.names = FALSE),
        lineage = rep(seq_along(tree_paths), lengths(tree_paths)))
      path_nodes[, child := shift(node, type = "lead"), by = lineage]
      tree_edges <- unique(path_nodes[!is.na(child),
        .(from = node, to = child)])
      tree_centers <- rowsum(geom_emb, tree_clusters) /
        as.vector(table(tree_clusters))
      tree_edges[, weight := sqrt(rowSums((
        tree_centers[as.character(from), , drop = FALSE] -
          tree_centers[as.character(to), , drop = FALSE])^2))]
      cluster_tree <- igraph::graph_from_data_frame(
        tree_edges, directed = FALSE)
      tree_degree <- igraph::degree(cluster_tree)
      tree_leaves <- names(tree_degree[tree_degree == 1L])
      tree_junctions <- names(tree_degree[tree_degree >= 3L])
      switch(as.character(length(tree_leaves) >= 2L &&
        length(tree_junctions) >= 1L &&
        as.character(selected_cluster) %in% tree_leaves), "TRUE" = {
          trunk_distance <- igraph::distances(
            cluster_tree, v = tree_leaves, to = tree_junctions,
            weights = igraph::edge_attr(cluster_tree, "weight"))
          trunk <- data.table(
            leaf = as.integer(tree_leaves),
            length = do.call(pmin, as.data.frame(trunk_distance))
          )[order(-length)]
          switch(as.character(trunk$leaf[1L] != selected_cluster &&
            trunk$length[1L] >=
              CONFIG$llr$long_trunk_orientation_min_ratio * trunk$length[2L]),
            "TRUE" = {
              source_cells <- which(tree_clusters == trunk$leaf[1L])
              source_center <- colMeans(
                geom_emb[source_cells, , drop = FALSE])
              trunk_time_origin_cell <- as.integer(source_cells[which.min(rowSums((
                geom_emb[source_cells, , drop = FALSE] -
                  matrix(source_center, nrow = length(source_cells),
                    ncol = ncol(geom_emb), byrow = TRUE))^2))])
              trunk_orientation_flip <- TRUE
              root_diagnostics$orientation_override <- list(
                distance_anchor_cell = start_cell,
                time_origin_cell = trunk_time_origin_cell,
                source_leaf = trunk$leaf[1L],
                trunk_ratio = trunk$length[1L] / trunk$length[2L])
              root_evidence_mode <- "connected_long_trunk_orientation"
              root_input_source <- paste0(
                root_input_source, ":cluster_tree_geometry")
              start_method <- paste0(start_method,
                ":connected_long_trunk_orientation")
              cat(sprintf("  [long trunk orientation %d -> %d] ",
                start_cell, trunk_time_origin_cell))
            })
        })
    })

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
  if (variant$time_method %in% c("raw_multibridge_terminal_direct",
                                 "tangent_multibridge_terminal_direct") &&
      max(bridge_component_labels) > 2L) {
    costs <- as.data.table(LineaLinkR:::ComputePseudotimeCostEdges(
      nn_idx, nn_dists
    ))
    costs[, metric_cost := raw_distance]
    if (identical(variant$time_method, "tangent_multibridge_terminal_direct")) {
      costs[, bridge := bridge_component_labels[from] !=
        bridge_component_labels[to]]
      stopifnot(sum(costs$bridge) == max(bridge_component_labels) - 1L)
      component_tangent <- do.call(rbind, lapply(
        split(seq_len(n_cells), bridge_component_labels),
        function(index) {
          local_embedding <- scale(geom_emb[index, , drop = FALSE],
                                   center = TRUE, scale = FALSE)
          eigen(crossprod(local_embedding), symmetric = TRUE)$vectors[, 1L]
        }
      ))
      bridge_rows <- which(costs$bridge)
      bridge_vector <- geom_emb[costs$to[bridge_rows], , drop = FALSE] -
        geom_emb[costs$from[bridge_rows], , drop = FALSE]
      source_tangent <- component_tangent[
        bridge_component_labels[costs$from[bridge_rows]], , drop = FALSE]
      target_tangent <- component_tangent[
        bridge_component_labels[costs$to[bridge_rows]], , drop = FALSE]
      costs$metric_cost[bridge_rows] <- (
        abs(rowSums(bridge_vector * source_tangent)) +
          abs(rowSums(bridge_vector * target_tangent))
      ) / 2
    }
    graph <- igraph::graph_from_data_frame(
      costs[, .(from, to, metric_cost)],
      directed = FALSE,
      vertices = data.frame(name = as.character(seq_len(n_cells)))
    )
    metric_distance <- as.numeric(igraph::distances(
      graph, v = as.character(start_cell), to = igraph::V(graph),
      weights = igraph::E(graph)$metric_cost,
      algorithm = "dijkstra"
    ))
    pseudotime <- metric_distance / max(metric_distance)
  } else if (variant$time_method %in% c("bridge_floor_cost",
                                 "bridge_terminal_distance",
                                 "bridge_terminal_direct_root",
                                 "raw_multibridge_terminal_direct",
                                 "tangent_multibridge_terminal_direct") &&
      max(bridge_component_labels) > 1L) {
    costs <- as.data.table(LineaLinkR:::ComputePseudotimeCostEdges(
      nn_idx, nn_dists
    ))
    costs[, bridge := bridge_component_labels[from] !=
      bridge_component_labels[to]]
    stopifnot(sum(costs$bridge) == max(bridge_component_labels) - 1L)
    costs[, floor_bridge := bridge &
      round(fuzzy_cost / raw_distance) == 1e6]
    if (any(costs$floor_bridge)) {
      unit_scale <- median(costs$fuzzy_cost[!costs$bridge]) /
        median(costs$raw_distance[!costs$bridge])
      costs[, calibrated_cost := fifelse(bridge,
        raw_distance * unit_scale, fuzzy_cost)]
      graph <- igraph::graph_from_data_frame(
        costs[, .(from, to, calibrated_cost)],
        directed = FALSE,
        vertices = data.frame(name = as.character(seq_len(n_cells)))
      )
      bridge_distance <- as.numeric(igraph::distances(
        graph, v = as.character(start_cell), to = igraph::V(graph),
        weights = igraph::E(graph)$calibrated_cost,
        algorithm = "dijkstra"
      ))
      pseudotime <- bridge_distance / max(bridge_distance)
    }
  }
  if (variant_id == "L1D0contSparsePearsonDiffusion") {
    diffusion_edges <- as.data.table(LineaLinkR:::ComputePseudotimeCostEdges(
      nn_idx, nn_dists
    ))
    diffusion_edges <- rbindlist(list(
      diffusion_edges[, .(i = from, j = to, distance = raw_distance)],
      diffusion_edges[, .(i = to, j = from, distance = raw_distance)]
    ))
    bandwidth <- diffusion_edges[, .(value = max(distance)), by = i]
    local_bandwidth <- bandwidth$value[
      match(seq_len(n_cells), bandwidth$i)
    ]
    diffusion_edges[, weight := exp(-distance^2 /
      (local_bandwidth[i] * local_bandwidth[j]))]
    diffusion_graph <- Matrix::sparseMatrix(
      i = diffusion_edges$i, j = diffusion_edges$j,
      x = diffusion_edges$weight, dims = c(n_cells, n_cells)
    )
    diffusion_degree <- as.numeric(Matrix::rowSums(diffusion_graph))
    degree_scale <- Matrix::Diagonal(x = 1 / sqrt(diffusion_degree))
    diffusion_operator <- degree_scale %*% diffusion_graph %*% degree_scale
    spectral <- RSpectra::eigs_sym(
      diffusion_operator,
      k = CONFIG$llr$sparse_diffusion_dims + 1L,
      which = "LA",
      opts = list(initvec = rep(1, n_cells))
    )
    eig_order <- order(spectral$values, decreasing = TRUE)
    eigenvalues <- spectral$values[eig_order]
    eigenvectors <- spectral$vectors[, eig_order, drop = FALSE]
    diffusion_embedding <- sweep(
      eigenvectors[, -1L, drop = FALSE], 1L,
      sqrt(diffusion_degree), "/"
    )
    diffusion_embedding <- sweep(
      diffusion_embedding, 2L,
      eigenvalues[-1L] / (1 - eigenvalues[-1L]), "*"
    )
    root_embedding <- matrix(
      diffusion_embedding[start_cell, ],
      nrow = n_cells, ncol = ncol(diffusion_embedding), byrow = TRUE
    )
    diffusion_distance <- sqrt(rowSums(
      (diffusion_embedding - root_embedding)^2
    ))
    pseudotime <- diffusion_distance / max(diffusion_distance)
  }
  if (!isTRUE(root_anchor_used)) {
    pseudotime <- LineaLinkR::OrientPseudotimeBlind(pseudotime, blind_emb)
  }
  if (variant_id %in% c("L1D0contTangentLocalRank",
                        "L1D0contStarLocalRank",
                        "L1D0contStarEarlyLocalRank",
                        "L1D0contStarEarlyPotencyRank",
                        "L1D0contStarCorePotencyRank",
                        "L1D0contRootBridgeLocalRank",
                        "L1D0contTypeCoherentAnchor",
                        "L1D0contLowFlowLeaf",
                        "L1D0contBranchSource",
                        "L1D0contSampleRankRoot",
                        "L1D0contSampleRankLowRoot",
                        "L1D0contSampleDominanceRoot",
                        "L1D0contTargeted",
                        "L1D0contTargetedInterval",
                        "L1D0contTargetedDonorAdaptiveInterval",
                        "L1D0contDonorStructuralInterval",
                        "L1D0contDualAttractorInterval",
                        "L1D0contDepthResidualInterval",
                        "L1D0contDepthResidualLeafInterval",
                        "L1D0contDepthResidualNoDirectInterval",
                        "L1D0contDepthResidualOrderedLeafInterval") &&
      max(bridge_component_labels) > 1L) {
    global_axis <- as.numeric(pca_embedding[, 1L]) *
      sign(stats::cor(pca_embedding[, 1L], pseudotime,
        method = "spearman"))
    component_time <- data.table(
      cell_index = seq_len(n_cells),
      component = bridge_component_labels,
      current = as.numeric(pseudotime),
      global_axis = global_axis
    )
    component_time[, local_axis := stats::prcomp(
      pca_embedding[cell_index, , drop = FALSE],
      center = TRUE, scale. = FALSE)$x[, 1L], by = component]
    if (variant_id %in% c("L1D0contStarEarlyPotencyRank",
                        "L1D0contStarCorePotencyRank") &&
        selected_root_potency_mode %in% c("detected_genes_high_state",
                                           "detected_genes_high_cellwise")) {
      component_time[, local_sign := -sign(stats::cor(
        local_axis, root_potency_vector[cell_index], method = "spearman")),
        by = component]
    } else {
      component_time[, local_sign := sign(stats::cor(
        local_axis, global_axis, method = "spearman")), by = component]
    }
    component_time[, local_time := sort(current)[data.table::frank(
      local_axis * local_sign, ties.method = "first")],
      by = component]
    component_time[component == bridge_component_labels[start_cell],
      local_time := current]
    switch(as.character(enable_component_axis_repair), "TRUE" = {
      component_time[, local_global_alignment := stats::cor(
        local_axis, global_axis, method = "spearman"), by = component]
      component_time[, local_graph_alignment := stats::cor(
        local_axis, current, method = "spearman"), by = component]
      component_time[
        component != bridge_component_labels[start_cell] &
          abs(local_global_alignment) <
            CONFIG$llr$component_axis_alignment_min,
        local_time := sort(current)[data.table::frank(
          global_axis, ties.method = "first")], by = component]
      component_time[
        component == bridge_component_labels[start_cell] &
          abs(local_global_alignment) > llr_cfg$root_axis_min_alignment &
          abs(local_graph_alignment) <
            CONFIG$llr$component_axis_alignment_min,
        local_time := sort(current)[data.table::frank(
          local_axis * sign(local_global_alignment),
          ties.method = "first")], by = component]
    })
    pseudotime <- component_time$local_time
  }
  switch(as.character(trunk_orientation_flip), "TRUE" = {
    pseudotime <- 1 - pseudotime
    start_cell <- trunk_time_origin_cell
    root_diagnostics$start_cell <- start_cell
  })
  bridge_clock_gap <- max(diff(sort(pseudotime)))
  bridge_clock_rank_applied <- main_candidate &&
    max(bridge_component_labels) > 1L &&
    bridge_clock_gap > CONFIG$llr$bridge_clock_max_gap
  switch(as.character(bridge_clock_rank_applied), "TRUE" = {
    pseudotime <- (rank(pseudotime, ties.method = "average") - 1) /
      (n_cells - 1L)
  })
  shared_induction_applied <- FALSE
  shared_induction_n_genes <- 0L
  switch(as.character(enable_shared_induction_clock &&
    max(bridge_component_labels) > 1L), "TRUE" = {
      component_membership <- Matrix::sparseMatrix(
        i = seq_len(n_cells), j = bridge_component_labels, x = 1,
        dims = c(n_cells, max(bridge_component_labels)))
      component_sizes <- as.numeric(Matrix::colSums(component_membership))
      centered_time <- pseudotime - (
        as.numeric(Matrix::crossprod(component_membership, pseudotime)) /
          component_sizes)[bridge_component_labels]
      clock_expression <- expr %*% Matrix::Diagonal(
        x = 10000 / Matrix::colSums(expr))
      clock_expression@x <- log1p(clock_expression@x)
      clock_covariance <- as.matrix(clock_expression %*% (
        Matrix::Diagonal(x = centered_time) %*% component_membership))
      clock_sum <- as.matrix(clock_expression %*% component_membership)
      clock_square_sum <- as.matrix((clock_expression^2) %*%
        component_membership)
      clock_variance <- pmax(clock_square_sum - sweep(
        clock_sum^2, 2L, component_sizes, "/"), 0)
      time_variance <- as.numeric(Matrix::crossprod(
        component_membership, centered_time^2))
      clock_correlation <- clock_covariance / sqrt(sweep(
        clock_variance, 2L, time_variance, "*"))
      clock_correlation[!is.finite(clock_correlation)] <- 0
      induction_genes <- rowSums(clock_correlation >
        CONFIG$llr$shared_induction_component_r) >= ceiling(
          CONFIG$llr$shared_induction_component_fraction *
            ncol(clock_correlation))
      shared_induction_n_genes <- sum(induction_genes)
      switch(as.character(shared_induction_n_genes >=
        CONFIG$llr$shared_induction_min_genes), "TRUE" = {
          clock_z <- t(scale(t(as.matrix(clock_expression[
            induction_genes, , drop = FALSE]))))
          clock_z[!is.finite(clock_z)] <- 0
          clock_score <- colMeans(clock_z)
          clock_score <- clock_score * sign(stats::cor(
            clock_score, pseudotime, method = "spearman"))
          pseudotime <- (clock_score - min(clock_score)) /
            diff(range(clock_score))
          shared_induction_applied <- TRUE
          cat(sprintf("  [shared induction clock: %d genes] ",
            shared_induction_n_genes))
        })
    })
  centroid_phase_applied <- FALSE
  switch(as.character(enable_centroid_phase &&
    !shared_induction_applied && max(bridge_component_labels) > 1L),
    "TRUE" = {
      component_membership <- Matrix::sparseMatrix(
        i = seq_len(n_cells), j = bridge_component_labels, x = 1,
        dims = c(n_cells, max(bridge_component_labels)))
      component_sizes <- as.numeric(Matrix::colSums(component_membership))
      component_centers <- sweep(as.matrix(Matrix::crossprod(
        component_membership,
        pca_embedding[, seq_len(min(2L, ncol(pca_embedding))),
          drop = FALSE])), 1L, component_sizes, "/")
      root_center <- component_centers[bridge_component_labels[start_cell], ]
      centroid_distance <- sqrt(rowSums((component_centers - matrix(
        root_center, nrow = nrow(component_centers),
        ncol = ncol(component_centers), byrow = TRUE))^2))
      component_stage <- centroid_distance / max(centroid_distance)
      component_mean_time <- as.numeric(Matrix::crossprod(
        component_membership, pseudotime)) / component_sizes
      phase_time <- component_stage[bridge_component_labels] + pseudotime -
        component_mean_time[bridge_component_labels]
      pseudotime <- (phase_time - min(phase_time)) / diff(range(phase_time))
      centroid_phase_applied <- TRUE
      cat(sprintf("  [centroid phase: %d components] ",
        length(component_stage)))
    })
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
  direction_diagnostics$bridge_clock_gap <- bridge_clock_gap
  direction_diagnostics$bridge_clock_rank_applied <- bridge_clock_rank_applied
  direction_diagnostics$shared_induction_applied <- shared_induction_applied
  direction_diagnostics$shared_induction_n_genes <- shared_induction_n_genes
  direction_diagnostics$centroid_phase_applied <- centroid_phase_applied
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
  source_topology_cl <- NULL
  transition_fragment_labels <- integer()
  if ("--root-topology-k" %in% arg_keys) {
    source_topology_cl <- lineage_topology$labels
    lineage_topology <- LineaLinkR::ClusterLineageTopology(
      geom_emb,
      k = min(as.integer(llr_cfg$lineage_topology_k), n_cells),
      nstart = as.integer(llr_cfg$lineage_topology_nstart %||% 10L),
      iter_max = 500L,
      seed = CONFIG$seed
    )
  }
  if (lineage_topology$empty_start_warnings > 0L) {
    cat(sprintf(
      "[audited empty starts=%d] ",
      lineage_topology$empty_start_warnings
    ))
  }
  lineage_cl <- lineage_topology$labels
  switch(as.character(enable_knn_topology && nzchar(lineage_state_column)), "TRUE" = {
    lineage_cl <- as.integer(factor(as.character(meta[[lineage_state_column]])))
  })
  start_cluster <- as.integer(lineage_cl[start_cell])
  lineage_cluster_source <- sprintf(
    "kmeans:%s:algorithm=%s:requested_k=%d:occupied_k=%d:purpose=lineage_topology",
    lineage_embedding_mode,
    lineage_topology$algorithm,
    lineage_topology$requested_k,
    lineage_topology$occupied_k
  )
  knn_topology_diag <- NULL
  if (use_knn_topology &&
      (!(donor_split_geometry_mode || donor_center_lr_mode) ||
       "--knn-topology" %in% arg_keys)) {
    # Variant A: paths, weights and fork placement all come from one tree
    # built on the NATIVE cell kNN graph (knn_conn.mtx; ConnectKnnComponents
    # bridges are in-memory only and never enter this tree).
    knn_tree <- LineaLinkR::BuildKnnClusterTree(knn_conn, lineage_cl, start_cluster)
    knn_paths <- LineaLinkR::ExtractLineagesFromKnnGraph(
      knn_conn, lineage_cl, start_cluster
    )
    lineage_paths <- knn_paths$paths
    lineage_weights <- LineaLinkR::ComputeLineageWeightsFromKnnGraph(
      knn_conn, geom_emb, lineage_cl, start_cluster, tree = knn_tree
    )
    knn_topology_diag <- list(
      n_tree_edges = nrow(knn_tree$edges),
      n_components = knn_tree$n_components,
      root_component = knn_tree$root_component,
      n_paths = length(lineage_paths),
      tree_edges = knn_tree$edges
    )
    cat(sprintf("  [knn-tree: %d edges, %d comps, %d paths] ",
                nrow(knn_tree$edges), knn_tree$n_components,
                length(lineage_paths)))
  } else {
    lineage_weights <- LineaLinkR::ComputeLineageWeightsTrunkAware(geom_emb, lineage_cl, start_cluster)
    if (!is.matrix(lineage_weights) || nrow(lineage_weights) != n_cells || ncol(lineage_weights) < 1L) {
      stop(sprintf("Invalid lineage weights for %s", ds_name))
    }
    n_lineages <- ncol(lineage_weights)
    lineage_paths <- LineaLinkR::ExtractLineages(
      geom_emb, lineage_cl, start_cluster
    )
  }
  if (enable_terminal_sibling_refine) {
    fine_cl <- LineaLinkR::ClusterLineageTopology(
      geom_emb, k = min(25L, n_cells),
      nstart = as.integer(llr_cfg$lineage_topology_nstart %||% 10L),
      iter_max = 500L, seed = CONFIG$seed
    )$labels
    path_nodes <- data.table(
      lineage = rep(seq_along(lineage_paths), lengths(lineage_paths)),
      node = unlist(lineage_paths, use.names = FALSE)
    )
    terminal_parent <- path_nodes[
      , .(coarse = node[.N], parent = node[.N - 1L]), by = lineage
    ][, .(coarse, parent)]
    knn_conn_t <- Matrix::t(knn_conn)
    sym_knn_conn <- (knn_conn + knn_conn_t +
                     abs(knn_conn - knn_conn_t)) / 2
    graph_edges <- as.data.table(summary(sym_knn_conn))
    graph_edges[, `:=`(
      coarse_i = lineage_cl[i], coarse_j = lineage_cl[j],
      fine_i = fine_cl[i], fine_j = fine_cl[j]
    )]
    graph_edges <- graph_edges[
      coarse_i == coarse_j & coarse_i %in% terminal_parent$coarse
    ]
    fragment_size <- data.table(coarse = lineage_cl, fine = fine_cl)[
      , .(n = .N), by = .(coarse, fine)
    ]
    fragment_volume <- graph_edges[
      , .(volume = sum(x)), by = .(coarse = coarse_i, fine = fine_i)
    ]
    fragment_inside <- graph_edges[fine_i == fine_j,
      .(inside = sum(x)), by = .(coarse = coarse_i, fine = fine_i)
    ]
    cluster_volume <- graph_edges[
      , .(cluster_volume = sum(x)), by = .(coarse = coarse_i)
    ]
    candidates <- merge(fragment_size, fragment_volume,
      by = c("coarse", "fine"))
    candidates <- merge(candidates, fragment_inside,
      by = c("coarse", "fine"), all.x = TRUE)
    candidates[is.na(inside), inside := 0]
    candidates <- merge(candidates, cluster_volume, by = "coarse")
    candidates[, conductance := (volume - inside) /
      pmin(volume, cluster_volume - volume)]

    state_counts <- data.table(coarse = lineage_cl,
      fine = fine_cl, cell_type = as.character(ct_labels))
    coarse_type <- state_counts[
      , .(coarse_n_type = .N), by = .(coarse, cell_type)
    ]
    fine_type <- state_counts[
      , .(fine_n_type = .N), by = .(coarse, fine, cell_type)
    ]
    type_grid <- merge(candidates[, .(coarse, fine)], coarse_type,
      by = "coarse", allow.cartesian = TRUE)
    type_grid <- merge(type_grid, fine_type,
      by = c("coarse", "fine", "cell_type"), all.x = TRUE)
    type_grid[is.na(fine_n_type), fine_n_type := 0L]
    split_majority <- type_grid[
      , .(after_majority = max(fine_n_type) +
          max(coarse_n_type - fine_n_type)),
      by = .(coarse, fine)
    ]
    coarse_majority <- coarse_type[
      , .(before_majority = max(coarse_n_type)), by = coarse
    ]
    coarse_size <- state_counts[, .(coarse_size = .N), by = coarse]
    candidates <- merge(candidates, split_majority,
      by = c("coarse", "fine"))
    candidates <- merge(candidates, coarse_majority, by = "coarse")
    candidates <- merge(candidates, coarse_size, by = "coarse")
    candidates[, `:=`(
      purity_gain = (after_majority - before_majority) / coarse_size,
      minor_fraction = pmin(n, coarse_size - n) / coarse_size
    )]
    selected <- candidates[
      n >= 20L & minor_fraction >= 0.05 &
        conductance <= 0.1 & purity_gain >= 0.05
    ][order(coarse, -purity_gain, conductance)][
      , .SD[1L], by = coarse
    ]
    selected <- merge(selected, terminal_parent, by = "coarse")
    if (nrow(selected) > 0L) {
      selected[, new_label := max(lineage_cl) + seq_len(.N)]
      cell_groups <- data.table(coarse = lineage_cl, fine = fine_cl)
      cell_groups[selected, on = .(coarse, fine),
        new_label := i.new_label]
      provisional_cl <- as.integer(ifelse(
        is.na(cell_groups$new_label), lineage_cl, cell_groups$new_label
      ))
      provisional_tree <- LineaLinkR::BuildKnnClusterTree(
        knn_conn, provisional_cl, provisional_cl[start_cell]
      )
      support <- provisional_tree$s_matrix
      new_other <- support[selected$new_label, , drop = FALSE]
      old_other <- support[selected$coarse, , drop = FALSE]
      new_other[cbind(seq_len(nrow(selected)), selected$coarse)] <- 0
      old_other[cbind(seq_len(nrow(selected)), selected$new_label)] <- 0
      selected <- selected[
        max.col(new_other, ties.method = "first") == parent &
        max.col(old_other, ties.method = "first") == parent &
        support[cbind(new_label, parent)] > 0
      ]
      if (nrow(selected) > 0L) {
        selected[, new_label := max(lineage_cl) + seq_len(.N)]
        cell_groups[, new_label := NULL]
        cell_groups[selected, on = .(coarse, fine),
          new_label := i.new_label]
        lineage_cl <- as.integer(ifelse(
          is.na(cell_groups$new_label), lineage_cl, cell_groups$new_label
        ))
        start_cluster <- as.integer(lineage_cl[start_cell])
        knn_tree <- LineaLinkR::BuildKnnClusterTree(
          knn_conn, lineage_cl, start_cluster
        )
        tree_edges <- as.data.table(knn_tree$edges)
        edge_keys <- paste(tree_edges$from, tree_edges$to, sep = ":")
        replacement <- selected[
          paste(pmin(coarse, new_label), pmax(coarse, new_label),
                sep = ":") %in% edge_keys &
          paste(pmin(coarse, parent), pmax(coarse, parent),
                sep = ":") %in% edge_keys
        ]
        sibling_edges <- replacement[
          , .(from = pmin(coarse, new_label),
              to = pmax(coarse, new_label))
        ]
        tree_edges <- tree_edges[!sibling_edges,
          on = .(from, to)]
        parent_edges <- replacement[
          , .(from = pmin(parent, new_label),
              to = pmax(parent, new_label),
              support = knn_tree$s_matrix[cbind(parent, new_label)])
        ]
        knn_tree$edges <- rbind(tree_edges, parent_edges)
        lineage_paths <- lapply(
          LineaLinkR:::rs_extract_lineages_from_edges(
            as.integer(knn_tree$edges$from),
            as.integer(knn_tree$edges$to),
            as.integer(length(unique(lineage_cl))),
            as.integer(start_cluster)
          ), as.integer
        )
        lineage_weights <- LineaLinkR::ComputeLineageWeightsFromKnnGraph(
          knn_conn, geom_emb, lineage_cl, start_cluster,
          tree = knn_tree
        )
        lineage_cluster_source <- paste0(
          lineage_cluster_source, "+terminal_state_siblings"
        )
        knn_topology_diag$n_tree_edges <- nrow(knn_tree$edges)
        knn_topology_diag$n_paths <- length(lineage_paths)
        knn_topology_diag$tree_edges <- knn_tree$edges
        knn_topology_diag$terminal_sibling_splits <- selected[
          , .(coarse, fine, new_label, parent,
              conductance, purity_gain)
        ]
      }
    }
  }
  if (enable_internal_transition_pocket) {
    path_nodes <- data.table(
      lineage = rep(seq_along(lineage_paths), lengths(lineage_paths)),
      position = sequence(lengths(lineage_paths)),
      node = unlist(lineage_paths, use.names = FALSE)
    )
    path_nodes[, parent := data.table::shift(node), by = lineage]
    tree_links <- unique(path_nodes[!is.na(parent), .(parent, node)])
    incoming <- tree_links[,
      .(n_parent = .N, parent = parent[1L]), by = node]
    outgoing <- tree_links[,
      .(n_child = .N, child = node[1L]), by = parent]
    setnames(outgoing, "parent", "node")
    pocket_nodes <- merge(incoming[n_parent == 1L],
      outgoing[n_child == 1L], by = "node")
    node_size <- data.table(node = lineage_cl)[, .(total_n = .N), by = node]
    pocket_nodes <- merge(pocket_nodes, node_size, by = "node")[total_n >= 3L]
    pocket_cells <- data.table(cell = seq_len(n_cells), node = lineage_cl)[
      pocket_nodes, on = "node", nomatch = 0L]
    pocket_fragments <- pocket_cells[, .(
      cell = cell,
      fragment = {
        set.seed(CONFIG$seed)
        stats::kmeans(geom_emb[cell, , drop = FALSE],
          centers = 3L, nstart = 1L, iter.max = 500L)$cluster
      }
    ), by = node]
    node_of_cell <- rep(NA_integer_, n_cells)
    fragment_of_cell <- rep(NA_integer_, n_cells)
    node_of_cell[pocket_fragments$cell] <- pocket_fragments$node
    fragment_of_cell[pocket_fragments$cell] <- pocket_fragments$fragment
    contacts <- as.data.table(summary(knn_conn))
    contacts[, `:=`(
      node = node_of_cell[i],
      fragment = fragment_of_cell[i],
      neighbor_node = lineage_cl[j]
    )]
    contacts <- merge(contacts[!is.na(node)],
      pocket_nodes[, .(node, parent, child)], by = "node")
    parent_contacts <- contacts[neighbor_node == parent,
      .(parent_links = sum(x)), by = .(node, fragment)]
    child_contacts <- contacts[neighbor_node == child,
      .(child_links = sum(x)), by = .(node, fragment)]
    fragment_stats <- pocket_fragments[, .(
      n = .N, pt_median = median(pseudotime[cell])),
      by = .(node, fragment)]
    fragment_stats <- merge(fragment_stats, parent_contacts,
      by = c("node", "fragment"), all.x = TRUE)
    fragment_stats <- merge(fragment_stats, child_contacts,
      by = c("node", "fragment"), all.x = TRUE)
    fragment_stats[is.na(parent_links), parent_links := 0]
    fragment_stats[is.na(child_links), child_links := 0]
    fragment_stats[, `:=`(
      parent_rate = parent_links / n,
      child_rate = child_links / n
    )]
    bridge <- fragment_stats[order(node, -child_rate), .SD[1L], by = node]
    isolated <- fragment_stats[order(node, child_rate), .SD[1L], by = node]
    pocket_geometry <- merge(bridge, isolated, by = "node",
      suffixes = c("_bridge", "_isolated"))
    pocket_geometry <- merge(pocket_geometry,
      fragment_stats[, .(min_fragment_n = min(n), total_n = sum(n)),
        by = node], by = "node")
    pocket_geometry <- merge(pocket_geometry,
      pocket_cells[, .(pt_iqr = IQR(pseudotime[cell])), by = node],
      by = "node")
    pocket_geometry[, `:=`(
      parent_ratio = parent_rate_isolated / parent_rate_bridge,
      child_ratio = child_rate_isolated / child_rate_bridge,
      time_separation_iqr =
        (pt_median_isolated - pt_median_bridge) / pt_iqr,
      min_fragment_fraction = min_fragment_n / total_n
    )]
    selected_pocket <- head(pocket_geometry[
      is.finite(parent_ratio) & is.finite(child_ratio) &
        is.finite(time_separation_iqr) &
        parent_ratio <= internal_transition_pocket_rule$boundary_depletion &
        child_ratio <= internal_transition_pocket_rule$boundary_depletion &
        time_separation_iqr >=
          internal_transition_pocket_rule$min_time_separation_iqr &
        min_fragment_n >= internal_transition_pocket_rule$min_fragment_n &
        min_fragment_fraction >=
          internal_transition_pocket_rule$min_fragment_fraction
    ][order(child_ratio, parent_ratio, -time_separation_iqr)], 1L)
    knn_topology_diag$internal_transition_pocket_candidates <-
      pocket_geometry[, .(node, parent_ratio, child_ratio,
        time_separation_iqr, min_fragment_n, min_fragment_fraction)]
    internal_transition_cluster <- switch(
      as.character(nrow(selected_pocket) > 0L),
      "TRUE" = selected_pocket$node[1L], "FALSE" = NULL)
  }
  if (!is.null(internal_transition_cluster)) {
    transition_cells <- which(lineage_cl == internal_transition_cluster)
    set.seed(CONFIG$seed)
    transition_groups <- stats::kmeans(
      geom_emb[transition_cells, , drop = FALSE],
      centers = 3L, nstart = 1L, iter.max = 500L
    )$cluster
    first_fragment <- max(lineage_cl) + 1L
    lineage_cl[transition_cells] <- first_fragment - 1L + transition_groups
    lineage_cl <- match(lineage_cl, sort(unique(lineage_cl)))
    fragment_labels <- lineage_cl[transition_cells][
      match(seq_len(3L), transition_groups)
    ]
    transition_fragment_labels <- fragment_labels[
      rep(enable_internal_transition_pocket, length(fragment_labels))]
    start_cluster <- as.integer(lineage_cl[start_cell])
    knn_tree <- LineaLinkR::BuildKnnClusterTree(
      knn_conn, lineage_cl, start_cluster
    )
    lineage_paths <- LineaLinkR::ExtractLineagesFromKnnGraph(
      knn_conn, lineage_cl, start_cluster
    )$paths
    lineage_weights <- LineaLinkR::ComputeLineageWeightsFromKnnGraph(
      knn_conn, geom_emb, lineage_cl, start_cluster, tree = knn_tree
    )
    knn_topology_diag$n_tree_edges <- nrow(knn_tree$edges)
    knn_topology_diag$n_paths <- length(lineage_paths)
    knn_topology_diag$tree_edges <- knn_tree$edges
    knn_topology_diag$internal_transition_split <- data.table(
      coarse = internal_transition_cluster,
      fragment = fragment_labels,
      n = tabulate(transition_groups, 3L)
    )
    lineage_cluster_source <- paste0(
      lineage_cluster_source, "+internal_transition_k3"
    )
  }
  if (!is.matrix(lineage_weights) || nrow(lineage_weights) != n_cells || ncol(lineage_weights) < 1L) {
    stop(sprintf("Invalid lineage weights for %s", ds_name))
  }
  n_lineages <- ncol(lineage_weights)
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

  if (variant_id %in% c("L1D0contBranchSource",
                        "L1D0contTargeted",
                        "L1D0contTargetedInterval",
                        "L1D0contTargetedDonorAdaptiveInterval",
                        "L1D0contDonorStructuralInterval",
                        "L1D0contDualAttractorInterval",
                        "L1D0contDepthResidualInterval",
                        "L1D0contDepthResidualLeafInterval",
                        "L1D0contDepthResidualNoDirectInterval",
                        "L1D0contDepthResidualOrderedLeafInterval") && n_lineages > 1L &&
      selected_root_potency_mode %in% c("detected_genes_high_state",
                                   "detected_genes_high_cellwise") &&
      max(bridge_component_labels) == 1L) {
    branch_cl <- lineage_cl
    branch_paths <- lineage_paths
    branch_weights <- lineage_weights
    if (!is.null(source_topology_cl)) {
      branch_cl <- source_topology_cl
      branch_start <- as.integer(branch_cl[start_cell])
      branch_paths <- LineaLinkR::ExtractLineages(
        geom_emb, branch_cl, branch_start
      )
      branch_weights <- LineaLinkR::ComputeLineageWeightsTrunkAware(
        geom_emb, branch_cl, branch_start
      )
    }
    leaf_clusters <- vapply(branch_paths, tail, integer(1L), 1L)
    cluster_size <- tabulate(branch_cl)
    cluster_source_score <- tapply(root_diagnostics$scores$root_score,
      branch_cl, mean)
    cluster_potency <- tapply(root_diagnostics$scores$potency,
      branch_cl, mean)
    source_mass <- cluster_size[leaf_clusters] *
      cluster_source_score[leaf_clusters]
    source_lineage <- which.max(source_mass)
    source_cluster <- leaf_clusters[source_lineage]
    source_ancestors <- branch_paths[[source_lineage]][
      -length(branch_paths[[source_lineage]])]
    valley_cluster <- source_ancestors[which.min(
      cluster_potency[source_ancestors])]
    source_rebound <- cluster_potency[source_cluster] -
      cluster_potency[valley_cluster]
    branch_source_active <-
      cluster_potency[source_cluster] >
        max(cluster_potency[leaf_clusters[-source_lineage]]) &&
      source_rebound > 0 &&
      median(pseudotime[branch_cl == source_cluster]) >
        median(pseudotime[branch_cl == valley_cluster])
    root_diagnostics$branch_source <- list(
      active = branch_source_active,
      lineage = source_lineage,
      source_cluster = source_cluster,
      valley_cluster = valley_cluster,
      potency_rebound = source_rebound,
      source_mass_fraction = source_mass[source_lineage] / sum(source_mass)
    )
    if (branch_source_active && !nzchar(Sys.getenv('LLR_DISABLE_BRANCH_SOURCE'))) {
      source_cells <- which(branch_cl == source_cluster)
      source_center <- colMeans(pca_embedding[source_cells, , drop = FALSE])
      source_distance <- rowSums((
        pca_embedding[source_cells, , drop = FALSE] -
          matrix(source_center, nrow = length(source_cells),
                 ncol = ncol(pca_embedding), byrow = TRUE)
      )^2)
      second_root <- source_cells[which.min(source_distance)]
      two_root_time <- LineaLinkR:::ComputePseudotimeFuzzyMulti(
        nn_idx, nn_dists, c(start_cell, second_root),
        rank_normalize = rank_norm
      )
      source_time <- pmin(two_root_time[, 1L], two_root_time[, 2L]) /
        max(two_root_time[, 1L])
      source_weight <- branch_weights[, source_lineage]
      pseudotime <- pseudotime * (1 - source_weight) +
        source_time * source_weight
      pseudotime <- pseudotime / max(pseudotime)
      root_diagnostics$branch_source$second_root <- second_root
      updated_direction <- LineaLinkR::AuditPseudotimeDirection(
        pseudotime,
        potency = root_potency_vector,
        embedding = blind_emb,
        n_pc = llr_cfg$root_direction_n_pc %||% 3L,
        endpoint_fraction = llr_cfg$root_direction_endpoint_fraction %||% 0.05,
        min_potency_delta = llr_cfg$root_direction_min_potency_delta %||% 0.10,
        min_blind_delta = llr_cfg$root_direction_min_blind_delta %||% 0.05
      )
      direction_diagnostics[names(updated_direction)] <- updated_direction
      direction_status <- direction_diagnostics$direction_status
      direction_margin <- abs(direction_diagnostics$blind_score_margin)
      root_diagnostics$direction_diagnostics <- direction_diagnostics
      cat(sprintf("  [branch source: lineage %d cluster %d cell %d] ",
                  source_lineage, source_cluster, second_root))
    }
  }

  state_clock_conn <- knn_conn
  switch(as.character((donor_split_geometry_mode &&
                       geometry_excludes_lr) || geometry_centers_lr), "TRUE" = {
    structural_dir <- file.path(result_dir, "structural_geometry")
    dir.create(structural_dir, recursive = TRUE, showWarnings = FALSE)
    structural_rows <- switch(as.character(geometry_centers_lr),
      "TRUE" = seq_len(nrow(expr)),
      "FALSE" = setdiff(seq_len(nrow(expr)), lr_rows))
    structural_matrix <- file.path(structural_dir,
                                   sprintf("%s_matrix.mtx", ds_name))
    writeMM(expr[structural_rows, , drop = FALSE], structural_matrix)
    writeLines(rownames(expr)[structural_rows],
      file.path(structural_dir, sprintf("%s_genes.txt", ds_name)),
      useBytes = TRUE)
    center_args <- character()
    switch(as.character(geometry_centers_lr), "TRUE" = {
      donor_file <- file.path(structural_dir, "donor.txt")
      lr_gene_file <- file.path(structural_dir, "lr_genes.txt")
      writeLines(as.character(meta[[sample_col]]), donor_file, useBytes = TRUE)
      writeLines(rownames(expr)[lr_rows], lr_gene_file, useBytes = TRUE)
      center_args <- c("no_umap", "30", "donor_center_lr",
                       donor_file, lr_gene_file)
    })
    structural_status <- sys::exec_wait(
      CONFIG$exe$python,
      args = c("-X", "utf8", knn_script, ds_name, structural_dir,
               structural_matrix, as.character(knn_k), center_args),
      std_out = file.path(CONFIG$dir$log,
                          sprintf("llr_%s_structure.log", ds_name)),
      std_err = file.path(CONFIG$dir$log,
                          sprintf("llr_%s_structure.err", ds_name)))
    stopifnot(structural_status == 0L)
    structural_conn <- as(readMM(file.path(structural_dir, "knn_conn.mtx")),
                          "CsparseMatrix")
    state_clock_conn <- structural_conn
    geom_emb <- as.matrix(fread(file.path(structural_dir, "pca.csv"),
                                header = FALSE))
    lineage_topology <- LineaLinkR::ClusterLineageTopology(
      geom_emb, k = selector_topology_k,
      nstart = as.integer(llr_cfg$lineage_topology_nstart %||% 10L),
      iter_max = 500L, seed = CONFIG$seed)
    lineage_cl <- lineage_topology$labels
    start_cluster <- as.integer(lineage_cl[start_cell])
    knn_tree <- LineaLinkR::BuildKnnClusterTree(
      structural_conn, lineage_cl, start_cluster)
    lineage_paths <- LineaLinkR::ExtractLineagesFromKnnGraph(
      structural_conn, lineage_cl, start_cluster)$paths
    lineage_weights <- LineaLinkR::ComputeLineageWeightsFromKnnGraph(
      structural_conn, geom_emb, lineage_cl, start_cluster, tree = knn_tree)
    n_lineages <- ncol(lineage_weights)
    sibling_map <- switch(as.character(n_lineages > 1L),
      "TRUE" = LineaLinkR::BuildSiblingContrastMap(lineage_paths),
      "FALSE" = data.frame(focal_lineage = integer(),
                           sibling_lineage = integer(),
                           fork_cluster = integer(),
                           shared_prefix_length = integer()))
    knn_topology_diag <- list(
      n_tree_edges = nrow(knn_tree$edges),
      n_components = knn_tree$n_components,
      root_component = knn_tree$root_component,
      n_paths = length(lineage_paths),
      tree_edges = knn_tree$edges)
    lineage_cluster_source <- sprintf(
      "%s:kmeans:pca:requested_k=%d:occupied_k=%d",
      c("donor_no_lr", "donor_center_lr")[1L + as.integer(geometry_centers_lr)],
      lineage_topology$requested_k, lineage_topology$occupied_k)
    cat(sprintf("  [donor geometry: %d LR genes %s in branch graph; %d paths] ",
                length(lr_rows),
                c("excluded", "centered")[1L + as.integer(geometry_centers_lr)],
                n_lineages))
  })

  contracted_spurs <- NULL
  if (variant_id == "L1D0contDonorStructuralInterval" &&
      run_tag != "L1D0contDonorBalanceNoSpurInterval" &&
      geometry_excludes_lr && n_lineages > 1L) {
    paths <- as.data.table(stack(setNames(lineage_paths,
                                          seq_along(lineage_paths))))
    setnames(paths, c("values", "ind"), c("node", "lineage"))
    paths[, lineage := as.integer(as.character(lineage))]
    paths[, parent := shift(node), by = lineage]
    edges <- unique(paths[!is.na(parent), .(parent, child = node)])
    terminal_nodes <- setdiff(edges$child, edges$parent)
    node_time <- data.table(node = lineage_cl, time = pseudotime)[
      , .(median_time = median(time)), by = node]
    spurs <- merge(edges[child %in% terminal_nodes], node_time,
                   by.x = "parent", by.y = "node")
    setnames(spurs, "median_time", "parent_time")
    spurs <- merge(spurs, node_time, by.x = "child", by.y = "node")
    setnames(spurs, "median_time", "child_time")
    spurs <- spurs[child_time <= parent_time]
    contracted_spurs <- spurs
    if (nrow(spurs) > 0L) {
      drop_lineages <- unique(paths[node %in% spurs$child &
                                      node %in% terminal_nodes, lineage])
      remaining <- paths[!lineage %in% drop_lineages]
      remaining_edges <- unique(remaining[!is.na(parent),
                                          .(parent, child = node)])
      contract_nodes <- remaining_edges[parent %in% spurs$parent,
        .(n_children = uniqueN(child)), by = parent][n_children == 1L, parent]
      remaining <- remaining[!node %in% contract_nodes]
      lineage_paths <- split(remaining$node, remaining$lineage)
      lineage_weights <- lineage_weights[,
        as.integer(names(lineage_paths)), drop = FALSE]
      n_lineages <- ncol(lineage_weights)
      sibling_map <- LineaLinkR::BuildSiblingContrastMap(lineage_paths)
      cat(sprintf("  [donor topology: %d backward tips, %d contracted nodes, %d paths]\n",
        nrow(spurs), length(contract_nodes), n_lineages))
    }
  }

  crossing_topology_diag <- NULL
  switch(as.character(main_candidate && component_ablation != "no_crossing_repair"), "TRUE" = {
    crossing_variance <- colSums(scale(geom_emb, center = TRUE, scale = FALSE)^2) /
      (nrow(geom_emb) - 1L)
    crossing_dims <- which.max(head(crossing_variance, -1L) /
                                tail(crossing_variance, -1L))
    crossing_emb <- geom_emb[, seq_len(crossing_dims), drop = FALSE]
    crossing_size <- tabulate(lineage_cl)
    crossing_center <- rowsum(crossing_emb, lineage_cl) / crossing_size
    crossing_time <- as.numeric(rowsum(pseudotime, lineage_cl) / crossing_size)
    crossing_rank <- rank(crossing_time, ties.method = "first")
    crossing_root <- lineage_cl[start_cell]
    crossing_rank[crossing_root] <- 0
    crossing_data <- split(as.data.table(crossing_emb),
      factor(lineage_cl, levels = seq_along(crossing_size)))
    crossing_cov <- lapply(crossing_data, cov)
    crossing_pairs <- CJ(from = seq_along(crossing_size),
                         to = seq_along(crossing_size))[
      crossing_rank[from] < crossing_rank[to]]
    crossing_delta <- crossing_center[crossing_pairs$to, , drop = FALSE] -
      crossing_center[crossing_pairs$from, , drop = FALSE]
    crossing_pooled <- Map("+", crossing_cov[crossing_pairs$from],
                          crossing_cov[crossing_pairs$to])
    crossing_scaled <- do.call(rbind, Map(solve, crossing_pooled,
                                          asplit(crossing_delta, 1L)))
    crossing_eigen <- lapply(crossing_cov, eigen, symmetric = TRUE)
    crossing_tangent <- do.call(rbind,
      lapply(lapply(crossing_eigen, "[[", "vectors"), "[", , 1L))
    crossing_pairs[, `:=`(
      distance = sqrt(rowSums(crossing_delta^2)),
      shape_distance = sqrt(rowSums(crossing_delta * crossing_scaled)),
      normal_error = sqrt(pmax(0, rowSums(crossing_delta^2) -
        (rowSums(crossing_delta * crossing_tangent[from, , drop = FALSE])^2 +
         rowSums(crossing_delta * crossing_tangent[to, , drop = FALSE])^2) / 2))
    )]
    crossing_edges <- crossing_pairs[order(shape_distance, from), .SD[1L], by = to]
    crossing_leaves <- setdiff(seq_along(crossing_size), crossing_edges$from)
    crossing_forks <- crossing_edges[, .N, by = from][N > 1L, from]
    crossing_ends <- crossing_edges[to %in% crossing_leaves][from %in% crossing_forks,
      .(to, old_from = from, old_error = normal_error)]
    crossing_joins <- merge(crossing_pairs[from %in% crossing_leaves],
                            crossing_ends, by = "to")
    # Keep the previously validated geometric repair as the reference topology.
    crossing_original <- rbindlist(lapply(lineage_paths, data.table), idcol = "lineage")
    setnames(crossing_original, "V1", "node")
    crossing_original[, original_from := shift(node), by = lineage]
    crossing_original <- unique(crossing_original[!is.na(original_from),
      .(to = node, original_from)])
    crossing_joins <- merge(crossing_joins, crossing_original, by = "to")
    crossing_original_delta <- crossing_center[crossing_joins$to, , drop = FALSE] -
      crossing_center[crossing_joins$original_from, , drop = FALSE]
    crossing_joins[, original_error := sqrt(pmax(0,
      rowSums(crossing_original_delta^2) -
      (rowSums(crossing_original_delta * crossing_tangent[original_from, , drop = FALSE])^2 +
       rowSums(crossing_original_delta * crossing_tangent[to, , drop = FALSE])^2) / 2))]
    crossing_joins <- crossing_joins[from != original_from][
      normal_error < old_error][normal_error < original_error]
    crossing_joins[, gain := original_error - normal_error]
    crossing_joins <- crossing_joins[order(-gain, distance), .SD[1L], by = to]
    crossing_joins <- crossing_joins[order(-gain, distance), .SD[1L], by = from]
    crossing_var_fraction <- sum(crossing_variance[seq_len(crossing_dims)]) /
      sum(crossing_variance)
    crossing_topology_diag <- list(
      retained_dims = crossing_dims,
      retained_variance_fraction = crossing_var_fraction,
      joins = crossing_joins)
    switch(as.character(nrow(crossing_joins) > 0L &&
      crossing_var_fraction >= llr_cfg$crossing_min_retained_variance), "TRUE" = {
      switch(as.character(component_ablation != "no_terminal_continuation"), "TRUE" = {
        crossing_edges[crossing_joins, on = "to", from := i.from]
      })
      lineage_paths <- LineaLinkR:::rs_extract_lineages_from_edges(
        as.integer(crossing_edges$from), as.integer(crossing_edges$to),
        as.integer(length(crossing_size)), as.integer(crossing_root))
      crossing_nodes <- rbindlist(lapply(lineage_paths, data.table), idcol = "lineage")
      setnames(crossing_nodes, "V1", "node")
      crossing_nodes[, next_node := shift(node, type = "lead"), by = lineage]
      crossing_segments <- crossing_nodes[!is.na(next_node)]
      crossing_direction <- crossing_center[crossing_segments$next_node, , drop = FALSE] -
        crossing_center[crossing_segments$node, , drop = FALSE]
      crossing_grid <- CJ(cell = seq_len(nrow(crossing_emb)),
                          segment = seq_len(nrow(crossing_segments)))
      crossing_residual <- crossing_emb[crossing_grid$cell, , drop = FALSE] -
        crossing_center[crossing_segments$node[crossing_grid$segment], , drop = FALSE]
      crossing_along <- rowSums(crossing_residual *
        crossing_direction[crossing_grid$segment, , drop = FALSE]) /
        rowSums(crossing_direction[crossing_grid$segment, , drop = FALSE]^2)
      crossing_along <- pmin(1, pmax(0, crossing_along))
      crossing_grid[, `:=`(lineage = crossing_segments$lineage[segment],
        distance2 = rowSums((crossing_residual -
          crossing_direction[segment, , drop = FALSE] * crossing_along)^2))]
      crossing_distance <- crossing_grid[, .(distance2 = min(distance2)), by = .(cell, lineage)]
      crossing_distance[, distance2 := distance2 - min(distance2), by = cell]
      crossing_distance <- as.matrix(dcast(crossing_distance, cell ~ lineage,
        value.var = "distance2")[, -1L])
      lineage_weights <- exp(-crossing_distance / 2)
      lineage_weights <- lineage_weights / rowSums(lineage_weights)
      switch(as.character(component_ablation == "no_curve_assignment"), "TRUE" = {
        lineage_weights <- LineaLinkR::ComputeLineageWeightsFromKnnGraph(
          knn_conn, geom_emb, lineage_cl, crossing_root,
          tree = list(edges = as.data.frame(crossing_edges[, .(from, to)])))
      })
      n_lineages <- ncol(lineage_weights)
      sibling_map <- LineaLinkR::BuildSiblingContrastMap(lineage_paths)
      cat(sprintf("  [crossing continuation: %d repaired edges, %d paths]\n",
        nrow(crossing_joins), n_lineages))
    })
  })
  source(file.path(.llr_work_root, "experiments/run/compute_state_clock.R"),
    local = TRUE, encoding = "UTF-8")

  # One globally anchored clock is shared across reference paths and their common trunk.
  pseudotime_per_lineage <- NULL
  if (n_lineages >= 2) {
    cat("  Step 4b: Per-lineage pseudotime ... shared across the reference forest\n")
  }

  # Step 5: static LR scoring.
  cat("  Step 5/7: LR scores (static) ... ")
  sp <- CONFIG$dataset_species[[ds_name]]
  if (is.null(sp)) sp <- "human"
  lr_scores <- NULL; lr_db <- NULL; expr_lr <- NULL
  lr_input_scale <- ifelse(ds_name %in% names(CONFIG$dataset_lr_input_scale),
    CONFIG$dataset_lr_input_scale[ds_name], "counts")
  if (!sp %in% c("human", "mouse", "zebrafish", "celegans")) {
    cat(sprintf("SKIP (no configured LR database for species: %s)\n", sp))
  } else {
    expr_lr <- PrepareLRExpression(expr, as.numeric(llr_cfg$lr_scale_factor), lr_input_scale)
    lr_db <- BuildLRDB(expr_lr, species = sp)
    if (length(lr_db$ligand_rows) < 1L) {
      cat(sprintf("SKIP (no LR pairs matched for species=%s)\n", sp))
      lr_db <- NULL
    } else {
      lr_scores <- LineaLinkR::ComputeLRScores(expr_lr, celltypes, lr_db, n_perms = 100L)
      cat(sprintf("done [%d pairs]\n", nrow(lr_scores)))
    }
  }

  program_model <- NULL
  program_profiles <- NULL
  switch(as.character(program_lineage), "TRUE" = {
    program_excluded_genes <- switch(as.character(!is.null(lr_db)),
      "TRUE" = rownames(expr_lr)[unique(unlist(c(lr_db$ligand_rows,
        lr_db$receptor_rows), use.names = FALSE))], character())
    program_embedding <- switch(as.character(persistent_rna_states),
      "TRUE" = pca_embedding, "FALSE" = geom_emb)
    program_model <- LineaLinkR::RefineLineageProgramTopologyFromCounts(
      expr, lineage_cl, pseudotime, lineage_paths, lineage_weights,
      program_embedding, exclude_genes = program_excluded_genes,
      expression_scale = as.numeric(llr_cfg$lr_scale_factor), knn_conn = knn_conn)
    switch(as.character(persistent_rna_states), "TRUE" = {
      state_condition <- rep("all", length(lineage_cl))
      switch(as.character(nzchar(state_condition_column)), "TRUE" = {
        state_condition <- meta[[state_condition_column]]
      })
      program_model <- LineaLinkR::RefinePersistentLineageStatesFromCounts(
        expr, lineage_cl, pseudotime, program_model, pca_embedding, knn_conn,
        meta[[state_time_column]], state_condition,
        expression_scale = as.numeric(llr_cfg$lr_scale_factor))
      lineage_cl <- program_model$cluster_labels
      n_lineage_cl <- uniqueN(lineage_cl)
    })
    lineage_paths <- program_model$lineage_paths
    lineage_weights <- program_model$lineage_weights
    program_profiles <- program_model$program_basis$program_profiles
    n_lineages <- length(lineage_paths)
    sibling_map <- LineaLinkR::BuildSiblingContrastMap(lineage_paths)
    lineage_cluster_source <- paste0(lineage_cluster_source, "+program_continuity")
    switch(as.character(persistent_rna_states), "TRUE" = {
      lineage_cluster_source <- paste0(lineage_cluster_source, "+persistent_rna_states")
    })
  })

  lineage_reference_weights <- lineage_weights
  lineage_assignment_weights <- lineage_weights
  switch(as.character(enable_local_branch_residual), "TRUE" = {
    lineage_assignment_weights <- switch(as.character(continuous_fork_priming),
      "TRUE" = LineaLinkR::ComputeLineageWeightsLocalBranchResidual(
        geom_emb, lineage_cl, lineage_paths, priming = TRUE),
      "FALSE" = LineaLinkR::ComputeLineageWeightsLocalBranchResidual(
        geom_emb, lineage_cl, lineage_paths))
    n_lineages <- ncol(lineage_assignment_weights)
    switch(as.character(enable_root_geometric_gate), "TRUE" = {
      root_cells <- lineage_cl == lineage_paths[[1L]][1L]
      first_children <- vapply(lineage_paths, "[", integer(1L), 2L)
      children <- sort(unique(first_children))
      child_membership <- outer(first_children, children, "==") * 1.0
      old_branch_mass <- lineage_reference_weights[root_cells, , drop = FALSE] %*%
        child_membership
      new_branch_mass <- lineage_assignment_weights[root_cells, , drop = FALSE] %*%
        child_membership
      child_index <- match(first_children, children)
      lineage_assignment_weights[root_cells, ] <-
        lineage_assignment_weights[root_cells, , drop = FALSE] *
        old_branch_mass[, child_index, drop = FALSE] /
        new_branch_mass[, child_index, drop = FALSE]
    })
    lineage_weights <- lineage_assignment_weights
  })
  switch(as.character(enable_separate_priming_weights), "TRUE" = {
    lineage_weights <- lineage_reference_weights
    lineage_priming_weights <- LineaLinkR::CombineLineageMemberships(
      lineage_reference_weights, lineage_assignment_weights)
  })
  switch(component_ablation, local_branch_reference =, local_point_reference = {
    lineage_weights <- lineage_assignment_weights
  }, supported_branch_reference = {
    lineage_weights <- LineaLinkR::CombineLineageMemberships(
      lineage_reference_weights, lineage_assignment_weights)
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
  })
  lineage_membership <- switch(as.character(program_lineage),
    "TRUE" = program_model$prediction, max.col(lineage_weights, ties.method = "first"))

  # ── Step 6: Temporal LR scoring ───────────────────────────────────
  source(file.path(.llr_work_root, "experiments/run/compute_communication.R"),
    local = TRUE, encoding = "UTF-8")

  # Map pair indices to gene names.
  pair_gene_map <- NULL
  if (!is.null(lr_scores) && !is.null(lr_db) && nrow(lr_scores) > 0) {
    pair_gene_map <- lr_db$pair_gene_map
  }

  resident_sender_events <- NULL
  switch(as.character(resident_sender && nrow(fork_contrasts) > 0L &&
    !is.null(attr(fork_contrasts, "event_units"))), "TRUE" = {
    resident_stage <- switch(as.character(nzchar(resident_stage_column)),
      "TRUE" = meta[[resident_stage_column]], rep("snapshot", n_cells))
    resident_time <- switch(as.character(nzchar(resident_time_column)),
      "TRUE" = meta[[resident_time_column]], NULL)
    resident_condition <- switch(as.character(nzchar(resident_condition_column)),
      "TRUE" = meta[[resident_condition_column]], rep("all", n_cells))
    resident_sender_events <- LineaLinkR::ComputeResidentSenderEventsFromCells(
      expr_lr, celltypes, lineage_cl, lineage_paths,
      attr(fork_contrasts, "event_units"), knn_conn, pair_gene_map,
      attr(temporal_lr$scores, "hill_scale"), stage = resident_stage,
      context_order = sort(unique(resident_stage)), embedding = geom_emb,
      observed_time = resident_time, condition = resident_condition,
      min_cells = llr_cfg$min_n_eff %||% 5, min_pct = llr_cfg$min_pct %||% 0.05)
  })

  # Save the result.
  exp_meta <- LlrExpectedMeta(use_knn_topology)
  result <- list(
    dataset          = ds_name,
    timestamp        = Sys.time(),
    result_timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    run_attempt      = run_attempt,
    package_version  = as.character(packageVersion("LineaLinkR")),
    cells            = cells,
    start_cell       = start_cell,
    pseudotime_origin_cell = ifelse(main_candidate, which.min(pseudotime), ifelse(trunk_orientation_flip,
      trunk_time_origin_cell,
      ifelse(enable_component_axis_repair, which.min(pseudotime), start_cell))),
    start_method     = start_method,
    benchmark_protocol = exp_meta$benchmark_protocol,
    result_revision  = exp_meta$result_revision,
    run_tag          = as.character(run_tag),
    event_anchor_rule = event_anchor_val,
    event_scoring_mode = event_scoring_val,
    event_expression_mean_mode = event_expression_mean_mode,
    lr_detection_intensity_rho = lr_detection_intensity_rho,
    event_activity_mode = event_activity_val,
    merge_fork_windows = isTRUE(enable_merged_windows),
    knn_graph_topology = isTRUE(use_knn_topology),
    knn_topology_diag = knn_topology_diag,
    contracted_spurs = contracted_spurs,
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
    geometry_donor_signal = geometry_donor_signal,
    geometry_excludes_lr = geometry_excludes_lr,
    geometry_centers_lr = geometry_centers_lr,
    n_cells          = n_cells,
    n_genes          = nrow(expr),
    n_types          = n_types,
    n_lineages       = n_lineages,
    celltype_column   = ct_col,
    celltype_labels   = unname(ct_labels),
    celltypes         = celltypes,
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
    branch_clock_calibration = branch_clock_calibration,
    state_clock     = state_clock,
    pseudotime_per_lineage = pseudotime_per_lineage,
    lineage_weights  = lineage_weights,
    lineage_membership = lineage_membership,
    lineage_program_model = program_model,
    lineage_crossing_diagnostics = crossing_topology_diag,
    lineage_paths    = lineage_paths,
    sibling_map      = sibling_map,
    fork_bins        = fork_bins,
    lr_scores        = lr_scores,
    lr_input_scale   = unname(lr_input_scale),
    temporal_lr      = temporal_lr,
    shape_temporal_lr = shape_temporal_lr,
    reference_temporal_lr = reference_temporal_lr,
    sample_temporal_lr = sample_temporal_lr,
    calendar_temporal_lr = calendar_temporal_lr,
    calendar_fork_dynamics = calendar_fork_dynamics,
    calendar_inference_status = calendar_inference_status,
    entropy_raw      = entropy_raw,
    h_tilde          = h_tilde,
    bif_df           = bif_df,
    fork_contrasts   = fork_contrasts,
    fork_evidence    = fork_evidence,
    ligand_nominations = attr(fork_evidence, "ligand_nominations"),
    resident_sender_events = resident_sender_events,
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
    input_extract_dir = LlrExtractDir(ds_name),
    knn_precompute_dir = knn_dir,
    pair_gene_map    = pair_gene_map,
    config_snapshot  = c(exp_meta, list(
      seed      = CONFIG$seed,
      n_threads = CONFIG$n_threads,
      n_bins    = llr_cfg$n_bins %||% 12L,
      no_fork_control = llr_cfg$no_fork_control,
      prefork_weighting = llr_cfg$prefork_weighting,
        activity_transform = llr_cfg$event_activity_transform,
        context_error_weighting = llr_cfg$context_error_weighting,
        late_bin_merging_enabled = enable_late_bin_merging,
      persistence_operator = llr_cfg$persistence_operator,
      persistence_sibling_scope = llr_cfg$persistence_sibling_scope,
      cluster_input_policy = cluster_input_policy,
      lineage_weight_mode = CONFIG$llr$lineage_weight_mode,
      lineage_curve_weight_mode = CONFIG$llr$lineage_curve_weight_mode,
      lineage_persistence_weight_mode = ifelse(use_assignment_layout,
        "reference_geometry",CONFIG$llr$lineage_curve_weight_mode),
      lineage_priming_mode = CONFIG$llr$lineage_priming_mode,
      lineage_assignment_mode = CONFIG$llr$lineage_assignment_mode,
      root_attractor_selection = root_diagnostics$attractor_selection,
      selected_root_potency_mode = selected_root_potency_mode,
      algorithm = sprintf("LineaLinkR v%s %s rev=%s potency=%s geom=%s topology_k=%d",
        packageVersion("LineaLinkR"), exp_meta$benchmark_protocol,
        exp_meta$result_revision, exp_meta$root_potency_mode,
        exp_meta$lineage_embedding, exp_meta$lineage_topology_k)
    ))
  )
  switch(as.character(enable_separate_priming_weights), "TRUE" = {
    result$lineage_reference_weights <- lineage_reference_weights
    result$lineage_assignment_weights <- lineage_assignment_weights
    result$lineage_curve_weights <- lineage_curve_weights
    result$lineage_shape_weights <- lineage_shape_weights
    result$lineage_priming_weights <- lineage_priming_weights
  })

  # Write beside the final path, validate the serialized object by reopening it,
  # then atomically publish it. A truncated connection never becomes the
  # trusted complete.rds path.
  cache_tmp <- file.path(
    result_dir,
    sprintf(".linealinkr_complete.rds.tmp-%s", Sys.getpid())
  )
  saveRDS(result, cache_tmp)
  saved_check <- readRDS(cache_tmp)
  if (!LlrRdsComplete(saved_check)) {
    stop(sprintf("Serialized result failed completeness validation: %s", cache_tmp))
  }
  if (!file.rename(cache_tmp, cache_file)) {
    stop(sprintf("Could not publish atomic result cache: %s", cache_file))
  }
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
    if (is.null(res)) stop(sprintf("Requested dataset failed: %s", ds))
    all_results[[ds]] <- ds
    rm(res)
    gc(verbose = FALSE)
  }

  # Summarize only datasets requested and completed in this execution.
  summary_file <- file.path(OUT, "linealinkr_summary.csv")
  summary_paths <- file.path(OUT, ds_names, "linealinkr_complete.rds")
  expected_summary_meta <- LlrExpectedMeta()
  summary <- rbindlist(lapply(summary_paths, function(path) {
    r <- readRDS(path)
    saved <- r$config_snapshot %||% list()
    saved$benchmark_protocol <- r$benchmark_protocol %||% saved$benchmark_protocol
    saved$result_revision <- r$result_revision %||% saved$result_revision
    expected_lr_scale <- unname(ifelse(
      r$dataset %in% names(CONFIG$dataset_lr_input_scale),
      CONFIG$dataset_lr_input_scale[r$dataset], "counts"))
    current <- LlrRdsComplete(r) &&
      identical(r$lr_input_scale %||% "counts", expected_lr_scale) &&
      identical(
        as.character(saved$benchmark_protocol),
        expected_summary_meta$benchmark_protocol
      ) &&
      identical(
        as.character(saved$result_revision),
        expected_summary_meta$result_revision
      )
    stopifnot(current)
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
  cat(sprintf("Completed requested datasets: %d\n", length(all_results)))
}

RunLineaLinkRBenchmark()
