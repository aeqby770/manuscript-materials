REV <- normalizePath(Sys.getenv("LLR_WORK_ROOT",getwd()),winslash="/",mustWork=TRUE)
BASELINE_ROOT <- file.path(REV,"inputs/comparison")
# Evaluate the selected LineaLinkR 0.6.5 result on the fixed ordering panel using
# the submission metrics (Spearman/Pearson/RMSE on normalized pseudotime,
# ARI/Branch_F1 on argmax lineage weights vs true_branch where present).
# Reuse saved predictions: the submission table plus completed native real-data metrics.
Sys.setenv(OMP_NUM_THREADS = as.integer(Sys.getenv("LLR_THREADS", "2")))
data.table::setDTthreads(as.integer(Sys.getenv("LLR_THREADS", "2")))
suppressPackageStartupMessages(library(data.table))

RUN_VARIANT <- c(commandArgs(trailingOnly = TRUE), "Final065")[[1L]]
RUN_ALGORITHM <- paste0("LLR_0.6.5_", RUN_VARIANT)
EXCLUDED_OWN <- c("LLR_0.6.4")
WORKSPACES <- c(CurrentRefit = "server_e/current_default_refit/main", ObservedGraph = "topology_representation_repair/ordering24",
  Final065_NoSignal = "final065_ablation/no_coherent_clock_signal",
  Final065 = "final065", L1D0fix = "linealinkr_workspace_L1D0",
  FeasibleClock = "feasible_clock",
  Unified = "unified",
  VerifyRootfix = "verify_rootfix",
  VerifyPreRootfix = "verify_prerootfix",
  VerifyNoKm = "verify_nokm",
  Main065Rootfix = "main065_final",
  ClockRepair = "clock_repair",
  DepthClock = "clock_depth_final",
  UnifiedNoCal = "unified_no_cal",
  AnchoredClock = "anchored_clock",
  OverlapClock = "state_clock",
  Revalidated20261007 = "main",
  RevalidatedInherited20261007 = "main",
  RevalidatedNoPhase20261007 = "ablation/no_component_phase_alignment",
  RevalidatedNoKMeans20261007 = "ablation/no_kmeans_refinement",
  RevalidatedGraphK1020261007 = "ablation/graph_k_10",
  RevalidatedGraphK3020261007 = "ablation/graph_k_30",
  L1D0contPrecise = "linealinkr_workspace_L1D0_fixcont_precise",
  L1D0contFateSolve = "linealinkr_workspace_L1D0_fixcont_fatesolve",
  L1D0contTargetedDual = "linealinkr_workspace_L1D0_fixcont_targeted_timepoint_root_majority_leaf_local_branch_residual_scaled_root_geometric_gate_separate_priming_weights",
  L1D0cont065Main = "linealinkr_workspace_L1D0_fixcont_targeted_interval_timepoint_root_majority_leaf_timepoint_dominant_early_leaf_local_branch_residual_scaled_root_geometric_gate_separate_priming_weights_dropout_aware_event_mean",
  CurrentMain = "llr065_preserve_selected_root_full",
  Main065Integrated = "llr065_preserve_selected_root_full_donor_center_lr_geometry_structural_event_rank_connected_long_trunk_orientation",
  Main065Axis = "llr065_preserve_selected_root_full_donor_center_lr_geometry_structural_event_rank_connected_long_trunk_orientation_component_axis_repair",
  Main065Induction = "llr065_preserve_selected_root_full_donor_center_lr_geometry_structural_event_rank_connected_long_trunk_orientation_component_axis_repair_shared_induction_clock",
  Main065Phase = "llr065_preserve_selected_root_full_donor_center_lr_geometry_structural_event_rank_connected_long_trunk_orientation_component_axis_repair_shared_induction_clock_centroid_phase",
  Main065PhasePca2 = "llr065_preserve_selected_root_full_donor_center_lr_geometry_structural_event_rank_connected_long_trunk_orientation_component_axis_repair_shared_induction_clock_centroid_phase_pca2",
  Main065Final = "llr065_preserve_selected_root_full_donor_center_lr_geometry_structural_event_rank_connected_long_trunk_orientation_component_axis_repair_shared_induction_clock_centroid_phase_pca2_late_bin_persistence",
  FateSolveVerified = "linealinkr_workspace_L1D0_fixcont_fatesolve")
WS <- switch(RUN_VARIANT,
  RetentionProbe = file.path(REV, "results/algorithm_revision/root_selection"),
  file.path(REV, "cache", WORKSPACES[[RUN_VARIANT]]))
PREDICTION_FILE <- switch(RUN_VARIANT,
  RetentionProbe = "fate_retention_prediction.rds",
  ClockRepair = , DepthClock = "linealinkr_ordering.rds", "linealinkr_complete.rds")
Sys.setenv(LLR_WORK_ROOT = REV, LLR_PIPELINE = file.path(REV, "experiments"), LLR_ASSETS = REV)
source(file.path(REV, "experiments/config/llr_bootstrap.R"), encoding = "UTF-8")
llr_source_config()
MAIN065_ROOTFIX_REVISION <- "linealinkr065_feasible_anchored_expression_clock_rootfix20261016"
MAIN065_PREROOTFIX_REVISION <- "linealinkr065_feasible_anchored_expression_clock"
EXPECTED_REVISION <- switch(RUN_VARIANT,
  CurrentRefit = "linealinkr065_common_parent_population_nomination",
  ObservedGraph = "linealinkr065_common_parent_population_nomination_observed_graph_context",
  Final065_NoSignal = "linealinkr065_common_parent_observed_late_bin_fallback_no_coherent_clock_signal",
  Final065 = "linealinkr065_common_parent_population_nomination_assigned_layout_reference",
  ClockRepair = "linealinkr065_expression_clock_repair20261016",
  DepthClock = "linealinkr065_depth_adjusted_clock20261016",
  Main065Rootfix = MAIN065_ROOTFIX_REVISION,
  Main065Prerootfix = MAIN065_PREROOTFIX_REVISION,
  L1D0contPrecise = paste0(CONFIG$llr$result_revision,
    "_root_time_from_start_fuzzy_f64_accumulation"),
  RetentionProbe = "linealinkr065_retention_probe_ilu0_gmres",
  FateSolveVerified = , L1D0contFateSolve = paste0(CONFIG$llr$result_revision,
    "_ilu0_gmres_candidate_availability"), NULL)
switch(as.character(nzchar(Sys.getenv("LLR_EVAL_EXPECTED_REVISION"))),
  "TRUE"=EXPECTED_REVISION <- Sys.getenv("LLR_EVAL_EXPECTED_REVISION"))
WS <- Sys.getenv("LLR_EVAL_WORKSPACE",WS)
# Only late LR ranking changed between these revisions; the complete saved
# ordering coordinates and lineage weights are reused without refitting.
COORDINATE_REFERENCE_REVISION <- switch(RUN_VARIANT,
  Final065 = "linealinkr065_final_common_parent_native_persistence",
  Final065_NoSignal = "linealinkr065_final_common_parent_native_persistence_no_coherent_clock_signal",
  NULL)
CMP_TOL <- 1e-10 # CSV floating-point round trips are ties at this scale.
REFD <- file.path(BASELINE_ROOT,"results/paper_tables/reference")
SUB <- file.path(BASELINE_ROOT,"results/paper_tables/table1_main")
COMPLETED_PT <- file.path(BASELINE_ROOT, "results/paper_tables/completed_real_ordering.csv")
OUT <- Sys.getenv("LLR_EVAL_OUT_DIR", file.path(REV, "results", "evaluation", "final065"))

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
DATASETS <- BENCHMARK_ORDERING_DATASETS
RUN_STATUS <- switch(RUN_VARIANT, FateSolveVerified = , RetentionProbe = , L1D0contFateSolve =
  fread(file.path(OUT, paste0("ordering_run_status_",
    switch(RUN_VARIANT, FateSolveVerified = "L1D0contFateSolve", RUN_VARIANT), ".tsv"))),
  data.table(dataset = character(), exit_code = integer()))
FAILED_DATASETS <- RUN_STATUS[exit_code != 0L, dataset]
PREDICTION_INPUTS <- data.table(dataset = DATASETS,
  prediction_file = file.path(WS, DATASETS, PREDICTION_FILE),
  expected_revision = rep_len(c(EXPECTED_REVISION, NA_character_)[1L], length(DATASETS)),
  solver = switch(RUN_VARIANT, FateSolveVerified = "gmres", NA_character_))
switch(RUN_VARIANT, Main065Axis = , Main065Induction = {
  unchanged_connected <- c("dentategyrus", "gastrulation_erythroid",
    "human_fetal_liver_organoid", "human_fetal_pancreas", "larry_invitro",
    "lung", "zebrafish_neural_dev", "zebrafish_tracerseq")
  PREDICTION_INPUTS[dataset %in% unchanged_connected,
    prediction_file := file.path(REV, "cache",
      WORKSPACES[["Main065Integrated"]], dataset, PREDICTION_FILE)]
})
switch(RUN_VARIANT, Main065Phase = , Main065PhasePca2 = {
  affected_disconnected <- c("sim_3br_500c", "sim_3br_2000c",
    "sim_3br_5000c", "sim_3br_2000c_batch",
    "sim_3br_2000c_drop30", "sim_3br_2000c_drop50",
    "simulated_complex_3branch")
  PREDICTION_INPUTS[!dataset %in% affected_disconnected,
    prediction_file := file.path(REV, "cache",
      WORKSPACES[["Main065Induction"]], dataset, PREDICTION_FILE)]
  unchanged_connected <- c("dentategyrus", "gastrulation_erythroid",
    "human_fetal_liver_organoid", "human_fetal_pancreas", "larry_invitro",
    "lung", "zebrafish_neural_dev", "zebrafish_tracerseq")
  PREDICTION_INPUTS[dataset %in% unchanged_connected,
    prediction_file := file.path(REV, "cache",
      WORKSPACES[["Main065Integrated"]], dataset, PREDICTION_FILE)]
})
switch(RUN_VARIANT, FateSolveVerified = {
  # Both solvers target the same absorption equation and pass the original
  # probability tolerance. Selection depends only on numerical convergence,
  # never on the reference ordering or a dataset's performance.
  PREDICTION_INPUTS[dataset %in% FAILED_DATASETS, `:=`(
    prediction_file = file.path(REV, "results/algorithm_revision/root_selection",
      dataset, "fate_direct_prediction.rds"),
    expected_revision = "linealinkr065_sparse_lu_absorption", solver = "sparse_lu")]
  FAILED_DATASETS <- character()
})

ARI <- function(pred, truth) {
  tbl <- table(pred, truth)
  a <- sum(choose(tbl, 2))
  b <- sum(choose(rowSums(tbl), 2)) - a
  c2 <- sum(choose(colSums(tbl), 2)) - a
  den <- (a + b + a + c2) / 2 - (a + b) * (a + c2) / choose(sum(tbl), 2)
  if (den == 0) 0 else (a - (a + b) * (a + c2) / choose(sum(tbl), 2)) / den
}
BranchF1 <- function(pred, truth) {
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

# Submission-protocol ground-truth lookup.
.extract_ordinal_time <- function(vals) {
  s <- as.character(vals)
  num <- suppressWarnings(as.numeric(
    sub(".*?(-?[0-9]+(?:\\.[0-9]+)?).*", "\\1", s, perl = TRUE)))
  num[!grepl("[0-9]", s)] <- NA_real_
  num
}
FindGroundTruth <- function(ref, ds_name) {
  registered_gt_cols <- c(
    lung = "time_point", larry_invitro = "time_info",
    human_fetal_liver_organoid = "timepoint",
    human_fetal_pancreas = "pcw", zebrafish_tracerseq = "time_hours")
  cid <- if ("cell_id" %in% names(ref)) "cell_id" else names(ref)[1L]
  gt_cols <- c(unname(registered_gt_cols[ds_name]),
               "true_pseudotime", "timepoint", "day", "age(days)", "stage")
  gt_cols <- gt_cols[!is.na(gt_cols)]
  for (col in unique(gt_cols)) {
    if (!col %in% names(ref)) next
    vals <- ref[[col]]
    numeric_vals <- suppressWarnings(as.numeric(vals))
    if (sum(!is.na(numeric_vals)) < 10) {
      extracted <- .extract_ordinal_time(vals)
      if (sum(!is.na(extracted)) >= 10 &&
          length(unique(extracted[!is.na(extracted)])) >= 2) {
        numeric_vals <- extracted
      }
    }
    valid <- !is.na(numeric_vals)
    if (sum(valid) > 10 && length(unique(numeric_vals[valid])) >= 2) {
      return(list(column = col,
                  cell_id = as.character(ref[[cid]]),
                  values = as.numeric(numeric_vals)))
    }
  }
  NULL
}

switch(RUN_VARIANT, RevalidatedInherited20261007 = {
  inherited_datasets <- c(setdiff(names(CONFIG$dataset_registry),
    c(names(CONFIG$datasets_revision_sim), "GSE135893_IPF")), "cellrank_stemnet")
  inherited_ablations <- c("full", "no_lr_module", "no_temporal_binning",
    "no_bifurcation_module", "no_lineage_specific_scoring", "no_kmeans_refinement",
    "bins_5", "bins_10", "bins_12", "bins_20", "graph_k_10", "graph_k_30")
  stopifnot(length(inherited_datasets) == 36L, !anyDuplicated(inherited_datasets))
  inherited_grid <- CJ(dataset = inherited_datasets, ablation_id = inherited_ablations)
  inherited_grid[, prediction_file := file.path(REV, "cache",
    paste0("llr065_revalidated_20261007_", ablation_id), dataset,
    "linealinkr_complete.rds")]
  inherited_grid[, execution_status := fifelse(file.exists(prediction_file),
    "completed", "not_completed")]
  inherited_grid[execution_status == "not_completed" &
    ((dataset == "paul15" & ablation_id == "graph_k_10") |
      (dataset == "cellagentchat_scc" & ablation_id %in% c("graph_k_10", "graph_k_30"))),
    execution_status := "failed_no_eligible_root"]
  inherited_metrics <- inherited_grid[execution_status == "completed", {
    x <- readRDS(prediction_file)
    stopifnot(x$dataset == dataset, x$package_version == "0.6.5")
    stats <- data.table(n_cells = length(x$cells), n_lineages = ncol(x$lineage_weights),
      root_cell = x$start_cell, n_bins = x$config_snapshot$n_bins,
      root_evidence_mode = x$root_evidence_mode, root_input_source = x$root_input_source,
      root_input_mode = ifelse(grepl("^observed_timepoint_", x$root_evidence_mode),
        "observed_timepoint_root", "blind_data_driven"),
      n_lr_rows = NROW(x$lr_scores), n_temporal_rows = NROW(x$temporal_lr$scores),
      n_bifurcation_points = sum(x$bif_df$reason == "crossing", na.rm = TRUE),
      n_fork_rows = NROW(x$fork_contrasts),
      n_evaluable = sum(x$fork_contrasts$evaluable %in% TRUE),
      n_ranked = sum(is.finite(x$fork_contrasts$event_ranking_score)),
      ground_truth = NA_character_, n_time_used = NA_integer_,
      Spearman = NA_real_, Pearson = NA_real_, RMSE = NA_real_,
      branch_scope = NA_character_, n_branch_used = NA_integer_,
      ARI = NA_real_, Branch_F1 = NA_real_)
    switch(as.character(dataset %in% DATASETS), "TRUE" = {
      ref <- fread(file.path(REFD, sprintf("ref_%s.csv", dataset)))
      gt <- FindGroundTruth(ref, dataset)
      switch(as.character(!is.null(gt)), "TRUE" = {
        matched <- merge(data.table(cell_id = gt$cell_id, truth = gt$values),
          data.table(cell_id = as.character(x$cells), prediction = x$pseudotime),
          by = "cell_id")[is.finite(truth) & is.finite(prediction)]
        p <- (matched$prediction - min(matched$prediction)) / diff(range(matched$prediction))
        t <- (matched$truth - min(matched$truth)) / diff(range(matched$truth))
        stats[, `:=`(ground_truth = gt$column, n_time_used = nrow(matched),
          Spearman = cor(p, t, method = "spearman"), Pearson = cor(p, t),
          RMSE = sqrt(mean((p - t)^2)))]
      })
      switch(as.character("true_branch" %in% names(ref)), "TRUE" = {
        branch_ref <- ref
        scope <- "all_cells_diagnostic"
        switch(as.character(dataset %in% BENCHMARK_BRANCH_DATASETS), "TRUE" = {
          threshold <- unname(BENCHMARK_BRANCH_TERMINAL_MIN_PSEUDOTIME[dataset])
          branch_ref <- ref[is.finite(as.numeric(true_pseudotime)) &
            as.numeric(true_pseudotime) >= threshold]
          scope <- "submission_terminal_cells"
        })
        cid <- c("cell_id", names(ref)[1L])[c("cell_id", names(ref)[1L]) %in% names(ref)][1L]
        bm <- merge(branch_ref[, .(cell_id = as.character(get(cid)), true_branch)],
          data.table(cell_id = as.character(x$cells),
            pred = max.col(x$lineage_weights, ties.method = "first")), by = "cell_id")
        stats[, `:=`(branch_scope = scope, n_branch_used = nrow(bm),
          ARI = ARI(bm$pred, bm$true_branch), Branch_F1 = BranchF1(bm$pred, bm$true_branch))]
      })
    })
    rm(x)
    gc(FALSE)
    stats
  }, by = .(dataset, ablation_id, prediction_file)]
  inherited_metrics <- merge(inherited_grid, inherited_metrics,
    by = c("dataset", "ablation_id", "prediction_file"), all.x = TRUE)
  inherited_output <- file.path(REV,
    "results/evaluation",
    "inherited_ablation_metrics.csv")
  fwrite(inherited_metrics, inherited_output, bom = FALSE)
  print(inherited_metrics[, .(n_datasets = .N, completed = sum(execution_status == "completed"),
    time_scored = sum(is.finite(Spearman)), branch_scored = sum(is.finite(Branch_F1))),
    by = ablation_id])
  quit(save = "no", status = 0L)
})

pt_rows <- list()
br_rows <- list()
br_saved_rows <- list()
stale_datasets <- character()
for (ds in DATASETS) {
  if (ds %in% FAILED_DATASETS) { cat(ds, "FAILED RUN\n"); next }
  rds <- PREDICTION_INPUTS[dataset == ds, prediction_file]
  if (!file.exists(rds)) { cat(ds, "NO WORKSPACE\n"); next }
  ref <- fread(file.path(REFD, sprintf("ref_%s.csv", ds)))
  x <- readRDS(rds)
  expected_revision <- PREDICTION_INPUTS[dataset == ds, expected_revision]
  if (!is.na(expected_revision) &&
      !x$config_snapshot$result_revision %in% c(expected_revision, COORDINATE_REFERENCE_REVISION)) {
    stale_datasets <- c(stale_datasets, ds)
    cat(ds, "STALE WORKSPACE\n")
    next
  }
  gt <- FindGroundTruth(ref, ds)
  cid <- if ("cell_id" %in% names(ref)) "cell_id" else names(ref)[1L]
  llr_cells <- as.character(x$cells)
  pred_lin <- max.col(x$lineage_weights, ties.method = "first")

  if (!is.null(gt)) {
    m <- merge(data.table(cell_id = gt$cell_id, gt_pt = gt$values),
               data.table(cell_id = llr_cells, pt = x$pseudotime),
               by = "cell_id")
    ok <- is.finite(m$pt) & is.finite(m$gt_pt)
    p01 <- (m$pt[ok] - min(m$pt[ok])) / diff(range(m$pt[ok]))
    t01 <- (m$gt_pt[ok] - min(m$gt_pt[ok])) / diff(range(m$gt_pt[ok]))
    pt_rows[[ds]] <- data.table(
      dataset = ds, algorithm = paste0("LineaLinkR_0.6.5_", RUN_VARIANT),
      ground_truth = gt$column, n_aligned = nrow(m), n_used = sum(ok),
      root_input_mode = ifelse(grepl("^observed_timepoint_", x$root_evidence_mode),
                               "observed_timepoint_root", "blind_data_driven"),
      Spearman = cor(p01, t01, method = "spearman"),
      Pearson = cor(p01, t01, method = "pearson"),
      RMSE = sqrt(mean((p01 - t01)^2)))
  }
  if ("true_branch" %in% names(ref)) {
    branch_ref <- ref
    branch_scope <- "all_cells_diagnostic"
    if (ds %in% BENCHMARK_BRANCH_DATASETS) {
      threshold <- unname(BENCHMARK_BRANCH_TERMINAL_MIN_PSEUDOTIME[ds])
      branch_ref <- ref[is.finite(as.numeric(true_pseudotime)) & as.numeric(true_pseudotime) >= threshold]
      branch_scope <- "submission_terminal_cells"
    }
    bm <- merge(branch_ref[, .(cell_id = as.character(get(cid)),
                        true_branch = as.factor(true_branch))],
                data.table(cell_id = llr_cells, pred = pred_lin),
                by = "cell_id")
    br_rows[[ds]] <- data.table(
      dataset = ds, algorithm = paste0("LineaLinkR_0.6.5_", RUN_VARIANT),
      n_lineages = ncol(x$lineage_weights), n_aligned = nrow(bm),
      branch_scope = branch_scope,
      ARI = ARI(bm$pred, bm$true_branch),
      Branch_F1 = BranchF1(bm$pred, bm$true_branch))
    slingshot_dir <- file.path(BASELINE_ROOT, "results/by_algorithm/slingshot", ds)
    saved_slingshot_cells <- fread(file.path(slingshot_dir, "slingshot_pseudotime.csv"))
    saved_slingshot_weights <- as.matrix(fread(file.path(slingshot_dir, "slingshot_weights.csv")))
    stopifnot(nrow(saved_slingshot_weights)==nrow(saved_slingshot_cells),
      all(is.finite(saved_slingshot_weights)),all(rowSums(saved_slingshot_weights)>0))
    saved_assignments <- data.table(cell_id=as.character(saved_slingshot_cells$cell_id),
      algorithm="slingshot",pred=max.col(saved_slingshot_weights,ties.method="first"))
    saved_aligned <- merge(branch_ref[, .(cell_id = as.character(get(cid)), true_branch)],
      saved_assignments, by = "cell_id")
    stopifnot(all(saved_aligned[, .N, by = algorithm]$N == nrow(bm)))
    br_saved_rows[[ds]] <- saved_aligned[, .(dataset = ds,
      n_aligned = .N, branch_scope = branch_scope,
      ARI = ARI(pred, true_branch), Branch_F1 = BranchF1(pred, true_branch),
      row_present = TRUE, metric_source = "reevaluated_completed_prediction"), by = algorithm]
  }
  cat(ds, "done\n")
}

pt_new <- rbindlist(pt_rows)
br_new <- rbindlist(br_rows)
fwrite(pt_new, file.path(OUT, paste0("ordering_pseudotime_", RUN_VARIANT, ".csv")))
fwrite(br_new, file.path(OUT, paste0("ordering_branch_", RUN_VARIANT, ".csv")))

# Reuse submission metrics and the recovered native competitor predictions.
# Every dataset x algorithm combination is retained so missing predictions are
# visible instead of being silently removed from a method average.
sub_pt <- fread(file.path(SUB, "table1_pseudotime_linealinkr064_submission.csv"))
pt_old <- sub_pt[, .(
  dataset,
  algorithm = fifelse(algorithm == "LineaLinkR", "LLR_0.6.4", algorithm),
  n_aligned, Spearman, Spearman_signed, Pearson, RMSE,
  root_input_mode, direction_policy, row_present = TRUE,
  metric_source = "reused_submission"
)]
completed_pt <- fread(COMPLETED_PT)
stopifnot(!anyDuplicated(completed_pt[, .(dataset, algorithm)]))
pt_old <- rbindlist(list(
  pt_old[!completed_pt, on = .(dataset, algorithm)], completed_pt),
  use.names = TRUE, fill = TRUE)
pt_new2 <- pt_new[, .(
  dataset, algorithm = RUN_ALGORITHM, n_aligned, Spearman,
  Spearman_signed = Spearman, Pearson, RMSE,
  root_input_mode, direction_policy = "signed",
  row_present = TRUE, metric_source = paste0("new_065_", RUN_VARIANT)
)]
reference_065_pt <- switch(RUN_VARIANT, FateSolveVerified = , RetentionProbe = , L1D0contFateSolve = {
  fread(file.path(OUT, "ordering_pseudotime_L1D0fix.csv"))[, .(
    dataset, algorithm = "LLR_0.6.5_L1D0fix", n_aligned, Spearman,
    Spearman_signed = Spearman, Pearson, RMSE,
    root_input_mode = "blind_data_driven", direction_policy = "signed",
    row_present = TRUE, metric_source = "saved_065_L1D0fix")]
}, pt_old[0L])
pt_old <- rbindlist(list(pt_old, reference_065_pt), use.names = TRUE)
pt_old<-pt_old[!algorithm%chin%EXCLUDED_OWN]
pt_methods <- unique(c(pt_old$algorithm, RUN_ALGORITHM))
pt_grid <- CJ(dataset = DATASETS, algorithm = pt_methods, unique = TRUE)
pt_all <- merge(pt_grid, rbind(pt_old, pt_new2, fill = TRUE),
                by = c("dataset", "algorithm"), all.x = TRUE, sort = FALSE)
pt_all[, output_status := fcase(
  algorithm == RUN_ALGORITHM & dataset %in% FAILED_DATASETS, "failed_prediction",
  algorithm == RUN_ALGORITHM & dataset %in% stale_datasets, "stale_prediction",
  is.na(row_present), "missing_prediction",
  !is.na(output_status) & output_status != "scored", output_status,
  is.finite(Spearman) & is.finite(RMSE), "scored",
  default = "nonfinite_output"
)]
pt_all[, metric_source := fifelse(
  is.na(metric_source), "not_available", metric_source
)]
fwrite(pt_all, file.path(OUT, paste0("ordering_pseudotime_complete_", RUN_VARIANT, ".csv")))
switch(RUN_VARIANT, FateSolveVerified = fwrite(merge(PREDICTION_INPUTS,
  pt_all[algorithm == RUN_ALGORITHM, .(dataset, output_status)], by = "dataset"),
  file.path(OUT, "ordering_solver_sources_FateSolveVerified.csv")))
fwrite(dcast(pt_all, dataset ~ algorithm, value.var = "Spearman"),
       file.path(OUT, paste0("ordering_pseudotime_wide",
         switch(RUN_VARIANT, L1D0fix = "", paste0("_", RUN_VARIANT)), ".csv")))

sub_ari <- fread(file.path(SUB, "table1_lineage_linealinkr064_submission.csv"))
sub_f1 <- fread(file.path(SUB, "table1_branch_f1_linealinkr064_submission.csv"))
br_old <- merge(
  sub_ari[, .(dataset, algorithm = fifelse(algorithm == "LineaLinkR", "LLR_0.6.4", algorithm),
              ARI, row_present = TRUE)],
  sub_f1[, .(dataset, algorithm = fifelse(algorithm == "LineaLinkR", "LLR_0.6.4", algorithm),
             Branch_F1)],
  by = c("dataset", "algorithm"), all = TRUE
)
br_old[, metric_source := "reused_submission"]
br_saved <- rbindlist(br_saved_rows)
br_old <- rbindlist(list(br_old[!br_saved, on = .(dataset, algorithm)], br_saved), fill = TRUE)
fwrite(br_saved, file.path(OUT, paste0("ordering_branch_reused_predictions_", RUN_VARIANT, ".csv")))
reference_065_br <- switch(RUN_VARIANT, FateSolveVerified = , RetentionProbe = , L1D0contFateSolve = {
  fread(file.path(OUT, "ordering_branch_L1D0fix.csv"))[, .(
    dataset, algorithm = "LLR_0.6.5_L1D0fix", ARI, Branch_F1,
    row_present = TRUE, metric_source = "saved_065_L1D0fix")]
}, br_old[0L])
br_old <- rbindlist(list(br_old, reference_065_br), use.names = TRUE)
br_new2 <- br_new[, .(dataset, algorithm = RUN_ALGORITHM,
                      ARI, Branch_F1, branch_scope, row_present = TRUE,
                      metric_source = paste0("new_065_", RUN_VARIANT))]
br_old<-br_old[!algorithm%chin%EXCLUDED_OWN]
br_saved<-br_saved[!algorithm%chin%EXCLUDED_OWN]
br_methods <- unique(c(br_old$algorithm, RUN_ALGORITHM))
br_grid <- CJ(dataset = DATASETS, algorithm = br_methods, unique = TRUE)
br_all <- merge(br_grid, rbind(br_old, br_new2, fill = TRUE),
                by = c("dataset", "algorithm"), all.x = TRUE, sort = FALSE)
br_all[, branch_status := fcase(
  algorithm == RUN_ALGORITHM & dataset %in% FAILED_DATASETS, "failed_prediction",
  algorithm == RUN_ALGORITHM & dataset %in% stale_datasets, "stale_prediction",
  is.na(row_present), "missing_prediction",
  is.finite(ARI) & is.finite(Branch_F1), "scored",
  default = "nonfinite_output"
)]
br_all[, metric_source := fifelse(is.na(metric_source), "not_available", metric_source)]
fwrite(br_all, file.path(OUT, paste0("ordering_branch_complete_", RUN_VARIANT, ".csv")))

all_methods <- unique(c(pt_methods, br_methods))
full_grid <- CJ(dataset = DATASETS, algorithm = all_methods, unique = TRUE)
full <- merge(full_grid, pt_all, by = c("dataset", "algorithm"), all.x = TRUE)
full <- merge(full,
              br_all[, .(dataset, algorithm, ARI, Branch_F1, branch_status,
                         branch_metric_source = metric_source)],
              by = c("dataset", "algorithm"), all.x = TRUE)
fwrite(full, file.path(OUT, paste0("ordering_complete_comparison_", RUN_VARIANT, ".csv")))

PairwisePT <- function(algo) {
  a <- pt_all[algorithm == RUN_ALGORITHM & output_status == "scored",
              .(dataset, llr_spearman = Spearman, llr_rmse = RMSE,
                llr_native_spearman = Spearman_signed)]
  b <- pt_all[algorithm == algo & output_status == "scored",
              .(dataset, competitor_spearman = Spearman, competitor_rmse = RMSE,
                competitor_native_spearman = Spearman_signed)]
  z <- merge(a, b, by = "dataset")
  # Preserve the inherited submission values, while also comparing both
  # methods under identical native-sign and direction-free conventions.
  comparison <- rbindlist(list(
    z[, .(dataset, metric = "Spearman", policy = "inherited_submission_direction_rules",
      llr = llr_spearman, competitor = competitor_spearman, higher_better = TRUE)],
    z[, .(dataset, metric = "RMSE", policy = "inherited_submission_direction_rules",
      llr = llr_rmse, competitor = competitor_rmse, higher_better = FALSE)],
    z[, .(dataset, metric = "Spearman_signed", policy = "native_signed_both_methods",
      llr = llr_native_spearman, competitor = competitor_native_spearman, higher_better = TRUE)],
    z[, .(dataset, metric = "Spearman_absolute", policy = "absolute_both_methods",
      llr = abs(llr_native_spearman), competitor = abs(competitor_native_spearman),
      higher_better = TRUE)]))
  comparison[, advantage := (llr - competitor) * fifelse(higher_better, 1, -1)]
  comparison[, .(
    comparator = algo, n_common = .N, llr_wins = sum(advantage > CMP_TOL),
    ties = sum(abs(advantage) <= CMP_TOL), competitor_wins = sum(advantage < -CMP_TOL),
    llr_mean = mean(llr), competitor_mean = mean(competitor),
    mean_delta_llr_minus_competitor = mean(llr - competitor)
  ), by = .(metric, policy)]
}
PairwiseBranch <- function(algo) {
  a <- br_all[algorithm == RUN_ALGORITHM & branch_status == "scored" &
                dataset %in% BENCHMARK_BRANCH_DATASETS,
              .(dataset, llr_ari = ARI, llr_f1 = Branch_F1)]
  b <- br_all[algorithm == algo & branch_status == "scored" &
                dataset %in% BENCHMARK_BRANCH_DATASETS,
              .(dataset, competitor_ari = ARI, competitor_f1 = Branch_F1)]
  z <- merge(a, b, by = "dataset")
  data.table(
    comparator = algo, metric = c("ARI", "Branch_F1"), n_common = nrow(z),
    llr_wins = c(sum(z$llr_ari > z$competitor_ari + CMP_TOL),
                 sum(z$llr_f1 > z$competitor_f1 + CMP_TOL)),
    ties = c(sum(abs(z$llr_ari - z$competitor_ari) <= CMP_TOL),
             sum(abs(z$llr_f1 - z$competitor_f1) <= CMP_TOL)),
    competitor_wins = c(sum(z$llr_ari < z$competitor_ari - CMP_TOL),
                         sum(z$llr_f1 < z$competitor_f1 - CMP_TOL)),
    llr_mean = c(mean(z$llr_ari), mean(z$llr_f1)),
    competitor_mean = c(mean(z$competitor_ari), mean(z$competitor_f1)),
    mean_delta_llr_minus_competitor = c(
      mean(z$llr_ari - z$competitor_ari),
      mean(z$llr_f1 - z$competitor_f1)))
}
pt_pairwise <- rbindlist(lapply(setdiff(pt_methods, RUN_ALGORITHM), PairwisePT),
                          fill = TRUE)
branch_pairwise <- rbindlist(lapply(setdiff(br_methods, RUN_ALGORITHM), PairwiseBranch),
                             fill = TRUE)
fwrite(pt_pairwise, file.path(OUT, paste0("ordering_pairwise_pseudotime_", RUN_VARIANT, ".csv")))
fwrite(branch_pairwise, file.path(OUT, paste0("ordering_pairwise_branch_", RUN_VARIANT, ".csv")))

# Separate a global/root-direction problem from within-branch ordering. This
# is diagnostic only: it does not alter the L1D0fix result or add a new model.
time_diag <- rbindlist(lapply(names(pt_rows), function(ds) {
  ref <- fread(file.path(REFD, sprintf("ref_%s.csv", ds)))
  if (!"true_branch" %in% names(ref)) return(NULL)
  x <- readRDS(PREDICTION_INPUTS[dataset == ds, prediction_file])
  gt <- FindGroundTruth(ref, ds)
  cid <- if ("cell_id" %in% names(ref)) "cell_id" else names(ref)[1L]
  m <- merge(
    ref[, .(cell_id = as.character(get(cid)),
            true_branch = as.factor(true_branch))],
    data.table(cell_id = as.character(x$cells), pt = as.numeric(x$pseudotime)),
    by = "cell_id"
  )
  mt <- merge(
    data.table(cell_id = gt$cell_id, gt_pt = gt$values),
    data.table(cell_id = as.character(x$cells), pt = as.numeric(x$pseudotime)),
    by = "cell_id"
  )
  ok <- is.finite(mt$pt) & is.finite(mt$gt_pt)
  m <- merge(m, data.table(cell_id = gt$cell_id, gt_pt = gt$values),
             by = "cell_id")
  within <- m[is.finite(pt) & is.finite(gt_pt),
              .(n_cells = .N, within_spearman = ifelse(
                  uniqueN(pt) > 1L & uniqueN(gt_pt) > 1L,
                  cor(pt, gt_pt, method = "spearman"), NA_real_)), by = true_branch]
  data.table(
    dataset = ds,
    global_spearman = cor(mt$pt[ok], mt$gt_pt[ok], method = "spearman"),
    median_within_branch_spearman = median(within$within_spearman, na.rm = TRUE),
    min_within_branch_spearman = as.numeric(quantile(within$within_spearman, 0, na.rm = TRUE)),
    max_within_branch_spearman = as.numeric(quantile(within$within_spearman, 1, na.rm = TRUE)),
    negative_within_branch_fraction = mean(within$within_spearman < 0, na.rm = TRUE),
    n_branches = nrow(within),
    global_vs_within_sign_mismatch = sign(cor(mt$pt[ok], mt$gt_pt[ok],
                                               method = "spearman")) !=
      sign(median(within$within_spearman, na.rm = TRUE))
  )
}), fill = TRUE)
fwrite(time_diag, file.path(OUT, paste0("ordering_time_diagnosis_", RUN_VARIANT, ".csv")))
best_signed <- pt_all[output_status == "scored" & !startsWith(algorithm, "LLR_") & is.finite(Spearman_signed),
  .(n_applicable_competitors = uniqueN(algorithm),
    strongest_signed_method = algorithm[which.max(Spearman_signed)],
    strongest_signed_spearman = max(Spearman_signed)), by = dataset]
best_rmse <- pt_all[output_status == "scored" & !startsWith(algorithm, "LLR_") & is.finite(RMSE),
  .(strongest_rmse_method = algorithm[which.min(RMSE)], strongest_rmse = min(RMSE)),
  by = dataset]
best_pt <- merge(best_signed, best_rmse, by = "dataset", all = TRUE)
main_pt <- merge(best_pt, pt_new[, .(dataset, n_aligned, Spearman, RMSE)],
  by = "dataset", all.x = TRUE)
main_pt[, `:=`(
  signed_gap = Spearman - strongest_signed_spearman,
  rmse_gap = RMSE - strongest_rmse
)]
fwrite(main_pt, file.path(OUT,
  paste0("ordering_vs_strongest_", RUN_VARIANT, ".csv")))
best_branch <- br_old[!startsWith(algorithm, "LLR_") & is.finite(ARI) & is.finite(Branch_F1),
  .(strongest_ari_method = algorithm[which.max(ARI)], strongest_ari = max(ARI),
    strongest_f1_method = algorithm[which.max(Branch_F1)], strongest_branch_f1 = max(Branch_F1)),
  by = dataset]
main_branch <- merge(best_branch, br_new[, .(dataset, ARI, Branch_F1, branch_scope)],
  by = "dataset", all.x = TRUE)
main_branch[, `:=`(
  ari_gap = ARI - strongest_ari,
  branch_f1_gap = Branch_F1 - strongest_branch_f1
)]
fwrite(main_branch, file.path(OUT,
  paste0("branch_vs_strongest_", RUN_VARIANT, ".csv")))
real_summary <- pt_all[dataset %in% completed_pt$dataset,
  .(n_datasets = .N, n_scored = sum(output_status == "scored"),
    n_failed = sum(output_status != "scored" & output_status != "missing_prediction"),
    n_missing = sum(output_status == "missing_prediction"),
    mean_Spearman_signed = mean(fifelse(output_status == "scored", Spearman_signed, NA_real_)),
    mean_Spearman_absolute = mean(fifelse(output_status == "scored", abs(Spearman_signed), NA_real_)),
    mean_RMSE = mean(fifelse(output_status == "scored", RMSE, NA_real_))),
  by = algorithm][order(-mean_Spearman_signed)]
fwrite(real_summary, file.path(OUT,
  paste0("ordering_real_summary_", RUN_VARIANT, ".csv")))
cat("Ordering evaluation saved:", nrow(pt_new), "datasets\n")
