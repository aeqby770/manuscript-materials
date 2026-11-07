WORK_ROOT <- normalizePath(Sys.getenv("LLR_WORK_ROOT",getwd()),winslash="/",mustWork=TRUE)
READONLY_REVISION_ROOT <- file.path(WORK_ROOT,"inputs/revision")
CONFIG <- list(
  revision = WORK_ROOT,
  scorecard = Sys.getenv("LLR_EVAL_RESULT_DIR", "results/evaluation/feasible_clock"),
  paper = file.path(WORK_ROOT,"results/paper"),
  workspace = Sys.getenv("LLR_OUTPUT_WORKSPACE", "feasible_clock"), library = "library",
  min_weight = 1e-8, min_n_eff = 5, tolerance = 1e-10,
  threads = as.integer(Sys.getenv("LLR_THREADS", "2"))
)
Sys.setenv(OMP_NUM_THREADS = CONFIG$threads)
library(data.table)
setDTthreads(CONFIG$threads)
scorer_library <- Sys.getenv("LLR_PACKAGE_LIBRARY",file.path(CONFIG$revision,CONFIG$library))
.libPaths(c(scorer_library,file.path(CONFIG$revision,CONFIG$library),.libPaths()))
library(LineaLinkR,lib.loc=scorer_library)
dataset <- commandArgs(trailingOnly = TRUE)[1L]
CONFIG$dataset <- dataset
base <- file.path(CONFIG$revision, CONFIG$scorecard)
out <- file.path(base, "dynamic_methods", dataset)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
x <- readRDS(file.path(CONFIG$revision, "cache", CONFIG$workspace, dataset, "linealinkr_complete.rds"))
detail <- fread(file.path(base, "lrsim_contrast_detail.csv"))[
  dataset == CONFIG$dataset][, `:=`(
    score = ranking_score, score_mode = "LineaLinkR_event_score")]
truth <- fread(file.path(READONLY_REVISION_ROOT, "sim", paste0(dataset, "_truth.csv")))
tra_dir <- file.path(CONFIG$revision, "results/by_algorithm/TraSig", dataset)
export_dir <- c(commandArgs(trailingOnly = TRUE)[-1L], tra_dir)[1L]
native <- fread(file.path(export_dir, "trasig_selected_task_pairs.csv"))
curves <- fread(file.path(export_dir, "trasig_temporal.csv"))
inputs <- fread(file.path(tra_dir, "trasig_cell_inputs.csv"))[match(x$cells, cell_id)]
stopifnot(identical(inputs$cell_id, x$cells))
pairs <- as.data.table(x$pair_gene_map)[,
  .(pair_idx, ligand = ligand_genes, receptor = receptor_genes)]
common <- merge(pairs, unique(native[, .(ligand, receptor)]), by = c("ligand", "receptor"))
fwrite(merge(truth, common[, .(ligand, receptor, in_common_resource = TRUE)],
  by = c("ligand", "receptor"), all.x = TRUE), file.path(out, "resource_coverage.csv"))

# Native LR-pair nomination: collapse paths/events without using the true branch.
# This task does not claim recovery of a specific fork or a transient phase.
llr_pairs <- detail[sender_direction %in% TRUE & evaluable %in% TRUE &
  duplicate_emission %in% c(FALSE, NA) & is.finite(score),
  .(score = max(score)), by = .(ligand, receptor)][, method := "LineaLinkR_event_max"]
# Same candidate filter, same windows and same time coordinate; only the
# sibling-branch contrast is removed; raw induction and the matched score remain distinct.
llr_focal_pairs <- detail[sender_direction %in% TRUE & evaluable %in% TRUE &
  duplicate_emission %in% c(FALSE, NA) & is.finite(post_fork_change),
  .(score = max(pmax(0, post_fork_change))), by = .(ligand, receptor)][
    , method := "LineaLinkR_focal_only"]
llr_matched_focal_pairs <- detail[sender_direction %in% TRUE & evaluable %in% TRUE &
  duplicate_emission %in% c(FALSE, NA) & is.finite(no_fork_ranking_score),
  .(score = max(no_fork_ranking_score)), by = .(ligand, receptor)][
    , method := "LineaLinkR_matched_no_fork"]
tra_pairs <- native[source_type == "SRC" & target_type == "TGT",
  .(score = max(observed_score)), by = .(ligand, receptor)][, method := "TraSig_native_max"]
static_pairs <- merge(as.data.table(x$lr_scores)[source == 1L & target == 2L], pairs,
  by = "pair_idx")[, .(ligand, receptor, score, method = "pooled_LR")]
pair_scores <- rbindlist(list(llr_pairs, llr_focal_pairs, llr_matched_focal_pairs, tra_pairs, static_pairs))
pair_grid <- truth[, .(method = c("LineaLinkR_event_max", "LineaLinkR_focal_only", "LineaLinkR_matched_no_fork", "TraSig_native_max", "pooled_LR")), by = .(ligand, receptor, true_positive, role)]
pair_grid <- merge(pair_grid, pair_scores, by = c("ligand", "receptor", "method"), all.x = TRUE)
pair_grid <- merge(pair_grid, common[, .(ligand, receptor, in_common_resource = TRUE)],
  by = c("ligand", "receptor"), all.x = TRUE)
pair_grid[, dataset := CONFIG$dataset]
fwrite(pair_grid, file.path(out, "native_pair_predictions.csv"))

# Controlled projection: retain each method's fitted time and curves, but use
# the same saved LineaLinkR cells, tree nodes, child weights and pre/post windows.
units <- attr(x$fork_contrasts, "event_units")$units
projection_source <- Sys.getenv("LLR_REUSE_PROJECTED_FROM","")
reuse_grid <- FALSE
if (nzchar(projection_source) && !is.null(x$upstream_result)) {
  parent <- readRDS(x$upstream_result)
  reuse_grid <- identical(x$cells,parent$cells) && identical(x$celltypes,parent$celltypes) &&
    identical(x$lineage_clusters,parent$lineage_clusters) && identical(x$lineage_paths,parent$lineage_paths) &&
    isTRUE(all.equal(x$lineage_weights,parent$lineage_weights,tolerance=1e-12)) &&
    isTRUE(all.equal(x$pseudotime,parent$pseudotime,tolerance=1e-12)) &&
    isTRUE(all.equal(units,attr(parent$fork_contrasts,"event_units")$units,tolerance=1e-12)) &&
    identical(as.data.frame(x$pair_gene_map),as.data.frame(parent$pair_gene_map))
  rm(parent)
}
old_grid <- file.path(projection_source,dataset,"projected_context_scores.csv")
fwrite(data.table(dataset=dataset,source=projection_source,matching_geometry=as.logical(reuse_grid),
  TraSig_grid_reused=reuse_grid&&file.exists(old_grid)),file.path(out,"projection_binding.csv"))
if (reuse_grid && file.exists(old_grid)) {
  grid <- fread(old_grid)
} else {
unit_table <- data.table(unit_id = seq_along(units), unit = units)
children <- unit_table[, .(child = unit[[1L]]$children, descendants = unit[[1L]]$desc), by = unit_id]
child_lineages <- children[, .(leaf = unlist(descendants)), by = .(unit_id, child)]
child_columns <- children[, .(child, lineage = seq_len(.N)), by = unit_id]
weights <- as.data.table(x$lineage_weights / rowSums(x$lineage_weights))
setnames(weights, as.character(seq_len(ncol(weights))))
weights[, cell_index := .I]
weights <- melt(weights, id.vars = "cell_index", variable.name = "leaf", value.name = "weight")
weights[, leaf := as.integer(as.character(leaf))]
weights <- merge(weights[weight > 0], child_lineages, by = "leaf", allow.cartesian = TRUE)
weights <- weights[, .(weight = sum(weight)), by = .(unit_id, child, cell_index)]
weights <- weights[weight >= CONFIG$min_weight]
cell_meta <- data.table(cell_index = seq_along(x$cells), node = x$lineage_clusters,
  cell_type = x$celltypes, path = inputs$cell_path, time = round(inputs$cell_time, 2))
pre_nodes <- unit_table[, .(node = unit[[1L]]$segment_pre_nodes), by = unit_id]
pre_nodes <- merge(children[, .(unit_id, child)], pre_nodes, by = "unit_id", allow.cartesian = TRUE)[, bin := 1L]
post_nodes <- unit_table[, .(child = unit[[1L]]$children, nodes = unit[[1L]]$segment_post_nodes), by = unit_id]
post_nodes <- post_nodes[, .(node = unlist(nodes)), by = .(unit_id, child)][, bin := 3L]
weights <- merge(weights, cell_meta, by = "cell_index")
weights <- merge(weights, rbindlist(list(pre_nodes, post_nodes)), by = c("unit_id", "child", "node"))
# The observed parent pool is shared across its children, matching the main estimand.
parent_mass <- weights[bin==1L,.(parent_mass=sum(weight)),by=.(unit_id,cell_index)]
shared_pre <- merge(pre_nodes,cell_meta,by="node",allow.cartesian=TRUE)
shared_pre <- merge(shared_pre,parent_mass,by=c("unit_id","cell_index"))
shared_pre <- merge(shared_pre,children[,.(n_children=.N),by=unit_id],by="unit_id")
shared_pre[,weight:=parent_mass/n_children]
weights <- rbindlist(list(weights[bin==3L],shared_pre[,names(weights),with=FALSE]),use.names=TRUE)

support <- weights[, .(mass = sum(weight), n_eff = sum(weight)^2 / sum(weight^2)),
  by = .(unit_id, child, bin, cell_type)]
sender_mix <- weights[, .(source_mass = sum(weight)), by = .(unit_id, child, bin, source = cell_type, source_path = path)]
sender_mix[, source_fraction := source_mass / sum(source_mass), by = .(unit_id, child, bin, source)]
receiver_mix <- weights[, .(receiver_mass = sum(weight)),
  by = .(unit_id, child, bin, target = cell_type, target_path = path, time)]
mixture <- merge(sender_mix, receiver_mix, by = c("unit_id", "child", "bin"), allow.cartesian = TRUE)
mixture[, path_pair := paste(source_path, target_path, sep = "_")]
curves <- merge(curves, common, by = c("ligand", "receptor"))
curve_values <- curves[, .(path_pair, time, pair_idx, score, source_cells_in_window, target_cells_in_window)]
stopifnot(!anyDuplicated(curve_values[, .(path_pair, time, pair_idx)]))
projected <- merge(mixture, curve_values, by = c("path_pair", "time"), all.x = TRUE, allow.cartesian = TRUE)
projected[, available := is.finite(score) & source_cells_in_window > 0 & target_cells_in_window > 0]
context <- c("unit_id", "child", "source", "target", "bin")
projected <- projected[, .(
  weighted_score = sum(fifelse(available, score, 0) * source_fraction * receiver_mass),
  supported_mass = sum(fifelse(available, source_fraction * receiver_mass, 0))
), by = c("pair_idx", context)]
grid <- unique(detail[pair_idx %in% common$pair_idx, .(unit_id, child, pair_idx, source, target)])[
  , .(bin = c(1L, 3L)), by = .(unit_id, child, pair_idx, source, target)]
grid <- merge(grid, child_columns, by = c("unit_id", "child"))
grid <- merge(grid, support[, .(unit_id, child, bin, source = cell_type, source_n_eff = n_eff)],
  by = c("unit_id", "child", "bin", "source"), all.x = TRUE)
grid <- merge(grid, support[, .(unit_id, child, bin, target = cell_type, receiver_n_eff = n_eff, mass)],
  by = c("unit_id", "child", "bin", "target"), all.x = TRUE)
grid <- merge(grid, projected, by = c("pair_idx", context), all.x = TRUE)
grid[, `:=`(score = weighted_score / mass, coverage = supported_mass / mass,
  w_eff = pmin(source_n_eff, receiver_n_eff), comp = fifelse(bin == 1L, "pre", "post"))]
grid[, row_valid := is.finite(w_eff) & w_eff >= CONFIG$min_n_eff &
  is.finite(coverage) & coverage >= 1 - CONFIG$tolerance]
grid[row_valid %in% FALSE, score := NA_real_]
}
fwrite(grid, file.path(out, "projected_context_scores.csv"))
events <- grid[, getFromNamespace(".event_unit_contrasts", "LineaLinkR")(
  .SD, units[[unit_id]], unit_id, "all", 1L, activity_floor = 0,
  activity_mode = "continuous")[, !"unit_id"], by = unit_id]
keys <- c("unit_id", "event_node", "child", "pair_idx", "source", "target")
metadata <- c(keys, "dataset", "ligand", "receptor", "source_label", "target_label",
  "sender_direction", "role", "true_positive", "donor_confounded", "truth_focal_branch",
  "truth_other_branch", "focal_set", "other_set", "duplicate_emission", "event_conflict",
  "unit_kind", "hit_target_key", "truth_row", "mirror_row", "resource_extra", "null_row",
  "false_fork_row", "prediction_key")
tra_detail <- merge(detail[, ..metadata], events, by = keys, all.x = TRUE)
tra_detail[, `:=`(score = sibling_contrast, score_mode = "TraSig_projected_contrast")]
tra_detail[is.na(reason), reason := "not_in_trasig_resource"]
tra_detail[is.na(evaluable), evaluable := FALSE]
event_predictions <- rbindlist(list(detail, tra_detail), fill = TRUE)
event_predictions[, common_resource := pair_idx %in% common$pair_idx]
fwrite(event_predictions, file.path(out, "projected_event_predictions.csv"))
print(pair_grid[, .(truth_pairs = .N, scored = sum(is.finite(score)),
  common_resource = sum(in_common_resource %in% TRUE)), by = method])
print(event_predictions[sender_direction %in% TRUE,
  .(events = .N, evaluable = sum(evaluable), true_evaluable = sum(truth_row & evaluable)), by = score_mode])

# Keep full truth denominators and pair/event tasks separate. Primary score
# comparisons use identical candidates; coverage outside that set is retained.
joint_pairs <- pair_grid[, .(jointly_scored = all(is.finite(score))), by = .(ligand, receptor)]
pair_ranked <- merge(pair_grid, joint_pairs, by = c("ligand", "receptor"))
pair_ranked <- rbindlist(list(
  copy(pair_ranked[is.finite(score)])[, evaluation_scope := "own_available"],
  copy(pair_ranked[jointly_scored %in% TRUE])[, evaluation_scope := "common_candidates"]
))
pair_ranked[, `:=`(truth_row = true_positive, task = "LR_pair_nomination",
  prediction_key = paste(ligand, receptor, sep = "|"))]
event_ranked <- event_predictions[sender_direction %in% TRUE &
  duplicate_emission %in% c(FALSE, NA) & resource_extra %in% c(FALSE, NA) &
  evaluable %in% TRUE & is.finite(score)]
setorder(event_ranked, score_mode, prediction_key, -score, lineage, pair_idx)
event_ranked <- event_ranked[, .SD[1L], by = .(score_mode, prediction_key)]
joint_events <- event_ranked[, .(jointly_scored = uniqueN(score_mode) == 2L), by = prediction_key]
event_ranked <- merge(event_ranked, joint_events, by = "prediction_key")
event_ranked <- rbindlist(list(
  copy(event_ranked)[, evaluation_scope := "own_available"],
  copy(event_ranked[jointly_scored %in% TRUE])[, evaluation_scope := "common_candidates"]
))
event_ranked[, `:=`(method = score_mode, task = "projected_fork_event")]
ranking_fields <- c("dataset", "task", "method", "evaluation_scope", "prediction_key", "truth_row", "score")
ranked <- rbindlist(list(pair_ranked[, ..ranking_fields], event_ranked[, ..ranking_fields]))
setorder(ranked, task, method, evaluation_scope, -score, prediction_key)
ranked[, candidate_rank := seq_len(.N), by = .(task, method, evaluation_scope)]
n_targets <- uniqueN(truth[true_positive %in% TRUE, .(ligand, receptor, focal_branch, other_branch)])
metrics <- ranked[, .(n_candidates = .N, n_positive = sum(truth_row), n_null = sum(!truth_row),
  n_targets = n_targets, precision_at_k = sum(truth_row[candidate_rank <= n_targets]) / n_targets),
  by = .(dataset, task, method, evaluation_scope)]
eligible <- metrics[n_positive >= 2L & n_null >= 2L, .(task, method, evaluation_scope)]
areas <- merge(ranked, eligible, by = c("task", "method", "evaluation_scope"))[, .(
  auprc = PRROC::pr.curve(scores.class0 = score[truth_row],
    scores.class1 = score[!truth_row], curve = FALSE)$auc.integral),
  by = .(task, method, evaluation_scope)]
metrics <- merge(metrics, areas, by = c("task", "method", "evaluation_scope"), all.x = TRUE)
metrics[n_targets == 0L, precision_at_k := NA_real_]
metrics[, inference_status := paste(unique(native$inference_status), collapse = "+")]
fwrite(ranked, file.path(out, "ranked_predictions.csv"))
fwrite(metrics, file.path(out, "metrics.csv"))
print(metrics[evaluation_scope == "common_candidates",
  .(task, method, n_candidates, n_positive, n_targets, auprc, precision_at_k)])

metric_files <- list.files(file.path(base, "dynamic_methods"),
  pattern = "^metrics[.]csv$", recursive = TRUE, full.names = TRUE)
all_metrics <- rbindlist(lapply(metric_files, fread), fill = TRUE)
dataset_scope <- fread(file.path(base, "lrsim_metrics.csv"))[, .(
  dataset, n_targets_total = n_targets)]
method_scope <- unique(all_metrics[, .(task, method, evaluation_scope)])
planned <- merge(dataset_scope[, join_key := 1L], method_scope[, join_key := 1L],
  by = "join_key", allow.cartesian = TRUE)[, join_key := NULL]
all_metrics <- merge(planned, all_metrics,
  by = c("dataset", "task", "method", "evaluation_scope"), all.x = TRUE)
all_metrics[, completed_native := dataset %in% basename(dirname(metric_files))]
all_metrics[is.na(n_targets), n_targets := n_targets_total]
all_metrics[completed_native & is.na(n_candidates), `:=`(
  n_candidates = 0L, n_positive = 0L, n_null = 0L)]
all_metrics[, evaluation_status := fcase(
  !completed_native, "missing_competitor_output",
  n_candidates == 0L, "no_evaluable_candidates",
  n_positive == 0L, "no_positive_targets_in_scored_candidates",
  default = "scored")]
all_metrics[is.na(inference_status), inference_status := fifelse(
  completed_native, "complete", "missing_competitor_output")]
all_metrics[, c("n_targets_total", "completed_native") := NULL]
fwrite(all_metrics, file.path(base, "dynamic_methods.csv"))
