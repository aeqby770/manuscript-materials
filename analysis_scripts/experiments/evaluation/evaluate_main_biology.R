CONFIG <- list(
  root = normalizePath(Sys.getenv("LLR_WORK_ROOT",getwd()),winslash="/",mustWork=TRUE),
  workspace = Sys.getenv("LLR_BIO_WS", "cache/final065"),
  output = Sys.getenv("LLR_BIO_OUT", "results/evaluation/biology_final065"),
  reporter_input = "results/evaluation/biology_main065_final/wnt_reporter_cells.csv",
  wnt = "GSE274389_signal_recorder",
  muscle = "GSE143437_muscle_regeneration",
  min_pct = 0.05, min_n_eff = 5,
  source = "FAPs", targets = c("Monocytes/Macrophages/Platelets",
    "Pro-inflammatory macrophages", "Resident Macrophages/APCs"),
  fixed_ligand = "Ccl2", fixed_receptor = "Ccr2",
  threads = parallel::detectCores()
)
Sys.setenv(OMP_NUM_THREADS = CONFIG$threads)
library(data.table)
library(Matrix)
setDTthreads(CONFIG$threads)
OUT <- file.path(CONFIG$root, CONFIG$output)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
switch(as.character(!"--muscle-only"%in%commandArgs(TRUE)),"TRUE"={
x <- readRDS(file.path(CONFIG$root, CONFIG$workspace, CONFIG$wnt, "linealinkr_complete.rds"))
meta <- fread(file.path(x$input_extract_dir, paste0(CONFIG$wnt, "_meta.csv")), encoding = "UTF-8")
reporter <- fread(file.path(CONFIG$root, CONFIG$reporter_input), encoding = "UTF-8")
stopifnot(identical(meta$barcode_id, x$cells), identical(reporter$barcode_id, x$cells))
types <- sort(unique(as.character(meta[[x$celltype_column]])))
celltypes <- match(as.character(meta[[x$celltype_column]]), types)
t96 <- meta$time_hours == 96
t120 <- meta$time_hours == 120
positive <- meta$wnt_recorded == "positive"
stopifnot(all(meta$wnt_recorded[t120] %in% c("positive", "negative")))
pmap <- as.data.table(x$pair_gene_map)
wnt_pairs <- pmap[pathway_name == "WNT" & startsWith(toupper(ligand_genes), "WNT")]
lw <- as.matrix(x$lineage_weights)

# Native temporal readout over the full resource and source-type universe.
temporal <- as.data.table(x$temporal_lr$scores)[pair_idx %in% wnt_pairs$pair_idx]
temporal[, fixed_score := fifelse(is.finite(score), score, 0)]
grid <- temporal[, .(native_temporal = sum(fixed_score) /
  (nrow(wnt_pairs) * length(types))), by = .(target, lineage, bin)]
bounds <- attr(x$temporal_lr$scores, "bin_bounds")[[1L]]
cell_bin <- findInterval(x$pseudotime, bounds, all.inside = TRUE, rightmost.closed = TRUE)
cell_lineage <- data.table(
  cell_index = rep(seq_along(x$cells), times = ncol(lw)),
  target = rep(celltypes, times = ncol(lw)),
  lineage = rep(seq_len(ncol(lw)), each = nrow(lw)),
  bin = rep(cell_bin, times = ncol(lw)), weight = as.vector(lw))
cell_lineage <- merge(cell_lineage, grid, by = c("target", "lineage", "bin"),
  all.x = TRUE, sort = FALSE)
cell_lineage[is.na(native_temporal), native_temporal := 0]
cell_scores <- cell_lineage[, .(native_temporal = sum(weight * native_temporal)), by = cell_index]
setorder(cell_scores, cell_index)

# Fork-masked aggregation (mechanism probe C1 variant V2): a cell only draws a
# lineage's binned score once that lineage has diverged (bin >= fork_bin).
# Uncommitted cells sitting on the shared trunk stop inheriting sibling-lineage
# early-bin signal; cells with no eligible lineage fall back to the
# dominant-lineage score. No outcome information is used.
fork_bin_l <- as.data.table(x$fork_bins)[order(lineage)]$fork_bin
maxw <- lw[cbind(seq_len(nrow(lw)), max.col(lw, ties.method = "first"))]
cell_lineage[, `:=`(lineage_eligible = bin >= fork_bin_l[lineage],
  is_max_weight = weight == maxw[cell_index])]
masked_scores <- cell_lineage[lineage_eligible == TRUE,
  .(masked = sum(weight * native_temporal) / sum(weight)), by = cell_index]
dominant_scores <- cell_lineage[is_max_weight == TRUE,
  .(dominant = mean(native_temporal)), by = cell_index]
cell_scores[masked_scores, native_temporal_fork_masked := i.masked, on = "cell_index"]
cell_scores[dominant_scores, .dominant_tmp := i.dominant, on = "cell_index"]
cell_scores[is.na(native_temporal_fork_masked),
  native_temporal_fork_masked := .dominant_tmp]
cell_scores[, .dominant_tmp := NULL]

# Formula-matched no-time/no-branch control; individual receiver expression.
genes <- readLines(file.path(x$input_extract_dir, paste0(CONFIG$wnt, "_genes.txt")), encoding = "UTF-8")
ligand_map <- wnt_pairs[, .(gene = unlist(strsplit(ligand_genes, "_", fixed = TRUE))), by = pair_idx]
receptor_map <- wnt_pairs[, .(gene = unlist(strsplit(receptor_genes, "_", fixed = TRUE))), by = pair_idx]
selected_genes <- unique(c(ligand_map$gene, receptor_map$gene))
counts <- as(readMM(file.path(x$input_extract_dir, paste0(CONFIG$wnt, "_matrix.mtx"))), "CsparseMatrix")
expression <- counts[match(selected_genes, genes), , drop = FALSE] %*%
  Diagonal(x = x$config_snapshot$lr_scale_factor / Matrix::colSums(counts))
expression@x <- log1p(expression@x)
rm(counts)
ligand_map[, gene_index := match(gene, selected_genes)]
receptor_map[, gene_index := match(gene, selected_genes)]
early_membership <- sparseMatrix(i = seq_len(sum(t96)), j = celltypes[t96], x = 1,
  dims = c(sum(t96), length(types)))
type_size <- Matrix::colSums(early_membership)
type_mean <- as.matrix(expression[, t96, drop = FALSE] %*% early_membership) /
  rep(type_size, each = length(selected_genes))
type_pct <- as.matrix((expression[, t96, drop = FALSE] > 0) %*% early_membership) /
  rep(type_size, each = length(selected_genes))
ligand_stats <- data.table(
  pair_idx = rep(ligand_map$pair_idx, times = length(types)),
  source = rep(seq_along(types), each = nrow(ligand_map)),
  value = as.vector(type_mean[ligand_map$gene_index, , drop = FALSE]),
  pct = as.vector(type_pct[ligand_map$gene_index, , drop = FALSE]))[
    , .(ligand_mean = min(value), ligand_pct = min(pct)), by = .(pair_idx, source)]
receptor_stats <- data.table(
  pair_idx = rep(receptor_map$pair_idx, times = sum(t96)),
  early_index = rep(seq_len(sum(t96)), each = nrow(receptor_map)),
  receptor_cell = as.vector(as.matrix(expression[receptor_map$gene_index, t96, drop = FALSE])),
  receptor_pct = as.vector(type_pct[receptor_map$gene_index, celltypes[t96], drop = FALSE]))[
    , .(receptor_cell = min(receptor_cell), receptor_pct = min(receptor_pct)),
    by = .(pair_idx, early_index)]
control <- merge(ligand_stats, receptor_stats, by = "pair_idx", allow.cartesian = TRUE)
control[, geo := sqrt(ligand_mean * receptor_cell)]
control[, score := fifelse(type_size[source] >= CONFIG$min_n_eff &
  type_size[celltypes[which(t96)][early_index]] >= CONFIG$min_n_eff &
  ligand_pct >= CONFIG$min_pct & receptor_pct >= CONFIG$min_pct,
  geo / (attr(x$temporal_lr$scores, "hill_scale") + geo), 0)]
control_cells <- control[, .(same_hill_no_time_no_branch = sum(score) /
  (nrow(wnt_pairs) * length(types))), by = early_index]
setorder(control_cells, early_index)
early <- data.table(barcode_id = x$cells[t96], target_type = types[celltypes[t96]],
  native_temporal = cell_scores$native_temporal[t96],
  same_hill_no_time_no_branch = control_cells$same_hill_no_time_no_branch,
  reporter = reporter$irfp_log1p_cptt[t96],
  native_temporal_fork_masked = cell_scores$native_temporal_fork_masked[t96])
stopifnot(nrow(early) == sum(t96), all(is.finite(as.matrix(early[, 3:6]))))
early_long <- melt(early, id.vars = c("barcode_id", "target_type", "reporter"),
  variable.name = "method", value.name = "prediction")
centered <- copy(early_long)
centered[, `:=`(prediction = prediction - mean(prediction),
  reporter = reporter - mean(reporter)), by = .(method, target_type)]
centered[, target_type := "within_type_pooled"]
early_long <- rbindlist(list(early_long, copy(early_long)[, target_type := "pooled"], centered))
early_metrics <- early_long[, .(n_cells = .N, n_distinct_prediction = uniqueN(prediction),
  spearman = cor(prediction, reporter, method = "spearman")), by = .(method, target_type)]
fwrite(early, file.path(OUT, "wnt_t96_cells.csv"))
fwrite(early_metrics, file.path(OUT, "wnt_t96_metrics.csv"))

# Transport t96 exposure to every t120 cell on the same native fate geometry.
early_w <- lw[t96, , drop = FALSE]
mass <- colSums(early_w)
effective <- mass^2 / colSums(early_w^2)
supported <- is.finite(effective) & effective >= CONFIG$min_n_eff
exposure <- crossprod(early_w, as.matrix(early[, 3:6])) / mass
exposure[!supported, ] <- 0
prediction <- lw[t120, , drop = FALSE] %*% exposure
colnames(prediction)[3L] <- "measured_t96_reporter_geometry_control"
late <- data.table(barcode_id = x$cells[t120], target_type = types[celltypes[t120]],
  recorded_positive = positive[t120], supported_mass = rowSums(lw[t120, supported, drop = FALSE]),
  as.data.table(prediction))
late_long <- melt(late, id.vars = c("barcode_id", "target_type", "recorded_positive", "supported_mass"),
  variable.name = "method", value.name = "prediction")
late_long <- rbindlist(list(late_long, copy(late_long)[, target_type := "pooled"]))
late_long[, prediction_rank := frank(prediction, ties.method = "average"), by = .(method, target_type)]
late_metrics <- late_long[, .(n_cells = .N, n_positive = sum(recorded_positive),
  supported_mass_mean = mean(supported_mass),
  auc = (sum(prediction_rank[recorded_positive]) - sum(recorded_positive) *
    (sum(recorded_positive) + 1) / 2) / (sum(recorded_positive) * sum(!recorded_positive))),
  by = .(method, target_type)]
fwrite(late, file.path(OUT, "wnt_t120_cells.csv"))
fwrite(late_metrics, file.path(OUT, "wnt_t120_metrics.csv"))

# Direct branch outcomes use observed descendant states, not future-fate mass
# assigned to cells still on a shared trunk.
fc <- as.data.table(x$fork_contrasts)
unit_list <- attr(fc,"event_units")$units
units <- data.table(unit_id=seq_along(unit_list),unit=unit_list)
descendants <- units[,.(event_node=unit[[1L]]$event_node,
  child=rep(unit[[1L]]$children,lengths(unit[[1L]]$desc)),
  lineage=unlist(unit[[1L]]$desc,use.names=FALSE)),by=unit_id]
children <- unique(descendants[,.(unit_id,event_node,child)])
children[,child_index:=.I]
descendants <- merge(descendants,children,by=c("unit_id","event_node","child"))
path_nodes <- data.table(lineage=rep(seq_along(x$lineage_paths),lengths(x$lineage_paths)),
  node=unlist(x$lineage_paths),position=sequence(lengths(x$lineage_paths)))
child_positions <- merge(descendants,
  path_nodes[,.(lineage,child=node,child_position=position)],by=c("lineage","child"))
observed_nodes <- unique(merge(child_positions,path_nodes,by="lineage",
  allow.cartesian=TRUE)[position>=child_position,.(child_index,node)])
observed_cells <- merge(data.table(cell_index=seq_along(x$cells),
  node=x$lineage_clusters),observed_nodes,by="node",allow.cartesian=TRUE)
child_w <- as.matrix(sparseMatrix(i=observed_cells$cell_index,
  j=observed_cells$child_index,x=1,dims=c(length(x$cells),nrow(children))))
status_membership <- sparseMatrix(i=seq_len(sum(t120)),
  j=(celltypes[t120]-1L)*2L+as.integer(positive[t120])+1L,x=1,
  dims=c(sum(t120),length(types)*2L))
endpoint_mass <- as.matrix(crossprod(child_w[t120,,drop=FALSE],status_membership))
endpoint_grid <- CJ(child_index=seq_len(nrow(children)),target=seq_along(types),status=1:2)
endpoint_grid[,weight_mass:=endpoint_mass[cbind(child_index,(target-1L)*2L+status)]]
endpoint <- dcast(endpoint_grid,child_index+target~status,value.var="weight_mass")
setnames(endpoint,c("1","2"),c("negative_mass","positive_mass"))
endpoint[,c("mass","n_eff"):=
  .(positive_mass+negative_mass,positive_mass+negative_mass)]
type_membership <- sparseMatrix(i=seq_len(sum(t120)),j=celltypes[t120],x=1,
  dims=c(sum(t120),length(types)))
reporter_sum <- as.matrix(crossprod(child_w[t120,,drop=FALSE],
  type_membership*reporter$irfp_log1p_cptt[t120]))
endpoint[,reporter_now:=reporter_sum[cbind(child_index,target)]/mass]
endpoint[,c("recorder_fraction","recorder_fraction_raw"):=.(
  (positive_mass/sum(positive[t120]))/
    (positive_mass/sum(positive[t120])+negative_mass/sum(!positive[t120])),
  positive_mass/mass)]
endpoint <- merge(endpoint,children,by="child_index")
focal <- endpoint[,.(unit_id,event_node,target,canonical_child=child,
  focal_fraction=recorder_fraction,focal_fraction_raw=recorder_fraction_raw,
  focal_n_eff=n_eff,focal_reporter_now=reporter_now)]
sister <- endpoint[,.(unit_id,event_node,target,canonical_sibling=child,
  sister_fraction=recorder_fraction,sister_fraction_raw=recorder_fraction_raw,
  sister_n_eff=n_eff,sister_reporter_now=reporter_now)]
panels <- merge(focal,sister,by=c("unit_id","event_node","target"),
  allow.cartesian=TRUE)[canonical_child<canonical_sibling]
panels[,c("target_type","outcome_observed","recorder_contrast",
  "recorder_contrast_raw","reporter_now_contrast","panel_id"):=.(
  types[target],focal_n_eff>=CONFIG$min_n_eff & sister_n_eff>=CONFIG$min_n_eff,
  focal_fraction-sister_fraction,focal_fraction_raw-sister_fraction_raw,
  focal_reporter_now-sister_reporter_now,
  paste(unit_id,target,canonical_child,canonical_sibling,sep=":"))]
wnt_events <- merge(fc,wnt_pairs[,.(pair_idx,ligand_genes,receptor_genes)],by="pair_idx")
event_predictions <- merge(wnt_events,panels[,.(unit_id,target,panel_id,
  canonical_child,canonical_sibling)],by=c("unit_id","target"),allow.cartesian=TRUE)
event_predictions <- event_predictions[child %in% c(canonical_child,canonical_sibling)]
event_predictions[,orientation:=fifelse(child==canonical_child,1,-1)]
source_pair <- event_predictions[evaluable %in% TRUE & is.finite(sibling_contrast) &
  is.finite(event_ranking_score) & is.finite(no_fork_ranking_score),.(
  n_supported_sides=uniqueN(child),
  native_contrast=mean(sibling_contrast*orientation),
  main_ranking_contrast=mean(event_ranking_score*orientation),
  matched_focal_ranking_contrast=mean(no_fork_ranking_score*orientation),
  focal_change_contrast=mean(post_fork_change*orientation),
  post_activity_contrast=mean(post_activity*orientation)),
  by=.(panel_id,source,pair_idx,ligand_genes,receptor_genes)]
own_scores <- source_pair[,.(n_native_source_pairs=.N,n_native_lr_pairs=uniqueN(pair_idx),
  native_contrast=mean(native_contrast)),by=panel_id]
# A scalar focal score needs both observed sides to form a paired contrast.
# The same two-side molecular contexts are used in every common-candidate arm.
common_pair <- source_pair[n_supported_sides==2L]
event_scores <- common_pair[,.(n_source_pairs=.N,n_lr_pairs=uniqueN(pair_idx),
  paired_native_contrast=mean(native_contrast),
  main_ranking_contrast=mean(main_ranking_contrast),
  matched_focal_ranking_contrast=mean(matched_focal_ranking_contrast),
  focal_change_contrast=mean(focal_change_contrast),
  post_activity_contrast=mean(post_activity_contrast)),by=panel_id]
panels <- merge(panels,own_scores,by="panel_id",all.x=TRUE)
panels <- merge(panels,event_scores,by="panel_id",all.x=TRUE)
source_pair <- merge(source_pair,panels[,.(panel_id,target_type,outcome_observed,
  focal_n_eff,sister_n_eff,reporter_now_contrast,recorder_contrast)],by="panel_id")
fwrite(source_pair,file.path(OUT,"wnt_branch_candidate_predictions.csv"))
methods <- c("native_contrast","paired_native_contrast","main_ranking_contrast",
  "matched_focal_ranking_contrast",
  "focal_change_contrast","post_activity_contrast")
panel_long <- melt(panels,id.vars=c("panel_id","outcome_observed"),
  measure.vars=methods,variable.name="method",value.name="prediction")
readouts <- melt(panels,id.vars=c("panel_id","outcome_observed"),
  measure.vars=c("reporter_now_contrast","recorder_contrast"),
  variable.name="readout",value.name="outcome")
panel_long <- merge(panel_long,readouts,by=c("panel_id","outcome_observed"),
  allow.cartesian=TRUE)
panel_long[,direction_correct:=sign(prediction)==sign(outcome)]
panel_metrics <- panel_long[,.(n_panels=.N,
  n_outcome_available=sum(outcome_observed & is.finite(outcome)),
  n_predicted=sum(outcome_observed & is.finite(prediction) & is.finite(outcome)),
  n_direction_correct=sum(outcome_observed & direction_correct,na.rm=TRUE),
  spearman=cor(prediction[outcome_observed],outcome[outcome_observed],
    method="spearman",use="complete.obs")),by=.(readout,method)]
fwrite(panels,file.path(OUT,"wnt_history_event_panels.csv"))
fwrite(panel_metrics,file.path(OUT,"wnt_history_event_metrics.csv"))
fwrite(wnt_events,file.path(OUT,"wnt_native_candidates.csv"))
print(early_metrics[target_type %in% c("pooled","within_type_pooled")])
print(late_metrics[target_type=="pooled"])
print(panel_metrics)

# Fixed FAP axis on the current native result; all receiver contexts stay visible.

})

y <- readRDS(file.path(CONFIG$root, CONFIG$workspace, CONFIG$muscle, "linealinkr_complete.rds"))
muscle_meta <- fread(file.path(y$input_extract_dir, paste0(CONFIG$muscle, "_meta.csv")), encoding = "UTF-8")
stopifnot(identical(muscle_meta$barcode_id, y$cells))
muscle_types <- sort(unique(as.character(muscle_meta[[y$celltype_column]])))
events <- merge(as.data.table(y$fork_contrasts), as.data.table(y$pair_gene_map), by = "pair_idx")
events[, `:=`(source_label = muscle_types[source], target_label = muscle_types[target])]
fap <- events[source_label == CONFIG$source]
fap[evaluable & is.finite(event_ranking_score),
  resource_rank := frank(-event_ranking_score, ties.method = "min"),
  by = .(target, event_node, child)]
fap_coverage <- fap[, .(n_candidate_rows = .N,
  n_evaluable = sum(evaluable & is.finite(sibling_contrast)),
  n_positive_main_score = sum(evaluable & event_ranking_score > 0, na.rm = TRUE)),
  by = target_label]
fixed <- fap[ligand_genes == CONFIG$fixed_ligand & receptor_genes == CONFIG$fixed_receptor]
pooled <- merge(as.data.table(y$lr_scores), as.data.table(y$pair_gene_map), by = "pair_idx")
pooled[, `:=`(source_label = muscle_types[source], target_label = muscle_types[target])]
pooled <- pooled[source_label == CONFIG$source]
pooled[, resource_rank := frank(-score, ties.method = "min"), by = target]
fixed_pooled <- pooled[ligand_genes == CONFIG$fixed_ligand & receptor_genes == CONFIG$fixed_receptor]
clock <- data.table(barcode_id = y$cells, sample_id = muscle_meta$sample_id,
  time_days = muscle_meta$time_days, cell_type = muscle_meta[[y$celltype_column]],
  pseudotime = y$pseudotime, node = y$lineage_clusters)
clock_summary <- clock[, .(n_cells = .N, median_time = median(pseudotime),
  time_rank_correlation = cor(pseudotime, time_days, method = "spearman")), by = cell_type]
fwrite(fap, file.path(OUT, "muscle_native_fap_candidates.csv"))
fwrite(fap_coverage, file.path(OUT, "muscle_native_fap_coverage.csv"))
fwrite(fixed, file.path(OUT, "muscle_fixed_ccl2_ccr2_events.csv"))
paired_rows <- as.data.table(y$fork_sample_contrasts)
pair_keys <- c("pair_idx","source","target","event_node","child","unit_id")
fwrite(merge(paired_rows,unique(fixed[,..pair_keys]),by=pair_keys),file.path(OUT,"muscle_ccl2_paired_samples.csv"))
fwrite(fixed_pooled, file.path(OUT, "muscle_fixed_ccl2_ccr2_pooled.csv"))
fwrite(clock, file.path(OUT, "muscle_clock_cells.csv"))
fwrite(clock_summary, file.path(OUT, "muscle_clock_metrics.csv"))
print(fap_coverage[target_label %in% CONFIG$targets])
print(fixed[target_label %in% CONFIG$targets, .(target_label, event_node, child,
  sibling_contrast, event_ranking_score, resource_rank, evaluable, reason)])
print(fixed_pooled[target_label %in% CONFIG$targets, .(target_label, score, resource_rank)])
print(clock_summary[cell_type %in% c(CONFIG$source, CONFIG$targets)])

switch(as.character(!"--muscle-only"%in%commandArgs(TRUE)),"TRUE"={
source(file.path(CONFIG$root,"experiments/evaluation/evaluate_wnt_local_reporter.R"),
  local=new.env(parent=globalenv()),encoding="UTF-8")

source(file.path(CONFIG$root,"experiments/evaluation/evaluate_muscle_event_calendar.R"),
  local=new.env(parent=globalenv()),encoding="UTF-8")

})
