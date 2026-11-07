CONFIG <- list(root=normalizePath(Sys.getenv("LLR_WORK_ROOT",getwd()),winslash="/",mustWork=TRUE),threads=parallel::detectCores(),
  workspace=Sys.getenv("LLR_BIO_WS","cache/final065"),
  reference_workspace=Sys.getenv("LLR_BIO_REFERENCE_WS",Sys.getenv("LLR_ANCHOR_WORKSPACE","")),
  output=Sys.getenv("LLR_BIO_OUT","results/evaluation/biology_final065"),
  dataset="GSE274389_signal_recorder",min_n_eff=5,
  reporter="results/evaluation/biology_main065_final/wnt_reporter_cells.csv",
  readout_weighting="shared")
Sys.setenv(OMP_NUM_THREADS=CONFIG$threads,OPENBLAS_NUM_THREADS=CONFIG$threads,MKL_NUM_THREADS=CONFIG$threads)
library(data.table)
library(Matrix)
setDTthreads(CONFIG$threads)
OUT <- file.path(CONFIG$root,CONFIG$output)
x <- readRDS(file.path(CONFIG$root,CONFIG$workspace,CONFIG$dataset,"linealinkr_complete.rds"))
reference <- switch(as.character(nzchar(CONFIG$reference_workspace)),
 "TRUE"=readRDS(file.path(CONFIG$root,"cache",CONFIG$reference_workspace,CONFIG$dataset,"linealinkr_complete.rds")),
 "FALSE"=x)
stopifnot(identical(reference$cells,x$cells),
 identical(reference$lineage_paths,x$lineage_paths),
 identical(reference$celltypes,x$celltypes),
 !is.null(reference$lineage_reference_weights) ||
  is.null(reference$config_snapshot$lineage_weight_mode) ||
  reference$config_snapshot$lineage_weight_mode=="reference_geometry")
reference_weights <- switch(as.character(is.null(reference$lineage_reference_weights)),
 "TRUE"=reference$lineage_weights,"FALSE"=reference$lineage_reference_weights)
meta <- fread(file.path(x$input_extract_dir,paste0(CONFIG$dataset,"_meta.csv")),encoding="UTF-8")
reporter <- fread(file.path(CONFIG$root,CONFIG$reporter),encoding="UTF-8")
stopifnot(identical(meta$barcode_id,x$cells),identical(reporter$barcode_id,x$cells))
cells <- data.table(cell_index=seq_along(x$cells),node=reference$lineage_clusters,
  target=x$celltypes,time_hours=meta$time_hours,reporter=reporter$irfp_log1p_cptt)
units <- data.table(unit_id=seq_along(attr(reference$fork_contrasts,"event_units")$units),
  unit=attr(reference$fork_contrasts,"event_units")$units)
desc <- units[,.(child=rep(unit[[1L]]$children,lengths(unit[[1L]]$desc)),
  lineage=unlist(unit[[1L]]$desc,use.names=FALSE)),by=unit_id]
children <- unique(desc[,.(unit_id,child)])
children[,child_index:=.I]
desc <- merge(desc,children,by=c("unit_id","child"))
membership <- sparseMatrix(i=desc$lineage,j=desc$child_index,x=1,
  dims=c(ncol(x$lineage_weights),nrow(children)))
child_weights <- as.matrix(reference_weights %*% membership) / rowSums(reference_weights)
pre_nodes <- merge(children,units[,.(node=unit[[1L]]$segment_pre_nodes),by=unit_id],
  by="unit_id",allow.cartesian=TRUE)
pre_nodes[,phase:="pre"]
post_nodes <- units[,.(child=rep(unit[[1L]]$children,lengths(unit[[1L]]$segment_post_nodes)),
  node=unlist(unit[[1L]]$segment_post_nodes,use.names=FALSE)),by=unit_id]
post_nodes <- merge(post_nodes,children,by=c("unit_id","child"))
post_nodes[,phase:="post"]
groups <- unique(rbindlist(list(pre_nodes,post_nodes),use.names=TRUE))
members <- merge(cells,groups,by="node",allow.cartesian=TRUE)
members[,weight:=child_weights[cbind(cell_index,child_index)]]
switch(as.character(identical(CONFIG$readout_weighting,"shared")),"TRUE"={
  parent <- members[phase=="pre",.(parent_mass=sum(weight),n_children=uniqueN(child)),by=.(unit_id,cell_index)]
  members <- merge(members,parent,by=c("unit_id","cell_index"),all.x=TRUE)
  members[phase=="pre",weight:=parent_mass/n_children]
})
members <- rbindlist(list(copy(members)[,readout_scope:="all_native_cells"],
  copy(members)[time_hours==120][,readout_scope:="t120"]))
means <- members[,.(mass=sum(weight),sum_sq=sum(weight^2),
  reporter_mean=sum(weight*reporter)/sum(weight),n_cells=.N,
  t120_fraction=sum(weight*(time_hours==120))/sum(weight)),
  by=.(readout_scope,unit_id,child,target,phase)]
means[,n_eff:=mass^2/sum_sq]
wide <- dcast(means,readout_scope+unit_id+child+target~phase,
  value.var=c("reporter_mean","n_eff","n_cells","t120_fraction"))
wide[,reporter_change:=reporter_mean_post-reporter_mean_pre]
panels <- fread(file.path(OUT,"wnt_history_event_panels.csv"))
identity <- panels[,.(panel_id,unit_id,target,canonical_child,canonical_sibling,target_type)]
local_focal <- wide[,.(readout_scope,unit_id,target,canonical_child=child,
  focal_pre=reporter_mean_pre,focal_post=reporter_mean_post,
  focal_change=reporter_change,focal_n_pre=n_eff_pre,focal_n_post=n_eff_post,
  focal_t120_pre=t120_fraction_pre,focal_t120_post=t120_fraction_post)]
local_sister <- wide[,.(readout_scope,unit_id,target,canonical_sibling=child,
  sibling_pre=reporter_mean_pre,sibling_post=reporter_mean_post,
  sibling_change=reporter_change,sibling_n_pre=n_eff_pre,sibling_n_post=n_eff_post,
  sibling_t120_pre=t120_fraction_pre,sibling_t120_post=t120_fraction_post)]
local <- merge(CJ(panel_id=identity$panel_id,readout_scope=c("all_native_cells","t120")),
  identity,by="panel_id")
local <- merge(local,local_focal,by=c("readout_scope","unit_id","target","canonical_child"),all.x=TRUE)
local <- merge(local,local_sister,by=c("readout_scope","unit_id","target","canonical_sibling"),all.x=TRUE)
local[,c("outcome_observed","relative_change"):=.(
  focal_n_pre>=CONFIG$min_n_eff & focal_n_post>=CONFIG$min_n_eff &
  sibling_n_pre>=CONFIG$min_n_eff & sibling_n_post>=CONFIG$min_n_eff,
  focal_change-sibling_change)]
pmap <- as.data.table(x$pair_gene_map)
wnt <- pmap[pathway_name=="WNT" & startsWith(toupper(ligand_genes),"WNT"),
  .(pair_idx,ligand_genes,receptor_genes)]
candidates <- merge(as.data.table(x$fork_contrasts),wnt,by="pair_idx")
candidates <- merge(candidates,identity[,.(panel_id,unit_id,target,canonical_child)],
  by=c("unit_id","target"))[child==canonical_child]
common <- fread(file.path(OUT,"wnt_branch_candidate_predictions.csv"))[
  n_supported_sides==2L,.(panel_id,source,pair_idx)]
candidates <- merge(candidates,common,by=c("panel_id","source","pair_idx"))
candidates <- candidates[evaluable %in% TRUE & is.finite(event_ranking_score) &
  is.finite(no_fork_ranking_score)]
predictions <- candidates[,.(n_candidates=.N,n_lr_pairs=uniqueN(pair_idx),
  full_local=mean(sibling_contrast),focal_local=mean(post_fork_change),
  full_ranking=mean(event_ranking_score),matched_focal_ranking=mean(no_fork_ranking_score),
  post_activity=mean(post_activity)),by=panel_id]
local <- merge(local,predictions,by="panel_id",all.x=TRUE)
pred <- melt(local,id.vars=c("panel_id","readout_scope","outcome_observed"),
  measure.vars=c("full_local","focal_local","full_ranking","matched_focal_ranking","post_activity"),
  variable.name="method",value.name="prediction")
truth <- melt(local,id.vars=c("panel_id","readout_scope","outcome_observed"),
  measure.vars=c("relative_change","focal_change"),variable.name="readout",value.name="outcome")
joined <- merge(pred,truth,by=c("panel_id","readout_scope","outcome_observed"),allow.cartesian=TRUE)
metrics <- joined[,.(n_panels=.N,
  n_outcome_available=sum(outcome_observed %in% TRUE & is.finite(outcome)),
  n_predicted=sum(outcome_observed %in% TRUE & is.finite(outcome) & is.finite(prediction)),
  n_direction_correct=sum(outcome_observed %in% TRUE & sign(outcome)==sign(prediction),na.rm=TRUE),
  spearman=ifelse(sum(outcome_observed %in% TRUE & is.finite(outcome) & is.finite(prediction))>=2L,
    cor(prediction[outcome_observed %in% TRUE],outcome[outcome_observed %in% TRUE],
      method="spearman",use="complete.obs"),NA_real_)),
  by=.(readout_scope,readout,method)]
fwrite(means,file.path(OUT,"wnt_local_reporter_groups.csv"))
fwrite(local,file.path(OUT,"wnt_local_reporter_panels.csv"))
fwrite(metrics,file.path(OUT,"wnt_local_reporter_metrics.csv"))
print(metrics[readout_scope=="t120"])
print(local[readout_scope=="t120" & outcome_observed %in% TRUE,
  .(panel_id,focal_change,sibling_change,relative_change,
    full_local,focal_local,full_ranking,matched_focal_ranking)])
