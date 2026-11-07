CONFIG <- list(root=normalizePath(Sys.getenv("LLR_WORK_ROOT",getwd()),winslash="/",mustWork=TRUE),threads=parallel::detectCores(),
  arms=c("full","no_coherent_clock_signal","no_sample_temporal_conditioning",
    "no_pulse_evidence","no_donor_local_conditioning","descendant_prefork",
    "positive_only_persistence","strict_common_late_window","legacy_event_baseline"),
  output="results/evaluation/final065",min_n_eff=5,
  datasets_file="experiments/config/communication13_datasets.txt",
  original_panels=c("1:4:5:13","4:4:3:12"))
Sys.setenv(OMP_NUM_THREADS=CONFIG$threads)
library(data.table)
.libPaths(c(file.path(CONFIG$root,"library"),.libPaths()))
setDTthreads(CONFIG$threads)
CONFIG$study <- Sys.getenv("LLR_ABLATION_STUDY","historical")
CONFIG$arms <- switch(CONFIG$study,historical=CONFIG$arms,
 current_key_ablation=c("full","descendant_prefork","positive_only_persistence",
  "no_pulse_evidence","legacy_matched"))
CONFIG$output <- switch(CONFIG$study,historical=CONFIG$output,
 current_key_ablation="results/evaluation/final065/key_ablation")
runs <- data.table(model=CONFIG$arms)
runs[,workspace:=switch(CONFIG$study,
 historical=fifelse(model=="full","final065",paste0("final065_ablation/",model)),
 current_key_ablation=fifelse(model=="full","final065",
  paste0("final065_ablation/",model)))]
runs[,biology:=switch(CONFIG$study,
 historical=fifelse(model=="full","biology_final065",paste0("biology_final065_ablation/",model)),
 current_key_ablation=fifelse(model=="full","biology_final065",
  paste0("biology_final065_ablation/",model)))]
CONFIG$refit_root<-Sys.getenv("LLR_ABLATION_REFIT_ROOT","")
runs[,validation_root:=file.path(CONFIG$root,"results/validation",workspace)]
switch(as.character(nzchar(CONFIG$refit_root)),"TRUE"={
 CONFIG$output<-file.path("results/evaluation",CONFIG$refit_root,"key_ablation")
 runs[,workspace:=file.path(CONFIG$refit_root,
  fifelse(model=="full","main",paste0("ablation/",model)))]
 runs[,biology:=file.path(CONFIG$refit_root,"biology",model)]
 runs[,validation_root:=file.path(CONFIG$root,"results/evaluation",
  CONFIG$refit_root,"validation",model)]
})
OUT <- file.path(CONFIG$root,CONFIG$output)
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
native <- rbindlist(setNames(lapply(file.path(CONFIG$root,"results/evaluation",runs$workspace,
  "lrsim_metrics.csv"),fread),runs$model),idcol="model",fill=TRUE)
dynamic <- rbindlist(setNames(lapply(file.path(CONFIG$root,"results/evaluation",runs$workspace,
  "dynamic_methods.csv"),fread),runs$model),idcol="model",fill=TRUE)
nom <- dynamic[task=="LR_pair_nomination" & method=="LineaLinkR_event_max" &
  evaluation_scope=="common_candidates",.(model,dataset,nomination_auprc=auprc,
    nomination_precision=precision_at_k,nomination_candidates=n_candidates,
    nomination_status=evaluation_status)]
recovery <- rbindlist(setNames(lapply(file.path(CONFIG$root,"results/evaluation",runs$workspace,
  "lrsim_target_recovery.csv"),fread),runs$model),idcol="model")
native[,c("n_targets_available","n_targets_correct_direction"):=NULL]
target_counts <- recovery[,.(n_targets_available=sum(correct_candidate_available %in% TRUE),
  n_targets_correct_direction=sum(correct_local_direction %in% TRUE)),by=.(model,dataset)]
fwrite(recovery,file.path(OUT,"ablation_target_recovery.csv"))
per_dataset <- merge(native,nom,by=c("model","dataset"),all.x=TRUE)
per_dataset <- merge(per_dataset,target_counts,by=c("model","dataset"),all.x=TRUE)
per_dataset[n_targets==0L & result_status=="available",
 c("n_targets_available","n_targets_correct_direction"):=0L]
per_dataset[,grep("^(n_recovered|recall)$",names(per_dataset),value=TRUE):=NULL]
per_dataset[,target_availability:=ifelse(n_targets>0,n_targets_available/n_targets,NA_real_)]
per_dataset[,budget_recall:=ifelse(n_targets>0,n_targets_in_budget/n_targets,NA_real_)]
per_dataset[,native_event_auprc:=auprc]
summary <- per_dataset[,.(n_datasets=.N,n_positive_datasets=sum(n_targets>0),
  n_nomination_scored=sum(n_targets>0 & is.finite(nomination_auprc)),
  n_native_AP_scored=sum(n_targets>0 & is.finite(native_event_auprc)),
  nomination_auprc=mean(nomination_auprc[n_targets>0],na.rm=TRUE),
  nomination_precision=mean(nomination_precision[n_targets>0],na.rm=TRUE),
  native_event_auprc=mean(native_event_auprc[n_targets>0],na.rm=TRUE),
  n_targets=sum(n_targets),n_targets_available=sum(n_targets_available),
  n_targets_positive=sum(n_targets_positive_score),n_targets_correct_direction=sum(n_targets_correct_direction,na.rm=TRUE),n_targets_in_budget=sum(n_targets_in_budget),
  n_null_scored=sum(n_null_scored[n_targets==0L]),
  n_null_positive=sum(null_rank_positive[n_targets==0L]),
  null_max_abs=max(null_rank_max_abs[n_targets==0L],na.rm=TRUE)),by=model]
fwrite(per_dataset,file.path(OUT,"ablation_per_dataset.csv"))
fwrite(summary,file.path(OUT,"ablation_summary.csv"))

# All-arm molecular intersection is an additional matched comparison; the
# native target/coverage denominators above remain unchanged.
datasets <- native[model=="full" & n_targets>0,dataset]
files <- merge(CJ(model=runs$model,dataset=datasets),runs,by="model")
files[,file_id:=as.character(.I)]
files[,path:=file.path(CONFIG$root,"results/evaluation",workspace,"dynamic_methods",dataset,
  "native_pair_predictions.csv")]
pairs <- rbindlist(setNames(lapply(files$path,fread),files$file_id),idcol="file_id")
pairs <- merge(pairs,files[,.(file_id,model)],by="file_id")
pairs <- pairs[method=="LineaLinkR_event_max" & in_common_resource %in% TRUE]
common <- pairs[is.finite(score),.(n_models=uniqueN(model)),by=.(dataset,ligand,receptor)][
  n_models==nrow(runs)]
matched <- merge(pairs,common[,!"n_models"],by=c("dataset","ligand","receptor"))
matched <- merge(matched,unique(native[,.(dataset,n_targets)]),by="dataset")
matched[,prediction_key:=paste(ligand,receptor,sep="|")]
setorder(matched,model,dataset,-score,prediction_key)
matched[,candidate_rank:=seq_len(.N),by=.(model,dataset)]
matched_metrics <- matched[,.(n_candidates=.N,n_positive=sum(true_positive),
  n_targets=unique(n_targets),precision_at_k=sum(true_positive[candidate_rank<=unique(n_targets)])/unique(n_targets)),
  by=.(model,dataset)]
eligible <- matched_metrics[n_positive>=2L & n_candidates-n_positive>=2L,.(model,dataset)]
matched_ap <- merge(matched,eligible,by=c("model","dataset"))[,.(
  auprc=PRROC::pr.curve(scores.class0=score[true_positive],scores.class1=score[!true_positive],
    curve=FALSE)$auc.integral),by=.(model,dataset)]
matched_metrics <- merge(matched_metrics,matched_ap,by=c("model","dataset"),all.x=TRUE)
fwrite(matched,file.path(OUT,"ablation_matched_pair_predictions.csv"))
fwrite(matched_metrics,file.path(OUT,"ablation_matched_pair_metrics.csv"))

# Readout weights are fixed to the observed common parent in every arm.
wnt_local <- rbindlist(setNames(lapply(file.path(CONFIG$root,"results/evaluation",runs$biology,
  "wnt_local_reporter_panels.csv"),fread),runs$model),idcol="model")
wnt_endpoint <- rbindlist(setNames(lapply(file.path(CONFIG$root,"results/evaluation",runs$biology,
  "wnt_history_event_metrics.csv"),fread),runs$model),idcol="model")
wnt_native_metrics <- rbindlist(setNames(lapply(file.path(CONFIG$root,"results/evaluation",runs$biology,
  "wnt_local_reporter_metrics.csv"),fread),runs$model),idcol="model")
fwrite(wnt_native_metrics,file.path(OUT,"wnt_ablation_native_local_metrics.csv"))
fwrite(wnt_endpoint,file.path(OUT,"wnt_ablation_endpoint_metrics.csv"))
identity <- unique(wnt_local[model=="full",.(panel_id,unit_id,target,canonical_child)])
outcomes <- wnt_local[model=="full" & readout_scope=="t120",
  .(panel_id,outcome_observed,relative_change,focal_change)]
outcome_check <- wnt_local[readout_scope=="t120",.(n_values=uniqueN(relative_change,na.rm=FALSE),
  n_status=uniqueN(outcome_observed,na.rm=FALSE)),by=panel_id]
stopifnot(all(outcome_check$n_values==1L),all(outcome_check$n_status==1L))
wnt_native <- rbindlist(setNames(lapply(file.path(CONFIG$root,"results/evaluation",runs$biology,
  "wnt_native_candidates.csv"),fread),runs$model),idcol="model",fill=TRUE)
pair_map_check <- wnt_native[,.(n_ligand=uniqueN(ligand_genes),n_receptor=uniqueN(receptor_genes)),by=pair_idx]
stopifnot(all(pair_map_check$n_ligand==1L),all(pair_map_check$n_receptor==1L))
two_sides <- rbindlist(setNames(lapply(file.path(CONFIG$root,"results/evaluation",runs$biology,
  "wnt_branch_candidate_predictions.csv"),fread),runs$model),idcol="model")[
  n_supported_sides==2L,.(model,panel_id,source,pair_idx)]
wnt_rows <- merge(wnt_native,identity,by=c("unit_id","target"))[
  child==canonical_child & evaluable %in% TRUE & is.finite(event_ranking_score) & is.finite(no_fork_ranking_score)]
wnt_rows <- merge(wnt_rows,two_sides,by=c("model","panel_id","source","pair_idx"))
wnt_common <- wnt_rows[,.(n_models=uniqueN(model)),by=.(panel_id,source,pair_idx)][n_models==nrow(runs)]
wnt_rows <- merge(wnt_rows,wnt_common[,!"n_models"],by=c("panel_id","source","pair_idx"))
wnt_pred <- wnt_rows[,.(n_candidates=.N,full_ranking=mean(event_ranking_score),
  matched_focal_ranking=mean(no_fork_ranking_score),full_local=mean(sibling_contrast),
  focal_local=mean(post_fork_change),lr_focal_pre=mean(focal_pre),lr_focal_post=mean(focal_post),
  lr_sibling_pre=mean(other_pre),lr_sibling_post=mean(other_post),
  lr_sibling_change=mean(sibling_change),
  lr_persistence=mean(persistent_positive_score),joint_pulse_weight=mean(joint_pulse_weight)),
  by=.(model,panel_id)]
wnt_fixed <- merge(CJ(model=CONFIG$arms,panel_id=identity$panel_id),outcomes,by="panel_id",all.x=TRUE)
wnt_fixed <- merge(wnt_fixed,wnt_pred,by=c("model","panel_id"),all.x=TRUE)
wnt_long <- melt(wnt_fixed,id.vars=c("model","panel_id","outcome_observed","relative_change"),
  measure.vars=c("full_ranking","matched_focal_ranking","full_local","focal_local"),
  variable.name="method",value.name="prediction")
wnt_metrics <- wnt_long[,.(n_planned=.N,n_observed=sum(outcome_observed %in% TRUE),
  n_predicted=sum(outcome_observed %in% TRUE & is.finite(prediction)),
  n_direction_correct=sum(outcome_observed %in% TRUE & sign(prediction)==sign(relative_change),na.rm=TRUE),
  spearman=cor(prediction[outcome_observed %in% TRUE],relative_change[outcome_observed %in% TRUE],
    method="spearman",use="complete.obs")),by=.(model,method)]
fwrite(wnt_fixed,file.path(OUT,"wnt_ablation_fixed_candidate_panels.csv"))
fwrite(wnt_fixed[panel_id %chin% CONFIG$original_panels],file.path(OUT,"wnt_original_two_errors_final.csv"))
fwrite(wnt_metrics,file.path(OUT,"wnt_ablation_fixed_candidate_metrics.csv"))

validation <- rbindlist(setNames(lapply(file.path(runs$validation_root,
  "tables/supp_muscle_depletion_validation.csv"),fread),runs$model),idcol="model")
entrain <- rbindlist(setNames(lapply(file.path(runs$validation_root,
  "tables/supp_muscle_depletion_ligand_level.csv"),fread),runs$model),idcol="model")
fwrite(validation,file.path(OUT,"muscle_ablation_fixed_candidate_metrics.csv"))
fwrite(entrain,file.path(OUT,"muscle_ablation_entrain_panel_metrics.csv"))
print(summary)
print(matched_metrics[,.(mean_AP=mean(auprc,na.rm=TRUE),mean_precision=mean(precision_at_k),n=.N),by=model])
print(wnt_metrics)
print(validation[method=="llr_full" & k==20L])

switch(CONFIG$study,current_key_ablation={
 comparison_files<-merge(CJ(model=runs[model!="full",model],
  dataset=readLines(file.path(CONFIG$root,CONFIG$datasets_file),encoding="UTF-8")),runs,by="model")
 upstream_match<-comparison_files[,{
  reference<-readRDS(file.path(CONFIG$root,"cache",runs[model=="full",workspace],dataset,"linealinkr_complete.rds"))
  candidate<-readRDS(file.path(CONFIG$root,"cache",workspace,dataset,"linealinkr_complete.rds"))
  .(n_reference_cells=reference$n_cells,n_candidate_cells=candidate$n_cells,
   component_ablation=as.character(candidate$config_snapshot$component_ablation)[1L],
   activity_transform=as.character(candidate$config_snapshot$activity_transform)[1L],
   context_error_weighting=as.logical(candidate$config_snapshot$context_error_weighting)[1L],
   late_bin_merging_enabled=as.logical(candidate$config_snapshot$late_bin_merging_enabled)[1L],
   roots_equal=identical(reference$start_cell,candidate$start_cell),
   paths_equal=identical(reference$lineage_paths,candidate$lineage_paths),
   cell_ids_equal=identical(reference$cells,candidate$cells),
   pseudotime_max_difference=max(abs(reference$pseudotime-candidate$pseudotime)),
   lineage_weight_max_difference=max(abs(reference$lineage_weights-candidate$lineage_weights)),
   molecular_candidates_equal=isTRUE(all.equal(as.data.frame(reference$pair_gene_map),
    as.data.frame(candidate$pair_gene_map),check.attributes=FALSE,tolerance=0)))
 },by=.(model,dataset,workspace)]
 fwrite(upstream_match,file.path(OUT,"matched_upstream_all13.csv"))
 corrected_controls<-upstream_match[model %chin% c("legacy_matched","positive_only_persistence")]
 stopifnot(nrow(corrected_controls)==26L,
  all(upstream_match$activity_transform=="log_exposure"),
  all(corrected_controls$component_ablation==fifelse(corrected_controls$model=="legacy_matched",
    "legacy_log_error_weighting","positive_only_persistence")),
  all(corrected_controls$context_error_weighting),all(corrected_controls$late_bin_merging_enabled))
 stopifnot(all(upstream_match$roots_equal),all(upstream_match$paths_equal),
  all(upstream_match$cell_ids_equal),all(upstream_match$molecular_candidates_equal),
  all(upstream_match$n_reference_cells==upstream_match$n_candidate_cells),
  max(upstream_match$pseudotime_max_difference)<1e-10,
  max(upstream_match$lineage_weight_max_difference)<1e-10)
 print(upstream_match[,.(n_datasets=.N,n_cells_equal=sum(n_reference_cells==n_candidate_cells),
  max_clock_difference=max(pseudotime_max_difference),
  max_weights_difference=max(lineage_weight_max_difference)),by=model])
})
