# Reproduce Table 1 and the focal-contact peak readout from distributed records.
# Usage: Rscript --vanilla reproduce_saved_results.R evaluation SourceData.zip
CONFIG <- list(
  paths=list(input=normalizePath(commandArgs(trailingOnly=TRUE)[1], winslash="/", mustWork=TRUE),
    source_data=normalizePath(commandArgs(trailingOnly=TRUE)[2], winslash="/", mustWork=TRUE)),
  threads=parallel::detectCores(),
  pair_methods=c("LineaLinkR_event_max","legacy_matched","no_pulse_evidence","geometric_layout",
    "LineaLinkR_focal_only","LineaLinkR_matched_no_fork",
    "TraSig_native_max","scComm_no_TF_activity","pooled_LR"),
  pair_labels=c("LineaLinkR","Matched score control","Persistence only","Geometric layout (G)",
    "Raw focal","Matched focal","Native TraSig","scComm (no TF)","Pooled LR"),
  event_methods=c("LineaLinkR_event_score","TraSig_projected_contrast","scComm_projected_native_kernel"),
  event_labels=c("LineaLinkR","TraSig projection","scComm projection"),
  tolerance=1e-12
)
Sys.setenv(OMP_NUM_THREADS=CONFIG$threads)
# R uses the library paths supplied by the reader for the declared dependencies.
library(data.table)
setDTthreads(CONFIG$threads)
OUT <- file.path(CONFIG$paths$input, "reproduced")
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

# Scores and truth labels are individual candidate records, not metric summaries.
score_files <- list.files(file.path(CONFIG$paths$input,"ranked_predictions"),
  pattern="^ranked_predictions[.]csv$", recursive=TRUE, full.names=TRUE)
scores <- rbindlist(lapply(score_files, fread))
common_fields <- c("dataset","task","method","prediction_key","truth_row","score")
component_pairs <- fread(file.path(CONFIG$paths$input,"component_tables","ablation_matched_pair_predictions.csv"))[
  model %chin% c("legacy_matched","no_pulse_evidence"),
  .(dataset,task="LR_pair_nomination",method=model,prediction_key,truth_row=true_positive,score)]
geometry_files <- list.files(file.path(CONFIG$paths$input,"component_tables","geometric_layout"),
  pattern="^ranked_predictions[.]csv$",recursive=TRUE,full.names=TRUE)
geometry_pairs <- rbindlist(lapply(geometry_files,fread))[
  evaluation_scope=="common_candidates" & task=="LR_pair_nomination" & method=="LineaLinkR_event_max",
  .(dataset,task,method="geometric_layout",prediction_key,truth_row,score)]
scores <- rbindlist(list(scores[evaluation_scope=="common_candidates",..common_fields],
  component_pairs,geometry_pairs))
stopifnot(all(is.finite(scores$score)),
  uniqueN(scores$dataset)==11L,
  !anyDuplicated(scores[,.(dataset,task,method,prediction_key)]))
setorder(scores,dataset,task,method,-score,prediction_key)
scores[,candidate_rank:=seq_len(.N),by=.(dataset,task,method)]

native_targets <- fread(file.path(CONFIG$paths$input,"task_tables/lrsim_target_recovery.csv"))[
  version=="Final065"]
stopifnot(nrow(native_targets)==42L,!anyDuplicated(native_targets$target_key))
budgets <- native_targets[,.(k_budget=.N,declared_k=unique(k_budget)),by=dataset]
stopifnot(all(budgets$k_budget==budgets$declared_k))
scores <- merge(scores,budgets[,.(dataset,k_budget)],by="dataset",all.x=TRUE)
stopifnot(all(scores[is.na(k_budget)]$dataset=="lrsim_shared_trend"))
scores[is.na(k_budget),k_budget:=0L]
counts <- scores[,.(n_candidates=.N,n_positive=sum(truth_row),n_null=sum(!truth_row),
  k_budget=unique(k_budget)),by=.(dataset,task,method)]
defined <- counts[n_positive>=2L & n_null>=2L,.(dataset,task,method)]
areas <- merge(scores,defined,by=c("dataset","task","method"))[,.(auprc=
  PRROC::pr.curve(scores.class0=score[truth_row],
    scores.class1=score[!truth_row],curve=FALSE)$auc.integral),
  by=.(dataset,task,method)]
precision <- scores[k_budget>0L,.(precision_at_k=
  sum(truth_row[candidate_rank<=k_budget[1L]])/k_budget[1L]),
  by=.(dataset,task,method)]
per_scenario <- merge(merge(counts,areas,by=c("dataset","task","method"),all.x=TRUE),
  precision,by=c("dataset","task","method"),all.x=TRUE)
per_scenario[,status:=fifelse(n_positive==0L,"zero_positive_reference","scored")]
pair <- per_scenario[task=="LR_pair_nomination",
  .(auprc=mean(auprc,na.rm=TRUE),precision_at_k=mean(precision_at_k,na.rm=TRUE),
    defined_inputs=sum(is.finite(auprc))),by=method][match(CONFIG$pair_methods,method)]
event_areas <- per_scenario[task=="projected_fork_event",
  .(auprc=mean(auprc,na.rm=TRUE),defined_inputs=sum(is.finite(auprc))),by=method]
stopifnot(all(pair$defined_inputs==10L),all(event_areas$defined_inputs==6L))

# Recompute each budget hit from its correct-context rank and declared K.
target_files <- list.files(file.path(CONFIG$paths$input,"projected_targets"),
  pattern="^projected_target_recovery[.]csv$",recursive=TRUE,full.names=TRUE)
targets <- rbindlist(lapply(target_files,fread))[evaluation_scope=="own_available"]
targets[,recomputed_hit:=is.finite(best_correct_rank)&best_correct_rank<=k_budget]
stopifnot(all(targets$recomputed_hit==targets$hit_in_budget),
  !anyDuplicated(targets[,.(method,target_key)]))
event_counts <- targets[,.(planned=.N,hits=sum(recomputed_hit)),by=method]
stopifnot(all(event_counts$planned==42L))
events <- merge(event_areas,event_counts,by="method")[match(CONFIG$event_methods,method)]
native_targets[,recomputed_hit:=is.finite(best_correct_rank)&best_correct_rank<=k_budget]
stopifnot(sum(native_targets$recomputed_hit)==events[method=="LineaLinkR_event_score"]$hits)

table1 <- rbindlist(list(
  pair[,.(Task="Molecular pair",Method=CONFIG$pair_labels,AUPRC=sprintf("%.3f",auprc),
    "Budget result"=sprintf("%.3f",precision_at_k))],
  events[,.(Task="Branch event",Method=CONFIG$event_labels,AUPRC=sprintf("%.3f",auprc),
    "Budget result"=paste0(hits,"/",planned))]
))
published <- as.data.table(read.csv(unz(CONFIG$paths$source_data,"Table1.csv"),
  check.names=FALSE,colClasses="character",fileEncoding="UTF-8"))
stopifnot(identical(as.data.frame(table1),as.data.frame(published)))
fwrite(per_scenario,file.path(OUT,"Table1_by_scenario.csv"),na="NA")
fwrite(targets,file.path(OUT,"Table1_event_recovery.csv"),na="NA")
fwrite(table1,file.path(OUT,"Table1_recomputed.csv"))

# Reconstruct physical focal fractions and their strict replicate-separated peaks.
CONTACT <- file.path(CONFIG$paths$input,"contact_validation")
gold <- fread(file.path(CONTACT,"ulipstic_all_mouse_child_contact.csv"))
stopifnot(nrow(gold)==54L,uniqueN(gold$sample_id)==6L,
  identical(sort(unique(gold$time_hours)),c(36L,50L,96L)))
ranges <- gold[,.(n_observed=sum(observable),
  focal_min=min(contact_focal),focal_max=max(contact_focal),
  relative_min=min(contact_relative),relative_max=max(contact_relative)),
  by=.(event_node,child,time_hours)]
profiles <- dcast(ranges,event_node+child~time_hours,
  value.var=c("n_observed","focal_min","focal_max","relative_min","relative_max"))
profiles[,complete:=n_observed_36==2L & n_observed_50==2L & n_observed_96==2L]
profiles[,focal_peak:=complete & focal_min_50>pmax(focal_max_36,focal_max_96)]
profiles[,relative_peak:=complete & relative_min_50>pmax(relative_max_36,relative_max_96)]

# The exported graph profiles use contact_focal, with the same state-row mapping.
contexts <- unique(gold[,.(event_node,child)])[order(event_node,child)]
contexts[,row:=.N-.I+1L]
plotted <- fread(file.path(CONTACT,"main6_all_mouse_contact_profiles.csv"))
recomputed_plot <- merge(gold,contexts,by=c("event_node","child"))[observable==TRUE,
  .(lower=min(contact_focal),upper=max(contact_focal),mean=mean(contact_focal),n=.N),
  by=.(row,time_hours)]
plot_check <- merge(plotted,recomputed_plot,by=c("row","time_hours"),all.x=TRUE,
  suffixes=c("_figure","_recomputed"))
stopifnot(all(plotted$readout=="focal_contact_fraction"),
  isTRUE(all.equal(as.matrix(plot_check[,.(lower_figure,upper_figure,mean_figure,n_figure)]),
    as.matrix(plot_check[,.(lower_recomputed,upper_recomputed,mean_recomputed,n_recomputed)]),
    tolerance=CONFIG$tolerance,check.attributes=FALSE)))

temporal <- fread(file.path(CONTACT,"ulipstic_all_contact_time_predictions.csv"))[
  model=="event_ranking_score"]
temporal <- merge(temporal,profiles[,.(event_node,child,focal_peak,relative_peak)],
  by=c("event_node","child"))
stopifnot(all(temporal$focal_pulse==temporal$focal_peak))
curves <- fread(file.path(CONTACT,"ulipstic_selected_query_calendar_support.csv"))[
  model=="event_ranking_score"]
curve_check <- merge(curves,gold[,.(event_node,child,sample_id,time_hours,contact_focal)],
  by=c("event_node","child","sample_id","time_hours"),all.x=TRUE,
  suffixes=c("_rna_record","_physical_record"))
stopifnot(isTRUE(all.equal(curve_check$contact_focal_rna_record,
  curve_check$contact_focal_physical_record,tolerance=CONFIG$tolerance)))
support <- curves[,.(n_rna_structural=sum(structural_valid),
  n_rna_expression=sum(expression_supported),
  n_supported_hours=uniqueN(time_hours[structural_valid&expression_supported])),
  by=.(event_node,child)]
peak_check <- merge(temporal[focal_peak==TRUE,.(event_node,child,pair_idx,interaction_name,
  source,target,focal_peak,focal_dynamic_status,focal_dynamic_pattern,
  matched_peak=focal_dynamic_status=="observed"&focal_dynamic_pattern=="pulse")],
  support,by=c("event_node","child"))
published_contact <- fread(file.path(CONTACT,"ulipstic_contact_time_summary.csv"))[model=="event_ranking_score"]
stopifnot(nrow(peak_check)==published_contact$n_strict_focal_pulse,
  sum(peak_check$matched_peak)==published_contact$focal_pulse_tp,
  sum(profiles$relative_peak)==published_contact$n_strict_relative_pulse)
cell_key <- fread(file.path(CONTACT,"ulipstic_cell_type_key.csv"))
peak_check <- merge(peak_check,cell_key,by.x="source",by.y="type_index")
setnames(peak_check,"cell_type","source_population")
peak_check <- merge(peak_check,cell_key,by.x="target",by.y="type_index")
setnames(peak_check,"cell_type","target_population")
setorder(peak_check,child)
fwrite(profiles,file.path(OUT,"S13_physical_peak_profiles.csv"),na="NA")
fwrite(peak_check,file.path(OUT,"S13_RNA_peak_check.csv"),na="NA")
writeLines(capture.output(sessionInfo()),file.path(OUT,"sessionInfo.txt"),useBytes=TRUE)
print(table1)
print(peak_check)
cat("Table 1 matches all ",nrow(table1)," SourceData rows; focal contact peaks matched ",
  sum(peak_check$matched_peak),"/",nrow(peak_check),
  "; strict relative-contact peaks ",sum(profiles$relative_peak),".\n",sep="")
