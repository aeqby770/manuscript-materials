# Evaluate returned events on the fixed generator tree and generating-time scale.
CONFIG <- list(root=normalizePath(Sys.getenv("LLR_WORK_ROOT",getwd()),winslash="/",mustWork=TRUE),
  workspace=Sys.getenv("LLR_SCORE_WORKSPACE","final065"),evaluation=Sys.getenv("LLR_SCORE_EVAL","final065"),
  curve_root=Sys.getenv("LLR_TIME_ACTIVITY_ROOT",""),output=Sys.getenv("LLR_LOCATION_OUT",""),
  python=Sys.getenv("LLR_PYTHON","python"),threads=parallel::detectCores(),min_effective_cells=5L)
Sys.setenv(OMP_NUM_THREADS=CONFIG$threads)
library(data.table)
setDTthreads(CONFIG$threads)
base <- file.path(CONFIG$root,"results/evaluation",CONFIG$evaluation)
OUT <- ifelse(nzchar(CONFIG$output),CONFIG$output,file.path(base,"event_location"))
curves_root <- ifelse(nzchar(CONFIG$curve_root),CONFIG$curve_root,file.path(base,"dynamic_methods"))
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
gen <- file.path(CONFIG$root,"experiments/simulation/generate_lr_fork_sim.py")
design_text <- system2(CONFIG$python,c("-X","utf8",shQuote(gen),"--describe"),stdout=TRUE)
design <- jsonlite::fromJSON(paste(design_text,collapse="\n"))
writeLines(design_text,file.path(OUT,"generator_design.json"),useBytes=TRUE)
datasets <- readLines(file.path(CONFIG$root,"experiments/config/lrsim11_datasets.txt"),encoding="UTF-8")
all_targets <- fread(file.path(base,"lrsim_target_recovery.csv"))
method_names <- c("LineaLinkR_event_score","TraSig_projected_contrast","scComm_projected_native_kernel")
target_plan <- unique(all_targets[,.(dataset,target_key,ligand,receptor,focal_branch,other_branch,k_budget)])
stopifnot(nrow(target_plan)==42L, !anyDuplicated(target_plan$target_key), length(datasets)==11L)
per_input <- lapply(datasets,function(ds) {
  panel <- file.path(base,"dynamic_methods",ds)
  x <- readRDS(file.path(CONFIG$root,"cache",CONFIG$workspace,ds,"linealinkr_complete.rds"))
  truth <- fread(file.path(CONFIG$root,"inputs/revision/sim",paste0(ds,"_truth.csv")))
  meta <- fread(file.path(CONFIG$root,"inputs/revision/sim",paste0(ds,"_cells.csv")))
  meta <- meta[match(x$cells,cell_id)]
  stopifnot(identical(meta$cell_id,x$cells))
  meta[,node:=as.integer(x$lineage_clusters)]
  scenario <- sub("^lrsim_","",ds)
  stopifnot(identical(as.character(truth$role),as.character(design$roles[[scenario]])))
  shape <- if(any(truth$role=="nested")) design$nested_segments else design$two_branch_segments
  segment <- data.table(label=as.character(shape[,1L]),lower=as.numeric(shape[,2L]),upper=as.numeric(shape[,3L]))
  positive <- truth[true_positive==TRUE,.(ligand,receptor,role,fork_pt,focal_branch,other_branch)]
  positive[,parent_segment:=segment$label[match(fork_pt,segment$upper)]]
  positive[,active_start:=fork_pt]
  positive[,active_end:=fifelse(role=="transient",design$pulse$end,1)]
  positive[,plateau_start:=fifelse(role=="transient",design$pulse$plateau_start,NA_real_)]
  positive[,plateau_end:=fifelse(role=="transient",design$pulse$plateau_end,NA_real_)]
  units <- attr(x$fork_contrasts,"event_units")$units
  ut <- data.table(unit_id=seq_along(units),unit=units)
  edges <- ut[,.(event_node=unit[[1L]]$event_node,child=unit[[1L]]$children,
    descendants=unit[[1L]]$desc,post_nodes=unit[[1L]]$segment_post_nodes),by=unit_id]
  counts <- meta[,.N,by=.(node,true_branch)]
  setorder(counts,node,-N,true_branch)
  nodes <- counts[,.(segment=true_branch[1L],purity=N[1L]/sum(N)),by=node]
  nodes <- merge(nodes,meta[,.(at=median(true_pseudotime)),by=node],by="node")
  edges <- merge(edges,nodes[,.(event_node=node,parent_label=segment,parent_purity=purity,parent_time=at)],by="event_node")
  edges <- merge(edges,nodes[,.(child=node,child_label=segment,child_purity=purity,child_time=at)],by="child")
  edges[,split_upper:=min(child_time),by=unit_id]
  edge_sisters <- merge(edges[,.(unit_id,child)],edges[,.(unit_id,sister=child,label=child_label)],by="unit_id",allow.cartesian=TRUE)[child!=sister,
    .(sister_node_labels=paste(sort(unique(label)),collapse="+")),by=.(unit_id,child)]
  edges <- merge(edges,edge_sisters,by=c("unit_id","child"))
  weights <- x$lineage_reference_weights/rowSums(x$lineage_reference_weights)
  leaves <- edges[,.(leaf=unlist(descendants)),by=.(unit_id,child)]
  incidence <- Matrix::sparseMatrix(i=leaves$leaf,j=match(paste(leaves$unit_id,leaves$child),
    paste(edges$unit_id,edges$child)),x=1,dims=c(ncol(weights),nrow(edges)))
  w <- as.data.table(as.matrix(weights%*%incidence));setnames(w,as.character(seq_len(nrow(edges))))
  w[,cell:=.I]
  w <- melt(w,id.vars="cell",variable.name="edge",value.name="mass",variable.factor=FALSE)
  w <- w[mass>=1e-8];w[,edge:=as.integer(edge)]
  bounds <- x$temporal_lr$intermediate$bin_bounds[[1L]]
  w[,c("unit_id","child","cell_type","true_time","node","bin"):=list(edges$unit_id[edge],edges$child[edge],
    meta$cell_type[cell],meta$true_pseudotime[cell],meta$node[cell],findInterval(x$pseudotime[cell],bounds,all.inside=TRUE))]
  bin_support <- w[,.(n_eff=sum(mass)^2/sum(mass^2)),by=.(unit_id,child,bin,cell_type)]
  supported <- dcast(bin_support,unit_id+child+bin~cell_type,value.var="n_eff")
  calibration <- w[cell_type=="TGT",.(true_center=weighted.mean(true_time,mass)),by=.(unit_id,child,bin)]
  calibration <- merge(calibration,supported,by=c("unit_id","child","bin"))[SRC>=CONFIG$min_effective_cells&TGT>=CONFIG$min_effective_cells]
  calibration[,mid:=(bounds[bin]+bounds[bin+1L])/2]
  post <- edges[,.(node=unlist(post_nodes)),by=.(unit_id,child)]
  post_times <- merge(w[cell_type=="TGT"],post,by=c("unit_id","child","node"))[
    ,.(post_start=min(true_time),post_end=max(true_time)),by=.(unit_id,child)]
  predictions <- fread(file.path(panel,"projected_event_predictions.csv"))
  original <- predictions[score_mode=="LineaLinkR_event_score" & sender_direction%in%TRUE & evaluable%in%TRUE & is.finite(score)]
  resource <- original[,.(score=max(score)),by=.(ligand,receptor)]
  resource <- merge(unique(as.data.table(x$pair_gene_map)[,.(ligand=ligand_genes,receptor=receptor_genes)]),
    resource,by=c("ligand","receptor"),all.x=TRUE)
  resource <- merge(resource,truth[,.(ligand,receptor,role,true_positive)],by=c("ligand","receptor"),all.x=TRUE)
  setorder(resource,-score,ligand,receptor,na.last=TRUE)
  resource[,rank:=fifelse(is.finite(score),.I,NA_integer_)]
  k <- nrow(positive)
  resource[,c("dataset","budget","truth_status","in_budget"):=list(ds,k,
    fifelse(is.na(role),"undefined_in_generator","defined"),is.finite(rank)&rank<=k)]
  candidates <- predictions[sender_direction%in%TRUE & evaluable%in%TRUE & is.finite(score) &
    duplicate_emission%in%c(FALSE,NA) & resource_extra%in%c(FALSE,NA)]
  setorder(candidates,score_mode,prediction_key,-score,lineage,pair_idx)
  candidates <- candidates[,.SD[1L],by=.(score_mode,prediction_key)]
  ranks <- fread(file.path(panel,"ranked_predictions.csv"))[task=="projected_fork_event" & evaluation_scope=="own_available",
    .(score_mode=method,prediction_key,candidate_rank,rank_score=score)]
  candidates <- merge(candidates,ranks,by=c("score_mode","prediction_key"))
  stopifnot(nrow(candidates)==nrow(ranks),all(abs(candidates$score-candidates$rank_score)<1e-10))
  ccols <- c("unit_id","child","score_mode","prediction_key","candidate_rank","pair_idx","ligand","receptor",
    "score","truth_row","hit_target_key","focal_set","other_set","fork_time","focal_peak_time","joint_pulse_weight")
  candidates <- merge(candidates[,..ccols],positive[,.(ligand,receptor,role,fork_pt,truth_focal_label=focal_branch,other_branch,parent_segment,active_start,active_end,plateau_start,plateau_end)],by=c("ligand","receptor"),all.x=TRUE)
  candidates <- merge(candidates,edges[,.(unit_id,child,event_node,parent_label,parent_purity,parent_time,child_label,child_purity,
    sister_node_labels,split_upper)],by=c("unit_id","child"))
  candidates <- merge(candidates,post_times,by=c("unit_id","child"),all.x=TRUE)
  candidates[,c("dataset","budget","in_budget"):=list(ds,k,candidate_rank<=k)]
  candidates[,strict_sister:=truth_row%in%TRUE & other_set==other_branch]
  candidates[,correct_parent_and_daughters:=is.finite(fork_pt) & parent_label==parent_segment &
    child_label==truth_focal_label & sister_node_labels==other_branch]
  candidates[,time_bracket:=is.finite(fork_pt)&parent_time<=fork_pt & split_upper>=fork_pt]
  candidates[,fork_bracket:=correct_parent_and_daughters%in%TRUE & time_bracket%in%TRUE]
  candidates[,window_iou:=fifelse(is.finite(active_start)&is.finite(post_start),
    pmax(0,pmin(post_end,active_end)-pmax(post_start,active_start))/
    (pmax(post_end,active_end)-pmin(post_start,active_start)),NA_real_)]
  candidates[,row_id:=.I]
  cal <- merge(candidates[,.(row_id,unit_id,child,fork_time,focal_peak_time)],calibration,by=c("unit_id","child"),allow.cartesian=TRUE)
  calp <- cal[,.(fork_time_cal=if(.N>=2L) approx(mid,true_center,xout=fork_time[1L])$y else NA_real_,
    native_peak_cal=if(.N>=3L && is.finite(focal_peak_time[1L])) approx(mid,true_center,xout=focal_peak_time[1L])$y else NA_real_),by=row_id]
  candidates <- merge(candidates,calp,by="row_id",all.x=TRUE)
  candidates[,fork_abs_error:=abs(fork_time_cal-fork_pt)]
  candidates[,peak_cal:=native_peak_cal]
  candidates[,activity_peak_cal:=NA_real_]
  if (any(positive$role=="transient")) {
    profiles <- fread(file.path(curves_root,ds,"common_temporal_activity_profiles.csv"))[is.finite(score)&is.finite(true_center)]
    activity <- profiles[,.(peak_cal=if(.N>=3L && diff(range(score))>1e-10) true_center[which.max(score)] else NA_real_),
      by=.(method,unit_id,child,pair_idx)]
    activity[,score_mode:=c(scComm_native_kernel="scComm_projected_native_kernel",
      TraSig_native_curve="TraSig_projected_contrast",LineaLinkR_G_activity="LineaLinkR_G_activity")[method]]
    candidates[activity,on=.(score_mode,unit_id,child,pair_idx),peak_cal:=i.peak_cal]
    reference_max <- activity[method=="LineaLinkR_G_activity"]
    reference_max[,score_mode:="LineaLinkR_event_score"]
    candidates[reference_max,on=.(score_mode,unit_id,child,pair_idx),activity_peak_cal:=i.peak_cal]
  }
  candidates[,peak_in_plateau:=is.finite(peak_cal)&is.finite(plateau_start)&peak_cal>=plateau_start&peak_cal<=plateau_end]
  candidates[,pulse_supported:=score_mode!="LineaLinkR_event_score" | (is.finite(joint_pulse_weight)&joint_pulse_weight>0)]
  selected <- candidates[truth_row%in%TRUE,.(reachable=TRUE,strict_reachable=any(strict_sister%in%TRUE),
    fork_reachable=any(fork_bracket%in%TRUE),budget_hit=any(in_budget),
    strict_budget_hit=any(in_budget&strict_sister%in%TRUE),fork_budget_hit=any(in_budget&fork_bracket%in%TRUE),
    timed_recommendation=any(in_budget & is.finite(peak_cal)&pulse_supported),
    activity_estimated=any(in_budget & is.finite(activity_peak_cal)),
    activity_peak_joint=any(in_budget & is.finite(activity_peak_cal) & is.finite(plateau_start) &
      activity_peak_cal>=plateau_start & activity_peak_cal<=plateau_end),
    branch_peak_joint=any(in_budget & peak_in_plateau & pulse_supported),
    strict_peak_joint=any(in_budget & strict_sister%in%TRUE & peak_in_plateau & pulse_supported),
    fork_peak_joint=any(in_budget & fork_bracket%in%TRUE & peak_in_plateau & pulse_supported)),
    by=.(score_mode,target_key=hit_target_key)]
  node_checks <- candidates[is.finite(fork_pt),.(
    local_topology_available=any(correct_parent_and_daughters%in%TRUE),
    local_topology_hit=any(in_budget & correct_parent_and_daughters%in%TRUE),
    generated_available=any(fork_bracket%in%TRUE),
    generated_hit=any(in_budget & fork_bracket%in%TRUE),
    generated_peak_hit=any(in_budget & fork_bracket%in%TRUE & peak_in_plateau & pulse_supported)),
    by=.(score_mode,ligand,receptor)]
  plan <- merge(target_plan[dataset==ds,.(dataset,target_key,ligand,receptor,k_budget,j=1L)],
    data.table(score_mode=method_names,j=1L),by="j",allow.cartesian=TRUE)[,j:=NULL]
  check <- merge(plan,selected,by=c("score_mode","target_key"),all.x=TRUE)
  check <- merge(check,node_checks,by=c("score_mode","ligand","receptor"),all.x=TRUE)
  flags <- c(setdiff(names(selected),c("score_mode","target_key")),setdiff(names(node_checks),c("score_mode","ligand","receptor")))
  check[,(flags):=lapply(.SD,function(a) a%in%TRUE),.SDcols=flags]
  check <- merge(check,positive[,.(ligand,receptor,role,fork_pt,active_start,active_end,plateau_start,plateau_end)],
    by=c("ligand","receptor"),all.x=TRUE)
  fc <- as.data.table(x$fork_contrasts)
  if (!"n_paired_samples"%in%names(fc)) {
    fc[,c("n_paired_samples","n_paired_late_samples","paired_sample_status"):=
      list(NA_integer_,NA_integer_,"archived_scoring")]
  }
  coverage <- fc[,.(events=.N,evaluable=sum(evaluable%in%TRUE),
    paired_evaluable=sum(n_paired_samples>0,na.rm=TRUE),single_sample=sum(n_paired_samples==1,na.rm=TRUE),
    paired_late=sum(n_paired_late_samples>0,na.rm=TRUE),sample_supports=sum(n_paired_samples,na.rm=TRUE)),
    by=paired_sample_status]
  list(candidates=candidates,targets=check,resource=resource,
    coverage=coverage[,dataset:=ds])
})
candidates <- rbindlist(lapply(per_input,"[[","candidates"),fill=TRUE)
targets <- rbindlist(lapply(per_input,"[[","targets"),fill=TRUE)
resources <- rbindlist(lapply(per_input,"[[","resource"),fill=TRUE)
coverage <- rbindlist(lapply(per_input,"[[","coverage"),fill=TRUE)
stopifnot(all(targets[,.N,by=score_mode]$N==42L),all(targets[role=="transient",.N,by=score_mode]$N==12L))
event_summary <- targets[,.(planned=.N,reachable=sum(reachable),budget_hits=sum(budget_hit),
  strict_reachable=sum(strict_reachable),strict_budget_hits=sum(strict_budget_hit),
  node_reachable=sum(local_topology_available),node_budget_hits=sum(local_topology_hit),
  fork_reachable=sum(generated_available),fork_budget_hits=sum(generated_hit),
  label_and_fork_hits=sum(fork_budget_hit)),by=score_mode]
time_summary <- targets[role=="transient",.(planned=.N,branch_budget_hits=sum(budget_hit),
  estimated_on_recommended=sum(timed_recommendation),branch_peak_joint=sum(branch_peak_joint),
  focal_activity_estimated=sum(activity_estimated),focal_activity_joint=sum(activity_peak_joint),
  strict_peak_joint=sum(strict_peak_joint),fork_peak_joint=sum(generated_peak_hit)),by=.(score_mode,dataset)]
budget_rows <- candidates[truth_row%in%TRUE & in_budget%in%TRUE][order(candidate_rank),
  .SD[1L],by=.(dataset,score_mode,hit_target_key)]
timing <- budget_rows[,.(branch_retrieved=.N,calibrated_anchors=sum(is.finite(fork_abs_error)),
  median_abs_fork_error=median(fork_abs_error[is.finite(fork_abs_error)]),
  bracket_contains_split=sum(time_bracket%in%TRUE),post_windows=sum(is.finite(window_iou)),
  median_post_window_iou=median(window_iou[is.finite(window_iou)])),by=score_mode]
fwrite(timing,file.path(OUT,"recommended_fork_and_window_summary.csv"))
fwrite(candidates,file.path(OUT,"event_node_window_records.csv"))
fwrite(targets,file.path(OUT,"target_recovery_by_definition.csv"))
fwrite(event_summary,file.path(OUT,"event_recovery_summary.csv"))
fwrite(time_summary,file.path(OUT,"joint_time_summary.csv"))
fwrite(resources,file.path(OUT,"full_resource_molecular_ranks.csv"))
fwrite(coverage,file.path(OUT,"sample_support.csv"))
print(event_summary)
print(time_summary)
