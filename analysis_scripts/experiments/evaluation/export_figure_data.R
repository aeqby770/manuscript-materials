CONFIG <- list(root=normalizePath(Sys.getenv("LLR_WORK_ROOT",getwd()),winslash="/",mustWork=TRUE),
  threads=parallel::detectCores(),evaluation="results/evaluation/final065",
  output="results/paper_tables/current/figure_data",datasets="experiments/config/lrsim11_datasets.txt")
Sys.setenv(OMP_NUM_THREADS=CONFIG$threads,OPENBLAS_NUM_THREADS=CONFIG$threads,MKL_NUM_THREADS=CONFIG$threads)
library(data.table)
setDTthreads(CONFIG$threads)
base <- file.path(CONFIG$root,CONFIG$evaluation)
out <- file.path(CONFIG$root,CONFIG$output)
dir.create(out,recursive=TRUE,showWarnings=FALSE)
datasets <- readLines(file.path(CONFIG$root,CONFIG$datasets),encoding="UTF-8")
label_map <- c(lrsim_donor_confounded="Donor confounding",
 lrsim_donor_ligand_only="Donor: ligand only",lrsim_donor_nested_fork="Donor + nested fork",
 lrsim_donor_receptor_only="Donor: receptor only",lrsim_donor_transient_focal="Donor + transient",
 lrsim_dropout="Dropout",lrsim_focal_branch="Focal branch",lrsim_nested_fork="Nested fork",
 lrsim_sampling_thinning="Thinning",lrsim_shared_trend="Shared trend (null)",
 lrsim_transient_focal="Transient")
datasets <- c(setdiff(datasets,"lrsim_shared_trend"),"lrsim_shared_trend")
dataset_labels <- data.table(dataset=datasets,label=unname(label_map[datasets]),input_id=paste0("I",seq_along(datasets)))
dataset_labels[,label:=paste(input_id,label)]
stopifnot(!anyNA(dataset_labels$label))
methods <- data.table(method=c("LineaLinkR_event_max","LineaLinkR_focal_only",
 "LineaLinkR_matched_no_fork","scComm_no_TF_activity","TraSig_native_max","pooled_LR"),
 label=c("LineaLinkR","Raw focal","Matched focal","scComm (no TF)","TraSig","Pooled LR"))
event_methods <- data.table(method=c("LineaLinkR_event_score","scComm_projected_native_kernel",
 "TraSig_projected_contrast"),label=c("LineaLinkR","scComm projection","TraSig projection"))
metrics <- fread(file.path(base,"dynamic_methods.csv"))[evaluation_scope=="common_candidates"]
nomination <- merge(metrics[task=="LR_pair_nomination"],methods,by="method")
nomination <- merge(nomination,dataset_labels,by="dataset",suffixes=c("_method","_dataset"))
nomination[,c("dataset","method"):=list(factor(label_dataset,levels=dataset_labels$label),
 factor(label_method,levels=methods$label))]
nomination <- dcast(nomination,dataset~method,value.var="auprc",drop=FALSE)
fwrite(nomination,file.path(out,"communication_pair_auprc.csv"))
projection <- merge(metrics[task=="projected_fork_event"],event_methods,by="method")
projection <- merge(projection,dataset_labels,by="dataset",suffixes=c("_method","_dataset"))
projection[,c("dataset","method"):=list(factor(label_dataset,levels=dataset_labels$label),
 factor(label_method,levels=event_methods$label))]
projection <- dcast(projection,dataset~method,value.var="auprc",drop=FALSE)
fwrite(projection,file.path(out,"communication_projected_auprc.csv"))
paths <- list.files(file.path(base,"dynamic_methods"),pattern="^projected_target_recovery[.]csv$",
 recursive=TRUE,full.names=TRUE)
recovery <- rbindlist(lapply(paths,fread),fill=TRUE)[evaluation_scope=="own_available"]
target_info <- fread(file.path(base,"lrsim_target_recovery.csv"))[,
 .(target_key,dataset,ligand,receptor,focal_branch,other_branch)]
target_info[,dataset_order:=match(dataset,datasets)]
setorder(target_info,dataset_order,target_key)
target_info[,row_id:=.I]
target_info <- merge(target_info,dataset_labels,by="dataset",all.x=TRUE)
setorder(target_info,row_id)
recovery <- merge(recovery,target_info[,.(target_key,row_id)],by="target_key")
recovery <- merge(recovery,event_methods,by="method")
recovery[,method:=factor(label,levels=event_methods$label)]
recovery[,rank_ratio:=best_correct_rank/k_budget]
recovery[correct_candidate_available%in%FALSE,rank_ratio:=NA_real_]
recovery[,event_status:=fcase(!correct_candidate_available,"No score",
 !correct_local_direction,"Opposite direction",hit_in_budget,"Within budget",
 default="Beyond budget")]
rank_matrix <- dcast(recovery,row_id~method,value.var="rank_ratio",drop=FALSE)
status_matrix <- dcast(recovery,row_id~method,value.var="event_status",drop=FALSE)
fwrite(rank_matrix,file.path(out,"event_rank_over_budget.csv"))
fwrite(status_matrix,file.path(out,"event_support_direction.csv"))
fwrite(target_info,file.path(out,"event_target_labels.csv"))
cat("FIGURE_DATA",nrow(nomination),"scenarios",nrow(target_info),"fixed_targets\n")

# Complete biological plotting data from the saved formal models.
bio_models <- data.table(dataset=c("GSE143437_muscle_regeneration","GSE274389_signal_recorder",
 "GSE291091_fap_depletion","GSE48968_BMDC_WT"),
 workspace=c("final065","final065","external_validation_fap_blind","external_validation_time_root_layout"),
 meta_kind=c("input","input","fap_posthoc","input"))
bio_models[,{
 x <- readRDS(file.path(CONFIG$root,"cache",workspace,dataset,"linealinkr_complete.rds"))
 knn_manifest <- readRDS(file.path(x$knn_precompute_dir,"knn_precompute_manifest.rds"))
 stopifnot(identical(knn_manifest$cells,x$cells))
 coordinates <- fread(file.path(x$knn_precompute_dir,"umap.csv"),header=FALSE)
 setnames(coordinates,c("umap1","umap2"))
 stopifnot(nrow(coordinates)==length(x$cells))
 meta <- switch(meta_kind,
  input=fread(file.path(x$input_extract_dir,paste0(dataset,"_meta.csv"))),
  fap_posthoc=fread(file.path(CONFIG$root,"results/validation/gse291091/cell_metadata.csv.gz")))
 id_column <- switch(meta_kind,input="barcode_id",fap_posthoc="cell_id")
 setnames(meta,id_column,"cell_id")
 meta <- meta[match(x$cells,cell_id)]
 stopifnot(identical(meta$cell_id,x$cells))
 cells <- cbind(coordinates,data.table(cell_id=x$cells,native_state=x$lineage_clusters,
  native_time=x$pseudotime),meta[,-"cell_id"])
 fwrite(cells,file.path(out,paste0(dataset,"_cells.csv.gz")))
 path_nodes <- data.table(path=rep(seq_along(x$lineage_paths),lengths(x$lineage_paths)),
  state=unlist(x$lineage_paths,use.names=FALSE))
 path_nodes[,position:=seq_len(.N),by=path]
 path_nodes[,parent:=shift(state),by=path]
 edges <- unique(path_nodes[!is.na(parent),.(from=parent,to=state)])
 ancestors <- merge(path_nodes[,.(path,ancestor=state,ancestor_position=position)],
  path_nodes[,.(path,descendant=state,descendant_position=position)],by="path",allow.cartesian=TRUE)[
   descendant_position>=ancestor_position,.(ancestor,descendant)]
 ancestors <- unique(ancestors)
 nodes <- cells[,.(n_owned=.N,mean_state_time=mean(native_time)),by=.(state=native_state)]
 descendant_counts <- merge(ancestors,nodes[,.(descendant=state,n_owned)],by="descendant")[
  ,.(n_descendant_cells=sum(n_owned)),by=.(state=ancestor)]
 nodes <- merge(nodes,descendant_counts,by="state")
 fwrite(nodes,file.path(out,paste0(dataset,"_nodes.csv")))
 fwrite(edges,file.path(out,paste0(dataset,"_edges.csv")))
 weights <- rbindlist(list(as.data.table(x$lineage_reference_weights)[,representation:="G"],
  as.data.table(x$lineage_assignment_weights)[,representation:="Q"],
  as.data.table(x$lineage_priming_weights)[,representation:="P"]))
 weights[,state:=rep(x$lineage_clusters,3L)]
 weights <- melt(weights,id.vars=c("representation","state"),variable.name="path",value.name="weight")
 state_weights <- weights[,.(weight=mean(weight)),by=.(state,representation,path)]
 fwrite(state_weights,file.path(out,paste0(dataset,"_state_GQP.csv")))
 cat("BIO_FIGURE_DATA",dataset,nrow(cells),"cells",nrow(nodes),"states\n")
 NULL
},by=dataset]
wnt_cells <- fread(file.path(out,"GSE274389_signal_recorder_cells.csv.gz"))
reporter <- fread(file.path(CONFIG$root,"results/evaluation/biology_main065_final/wnt_reporter_cells.csv"))
setnames(reporter,"barcode_id","cell_id")
stopifnot(setequal(wnt_cells$cell_id,reporter$cell_id))
wnt_cells <- merge(wnt_cells,reporter,by="cell_id",all.x=TRUE)
fwrite(wnt_cells,file.path(out,"GSE274389_signal_recorder_cells.csv.gz"))

# Fixed WNT readouts and complete paired-context support.
wnt_local <- fread(file.path(CONFIG$root,"results/evaluation/biology_final065/wnt_local_reporter_panels.csv"))[
 readout_scope=="t120"]
wnt_history <- fread(file.path(CONFIG$root,"results/evaluation/biology_final065/wnt_history_event_panels.csv"))
wnt_local[,context:=paste0("F",unit_id,":S",canonical_child,"/S",canonical_sibling," r",target_type+1L)]
wnt_local[,fork:=paste0("F",unit_id,":S",canonical_child,"/S",canonical_sibling)]
wnt_local[,receiver:=paste0("r",target_type+1L)]
wnt_local[,support:=fcase(!outcome_observed,"No paired observation",
 !is.finite(full_ranking),"No prediction",default="Scored")]
wnt_support <- dcast(wnt_local,receiver~fork,value.var="support")
fwrite(wnt_support,file.path(out,"wnt_complete_support_matrix.csv"))
wnt_observed <- wnt_local[outcome_observed%in%TRUE & is.finite(full_ranking)]
wnt_means <- wnt_observed[,.(context,parent=focal_pre,focal=focal_post,sisters=sibling_post)]
fwrite(wnt_means,file.path(out,"wnt_native_reporter_mean_matrix.csv"))
wnt_history <- wnt_history[match(wnt_observed$panel_id,panel_id)]
stopifnot(identical(wnt_history$panel_id,wnt_observed$panel_id))
wnt_directions <- data.table(context=wnt_observed$context,
 full=sign(wnt_history$main_ranking_contrast),
 matched_focal=sign(wnt_history$matched_focal_ranking_contrast),
 current_reporter=sign(wnt_history$reporter_now_contrast),
 recorded_history=sign(wnt_history$recorder_contrast))
fwrite(wnt_directions,file.path(out,"wnt_matched_endpoint_directions.csv"))
# Gene expression is prepared on the same original atlas cells.
muscle_extract <- file.path(CONFIG$root,"cache/final065/GSE143437_muscle_regeneration/extract")
muscle_genes <- readLines(file.path(muscle_extract,"GSE143437_muscle_regeneration_genes.txt"),encoding="UTF-8")
muscle_ids <- readLines(file.path(muscle_extract,"GSE143437_muscle_regeneration_cells.txt"),encoding="UTF-8")
library(Matrix)
muscle_counts <- as(readMM(file.path(muscle_extract,"GSE143437_muscle_regeneration_matrix.mtx")),"CsparseMatrix")
muscle_depth <- Matrix::colSums(muscle_counts)
muscle_gene_values <- as.matrix(log1p(muscle_counts[match(c("Ccl2","Ccr2"),muscle_genes),,drop=FALSE]%*%Diagonal(x=10000/muscle_depth)))
muscle_gene_frame <- data.table(cell_id=muscle_ids,Ccl2=muscle_gene_values[1L,],
 Ccr2=muscle_gene_values[2L,])
muscle_cells <- fread(file.path(out,"GSE143437_muscle_regeneration_cells.csv.gz"))
stopifnot(setequal(muscle_gene_frame$cell_id,muscle_cells$cell_id))
muscle_cells <- merge(muscle_cells,muscle_gene_frame,by="cell_id",all.x=TRUE)
fwrite(muscle_cells,file.path(out,"GSE143437_muscle_regeneration_cells.csv.gz"))

# Complete native biological panels; no additional fitting or cell filtering.
CONFIG$bio <- list(muscle="results/evaluation/biology_final065",
 response="results/validation/final065/tables",depletion="results/validation/gse291091/lineage_response",
 ifn="results/validation/gse48968/external_validation_time_root_layout",
 tables="results/paper_tables/current")
events <- fread(file.path(CONFIG$root,CONFIG$bio$muscle,"muscle_fixed_ccl2_ccr2_events.csv"))
events <- events[evaluable%in%TRUE]
events[,context:=paste0("F",unit_id,":S",child," | ",target_label)]
events[,row_id:=.I]
fwrite(events[,.(row_id,unit_id,child,target,context,theta,late_contrast,persistent_positive_score,late_change_contrast,event_ranking_score,
 resource_rank,focal_pre,focal_post,other_post)],file.path(out,"muscle_ccl2_complete_event_readout.csv"))
muscle_day <- fread(file.path(CONFIG$root,CONFIG$bio$muscle,"muscle_ccl2_observed_day_distribution.csv"))
muscle_day <- merge(events[,.(row_id,unit_id,child,type=target,context)],
 muscle_day[role=="receiver"],by=c("unit_id","child","type"),allow.cartesian=TRUE)
muscle_day[,column:=paste0(phase," d",time_days)]
muscle_day <- dcast(muscle_day,row_id+context~column,value.var="fraction",fill=0)
setorder(muscle_day,row_id)
fwrite(muscle_day,file.path(out,"muscle_ccl2_receiver_day_fractions.csv"))
depletion_cells <- fread(file.path(out,"GSE291091_fap_depletion_cells.csv.gz"))
depletion_counts <- depletion_cells[,.(n_cells=.N),by=.(sample_id,condition,cell_type,native_state)]
depletion_counts[,fraction_of_sample:=n_cells/sum(n_cells),by=sample_id]
fwrite(depletion_counts,file.path(out,"depletion_native_sample_composition.csv"))
depletion_types <- depletion_counts[,.(fraction=sum(fraction_of_sample)),by=.(cell_type,sample_id)]
depletion_types <- dcast(depletion_types,cell_type~sample_id,value.var="fraction",fill=0)
fwrite(depletion_types,file.path(out,"depletion_all_celltype_fractions.csv"))
depletion_mnp <- fread(file.path(CONFIG$root,CONFIG$bio$depletion,"native_MNP_response.csv"))
depletion_mnp <- depletion_mnp[,.(fraction=sum(fraction_of_MNP)),by=.(native_node,sample_id)]
depletion_mnp[,context:=paste0("S",native_node)]
depletion_mnp <- dcast(depletion_mnp,context~sample_id,value.var="fraction",fill=0)
fwrite(depletion_mnp,file.path(out,"depletion_native_MNP_fractions.csv"))
response <- fread(file.path(CONFIG$root,CONFIG$bio$response,"supp_muscle_depletion_validation.csv"))
response[,panel:=paste0(context," (",n_candidates,") K",k)]
response75 <- fread(file.path(CONFIG$root,CONFIG$bio$response,"supp_muscle_depletion_ligand_level.csv"))
response75[,panel:=paste0("All MNP (75) K",k)]
response <- rbindlist(list(response,response75),use.names=TRUE,fill=TRUE)
response[,method:=factor(method,levels=c("llr_full","llr_nofork","llr_unadjusted","llr_previous","timetalk","nichenet","entrain"))]
response_matrix <- dcast(response,panel~method,value.var="enrichment")
fwrite(response_matrix,file.path(out,"depletion_response_enrichment_matrix.csv"))
fwrite(response,file.path(out,"depletion_response_complete_values.csv"))
functional <- fread(file.path(CONFIG$root,CONFIG$bio$tables,"supp_gse48968_functional_gene_matrix.csv"))
functional[,module:=gsub("\n"," ",module,fixed=TRUE)]
functional[,condition:=factor(condition,levels=c("WT","WT replicate","Isolated","Ifnar1 KO","Stat1 KO","Tnfr KO"))]
functional[,module:=factor(module,levels=c("Core antiviral","Peaked inflammatory","Sustained inflammatory"))]
functional <- dcast(functional,module+GENE~condition,value.var="activation_fraction")
setorder(functional,module,GENE)
fwrite(functional,file.path(out,"bmdc_all_functional_genes_matrix.csv"))
state_program <- fread(file.path(CONFIG$root,CONFIG$bio$tables,"supp_gse48968_state_program_matrix.csv"))
state_program[,module:=gsub("\n"," ",module,fixed=TRUE)]
state_program <- dcast(state_program,state~module,value.var="detection")
fwrite(state_program,file.path(out,"bmdc_native_state_programs.csv"))
population_ifn <- fread(file.path(CONFIG$root,CONFIG$bio$ifn,"population_ifn_profiles.csv"))
fwrite(population_ifn,file.path(out,"bmdc_observed_IFNB1_profiles.csv"))
calendar_ifn <- fread(file.path(CONFIG$root,CONFIG$bio$ifn,"lps_calendar_profiles.csv"))[ligand_genes=="IFNB1"&receptor_genes=="IFNAR1_IFNAR2"]
fwrite(calendar_ifn,file.path(out,"bmdc_IFNB1_IFNAR_activity.csv"))
axis_rank <- fread(file.path(CONFIG$root,CONFIG$bio$tables,"gse48968_axis_retrieval.csv"))[
 method%chin%c("LineaLinkR_G_amplitude_Q_layout","TraSig_native_dot_same_resource","maximum_observed_mean_LR_product","pooled_mean_LR_product")]
axis_rank[,label:=factor(method,levels=c("LineaLinkR_G_amplitude_Q_layout","TraSig_native_dot_same_resource",
 "maximum_observed_mean_LR_product","pooled_mean_LR_product"),labels=c("Full","TraSig","Peak LR product","Pooled LR product"))]
fwrite(axis_rank,file.path(out,"bmdc_known_axis_common_ranks.csv"))
cat("BIOLOGICAL_PANEL_EXPORT",nrow(events),"muscle events",nrow(functional),"functional genes\n")

# Saved-model mechanism example: select the highest ranked native pulse event.
CONFIG$method_example <- list(dataset="lrsim_donor_transient_focal",
 workspace="cache/final065/lrsim_donor_transient_focal",
 predictions="results/evaluation/final065/dynamic_methods/lrsim_donor_transient_focal/projected_event_predictions.csv")
example <- readRDS(file.path(CONFIG$root,CONFIG$method_example$workspace,"linealinkr_complete.rds"))
example_events <- fread(file.path(CONFIG$root,CONFIG$method_example$predictions))[
 score_mode=="LineaLinkR_event_score"&evaluable%in%TRUE&is.finite(event_ranking_score)&focal_pulse_weight>0]
setorder(example_events,-event_ranking_score)
example_event <- example_events[1L]
example_unit <- attr(example$fork_contrasts,"event_units")$units[[example_event$unit_id]]
example_children <- data.table(child=example_unit$children,leaves=example_unit$desc,
 post_states=example_unit$segment_post_nodes)
example_child_paths <- example_children[,.(lineage=unlist(leaves)),by=child]
example_windows <- example_children[,.(state=unlist(post_states)),by=child]
example_windows[,role:=fifelse(child==example_event$child,"Focal window","Sister windows")]
example_windows <- rbind(example_windows[,.(state,role)],
 data.table(state=example_unit$segment_pre_nodes,role="Shared parent"),use.names=TRUE)
example_paths <- data.table(path=rep(seq_along(example$lineage_paths),lengths(example$lineage_paths)),
 state=unlist(example$lineage_paths))
example_paths[,position:=seq_len(.N),by=path]
example_paths[,parent:=shift(state),by=path]
example_edges <- unique(example_paths[!is.na(parent),.(from=parent,to=state)])
example_nodes <- data.table(state=example$lineage_clusters)[,.(n_owned=.N),by=state]
example_nodes <- merge(example_nodes,example_windows,by="state",all.x=TRUE)
example_nodes[is.na(role),role:="Outside window"]
example_manifest <- readRDS(file.path(example$knn_precompute_dir,"knn_precompute_manifest.rds"))
stopifnot(identical(example_manifest$cells,example$cells))
example_umap <- fread(file.path(example$knn_precompute_dir,"umap.csv"),header=FALSE)
setnames(example_umap,c("umap1","umap2"))
example_cells <- cbind(example_umap,data.table(cell_id=example$cells,native_state=example$lineage_clusters,
 native_time=example$pseudotime,cell_type=as.character(example$celltype_labels)))
example_cells <- merge(example_cells,example_nodes[,.(native_state=state,role)],by="native_state",all.x=TRUE)
example_weights <- rbindlist(list(as.data.table(example$lineage_reference_weights)[,representation:="G"],
 as.data.table(example$lineage_assignment_weights)[,representation:="Q"],
 as.data.table(example$lineage_priming_weights)[,representation:="P"]))
example_weights[,state:=rep(example$lineage_clusters,3L)]
example_weights <- melt(example_weights,id.vars=c("representation","state"),
 variable.name="path",value.name="weight")
example_state_weights <- example_weights[,.(weight=mean(weight)),by=.(state,representation,path)]
example_curves <- rbindlist(list(as.data.table(example$reference_temporal_lr$scores)[,representation:="G"],
 as.data.table(example$shape_temporal_lr$scores)[,representation:="Q"]))
example_curves <- example_curves[pair_idx==example_event$pair_idx&source==example_event$source&target==example_event$target]
example_bounds <- example$temporal_lr$intermediate$bin_bounds[[1L]]
stopifnot(identical(example_bounds,example$reference_temporal_lr$intermediate$bin_bounds[[1L]]),
 identical(example_bounds,example$shape_temporal_lr$intermediate$bin_bounds[[1L]]))
example_bins <- data.table(bin=seq_len(length(example_bounds)-1L),
 lower=head(example_bounds,-1L),upper=tail(example_bounds,-1L))
example_bins[,mid:=(lower+upper)/2]
example_curves <- merge(merge(example_curves,example_bins,by="bin"),example_child_paths,by="lineage")
example_curves[,role:=fifelse(child==example_event$child,"Focal window","Sister windows")]
example_curves[,activity:=fifelse(structural_valid&expression_supported,score,NA_real_)]
example_point <- as.data.table(example$fork_contrasts)[pair_idx==example_event$pair_idx&
 unit_id==example_event$unit_id&source==example_event$source&target==example_event$target]
stopifnot(max(example_point$focal_pre)-min(example_point$focal_pre)<1e-12)
example_point[,role:=fifelse(child==example_event$child,"Focal window","Sister windows")]
fwrite(example_event,file.path(out,"method_example_event.csv"))
fwrite(example_cells,file.path(out,"method_example_cells.csv.gz"))
fwrite(example_nodes,file.path(out,"method_example_nodes.csv"))
fwrite(example_edges,file.path(out,"method_example_edges.csv"))
fwrite(example_windows,file.path(out,"method_example_windows.csv"))
fwrite(example_weights,file.path(out,"method_example_cell_weights.csv.gz"))
fwrite(example_state_weights,file.path(out,"method_example_state_weights.csv"))
fwrite(example_curves,file.path(out,"method_example_temporal_profiles.csv"))
fwrite(example_bins,file.path(out,"method_example_Q_bins.csv"))
fwrite(example_point,file.path(out,"method_example_shared_parent_points.csv"))
cat("NATIVE_METHOD_EXAMPLE",nrow(example_cells),"cells;",example_event$ligand,example_event$receptor,
 "; parent",example_event$event_node,"child",example_event$child,"; mean mode",example$event_expression_mean_mode,"\n")

# Case-study display data retain every atlas cell and biological sample.
CONFIG$case_display <- list(muscle="cache/final065/GSE143437_muscle_regeneration/linealinkr_complete.rds",
 events="results/evaluation/biology_final065/muscle_fixed_ccl2_ccr2_events.csv",
 depletion="results/validation/gse291091/lineage_response/native_MNP_response.csv")
muscle_model <- readRDS(file.path(CONFIG$root,CONFIG$case_display$muscle))
muscle_events <- fread(file.path(CONFIG$root,CONFIG$case_display$events))
muscle_focus <- muscle_events[evaluable%in%TRUE&target_label=="Resident Macrophages/APCs"&theta>0]
muscle_focus <- muscle_focus[which.max(event_ranking_score)]
muscle_unit <- attr(muscle_model$fork_contrasts,"event_units")$units[[muscle_focus$unit_id]]
muscle_children <- data.table(child=muscle_unit$children,states=muscle_unit$segment_post_nodes)
muscle_windows <- muscle_children[,.(state=unlist(states)),by=child]
muscle_windows[,role:=fifelse(child==muscle_focus$child,"Focal window","Sister windows")]
muscle_windows <- rbind(muscle_windows[,.(state,role)],
 data.table(state=muscle_unit$segment_pre_nodes,role="Shared parent"))
muscle_cells <- fread(file.path(out,"GSE143437_muscle_regeneration_cells.csv.gz"))
muscle_cells <- merge(muscle_cells,muscle_windows[,.(native_state=state,role)],by="native_state",all.x=TRUE)
muscle_cells[is.na(role),role:="Outside window"]
gene_display <- melt(muscle_cells,id.vars=c("cell_id","cell_type"),
 measure.vars=c("Ccl2","Ccr2"),variable.name="gene",value.name="log_expression")
gene_display <- gene_display[,.(mean_log_expression=mean(log_expression),
 detected_fraction=mean(log_expression>0),n_cells=.N),by=.(cell_type,gene)]
muscle_nodes <- fread(file.path(out,"GSE143437_muscle_regeneration_nodes.csv"))
muscle_nodes <- merge(muscle_nodes,muscle_windows,by="state",all.x=TRUE)
muscle_nodes[is.na(role),role:="Outside window"]
muscle_selected <- muscle_events[unit_id==muscle_focus$unit_id&
 source==muscle_focus$source&target==muscle_focus$target]
fwrite(muscle_cells,file.path(out,"muscle_case_cells.csv.gz"))
fwrite(muscle_nodes,file.path(out,"muscle_case_nodes.csv"))
fwrite(muscle_focus,file.path(out,"muscle_case_focus_event.csv"))
fwrite(muscle_selected,file.path(out,"muscle_case_all_siblings.csv"))
fwrite(gene_display,file.path(out,"muscle_case_celltype_gene_expression.csv"))
mnp <- fread(file.path(CONFIG$root,CONFIG$case_display$depletion))
mnp_phenotypes <- mnp[,.(fraction_of_MNP=sum(fraction_of_MNP),n_cells=sum(N)),
 by=.(sample_id,condition,cell_type)]
fwrite(mnp_phenotypes,file.path(out,"depletion_MNP_phenotype_by_sample.csv"))
cat("CASE_DISPLAY",nrow(muscle_cells),"muscle cells; native focal",muscle_focus$child,
 "; MNP types",uniqueN(mnp_phenotypes$cell_type),"; samples",uniqueN(mnp_phenotypes$sample_id),"\n")
print(unique(gene_display$cell_type))
print(mnp_phenotypes[,.(cells=sum(n_cells)),by=sample_id])

CONFIG$case_display$response_predictions <- "results/validation/final065/predictions"
differences102 <- fread(file.path(CONFIG$root,CONFIG$case_display$response_predictions,"method_difference_assignments.csv"))[
 k==10L&context%chin%c("All MNP","APC")&comparison=="Full minus NicheNet"]
differences102[,panel:=fifelse(context=="All MNP","mnp102_nichenet","apc90_nichenet")]
differences75 <- fread(file.path(CONFIG$root,CONFIG$case_display$response_predictions,"entrain_common_assignments.csv"))[k==10L]
differences75 <- dcast(differences75,k+assignment~method,value.var="enrichment")
plot_differences75 <- rbindlist(list(
 differences75[,.(assignment,panel="mnp75_entrain",difference=llr_full-entrain)],
 differences75[,.(assignment,panel="mnp75_nichenet",difference=llr_full-nichenet)]))
plot_differences <- rbind(differences102[,.(assignment,panel,difference)],plot_differences75)
fwrite(plot_differences,file.path(out,"depletion_four_sample_method_differences.csv"))
print(plot_differences[assignment==1L])

# RNA snapshots at the biological-sample level for the fixed muscle axis.
muscle_cells <- fread(file.path(out,"muscle_case_cells.csv.gz"))
muscle_sample_genes <- melt(muscle_cells,id.vars=c("cell_id","sample_id","time_days","cell_type"),
 measure.vars=c("Ccl2","Ccr2"),variable.name="gene",value.name="log_expression")
muscle_sample_genes <- muscle_sample_genes[,.(mean_log_expression=mean(log_expression),
 detected_fraction=mean(log_expression>0),n_cells=.N),by=.(sample_id,time_days,cell_type,gene)]
fwrite(muscle_sample_genes,file.path(out,"muscle_case_sample_gene_expression.csv"))
cat("MUSCLE_BIOLOGICAL_SAMPLES",uniqueN(muscle_sample_genes$sample_id),"\n")

# Recorded class composition within each unsupervised depletion state.
mnp_classes <- fread(file.path(CONFIG$root,"results/validation/gse291091/lineage_response/native_MNP_response.csv"))
mnp_classes <- mnp_classes[,.(n_cells=sum(N)),by=.(native_node,cell_type)]
mnp_classes <- merge(CJ(native_node=sort(unique(mnp_classes$native_node)),cell_type=sort(unique(mnp_classes$cell_type))),mnp_classes,by=c("native_node","cell_type"),all.x=TRUE)
mnp_classes[is.na(n_cells),n_cells:=0L]
mnp_classes[,fraction:=n_cells/fifelse(sum(n_cells)>0,sum(n_cells),NA_real_),by=native_node]
mnp_classes <- dcast(mnp_classes,native_node~cell_type,value.var="fraction",fill=NA_real_)
fwrite(mnp_classes,file.path(out,"depletion_native_MNP_class_composition.csv"))

# Matched current and recorded WNT readouts on the saved native contexts.
CONFIG$wnt_display <- list(local="results/evaluation/biology_final065/wnt_local_reporter_panels.csv",
 history="results/evaluation/biology_final065/wnt_history_event_panels.csv")
wnt_local <- fread(file.path(CONFIG$root,CONFIG$wnt_display$local))[readout_scope=="t120"]
wnt_history <- fread(file.path(CONFIG$root,CONFIG$wnt_display$history))
wnt_matched <- wnt_local[outcome_observed%in%TRUE&is.finite(full_ranking)]
wnt_matched[,context:=paste0("F",unit_id,":S",canonical_child,"/S",canonical_sibling," r",target_type+1L)]
wnt_history_matched <- wnt_history[match(wnt_matched$panel_id,panel_id)]
stopifnot(identical(wnt_history_matched$panel_id,wnt_matched$panel_id))
wnt_endpoint_contrasts <- data.table(context=wnt_matched$context,
 native_score=wnt_history_matched$main_ranking_contrast,
 matched_score=wnt_history_matched$matched_focal_ranking_contrast,
 current_reporter=wnt_history_matched$reporter_now_contrast,
 recorded_history=wnt_history_matched$recorder_contrast)
wnt_endpoint_contrasts[,context_id:=.I]
fwrite(wnt_endpoint_contrasts,file.path(out,"wnt_matched_continuous_contrasts.csv"))
fwrite(unique(wnt_history[,.(unit_id,event_node)]),file.path(out,"wnt_native_fork_nodes.csv"))
wnt_cells <- fread(file.path(out,"GSE274389_signal_recorder_cells.csv.gz"))
wnt_current_nodes <- wnt_cells[time_hours==120,.(mean_irfp=mean(irfp_log1p_cptt),n_reporter_cells=.N),
 by=.(state=native_state)]
fwrite(wnt_current_nodes,file.path(out,"wnt_native_current_reporter_means.csv"))
cat("WNT_CONTINUOUS_CONTEXTS",nrow(wnt_endpoint_contrasts),"\n")
print(wnt_endpoint_contrasts[,.(current=sum(sign(native_score)==sign(current_reporter)),
 history=sum(sign(native_score)==sign(recorded_history)))])


# Native PRROC curves from the exact saved common-candidate rankings.
CONFIG$pr_display <- list(evaluation="results/evaluation/final065")
.libPaths(c(file.path(CONFIG$root,"library"),.libPaths()))
library(PRROC)
ranking_files <- data.table(file=list.files(file.path(CONFIG$root,CONFIG$pr_display$evaluation,"dynamic_methods"),
 pattern="^ranked_predictions[.]csv$",recursive=TRUE,full.names=TRUE))
ranked_display <- ranking_files[,fread(file,encoding="UTF-8"),by=file][evaluation_scope=="common_candidates"]
ranked_display[,n_positive:=sum(truth_row),by=.(dataset,task,method)]
pr_display <- ranked_display[n_positive>0,{
 pr <- PRROC::pr.curve(scores.class0=score[truth_row],scores.class1=score[!truth_row],curve=TRUE)
 data.table(recall=pr$curve[,1L],precision=pr$curve[,2L],recomputed_auprc=pr$auc.integral)
},by=.(dataset,task,method)]
pr_coverage <- ranked_display[,.(n_candidates=.N,n_positive=sum(truth_row)),by=.(dataset,task,method)]
saved_pr_metrics <- fread(file.path(CONFIG$root,CONFIG$pr_display$evaluation,"dynamic_methods.csv"))[
 evaluation_scope=="common_candidates",.(dataset,task,method,auprc)]
pr_coverage <- merge(pr_coverage,saved_pr_metrics,by=c("dataset","task","method"),all.x=TRUE)
pr_values <- unique(pr_display[,.(dataset,task,method,recomputed_auprc)])
pr_values <- merge(pr_values,saved_pr_metrics,by=c("dataset","task","method"))
pr_difference <- max(abs(pr_values$recomputed_auprc-pr_values$auprc),na.rm=TRUE)
stopifnot(pr_difference<1e-10,uniqueN(ranked_display$dataset)==11L)
fwrite(pr_display[,.(dataset,task,method,recall,precision)],file.path(out,"communication_native_pr_curves.csv"))
fwrite(pr_coverage,file.path(out,"communication_native_pr_coverage.csv"))
null_scores <- fread(file.path(CONFIG$root,CONFIG$pr_display$evaluation,"lrsim_contrast_detail.csv"))[
 dataset=="lrsim_shared_trend"&sender_direction%in%TRUE&null_row%in%TRUE&duplicate_emission%in%FALSE]
null_summary <- fread(file.path(base,"lrsim_metrics.csv"))[dataset=="lrsim_shared_trend"]
stopifnot(nrow(null_scores)==null_summary$n_null_scored,sum(null_scores$ranking_score>0)==null_summary$null_rank_positive)
fwrite(null_scores,file.path(out,"communication_complete_native_null.csv"))
cat("NATIVE_PR_FIGURE_DATA",nrow(pr_values),"defined curves;",
 sum(pr_coverage$n_positive==0),"undefined groups; max AUPRC difference",pr_difference,"\n")

# Native windows for the already-recorded best correct predictions.
window_files <- data.table(dataset=unique(target_info$dataset))
window_pred <- window_files[,{
 p <- fread(file.path(base,"dynamic_methods",.BY$dataset,"projected_event_predictions.csv"))
 p[score_mode=="LineaLinkR_event_score"&duplicate_emission%in%FALSE,
 .(prediction_key,event_node,child,unit_id,fork_time,dynamic_time_start,dynamic_time_end,
 persistence_phase_start,persistence_phase_end,dynamic_status,persistence_support,score)]
},by=dataset]
setorder(window_pred,-score)
window_pred <- unique(window_pred,by=c("dataset","prediction_key"))
window_recovery <- fread(file.path(base,"lrsim_target_recovery.csv"))
native_windows <- merge(window_recovery,window_pred,
 by.x=c("dataset","best_correct_prediction"),by.y=c("dataset","prediction_key"),all.x=TRUE)
native_windows <- merge(fread(file.path(out,"event_target_labels.csv"))[,.(target_key,row_id,input_id)],
 native_windows,by="target_key",all.x=TRUE)
setorder(native_windows,row_id)
stopifnot(nrow(native_windows)==42L,all(is.finite(native_windows$child)),
 max(abs(native_windows$ranking_score-native_windows$score))<1e-12)
fwrite(native_windows,file.path(out,"communication_target_native_windows.csv"))
