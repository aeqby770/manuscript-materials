CONFIG <- list(root=normalizePath(Sys.getenv("LLR_WORK_ROOT",getwd()),winslash="/",mustWork=TRUE),
 workspace=Sys.getenv("LLR_GSE_WS","external_validation_time_root_layout"),
 output="results/validation/gse48968", dataset="GSE48968_BMDC_WT", threads=parallel::detectCores())
Sys.setenv(OMP_NUM_THREADS=CONFIG$threads,OPENBLAS_NUM_THREADS=CONFIG$threads,MKL_NUM_THREADS=CONFIG$threads)
library(data.table)
setDTthreads(CONFIG$threads)
OUT <- file.path(CONFIG$root,CONFIG$output,CONFIG$workspace)
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
x <- readRDS(file.path(CONFIG$root,"cache",CONFIG$workspace,CONFIG$dataset,"linealinkr_complete.rds"))
meta <- fread(file.path(x$input_extract_dir,paste0(CONFIG$dataset,"_meta.csv")),encoding="UTF-8")
stopifnot(identical(meta$barcode_id,x$cells))
pmap <- as.data.table(x$pair_gene_map)
fwrite(pmap,file.path(OUT,"all_lr_pairs.csv"))
event <- merge(as.data.table(x$fork_contrasts),pmap,by="pair_idx",all.x=TRUE)
fwrite(event,file.path(OUT,"all_native_events.csv.gz"))
scores <- event[evaluable==TRUE & is.finite(event_ranking_score),
 .(score=max(event_ranking_score),n_evaluable_events=.N),by=.(pair_idx)]
scores <- merge(pmap,scores,by="pair_idx",all.x=TRUE)
scores[,rank:=frank(-score,ties.method="min",na.last="keep")]
scores[,n_ranked_pairs:=sum(is.finite(score))]
scores[,covered:=is.finite(score)]
fwrite(scores,file.path(OUT,"native_pair_ranking.csv"))
axis <- scores[ligand_genes=="IFNB1" & receptor_genes=="IFNAR1_IFNAR2"]
stopifnot(nrow(axis)==1L)
# This is a molecular retrieval readout. Untested resource pairs are not negatives.
axis[,`:=`(hit_k10=rank<=10L,hit_k20=rank<=20L,truth_source="Ifnar1 knockout, isolated culture, IFNB stimulation and timed secretion blockade; Shalek 2014")]
fwrite(axis,file.path(OUT,"validated_axis_retrieval.csv"))
calendar <- copy(as.data.table(x$calendar_temporal_lr$scores))
calendar[,native_node:=x$calendar_temporal_lr$profile_nodes[lineage]]
population_lineage <- which(lengths(x$calendar_temporal_lr$lineage_paths)==1L)
calendar <- merge(calendar,pmap,by="pair_idx",all.x=TRUE)
calendar[,`:=`(experimental_block=biological_sample,sampling_unit=x$calendar_temporal_lr$sampling_unit)]
calendar[,biological_sample:=NULL]
fwrite(calendar,file.path(OUT,"all_calendar_profiles.csv.gz"))
lps <- calendar[lineage%in%population_lineage & structural_valid%in%TRUE &
 (grepl("^LPS_[1246]h$",experimental_block)|experimental_block=="Unstimulated")]
lps_profile <- lps[,.(signed_score=sum(score*weight_mass_cells)/sum(weight_mass_cells),
 positive_activity=sum(score_positive*weight_mass_cells)/sum(weight_mass_cells),
 n_paths=.N,effective_cells=sum(n_eff_cells)),by=.(pair_idx,ligand_genes,receptor_genes,observed_time)]
setorder(lps_profile,pair_idx,observed_time)
fwrite(lps_profile,file.path(OUT,"lps_calendar_profiles.csv"))
dynamics <- merge(as.data.table(x$calendar_fork_dynamics),pmap,by="pair_idx",all.x=TRUE)
fwrite(dynamics,file.path(OUT,"all_calendar_event_dynamics.csv.gz"))
axis_dynamics <- dynamics[ligand_genes=="IFNB1" & receptor_genes=="IFNAR1_IFNAR2"]
axis_event_rows <- as.data.table(x$fork_contrasts)[pair_idx==axis$pair_idx]
lps_calendar <- x$calendar_temporal_lr
lps_calendar$scores <- as.data.table(x$calendar_temporal_lr$scores)[pair_idx==axis$pair_idx &
 (grepl("^LPS_[1246]h$",biological_sample)|biological_sample=="Unstimulated")]
pam_calendar <- x$calendar_temporal_lr
pam_calendar$scores <- as.data.table(x$calendar_temporal_lr$scores)[pair_idx==axis$pair_idx &
 (grepl("^PAM_[1246]h$",biological_sample)|biological_sample=="Unstimulated")]
pic_calendar <- x$calendar_temporal_lr
pic_calendar$scores <- as.data.table(x$calendar_temporal_lr$scores)[pair_idx==axis$pair_idx &
 (grepl("^PIC_[1246]h$",biological_sample)|biological_sample=="Unstimulated")]
.libPaths(c(file.path(CONFIG$root,"library"),.libPaths()))
lps_dynamics <- as.data.table(LineaLinkR::SummarizeForkDynamics(axis_event_rows,lps_calendar,x$lineage_paths,
 sample_labels=lps_calendar$scores$biological_sample,time_coordinate="observed_calendar"))
lps_dynamics[,stimulation:="LPS"]
pam_dynamics <- as.data.table(LineaLinkR::SummarizeForkDynamics(axis_event_rows,pam_calendar,x$lineage_paths,
 sample_labels=pam_calendar$scores$biological_sample,time_coordinate="observed_calendar"))
pam_dynamics[,stimulation:="PAM"]
pic_dynamics <- as.data.table(LineaLinkR::SummarizeForkDynamics(axis_event_rows,pic_calendar,x$lineage_paths,
 sample_labels=pic_calendar$scores$biological_sample,time_coordinate="observed_calendar"))
pic_dynamics[,stimulation:="PIC"]
descriptive <- rbindlist(list(lps_dynamics,pam_dynamics,pic_dynamics),use.names=TRUE)
descriptive[,inference_unit:="descriptive_experimental_time_profiles"]
fwrite(descriptive,file.path(OUT,"validated_axis_stimulus_specific_dynamics.csv"))


population_temporal <- x$calendar_temporal_lr
population_temporal$scores <- copy(as.data.table(population_temporal$scores)[lineage%in%population_lineage])
population_temporal$scores[,lineage:=1L]
population_temporal$intermediate$bin_bounds <- list(
 x$calendar_temporal_lr$intermediate$bin_bounds[[population_lineage]])
population_temporal$lineage_paths <- list(c(0L,1L))
population_temporal$observation_scope <- "observed_celltype_population"
population_events <- as.data.frame(unique(as.data.table(x$fork_contrasts)[,
 .(pair_idx,source,target,unit_id=0L,event_node=0L,child=1L,
 sibling_children="2",evaluable=TRUE,sibling_contrast=0)]))
attr(population_events,"event_units") <- list(
 activity_transform=attr(x$fork_contrasts,"event_units")$activity_transform)
population <- as.data.table(LineaLinkR::SummarizeForkDynamics(
 population_events,population_temporal,population_temporal$lineage_paths,
 sample_labels=population_temporal$scores$biological_sample,
 time_coordinate="observed_calendar",background="focal_only",
 condition_labels=population_temporal$scores$condition,reference_conditions="none"))
population <- merge(population,pmap,by="pair_idx")
population[,c("unit_id","event_node","child","sibling_children","evaluable","sibling_contrast"):=NULL]
fwrite(population,file.path(OUT,"population_stimulus_dynamics.csv"))
population_axis <- population[pair_idx==axis$pair_idx]
fwrite(population_axis,file.path(OUT,"validated_axis_population_timing.csv"))

timing <- data.table(dataset=CONFIG$dataset,workspace=CONFIG$workspace,
 root_time_hours=meta$time_hours[x$start_cell],root_prior=x$root_evidence_mode,
 sampled_times=paste(sort(unique(meta$time_hours)),collapse=","),
 perturbation_supported_window="before 1 h; 0 h block strong, 1 h block intermediate, 2 h block mostly preserves antiviral activation",
 n_axis_events=nrow(axis_event_rows),n_evaluable_axis_events=sum(axis_event_rows$evaluable),
 n_condition_specific_event_profiles=nrow(axis_dynamics),
 lps_population_transition_lower=population_axis[dynamic_condition=="LPS",dynamic_transition_lower],
 lps_population_transition_upper=population_axis[dynamic_condition=="LPS",dynamic_transition_upper],
 pic_population_transition_lower=population_axis[dynamic_condition=="PIC",dynamic_transition_lower],
 pic_population_transition_upper=population_axis[dynamic_condition=="PIC",dynamic_transition_upper],
 n_condition_specific_focal_peaks=sum(is.finite(axis_dynamics$focal_peak_time)),
 n_lps_descriptive_focal_peaks=sum(is.finite(lps_dynamics$focal_peak_time)),
 lps_descriptive_peak_hours=paste(sort(unique(na.omit(lps_dynamics$focal_peak_time))),collapse=","),
 pic_descriptive_peak_hours=paste(sort(unique(na.omit(pic_dynamics$focal_peak_time))),collapse=","),
 n_condition_specific_focal_changes=sum(is.finite(axis_dynamics$focal_dynamic_change_time)),
 native_timing_status="reported without converting RNA activity peaks into causal secretion times",
 inference_unit=x$fork_inference_unit,calendar_sampling_unit=x$calendar_temporal_lr$sampling_unit)
fwrite(timing,file.path(OUT,"timing_validation.csv"))
cell_state <- data.table(gsm=x$cells,cell_index=seq_along(x$cells),
 state=x$lineage_clusters,stimulation=meta$stimulation,time_hours=meta$time_hours,
 pseudotime=x$pseudotime,reference_lineage=max.col(x$lineage_reference_weights,ties.method="first"),
 assignment_lineage=max.col(x$lineage_assignment_weights,ties.method="first"))
fwrite(cell_state,file.path(OUT,"fitted_cell_states.csv.gz"))
modules <- fread(file.path(CONFIG$root,CONFIG$output,"all_cell_module_detection.csv.gz"))
modules <- merge(cell_state,modules[,.(gsm,CLUSTER,fraction_detected,n_genes)],by="gsm")
programs <- modules[,.(n_cells=.N,mean_fraction_detected=mean(fraction_detected),median_fraction_detected=median(fraction_detected),
 mean_time_hours=mean(time_hours)),by=.(state,stimulation,CLUSTER)]
fwrite(programs,file.path(OUT,"state_programs.csv"))
state_composition <- cell_state[,.(n_cells=.N,mean_pseudotime=mean(pseudotime)),by=.(state,stimulation,time_hours)]
fwrite(state_composition,file.path(OUT,"state_experimental_composition.csv"))
# Functional annotation of every native IFN event, using the author's fixed
# antiviral module. These overlapping RNA states are not independent replicates.
native_events <- unique(event[pair_idx==axis$pair_idx,
 .(unit_id,event_node,child,sibling_children,event_ranking_score,evaluable)])
native_events[,query_id:=.I]
fwrite(native_events[,.(unit_id,event_node,child,sibling_children,evaluable,event_ranking_score,event_row=query_id)],file.path(OUT,"primary_native_ifn_queries.csv"))
paths <- data.table(lineage=rep(seq_along(x$lineage_paths),lengths(x$lineage_paths)),
 node=as.integer(unlist(x$lineage_paths,use.names=FALSE)))
paths[,position:=seq_len(.N),by=lineage]
descendants <- merge(paths[,.(lineage,side_child=node,start_position=position)],
 paths,by="lineage",allow.cartesian=TRUE)
descendants <- unique(descendants[position>=start_position,.(side_child,node)])
focal_panels <- native_events[,.(query_id,unit_id,event_node,child,event_ranking_score,evaluable,
 side_child=child,side="focal")]
sister_panels <- native_events[,.(side_child=as.integer(unlist(strsplit(sibling_children,",",fixed=TRUE)))),
 by=.(query_id,unit_id,event_node,child,event_ranking_score,evaluable)]
sister_panels[,side:="sisters"]
panels <- merge(rbindlist(list(focal_panels,sister_panels)),descendants,by="side_child",allow.cartesian=TRUE)
annotated <- merge(panels,modules[CLUSTER=="Id",.(node=state,gsm,fraction_detected,stimulation,time_hours)],
 by="node",allow.cartesian=TRUE)
annotation <- annotated[,.(n_cells=.N,mean_core_gene_detection=mean(fraction_detected),
 mean_time_hours=mean(time_hours),fraction_LPS=mean(stimulation=="LPS"),
 fraction_PAM=mean(stimulation=="PAM"),fraction_PIC=mean(stimulation=="PIC")),
 by=.(unit_id,event_node,child,event_ranking_score,evaluable,side)]
fwrite(annotation,file.path(OUT,"all_native_ifn_event_program_annotation.csv"))

clock_check <- cell_state[,.(spearman=cor(pseudotime,time_hours,method="spearman"),n_cells=.N,n_times=uniqueN(time_hours)),by=stimulation]
fwrite(clock_check,file.path(OUT,"observed_time_ordering.csv"))
print(axis[,.(ligand_genes,receptor_genes,score,rank,n_ranked_pairs,hit_k10,hit_k20)])
print(timing)
