CONFIG<-list(root=normalizePath(Sys.getenv("LLR_WORK_ROOT",getwd()),winslash="/",mustWork=TRUE),threads=parallel::detectCores(),
 workspace=Sys.getenv("LLR_CALENDAR_INPUT_WS","cache/final065"),dataset="GSE143437_muscle_regeneration",
 output=Sys.getenv("LLR_CALENDAR_OUT","results/evaluation/biology_final065"),time_column="time_days",sample_column="sample_id",
 min_n_eff=5,min_pct=0.05,lr_scale=10000,
 scope="All native events and complete muscle data; development cohort, no external timing truth",
 time_object="Observed sampling day, separate from state pseudotime",
 unit="Biological sample, each sampled once")
Sys.setenv(OMP_NUM_THREADS=CONFIG$threads)
library(data.table)
library(Matrix)
.libPaths(c(file.path(CONFIG$root,"library"),.libPaths()))
library(LineaLinkR)
setDTthreads(CONFIG$threads)
OUT<-file.path(CONFIG$root,CONFIG$output)
e<-asNamespace("LineaLinkR")

reference<-readRDS(file.path(CONFIG$root,CONFIG$workspace,CONFIG$dataset,"calendar_assignment_reference.rds"))
x<-readRDS(file.path(CONFIG$root,CONFIG$workspace,CONFIG$dataset,"linealinkr_complete.rds"))
meta<-fread(file.path(x$input_extract_dir,paste0(CONFIG$dataset,"_meta.csv")),encoding="UTF-8")
stopifnot(identical(meta$barcode_id,x$cells))
times<-sort(unique(meta[[CONFIG$time_column]]))
calendar<-x$calendar_temporal_lr
curves<-as.data.table(calendar$scores)
samples<-unique(curves[,.(biological_sample,observed_time)])
stopifnot(identical(sort(unique(curves$observed_time)),times),
 samples[,all(.N==1L),by=biological_sample]$V1,
 setequal(samples$biological_sample,as.character(meta[[CONFIG$sample_column]])))
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
saveRDS(calendar,file.path(OUT,"muscle_calendar_lr.rds"))
message("Calendar LR complete: ",nrow(curves))
repaired<-x$calendar_fork_dynamics
before<-reference$dynamics
stopifnot(identical(calendar$observation_scope,"observed_native_branch_subtrees"))
diagnostic<-rbindlist(list(
 data.table(model="previous_assigned_calendar",status=before$dynamic_status,pattern=before$dynamic_pattern),
 data.table(model="observed_subtree_calendar",status=repaired$dynamic_status,pattern=repaired$dynamic_pattern)),
 use.names=TRUE)[,.N,by=.(model,status,pattern)]
summary<-data.table(dataset=CONFIG$dataset,n_samples=nrow(samples),n_observed_times=length(times),
 n_native_events=nrow(repaired),n_state_evaluable=sum(repaired$evaluable%in%TRUE),
 n_calendar_before=sum(before$dynamic_status=="observed"),
 n_calendar_observed=sum(repaired$dynamic_status=="observed"),
 n_focal_calendar_observed=sum(repaired$focal_dynamic_status=="observed"),
 n_calendar_beyond_state_support=sum(repaired$dynamic_status=="observed"&!(repaired$evaluable%in%TRUE)),
 n_calendar_varying=sum(repaired$dynamic_status=="observed"&repaired$dynamic_model!="constant"),
 n_calendar_pulses=sum(repaired$dynamic_pattern=="pulse",na.rm=TRUE))
saveRDS(list(config=CONFIG,calendar_fork_dynamics=repaired,summary=summary),
 file.path(OUT,"muscle_calendar_dynamics.rds"))
fwrite(diagnostic,file.path(OUT,"muscle_calendar_model_counts.csv"))
fwrite(summary,file.path(OUT,"muscle_calendar_dynamics_summary.csv"))
pmap<-as.data.table(x$pair_gene_map)[,.(pair_idx,ligand=ligand_genes,receptor=receptor_genes)]
calendar_output<-merge(as.data.table(repaired)[,c("pair_idx","source","target","event_node","child",
 "evaluable",grep("dynamic_|pulse_bic_gain|focal_peak_time",names(repaired),value=TRUE)),with=FALSE],
 pmap,by="pair_idx")
fwrite(calendar_output,file.path(OUT,"muscle_all_calendar_event_profiles.csv"))
print(summary)
print(diagnostic)
