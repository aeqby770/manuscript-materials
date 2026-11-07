CONFIG <- list(root=normalizePath(Sys.getenv("LLR_WORK_ROOT",getwd()),winslash="/",mustWork=TRUE),threads=parallel::detectCores(),
  workspace=Sys.getenv("LLR_BIO_WS","cache/final065"),
  output=Sys.getenv("LLR_BIO_OUT","results/evaluation/biology_final065"),
  dataset="GSE143437_muscle_regeneration",early_end=3,min_n_eff=5)
Sys.setenv(OMP_NUM_THREADS=CONFIG$threads)
library(data.table)
library(Matrix)
setDTthreads(CONFIG$threads)
OUT <- file.path(CONFIG$root,CONFIG$output)
x <- readRDS(file.path(CONFIG$root,CONFIG$workspace,CONFIG$dataset,"linealinkr_complete.rds"))
meta <- fread(file.path(x$input_extract_dir,paste0(CONFIG$dataset,"_meta.csv")),encoding="UTF-8")
fixed <- fread(file.path(OUT,"muscle_fixed_ccl2_ccr2_events.csv"))
stopifnot(identical(meta$barcode_id,x$cells))
units <- data.table(unit_id=seq_along(attr(x$fork_contrasts,"event_units")$units),
  unit=attr(x$fork_contrasts,"event_units")$units)
desc <- units[,.(child=rep(unit[[1L]]$children,lengths(unit[[1L]]$desc)),
  lineage=unlist(unit[[1L]]$desc,use.names=FALSE)),by=unit_id]
children <- unique(desc[,.(unit_id,child)])
children[,child_index:=.I]
desc <- merge(desc,children,by=c("unit_id","child"))
membership <- sparseMatrix(i=desc$lineage,j=desc$child_index,x=1,
  dims=c(ncol(x$lineage_weights),nrow(children)))
weights <- as.matrix(x$lineage_weights %*% membership)/rowSums(x$lineage_weights)
pre <- merge(children,units[,.(node=unit[[1L]]$segment_pre_nodes),by=unit_id],
  by="unit_id",allow.cartesian=TRUE)
pre[,phase:="pre"]
post <- units[,.(child=rep(unit[[1L]]$children,lengths(unit[[1L]]$segment_post_nodes)),
  node=unlist(unit[[1L]]$segment_post_nodes,use.names=FALSE)),by=unit_id]
post <- merge(post,children,by=c("unit_id","child"))
post[,phase:="post"]
groups <- unique(rbindlist(list(pre,post),use.names=TRUE))
cells <- data.table(cell_index=seq_along(x$cells),node=x$lineage_clusters,
  type=x$celltypes,time_days=meta$time_days,sample_id=meta$sample_id)
members <- merge(cells,groups,by="node",allow.cartesian=TRUE)
members[,weight:=weights[cbind(cell_index,child_index)]]
switch(as.character(identical(attr(x$fork_contrasts,"event_units")$prefork_weighting,"shared")),"TRUE"={
  parent <- members[phase=="pre",.(parent_mass=sum(weight),n_children=uniqueN(child)),by=.(unit_id,cell_index)]
  members <- merge(members,parent,by=c("unit_id","cell_index"),all.x=TRUE)
  members[phase=="pre",weight:=parent_mass/n_children]
})
populations <- unique(rbindlist(list(
  unique(fixed[,.(unit_id,child,type=source,role="source")]),
  unique(fixed[,.(unit_id,child,type=target,role="receiver")]))))
members <- merge(members,populations,by=c("unit_id","child","type"),allow.cartesian=TRUE)
calendar <- members[,.(n_cells=.N,mass=sum(weight),n_eff=sum(weight)^2/sum(weight^2),
  mean_day=sum(weight*time_days)/sum(weight),
  early_fraction=sum(weight*(time_days<=CONFIG$early_end))/sum(weight),
  n_samples=uniqueN(sample_id)),by=.(unit_id,child,type,role,phase)]
calendar_day <- members[,.(mass=sum(weight)),by=.(unit_id,child,type,role,phase,time_days)]
calendar_day[,fraction:=mass/sum(mass),by=.(unit_id,child,type,role,phase)]
wide <- dcast(calendar,unit_id+child+type+role~phase,
  value.var=c("n_eff","mean_day","early_fraction","n_samples"))
source <- wide[role=="source",!"role"]
setnames(source,"type","source")
setnames(source,setdiff(names(source),c("unit_id","child","source")),
  paste0("source_",setdiff(names(source),c("unit_id","child","source"))))
target <- wide[role=="receiver",!"role"]
setnames(target,"type","target")
setnames(target,setdiff(names(target),c("unit_id","child","target")),
  paste0("receiver_",setdiff(names(target),c("unit_id","child","target"))))
result <- merge(fixed,source,by=c("unit_id","child","source"),all.x=TRUE)
result <- merge(result,target,by=c("unit_id","child","target"),all.x=TRUE)
fwrite(result,file.path(OUT,"muscle_ccl2_observed_time_windows.csv"))
fwrite(calendar_day,file.path(OUT,"muscle_ccl2_observed_day_distribution.csv"))
print(result[evaluable %in% TRUE,.(target_label,event_node,child,event_ranking_score,
  source_mean_day_pre,source_mean_day_post,receiver_mean_day_pre,receiver_mean_day_post,
  source_early_fraction_pre,source_early_fraction_post,
  receiver_early_fraction_pre,receiver_early_fraction_post)])
