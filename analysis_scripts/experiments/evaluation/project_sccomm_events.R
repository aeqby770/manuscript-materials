CONFIG <- list(root=normalizePath(Sys.getenv("LLR_WORK_ROOT",getwd()),winslash="/",mustWork=TRUE),
  threads=parallel::detectCores(),dataset=commandArgs(TRUE)[1L],tolerance=1e-10,
  base=file.path("results/evaluation",Sys.getenv("LLR_SCORE_EVAL","final065")),
  workspace=file.path("cache",Sys.getenv("LLR_SCORE_WORKSPACE","final065")),
  native="results/evaluation/final065/sccomm/no_TF_activity",
  min_weight=1e-8,min_n_eff=5,score_mode="scComm_projected_native_kernel")
Sys.setenv(OMP_NUM_THREADS=CONFIG$threads,OPENBLAS_NUM_THREADS=CONFIG$threads,
  MKL_NUM_THREADS=CONFIG$threads)
library(data.table)
library(Matrix)
.libPaths(c(Sys.getenv("LLR_PACKAGE_LIBRARY",file.path(CONFIG$root,"library")),file.path(CONFIG$root,"library"),.libPaths()))
library(LineaLinkR)
setDTthreads(CONFIG$threads)
base <- file.path(CONFIG$root,CONFIG$base)
out <- file.path(base,"dynamic_methods",CONFIG$dataset)
x <- readRDS(file.path(CONFIG$root,CONFIG$workspace,CONFIG$dataset,"linealinkr_complete.rds"))
projection_source <- Sys.getenv("LLR_REUSE_PROJECTED_FROM","")
binding <- fread(file.path(out,"projection_binding.csv"))
old_grid <- file.path(projection_source,CONFIG$dataset,"sccomm_projected_context_scores.csv")
reused_sc <- nzchar(projection_source) && identical(binding$source,projection_source) &&
  isTRUE(binding$matching_geometry) && file.exists(old_grid) && !("--time-profiles"%in%commandArgs(TRUE))
all_predictions <- fread(file.path(out,"projected_event_predictions.csv"))[score_mode!=CONFIG$score_mode]
detail <- all_predictions[score_mode=="LineaLinkR_event_score"]
units <- attr(x$fork_contrasts,"event_units")$units
if (reused_sc) {
  grid <- fread(old_grid)
  equivalence <- fread(file.path(projection_source,CONFIG$dataset,"sccomm_projection_equivalence.csv"))
  stopifnot(max(equivalence$absolute_error)<CONFIG$tolerance)
  fwrite(equivalence,file.path(out,"sccomm_projection_equivalence.csv"))
} else {
native <- readRDS(file.path(CONFIG$root,CONFIG$native,CONFIG$dataset,"sccomm_native.rds"))
stopifnot(identical(colnames(native$ccires$sender_lr),x$cells))
extract_dir <- file.path(CONFIG$root,"cache/final065",CONFIG$dataset,"extract")
genes <- readLines(file.path(extract_dir,paste0(CONFIG$dataset,"_genes.txt")),encoding="UTF-8")
cells <- readLines(file.path(extract_dir,paste0(CONFIG$dataset,"_cells.txt")),encoding="UTF-8")
counts <- as(readMM(file.path(extract_dir,paste0(CONFIG$dataset,"_matrix.mtx"))),"CsparseMatrix")
counts <- counts[,match(x$cells,cells),drop=FALSE]
expression <- log1p(counts %*% Diagonal(x=10000/Matrix::colSums(counts)))
pair_keys <- rownames(native$ccires$sender_lr)
pairs <- data.table(pair_key=pair_keys)
pairs[,c("ligand","receptor"):=tstrsplit(pair_key,"---",fixed=TRUE)]
pairs <- merge(pairs,as.data.table(x$pair_gene_map)[,.(pair_idx,ligand=ligand_genes,receptor=receptor_genes)],
  by=c("ligand","receptor"))
pairs <- pairs[match(pair_keys,pair_key)]
stopifnot(identical(pairs$pair_key,pair_keys))
ligand <- as.matrix(expression[match(pairs$ligand,genes),,drop=FALSE])
receptor <- as.matrix(expression[match(pairs$receptor,genes),,drop=FALSE])
u <- sqrt(ligand)
v <- sqrt(receptor)
n_cells <- length(x$cells)
n_positive <- rowSums(u>0)*rowSums(v>0)
nonzero_mean <- rowSums(u)*rowSums(v)/n_positive
nonzero_sd <- sqrt(pmax(0,(rowSums(ligand)*rowSums(receptor)-
  n_positive*nonzero_mean^2)/(n_positive-1)))
fraction <- n_positive/n_cells^2
kernel_mean <- nonzero_mean*fraction
kernel_sd <- sqrt((nonzero_sd^2+(1-fraction)*nonzero_mean^2)*fraction)
gene_sd <- sqrt(matrixStats::rowSds(ligand)*matrixStats::rowSds(receptor))
source_index <- rep(seq_len(n_cells),times=n_cells)
target_index <- rep(seq_len(n_cells),each=n_cells)
# Exact deterministic CalcCCI kernel from saved fitted weights, with no fitting
# or permutation. Standardization and clipping follow the authors' code.
kernel <- sweep(sweep(u[,source_index]*v[,target_index],
  1L,kernel_mean,"-"),1L,kernel_sd,"/")
kernel <- pmin(pmax(kernel,0),5)
kernel[gene_sd<0.1 | fraction<=0.0001,] <- 0
stopifnot(all(is.finite(kernel)))
labels <- as.character(x$celltype_labels)
type_map <- unique(data.table(cell_type=x$celltypes,label=labels))
source_type <- type_map[label=="SRC",cell_type]
target_type <- type_map[label=="TGT",cell_type]
stopifnot(length(source_type)==1L,length(target_type)==1L)
native_weight <- native$weights$totalweight["SRC","TGT",]
global_mean <- rowMeans(kernel[,labels[source_index]=="SRC" & labels[target_index]=="TGT",drop=FALSE])*native_weight
native_mean <- as.numeric(native$ccires$lrscore["SRC","TGT",])
equivalence <- data.table(dataset=CONFIG$dataset,pair_key=pair_keys,
  native_score=native_mean,reaggregated_score=global_mean,
  absolute_error=abs(native_mean-global_mean))
stopifnot(max(equivalence$absolute_error)<CONFIG$tolerance)
fwrite(equivalence,file.path(out,"sccomm_projection_equivalence.csv"))
switch(as.character("--time-profiles"%in%commandArgs(TRUE)), "TRUE" =
  source(file.path(CONFIG$root,"experiments/evaluation/export_common_temporal_curves.R"),
    local=TRUE,encoding="UTF-8"))

unit_table <- data.table(unit_id=seq_along(units),unit=units)
children <- unit_table[,.(child=unit[[1L]]$children,descendants=unit[[1L]]$desc),by=unit_id]
child_lineages <- children[,.(leaf=unlist(descendants)),by=.(unit_id,child)]
child_columns <- children[,.(child,lineage=seq_len(.N)),by=unit_id]
weights <- as.data.table(x$lineage_weights/rowSums(x$lineage_weights))
setnames(weights,as.character(seq_len(ncol(weights))))
weights[,cell_index:=.I]
weights <- melt(weights,id.vars="cell_index",variable.name="leaf",value.name="weight")
weights[,leaf:=as.integer(as.character(leaf))]
weights <- merge(weights[weight>0],child_lineages,by="leaf",allow.cartesian=TRUE)
weights <- weights[,.(weight=sum(weight)),by=.(unit_id,child,cell_index)][weight>=CONFIG$min_weight]
cell_meta <- data.table(cell_index=seq_along(x$cells),node=x$lineage_clusters,cell_type=x$celltypes)
pre_nodes <- unit_table[,.(node=unit[[1L]]$segment_pre_nodes),by=unit_id]
pre_nodes <- merge(children[,.(unit_id,child)],pre_nodes,by="unit_id",allow.cartesian=TRUE)[,bin:=1L]
post_nodes <- unit_table[,.(child=unit[[1L]]$children,nodes=unit[[1L]]$segment_post_nodes),by=unit_id]
post_nodes <- post_nodes[,.(node=unlist(nodes)),by=.(unit_id,child)][,bin:=3L]
weights <- merge(weights,cell_meta,by="cell_index")
weights <- merge(weights,rbindlist(list(pre_nodes,post_nodes)),by=c("unit_id","child","node"))
parent_mass <- weights[bin==1L,.(parent_mass=sum(weight)),by=.(unit_id,cell_index)]
shared_pre <- merge(pre_nodes,cell_meta,by="node",allow.cartesian=TRUE)
shared_pre <- merge(shared_pre,parent_mass,by=c("unit_id","cell_index"))
shared_pre <- merge(shared_pre,children[,.(n_children=.N),by=unit_id],by="unit_id")
shared_pre[,weight:=parent_mass/n_children]
weights <- rbindlist(list(weights[bin==3L],shared_pre[,names(weights),with=FALSE]))
support <- weights[,.(mass=sum(weight),n_eff=sum(weight)^2/sum(weight^2)),
  by=.(unit_id,child,bin,cell_type)]
windows <- unique(weights[,.(unit_id,child,bin)])[,window_id:=.I]
weights <- merge(weights,windows,by=c("unit_id","child","bin"))
weights[,normalized_weight:=weight/sum(weight),by=.(window_id,cell_type)]
sender <- weights[cell_type==source_type]
receiver <- weights[cell_type==target_type]
source_weights <- sparseMatrix(i=sender$cell_index,j=sender$window_id,x=sender$normalized_weight,
  dims=c(n_cells,nrow(windows)))
target_weights <- sparseMatrix(i=receiver$cell_index,j=receiver$window_id,x=receiver$normalized_weight,
  dims=c(n_cells,nrow(windows)))
product_weights <- as.matrix(source_weights)[source_index,,drop=FALSE]*
  as.matrix(target_weights)[target_index,,drop=FALSE]
projected <- (kernel %*% product_weights)*native_weight
dimnames(projected) <- list(pair_keys,as.character(windows$window_id))
projected <- melt(as.data.table(projected,keep.rownames="pair_key"),id.vars="pair_key",
  variable.name="window_id",value.name="score")
projected[,window_id:=as.integer(as.character(window_id))]
projected <- merge(merge(projected,pairs,by="pair_key"),windows,by="window_id")
projected[,c("source","target"):=list(source_type,target_type)]
grid <- merge(projected,child_columns,by=c("unit_id","child"))
grid <- merge(grid,support[cell_type==source_type,
  .(unit_id,child,bin,source_n_eff=n_eff)],by=c("unit_id","child","bin"),all.x=TRUE)
grid <- merge(grid,support[cell_type==target_type,
  .(unit_id,child,bin,receiver_n_eff=n_eff,mass)],by=c("unit_id","child","bin"),all.x=TRUE)
grid[,c("w_eff","comp"):=list(pmin(source_n_eff,receiver_n_eff),fifelse(bin==1L,"pre","post"))]
grid[,row_valid:=is.finite(w_eff)&w_eff>=CONFIG$min_n_eff]
grid[row_valid%in%FALSE,score:=NA_real_]
}
binding[,scComm_grid_reused:=reused_sc]
fwrite(binding,file.path(out,"projection_binding.csv"))
fwrite(grid,file.path(out,"sccomm_projected_context_scores.csv"))
events <- grid[,getFromNamespace(".event_unit_contrasts","LineaLinkR")(
  .SD,units[[unit_id]],unit_id,"all",1L,activity_floor=0,
  activity_mode="continuous")[,!"unit_id"],by=unit_id]
keys <- c("unit_id","event_node","child","pair_idx","source","target")
metadata <- c(keys,"dataset","ligand","receptor","source_label","target_label",
  "sender_direction","role","true_positive","donor_confounded","truth_focal_branch",
  "truth_other_branch","focal_set","other_set","duplicate_emission","event_conflict",
  "unit_kind","hit_target_key","truth_row","mirror_row","resource_extra","null_row",
  "false_fork_row","prediction_key")
sc_detail <- merge(detail[sender_direction%in%TRUE,..metadata],events,by=keys,all.x=TRUE)
sc_detail[,c("score","score_mode"):=list(sibling_contrast,CONFIG$score_mode)]
sc_detail[is.na(evaluable),evaluable:=FALSE]
sc_detail[is.na(reason),reason:="no_native_common_window"]
all_predictions <- rbindlist(list(all_predictions,sc_detail),fill=TRUE)
fwrite(all_predictions,file.path(out,"projected_event_predictions.csv"))
n_methods <- uniqueN(all_predictions$score_mode)
ranked <- all_predictions[sender_direction%in%TRUE & duplicate_emission%in%c(FALSE,NA) &
  resource_extra%in%c(FALSE,NA) & evaluable%in%TRUE & is.finite(score)]
setorder(ranked,score_mode,prediction_key,-score,lineage,pair_idx)
ranked <- ranked[,.SD[1L],by=.(score_mode,prediction_key)]
joint <- ranked[,.(jointly_scored=uniqueN(score_mode)==n_methods),by=prediction_key]
ranked <- merge(ranked,joint,by="prediction_key")
ranked <- rbindlist(list(copy(ranked)[,evaluation_scope:="own_available"],
  copy(ranked[jointly_scored%in%TRUE])[,evaluation_scope:="common_candidates"]))
ranked[,c("method","task"):=list(score_mode,"projected_fork_event")]
setorder(ranked,method,evaluation_scope,-score,prediction_key)
ranked[,candidate_rank:=seq_len(.N),by=.(method,evaluation_scope)]
targets <- fread(file.path(base,"lrsim_target_recovery.csv"))[dataset==CONFIG$dataset,.(dataset,target_key,k_budget)]
n_targets <- nrow(targets)
metrics <- ranked[,.(n_candidates=.N,n_positive=sum(truth_row),n_null=sum(!truth_row),
  n_targets=n_targets,precision_at_k=sum(truth_row[candidate_rank<=n_targets])/n_targets,
  recall_at_k=uniqueN(hit_target_key[truth_row & candidate_rank<=n_targets])/n_targets),
  by=.(dataset,task,method,evaluation_scope)]
eligible <- metrics[n_positive>=2L & n_null>=2L,.(method,evaluation_scope)]
areas <- merge(ranked,eligible,by=c("method","evaluation_scope"))[,.(auprc=PRROC::pr.curve(
  scores.class0=score[truth_row],scores.class1=score[!truth_row],curve=FALSE)$auc.integral),
  by=.(method,evaluation_scope)]
metrics <- merge(metrics,areas,by=c("method","evaluation_scope"),all.x=TRUE)
metrics[n_targets==0L,c("precision_at_k","recall_at_k"):=list(NA_real_,NA_real_)]
metrics[,inference_status:="saved_native_scores_controlled_common_window_projection"]
previous_metrics <- fread(file.path(out,"metrics.csv"))[task!="projected_fork_event"]
fwrite(rbindlist(list(previous_metrics,metrics),fill=TRUE),file.path(out,"metrics.csv"))
rank_fields <- c("dataset","task","method","evaluation_scope","prediction_key","truth_row","score","candidate_rank")
previous_ranked <- fread(file.path(out,"ranked_predictions.csv"))[task!="projected_fork_event"]
fwrite(rbindlist(list(previous_ranked,ranked[,..rank_fields]),fill=TRUE),file.path(out,"ranked_predictions.csv"))
correct <- ranked[truth_row%in%TRUE,.(best_correct_rank=sort(candidate_rank)[1L],
  local_contrast=sibling_contrast[which.min(candidate_rank)],
  ranking_score=score[which.min(candidate_rank)]),
  by=.(method,evaluation_scope,target_key=hit_target_key)]
recovery_methods <- CJ(method=unique(all_predictions$score_mode),
  evaluation_scope=c("own_available","common_candidates"))[,join_key:=1L]
planned <- merge(targets[,join_key:=1L],recovery_methods,
  by="join_key",allow.cartesian=TRUE)[,join_key:=NULL]
recovery <- merge(planned,correct,by=c("method","evaluation_scope","target_key"),all.x=TRUE)
recovery[,c("correct_candidate_available","correct_local_direction","hit_in_budget"):=list(
  is.finite(best_correct_rank),is.finite(local_contrast)&local_contrast>0,
  is.finite(best_correct_rank)&best_correct_rank<=k_budget)]
fwrite(recovery,file.path(out,"projected_target_recovery.csv"))
metric_files <- list.files(file.path(base,"dynamic_methods"),pattern="^metrics[.]csv$",
  recursive=TRUE,full.names=TRUE)
all_metrics <- rbindlist(lapply(metric_files,fread),fill=TRUE)
scope <- fread(file.path(base,"lrsim_metrics.csv"))[,.(dataset,n_targets_total=n_targets)]
method_scope <- unique(all_metrics[,.(task,method,evaluation_scope)])
complete <- merge(scope[,join_key:=1L],method_scope[,join_key:=1L],by="join_key",
  allow.cartesian=TRUE)[,join_key:=NULL]
all_metrics <- merge(complete,all_metrics,by=c("dataset","task","method","evaluation_scope"),all.x=TRUE)
all_metrics[is.na(n_targets),n_targets:=n_targets_total]
all_metrics[,completed_projection:=dataset%chin%basename(dirname(metric_files))]
all_metrics[completed_projection & is.na(n_candidates),
  c("n_candidates","n_positive","n_null"):=list(0L,0L,0L)]
all_metrics[,evaluation_status:=fcase(!completed_projection,"missing_saved_projection",
  n_candidates==0L & evaluation_scope=="common_candidates","no_jointly_evaluable_candidates",
  n_candidates==0L,"no_evaluable_native_output",
  n_positive==0L,"no_positive_targets_in_scored_candidates",default="scored")]
all_metrics[is.na(inference_status),inference_status:=fifelse(completed_projection,
  "saved_native_output_unavailable_on_scope","missing_native_output")]
all_metrics[,c("n_targets_total","completed_projection"):=NULL]
fwrite(all_metrics,file.path(base,"dynamic_methods.csv"))
cat("NATIVE_KERNEL_MAX_ERROR",max(equivalence$absolute_error),"\n")
print(metrics[evaluation_scope=="common_candidates",.(method,n_candidates,n_positive,n_targets,auprc,precision_at_k,recall_at_k)])
print(recovery[evaluation_scope=="own_available",.(n_targets=.N,
  n_available=sum(correct_candidate_available),n_direction=sum(correct_local_direction),
  n_in_budget=sum(hit_in_budget)),by=method])
