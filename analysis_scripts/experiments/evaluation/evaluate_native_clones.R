CONFIG <- list(paths=list(root=normalizePath(commandArgs(TRUE)[1L],winslash="/",mustWork=TRUE),
  model="cache/final065/larry_invitro/linealinkr_complete.rds",
  metadata="inputs/external_validation/LARRY/stateFate_inVitro_metadata.txt.gz",
  clone_matrix="inputs/external_validation/LARRY/stateFate_inVitro_clone_matrix.mtx.gz",
  native_probabilities="results/evaluation/final065/saved_native/native_CellRank2_larry_probabilities.csv.gz",
  baseline_scores="results/evaluation/calendar_lineage_validation/larry_invitro_all_clone_composition_scores.csv",
  output="results/evaluation/final065/saved_native"),
  threads=parallel::detectCores(),terminal_day=6,epsilon=1e-12)
Sys.setenv(OMP_NUM_THREADS=CONFIG$threads,OPENBLAS_NUM_THREADS=CONFIG$threads,MKL_NUM_THREADS=CONFIG$threads)
library(data.table)
library(Matrix)
setDTthreads(CONFIG$threads)
root <- CONFIG$paths$root
out <- file.path(root,CONFIG$paths$output)
x <- readRDS(file.path(root,CONFIG$paths$model))
meta <- fread(file.path(root,CONFIG$paths$metadata),encoding="UTF-8")
meta[,cell_id:=paste(Library,gsub("-","",get("Cell barcode"),fixed=TRUE),sep=":")]
meta[,raw_row:=seq_len(.N)]
meta[,context:=as.character(get("Starting population"))]
meta[,day:=as.numeric(get("Time point"))]
meta[,fate:=as.character(get("Cell type annotation"))]
clones <- as(readMM(gzfile(file.path(root,CONFIG$paths$clone_matrix))),"CsparseMatrix")
stopifnot(nrow(clones)==nrow(meta),all(Matrix::rowSums(clones)%in%c(0,1)),!anyDuplicated(meta$cell_id))
members <- as.data.table(summary(clones))
setnames(members,c("i","j","x"),c("raw_row","clone","membership"))
trace <- merge(members,meta[,.(raw_row,context,day,fate)],by="raw_row")
labels <- sort(unique(meta$fate))
contexts <- sort(unique(meta$context))
future <- trace[day==CONFIG$terminal_day,.(n=sum(membership)),by=.(context,clone,fate)]
future[,context:=factor(context,levels=contexts)]
future[,fate:=factor(fate,levels=labels)]
truth <- dcast(future,context+clone~fate,value.var="n",fill=0,drop=FALSE)
truth[,context:=as.character(context)]
truth[,n_future_cells:=rowSums(.SD),.SDcols=labels]
truth[,labels] <- as.data.frame(as.matrix(truth[,labels,with=FALSE])/pmax(truth$n_future_cells,1))
fitted <- meta[match(x$cells,cell_id),.(cell_id,raw_row,context,day,fate)]
fitted[,index:=seq_len(.N)]
stopifnot(identical(fitted$cell_id,x$cells),!anyNA(fitted$raw_row))
n <- fread(file.path(root,CONFIG$paths$native_probabilities),encoding="UTF-8",header=TRUE)
stopifnot(setequal(n$cell_id,x$cells),!anyDuplicated(n$cell_id))
n <- n[match(x$cells,cell_id)]
native <- as.matrix(n[,!"cell_id"])
reference <- as.matrix(switch(as.character(is.null(x$lineage_reference_weights)),
  "TRUE"=x$lineage_weights,"FALSE"=x$lineage_reference_weights))
raw_mass_error <- data.table(method=c("LLR_0.6.5_reference","CellRank2"),
  n_saved_cells=c(nrow(reference),nrow(native)),
  maximum_raw_row_mass_error=c(max(abs(rowSums(reference)-1)),max(abs(rowSums(native)-1))))
reference <- reference/rowSums(reference)
native <- native/rowSums(native)
stopifnot(all(is.finite(native)),all(native>=0),all(is.finite(reference)))
onehot <- outer(fitted$fate,labels,"==")*1
colnames(onehot) <- labels
queries <- fitted[day<CONFIG$terminal_day]
terminal_contexts <- unique(fitted[day==CONFIG$terminal_day,context])
design <- CJ(context=terminal_contexts,method=c("LLR_0.6.5_reference","CellRank2"))
predictions <- design[,{
  w <- switch(.BY$method,LLR_0.6.5_reference=reference,CellRank2=native)
  late <- which(fitted$context==.BY$context&fitted$day==CONFIG$terminal_day)
  early <- which(fitted$context==.BY$context&fitted$day<CONFIG$terminal_day)
  mass <- colSums(w[late,,drop=FALSE])
  stopifnot(all(mass>0))
  decoder <- crossprod(w[late,,drop=FALSE],onehot[late,,drop=FALSE])/mass
  colnames(decoder) <- labels
  value <- as.data.table(w[early,,drop=FALSE]%*%decoder)
  cbind(fitted[early,.(cell_id,raw_row,day)],value)
},by=.(context,method)]
plan <- merge(predictions,members[,.(raw_row,clone)],by="raw_row",all.x=TRUE)
coverage <- plan[,.(n_planned_cells=.N,n_traced_cells=sum(!is.na(clone)),
  n_finite_cells=sum(is.finite(get(labels[1L])))),by=.(context,day,method)]
traced <- plan[!is.na(clone)]
composition <- traced[,as.list(colMeans(as.matrix(.SD))),by=.(method,context,day,clone),.SDcols=labels]
composition <- merge(composition,traced[,.(n_current_cells=.N),by=.(method,context,day,clone)],
  by=c("method","context","day","clone"))
composition <- merge(composition,truth,by=c("context","clone"),all.x=TRUE,suffixes=c("","_truth"))
p <- as.matrix(composition[,labels,with=FALSE])
q <- as.matrix(composition[,paste0(labels,"_truth"),with=FALSE])
composition[,truth_available:=is.finite(n_future_cells)&n_future_cells>0]
composition[,prediction_available:=rowSums(is.finite(p))==length(labels)&rowSums(p)>0]
composition[,soft_brier:=rowSums((p-q)^2)]
composition[,descendant_brier:=rowSums(p^2)-2*rowSums(p*q)+1]
composition[,log_loss:=-rowSums(q*log(pmax(p,CONFIG$epsilon)))]
composition[truth_available%in%FALSE,c("soft_brier","descendant_brier","log_loss"):=list(NA_real_,NA_real_,NA_real_)]
summary <- composition[,.(n_planned_clones=.N,n_truth_available=sum(truth_available),
  n_scored=sum(truth_available&prediction_available),
  soft_brier=mean(soft_brier[truth_available&prediction_available]),
  descendant_brier=mean(descendant_brier[truth_available&prediction_available]),
  log_loss=mean(log_loss[truth_available&prediction_available])),by=.(method,context,day)]
summary[,dataset:="larry_invitro"]
summary <- merge(summary,raw_mass_error,by="method")
summary[,c("task","unit","reference_protocol","native_probability_policy","fit_status","ground_truth_type"):=list(
  "clone_fate_composition","clone_within_context_and_current_day",
  "same_context_day6_RNA_labels_weighted_terminal_average;no_clone_ID_in_decoder",
  "same_row_sum_normalization_as_existing_LLR_evaluation;raw_masses_saved",
  "read_saved_native_fates;no_competitor_refit",
  "observed_future_clonal_composition;not_individual_cell_fate")]
old <- fread(file.path(root,CONFIG$paths$baseline_scores))
baseline <- merge(old[model=="native_reference",.(context,day,clone,old_soft_brier=soft_brier,old_log_loss=log_loss)],
  composition[method=="LLR_0.6.5_reference",.(context,day,clone,soft_brier,log_loss)],
  by=c("context","day","clone"),all=TRUE)
baseline_check <- baseline[,.(n_old=sum(is.finite(old_soft_brier)),n_current=sum(is.finite(soft_brier)),
  maximum_brier_difference=max(abs(old_soft_brier-soft_brier),na.rm=TRUE),
  maximum_log_loss_difference=max(abs(old_log_loss-log_loss),na.rm=TRUE))]
fwrite(summary,file.path(out,"classic_native_clonal_metrics.csv"))
fwrite(coverage,file.path(out,"native_clonal_input_coverage.csv"))
fwrite(composition,file.path(out,"native_clonal_composition_scores.csv.gz"))
fwrite(baseline_check,file.path(out,"native_clonal_existing_protocol_readback.csv"))
print(baseline_check)
print(summary[,.(method,context,day,n_planned_clones,n_truth_available,n_scored,soft_brier,log_loss)])
