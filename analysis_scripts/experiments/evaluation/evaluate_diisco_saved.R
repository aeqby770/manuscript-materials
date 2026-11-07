CONFIG <- list(root=normalizePath(Sys.getenv("LLR_WORK_ROOT",getwd()),winslash="/",mustWork=TRUE),
  threads=parallel::detectCores(),native="results/by_algorithm/DIISCO",
  output="results/evaluation/final065/diisco")
Sys.setenv(OMP_NUM_THREADS=CONFIG$threads,OPENBLAS_NUM_THREADS=CONFIG$threads,MKL_NUM_THREADS=CONFIG$threads)
library(data.table)
setDTthreads(CONFIG$threads)
paths <- list.files(file.path(CONFIG$root,CONFIG$native),pattern="^coverage[.]csv$",
  recursive=TRUE,full.names=TRUE)
paths <- paths[!grepl("[/\\\\]failed_",paths)]
coverage <- rbindlist(lapply(paths,fread),idcol="file_id",fill=TRUE)
coverage[,coverage_path:=paths[as.integer(file_id)]]
coverage[,unit_dir:=dirname(coverage_path)]
fitted <- coverage[grepl("^scored",result_status),{
  network <- fread(file.path(unit_dir,"interaction_posterior_mean.csv.gz"))
  latent <- fread(file.path(unit_dir,"latent_proportions.csv"))
  observed <- fread(file.path(unit_dir,"input_sample_proportions.csv"))
  setnames(latent,1L,"sample")
  setnames(observed,1L,"sample")
  network[,sample:=as.character(sample)]
  latent[,sample:=as.character(sample)]
  observed[,c("sample","time","n_cells"):=list(as.character(sample),as.numeric(time),as.integer(n_cells))]
  state <- melt(latent,id.vars="sample",variable.name="sender",value.name="latent_f")
  product <- merge(network,state,by=c("sample","sender"))
  # The authors' Normal emission mean is W f. The saved F_mean alone is
  # a latent state; its historical filename does not make it a proportion.
  prediction <- product[,.(fitted_mean=sum(posterior_mean*latent_f)),by=.(sample,receiver)]
  truth <- melt(observed,id.vars=c("sample","time","n_cells"),
    variable.name="receiver",value.name="observed_proportion")
  truth[,observed_proportion:=as.numeric(observed_proportion)]
  matched <- merge(prediction,truth,by=c("sample","receiver"))
  stopifnot(nrow(matched)==nrow(truth),all(is.finite(matched$fitted_mean)))
  matched[,constant_mean:=mean(observed_proportion),by=receiver]
  matched
},by=.(dataset,unit_dir)]
summary <- fitted[,.(n_samples=uniqueN(sample),n_cell_types=uniqueN(receiver),
  n_timepoints=uniqueN(time),n_rows=.N,
  fitted_rmse=sqrt(mean((fitted_mean-observed_proportion)^2)),
  fitted_mae=mean(abs(fitted_mean-observed_proportion)),
  constant_mean_rmse=sqrt(mean((constant_mean-observed_proportion)^2)),
  r2_over_cell_type_mean=1-sum((fitted_mean-observed_proportion)^2)/
    sum((constant_mean-observed_proportion)^2),
  n_fitted_means_outside_unit_interval=sum(fitted_mean<0|fitted_mean>1)),
  by=.(dataset,unit_dir)]
summary[,c("task","interpretation"):=list("observed_cell_type_proportion_reconstruction",
  "in_sample_factorized_variational_mean_W_times_f; not_LR_or_branch_accuracy")]
out <- file.path(CONFIG$root,CONFIG$output)
dir.create(out,recursive=TRUE,showWarnings=FALSE)
fwrite(fitted,file.path(out,"native_observation_reconstruction.csv"))
fwrite(summary,file.path(out,"native_reconstruction_metrics.csv"))
print(summary[,!"unit_dir"])
