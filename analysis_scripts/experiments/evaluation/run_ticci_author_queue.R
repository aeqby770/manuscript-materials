CONFIG <- list(paths=list(root=normalizePath(commandArgs(TRUE)[1L],winslash="/",mustWork=TRUE),
 manifest="experiments/evaluation/ticci_author_protocol_inputs.csv",
 output="results/by_algorithm/TICCI_author_functional_entropy",
 python="E:/Python/biotrade/Scripts/python.exe",R="D:/R/R-4.5.3/bin/Rscript.exe"),
 threads=as.integer(Sys.getenv("LLR_NATIVE_THREADS",parallel::detectCores())),plot=list())
Sys.setenv(OMP_NUM_THREADS=CONFIG$threads,OPENBLAS_NUM_THREADS=CONFIG$threads,
 MKL_NUM_THREADS=CONFIG$threads,PYTHONIOENCODING="utf-8",
 LLR_TICCI_ROOT=file.path(CONFIG$paths$root,"inputs/comparison/external_methods/TICCI"),LLR_WORK_ROOT=CONFIG$paths$root)
.libPaths(c(file.path(CONFIG$paths$root,"library"),.libPaths()))
library(data.table)
setDTthreads(CONFIG$threads)
root <- CONFIG$paths$root
jobs <- fread(file.path(root,CONFIG$paths$manifest))
scope <- commandArgs(TRUE)[2L]
jobs <- switch(scope,
 simulations=jobs[entropy_mode=="expression_gene_groups"],
 dentategyrus=jobs[dataset=="dentategyrus"],
 mouse_remaining=jobs[species=="mouse"&dataset!="dentategyrus"],
 human=jobs[species=="human"],
 jobs[dataset==scope])
jobs<-jobs[!file.exists(file.path(root,CONFIG$paths$output,dataset,"ticci_author_result.h5ad"))|
 !file.exists(file.path(root,CONFIG$paths$output,dataset,"ticci_pseudotime.csv"))]
results <- jobs[,{
 out <- file.path(root,CONFIG$paths$output,dataset)
 dir.create(out,recursive=TRUE,showWarnings=FALSE)
 previous<-fread(file.path(out,"run_status.csv"))
 started<-as.POSIXct(gsub("T"," ",as.character(previous$started[1L]),fixed=TRUE),format="%Y-%m-%d %H:%M:%OS",tz="Asia/Shanghai")
 input_file<-fifelse(grepl("^[A-Za-z]:",input_h5ad),input_h5ad,file.path(root,input_h5ad))
 fwrite(data.table(dataset,species,input_h5ad,status="preparing_full_native_input",
  started=as.character(started)),file.path(out,"run_status.csv"))
 prepare_exit<-switch(as.character(file.exists(file.path(out,"input_binding.csv"))),
 "TRUE"=0L,"FALSE"=system2(CONFIG$paths$python,c("-X","utf8",
  shQuote(file.path(root,"experiments/evaluation/prepare_ticci_author_input.py")),shQuote(root),shQuote(dataset)),
  stdout=file.path(out,"prepare.log"),stderr=file.path(out,"prepare.err")))
 entropy_exit<-switch(as.character(prepare_exit==0L),"TRUE"={
  switch(as.character(file.exists(file.path(out,"ticci_functional_entropy.csv"))),
   "TRUE"=0L,"FALSE"=system2(CONFIG$paths$R,c("--encoding=UTF-8",
    shQuote(file.path(root,"experiments/evaluation/run_ticci_author_entropy.R")),
    shQuote(root),shQuote(dataset)),stdout=file.path(out,"entropy.log"),stderr=file.path(out,"entropy.err")))
 },"FALSE"=NA_integer_)
 cellchat_exit<-switch(as.character(prepare_exit==0L&entropy_exit%in%0L),"TRUE"={
  switch(cellchat_mode,
   empty_molecular_resource={
    fwrite(data.table(source=character(),target=character(),interaction_name_2=character(),prob=numeric()),
     file.path(out,"ticci_cellchat_interactions.csv"))
    0L
   },
   reuse_saved_calls=0L,
   new_author_CellChat=switch(as.character(all(file.exists(file.path(out,c("ticci_cellchat_interactions.csv","ticci_cellchat.rds"))))),
    "TRUE"=0L,"FALSE"=system2(CONFIG$paths$R,c("--encoding=UTF-8",
     shQuote(file.path(root,"experiments/evaluation/run_ticci_author_cellchat.R")),
     shQuote(root),shQuote(dataset)),stdout=file.path(out,"cellchat.log"),stderr=file.path(out,"cellchat.err"))))
 },"FALSE"=NA_integer_)
 trajectory_exit <- switch(as.character(prepare_exit==0L&entropy_exit%in%0L&cellchat_exit%in%0L),"TRUE"={
  switch(cellchat_mode,new_author_CellChat={
   fwrite(data.table(dataset,species,input_h5ad,status="binding_saved_CellChat_to_full_input",
    started=as.character(started)),file.path(out,"run_status.csv"))
   library(CellChat)
   library(Matrix)
   native_cellchat<-readRDS(file.path(out,"ticci_cellchat.rds"))
   full_cells<-readLines(file.path(out,"entropy_cells.txt"),encoding="UTF-8")
   original_genes<-readLines(file.path(out,"entropy_genes.txt"),encoding="UTF-8")
   gene_binding<-fread(file.path(out,"ticci_gene_id_binding.csv"))
   input_groups<-fread(file.path(out,"ticci_cellchat_groups.csv"),colClasses="character")[match(full_cells,cell_id)]
   input_shape<-rhdf5::h5read(file.path(out,"entropy_input.h5"),"shape")
   original_matrix<-new("dgCMatrix",
    i=as.integer(rhdf5::h5read(file.path(out,"entropy_input.h5"),"indices")),
    p=as.integer(rhdf5::h5read(file.path(out,"entropy_input.h5"),"indptr")),
    x=as.numeric(rhdf5::h5read(file.path(out,"entropy_input.h5"),"data")),Dim=as.integer(rev(input_shape)))
   original_matrix<-original_matrix%*%Diagonal(x=10000/colSums(original_matrix))
   original_matrix@x<-log1p(original_matrix@x)
   stopifnot(identical(dim(original_matrix),dim(native_cellchat@data)),
    identical(original_genes,gene_binding$original_gene_id),
    identical(rownames(native_cellchat@data),gene_binding$resource_symbol),
    identical(input_groups$cell_id,full_cells),
    identical(as.character(native_cellchat@idents),input_groups$Classification),
    identical(as.character(native_cellchat@meta$Classification),input_groups$Classification))
   matrix_error<-max(c(0,abs((native_cellchat@data-original_matrix)@x)))
   stopifnot(matrix_error==0)
   native_names<-as.character(colnames(native_cellchat@data))
   stopifnot(length(native_names)%in%c(0L,length(full_cells)))
   native_names<-c(native_names,rep(NA_character_,length(full_cells)-length(native_names)))
   fwrite(data.table(native_column_index=seq_along(full_cells),native_cell_id=native_names,
    native_metadata_row_id=rownames(native_cellchat@meta),original_cell_id=full_cells,
    Classification=input_groups$Classification),file.path(out,"ticci_cellchat_cell_id_binding.csv"))
   fwrite(data.table(dataset,n_cells=length(full_cells),n_original_gene_rows=length(original_genes),
    n_native_cell_names=sum(!is.na(native_names)),native_cell_names_equal_original=identical(native_names,full_cells),
    max_absolute_expression_error=matrix_error,all_population_labels_exact=TRUE,
    binding_mode="unchanged_native_column_order;complete_input_matrix_and_population_labels_exact",
    native_result_modified=FALSE),file.path(out,"ticci_cellchat_binding_summary.csv"))
   rm(native_cellchat,original_matrix,gene_binding,input_groups,native_names);gc()
  },reuse_saved_calls=invisible(NULL),empty_molecular_resource=invisible(NULL))
  fwrite(data.table(dataset,species,input_h5ad,status="author_lineage_and_DPT",
   started=as.character(started)),file.path(out,"run_status.csv"))
  system2(CONFIG$paths$python,c("-X","utf8",
   shQuote(file.path(root,"experiments/evaluation/run_ticci_author_trajectory.py")),
   shQuote(input_file),shQuote(out)),
   stdout=file.path(out,"trajectory.log"),stderr=file.path(out,"trajectory.err"))
 },"FALSE"=NA_integer_)
 status <- fcase(prepare_exit!=0L,"input_binding_or_preparation_failed",
  is.na(entropy_exit)|entropy_exit!=0L,"functional_entropy_failed",
  is.na(cellchat_exit)|cellchat_exit!=0L,"author_CellChat_failed",
  is.na(trajectory_exit)|trajectory_exit!=0L,"lineage_or_DPT_failed",
  default="completed_author_functional_entropy_fixed_graph_protocol")
 progress <- data.table(dataset,species,input_h5ad,status,
  prepare_exit,entropy_exit,cellchat_exit,trajectory_exit,started=as.character(started),
  finished=as.character(Sys.time()),seconds=as.numeric(difftime(Sys.time(),started,units="secs")))
 fwrite(progress,file.path(out,"run_status.csv"))
 print(progress)
 stopifnot(status=="completed_author_functional_entropy_fixed_graph_protocol")
 progress[,!"dataset"]
},by=dataset]
fwrite(results,file.path(root,CONFIG$paths$output,paste0("queue_",scope,".csv")))
stopifnot(all(results$status=="completed_author_functional_entropy_fixed_graph_protocol"))
