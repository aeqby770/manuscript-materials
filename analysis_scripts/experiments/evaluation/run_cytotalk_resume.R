CONFIG <- list(root=Sys.getenv("LLR_WORK_ROOT"),threads=as.integer(Sys.getenv("LLR_NATIVE_THREADS","4")),seed=12345L,dataset="paul15",ntrial=1000L)
Sys.setenv(OMP_NUM_THREADS=CONFIG$threads,OPENBLAS_NUM_THREADS=CONFIG$threads,MKL_NUM_THREADS=CONFIG$threads)
library(data.table)
library(Matrix)
library(Rcpp)
library(infotheo)
library(CytoTalk)
data(list=c("pcg_mouse","lrp_mouse"),package="CytoTalk")
setDTthreads(CONFIG$threads)
root <- CONFIG$root
reticulate::use_python("E:/Python/biotrade/Scripts/python.exe",required=TRUE)
reticulate::py_run_string(paste0("import sys;sys.path.insert(0,",jsonlite::toJSON(file.path(root,"tools/competitors/dependencies/CytoTalk"),auto_unbox=TRUE),");import pcst_fast;print(pcst_fast.__file__)"))
work <- file.path(root,"results/by_algorithm/CytoTalk",CONFIG$dataset)
dir.create(work,recursive=TRUE,showWarnings=FALSE)
cache <- file.path(work,"native_type_cache")
dir.create(cache,recursive=TRUE,showWarnings=FALSE)
proof <- fread(file.path(root,"results/evaluation/final065/saved_native/cytotalk_mi_equivalence.csv"))
stopifnot(proof$identical_MI,proof$identical_ARACNE,proof$max_abs_error==0)
sourceCpp(file.path(root,"experiments/evaluation/cytotalk_mi_native.cpp"),cacheDir=file.path(root,"build/cytotalk_cpp"))
cellchat <- readRDS(file.path(root,"inputs/comparison/recovered_native/by_algorithm/CellChat/paul15/cellchat.rds"))
meta <- fread(file.path(root,"inputs/comparison/recovered_native/input_metadata/paul15/meta.csv"))
stopifnot(identical(colnames(cellchat@data),meta$cell_id),identical(as.character(cellchat@meta$labels),meta$paul15_clusters))
lst <- list(mat=cellchat@data,cell_types=as.character(cellchat@meta$labels))
mat_pem_shared <- CytoTalk::pem(lst)
old <- readRDS(file.path(root,"inputs/comparison/recovered_native/by_algorithm/CytoTalk/paul15/X13Baso__X10GMP/cytotalk_result.rds"))
stopifnot(identical(mat_pem_shared,old$pem))
rm(cellchat,old)
all_types <- names(sort(table(lst$cell_types),decreasing=TRUE))
pairs <- as.data.table(t(combn(all_types,2)))
setnames(pairs,c("population_a","population_b"))
pairs[,pair_order:=.I]
pairs[,pair_key:=paste(pmin(make.names(population_a),make.names(population_b)),pmax(make.names(population_a),make.names(population_b)),sep="__")]
plan <- fread(file.path(root,"results/evaluation/final065/saved_native/native_CytoTalk_pair_coverage.csv"))
missing <- plan[dataset==CONFIG$dataset & planned_pair & native_model_status!="valid_native_network",.(pair_key)]
queue <- merge(pairs,missing,by="pair_key",sort=FALSE)
setorder(queue,pair_order)
native_base <- file.path(root,"inputs/comparison/recovered_native/by_algorithm/CytoTalk",CONFIG$dataset)
queue[,pair_dir:=file.path(native_base,paste(make.names(population_a),make.names(population_b),sep="__"))]
queue <- queue[is.na(file.info(file.path(pair_dir,"cytotalk_result.rds"))$size)]
fwrite(data.table(dataset=CONFIG$dataset,n_input_cells=ncol(lst$mat),n_input_genes=nrow(lst$mat),planned_pairs=nrow(pairs),
 saved_pairs_before=sum(plan[dataset==CONFIG$dataset]$saved_network),missing_pairs=nrow(queue),seed=CONFIG$seed,
 ntrial=CONFIG$ntrial,input_source="saved_CellChat_full_expression_exact_PEM_match",MI_engine="same_native_infotheo_entropyR_batched_Cpp",
 random_nonselftalk="original_per_pair_calls_preserved",
 relevance_boundary="zero_degree_inverse0;positive_degree_unchanged;all_nodes_retained",
 relevance_equivalence="7MEP524_output_bitwise;18Eos1604_positive_block_bitwise;full1615_finite"),file.path(work,"input_binding.csv"))
types <- data.table(cell_type=unique(c(queue$population_a,queue$population_b)))
types[,cache_path:=file.path(cache,paste0(make.names(cell_type),".rds"))]
cache_summary <- types[,{
 started <- proc.time()[["elapsed"]]
 switch(as.character(file.exists(cache_path)),
  "FALSE"={
   mat_type <- CytoTalk::extract_group(cell_type,lst)
   filtered <- CytoTalk:::subset_rownames(CytoTalk:::subset_non_zero_old(mat_type,0.2),pcg_mouse)
   discrete <- CytoTalk::discretize_sparse(t(filtered))
   mi <- cytotalk_mi_native(as.matrix(discrete),getNativeSymbolInfo("entropyR",PACKAGE="infotheo")$address)
   dimnames(mi) <- list(colnames(discrete),colnames(discrete))
   intra <- Matrix(parmigene::aracne.m(CytoTalk:::zero_diag(mi)))
   saveRDS(intra,cache_path)
   cat("TYPE_CACHE",cell_type,ncol(mat_type),nrow(intra),"seconds",proc.time()[["elapsed"]]-started,"\n")
   flush.console()
  })
 list(seconds=proc.time()[["elapsed"]]-started)
},by=.(cell_type,cache_path)]
switch(as.character(file.exists(file.path(work,"native_type_cache_summary.csv"))),
 "FALSE"=fwrite(cache_summary,file.path(work,"native_type_cache_summary.csv")))
# Official master/feature_RcallPy/pipeline retain mpower(D,-0.5) at degree zero.
# The diagonal generalized inverse keeps every gene; positive-degree arithmetic is unchanged.
# Verified: full 524-gene 7MEP output and positive 1604-gene 18Eos block are bitwise identical.
gene_relevance_boundary <- CytoTalk::gene_relevance
relevance_statements <- as.list(body(gene_relevance_boundary))
relevance_index <- which(vapply(relevance_statements,deparse1,character(1L))=="mat_nsq <- corpcor::mpower(mat_nsq, -0.5)")
stopifnot(length(relevance_index)==1L)
relevance_statements[[relevance_index]] <- quote({
 degree <- diag(mat_nsq)
 mat_nsq <- diag(ifelse(degree>0,degree^(-0.5),0),nrow=length(degree))
})
body(gene_relevance_boundary) <- as.call(relevance_statements)
types[,relevance_path:=file.path(cache,paste0(make.names(cell_type),"_native_relevance.rds"))]
relevance_summary <- types[,{
 started<-proc.time()[["elapsed"]]
 switch(as.character(file.exists(relevance_path)),"FALSE"={
  relevance<-gene_relevance_boundary(readRDS(cache_path),lrp_mouse)
  saveRDS(relevance,relevance_path)
  stopifnot(identical(relevance,readRDS(relevance_path)))
  cat("NATIVE_RELEVANCE_CACHE",cell_type,proc.time()[["elapsed"]]-started,"seconds\n");flush.console()
 })
 relevance_before <- readRDS(relevance_path)
 switch(as.character(any(!is.finite(relevance_before))),"TRUE"={
  relevance <- gene_relevance_boundary(readRDS(cache_path),lrp_mouse)
  stopifnot(all(is.finite(relevance)),identical(names(relevance),names(relevance_before)))
  saveRDS(relevance,relevance_path)
  cat("ZERO_DEGREE_BOUNDARY_REPAIRED",cell_type,sum(!is.finite(relevance_before)),"to",sum(!is.finite(relevance)),"\n");flush.console()
 })
 list(seconds=proc.time()[["elapsed"]]-started,n_nonfinite_before=sum(!is.finite(relevance_before)))
},by=.(cell_type,cache_path,relevance_path)]
switch(as.character(file.exists(file.path(work,"native_relevance_cache_summary.csv"))),
 "FALSE"=fwrite(relevance_summary,file.path(work,"native_relevance_cache_summary.csv")))
type_validity<-types[,{
 intra<-readRDS(cache_path)
 relevance<-readRDS(relevance_path)
 zero<-Matrix::rowSums(intra)==0
 list(n_cells=sum(lst$cell_types==cell_type),n_genes=nrow(intra),n_zero_degree=sum(zero),
  zero_degree_gene_ids=paste(rownames(intra)[zero],collapse=";"),
  n_nonfinite_relevance=sum(!is.finite(relevance)))
},by=cell_type]
fwrite(type_validity,file.path(work,"native_type_validity.csv"))
stopifnot(all(type_validity$n_nonfinite_relevance==0L))
failures <- fread(file.path(work,"cytotalk_pair_failures.csv"))
cat("VALID_MISSING_PAIRS",nrow(queue),"BOUNDARY_HISTORY_PAIRS",nrow(failures),"\n");flush.console()
native_runner <- CytoTalk::run_cytotalk
env <- new.env(parent=environment(native_runner))
env$mat_pem_shared <- mat_pem_shared
env$cache_paths <- setNames(types$cache_path,types$cell_type)
env$relevance_paths <- setNames(types$relevance_path,types$cell_type)
native_integrate <- CytoTalk::integrate_network
environment(native_integrate)<-env
integrate_statements<-as.list(body(native_integrate))
integrate_text<-vapply(integrate_statements,deparse1,character(1L))
integrate_statements[[which(startsWith(integrate_text,"vec_gr_a <-"))]]<-quote(vec_gr_a<-readRDS(relevance_paths[[cell_type_a]]))
integrate_statements[[which(startsWith(integrate_text,"vec_gr_b <-"))]]<-quote(vec_gr_b<-readRDS(relevance_paths[[cell_type_b]]))
body(native_integrate)<-as.call(integrate_statements)
env$integrate_network<-native_integrate
environment(native_runner) <- env
native_graph <- CytoTalk::graph_pathway
graph_statements<-as.list(body(native_graph))
graph_text<-vapply(graph_statements,deparse1,character(1L))
graph_statements[[which(startsWith(graph_text,"color <-"))]]<-quote(color<-grDevices::hsv(
 ifelse(ew(df_node$node),0.02,0.55),
 pmin(1,pmax(0,ifelse(is.finite(df_node$pem),df_node$pem,0.5))),1))
body(native_graph)<-as.call(graph_statements)
env$graph_pathway<-native_graph
statements <- as.list(body(native_runner))
text <- vapply(statements,deparse1,character(1L))
statements[[which(startsWith(text,"lst_pcst <-"))]]<-quote({
 stopifnot(all(is.finite(lst_net$nodes$prize)),all(is.finite(lst_net$edges$cost)))
 # Existing partial PCST files for failed Eos pairs used undefined prizes.
 # Recompute these affected pairs; completed valid RDS pairs never enter this queue.
 lst_pcst<-run_pcst(lst_net,beta_max,omega_min,omega_max)})
statements[[which(startsWith(text,"mat_pem <-"))]] <- quote(mat_pem<-mat_pem_shared)
idx <- which(startsWith(text,"mat_disc_")|startsWith(text,"mat_mi_")|startsWith(text,"dimnames(mat_mi_"))
statements[idx] <- rep(list(quote(NULL)),length(idx))
statements[[which(startsWith(text,"mat_intra_a <-"))]] <- quote(mat_intra_a<-readRDS(cache_paths[[cell_type_a]]))
statements[[which(startsWith(text,"mat_intra_b <-"))]] <- quote(mat_intra_b<-readRDS(cache_paths[[cell_type_b]]))
body(native_runner) <- as.call(statements)
set.seed(CONFIG$seed)
switch(as.character(file.exists(file.path(work,"native_missing_pair_RNG.rds"))),
 "TRUE"={.Random.seed<-readRDS(file.path(work,"native_missing_pair_RNG.rds"))})
progress <- file.path(work,"cytotalk_summary.csv")
switch(as.character(file.exists(progress)),"FALSE"=fwrite(data.table(dataset=character(),cell_type_a=character(),cell_type_b=character(),pair_key=character(),
 n_input_cells=integer(),n_pair_cells=integer(),n_native_edges=integer(),ntrial=integer(),n_native_pathways=integer(),n_finite_pathway_pvalues=integer(),seconds=numeric(),status=character()),progress))
result_summary <- queue[,{
 dir.create(pair_dir,recursive=TRUE,showWarnings=FALSE)
 started <- proc.time()[["elapsed"]]
 result <- native_runner(lst_scrna=lst,cell_type_a=population_a,cell_type_b=population_b,
  cutoff_a=0.2,cutoff_b=0.2,pcg=pcg_mouse,lrp=lrp_mouse,ntrial=CONFIG$ntrial,
  cores=1L,echo=TRUE,dir_out=pair_dir)
 pvals<-as.data.table(result$pathways$df_pval)
 pcols<-intersect(names(pvals),c("pval_prize","pval_cost","pval_potential"))
 n_pathways<-length(result$pathways$raw)
 n_finite_p<-sum(is.finite(as.matrix(pvals[,.SD,.SDcols=pcols])))
 native_status<-fifelse(n_pathways==0L,"completed_native_network_no_pathways",
  fifelse(n_finite_p==3L*n_pathways,"completed_native","partial_native_pathway_significance"))
 result$implementation_repair <- "zero_degree_generalized_inverse;all_genes_retained;positive_degree_formula_unchanged"
 result$input_binding<-list(n_global_cells=ncol(lst$mat),n_pair_cells=sum(lst$cell_types %in% c(population_a,population_b)),
  n_input_genes=nrow(lst$mat),network_inputs_finite=TRUE,ntrial_requested=CONFIG$ntrial)
 result$native_phase_status<-native_status
 saveRDS(result,file.path(pair_dir,"cytotalk_result.rds"))
 fwrite(as.data.table(result$pcst$final_network),file.path(pair_dir,"FinalNetwork.txt"),sep="\t")
 status <- data.table(dataset=CONFIG$dataset,cell_type_a=population_a,cell_type_b=population_b,
  pair_key=pair_key,n_input_cells=ncol(lst$mat),n_pair_cells=sum(lst$cell_types %in% c(population_a,population_b)),n_native_edges=nrow(result$pcst$final_network),
  ntrial=CONFIG$ntrial,n_native_pathways=n_pathways,n_finite_pathway_pvalues=n_finite_p,
  seconds=proc.time()[["elapsed"]]-started,status=native_status)
 fwrite(status,progress,append=TRUE)
 failures[failures$pair_key==status$pair_key,"current_result_status":=native_status]
 fwrite(failures,file.path(work,"cytotalk_pair_failures.csv"))
 saveRDS(.Random.seed,file.path(work,"native_missing_pair_RNG.rds"))
 cat("PAIR_FINISHED",pair_key,status$n_native_edges,status$seconds,"\n");flush.console()
 list(status=as.character(status$status))
},by=.(pair_key,population_a,population_b,pair_dir)]
cat("COMPLETE_MISSING_PAIRS",nrow(result_summary),"\n")
