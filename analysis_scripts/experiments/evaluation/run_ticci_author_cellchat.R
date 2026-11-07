CONFIG <- list(root=Sys.getenv("LLR_WORK_ROOT"),dataset=commandArgs(TRUE)[2L],
 threads=as.integer(Sys.getenv("LLR_NATIVE_THREADS",parallel::detectCores())),
 seed=12345L,K_max=10L,gap_B=50L,workers=parallel::detectCores(logical=FALSE),
 kernel=file.path(Sys.getenv("LLR_WORK_ROOT"),"tools/competitors/optimized/TICCI/ticci_kmeans.R"))
Sys.setenv(OMP_NUM_THREADS=CONFIG$threads,OPENBLAS_NUM_THREADS=CONFIG$threads,MKL_NUM_THREADS=CONFIG$threads)
.libPaths(c(file.path(CONFIG$root,"library"),.libPaths()))
library(data.table)
library(Matrix)
library(rhdf5)
library(cluster)
library(CellChat)
setDTthreads(CONFIG$threads)
out <- file.path(CONFIG$root,"results/by_algorithm/TICCI_author_functional_entropy",CONFIG$dataset)
h5file <- file.path(out,"entropy_input.h5")
shape <- h5read(h5file,"shape")
expr <- new("dgCMatrix",i=as.integer(h5read(h5file,"indices")),
 p=as.integer(h5read(h5file,"indptr")),x=as.numeric(h5read(h5file,"data")),
 Dim=as.integer(rev(shape)))
genes <- readLines(file.path(out,"entropy_genes.txt"),encoding="UTF-8")
cells <- readLines(file.path(out,"entropy_cells.txt"),encoding="UTF-8")
dimnames(expr) <- list(genes,cells)
expr <- expr%*%Diagonal(x=10000/colSums(expr))
expr@x <- log1p(expr@x)
groups <- fread(file.path(out,"ticci_cellchat_groups.csv"),colClasses="character")
group_complete <- nrow(groups)==ncol(expr)&&setequal(groups$cell_id,cells)
switch(as.character(group_complete),"FALSE"={
 X <- t(as.matrix(expr))
 native_status <- fread(file.path(out,"run_status.csv"),colClasses="character")
 native_status[,":="(status="author_Gap_K10_B50_parallel_Hartigan_exact",
  gap_reference_completed=0L,gap_reference_total=CONFIG$gap_B)]
 fwrite(native_status,file.path(out,"run_status.csv"))
 source(CONFIG$kernel,encoding="UTF-8")
 gap_fast <- cluster::clusGap
 environment(gap_fast) <- new.env(parent=environment(gap_fast))
 environment(gap_fast)$gap_progress_file <- file.path(out,"ticci_author_gap_partial.rds")
 environment(gap_fast)$gap_status_file <- file.path(out,"run_status.csv")
 environment(gap_fast)$gap_status <- native_status
 environment(gap_fast)$gap_workers <- CONFIG$workers
 gap_body <- as.list(body(gap_fast))
 gap_declarations <- vapply(gap_body,deparse1,character(1L))
 wk <- which(startsWith(gap_declarations,"W.k <- function"))
 stopifnot(length(wk)==1L)
 gap_body[[wk]] <- quote(W.k <- function(X,kk) {
  clus <- if(kk>1) FUNcluster(X,kk,...)$cluster else rep.int(1L,nrow(X))
  switch(as.character(kk==1L),"TRUE"={
   distance_cache <<- NULL
   environment(FUNcluster)$.rowmajor <- NULL
   environment(FUNcluster)$.totss <- NA_real_
   gc()
   data.table::set(gap_status,j="gap_active_reference",
    value=get0("b",envir=parent.frame(),inherits=FALSE,ifnotfound=0L))
   data.table::set(gap_status,j="gap_phase",value="full_gene_pairwise_distances")
   data.table::fwrite(gap_status,gap_status_file)
   norms <- rowSums(X*X)
   pyX <- reticulate::r_to_py(X)
   gram <- reticulate::py_to_r(blas_numpy$matmul(pyX,blas_numpy$transpose(pyX)))
   distance_cache <<- sqrt(.Call(environment(FUNcluster)$.geometry,gram,norms))^d.power
   diag(distance_cache) <<- 0
   rm(pyX,gram,norms)
   gc()
   environment(FUNcluster)$.totss <- sum(scale(X,scale=FALSE)^2)
   data.table::set(gap_status,j="gap_phase",value="parallel_population_clustering")
   data.table::fwrite(gap_status,gap_status_file)
   environment(FUNcluster)$.prepare(X,K.max,gap_workers)
  })
  0.25*sum(vapply(split(ii,clus),function(I) sum(distance_cache[I,I,drop=FALSE])/length(I),0))
 })
 gap_body <- append(gap_body,list(quote(distance_cache<-NULL),
 quote(blas_numpy<-reticulate::import("numpy",convert=FALSE))),after=wk-1L)
 body(gap_fast) <- as.call(gap_body)
 reticulate::use_python("E:/Python/biotrade/Scripts/python.exe",required=TRUE)
 gap_text <- paste(deparse(body(gap_fast),width.cutoff=500L),collapse="\n")
 gap_text <- sub("xs <- xs %*% V.sx",
  "xs <- reticulate::py_to_r(blas_numpy$matmul(reticulate::r_to_py(xs),reticulate::r_to_py(V.sx)))",
  gap_text,fixed=TRUE)
 gap_text <- sub("tcrossprod(z1, V.sx)",
  "reticulate::py_to_r(blas_numpy$matmul(reticulate::r_to_py(z1),blas_numpy$transpose(reticulate::r_to_py(V.sx))))",
  gap_text,fixed=TRUE)
 gap_text <- sub("V.sx <- svd(xs, nu = 0)$v",
  paste0("V.sx <- switch(as.character(file.exists('",gsub("\\\\","/",file.path(out,"ticci_gap_native_svd.rds")),"')), ",
   "'TRUE'=readRDS('",gsub("\\\\","/",file.path(out,"ticci_gap_native_svd.rds")),"'), ",
   "'FALSE'={v <- reticulate::py_to_r(blas_numpy$transpose(reticulate::py_get_item(blas_numpy$linalg$svd(reticulate::r_to_py(xs),full_matrices=FALSE),2L)));saveRDS(v,'",gsub("\\\\","/",file.path(out,"ticci_gap_native_svd.rds")),
   "',compress=FALSE);v})"),
  gap_text,fixed=TRUE)
 gap_text <- sub("xs <- scale(x, center = TRUE, scale = FALSE)",
  "distance_cache <- NULL;environment(FUNcluster)$.rowmajor <- NULL;gc();xs <- scale(x, center=TRUE,scale=FALSE)",
  gap_text,fixed=TRUE)
 body(gap_fast) <- parse(text=gap_text)[[1L]]
 gap_parts <- as.list(body(gap_fast))
 gap_declarations <- vapply(gap_parts,deparse1,character(1L))
 bootstrap_index <- which(startsWith(gap_declarations,"for (b in 1:B)"))
 stopifnot(length(bootstrap_index)==1L)
 bootstrap_parts <- as.list(gap_parts[[bootstrap_index]][[4L]])
 bootstrap_parts <- append(bootstrap_parts,list(quote({
  gap_reference_started <- proc.time()[["elapsed"]]
  data.table::set(gap_status,j="gap_active_reference",value=b)
  data.table::set(gap_status,j="gap_phase",value="reference_generation")
  data.table::fwrite(gap_status,gap_status_file)
 })),after=1L)
 bootstrap_parts <- append(bootstrap_parts,list(quote({
  gap_reference_seconds <- proc.time()[["elapsed"]]-gap_reference_started
  saveRDS(list(completed_reference=b,logWks=logWks,
   rng_state=get(".Random.seed",envir=.GlobalEnv,inherits=FALSE),
   K_max=K.max,B=B,n_cells=nrow(x),n_genes=ncol(x)),
   gap_progress_file)
  data.table::set(gap_status,j="status",value="author_Gap_K10_B50_parallel_Hartigan_exact")
  data.table::set(gap_status,j="gap_phase",value="reference_complete")
  data.table::set(gap_status,j="gap_reference_completed",value=b)
  data.table::set(gap_status,j="gap_reference_total",value=B)
  data.table::set(gap_status,j="gap_last_reference_seconds",value=gap_reference_seconds)
  data.table::fwrite(gap_status,gap_status_file)
  cat("TICCI_GAP_REFERENCE_COMPLETE",b,B,"seconds",gap_reference_seconds,"\n")
 })))
 gap_parts[[bootstrap_index]][[4L]] <- as.call(bootstrap_parts)
 gap_parts[[bootstrap_index]][[3L]] <- quote(gap_b_sequence)
 gap_parts <- append(gap_parts,list(quote({
  gap_partial <- switch(as.character(file.exists(gap_progress_file)),
   "TRUE"=readRDS(gap_progress_file),"FALSE"=NULL)
  gap_b_sequence <- seq_len(B)
  switch(as.character(is.null(gap_partial)),"FALSE"={
   stopifnot(gap_partial$K_max==K.max,gap_partial$B==B,
    gap_partial$n_cells==nrow(x),gap_partial$n_genes==ncol(x))
   logWks <- gap_partial$logWks
   assign(".Random.seed",gap_partial$rng_state,envir=.GlobalEnv)
   gap_b_sequence <- gap_b_sequence[gap_b_sequence>gap_partial$completed_reference]
  })
 })),after=bootstrap_index-1L)
 body(gap_fast) <- as.call(gap_parts)
 set.seed(CONFIG$seed)
 gap_state <- switch(as.character(file.exists(file.path(out,"ticci_author_gap.rds"))),
  "TRUE"=readRDS(file.path(out,"ticci_author_gap.rds")),
  "FALSE"={
   started <- proc.time()[["elapsed"]]
   gaps <- gap_fast(X,FUNcluster=kmeans_cached,K.max=CONFIG$K_max,B=CONFIG$gap_B,verbose=TRUE)
   result <- list(gaps=gaps,rng_state=.Random.seed,seconds=proc.time()[["elapsed"]]-started)
   saveRDS(result,file.path(out,"ticci_author_gap.rds"))
   result
  })
 ticci_finish_parallel()
 .Random.seed <- gap_state$rng_state
 best_k <- which.max(gap_state$gaps$Tab[,"gap"])
 ticci_kmeans_state$.totss <- sum(scale(X,scale=FALSE)^2)
 ticci_kmeans_state$.rowmajor <- t(X)
 grouped <- kmeans_fast(X,centers=best_k)
 groups <- data.table(cell_id=cells,Classification=as.character(grouped$cluster))
 fwrite(groups,file.path(out,"ticci_cellchat_groups.csv"))
 fwrite(data.table(k=seq_len(nrow(gap_state$gaps$Tab)),as.data.table(gap_state$gaps$Tab)),
  file.path(out,"ticci_author_gap_statistics.csv"))
 fwrite(data.table(n_cells=nrow(X),n_input_genes=ncol(X),K_max=CONFIG$K_max,
  B=CONFIG$gap_B,best_k=best_k,final_kmeans_ifault=grouped$ifault,
  native_RNG_sequence_preserved=TRUE,dispersion_engine="cached_exact_full_gene_Euclidean_distances",
  svd_engine="DGESDD_NumPy_OpenBLAS_ILP64",n_SVD_dimensions=min(nrow(X),ncol(X)),
  clustering_engine="R453_AS136_contiguous_shared_matrix_exact_order",
  clustering_workers=CONFIG$workers,random_center_order="unchanged_native_R_sequential",
  cell_or_gene_subsetting=FALSE),file.path(out,"ticci_author_grouping_summary.csv"))
 ticci_kmeans_state$.rowmajor <- NULL
 rm(X,gap_state,grouped)
 gc()
})
groups <- groups[match(cells,cell_id)]
stopifnot(identical(groups$cell_id,cells),nrow(groups)==ncol(expr),!anyNA(groups$Classification))
gene_file <- file.path(out,"ticci_gene_id_binding.csv")
switch(as.character(file.exists(gene_file)),"TRUE"={
 gene_binding <- fread(gene_file)
 stopifnot(identical(genes,gene_binding$original_gene_id))
 rownames(expr) <- gene_binding$resource_symbol
})
native_status <- fread(file.path(out,"run_status.csv"),colClasses="character")
native_status[,status:="author_CellChat"]
fwrite(native_status,file.path(out,"run_status.csv"))
meta <- data.frame(Classification=factor(groups$Classification),row.names=cells)
set.seed(CONFIG$seed)
cellchat <- createCellChat(object=expr,meta=meta,group.by="Classification")
species <- fread(file.path(out,"input_binding.csv"))$species
cellchat@DB <- list(mouse=CellChatDB.mouse,human=CellChatDB.human,zebrafish=CellChatDB.zebrafish)[[species]]
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat,raw.use=TRUE)
cellchat <- filterCommunication(cellchat,min.cells=10)
stopifnot(ncol(cellchat@data)==length(cells))
interactions <- subsetCommunication(cellchat)
fwrite(interactions,file.path(out,"ticci_cellchat_interactions.csv"))
saveRDS(cellchat,file.path(out,"ticci_cellchat.rds"))
fwrite(data.table(dataset=CONFIG$dataset,n_cells=ncol(cellchat@data),
 n_groups=nlevels(cellchat@idents),n_interactions=nrow(interactions),
 group_mode="author_kmeans_and_Gap",group_Kmax=CONFIG$K_max,group_bootstrap_B=CONFIG$gap_B,
 input_genes_retained=length(genes),cellchat_version=as.character(packageVersion("CellChat"))),
 file.path(out,"ticci_cellchat_summary.csv"))
cat("TICCI_CELLCHAT_COMPLETED",CONFIG$dataset,nrow(interactions),"interactions\n")
