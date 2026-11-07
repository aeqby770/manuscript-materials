CONFIG <- list(root="D:/biotrade/LineaLinkR_0.6.5", threads=parallel::detectCores(), seed=12345, expression_pct=0.10, sender_pct=0.01, lfc_cutoff=0.25, covariance_cutoff=0.05, timetalk_root="D:/biotrade/0.64LLR/revision_20260919/results/by_algorithm/TimeTalk")
args <- commandArgs(TRUE)
CONFIG$dataset <- args[1L]
CONFIG$species <- args[2L]
CONFIG$input <- args[3L]
CONFIG$clock_mode <- args[4L]
CONFIG$clock_column <- args[5L]
CONFIG$method <- args[6L]
CONFIG$layout <- args[7L]
run_methods <- switch(CONFIG$method,both=c(TRUE,TRUE),NicheNet=c(TRUE,FALSE),ENTRAIN=c(FALSE,TRUE))
Sys.setenv(OMP_NUM_THREADS=CONFIG$threads,OPENBLAS_NUM_THREADS=CONFIG$threads,MKL_NUM_THREADS=CONFIG$threads)
.libPaths(c(file.path(CONFIG$root,"library"),.libPaths()))
library(data.table)
library(Matrix)
library(Seurat)
library(nichenetr)
library(dplyr)
library(tidyr)
setDTthreads(CONFIG$threads)
options(Seurat.object.assay.version="v3")
vendor <- file.path(CONFIG$root,"tools/competitors")
niche_out <- file.path(CONFIG$root,"results/by_algorithm/NicheNet",CONFIG$dataset)
entrain_out <- file.path(CONFIG$root,"results/by_algorithm/ENTRAIN",CONFIG$dataset)
dir.create(niche_out,recursive=TRUE,showWarnings=FALSE)
dir.create(entrain_out,recursive=TRUE,showWarnings=FALSE)
matrix_file <- switch(CONFIG$layout,prefixed=paste0(CONFIG$dataset,"_matrix.mtx"),plain="expr.mtx")
gene_file <- switch(CONFIG$layout,prefixed=paste0(CONFIG$dataset,"_genes.txt"),plain="genes.txt")
cell_file <- switch(CONFIG$layout,prefixed=paste0(CONFIG$dataset,"_cells.txt"),plain="cells.txt")
meta_file <- switch(CONFIG$layout,prefixed=paste0(CONFIG$dataset,"_meta.csv"),plain="meta.csv")
expr <- as(readMM(file.path(CONFIG$input,matrix_file)),"CsparseMatrix")
genes <- readLines(file.path(CONFIG$input,gene_file),encoding="UTF-8")
cells <- readLines(file.path(CONFIG$input,cell_file),encoding="UTF-8")
dimnames(expr) <- list(genes,cells)
meta <- fread(file.path(CONFIG$input,meta_file))
setnames(meta,1L,"barcode_id")
meta <- meta[match(cells,barcode_id)]
stopifnot(identical(meta$barcode_id,cells),!anyDuplicated(genes),!anyDuplicated(cells))
group_column <- switch(as.character(length(args)>=8L), "TRUE"=args[[8L]],
 "FALSE"=c(intersect(c("cell_type","clusters","celltype","cluster","paul15_clusters","ct","state_info"),names(meta)),".no_cell_type")[1L])
meta[, .no_cell_type:=NA_character_]
group <- as.character(meta[[group_column]])
input_is_counts <- all(expr@x>=0) && all(expr@x==floor(expr@x))
normalization <- switch(as.character(input_is_counts),"TRUE"="library_size_10000_log1p","FALSE"="provided_expression")
expr <- switch(as.character(input_is_counts),"TRUE"=NormalizeData(expr,normalization.method="LogNormalize",scale.factor=10000,verbose=FALSE),"FALSE"=expr)
print(data.table(dataset=CONFIG$dataset,n_cells=ncol(expr),n_genes=nrow(expr),normalization=normalization,group_column=group_column,n_groups=uniqueN(group)))
prior_filename <- switch(CONFIG$species,human="ligand_target_matrix_nsga2r_final.rds",mouse="ligand_target_matrix_nsga2r_final_mouse.rds")
ligand_target_matrix <- readRDS(file.path(vendor,"resources/NicheNet",prior_filename))
lr_network <- as.data.table(readRDS(file.path(vendor,"resources/NicheNet",paste0("lr_network_",CONFIG$species,"_21122021.rds"))))
switch(as.character(run_methods[1L]),"TRUE"={
optimized <- new.env(parent=asNamespace("nichenetr"))
invisible(lapply(list.files(file.path(vendor,"optimized/NicheNet/R"),pattern="\\.R$",full.names=TRUE),sys.source,envir=optimized))
meta_frame <- as.data.frame(meta)
rownames(meta_frame) <- cells
obj <- CreateSeuratObject(counts=expr,meta.data=meta_frame,min.cells=0,min.features=0)
seurat_expr <- expr
rownames(seurat_expr) <- rownames(obj)
obj[["RNA"]]@data <- seurat_expr
Idents(obj) <- factor(group)
stopifnot(nrow(obj)==length(genes),identical(colnames(obj),cells))
feature_names <- data.table(input_gene=genes,seurat_gene=rownames(obj))
fwrite(feature_names[input_gene!=seurat_gene],file.path(niche_out,"feature_api_names.csv"))
before <- proc.time()[["elapsed"]]
markers <- as.data.table(FindAllMarkers(obj,only.pos=FALSE,min.pct=CONFIG$expression_pct,logfc.threshold=CONFIG$lfc_cutoff,test.use="wilcox",max.cells.per.ident=Inf,random.seed=CONFIG$seed,return.thresh=1,verbose=TRUE))
markers[,gene:=genes[match(gene,rownames(obj))]]
marker_seconds <- proc.time()[["elapsed"]]-before
fwrite(markers,file.path(niche_out,"receiver_markers.csv.gz"))
membership <- sparseMatrix(i=seq_along(group),j=match(group,levels(Idents(obj))),x=1,dims=c(length(group),nlevels(Idents(obj))))
fraction <- sweep(as.matrix((expr>0)%*%membership),2,as.numeric(table(Idents(obj))),"/")
dimnames(fraction) <- list(genes,levels(Idents(obj)))
expressed_index <- which(fraction>=CONFIG$expression_pct,arr.ind=TRUE)
expressed <- data.table(gene=genes[expressed_index[,1]],receiver=colnames(fraction)[expressed_index[,2]])
background <- expressed[gene %in% rownames(ligand_target_matrix)]
genesets <- merge(markers[p_val_adj<=0.05 & abs(avg_log2FC)>=CONFIG$lfc_cutoff,.(gene,receiver=as.character(cluster))],background,by=c("gene","receiver"))
potential_network <- merge(lr_network[from %in% unique(expressed$gene) & from %in% colnames(ligand_target_matrix)],expressed,by.x="to",by.y="gene",allow.cartesian=TRUE)
potential <- unique(potential_network[,.(receiver,ligand=from)])
coverage <- data.table(receiver=levels(Idents(obj)),n_cells=as.numeric(table(Idents(obj))))
coverage <- merge(coverage,background[,.(n_background=.N),by=receiver],by="receiver",all.x=TRUE)
coverage <- merge(coverage,genesets[,.(n_response_genes=.N),by=receiver],by="receiver",all.x=TRUE)
coverage <- merge(coverage,potential[,.(n_ligands=.N),by=receiver],by="receiver",all.x=TRUE)
coverage[is.na(n_response_genes),n_response_genes:=0L]
coverage[is.na(n_ligands),n_ligands:=0L]
coverage[is.na(n_background),n_background:=0L]
coverage[,result_status:=fcase(n_cells<3,"not_applicable_too_few_receiver_cells",n_background==0,"not_applicable_no_background_prior_overlap",n_response_genes==0,"not_applicable_no_receiver_response_genes",n_response_genes>=n_background,"not_applicable_no_negative_background",n_ligands==0,"not_applicable_no_candidate_ligands",default="ready")]
valid <- coverage[result_status=="ready",receiver]
gene_lists <- split(genesets$gene,genesets$receiver)[valid]
background_lists <- split(background$gene,background$receiver)[valid]
ligand_lists <- split(potential$ligand,potential$receiver)[valid]
before <- proc.time()[["elapsed"]]
activities <- rbindlist(mapply(optimized$predict_ligand_activities,geneset=gene_lists,background_expressed_genes=background_lists,potential_ligands=ligand_lists,MoreArgs=list(ligand_target_matrix=ligand_target_matrix),SIMPLIFY=FALSE),idcol="receiver")
activity_seconds <- proc.time()[["elapsed"]]-before
activities[,dataset:=CONFIG$dataset]
fwrite(activities,file.path(niche_out,"ligand_activities.csv"))
coverage[result_status=="ready",result_status:="scored"]
coverage[,`:=`(dataset=CONFIG$dataset,method="NicheNet",n_input_cells=ncol(expr),n_input_genes=nrow(expr),normalization=normalization,marker_seconds=marker_seconds,activity_seconds=activity_seconds,contrast="receiver_vs_all_other_cell_types",ground_truth_accuracy_available=FALSE)]
fwrite(coverage,file.path(niche_out,"coverage.csv"))
saveRDS(list(activities=activities,genesets=gene_lists,background=background_lists,potential_ligands=ligand_lists,config=CONFIG),file.path(niche_out,"nichenet_native.rds"))
print(coverage)
})
switch(as.character(run_methods[2L]),"TRUE"={
clock <- switch(CONFIG$clock_mode,
    feasible_clock={clock_object<-readRDS(file.path(CONFIG$root,"cache/feasible_clock",CONFIG$dataset,"linealinkr_complete.rds"));clock_object$pseudotime[match(cells,clock_object$cells)]},
    observed=as.numeric(meta[[CONFIG$clock_column]]),
    TimeTalk={clock_table<-fread(file.path(CONFIG$timetalk_root,CONFIG$dataset,"timetalk_cell_inputs.csv"));stopifnot(setequal(cells,clock_table$cell_id));clock_table$pseudotime[match(cells,clock_table$cell_id)]},
    TICCI={clock_table<-fread(file.path(CONFIG$root,"inputs/comparison/results/by_algorithm/TICCI",CONFIG$dataset,"ticci_pseudotime.csv"));clock_table$pseudotime[match(cells,clock_table$cell_id)]})
fwrite(data.table(cell_id=cells,pseudotime=clock,cell_type=group,clock_source=CONFIG$clock_mode),file.path(entrain_out,"input_clock.csv.gz"))
clock_valid <- is.finite(clock)
entrain_expr <- expr[,clock_valid,drop=FALSE]
entrain_clock <- clock[clock_valid]
entrain <- new.env(parent=globalenv())
sys.source(file.path(vendor,"optimized/ENTRAIN/R/entrain_fitting.R"),envir=entrain)
expressed_ligands <- intersect(genes[Matrix::rowMeans(expr>0)>=CONFIG$sender_pct],colnames(ligand_target_matrix))
receiver_genes <- genes[Matrix::rowMeans(entrain_expr>0)>=CONFIG$expression_pct]
active_network <- lr_network[from %in% expressed_ligands & to %in% receiver_genes]
active_ligands <- unique(active_network$from)
background_genes <- entrain$get_background_genes_path(entrain_expr,pct=CONFIG$expression_pct)
before <- proc.time()[["elapsed"]]
training <- entrain$get_branch_genes_covar(as.matrix(entrain_expr[background_genes,,drop=FALSE]),entrain_clock)
cutoff <- stats::quantile(training$Covariances,CONFIG$covariance_cutoff)
training_top <- training[training$Covariances>cutoff,,drop=FALSE]
set.seed(CONFIG$seed)
fit <- entrain$get_ligand_trajectory_scores_regression(training_top,ligand_target_matrix[,active_ligands,drop=FALSE])
entrain_seconds <- proc.time()[["elapsed"]]-before
importance <- sort(fit$model$importance[,"IncNodePurity"],decreasing=TRUE)
fwrite(data.table(dataset=CONFIG$dataset,ligand=names(importance),importance=importance),file.path(entrain_out,"ligand_importance.csv"))
fwrite(as.data.table(training,keep.rownames="gene"),file.path(entrain_out,"pseudotime_gene_covariance.csv.gz"))
saveRDS(list(model=fit$model,importance=importance,active_network=active_network,training_genes=rownames(training_top),clock_source=CONFIG$clock_mode,config=CONFIG),file.path(entrain_out,"entrain_native.rds"))
entrain_coverage <- data.table(dataset=CONFIG$dataset,method="ENTRAIN",n_input_cells=ncol(expr),n_used_cells=ncol(entrain_expr),n_nonfinite_clock=sum(!clock_valid),n_input_genes=nrow(expr),n_ligands=length(importance),n_trees=fit$model$ntree,normalization=normalization,clock_source=CONFIG$clock_mode,result_status=fifelse(all(clock_valid),"scored","scored_native_nonfinite_clock_exclusion"),core_seconds=entrain_seconds,ground_truth_accuracy_available=FALSE)
fwrite(entrain_coverage,file.path(entrain_out,"coverage.csv"))
print(entrain_coverage)
})
