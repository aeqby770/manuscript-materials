CONFIG<-list(root=normalizePath(Sys.getenv("LLR_WORK_ROOT",getwd()),winslash="/",mustWork=TRUE),
 input="results/by_algorithm/GSE48968_BMDC_WT/log_tpm_input",
 output="results/by_algorithm/TimeTalk/GSE48968_BMDC_WT",
 pairs="results/by_algorithm/GSE48968_BMDC_WT/all_lr_pairs.csv",
 common_pairs="results/by_algorithm/GSE48968_BMDC_WT/same_resource_pair_rankings.csv",
 common_ligands="results/validation/gse48968/all_method_ligand_rankings.csv",
 threads=parallel::detectCores(),workers=parallel::detectCores(),seed=12345L)
Sys.setenv(OMP_NUM_THREADS=CONFIG$threads,OPENBLAS_NUM_THREADS=CONFIG$threads,MKL_NUM_THREADS=CONFIG$threads)
.libPaths(c(file.path(CONFIG$root,"library"),.libPaths()))
library(data.table)
library(Matrix)
library(Seurat)
library(monocle3)
library(TimeTalk)
library(RTN)
library(SeuratWrappers)
library(dorothea)
library(dplyr)
library(tidyr)
library(tibble)
library(future)
library(future.apply)
library(snow)
library(lmtest)
library(limma)
setDTthreads(CONFIG$threads)
options(Seurat.object.assay.version="v3",future.globals.maxSize=16*1024^3)
set.seed(CONFIG$seed)
IN<-file.path(CONFIG$root,CONFIG$input)
OUT<-file.path(CONFIG$root,CONFIG$output)
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
pairs<-fread(file.path(CONFIG$root,CONFIG$pairs))
pair_common<-unique(fread(file.path(CONFIG$root,CONFIG$common_pairs))[
 in_common_evaluable_panel%in%TRUE,.(pair_idx,ligand_genes,receptor_genes)])
ligand_common<-unique(fread(file.path(CONFIG$root,CONFIG$common_ligands))[
 is.finite(rank_common),.(ligand)])
stopifnot(nrow(pair_common)==420L,nrow(ligand_common)==138L)
COMPLEX_PARTS<-unique(rbindlist(list(
 pairs[,.(partner=ligand_genes)],pairs[,.(partner=receptor_genes)])))
COMPLEX_PARTS<-COMPLEX_PARTS[,.(gene=unlist(strsplit(partner,"_",fixed=TRUE))),by=partner]
engine<-new.env(parent=asNamespace("TimeTalk"))
engine$CONFIG<-CONFIG
engine$COMPLEX_PARTS<-COMPLEX_PARTS
engine$data.table<-data.table::data.table
engine$setorder<-data.table::setorder
engine$dcast<-data.table::dcast
engine$pull<-dplyr::pull
engine$rownames_to_column<-tibble::rownames_to_column
engine$plan<-future::plan
engine$future_lapply<-future.apply::future_lapply
engine$grangertest<-lmtest::grangertest
engine$RunALRA<-SeuratWrappers::RunALRA
engine$stopCluster<-snow::stopCluster
variable_run<-TimeTalk::GetVariableGeneInTrajectory
environment(variable_run)<-engine
tt_text<-paste(deparse(TimeTalk::RunTimeTalk,width.cutoff=500L),collapse="\n")
# Enable data.table semantics inside the cloned native function.
tt_text<-sub("{","{\n.datatable.aware <- TRUE",tt_text,fixed=TRUE)
stopifnot(grepl("LRpairs <- LRpairs[-index.remove]",tt_text,fixed=TRUE))
# Preserve eligible genes when the removal set is empty.
tt_text<-sub("LRpairs <- LRpairs[-index.remove]","LRpairs <- LRpairs[setdiff(seq_along(LRpairs),index.remove)]",tt_text,fixed=TRUE)
tt_text<-sub("Lgenelist <- Lgenelist[-index.remove]","Lgenelist <- Lgenelist[setdiff(seq_along(Lgenelist),index.remove)]",tt_text,fixed=TRUE)
tt_text<-sub("Rgenelist <- Rgenelist[-index.remove]","Rgenelist <- Rgenelist[setdiff(seq_along(Rgenelist),index.remove)]",tt_text,fixed=TRUE)
tt_text<-sub("gene_symbols <- rownames(tmp.data)","gene_symbols <- unique(c(rownames(tmp.data),COMPLEX_PARTS$partner))",tt_text,fixed=TRUE)
tt_text<-sub("PCC = unlist(tmp.res.list)","PCC = unlist(tmp.res.list.PCC)",tt_text,fixed=TRUE)
tt_text<-sub("SCC = unlist(tmp.res.list)","SCC = unlist(tmp.res.list.SCC)",tt_text,fixed=TRUE)
# Average measured subunits first, then take the limiting subunit before native scaling.
complex_code<-"\nprofile_parts_1<-data.table(\n partner=rep(COMPLEX_PARTS$partner,each=ncol(inter.tmp.mat.1$interpolatedVals)),\n grid=rep(seq_len(ncol(inter.tmp.mat.1$interpolatedVals)),times=nrow(COMPLEX_PARTS)),\n value=as.vector(t(inter.tmp.mat.1$interpolatedVals[COMPLEX_PARTS$gene,,drop=FALSE])),\n error=as.vector(t(inter.tmp.mat.1$error[COMPLEX_PARTS$gene,,drop=FALSE])))\nsetorder(profile_parts_1,partner,grid,value)\nprofile_parts_1<-unique(profile_parts_1,by=c(\"partner\",\"grid\"))\ncompound_values_1<-dcast(profile_parts_1,partner~grid,value.var=\"value\")\ncompound_errors_1<-dcast(profile_parts_1,partner~grid,value.var=\"error\")\ncompound_names_1<-setdiff(compound_values_1$partner,rownames(inter.tmp.mat.1$interpolatedVals))\ncompound_matrix_1<-as.matrix(compound_values_1[match(compound_names_1,partner),-\"partner\"])\ncompound_error_matrix_1<-as.matrix(compound_errors_1[match(compound_names_1,partner),-\"partner\"])\ndimnames(compound_matrix_1)<-list(compound_names_1,colnames(inter.tmp.mat.1$interpolatedVals))\ndimnames(compound_error_matrix_1)<-list(compound_names_1,colnames(inter.tmp.mat.1$error))\ninter.tmp.mat.1$interpolatedVals<-rbind(inter.tmp.mat.1$interpolatedVals,compound_matrix_1)\ninter.tmp.mat.1$error<-rbind(inter.tmp.mat.1$error,compound_error_matrix_1)\n\n\nprofile_parts_2<-data.table(\n partner=rep(COMPLEX_PARTS$partner,each=ncol(inter.tmp.mat.2$interpolatedVals)),\n grid=rep(seq_len(ncol(inter.tmp.mat.2$interpolatedVals)),times=nrow(COMPLEX_PARTS)),\n value=as.vector(t(inter.tmp.mat.2$interpolatedVals[COMPLEX_PARTS$gene,,drop=FALSE])),\n error=as.vector(t(inter.tmp.mat.2$error[COMPLEX_PARTS$gene,,drop=FALSE])))\nsetorder(profile_parts_2,partner,grid,value)\nprofile_parts_2<-unique(profile_parts_2,by=c(\"partner\",\"grid\"))\ncompound_values_2<-dcast(profile_parts_2,partner~grid,value.var=\"value\")\ncompound_errors_2<-dcast(profile_parts_2,partner~grid,value.var=\"error\")\ncompound_names_2<-setdiff(compound_values_2$partner,rownames(inter.tmp.mat.2$interpolatedVals))\ncompound_matrix_2<-as.matrix(compound_values_2[match(compound_names_2,partner),-\"partner\"])\ncompound_error_matrix_2<-as.matrix(compound_errors_2[match(compound_names_2,partner),-\"partner\"])\ndimnames(compound_matrix_2)<-list(compound_names_2,colnames(inter.tmp.mat.2$interpolatedVals))\ndimnames(compound_error_matrix_2)<-list(compound_names_2,colnames(inter.tmp.mat.2$error))\ninter.tmp.mat.2$interpolatedVals<-rbind(inter.tmp.mat.2$interpolatedVals,compound_matrix_2)\ninter.tmp.mat.2$error<-rbind(inter.tmp.mat.2$error,compound_error_matrix_2)\n"
needle<-"inter.tmp.mat.1 <- cellAlign::scaleInterpolate(inter.tmp.mat.1)"
stopifnot(grepl(needle,tt_text,fixed=TRUE))
tt_text<-sub(needle,paste0(complex_code,"\n",needle),tt_text,fixed=TRUE)
needle<-"if (length(tmp.TF.gene) == 0L)"
tt_text<-sub(needle,paste0("tmp.TF.gene <- intersect(tmp.TF.gene,rownames(inter.tmp.mat.2))\n",needle),tt_text,fixed=TRUE)
# The autocrine fields contain the same BMDC cells and observed-time coordinate.
interp1<-'inter.tmp.mat.1 <- cellAlign::interWeights(expDataBatch = tmp.mat.1, trajCond = tmp.mat.pseudotime.1, winSz = tmp.winsz, numPts = numPts)'
interp2<-'inter.tmp.mat.2 <- cellAlign::interWeights(expDataBatch = tmp.mat.2, trajCond = tmp.mat.pseudotime.2, winSz = tmp.winsz, numPts = numPts)'
stopifnot(grepl(interp1,tt_text,fixed=TRUE),grepl(interp2,tt_text,fixed=TRUE))
profile_text<-paste0(
 'stopifnot(identical(tmp.cell.meta.1$cell_id,tmp.cell.meta.2$cell_id),identical(tmp.mat.pseudotime.1,tmp.mat.pseudotime.2))\n',
 'profile_file<-file.path(CONFIG$root,CONFIG$output,paste0("TimeTalk_",tmp.orig.ident,"_interpolation.rds"))\n',
 'inter.tmp.mat.1<-switch(as.character("--reuse-profiles"%in%commandArgs(TRUE)&&file.exists(profile_file)),"TRUE"=readRDS(profile_file),"FALSE"=cellAlign::interWeights(expDataBatch=tmp.mat.1,trajCond=tmp.mat.pseudotime.1,winSz=tmp.winsz,numPts=numPts))\n',
 'saveRDS(inter.tmp.mat.1,profile_file)')
tt_text<-sub(interp1,profile_text,tt_text,fixed=TRUE)
tt_text<-sub(interp2,"inter.tmp.mat.2 <- inter.tmp.mat.1",tt_text,fixed=TRUE)
tt_run<-eval(parse(text=tt_text),envir=engine)
writeLines(tt_text,file.path(OUT,"TimeTalk_bound_core.R"),useBytes=TRUE)
switch(as.character(!any(c("--score-only","--mra-only","--aggregate-only")%in%commandArgs(TRUE))),"TRUE"={
 expr<-as(readMM(file.path(IN,"expr.mtx")),"CsparseMatrix")
 genes<-readLines(file.path(IN,"genes.txt"),encoding="UTF-8")
 cells<-readLines(file.path(IN,"cells.txt"),encoding="UTF-8")
 stopifnot(!anyDuplicated(toupper(genes)))
 dimnames(expr)<-list(toupper(genes),cells)
 meta<-as.data.frame(fread(file.path(IN,"meta.csv")))
 rownames(meta)<-meta$barcode_id
 meta<-meta[cells,,drop=FALSE]
 stopifnot(identical(rownames(meta),colnames(expr)),all(COMPLEX_PARTS$gene%in%rownames(expr)))
 linear_tpm<-expr
 linear_tpm@x<-expm1(linear_tpm@x)*100
 obj<-CreateSeuratObject(counts=linear_tpm,meta.data=meta,min.cells=0,min.features=0)
 gene_binding<-data.table(input_gene=rownames(expr),RNA_feature=rownames(obj))
 fwrite(gene_binding,file.path(OUT,"TimeTalk_feature_binding.csv"))
 rownames(expr)<-rownames(obj)
 obj[["RNA"]]@data<-expr
 obj$CellType<-"BMDC"
 Idents(obj)<-"CellType"
 obj<-FindVariableFeatures(obj,nfeatures=2000,verbose=FALSE)
 obj<-ScaleData(obj,features=VariableFeatures(obj),verbose=FALSE)
 obj<-RunPCA(obj,npcs=50,seed.use=CONFIG$seed,verbose=FALSE)
 obj<-RunUMAP(obj,dims=1:30,umap.method="uwot",seed.use=CONFIG$seed,verbose=FALSE)
 count_data<-Seurat::GetAssayData(obj,assay="RNA",layer="counts")
 gene_data<-data.frame(gene_short_name=rownames(count_data),row.names=rownames(count_data))
 cds<-monocle3::new_cell_data_set(count_data,cell_metadata=meta,gene_metadata=gene_data)
 SingleCellExperiment::reducedDim(cds,"UMAP")<-Seurat::Embeddings(obj,"umap")
 cds@principal_graph_aux[["UMAP"]]$pseudotime<-setNames(as.numeric(meta$time_hours),cells)
 saveRDS(list(expr=expr,meta=meta,obj=obj,cds=cds,gene_binding=gene_binding),
  file.path(OUT,"TimeTalk_RNA_input.rds"))
})
switch(as.character(!any(c("--score-only","--aggregate-only")%in%commandArgs(TRUE))),"TRUE"={
 list2env(readRDS(file.path(OUT,"TimeTalk_RNA_input.rds")),envir=environment())
 condition_group<-factor(paste(meta$stimulation,meta$time_hours,sep="_"))
 design<-model.matrix(~condition_group)
 fit<-limma::eBayes(limma::lmFit(as.matrix(expr),design))
 association<-as.data.table(limma::topTable(fit,coef=seq.int(2L,ncol(design)),
  number=Inf,sort.by="none"),keep.rownames="gene")
 fwrite(association,file.path(OUT,"TimeTalk_time_response_phenotype.csv"))
 hits<-association[adj.P.Val<.05&is.finite(F)]
 phenotype<-setNames(hits$F,hits$gene)
 utils::data("dorothea_mm",package="dorothea")
 gene_tf<-intersect(toupper(unique(dorothea_mm$tf)),rownames(expr))
 stopifnot(length(gene_tf)>0L,length(phenotype)>0L)
 cat("BMDC_TRN_INPUT",ncol(expr),"cells",nrow(expr),"genes",length(gene_tf),"TFs",length(phenotype),"phenotype genes\n")
 trn_text<-paste(deparse(TimeTalk::GetCellTypeSpecificTRN,width.cutoff=500L),collapse="\n")
 needle<-'options(cluster = snow::makeCluster(spec = 10, "SOCK"))'
 stopifnot(grepl(needle,trn_text,fixed=TRUE))
 trn_text<-sub(needle,paste0(needle,"\nparallel::clusterSetRNGStream(getOption('cluster'),CONFIG$seed)"),trn_text,fixed=TRUE)
 trn_run<-eval(parse(text=trn_text),envir=engine)
 rtna<-trn_run(obj,tmp.ident="BMDC",gene.TF=gene_tf,tmp.phenotype=phenotype)
 saveRDS(rtna,file.path(OUT,"TimeTalk_RTN_MRA.rds"))
 mra<-as.data.table(RTN::tna.get(rtna,what="mra",ntop=-1))
 fwrite(mra,file.path(OUT,"TimeTalk_all_MRA_TFs.csv"))
 mra_sig<-mra[Pvalue<.05&Regulon%chin%rownames(expr),.(Regulon,Pvalue,group="BMDC")]
 saveRDS(list(obj=obj,cds=cds,meta=meta,gene_binding=gene_binding,mra_sig=mra_sig,n_mra_tested=nrow(mra),
  n_genes=nrow(expr),n_cells=ncol(expr)),file.path(OUT,"TimeTalk_score_input.rds"))
 cat("BMDC_MRA_COMPLETE",nrow(mra_sig),"positive TFs\n")
})
list2env(readRDS(file.path(OUT,"TimeTalk_score_input.rds")),envir=environment())
COMPLEX_PARTS[,gene:=gene_binding$RNA_feature[match(gene,gene_binding$input_gene)]]
stopifnot(!anyNA(COMPLEX_PARTS$gene))
engine$COMPLEX_PARTS<-COMPLEX_PARTS
LRpairs.df<-as.data.frame(unique(pairs[,.(lr_pair=paste(ligand_genes,receptor_genes,sep="-"),
 ligand_gene_symbol=ligand_genes,receptor_gene_symbol=receptor_genes)]))
switch(as.character(!"--aggregate-only"%in%commandArgs(TRUE)),"TRUE"={
conditions<-data.table(condition=sort(unique(as.character(meta$stimulation[meta$time_hours>0]))))
result<-conditions[,{
 obj$orig.ident<-as.character(meta$stimulation)
 obj$orig.ident[meta$time_hours==0]<-condition
 native<-as.data.table(tt_run(cds,obj,tmp.ident.1="BMDC",tmp.ident.2="BMDC",
  LRpairs.df=LRpairs.df,tmp.mra.res=as.data.frame(mra_sig),tmp.orig.ident=condition,
  tmp.winsz=.1,tmp.lags=1L,numPts=200L,tmp.cores=CONFIG$workers,
  tmp.SCC.cutoff=.2,tmp.granger.cutoff=.01))
 fwrite(native,file.path(OUT,paste0("TimeTalk_",condition,"_TF_results.csv.gz")))
 native
},by=condition]
fwrite(result,file.path(OUT,"timetalk_result.csv.gz"))
})
result<-fread(file.path(OUT,"timetalk_result.csv.gz"))
condition_strength<-result[,.(n_TF_rows=.N,n_pass=sum(category=="PASS"),
 SCC_ens=max(fifelse(category=="PASS",SCC,-Inf))),by=.(condition,L,R)]
condition_strength[,score:=fifelse(n_pass>0&abs(SCC_ens)>.8,abs(SCC_ens),0)]
fwrite(condition_strength,file.path(OUT,"TimeTalk_condition_pair_strength.csv"))
strength<-condition_strength[,.(n_TF_rows=sum(n_TF_rows),n_pass=sum(n_pass),score=max(score)),by=.(L,R)]
all_scores<-merge(pairs,strength,by.x=c("ligand_genes","receptor_genes"),by.y=c("L","R"),all.x=TRUE)
all_scores[,c("score","n_TF_rows","n_pass"):=.(fcoalesce(score,0),fcoalesce(n_TF_rows,0L),fcoalesce(n_pass,0L))]
common<-merge(pair_common,all_scores[,.(pair_idx,score,n_TF_rows,n_pass)],by="pair_idx")
stopifnot(nrow(common)==420L,all(is.finite(common$score)))
common[,c("rank_lower","rank_upper"):=.(frank(-score,ties.method="min"),frank(-score,ties.method="max"))]
ligands<-all_scores[,.(score=max(score)),by=.(ligand=ligand_genes)]
ligands<-merge(ligand_common,ligands,by="ligand")
stopifnot(nrow(ligands)==138L,all(is.finite(ligands$score)))
ligands[,c("rank_lower","rank_upper"):=.(frank(-score,ties.method="min"),frank(-score,ties.method="max"))]
fwrite(common,file.path(OUT,"TimeTalk_common_pair_ranking.csv"))
fwrite(ligands,file.path(OUT,"TimeTalk_common_ligand_ranking.csv"))
coverage<-data.table(dataset="GSE48968_BMDC_WT",method="TimeTalk",n_cells=n_cells,n_genes=n_genes,
 n_cell_types=1L,n_mra_tested=n_mra_tested,n_positive_MRA_TFs=nrow(mra_sig),
 n_input_molecular_pairs=nrow(pairs),n_common_pairs=nrow(common),n_common_ligands=nrow(ligands),n_native_TF_rows=nrow(result),
 result_status="completed_BMDC_MRA_autocrine",time_input="observed_collection_hours",
 phenotype="omnibus_stimulus_time_RNA_association;limma_F;FDR0.05",
 molecular_complexes="minimum_of_interpolated_subunit_means_before_native_scaling",
 TF_selection="BMDC_wide_MRA_Pvalue_below_0.05",baseline="shared_zero_hour_cells_per_stimulus",
 score="native_PASS_maximum_signed_SCC_within_condition_then_absolute_threshold_0.8;maximum_strength_across_conditions",seed=CONFIG$seed,
 TimeTalk_version=as.character(packageVersion("TimeTalk")))
fwrite(coverage,file.path(OUT,"coverage.csv"))
print(common[ligand_genes=="IFNB1"&receptor_genes=="IFNAR1_IFNAR2"])
print(ligands[ligand=="IFNB1"])
print(coverage)
