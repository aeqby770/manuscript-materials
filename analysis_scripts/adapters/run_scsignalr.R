suppressPackageStartupMessages({library(Matrix);library(data.table);library(SingleCellSignalR);library(SingleCellExperiment)})

# path bootstrap (relative; set LLR_ASSETS or LLR_PIPELINE to override)
.llr_cmd <- commandArgs(trailingOnly = FALSE)
.llr_file <- sub("^--file=", "", grep("^--file=", .llr_cmd, value = TRUE)[1])
if (!length(.llr_file) || is.na(.llr_file) || !nzchar(.llr_file)) {
  stop("Run this script via Rscript so --file= is set")
}
.llr_here <- dirname(normalizePath(.llr_file, winslash = "/", mustWork = TRUE))
source(file.path(.llr_here, "../config/llr_bootstrap.R"), local = FALSE)

llr_source_config()
llr_source_rel("adapters/_mem_guard.R")
args<-commandArgs(trailingOnly=TRUE);ds<-args[1];outdir<-args[2]
if(is.na(ds))stop("Usage: Rscript run_scsignalr.R <dataset> <outdir>")
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
rd<-file.path(CONFIG$dir$r_export,ds)
if(!file.exists(file.path(rd,"expr_lognorm.mtx")))
  sys::exec_wait(CONFIG$exe$rscript,
    c(file.path(llr_pipeline_dir(), "adapters/00_export_R.R"),ds))
expr_s<-readMM(file.path(rd,"expr_lognorm.mtx"))
meta<-fread(file.path(rd,"meta.csv"))
genes<-readLines(file.path(rd,"genes.txt"))
cells<-readLines(file.path(rd,"cells.txt"))
rownames(expr_s)<-genes
colnames(expr_s)<-cells[seq_len(ncol(expr_s))]
gene_var<-Matrix::rowMeans(expr_s^2)-Matrix::rowMeans(expr_s)^2
keep<-which(is.finite(gene_var)&gene_var>1e-10)
top_n<-min(2000L,length(keep))
keep<-keep[order(gene_var[keep],decreasing=TRUE)[seq_len(top_n)]]
mem_require_under_gb(
  mem_dense_gb(ncol(expr_s),length(keep),copies=2),
  what=paste0("SingleCellSignalR/",ds)
)
expr<-as.matrix(expr_s[keep,,drop=FALSE])
rm(expr_s)
gc()
ct_col<-NA_character_
for(cn in c("clusters","cluster","cell_type","celltype","seurat_clusters","state_info","ct","ct_subtype","paul15_clusters","louvain","leiden","cell_state")){
  if(!cn%in%names(meta))next
  v<-as.character(meta[[cn]]);v[is.na(v)]<-""
  if(length(unique(v[v!=""]))>=2L){ct_col<-cn;break}
}
if(is.na(ct_col))stop("SingleCellSignalR: no cluster column with >= 2 groups")
ct<-meta[[ct_col]]
pop<-as.character(ct)
sp<-CONFIG$dataset_species[[ds]]
if(is.null(sp))sp<-"human"
if(!sp%in%c("human","mouse")){
  stop("SingleCellSignalR: species '",sp,"' is not supported (LRdb is human-only); use human/mouse datasets only")
}
if(CONFIG$dataset_species[[ds]]=="mouse"){
  data(geneinfo_2022, package="nichenetr", envir=environment())
  mapped<-nichenetr::convert_mouse_to_human_symbols(rownames(expr), version=2)
  keep<-!is.na(mapped)
  expr<-rowsum(expr[keep,,drop=FALSE], group=mapped[keep], reorder=FALSE)
}
scsr<-SCSRNoNet(counts=expr,populations=pop,normalize=FALSE,method="log-only",log.transformed=TRUE,species="hsapiens",min.LR.found=80)
# Run the actual L-R inference. performInferences() infers paracrine interactions
# across ALL ordered population pairs and autocrine within each population
# (vendor defaults: subsample.size=50, n.resample=50, max.pval=0.05, min.logFC=1).
scsr<-performInferences(scsr,autocrine=TRUE,paracrine=TRUE,verbose=FALSE)
saveRDS(scsr,file.path(outdir,"scsignalr_scsr.rds"))
# Flatten paracrine + autocrine interaction tables into one CSV.
para<-scsr@paracrines
auto<-scsr@autocrines
flat<-list()
if(length(para)>0) for(nm in names(para)){d<-para[[nm]];if(!is.null(d)&&nrow(d)>0){d$pair<-nm;d$kind<-"paracrine";flat[[length(flat)+1]]<-d}}
if(length(auto)>0) for(nm in names(auto)){d<-auto[[nm]];if(!is.null(d)&&nrow(d)>0){d$pair<-nm;d$kind<-"autocrine";flat[[length(flat)+1]]<-d}}
if(length(flat)>0){
  res<-data.table::rbindlist(flat,use.names=TRUE,fill=TRUE)
  data.table::fwrite(res,file.path(outdir,"scsignalr_interactions.csv"))
  cat("OK:",ncol(expr),"cells,",nrow(res),"LR interactions,",length(para),"paracrine pairs\n")
}else{
  cat("FAIL: 0 L-R interactions inferred\n", file=stderr())
  quit(save="no", status=1)
}
