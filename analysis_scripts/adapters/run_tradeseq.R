suppressPackageStartupMessages({library(Matrix);library(data.table);library(tradeSeq)})

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
smoke_root <- Sys.getenv("LLR_SMOKE_ROOT", unset = NA_character_)
if (!is.na(smoke_root) && nzchar(smoke_root)) CONFIG$dir$algorithm_out <- file.path(smoke_root, "by_algorithm")
args<-commandArgs(trailingOnly=TRUE);ds<-args[1];outdir<-args[2]
if(is.na(ds))stop("Usage: Rscript run_tradeseq.R <dataset> <outdir>")
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
rd<-file.path(CONFIG$dir$r_export,ds)
if(!file.exists(file.path(rd,"expr.mtx")))
  sys::exec_wait(CONFIG$exe$rscript,
    c(file.path(llr_pipeline_dir(), "adapters/00_export_R.R"),ds))
cnt<-as(readMM(file.path(rd,"expr.mtx")),"CsparseMatrix")
meta<-fread(file.path(rd,"meta.csv"))
genes<-readLines(file.path(rd,"genes.txt"))
cells<-readLines(file.path(rd,"cells.txt"))
rownames(cnt)<-genes
colnames(cnt)<-cells
# Speed-only QC (tradeSeq vignette): drop lowly-expressed genes before fitGAM.
# Genes detected in too few cells cannot yield a reliable NB-GAM fit, so this
# is standard preprocessing and does not change the model/assumptions.
keep_min<-max(10L, ceiling(0.01*ncol(cnt)))
gene_detect<-rowSums(cnt>0)
n_before<-nrow(cnt)
cnt<-cnt[gene_detect>=keep_min, , drop=FALSE]
cat(sprintf("tradeSeq: kept %d/%d genes expressed in >= %d cells\n", nrow(cnt), n_before, keep_min))
# Large panels: top-2000 HVG by lognorm variance (same ranking as run_slingshot.R).
# Small panels (<=5000 expressed genes): keep all expressed genes (closer to official vignette).
HVG_TRIGGER<-5000L; HVG_TOP<-2000L
n_after_expr<-nrow(cnt)
if(n_after_expr>HVG_TRIGGER){
  lognorm<-as(readMM(file.path(rd,"expr_lognorm.mtx")),"CsparseMatrix")
  rownames(lognorm)<-genes; colnames(lognorm)<-cells
  lognorm<-lognorm[rownames(cnt), , drop=FALSE]
  gene_var<-Matrix::rowMeans(lognorm^2)-Matrix::rowMeans(lognorm)^2
  hvg_keep<-which(is.finite(gene_var)&gene_var>1e-10)
  top_n<-min(HVG_TOP, length(hvg_keep))
  hvg_keep<-hvg_keep[order(gene_var[hvg_keep], decreasing=TRUE)[seq_len(top_n)]]
  cnt<-cnt[hvg_keep, , drop=FALSE]
  cat(sprintf("tradeSeq: HVG %d/%d genes (top %d lognorm variance, slingshot-aligned)\n",
              nrow(cnt), n_after_expr, top_n))
}
# fitGAM densifies the (filtered) counts and copies slices across parallel workers.
mem_require_under_gb(mem_dense_gb(ncol(cnt), nrow(cnt), copies = 2), what = paste0("tradeSeq/", ds))
# Need pseudotime + cellWeights — reuse slingshot output if available
pt_file<-file.path(CONFIG$dir$algorithm_out,"slingshot",ds,"slingshot_pseudotime.csv")
cw_file<-file.path(CONFIG$dir$algorithm_out,"slingshot",ds,"slingshot_weights.csv")
if(!file.exists(pt_file) || !file.exists(cw_file)) stop("tradeSeq requires slingshot pseudotime/weights for ", ds)
pt<-as.matrix(fread(pt_file, header=TRUE))
cw<-as.matrix(fread(cw_file, header=TRUE))
if(nrow(pt) != length(cells) || nrow(cw) != length(cells)) {
  stop(sprintf("tradeSeq: slingshot nrow mismatch for %s (pt=%d cw=%d cells=%d)", ds, nrow(pt), nrow(cw), length(cells)))
}
rownames(pt)<-cells
rownames(cw)<-cells
sce<-fitGAM(counts=cnt,pseudotime=pt,cellWeights=cw,nknots=3,verbose=FALSE,parallel=TRUE,BPPARAM=BiocParallel::SnowParam(workers=CONFIG$n_threads,progressbar=FALSE))
assoc<-associationTest(sce)
saveRDS(sce,file.path(outdir,"tradeseq_sce.rds"))
fwrite(as.data.table(assoc, keep.rownames="gene"),file.path(outdir,"tradeseq_association_test.csv"))
cat("OK\n")
