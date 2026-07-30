suppressPackageStartupMessages({library(Matrix);library(data.table);library(CellChat)})

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
args <- commandArgs(trailingOnly=TRUE); ds <- args[1]; outdir <- args[2]
if (is.na(ds)) stop("Usage: Rscript run_cellchat.R <dataset> <outdir>")
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
rd <- file.path(CONFIG$dir$r_export, ds)
if (!file.exists(file.path(rd, "expr_lognorm.mtx")))
  sys::exec_wait(CONFIG$exe$rscript,
    c(file.path(llr_pipeline_dir(), "adapters/00_export_R.R"), ds))
expr_s <- readMM(file.path(rd, "expr_lognorm.mtx"))
genes <- readLines(file.path(rd, "genes.txt"))
cells <- readLines(file.path(rd, "cells.txt"))
rownames(expr_s) <- genes
colnames(expr_s) <- cells[seq_len(ncol(expr_s))]
cells <- colnames(expr_s)
meta <- fread(file.path(rd, "meta.csv"))
if (!"cell_id" %in% names(meta)) meta$cell_id <- as.character(meta[[1L]])
idx <- match(cells, meta$cell_id)
if (any(is.na(idx))) stop("CellChat: cells not found in meta barcode column")
meta <- meta[idx]
pick_col <- function(m) {
  for (cn in c("cell_type", "celltype", "CellType", "clusters", "cluster", "seurat_clusters", "state_info", "cell_ontology_class")) {
    if (!cn %in% names(m)) next
    v <- as.character(m[[cn]])
    v[is.na(v)] <- ""
    if (length(unique(v[v != ""])) >= 2L) return(cn)
  }
  NULL
}
ct_col <- pick_col(meta)
if (is.null(ct_col)) stop("CellChat: no cluster column with >= 2 groups")
sp <- CONFIG$dataset_species[[ds]]
if (is.null(sp)) sp <- "human"
if (!sp %in% c("human", "mouse")) {
  stop("CellChat: species '", sp, "' has no official CellChatDB; use human/mouse datasets only")
}
db <- if (sp == "mouse") CellChatDB.mouse else CellChatDB.human
lr_genes <- unique(unlist(strsplit(unique(c(db$interaction$ligand, db$interaction$receptor)), "_")))
lr_idx <- which(genes %in% lr_genes)
n_cells <- ncol(expr_s)
n_genes <- nrow(expr_s)
if (as.numeric(n_cells) * as.numeric(n_genes) <= 3e7) {
  keep <- seq_len(n_genes)
} else {
  gene_var <- Matrix::rowMeans(expr_s^2) - Matrix::rowMeans(expr_s)^2
  vkeep <- which(is.finite(gene_var) & gene_var > 1e-10)
  top_n <- min(2000L, length(vkeep))
  keep <- vkeep[order(gene_var[vkeep], decreasing = TRUE)[seq_len(top_n)]]
  keep <- unique(c(keep, lr_idx))
}
mem_require_under_gb(
  mem_dense_gb(n_cells, length(keep), copies = 2),
  what = paste0("CellChat/", ds)
)
expr <- as.matrix(expr_s[keep, , drop = FALSE])
rm(expr_s)
gc()
ct <- as.character(meta[[ct_col]])
ct[is.na(ct) | ct == ""] <- "unknown"
# CellChat rejects labels containing "0" (e.g. "0", "10", "20"). Prefix numeric labels.
if (any(grepl("^\\d+$", ct))) ct <- paste0("C", ct)
meta2 <- data.frame(labels = ct, row.names = cells)
cc <- createCellChat(object = expr, meta = meta2, group.by = "labels")
cc <- setIdent(cc, ident.use = "labels")
cc@DB <- db
cc <- subsetData(cc)
cc <- identifyOverExpressedGenes(cc)
cc <- identifyOverExpressedInteractions(cc)
# projectData removed in CellChat v2.2+; PPI projection is optional smoothing.
# Skipping it does not affect computeCommunProb downstream.
if (exists("projectData", where = asNamespace("CellChat"))) {
  cc <- projectData(cc, if (sp == "mouse") PPI.mouse else PPI.human)
}
cc <- computeCommunProb(cc, type = "truncatedMean", trim = 0.1)
cc <- filterCommunication(cc, min.cells = 10)
cc <- computeCommunProbPathway(cc)
cc <- aggregateNet(cc)
saveRDS(cc, file.path(outdir, "cellchat.rds"))
cat("OK:", length(cc@netP$pathways), "pathways,", nrow(cc@LR$LRsig), "sig LR pairs\n")
