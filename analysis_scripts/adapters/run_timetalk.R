
# path bootstrap (relative; set LLR_ASSETS or LLR_PIPELINE to override)
.llr_cmd <- commandArgs(trailingOnly = FALSE)
.llr_file <- sub("^--file=", "", grep("^--file=", .llr_cmd, value = TRUE)[1])
if (!length(.llr_file) || is.na(.llr_file) || !nzchar(.llr_file)) {
  stop("Run this script via Rscript so --file= is set")
}
.llr_here <- dirname(normalizePath(.llr_file, winslash = "/", mustWork = TRUE))
source(file.path(.llr_here, "../config/llr_bootstrap.R"), local = FALSE)

set.seed(12345)
Sys.setenv(OMP_NUM_THREADS = parallel::detectCores())
data.table::setDTthreads(parallel::detectCores())
llr_source_config()

CONFIG_TIMETALK <- list(
  assay = "RNA",
  minimal_branch_len = 5L,
  cluster_k = 100L,
  cores = parallel::detectCores(),
  q_cutoff = 0.05,
  winsz = 0.1,
  num_pts = 200L,
  scc_cutoff = 0.2,
  granger_cutoff = 1e-2,
  min_pair_cells = 10L,
  lags = 1L
)

library(Matrix)
library(data.table)
library(Seurat)
library(monocle3)
library(TimeTalk)
library(RTN)
source(file.path(CONFIG$dir$algorithms, "TimeTalk", "R", "analysis.R"))
library(dorothea)
library(dplyr)
library(tidyr)
library(tibble)
library(CellChat)
library(SeuratWrappers)
library(snow)
if (requireNamespace("future", quietly = TRUE)) {
  options(future.globals.maxSize = 8 * 1024^3)
}
library(future.apply)
llr_source_rel("adapters/_mem_guard.R")

pick_col <- function(m, candidates) {
  for (cn in candidates) {
    if (!cn %in% names(m)) next
    v <- as.character(m[[cn]])
    v[is.na(v)] <- ""
    if (length(unique(v[v != ""])) >= 2L) return(cn)
  }
  NULL
}

pick_timepoint_col <- function(m) {
  pick_col(m, c("timepoint", "time_point", "Timepoint", "time_info", "stage", "theiler",
                "age(days)", "age_days", "age"))
}

align_meta_cells <- function(meta, cells) {
  if ("cell_id" %in% names(meta)) {
    bc <- as.character(meta$cell_id)
  } else {
    cols <- names(meta)
    if (!nzchar(cols[1L])) meta[[1L]] <- as.character(meta[[1L]])
    scores <- vapply(cols, function(cn) sum(cells %in% as.character(meta[[cn]])), integer(1L))
    best <- cols[which.max(scores)]
    if (max(scores) < length(cells) * 0.9) {
      stop("TimeTalk: cells not found in meta barcode column")
    }
    bc <- as.character(meta[[best]])
  }
  meta$cell_id <- bc
  idx <- match(cells, meta$cell_id)
  if (any(is.na(idx))) {
    stop(sprintf("TimeTalk: %d cells not in meta", sum(is.na(idx))))
  }
  meta[idx]
}

norm_mra_block <- function(mra, gene_use, ct_name) {
  if (is.null(mra) || !is.data.frame(mra) || nrow(mra) == 0) return(NULL)
  reg_col <- intersect(c("Regulon", "regulon", "Regulator", "TF"), names(mra))[1]
  pval_col <- intersect(c("Pvalue", "pvalue", "P.Value", "p.value"), names(mra))[1]
  if (is.na(reg_col) || is.na(pval_col)) return(NULL)
  out <- data.frame(
    Regulon = as.character(mra[[reg_col]]),
    Pvalue = as.numeric(mra[[pval_col]]),
    group = ct_name,
    stringsAsFactors = FALSE
  )
  out <- out[!is.na(out$Pvalue) & out$Pvalue < 0.05, , drop = FALSE]
  if (!length(gene_use)) return(out)
  out[out$Regulon %in% gene_use, , drop = FALSE]
}

args <- commandArgs(trailingOnly = TRUE)
ds <- args[1]
outdir <- args[2]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
marker <- file.path(outdir, "timetalk_done.txt")

# Species-aware databases: mouse datasets must use the mouse DoRothEA regulons and
# the mouse CellChat L-R interactions, not the human ones. Using the wrong species
# database mismatches symbols and removes most available pairs.
sp <- CONFIG$dataset_species[[ds]]; if (is.null(sp)) sp <- "human"
if (identical(sp, "mouse")) {
  utils::data("dorothea_mm", package = "dorothea"); dorothea_tbl <- get("dorothea_mm")
  cc_db <- CellChat::CellChatDB.mouse
} else {
  utils::data("dorothea_hs", package = "dorothea"); dorothea_tbl <- get("dorothea_hs")
  cc_db <- CellChat::CellChatDB.human
}

rd <- file.path(CONFIG$dir$r_export, ds)
if (!file.exists(file.path(rd, "expr.mtx"))) {
  sys::exec_wait(
    CONFIG$exe$rscript,
    c(file.path(llr_pipeline_dir(), "adapters/00_export_R.R"), ds)
  )
}
cells <- readLines(file.path(rd, "cells.txt"))
genes <- readLines(file.path(rd, "genes.txt"))
meta <- align_meta_cells(fread(file.path(rd, "meta.csv")), cells)
tp_col <- pick_timepoint_col(meta)
if (is.null(tp_col)) {
  stop("TimeTalk: no timepoint column with >= 2 distinct values")
}
ct_col <- pick_col(meta, c("cell_type", "celltype", "clusters", "cluster", "seurat_clusters", "state_info"))
if (is.null(ct_col)) stop("TimeTalk: no cell type column with >= 2 groups")
ct <- as.character(meta[[ct_col]])

cnt_s <- readMM(file.path(rd, "expr.mtx"))
n_cells <- ncol(cnt_s)
n_genes <- nrow(cnt_s)
if (n_genes > 3000L) {
  mu <- Matrix::rowMeans(cnt_s)
  va <- Matrix::rowMeans(cnt_s^2) - mu^2
  va[!is.finite(va)] <- 0
  keep <- order(va, decreasing = TRUE)[seq_len(min(2000L, n_genes))]
  cnt_s <- cnt_s[keep, , drop = FALSE]
  genes <- genes[keep]
  n_genes <- length(genes)
}
mem_require_under_gb(
  mem_dense_gb(n_cells, n_genes, copies = 2),
  what = paste0("TimeTalk/", ds)
)
peak_main_gb <- mem_dense_gb(n_cells, n_genes, copies = 2)
# snow duplicates the full Seurat object per worker; stay single-threaded on panel sizes.
tt_cores <- 1L
Sys.setenv(OMP_NUM_THREADS = tt_cores)
data.table::setDTthreads(tt_cores)
cat(sprintf("TimeTalk/%s: %d cells x %d genes, cores=%d\n", ds, n_cells, n_genes, tt_cores))
cnt <- as.matrix(cnt_s)
rm(cnt_s)
gc()
rownames(cnt) <- genes
colnames(cnt) <- cells

seu <- CreateSeuratObject(counts = cnt)
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu, nfeatures = 2000)
seu <- ScaleData(seu)
seu <- RunPCA(seu)
seu <- RunUMAP(seu, dims = 1:30, umap.method = "uwot")
seu <- AddMetaData(seu, metadata = data.frame(
  CellType = ct,
  orig.ident = ds,
  row.names = cells
))
Idents(seu) <- "CellType"

cds <- RunMonocle3(tmp.seu = seu, tmp.assay = CONFIG_TIMETALK$assay,
  tmp.minimal_branch_len = CONFIG_TIMETALK$minimal_branch_len,
  tmp.cluster_k = CONFIG_TIMETALK$cluster_k)
time_order <- sort(unique(as.character(meta[[tp_col]])))
time_order <- time_order[!is.na(time_order) & time_order != ""]
root_cells <- cells[as.character(meta[[tp_col]]) == time_order[1L]]
if (!length(root_cells)) stop("TimeTalk: no root cells at first timepoint")
cds <- monocle3::order_cells(cds, root_cells = root_cells)
seu$Pseudotime <- monocle3::pseudotime(cds)

phenotype <- GetVariableGeneInTrajectory(tmp.cds = cds,
  tmp.cores = tt_cores,
  tmp.neighbor_graph = "principal_graph",
  tmp.q_value_cutoff = CONFIG_TIMETALK$q_cutoff)

gene_tf <- intersect(unique(dorothea_tbl$tf), rownames(cnt))

meta_ct <- data.table(CellType = ct, timepoint = as.character(meta[[tp_col]]))
ct_time <- meta_ct[, .N, by = .(CellType, timepoint)]
# Full temporal cell-type matrix (no top-2 restriction): every early-timepoint
# cell type as sender x every late-timepoint cell type as receiver, preserving
# TimeTalk's early->late direction.
sender_cts <- ct_time[timepoint == time_order[1], sort(unique(as.character(CellType)))]
receiver_cts <- ct_time[timepoint == time_order[length(time_order)], sort(unique(as.character(CellType)))]
if (length(sender_cts) < 1L || length(receiver_cts) < 1L) {
  stop("TimeTalk: no sender/receiver cell types at first/last timepoint")
}
involved_ct <- union(sender_cts, receiver_cts)

rtna_list <- stats::setNames(
  lapply(involved_ct, function(ct_name) {
    tryCatch(
      GetCellTypeSpecificTRN(tmp.seu = seu, tmp.ident = ct_name,
        gene.TF = gene_tf, tmp.phenotype = phenotype),
      error = function(e) {
        message("TimeTalk TRN skip ", ct_name, ": ", conditionMessage(e))
        NULL
      }
    )
  }),
  involved_ct
)
rtna_list <- rtna_list[!vapply(rtna_list, is.null, logical(1L))]
if (!length(rtna_list)) stop("TimeTalk: no cell-type TRN could be built")

all_markers <- FindAllMarkers(seu, assay = "RNA", only.pos = TRUE)

mra_list <- lapply(names(rtna_list), function(ct_name) {
  rtna <- rtna_list[[ct_name]]
  mra <- RTN::tna.get(rtna, what = "mra", ntop = -1)
  gene_use <- all_markers %>%
    dplyr::filter(cluster == ct_name, p_val_adj < 0.05) %>%
    dplyr::pull(gene)
  norm_mra_block(mra, gene_use, ct_name)
})
mra_list <- mra_list[!vapply(mra_list, is.null, logical(1L))]
if (!length(mra_list)) stop("TimeTalk: no MRA results after TRN filtering")
mra_res <- data.table::rbindlist(mra_list, use.names = TRUE)
setDF(mra_res)

lr_inter <- cc_db$interaction
lr_list <- lapply(seq_len(nrow(lr_inter)), function(i) {
  l <- lr_inter$ligand[i]
  rs <- strsplit(lr_inter$receptor[i], "_")[[1]]
  lapply(rs, function(r) {
    data.frame(lr_pair = paste(l, r, sep = "_"),
      ligand_gene_symbol = l, receptor_gene_symbol = r,
      stringsAsFactors = FALSE)
  })
})
lr_df <- unique(data.table::rbindlist(unlist(lr_list, recursive = FALSE), use.names = TRUE))
setDF(lr_df)

# Per-pair disk cache (resumable, identical result): skip already-computed pairs.
pair_cache <- file.path(outdir, "timetalk_pairs")
dir.create(pair_cache, recursive = TRUE, showWarnings = FALSE)
# RunTimeTalk may use future multisession; force sequential + large export cap.
if (requireNamespace("future", quietly = TRUE)) {
  options(future.globals.maxSize = 8 * 1024^3)
  future::plan(future::sequential)
}
all_res <- list()
for (s_ct in sender_cts) {
  for (r_ct in receiver_cts) {
    if (s_ct == r_ct) next
    n_send <- sum(ct == s_ct & as.character(meta[[tp_col]]) == time_order[1L])
    n_recv <- sum(ct == r_ct & as.character(meta[[tp_col]]) == time_order[length(time_order)])
    if (n_send < CONFIG_TIMETALK$min_pair_cells || n_recv < CONFIG_TIMETALK$min_pair_cells) next
    pair_csv <- file.path(pair_cache, paste0(make.names(s_ct), "__", make.names(r_ct), ".csv"))
    if (file.exists(pair_csv) && file.info(pair_csv)$size > 0) {
      all_res[[length(all_res) + 1L]] <- fread(pair_csv)
      next
    }
    res <- tryCatch(
      RunTimeTalk(
        tmp.cds = cds, tmp.seu = seu, tmp.orig.ident = ds,
        tmp.ident.1 = s_ct, tmp.ident.2 = r_ct,
        LRpairs.df = lr_df, tmp.mra.res = mra_res,
        tmp.winsz = CONFIG_TIMETALK$winsz, tmp.lags = CONFIG_TIMETALK$lags,
        numPts = CONFIG_TIMETALK$num_pts,
        tmp.SCC.cutoff = CONFIG_TIMETALK$scc_cutoff,
        tmp.granger.cutoff = CONFIG_TIMETALK$granger_cutoff,
        tmp.cores = tt_cores),
      error = function(e) {
        message("TimeTalk pair skip ", s_ct, " -> ", r_ct, ": ", conditionMessage(e))
        NULL
      }
    )
    if (is.null(res)) next
    rdt <- tryCatch(as.data.table(res), error = function(e) NULL)
    if (is.null(rdt) || nrow(rdt) == 0L) next
    rdt[, sender := s_ct]
    rdt[, receiver := r_ct]
    fwrite(rdt, pair_csv)
    all_res[[length(all_res) + 1L]] <- rdt
  }
}
if (length(all_res) == 0L) {
  stop("TimeTalk: no sender->receiver pair produced results")
}
result <- data.table::rbindlist(all_res, use.names = TRUE, fill = TRUE)
fwrite(result, file.path(outdir, "timetalk_result.csv"))
writeLines(as.character(Sys.time()), marker)
cat("OK:", length(all_res), "sender->receiver pairs\n")
