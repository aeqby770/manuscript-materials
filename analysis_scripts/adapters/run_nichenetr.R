#!/usr/bin/env Rscript
# NicheNet ligand activity prediction (sender -> receiver cell groups).
suppressPackageStartupMessages({
  library(Matrix)

# path bootstrap (relative; set LLR_ASSETS or LLR_PIPELINE to override)
.llr_cmd <- commandArgs(trailingOnly = FALSE)
.llr_file <- sub("^--file=", "", grep("^--file=", .llr_cmd, value = TRUE)[1])
if (!length(.llr_file) || is.na(.llr_file) || !nzchar(.llr_file)) {
  stop("Run this script via Rscript so --file= is set")
}
.llr_here <- dirname(normalizePath(.llr_file, winslash = "/", mustWork = TRUE))
source(file.path(.llr_here, "../config/llr_bootstrap.R"), local = FALSE)

  library(data.table)
  library(dplyr)
  library(nichenetr)
})
llr_source_config()
llr_source_rel("adapters/_mem_guard.R")

args <- commandArgs(trailingOnly = TRUE)
ds <- args[1]
out_file <- args[2]
if (is.na(ds)) stop("Usage: Rscript run_nichenetr.R <dataset> <output.csv>")

sp <- CONFIG$dataset_species[[ds]]
if (is.null(sp)) sp <- "human"
if (!sp %in% c("human", "mouse")) {
  stop("NicheNet requires human or mouse species, got: ", sp)
}

rd <- file.path(CONFIG$dir$r_export, ds)
if (!file.exists(file.path(rd, "expr_lognorm.mtx"))) {
  sys::exec_wait(
    CONFIG$exe$rscript,
    c(file.path(llr_pipeline_dir(), "adapters/00_export_R.R"), ds)
  )
}
expr_s <- readMM(file.path(rd, "expr_lognorm.mtx"))
genes <- readLines(file.path(rd, "genes.txt"))
cells <- readLines(file.path(rd, "cells.txt"))
rownames(expr_s) <- genes
colnames(expr_s) <- cells[seq_len(ncol(expr_s))]
meta <- fread(file.path(rd, "meta.csv"))
if (!"cell_id" %in% names(meta)) meta$cell_id <- as.character(meta[[1L]])
idx <- match(cells, meta$cell_id)
if (any(is.na(idx))) stop("NicheNet: cells not found in meta barcode column")
meta <- meta[idx]
gene_var <- Matrix::rowMeans(expr_s^2) - Matrix::rowMeans(expr_s)^2
keep <- which(is.finite(gene_var) & gene_var > 1e-10)
top_n <- min(2000L, length(keep))
keep <- keep[order(gene_var[keep], decreasing = TRUE)[seq_len(top_n)]]
mem_require_under_gb(
  mem_dense_gb(ncol(expr_s), length(keep), copies = 2),
  what = paste0("NicheNet/", ds)
)
expr <- as.matrix(expr_s[keep, , drop = FALSE])
rm(expr_s)
gc()

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
if (is.null(ct_col)) {
  stop("NicheNet: no cluster column found")
}
ct <- as.character(meta[[ct_col]])
ct[is.na(ct)] <- "unknown"
tab <- sort(table(ct), decreasing = TRUE)
if (length(tab) < 2L) {
  stop("NicheNet: fewer than 2 cell groups")
}
cell_groups <- names(tab)

if (sp == "mouse") {
  data("geneinfo_2022", package = "nichenetr", envir = environment())
  mapped <- nichenetr::convert_mouse_to_human_symbols(rownames(expr), version = 2)
  keep <- !is.na(mapped)
  expr <- rowsum(expr[keep, , drop = FALSE], group = mapped[keep], reorder = FALSE)
}

geneset_oi <- rownames(expr)[Matrix::rowSums(expr > 0) >= 10]
if (length(geneset_oi) < 50L) {
  stop("NicheNet: too few expressed genes (", length(geneset_oi), " < 50)")
}

cache_dir <- file.path(CONFIG$dir$cache, "nichenetr")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
fetch_rds <- function(name, url) {
  path <- file.path(cache_dir, name)
  if (!file.exists(path)) download.file(url, path, mode = "wb", quiet = TRUE)
  readRDS(path)
}
ligand_target_matrix <- fetch_rds(
  "ligand_target_matrix_nsga2r_final.rds",
  "https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds"
)
lr_network <- fetch_rds(
  "lr_network_human_21122021.rds",
  "https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"
)
weighted_networks <- fetch_rds(
  "weighted_networks_nsga2r_final.rds",
  "https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final.rds"
)

ligands <- unique(lr_network$from)
receptors <- unique(lr_network$to)

# Full sender x receiver matrix (no top-2 restriction). Each ordered pair runs
# the standard NicheNet ligand-activity prediction; a pair with too few
# expressed ligands/receptors/targets is skipped, not the whole run.
compute_pair <- function(sender, receiver) {
  sender_cells <- colnames(expr)[ct == sender]
  receiver_cells <- colnames(expr)[ct == receiver]
  geneset_sender <- geneset_oi[Matrix::rowSums(expr[geneset_oi, sender_cells, drop = FALSE] > 0) >= 5]
  geneset_receiver <- geneset_oi[Matrix::rowSums(expr[geneset_oi, receiver_cells, drop = FALSE] > 0) >= 5]
  if (length(geneset_receiver) < 10L) return(NULL)
  receiver_logfc <- Matrix::rowMeans(expr[geneset_receiver, receiver_cells, drop = FALSE]) -
    Matrix::rowMeans(expr[geneset_receiver, sender_cells, drop = FALSE])
  target_genes <- head(names(sort(receiver_logfc[receiver_logfc > 0], decreasing = TRUE)), 250L)
  if (length(target_genes) < 10L) return(NULL)
  expressed_ligands <- ligands[ligands %in% geneset_sender]
  expressed_receptors <- receptors[receptors %in% geneset_receiver]
  if (length(expressed_ligands) < 3L || length(expressed_receptors) < 3L) return(NULL)
  potential_ligands <- unique(lr_network$from[lr_network$from %in% expressed_ligands & lr_network$to %in% expressed_receptors])
  if (length(potential_ligands) < 3L) return(NULL)
  act <- predict_ligand_activities(
    geneset = target_genes,
    background_expressed_genes = geneset_receiver,
    ligand_target_matrix = ligand_target_matrix,
    potential_ligands = potential_ligands
  )
  setDT(act)
  act[, sender := sender]
  act[, receiver := receiver]
  act[, target_gene_count := length(target_genes)]
  act
}
# Per-pair disk cache (resumable, identical result): skip already-computed pairs.
pair_cache <- file.path(dirname(out_file), "nichenetr_pairs")
dir.create(pair_cache, recursive = TRUE, showWarnings = FALSE)
all_act <- list()
for (sender in cell_groups) {
  for (receiver in cell_groups) {
    if (sender == receiver) next
    pair_csv <- file.path(pair_cache, paste0(make.names(sender), "__", make.names(receiver), ".csv"))
    if (file.exists(pair_csv) && file.info(pair_csv)$size > 0) {
      all_act[[length(all_act) + 1L]] <- fread(pair_csv)
      next
    }
    res <- compute_pair(sender, receiver)
    if (!is.null(res)) {
      fwrite(res, pair_csv)
      all_act[[length(all_act) + 1L]] <- res
    }
  }
}
if (length(all_act) == 0L) {
  cat("FAIL: no sender->receiver pair produced ligand activities\n", file = stderr())
  quit(save = "no", status = 1)
}
activities <- rbindlist(all_act, use.names = TRUE, fill = TRUE)
activities[, dataset := ds]
fwrite(activities, out_file)
cat("OK:", nrow(activities), "ligand activities across", length(cell_groups), "groups ->", out_file, "\n")
