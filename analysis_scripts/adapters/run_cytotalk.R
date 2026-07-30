#!/usr/bin/env Rscript
# CytoTalk intercellular pathway inference between two cell types.
suppressPackageStartupMessages(library(data.table))

# path bootstrap (relative; set LLR_ASSETS or LLR_PIPELINE to override)
.llr_cmd <- commandArgs(trailingOnly = FALSE)
.llr_file <- sub("^--file=", "", grep("^--file=", .llr_cmd, value = TRUE)[1])
if (!length(.llr_file) || is.na(.llr_file) || !nzchar(.llr_file)) {
  stop("Run this script via Rscript so --file= is set")
}
.llr_here <- dirname(normalizePath(.llr_file, winslash = "/", mustWork = TRUE))
source(file.path(.llr_here, "../config/llr_bootstrap.R"), local = FALSE)

llr_source_config()

args <- commandArgs(trailingOnly = TRUE)
ds <- args[1]
outdir <- args[2]
if (is.na(ds)) stop("Usage: Rscript run_cytotalk.R <dataset> <outdir>")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
marker <- file.path(outdir, "cytotalk_done.txt")

vendor <- file.path(CONFIG$dir$algorithms, "CytoTalk")
if (!dir.exists(vendor)) stop("CytoTalk vendor missing: ", vendor)
for (pkg in c("parmigene", "DiagrammeRsvg", "vroom", "Matrix")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("CytoTalk dependency missing: ", pkg)
  }
}
suppressPackageStartupMessages(devtools::load_all(vendor, quiet = TRUE))
reticulate::use_python(CONFIG$exe$python, required = TRUE)
if (dir.exists(CONFIG$exe$mingw_bin)) {
  reticulate::py_run_string(sprintf(
    "import os\nos.add_dll_directory(r'%s')",
    CONFIG$exe$mingw_bin
  ))
}

sp <- CONFIG$dataset_species[[ds]]
if (is.null(sp)) sp <- "human"
if (!sp %in% c("human", "mouse")) {
  stop("CytoTalk requires human or mouse species, got: ", sp)
}
pcg <- if (sp == "mouse") CytoTalk::pcg_mouse else CytoTalk::pcg_human
lrp <- if (sp == "mouse") CytoTalk::lrp_mouse else CytoTalk::lrp_human

rd <- file.path(CONFIG$dir$r_export, ds)
if (!file.exists(file.path(rd, "expr_lognorm.mtx"))) {
  sys::exec_wait(
    CONFIG$exe$rscript,
    c(file.path(llr_pipeline_dir(), "adapters/00_export_R.R"), ds)
  )
}
expr <- Matrix::readMM(file.path(rd, "expr_lognorm.mtx"))
genes <- readLines(file.path(rd, "genes.txt"))
cells <- readLines(file.path(rd, "cells.txt"))
rownames(expr) <- genes
colnames(expr) <- cells[seq_len(ncol(expr))]
meta <- fread(file.path(rd, "meta.csv"))
if (!"cell_id" %in% names(meta)) meta$cell_id <- as.character(meta[[1L]])
idx <- match(cells, meta$cell_id)
if (any(is.na(idx))) stop("CytoTalk: cells not found in meta barcode column")
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
if (is.null(ct_col)) {
  stop("CytoTalk: no cell type column found")
}
ct <- as.character(meta[[ct_col]])
ct[is.na(ct)] <- "unknown"
tab <- sort(table(ct), decreasing = TRUE)
all_types <- names(tab)
if (length(all_types) < 2L) {
  stop("CytoTalk: fewer than 2 cell types")
}

lst <- list(
  mat = expr,
  cell_types = ct
)

cyto_cores <- 1L

# Full cell-type pair matrix (no top-2 restriction): every unordered pair.
# Each pair runs the official run_cytotalk with vendor ntrial=1000. Pair-level
# errors stop the dataset; completed pairs remain available in the RDS cache.
summ <- list()
for (ia in seq_len(length(all_types) - 1L)) {
  for (ib in seq.int(ia + 1L, length(all_types))) {
    a <- all_types[ia]
    b <- all_types[ib]
    pair_dir <- file.path(outdir, paste0(make.names(a), "__", make.names(b)))
    dir.create(pair_dir, recursive = TRUE, showWarnings = FALSE)
    # Per-pair disk cache (resumable, identical result): skip already-computed pairs.
    res_rds <- file.path(pair_dir, "cytotalk_result.rds")
    net_txt <- file.path(pair_dir, "FinalNetwork.txt")
    if (file.exists(res_rds) && file.info(res_rds)$size > 100) {
      summ[[length(summ) + 1L]] <- data.table(
        dataset = ds, cell_type_a = a, cell_type_b = b, species = sp, status = "CACHED"
      )
      next
    }
    if (file.exists(net_txt) && file.info(net_txt)$size > 0) {
      summ[[length(summ) + 1L]] <- data.table(
        dataset = ds, cell_type_a = a, cell_type_b = b, species = sp, status = "OK"
      )
      next
    }
    res <- tryCatch(
      CytoTalk::run_cytotalk(
        lst_scrna = lst,
        cell_type_a = a,
        cell_type_b = b,
        cutoff_a = 0.2,
        cutoff_b = 0.2,
        pcg = pcg,
        lrp = lrp,
        ntrial = 1000L,
        cores = cyto_cores,
        echo = FALSE,
        dir_out = pair_dir
      ),
      error = function(e) e
    )
    if (inherits(res, "error")) {
      msg <- conditionMessage(res)
      if (file.exists(net_txt) && file.info(net_txt)$size > 0) {
        summ[[length(summ) + 1L]] <- data.table(
          dataset = ds, cell_type_a = a, cell_type_b = b, species = sp,
          status = "OK", note = paste0("partial:", msg)
        )
      } else {
        summ[[length(summ) + 1L]] <- data.table(
          dataset = ds, cell_type_a = a, cell_type_b = b, species = sp,
          status = "EMPTY", note = msg
        )
      }
      next
    }
    saveRDS(res, file.path(pair_dir, "cytotalk_result.rds"))
    if (!is.null(res$pcst$final_network)) {
      fwrite(as.data.table(res$pcst$final_network), file.path(pair_dir, "FinalNetwork.txt"), sep = "\t")
    }
    if (!is.null(res$pcst$ks_test_pval)) {
      fwrite(as.data.table(res$pcst$ks_test_pval), file.path(pair_dir, "PCSTScores.txt"), sep = "\t")
    }
    if (!is.null(res$pathways$df_pval)) {
      fwrite(as.data.table(res$pathways$df_pval), file.path(pair_dir, "PathwayAnalysis.txt"), sep = "\t")
    }
    summ[[length(summ) + 1L]] <- data.table(
      dataset = ds, cell_type_a = a, cell_type_b = b, species = sp, status = "OK"
    )
  }
}
summary_dt <- data.table::rbindlist(summ, use.names = TRUE, fill = TRUE)
fwrite(summary_dt, file.path(outdir, "cytotalk_summary.csv"))
n_ok <- sum(summary_dt$status %in% c("OK", "CACHED"))
if (n_ok == 0L) {
  cat("FAIL: no cell-type pair produced CytoTalk output\n", file = stderr())
  quit(save = "no", status = 1)
}
writeLines(as.character(Sys.time()), marker)
cat("OK:", n_ok, "of", nrow(summary_dt), "cell-type pairs\n")
