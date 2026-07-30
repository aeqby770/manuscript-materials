
# path bootstrap (relative; set LLR_ASSETS or LLR_PIPELINE to override)
.llr_cmd <- commandArgs(trailingOnly = FALSE)
.llr_file <- sub("^--file=", "", grep("^--file=", .llr_cmd, value = TRUE)[1])
if (!length(.llr_file) || is.na(.llr_file) || !nzchar(.llr_file)) {
  stop("Run this script via Rscript so --file= is set")
}
.llr_here <- dirname(normalizePath(.llr_file, winslash = "/", mustWork = TRUE))
source(file.path(.llr_here, "../config/llr_bootstrap.R"), local = FALSE)

# Shared export and label helpers for official benchmark adapters.
llr_source_config()

EnsureExport <- function(ds_name, need = c("expr.mtx", "meta.csv", "cells.txt")) {
  rd <- file.path(CONFIG$dir$r_export, ds_name)
  export_status <- sys::exec_wait(
    CONFIG$exe$rscript,
    c("--encoding=UTF-8", file.path(llr_pipeline_dir(), "adapters/00_export_R.R"), ds_name)
  )
  if (export_status != 0L) {
    stop("Export validation failed for ", ds_name)
  }
  missing <- need[!file.exists(file.path(rd, need))]
  if (length(missing)) {
    stop("Export incomplete for ", ds_name, ": ", paste(missing, collapse = ", "))
  }
  rd
}

BenchmarkTrack <- function() {
  track <- Sys.getenv("LLR_BENCHMARK_TRACK", unset = "annotation_free")
  match.arg(track, c("annotation_free", "common_root", "legacy_annotated"))
}

PickClusterLabels <- function(meta, cells) {
  meta <- as.data.frame(meta)
  # 00_export_R: cells.txt = obs_names; meta.csv rows follow the same h5ad order.
  if (nrow(meta) == length(cells)) {
    idx <- match(cells, meta$cell_id)
    if (!any(is.na(idx))) {
      meta <- meta[idx, , drop = FALSE]
    } else {
      idx <- match(cells, as.character(meta[[1L]]))
      if (!any(is.na(idx))) {
        meta <- meta[idx, , drop = FALSE]
      } else {
        meta <- meta[seq_along(cells), , drop = FALSE]
      }
    }
  } else if ("cell_id" %in% names(meta)) {
    idx <- match(cells, meta$cell_id)
    if (any(is.na(idx))) {
      stop("PickClusterLabels: ", sum(is.na(idx)), " cells missing from meta$cell_id")
    }
    meta <- meta[idx, , drop = FALSE]
  } else {
    stop("PickClusterLabels: cannot align meta (", nrow(meta), " rows) to cells (", length(cells), ")")
  }
  cluster_columns <- if (BenchmarkTrack() == "legacy_annotated") {
    c("clusters", "cluster", "cell_type", "seurat_clusters")
  } else {
    "benchmark_clusters"
  }
  for (cn in cluster_columns) {
    if (!cn %in% names(meta)) next
    v <- as.character(meta[[cn]])
    v[is.na(v)] <- ""
    if (length(unique(v[v != ""])) >= 2L) return(v)
  }
  NULL
}

ComputeFallbackClusters <- function(expr, cells, n_clusters = 5L) {
  # Derive cluster labels from PCA + kmeans when no annotation columns exist.
  # expr: genes x cells sparse/dense matrix; cells: character vector of cell names.
  # Returns character vector of cluster labels (same length as cells).
  if (ncol(expr) < 10L) n_clusters <- max(2L, ncol(expr) %/% 3L)
  # Filter zero-variance genes to avoid prcomp scale error and reduce memory
  gene_var <- if (inherits(expr, "sparseMatrix")) {
    # Compute row variance on sparse matrix without dense conversion
    Matrix::rowMeans(expr^2) - Matrix::rowMeans(expr)^2
  } else {
    apply(expr, 1, var)
  }
  keep <- which(is.finite(gene_var) & gene_var > 1e-10)
  if (length(keep) < 2L) {
    cl <- rep("0", length(cells))
    names(cl) <- cells
    return(cl)
  }
  # Limit to top variable genes (cap memory)
  top_n <- min(3000L, length(keep))
  if (length(keep) > top_n) {
    keep <- keep[order(gene_var[keep], decreasing = TRUE)[seq_len(top_n)]]
  }
  expr_filt <- expr[keep, , drop = FALSE]
  expr_t <- t(as.matrix(expr_filt))
  n_pcs <- min(30L, ncol(expr_t) - 1L, nrow(expr_t) - 1L)
  if (n_pcs < 2L) {
    cl <- rep("0", length(cells))
    names(cl) <- cells
    return(cl)
  }
  pc <- stats::prcomp(expr_t, rank. = n_pcs, scale. = TRUE, center = TRUE)$x
  km <- stats::kmeans(pc, centers = min(n_clusters, nrow(pc)), nstart = 10, iter.max = 50)
  cl <- as.character(km$cluster)
  names(cl) <- cells
  cl
}

PickMonocle3RootCells <- function(cds, meta) {
  # Official monocle3 requires user-specified root_cells; benchmark protocol:
  # use cluster label '0' / 'root' / minimum numeric cluster if present in meta.
  if (BenchmarkTrack() == "common_root") {
    root_cell <- Sys.getenv("LLR_BENCHMARK_ROOT_CELL")
    if (!root_cell %in% colnames(cds)) {
      stop("monocle3: common-root cell is absent from input: ", root_cell)
    }
    return(root_cell)
  }
  ct_col <- NULL
  cluster_columns <- if (BenchmarkTrack() == "legacy_annotated") {
    c("clusters", "cluster", "cell_type", "celltype", "seurat_clusters")
  } else {
    "benchmark_clusters"
  }
  for (cn in cluster_columns) {
    if (!cn %in% names(meta)) next
    v <- as.character(meta[[cn]])
    valid <- !is.na(v) & nzchar(v)
    if (!any(valid)) next
    ct_col <- cn
    break
  }
  if (is.null(ct_col)) {
    # No annotation in meta — use monocle3's internal cluster_cells() result.
    # monocle3 stores clusters in colData as "<reduction>_cluster" (e.g. "UMAP_cluster")
    # and in cds@clusters$UMAP$clusters.
    cl_cds <- NULL
    cd <- colData(cds)
    for (cluster_col in c("UMAP_cluster", "cluster", "Lattice_cluster")) {
      if (cluster_col %in% colnames(cd)) {
        cl_cds <- as.character(cd[[cluster_col]])
        if (length(cl_cds) > 0L && any(!is.na(cl_cds))) break
        cl_cds <- NULL
      }
    }
    if (is.null(cl_cds) && !is.null(cds@clusters) && "UMAP" %in% names(cds@clusters)) {
      obj <- cds@clusters$UMAP$clusters
      if (!is.null(obj) && length(obj) > 0L) cl_cds <- as.character(obj)
    }
    if (is.null(cl_cds) || length(cl_cds) == 0L) {
      stop("monocle3: need cluster column in meta or cds to pick root_cells (official order_cells requirement)")
    }
    cl <- cl_cds
  } else {
    cl <- as.character(meta[[ct_col]])
  }
  root_cl <- NULL
  for (candid in c("0", "root", "Root", "STEM", "stem", "branch_0", "early", "Early")) {
    if (candid %in% cl) { root_cl <- candid; break }
  }
  if (is.null(root_cl)) {
    num <- suppressWarnings(as.numeric(cl))
    if (any(is.finite(num))) root_cl <- as.character(min(num, na.rm = TRUE))
  }
  if (is.null(root_cl) && "timepoint" %in% names(meta)) {
    tp <- as.character(meta$timepoint)
    valid <- !is.na(tp) & nzchar(tp)
    if (any(valid)) {
      tp_valid <- tp[valid]
      root_tp <- sort(unique(tp_valid))[1]
      cells <- colnames(cds)[tp == root_tp & valid]
      if (length(cells) >= 1L) return(cells)
    }
  }
  # fallback: if still unknown, use the most frequent cluster label
  if (is.null(root_cl)) {
    tab <- sort(table(cl), decreasing = TRUE)
    if (length(tab) >= 1L) root_cl <- names(tab)[1]
  }
  if (is.null(root_cl)) {
    stop("monocle3: cannot determine root cluster from meta$", ct_col)
  }
  cells <- colnames(cds)[cl == root_cl]
  if (length(cells) < 1L) stop("monocle3: root cluster '", root_cl, "' has zero cells")
  cells
}
