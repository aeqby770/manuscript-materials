
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

ADAPTER_CONFIG <- list(
  pipeline_config = file.path(llr_pipeline_dir(), "config/00_config.R"),
  gsva_min_size = 5L,
  gsva_max_size = 200L,
  gsva_kcdf = "Gaussian",
  n_pcs = 3L,
  db_species = c(human = "HS", mouse = "MM"),
  msigdb_species = c(human = "Homo sapiens", mouse = "Mus musculus"),
  msigdb_collection = c(human = "H", mouse = "MH"),
  implementation_class = "OFFICIAL_TEMPORA_CREATEOBJECT_BUILDTRAJECTORY_GSVA2_COMPAT",
  risk_class = "MEDIUM"
)

library(Matrix)
library(data.table)
library(Tempora)
library(GSVA)
library(msigdbr)

source(ADAPTER_CONFIG$pipeline_config)
llr_source_rel("adapters/_benchmark_io.R")

args <- commandArgs(trailingOnly = TRUE)
ds <- args[[1]]
outdir <- args[[2]]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

rd <- EnsureExport(ds, need = c("expr_lognorm.mtx", "meta.csv", "genes.txt", "cells.txt"))
expr <- as.matrix(readMM(file.path(rd, "expr_lognorm.mtx")))
genes <- fread(file.path(rd, "genes.txt"), header = FALSE)[["V1"]]
cells <- fread(file.path(rd, "cells.txt"), header = FALSE)[["V1"]]
rownames(expr) <- genes
colnames(expr) <- cells

meta <- fread(file.path(rd, "meta.csv"))
setDF(meta)
idx <- match(cells, meta$cell_id)
if (any(is.na(idx))) idx <- match(cells, as.character(meta[[1L]]))
meta <- meta[idx, , drop = FALSE]
rownames(meta) <- cells

cluster_labels <- PickClusterLabels(meta, cells)
if (is.null(cluster_labels)) {
  cluster_labels <- ComputeFallbackClusters(expr, cells)
  cat("Tempora: using fallback PCA+kmeans clusters\n")
}
cluster_factor <- factor(cluster_labels)

# Timepoint resolution: check timepoint first, then alternative time columns
tp <- NULL
tp_col <- NULL
if ("timepoint" %in% names(meta)) {
  v <- as.character(meta$timepoint)
  if (any(!is.na(v) & nzchar(v))) { tp <- v; tp_col <- "timepoint" }
}
if (is.null(tp)) {
  for (alt in c("time_info", "stage", "day", "age(days)", "theiler", "time_point")) {
    if (!alt %in% names(meta)) next
    v <- as.character(meta[[alt]])
    if (any(!is.na(v) & nzchar(v))) { tp <- v; tp_col <- alt; break }
  }
}
if (is.null(tp)) {
  stop("Tempora: no timepoint or alternative time column (time_info/stage/day/age/theiler/time_point) available; dataset is incompatible with Tempora\n")
}
tp_valid <- !is.na(tp) & nzchar(tp)
if (!any(tp_valid)) {
  stop("Tempora: timepoint column '", tp_col, "' has no non-missing values; dataset is incompatible with Tempora\n")
}
time_values <- unique(tp[tp_valid])
time_numeric <- as.numeric(gsub("[^0-9.]+", "", as.character(time_values)))
time_order <- time_values[order(
  fifelse(is.na(time_numeric), seq_along(time_numeric), time_numeric),
  as.character(time_values)
)]

meta_tempora <- data.frame(
  Timepoints = tp,
  Clusters = as.integer(cluster_factor),
  row.names = cells
)
cluster_map <- data.table(
  cluster_id = seq_along(levels(cluster_factor)),
  cluster_label = levels(cluster_factor)
)

tempora_obj <- CreateTemporaObject(expr, meta_tempora, time_order)
cluster_levels <- sort(unique(tempora_obj@meta.data[["Clusters"]]))
cluster_design <- Matrix::sparse.model.matrix(~ factor(tempora_obj@meta.data[["Clusters"]], levels = cluster_levels) + 0)
cluster_design <- cluster_design %*% Matrix::Diagonal(x = 1 / Matrix::colSums(cluster_design))
expr_by_cluster <- expr %*% cluster_design
colnames(expr_by_cluster) <- cluster_levels
rownames(expr_by_cluster) <- rownames(expr)

species_key <- CONFIG$dataset_species[[ds]]
gene_set_table <- msigdbr(
  db_species = ADAPTER_CONFIG$db_species[[species_key]],
  species = ADAPTER_CONFIG$msigdb_species[[species_key]],
  collection = ADAPTER_CONFIG$msigdb_collection[[species_key]]
)
gene_sets <- split(gene_set_table[["gene_symbol"]], gene_set_table[["gs_name"]])
gsva_param <- GSVA::gsvaParam(
  as.matrix(expr_by_cluster),
  gene_sets,
  minSize = ADAPTER_CONFIG$gsva_min_size,
  maxSize = ADAPTER_CONFIG$gsva_max_size,
  kcdf = ADAPTER_CONFIG$gsva_kcdf,
  verbose = FALSE
)

tempora_obj@cluster.pathways <- GSVA::gsva(gsva_param)
tempora_obj@cluster.pathways.dr <- prcomp(t(tempora_obj@cluster.pathways), scale = TRUE, center = TRUE)
tempora_obj <- BuildTrajectory(
  tempora_obj,
  n_pcs = min(ADAPTER_CONFIG$n_pcs, ncol(tempora_obj@cluster.pathways.dr[["rotation"]]))
)

saveRDS(tempora_obj, file.path(outdir, "tempora_result.rds"))
fwrite(as.data.table(tempora_obj@trajectory), file.path(outdir, "tempora_trajectory_edges.csv"))
fwrite(as.data.table(tempora_obj@cluster.metadata), file.path(outdir, "tempora_cluster_metadata.csv"))
fwrite(cluster_map, file.path(outdir, "tempora_cluster_map.csv"))
fwrite(as.data.table(tempora_obj@cluster.pathways, keep.rownames = "pathway"), file.path(outdir, "tempora_pathway_scores.csv"))
fwrite(
  data.table(
    algorithm = "Tempora",
    dataset = ds,
    n_cells = ncol(expr),
    n_genes = nrow(expr),
    n_clusters = uniqueN(meta_tempora$Clusters),
    n_timepoints = uniqueN(meta_tempora$Timepoints),
    n_pathways = nrow(tempora_obj@cluster.pathways),
    n_edges = nrow(tempora_obj@trajectory),
    n_pcs = tempora_obj@n.pcs,
    implementation_class = ADAPTER_CONFIG$implementation_class,
    misrepresentation_risk = ADAPTER_CONFIG$risk_class,
    metric_layer = "trajectory_topology_trace"
  ),
  file.path(outdir, "tempora_summary.csv")
)
