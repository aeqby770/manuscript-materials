# One LineaLinkR pseudotime/lineage core replicate using a freshly computed kNN graph.
# path bootstrap (relative; set LLR_ASSETS to override)
.llr_cmd <- commandArgs(trailingOnly = FALSE)
.llr_file <- sub("^--file=", "", grep("^--file=", .llr_cmd, value = TRUE)[1])
if (!length(.llr_file) || is.na(.llr_file) || !nzchar(.llr_file)) {
  stop("Run via Rscript so --file= is set")
}
.llr_here <- dirname(normalizePath(.llr_file, winslash = "/", mustWork = TRUE))
source(file.path(.llr_here, "../config/llr_bootstrap.R"), local = FALSE)
.llr_assets <- llr_assets_dir()

.libPaths(c(
  Sys.getenv("LLR064_RLIB", unset = "F:/biotrade/0.64LLR/rlib"),
  .libPaths()
))
set.seed(12345)
Sys.setenv(OMP_NUM_THREADS = parallel::detectCores())
data.table::setDTthreads(parallel::detectCores())

library(data.table)
library(Matrix)
library(LineaLinkR)
llr_source_config()
PROJECT_CONFIG <- CONFIG

CONFIG <- list(
  paths = list(
    workspace = file.path(.llr_assets, "cache/linealinkr_workspace")
  ),
  seed = 12345L,
  knn_neighbors = 50L,
  n_bins = 12L,
  min_n_eff = 5.0,
  fork_window_bins = 2L,
  fork_activity_floor = 1e-4,
  lr_scale_factor = 10000,
  return_complete_grid = identical(
    Sys.getenv("LLR_RUNTIME_COMPLETE_GRID", unset = "true"),
    "true"
  )
)

args <- commandArgs(trailingOnly = TRUE)
dataset <- args[[1L]]
replicate_dir <- args[[2L]]
result_path <- file.path(CONFIG$paths$workspace, dataset, "linealinkr_complete.rds")
result <- readRDS(result_path)
extract_dir <- file.path(CONFIG$paths$workspace, dataset, "extract")
genes <- readLines(file.path(extract_dir, sprintf("%s_genes.txt", dataset)), warn = FALSE)
cells <- readLines(file.path(extract_dir, sprintf("%s_cells.txt", dataset)), warn = FALSE)
stopifnot(
  identical(as.character(packageVersion("LineaLinkR")), PROJECT_CONFIG$llr$package_version),
  identical(as.character(result$package_version), PROJECT_CONFIG$llr$package_version),
  identical(as.character(result$benchmark_protocol), PROJECT_CONFIG$llr$benchmark_protocol),
  identical(as.character(result$result_revision), PROJECT_CONFIG$llr$result_revision),
  identical(as.character(result$cells), as.character(cells))
)
meta <- fread(file.path(extract_dir, sprintf("%s_meta.csv", dataset)), encoding = "UTF-8")
meta <- meta[match(cells, meta[[1L]])]
celltypes <- as.integer(as.factor(as.character(meta[[result$celltype_column]])))
expr <- as(
  Matrix::readMM(file.path(extract_dir, sprintf("%s_matrix.mtx", dataset))),
  "CsparseMatrix"
)
rownames(expr) <- genes
colnames(expr) <- cells
n_cells <- length(result$cells)
knn_dist <- as(Matrix::readMM(file.path(replicate_dir, "knn_dist.mtx")), "CsparseMatrix")
triplets <- as.data.table(Matrix::summary(knn_dist))
setorder(triplets, i, x)
triplets[, neighbor_rank := seq_len(.N), by = i]
triplets <- triplets[neighbor_rank <= min(CONFIG$knn_neighbors, n_cells - 1L)]
nn_idx <- as.matrix(dcast(triplets, i ~ neighbor_rank, value.var = "j", fill = 0L))[, -1L, drop = FALSE]
nn_dists <- as.matrix(dcast(triplets, i ~ neighbor_rank, value.var = "x", fill = Inf))[, -1L, drop = FALSE]
pca <- as.matrix(fread(file.path(replicate_dir, "pca.csv"), header = FALSE))
cluster_source <- as.character(result$lineage_cluster_source)
clusters <- if (startsWith(cluster_source, "metadata:")) {
  cluster_column <- sub("^metadata:", "", cluster_source)
  as.integer(as.factor(as.character(meta[[cluster_column]])))
} else {
  n_clusters <- as.integer(sub(".*requested_k=([0-9]+).*", "\\1", cluster_source))
  set.seed(CONFIG$seed)
  kmeans(pca, centers = n_clusters, nstart = 10L, iter.max = 50L)$cluster
}
potency <- log1p(as.numeric(Matrix::colSums(expr > 0)))

t0 <- proc.time()[["elapsed"]]
connected <- LineaLinkR::ConnectKnnComponents(nn_idx, nn_dists, pca)
t_connected <- proc.time()[["elapsed"]]
root_selection <- LineaLinkR::InferStartCellAttractor(
  pca,
  clusters,
  connected$nn_idx,
  connected$nn_dists,
  potency = potency,
  top_fraction = 0.05,
  n_attractors = 0L,
  max_iter = 300L,
  tol = 1e-5,
  weights = c(0.25, 0.20, 0.30, 0.15, 0.10),
  seed = CONFIG$seed
)
start_cell <- as.integer(root_selection$start_cell)
start_cluster <- clusters[[start_cell]]
t_root <- proc.time()[["elapsed"]]
bridged <- LineaLinkR::AugmentKnnWithMstBridges(
  connected$nn_idx,
  connected$nn_dists,
  pca,
  clusters,
  start_cluster
)
pseudotime <- LineaLinkR::ComputePseudotimeFuzzy(
  bridged$nn_idx,
  bridged$nn_dists,
  start_cell = start_cell,
  rank_normalize = FALSE
)
t_pseudotime <- proc.time()[["elapsed"]]
lineage_weights <- LineaLinkR::ComputeLineageWeightsTrunkAware(
  pca,
  clusters,
  start_cluster = start_cluster
)
lineage_paths <- LineaLinkR::ExtractLineages(pca, clusters, start_cluster)
sibling_map <- LineaLinkR::BuildSiblingContrastMap(lineage_paths)
t_lineage <- proc.time()[["elapsed"]]

expr_lr <- expr %*% Matrix::Diagonal(
  x = as.numeric(CONFIG$lr_scale_factor) / Matrix::colSums(expr)
)
expr_lr@x <- log1p(expr_lr@x)
ligand_parts <- strsplit(result$pair_gene_map$ligand_genes, "_", fixed = TRUE)
receptor_parts <- strsplit(result$pair_gene_map$receptor_genes, "_", fixed = TRUE)
lr_db <- list(
  ligand_rows = lapply(ligand_parts, match, table = genes),
  receptor_rows = lapply(receptor_parts, match, table = genes)
)
lr_scores <- LineaLinkR::ComputeLRScores(
  expr_lr,
  celltypes,
  lr_db,
  n_perms = 100L
)
t_static <- proc.time()[["elapsed"]]

temporal_lr <- LineaLinkR::ComputeTemporalLRWithIntermediate(
  expr_lr,
  celltypes,
  lr_db,
  pseudotime,
  lineage_weights,
  n_bins = CONFIG$n_bins,
  min_n_eff = CONFIG$min_n_eff,
  return_complete_grid = CONFIG$return_complete_grid
)
t_temporal_core <- proc.time()[["elapsed"]]
entropy_raw <- LineaLinkR::ComputeLineageEntropyRaw(temporal_lr$intermediate)
h_tilde <- LineaLinkR::SmoothMA3(entropy_raw)
t_entropy <- proc.time()[["elapsed"]]
bif_df <- LineaLinkR::FindBifurcationBins(h_tilde)
t_bifurcation <- proc.time()[["elapsed"]]
fork_contrasts <- LineaLinkR::ComputeForkContrasts(
  temporal_lr$scores,
  bif_df,
  sibling_map,
  window_bins = CONFIG$fork_window_bins,
  activity_floor = CONFIG$fork_activity_floor,
  min_n_eff = CONFIG$min_n_eff
)
t_fork_contrasts <- proc.time()[["elapsed"]]
fork_evidence <- LineaLinkR::ComposeForkEvidence(fork_contrasts)
t_temporal <- proc.time()[["elapsed"]]

runtime_result <- list(
  dataset = dataset,
  n_cells = n_cells,
  start_cell = start_cell,
  start_cluster = start_cluster,
  root_method = root_selection$method,
  lineage_cluster_source = cluster_source,
  pseudotime = pseudotime,
  lineage_weights = lineage_weights,
  lr_scores = lr_scores,
  temporal_lr = temporal_lr,
  bif_df = bif_df,
  fork_contrasts = fork_contrasts,
  fork_evidence = fork_evidence
)
saveRDS(runtime_result, file.path(replicate_dir, "runtime_result.rds"))
t_serialized <- proc.time()[["elapsed"]]

cat(sprintf(
  paste0(
    "%s,%d,%.9f,%.9f,%.9f,%.9f,%.9f,",
    "%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%s,%s,%s,%s\n"
  ),
  dataset,
  n_cells,
  t_connected - t0,
  t_root - t_connected,
  t_pseudotime - t_root,
  t_lineage - t_pseudotime,
  t_static - t_lineage,
  t_temporal_core - t_static,
  t_entropy - t_temporal_core,
  t_bifurcation - t_entropy,
  t_fork_contrasts - t_bifurcation,
  t_temporal - t_fork_contrasts,
  t_temporal - t_static,
  t_serialized - t_temporal,
  t_lineage - t0,
  t_temporal - t0,
  as.character(result$package_version),
  as.character(result$benchmark_protocol),
  as.character(result$result_revision),
  normalizePath(result_path, winslash = "/", mustWork = TRUE)
))
