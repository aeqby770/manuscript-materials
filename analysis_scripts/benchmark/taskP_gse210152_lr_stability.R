# GSE210152 temporal-LR stability under atlas-coverage and patient perturbations.

# path bootstrap (relative; set LLR064_WORK_ROOT to override)
.llr_cmd <- commandArgs(trailingOnly = FALSE)
.llr_file <- sub("^--file=", "", grep("^--file=", .llr_cmd, value = TRUE)[1])
if (!length(.llr_file) || is.na(.llr_file) || !nzchar(.llr_file)) {
  stop("Run via Rscript so --file= is set")
}
.llr_here <- dirname(normalizePath(.llr_file, winslash = "/", mustWork = TRUE))
source(file.path(.llr_here, "../config/llr_bootstrap.R"), local = FALSE)
.llr_assets <- llr_assets_dir()
.llr_config_env <- new.env(parent = globalenv())
sys.source(file.path(.llr_assets, "pipeline", "config", "00_config.R"), envir = .llr_config_env)
PROJECT_CONFIG <- .llr_config_env$CONFIG
rm(.llr_config_env)
.llr_package_library <- Sys.getenv(
  "LLR064_RLIB",
  unset = file.path(.llr_assets, "rlib")
)
.libPaths(c(.llr_package_library, .libPaths()))

set.seed(12345)
Sys.setenv(OMP_NUM_THREADS = parallel::detectCores())
data.table::setDTthreads(parallel::detectCores())

CONFIG <- list(
  paths = list(
    rds = file.path(.llr_assets, "cache/linealinkr_workspace/GSE210152_athero/linealinkr_complete.rds"),
    matrix = file.path(.llr_assets, "cache/linealinkr_workspace/GSE210152_athero/extract/GSE210152_athero_matrix.mtx"),
    genes = file.path(.llr_assets, "cache/linealinkr_workspace/GSE210152_athero/extract/GSE210152_athero_genes.txt"),
    cells = file.path(.llr_assets, "cache/linealinkr_workspace/GSE210152_athero/extract/GSE210152_athero_cells.txt"),
    meta = file.path(.llr_assets, "cache/linealinkr_workspace/GSE210152_athero/extract/GSE210152_athero_meta.csv"),
    mapping = file.path(.llr_assets, "results/biology/GSE210152/gse210152_atlas_cell_annotations.csv.gz"),
    out = file.path(.llr_assets, "results/biology/GSE210152")
  ),
  threads = parallel::detectCores(),
  seed = 12345,
  package_version = PROJECT_CONFIG$llr$package_version,
  benchmark_protocol = PROJECT_CONFIG$llr$benchmark_protocol,
  result_revision = PROJECT_CONFIG$llr$result_revision,
  n_bins = 12L,
  top_pool = 100L,
  top_set = 30L,
  lr_scale_factor = 10000
)

dir.create(CONFIG$paths$out, recursive = TRUE, showWarnings = FALSE)

result <- readRDS(CONFIG$paths$rds)
stopifnot(
  identical(as.character(packageVersion("LineaLinkR")), CONFIG$package_version),
  identical(as.character(result$package_version), CONFIG$package_version),
  identical(
    as.character(result$benchmark_protocol),
    CONFIG$benchmark_protocol
  ),
  identical(
    as.character(result$result_revision),
    CONFIG$result_revision
  )
)
expr <- as(Matrix::readMM(CONFIG$paths$matrix), "CsparseMatrix")
genes <- readLines(CONFIG$paths$genes, encoding = "UTF-8")
cells <- readLines(CONFIG$paths$cells, encoding = "UTF-8")
meta <- data.table::fread(CONFIG$paths$meta, encoding = "UTF-8")
mapping <- data.table::fread(CONFIG$paths$mapping, encoding = "UTF-8")

rownames(expr) <- genes
colnames(expr) <- cells
meta <- meta[match(cells, barcode_id)]
mapping <- mapping[match(cells, barcode_id)]
stopifnot(identical(cells, result$cells))
stopifnot(all(meta$barcode_id == cells), all(mapping$barcode_id == cells))

celltypes <- as.integer(as.factor(as.character(mapping$clusters)))
expr_lr <- expr %*% Matrix::Diagonal(x = CONFIG$lr_scale_factor / Matrix::colSums(expr))
expr_lr@x <- log1p(expr_lr@x)

pair_map <- data.table::as.data.table(result$pair_gene_map)
ligand_rows <- lapply(strsplit(pair_map$ligand_genes, "_", fixed = TRUE), match, table = genes)
receptor_rows <- lapply(strsplit(pair_map$receptor_genes, "_", fixed = TRUE), match, table = genes)
lr_db <- list(ligand_rows = ligand_rows, receptor_rows = receptor_rows)

key_cols <- c("pair_idx", "source", "target", "lineage")
full_scores <- data.table::as.data.table(result$temporal_lr$scores)
full_rank <- full_scores[
  is.finite(score) & structural_valid & expression_supported,
  .(
    temporal_effect = max(score) - min(score),
    mean_temporal_score = mean(score),
    n_valid_bins = .N
  ),
  by = key_cols
][n_valid_bins >= 2L][order(-temporal_effect, -mean_temporal_score)]
full_rank <- pair_map[full_rank, on = "pair_idx"]
full_rank[, pair_name := paste(ligand_genes, receptor_genes, sep = "–")]
full_rank[, trajectory_key := paste(pair_idx, source, target, lineage, sep = "|")]
data.table::fwrite(full_rank, file.path(CONFIG$paths$out, "gse210152_full_temporal_lr_rank.csv"))

keep_no_insufficient <- !mapping$clusters %in% c(2L, 3L)
temporal_no_insufficient <- LineaLinkR::ComputeTemporalLRWithIntermediate(
  expr_lr[, keep_no_insufficient, drop = FALSE],
  celltypes[keep_no_insufficient],
  lr_db,
  result$pseudotime[keep_no_insufficient],
  result$lineage_weights[keep_no_insufficient, , drop = FALSE],
  n_bins = CONFIG$n_bins
)
rank_no_insufficient <- data.table::as.data.table(temporal_no_insufficient$scores)[
  is.finite(score) & structural_valid & expression_supported,
  .(temporal_effect = max(score) - min(score), mean_temporal_score = mean(score), n_valid_bins = .N),
  by = key_cols
][n_valid_bins >= 2L][order(-temporal_effect, -mean_temporal_score)]
rank_no_insufficient <- pair_map[rank_no_insufficient, on = "pair_idx"]
rank_no_insufficient[, pair_name := paste(ligand_genes, receptor_genes, sep = "–")]
rank_no_insufficient[, trajectory_key := paste(pair_idx, source, target, lineage, sep = "|")]
data.table::fwrite(rank_no_insufficient, file.path(CONFIG$paths$out, "gse210152_exclude_leiden2_3_temporal_lr_rank.csv"))

keep_no_patient_cluster <- mapping$clusters != 3L
temporal_no_patient_cluster <- LineaLinkR::ComputeTemporalLRWithIntermediate(
  expr_lr[, keep_no_patient_cluster, drop = FALSE],
  celltypes[keep_no_patient_cluster],
  lr_db,
  result$pseudotime[keep_no_patient_cluster],
  result$lineage_weights[keep_no_patient_cluster, , drop = FALSE],
  n_bins = CONFIG$n_bins
)
rank_no_patient_cluster <- data.table::as.data.table(temporal_no_patient_cluster$scores)[
  is.finite(score) & structural_valid & expression_supported,
  .(temporal_effect = max(score) - min(score), mean_temporal_score = mean(score), n_valid_bins = .N),
  by = key_cols
][n_valid_bins >= 2L][order(-temporal_effect, -mean_temporal_score)]
rank_no_patient_cluster <- pair_map[rank_no_patient_cluster, on = "pair_idx"]
rank_no_patient_cluster[, pair_name := paste(ligand_genes, receptor_genes, sep = "–")]
rank_no_patient_cluster[, trajectory_key := paste(pair_idx, source, target, lineage, sep = "|")]
data.table::fwrite(rank_no_patient_cluster, file.path(CONFIG$paths$out, "gse210152_exclude_leiden3_temporal_lr_rank.csv"))

keep_exact_mapped <- mapping$atlas_mapped
temporal_exact_mapped <- LineaLinkR::ComputeTemporalLRWithIntermediate(
  expr_lr[, keep_exact_mapped, drop = FALSE],
  celltypes[keep_exact_mapped],
  lr_db,
  result$pseudotime[keep_exact_mapped],
  result$lineage_weights[keep_exact_mapped, , drop = FALSE],
  n_bins = CONFIG$n_bins
)
rank_exact_mapped <- data.table::as.data.table(temporal_exact_mapped$scores)[
  is.finite(score) & structural_valid & expression_supported,
  .(temporal_effect = max(score) - min(score), mean_temporal_score = mean(score), n_valid_bins = .N),
  by = key_cols
][n_valid_bins >= 2L][order(-temporal_effect, -mean_temporal_score)]
rank_exact_mapped <- pair_map[rank_exact_mapped, on = "pair_idx"]
rank_exact_mapped[, pair_name := paste(ligand_genes, receptor_genes, sep = "–")]
rank_exact_mapped[, trajectory_key := paste(pair_idx, source, target, lineage, sep = "|")]
data.table::fwrite(rank_exact_mapped, file.path(CONFIG$paths$out, "gse210152_exact_mapped_temporal_lr_rank.csv"))

common_no_insufficient <- merge(
  full_rank[, .(trajectory_key, full_effect = temporal_effect)],
  rank_no_insufficient[, .(trajectory_key, subset_effect = temporal_effect)],
  by = "trajectory_key"
)
common_no_patient_cluster <- merge(
  full_rank[, .(trajectory_key, full_effect = temporal_effect)],
  rank_no_patient_cluster[, .(trajectory_key, subset_effect = temporal_effect)],
  by = "trajectory_key"
)
common_exact_mapped <- merge(
  full_rank[, .(trajectory_key, full_effect = temporal_effect)],
  rank_exact_mapped[, .(trajectory_key, subset_effect = temporal_effect)],
  by = "trajectory_key"
)

top_full <- head(full_rank$trajectory_key, CONFIG$top_set)
top_no_insufficient <- head(rank_no_insufficient$trajectory_key, CONFIG$top_set)
top_no_patient_cluster <- head(rank_no_patient_cluster$trajectory_key, CONFIG$top_set)
top_exact_mapped <- head(rank_exact_mapped$trajectory_key, CONFIG$top_set)

subset_summary <- data.table::data.table(
  analysis = c("full", "exclude_Leiden2_3", "exclude_Leiden3", "exact_atlas_mapped_only"),
  n_cells = c(length(cells), sum(keep_no_insufficient), sum(keep_no_patient_cluster), sum(keep_exact_mapped)),
  n_temporal_trajectories = c(nrow(full_rank), nrow(rank_no_insufficient), nrow(rank_no_patient_cluster), nrow(rank_exact_mapped)),
  n_common_with_full = c(nrow(full_rank), nrow(common_no_insufficient), nrow(common_no_patient_cluster), nrow(common_exact_mapped)),
  effect_spearman_vs_full = c(
    1,
    stats::cor(common_no_insufficient$full_effect, common_no_insufficient$subset_effect, method = "spearman"),
    stats::cor(common_no_patient_cluster$full_effect, common_no_patient_cluster$subset_effect, method = "spearman"),
    stats::cor(common_exact_mapped$full_effect, common_exact_mapped$subset_effect, method = "spearman")
  ),
  top30_jaccard_vs_full = c(
    1,
    length(intersect(top_full, top_no_insufficient)) / length(union(top_full, top_no_insufficient)),
    length(intersect(top_full, top_no_patient_cluster)) / length(union(top_full, top_no_patient_cluster)),
    length(intersect(top_full, top_exact_mapped)) / length(union(top_full, top_exact_mapped))
  )
)
data.table::fwrite(subset_summary, file.path(CONFIG$paths$out, "gse210152_subset_temporal_lr_stability.csv"))
data.table::fwrite(
  subset_summary,
  file.path(CONFIG$paths$out, "gse210152_pooled_stability_064.csv")
)

pool_keys <- full_rank[seq_len(min(CONFIG$top_pool, nrow(full_rank))), ..key_cols]
pool_scores <- merge(full_scores, pool_keys, by = key_cols)
patient_effects <- LineaLinkR::ComputeTemporalEffectsBySample(
  scores_df = as.data.frame(pool_scores),
  block_labels = meta$patient,
  pseudotime = result$pseudotime,
  lineage_weights = result$lineage_weights,
  celltypes = celltypes,
  expr = expr_lr,
  lr_db = lr_db,
  temporal_null_per_block = FALSE,
  n_null_median_draws = 199L,
  n_bins = CONFIG$n_bins,
  min_block_cells = 10L
)
patient_meta <- LineaLinkR::MetaAnalyzeTemporalEffects(
  patient_effects,
  effect_col = "temporal_effect",
  se_col = "effect_se",
  method = "fixed",
  include_lodo = FALSE
)
patient_meta_random <- LineaLinkR::MetaAnalyzeTemporalEffects(
  patient_effects,
  effect_col = "temporal_effect",
  se_col = "effect_se",
  method = "random",
  include_lodo = FALSE
)
patient_lopo <- LineaLinkR::LeaveOneSampleOut(
  patient_effects,
  effect_col = "temporal_effect",
  se_col = "effect_se"
)
patient_effects <- pair_map[data.table::as.data.table(patient_effects), on = "pair_idx"]
patient_meta <- pair_map[data.table::as.data.table(patient_meta), on = "pair_idx"]
patient_meta_random <- pair_map[data.table::as.data.table(patient_meta_random), on = "pair_idx"]
patient_lopo <- pair_map[data.table::as.data.table(patient_lopo), on = "pair_idx"]
data.table::fwrite(patient_effects, file.path(CONFIG$paths$out, "gse210152_patient_temporal_effects.csv"))
data.table::fwrite(patient_meta, file.path(CONFIG$paths$out, "gse210152_patient_meta_effects.csv"))
data.table::fwrite(patient_meta_random, file.path(CONFIG$paths$out, "gse210152_patient_meta_random_effects.csv"))
data.table::fwrite(patient_lopo, file.path(CONFIG$paths$out, "gse210152_patient_lopo_effects.csv"))

print(subset_summary)
print(dim(patient_effects))
print(dim(patient_meta))
print(dim(patient_meta_random))
print(dim(patient_lopo))
