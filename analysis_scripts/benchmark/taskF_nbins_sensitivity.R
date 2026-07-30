# Task F: n_bins sensitivity on pancreas (top-pair stability + role distribution)
# path bootstrap (relative; set LLR_ASSETS to override)
.llr_cmd <- commandArgs(trailingOnly = FALSE)
.llr_file <- sub("^--file=", "", grep("^--file=", .llr_cmd, value = TRUE)[1])
if (!length(.llr_file) || is.na(.llr_file) || !nzchar(.llr_file)) {
  stop("Run via Rscript so --file= is set")
}
.llr_here <- dirname(normalizePath(.llr_file, winslash = "/", mustWork = TRUE))
source(file.path(.llr_here, "../config/llr_bootstrap.R"), local = FALSE)
.llr_assets <- llr_assets_dir()
llr_source_config()
PROJECT_CONFIG <- CONFIG

set.seed(12345)
Sys.setenv(OMP_NUM_THREADS = as.character(parallel::detectCores()))
data.table::setDTthreads(parallel::detectCores())
.libPaths(c(file.path(.llr_assets, "rlib"), .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(LineaLinkR)
  library(CellChat)
})

CONFIG <- list(
  paths = list(
    workspace = file.path(.llr_assets, "cache/linealinkr_workspace"),
    out = file.path(.llr_assets, "results/paper_tables/supplement/tables/sensitivity_sweep")
  ),
  dataset = "pancreas",
  n_bins_grid = c(8L, 12L, 16L, 20L),
  lr_scale_factor = 10000,
  seed = 12345,
  top_n = 30L
)

dir.create(CONFIG$paths$out, recursive = TRUE, showWarnings = FALSE)

BuildLRDB <- function(expr, species = "human") {
  if (!species %in% c("human", "mouse")) {
    return(list(ligand_rows = list(), receptor_rows = list()))
  }
  cc_db <- if (species == "mouse") CellChat::CellChatDB.mouse else CellChat::CellChatDB.human
  lr_pairs <- cc_db$interaction
  ligand_rows <- list()
  receptor_rows <- list()
  gene_names <- rownames(expr)
  for (i in seq_len(nrow(lr_pairs))) {
    lig <- strsplit(lr_pairs$ligand[i], "_")[[1]]
    rec <- strsplit(lr_pairs$receptor[i], "_")[[1]]
    lig_idx <- match(lig, gene_names)
    rec_idx <- match(rec, gene_names)
    if (any(is.na(lig_idx)) || any(is.na(rec_idx))) next
    ligand_rows[[length(ligand_rows) + 1]] <- lig_idx
    receptor_rows[[length(receptor_rows) + 1]] <- rec_idx
  }
  list(ligand_rows = ligand_rows, receptor_rows = receptor_rows)
}
ds <- CONFIG$dataset
rds <- readRDS(file.path(CONFIG$paths$workspace, ds, "linealinkr_complete.rds"))
if (!identical(as.character(rds$package_version), PROJECT_CONFIG$llr$package_version) ||
    !identical(as.character(rds$benchmark_protocol), PROJECT_CONFIG$llr$benchmark_protocol) ||
    !identical(as.character(rds$result_revision), PROJECT_CONFIG$llr$result_revision)) {
  stop("Pancreas RDS does not match the active 0.6.4 result contract.", call. = FALSE)
}
extract <- file.path(CONFIG$paths$workspace, ds, "extract")
genes <- readLines(file.path(extract, paste0(ds, "_genes.txt")), warn = FALSE)
cells <- readLines(file.path(extract, paste0(ds, "_cells.txt")), warn = FALSE)
expr <- Matrix::readMM(file.path(extract, paste0(ds, "_matrix.mtx")))
if (nrow(expr) != length(genes) && ncol(expr) == length(genes)) expr <- t(expr)
expr <- as(expr, "CsparseMatrix")
rownames(expr) <- genes
colnames(expr) <- cells

meta <- data.table::fread(file.path(extract, paste0(ds, "_meta.csv")), encoding = "UTF-8")
meta <- meta[match(cells, as.character(barcode_id))]
if (anyNA(meta$barcode_id)) stop("Pancreas metadata does not align with the active RDS cells.", call. = FALSE)
celltypes <- as.integer(as.factor(as.character(meta$clusters)))
pseudotime <- as.numeric(rds$pseudotime)
lineage_weights <- rds$lineage_weights

expr_lr <- expr %*% Matrix::Diagonal(x = CONFIG$lr_scale_factor / Matrix::colSums(expr))
expr_lr <- as(expr_lr, "dgCMatrix")
expr_lr@x <- log1p(expr_lr@x)
lr_db <- BuildLRDB(expr_lr, species = "mouse")
pair_map <- as.data.table(rds$pair_gene_map)

pair_tables <- list()
summary_rows <- list()

for (nb in CONFIG$n_bins_grid) {
  message("n_bins=", nb)
  temporal_lr <- LineaLinkR::ComputeTemporalLRWithIntermediate(
    expr_lr, celltypes, lr_db, pseudotime, lineage_weights, n_bins = as.integer(nb)
  )
  sc_all <- as.data.table(temporal_lr$scores)
  sc <- sc_all[
    structural_valid %in% TRUE & expression_supported %in% TRUE &
      is.finite(score)
  ]
  sc <- merge(sc, pair_map, by = "pair_idx", all.x = TRUE)
  agg <- sc[, .(
    temporal_mean_score = mean(score, na.rm = TRUE),
    temporal_max_score = max(score, na.rm = TRUE),
    n_bins_detected = uniqueN(bin),
    n_lineages_detected = uniqueN(lineage)
  ), by = .(ligand_genes, receptor_genes)]
  setorder(agg, -temporal_mean_score)
  agg[, pair_rank := seq_len(.N)]
  agg[, n_bins_param := nb]
  pair_tables[[as.character(nb)]] <- agg

  entropy_raw <- LineaLinkR::ComputeLineageEntropyRaw(temporal_lr$intermediate)
  h_tilde <- LineaLinkR::SmoothMA3(entropy_raw)
  bif_df <- LineaLinkR::FindBifurcationBins(h_tilde)
  n_bif <- sum(bif_df$reason == "crossing", na.rm = TRUE)
  top <- agg[seq_len(min(CONFIG$top_n, .N))]
  summary_rows[[length(summary_rows) + 1L]] <- data.table(
    n_bins = nb,
    n_temporal_rows = nrow(sc_all),
    n_valid_temporal_rows = nrow(sc),
    n_pairs = nrow(agg),
    n_bifurcation_lineages = n_bif,
    top1 = paste(top$ligand_genes[1], top$receptor_genes[1], sep = "_"),
    top1_score = top$temporal_mean_score[1]
  )
}

pair_all <- rbindlist(pair_tables, fill = TRUE)
fwrite(pair_all, file.path(CONFIG$paths$out, "pancreas_nbins_pair_scores.csv"))
sum_dt <- rbindlist(summary_rows)
# Jaccard of top-N pairs vs default 12
def <- pair_tables[["12"]][seq_len(min(CONFIG$top_n, .N)), paste(ligand_genes, receptor_genes, sep = "___")]
jacc <- vapply(names(pair_tables), function(nm) {
  cur <- pair_tables[[nm]][seq_len(min(CONFIG$top_n, .N)), paste(ligand_genes, receptor_genes, sep = "___")]
  inter <- length(intersect(def, cur))
  uni <- length(union(def, cur))
  if (uni == 0) NA_real_ else inter / uni
}, numeric(1))
sum_dt[, jaccard_top30_vs_12 := jacc[as.character(n_bins)]]
fwrite(sum_dt, file.path(CONFIG$paths$out, "pancreas_nbins_sensitivity_summary.csv"))

dentategyrus <- readRDS(file.path(
  CONFIG$paths$workspace, "dentategyrus", "linealinkr_complete.rds"
))
hill_grid <- data.table(hill_k = c(0.5, 1, 2, 5))
hill_grid[, join_key := 1L]
pan_hill <- merge(
  as.data.table(rds$temporal_lr$scores)[
    structural_valid %in% TRUE & expression_supported %in% TRUE &
      is.finite(score),
    .(pair_idx, score, join_key = 1L)
  ],
  hill_grid,
  by = "join_key",
  allow.cartesian = TRUE
)[, join_key := NULL]
dent_hill <- merge(
  as.data.table(dentategyrus$temporal_lr$scores)[
    structural_valid %in% TRUE & expression_supported %in% TRUE &
      is.finite(score),
    .(pair_idx, score, join_key = 1L)
  ],
  hill_grid,
  by = "join_key",
  allow.cartesian = TRUE
)[, join_key := NULL]
pan_hill[, score_k := score / (score + hill_k * (1 - score))]
dent_hill[, score_k := score / (score + hill_k * (1 - score))]
pan_hill_pairs <- merge(
  pan_hill[, .(mean_score = mean(score_k)), by = .(hill_k, pair_idx)],
  unique(as.data.table(rds$pair_gene_map)),
  by = "pair_idx", all.x = TRUE, sort = FALSE
)
dent_hill_pairs <- merge(
  dent_hill[, .(mean_score = mean(score_k)), by = .(hill_k, pair_idx)],
  unique(as.data.table(dentategyrus$pair_gene_map)),
  by = "pair_idx", all.x = TRUE, sort = FALSE
)
setorder(pan_hill_pairs, hill_k, -mean_score, pair_idx)
setorder(dent_hill_pairs, hill_k, -mean_score, pair_idx)
pan_hill_pairs[, pair_rank := seq_len(.N), by = hill_k]
dent_hill_pairs[, pair_rank := seq_len(.N), by = hill_k]
pan_hill_pairs[, pair_key := paste(
  ligand_genes, receptor_genes, sep = "___"
)]
dent_hill_pairs[, pair_key := paste(
  ligand_genes, receptor_genes, sep = "___"
)]
pan_hill_default <- pan_hill_pairs[hill_k == 1 & pair_rank <= 30, pair_key]
dent_hill_default <- dent_hill_pairs[hill_k == 1 & pair_rank <= 30, pair_key]
pan_hill_summary <- pan_hill_pairs[
  pair_rank <= 30,
  .(
    top30_jaccard_vs_k1 = length(intersect(pair_key, pan_hill_default)) /
      length(union(pair_key, pan_hill_default)),
    top_pair = pair_key[which.min(pair_rank)],
    top_pair_score = mean_score[which.min(pair_rank)]
  ),
  by = hill_k
]
dent_hill_summary <- dent_hill_pairs[
  pair_rank <= 30,
  .(
    top30_jaccard_vs_k1 = length(intersect(pair_key, dent_hill_default)) /
      length(union(pair_key, dent_hill_default)),
    top_pair = pair_key[which.min(pair_rank)],
    top_pair_score = mean_score[which.min(pair_rank)]
  ),
  by = hill_k
]
pan_hill_summary[, dataset := "pancreas"]
dent_hill_summary[, dataset := "dentategyrus"]
hill_summary <- rbindlist(list(pan_hill_summary, dent_hill_summary), fill = TRUE)
hill_summary[, transformation := "score/(score+k*(1-score))"]
setcolorder(hill_summary, c(
  "dataset", "hill_k", "top30_jaccard_vs_k1", "top_pair",
  "top_pair_score", "transformation"
))
fwrite(
  rbindlist(list(
    copy(pan_hill_pairs)[, dataset := "pancreas"],
    copy(dent_hill_pairs)[, dataset := "dentategyrus"]
  ), fill = TRUE),
  file.path(CONFIG$paths$out, "developmental_hill_pair_scores_064.csv")
)
fwrite(
  hill_summary,
  file.path(CONFIG$paths$out, "developmental_hill_sensitivity_064.csv")
)

pan_fork_w1 <- LineaLinkR::ComputeForkContrasts(
  rds$temporal_lr$scores, rds$bif_df, rds$sibling_map,
  window_bins = 1L, activity_floor = 1e-4, min_n_eff = 5
)
pan_fork_w2 <- LineaLinkR::ComputeForkContrasts(
  rds$temporal_lr$scores, rds$bif_df, rds$sibling_map,
  window_bins = 2L, activity_floor = 1e-4, min_n_eff = 5
)
pan_fork_w3 <- LineaLinkR::ComputeForkContrasts(
  rds$temporal_lr$scores, rds$bif_df, rds$sibling_map,
  window_bins = 3L, activity_floor = 1e-4, min_n_eff = 5
)
dent_fork_w1 <- LineaLinkR::ComputeForkContrasts(
  dentategyrus$temporal_lr$scores,
  dentategyrus$bif_df,
  dentategyrus$sibling_map,
  window_bins = 1L, activity_floor = 1e-4, min_n_eff = 5
)
dent_fork_w2 <- LineaLinkR::ComputeForkContrasts(
  dentategyrus$temporal_lr$scores,
  dentategyrus$bif_df,
  dentategyrus$sibling_map,
  window_bins = 2L, activity_floor = 1e-4, min_n_eff = 5
)
dent_fork_w3 <- LineaLinkR::ComputeForkContrasts(
  dentategyrus$temporal_lr$scores,
  dentategyrus$bif_df,
  dentategyrus$sibling_map,
  window_bins = 3L, activity_floor = 1e-4, min_n_eff = 5
)
pan_fork <- rbindlist(list(
  copy(as.data.table(pan_fork_w1))[, window_bins := 1L],
  copy(as.data.table(pan_fork_w2))[, window_bins := 2L],
  copy(as.data.table(pan_fork_w3))[, window_bins := 3L]
), fill = TRUE)
dent_fork <- rbindlist(list(
  copy(as.data.table(dent_fork_w1))[, window_bins := 1L],
  copy(as.data.table(dent_fork_w2))[, window_bins := 2L],
  copy(as.data.table(dent_fork_w3))[, window_bins := 3L]
), fill = TRUE)
pan_fork[, trajectory_key := paste(
  pair_idx, source, target, lineage, sep = "--"
)]
dent_fork[, trajectory_key := paste(
  pair_idx, source, target, lineage, sep = "--"
)]
pan_fork_top <- pan_fork[
  evaluable %in% TRUE & is.finite(sibling_contrast)
][order(window_bins, -abs(sibling_contrast))][
  , head(.SD, 30L), by = window_bins
]
dent_fork_top <- dent_fork[
  evaluable %in% TRUE & is.finite(sibling_contrast)
][order(window_bins, -abs(sibling_contrast))][
  , head(.SD, 30L), by = window_bins
]
pan_fork_default <- pan_fork_top[window_bins == 2L, trajectory_key]
dent_fork_default <- dent_fork_top[window_bins == 2L, trajectory_key]
pan_fork_summary <- pan_fork[
  , .(
    n_rows = .N,
    n_evaluable = sum(evaluable %in% TRUE),
    median_abs_sibling_contrast = median(
      abs(sibling_contrast[evaluable %in% TRUE]), na.rm = TRUE
    )
  ),
  by = window_bins
]
dent_fork_summary <- dent_fork[
  , .(
    n_rows = .N,
    n_evaluable = sum(evaluable %in% TRUE),
    median_abs_sibling_contrast = median(
      abs(sibling_contrast[evaluable %in% TRUE]), na.rm = TRUE
    )
  ),
  by = window_bins
]
pan_fork_summary[, top30_jaccard_vs_w2 := vapply(
  window_bins,
  function(w) {
    current <- pan_fork_top[window_bins == w, trajectory_key]
    if (!length(union(current, pan_fork_default))) NA_real_ else
      length(intersect(current, pan_fork_default)) /
      length(union(current, pan_fork_default))
  },
  numeric(1L)
)]
dent_fork_summary[, top30_jaccard_vs_w2 := vapply(
  window_bins,
  function(w) {
    current <- dent_fork_top[window_bins == w, trajectory_key]
    if (!length(union(current, dent_fork_default))) NA_real_ else
      length(intersect(current, dent_fork_default)) /
      length(union(current, dent_fork_default))
  },
  numeric(1L)
)]
pan_fork_summary[, dataset := "pancreas"]
dent_fork_summary[, dataset := "dentategyrus"]
fork_summary <- rbindlist(list(
  pan_fork_summary, dent_fork_summary
), fill = TRUE)
setcolorder(fork_summary, c(
  "dataset", "window_bins", "n_rows", "n_evaluable",
  "median_abs_sibling_contrast", "top30_jaccard_vs_w2"
))
fwrite(
  fork_summary,
  file.path(CONFIG$paths$out, "developmental_fork_window_sensitivity_064.csv")
)

print(sum_dt)
print(hill_summary)
print(fork_summary)
message("done")
