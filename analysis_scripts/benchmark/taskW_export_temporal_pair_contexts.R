.llr_cmd <- commandArgs(trailingOnly = FALSE)
.llr_file <- sub("^--file=", "", grep("^--file=", .llr_cmd, value = TRUE)[1L])
stopifnot(length(.llr_file) == 1L, !is.na(.llr_file), nzchar(.llr_file))
.llr_here <- dirname(normalizePath(.llr_file, winslash = "/", mustWork = TRUE))
source(file.path(.llr_here, "../config/llr_bootstrap.R"), local = TRUE, encoding = "UTF-8")
.llr_assets <- llr_assets_dir()
.llr_config_env <- new.env(parent = globalenv())
sys.source(file.path(.llr_assets, "pipeline", "config", "00_config.R"), envir = .llr_config_env)
PROJECT_CONFIG <- .llr_config_env$CONFIG
rm(.llr_config_env)

suppressPackageStartupMessages(library(data.table))
set.seed(12345)
Sys.setenv(OMP_NUM_THREADS = parallel::detectCores())
data.table::setDTthreads(parallel::detectCores())

CONFIG <- list(
  paths = list(
    workspace = PROJECT_CONFIG$dir$linealinkr,
    external = file.path(
      PROJECT_CONFIG$dir$competitor_baseline,
      "paper_tables", "supplement", "tables", "external_cci_comparison"
    ),
    output = file.path(.llr_assets, "results", "biology", "temporal_context"),
    pancreas_temporal_only = file.path(
      .llr_assets, "results", "biology", "pancreas", "pancreas_temporal_only_064.csv"
    ),
    dentategyrus_temporal_only = file.path(
      .llr_assets, "results", "biology", "dentategyrus",
      "dentategyrus_temporal_only_064.csv"
    )
  ),
  package_version = PROJECT_CONFIG$llr$package_version,
  benchmark_protocol = PROJECT_CONFIG$llr$benchmark_protocol,
  result_revision = PROJECT_CONFIG$llr$result_revision,
  temporal_quantile = 0.90,
  pancreas_n = 4L,
  dentategyrus_n = 2L,
  threads = parallel::detectCores(),
  seed = 12345L
)

dir.create(CONFIG$paths$output, recursive = TRUE, showWarnings = FALSE)

pancreas <- readRDS(file.path(
  CONFIG$paths$workspace, "pancreas", "linealinkr_complete.rds"
))
dentategyrus <- readRDS(file.path(
  CONFIG$paths$workspace, "dentategyrus", "linealinkr_complete.rds"
))
stopifnot(
  identical(pancreas$package_version, CONFIG$package_version),
  identical(dentategyrus$package_version, CONFIG$package_version),
  identical(pancreas$benchmark_protocol, CONFIG$benchmark_protocol),
  identical(dentategyrus$benchmark_protocol, CONFIG$benchmark_protocol),
  identical(pancreas$result_revision, CONFIG$result_revision),
  identical(dentategyrus$result_revision, CONFIG$result_revision)
)

pan_scores <- merge(
  as.data.table(pancreas$temporal_lr$scores),
  unique(as.data.table(pancreas$pair_gene_map)),
  by = "pair_idx", all.x = TRUE, sort = FALSE
)
dent_scores <- merge(
  as.data.table(dentategyrus$temporal_lr$scores),
  unique(as.data.table(dentategyrus$pair_gene_map)),
  by = "pair_idx", all.x = TRUE, sort = FALSE
)
pan_scores[, pair_key := tolower(paste(
  ligand_genes, receptor_genes, sep = "___"
))]
dent_scores[, pair_key := tolower(paste(
  ligand_genes, receptor_genes, sep = "___"
))]
pan_programs <- pan_scores[
  structural_valid %in% TRUE & expression_supported %in% TRUE &
    is.finite(score),
  .(
    temporal_mean_score = mean(score),
    temporal_max_score = max(score),
    n_bins_detected = uniqueN(bin),
    n_lineages_detected = uniqueN(lineage),
    n_source_target = uniqueN(paste(source, target, sep = "--"))
  ),
  by = .(ligand_genes, receptor_genes)
]
dent_programs <- dent_scores[
  structural_valid %in% TRUE & expression_supported %in% TRUE &
    is.finite(score),
  .(
    temporal_mean_score = mean(score),
    temporal_max_score = max(score),
    n_bins_detected = uniqueN(bin),
    n_lineages_detected = uniqueN(lineage),
    n_source_target = uniqueN(paste(source, target, sep = "--"))
  ),
  by = .(ligand_genes, receptor_genes)
]
pan_programs[, `:=`(
  pair_name = paste(ligand_genes, receptor_genes, sep = "_"),
  pair_key = tolower(paste(ligand_genes, receptor_genes, sep = "___"))
)]
dent_programs[, `:=`(
  pair_name = paste(ligand_genes, receptor_genes, sep = "_"),
  pair_key = tolower(paste(ligand_genes, receptor_genes, sep = "___"))
)]
setorder(pan_programs, -temporal_mean_score, pair_key)
setorder(dent_programs, -temporal_mean_score, pair_key)
pan_programs[, temporal_pair_rank := seq_len(.N)]
dent_programs[, temporal_pair_rank := seq_len(.N)]
pan_programs[, temporal_top10pct := temporal_pair_rank <= ceiling(.N * 0.10)]
dent_programs[, temporal_top10pct := temporal_pair_rank <= ceiling(.N * 0.10)]

pan_external <- fread(file.path(
  CONFIG$paths$external, "pancreas_llr_vs_external_pair_table.csv"
))[, .(
  pair_key,
  cellchat_supported = cellchat_supported %in% TRUE,
  liana_supported = liana_supported %in% TRUE,
  liana_top10pct = liana_top10pct %in% TRUE,
  cellchat_mean_prob,
  liana_pair_rank
)]
dent_external <- fread(file.path(
  CONFIG$paths$external, "dentategyrus_llr_vs_external_pair_table.csv"
))[, .(
  pair_key,
  cellchat_supported = cellchat_supported %in% TRUE,
  liana_supported = liana_supported %in% TRUE,
  liana_top10pct = liana_top10pct %in% TRUE,
  cellchat_mean_prob,
  liana_pair_rank
)]
pan_programs <- merge(
  pan_programs, pan_external, by = "pair_key", all.x = TRUE, sort = FALSE
)
dent_programs <- merge(
  dent_programs, dent_external, by = "pair_key", all.x = TRUE, sort = FALSE
)
pan_programs[, `:=`(
  cellchat_supported = fifelse(is.na(cellchat_supported), FALSE, cellchat_supported),
  liana_supported = fifelse(is.na(liana_supported), FALSE, liana_supported),
  liana_top10pct = fifelse(is.na(liana_top10pct), FALSE, liana_top10pct)
)]
dent_programs[, `:=`(
  cellchat_supported = fifelse(is.na(cellchat_supported), FALSE, cellchat_supported),
  liana_supported = fifelse(is.na(liana_supported), FALSE, liana_supported),
  liana_top10pct = fifelse(is.na(liana_top10pct), FALSE, liana_top10pct)
)]
pan_programs[, temporal_only_vs_external :=
  temporal_top10pct & !cellchat_supported & !liana_top10pct]
dent_programs[, temporal_only_vs_external :=
  temporal_top10pct & !cellchat_supported & !liana_top10pct]
setorder(pan_programs, -temporal_mean_score, pair_key)
setorder(dent_programs, -temporal_mean_score, pair_key)

pan_selected <- pan_programs[
  temporal_only_vs_external %in% TRUE
][seq_len(min(CONFIG$pancreas_n, .N))]
dent_selected <- dent_programs[
  temporal_only_vs_external %in% TRUE
][seq_len(min(CONFIG$dentategyrus_n, .N))]
stopifnot(nrow(pan_selected) == 4L, nrow(dent_selected) == 2L)

pan_context <- merge(
  pan_scores,
  pan_selected[, .(pair_key, pair_name)],
  by = "pair_key", all = FALSE, sort = FALSE
)[
  structural_valid %in% TRUE & expression_supported %in% TRUE &
    is.finite(score)
][
  order(-score, bin, lineage),
  .SD[1L],
  by = pair_key
]
dent_context <- merge(
  dent_scores,
  dent_selected[, .(pair_key, pair_name)],
  by = "pair_key", all = FALSE, sort = FALSE
)[
  structural_valid %in% TRUE & expression_supported %in% TRUE &
    is.finite(score)
][
  order(-score, bin, lineage),
  .SD[1L],
  by = pair_key
]

pan_levels <- c(
  "Alpha", "Beta", "Delta", "Ductal", "Epsilon",
  "Ngn3 high EP", "Ngn3 low EP", "Pre-endocrine"
)
dent_levels <- c(
  "Astrocytes", "Cajal Retzius", "Cck-Tox", "Endothelial", "GABA",
  "Granule immature", "Granule mature", "Microglia", "Mossy", "Neuroblast",
  "OL", "OPC", "Radial Glia-like", "nIPC"
)
pan_context[, `:=`(
  dataset = "pancreas",
  source_label = pan_levels[source],
  target_label = pan_levels[target]
)]
dent_context[, `:=`(
  dataset = "dentategyrus",
  source_label = dent_levels[source],
  target_label = dent_levels[target]
)]
context_columns <- c(
  "dataset", "pair_name", "pair_idx", "ligand_genes", "receptor_genes",
  "source", "target", "source_label", "target_label", "lineage", "bin",
  "score"
)
contexts <- rbindlist(list(
  pan_context[, ..context_columns], dent_context[, ..context_columns]
), fill = TRUE)
setnames(contexts, c("bin", "score"), c("peak_score_bin", "peak_score"))
setorder(contexts, dataset, -peak_score, pair_name)

pan_weights <- as.data.table(pancreas$lineage_weights)
setnames(pan_weights, paste0("L", seq_len(ncol(pan_weights))))
pan_weights[, `:=`(
  cell_id = as.character(pancreas$cells),
  pseudotime = as.numeric(pancreas$pseudotime)
)]
setcolorder(
  pan_weights,
  c("cell_id", "pseudotime", paste0("L", seq_len(pancreas$n_lineages)))
)
dent_weights <- as.data.table(dentategyrus$lineage_weights)
setnames(dent_weights, paste0("L", seq_len(ncol(dent_weights))))
dent_weights[, `:=`(
  cell_id = as.character(dentategyrus$cells),
  pseudotime = as.numeric(dentategyrus$pseudotime)
)]
setcolorder(
  dent_weights,
  c(
    "cell_id", "pseudotime",
    paste0("L", seq_len(dentategyrus$n_lineages))
  )
)

fwrite(
  pan_programs,
  CONFIG$paths$pancreas_temporal_only
)
fwrite(
  dent_programs,
  CONFIG$paths$dentategyrus_temporal_only
)
fwrite(contexts, file.path(CONFIG$paths$output, "temporal_only_top_lineage_contexts.csv"))
fwrite(
  pan_weights,
  file.path(CONFIG$paths$output, "pancreas_cell_pseudotime_lineage_weights.csv")
)
fwrite(
  dent_weights,
  file.path(
    CONFIG$paths$output,
    "dentategyrus_cell_pseudotime_lineage_weights.csv"
  )
)

print(contexts)
