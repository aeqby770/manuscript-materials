#!/usr/bin/env Rscript

.llr_cmd <- commandArgs(trailingOnly = FALSE)
.llr_file <- sub("^--file=", "", grep("^--file=", .llr_cmd, value = TRUE)[1])
.llr_here <- dirname(normalizePath(.llr_file, winslash = "/", mustWork = TRUE))
source(file.path(.llr_here, "../config/llr_bootstrap.R"), local = FALSE)
.llr_assets <- llr_assets_dir()
.libPaths(c(
  Sys.getenv("LLR064_RLIB", unset = file.path(.llr_assets, "rlib")),
  .libPaths()
))

library(data.table)
llr_source_config()

EXPORT_CONFIG <- list(
  paths = list(
    rds = file.path(
      .llr_assets,
      "cache/linealinkr_workspace/",
      "GSE274389_signal_recorder/linealinkr_complete.rds"
    ),
    cell_output = file.path(
      .llr_assets,
      "results/biology/gse274389/",
      "gse274389_linealinkr_cells.csv"
    ),
    wnt_output = file.path(
      .llr_assets,
      "results/biology/gse274389/",
      "gse274389_wnt_temporal_scores.csv"
    ),
    root_output = file.path(
      .llr_assets,
      "results/biology/gse274389/",
      "gse274389_root_candidates.csv"
    )
  ),
  package_version = CONFIG$llr$package_version,
  benchmark_protocol = CONFIG$llr$benchmark_protocol,
  result_revision = CONFIG$llr$result_revision,
  threads = parallel::detectCores(),
  seed = 12345L
)

set.seed(EXPORT_CONFIG$seed)
Sys.setenv(OMP_NUM_THREADS = EXPORT_CONFIG$threads)
data.table::setDTthreads(EXPORT_CONFIG$threads)

x <- readRDS(EXPORT_CONFIG$paths$rds)
stopifnot(
  identical(x$package_version, EXPORT_CONFIG$package_version),
  identical(x$benchmark_protocol, EXPORT_CONFIG$benchmark_protocol),
  identical(x$result_revision, EXPORT_CONFIG$result_revision),
  grepl("^(blind|anchored)", x$root_evidence_mode),
  grepl("expression_detected_genes_pca_knn_no_gt", x$root_input_source),
  !grepl(
    "ground_truth|true_pseudotime|timepoint|time_hours|day|stage|clone|fate",
    x$root_input_source,
    ignore.case = TRUE
  ),
  length(x$pseudotime_alternative) %in% c(0L, x$n_cells),
  is.null(x$root_alternative) || (
    is.data.frame(x$root_alternative) && nrow(x$root_alternative) == 1L
  ),
  identical(
    length(x$pseudotime_alternative) == x$n_cells,
    is.data.frame(x$root_alternative) && nrow(x$root_alternative) == 1L
  )
)
lineage_weights <- as.matrix(x$lineage_weights)
dominant_lineage <- max.col(lineage_weights, ties.method = "first")
dominant_weight <- lineage_weights[cbind(seq_len(nrow(lineage_weights)), dominant_lineage)]

cell_output <- data.table(
  cell = as.character(x$cells),
  llr_pseudotime = as.numeric(x$pseudotime),
  llr_pseudotime_alternative = c(
    as.numeric(x$pseudotime_alternative),
    rep(NA_real_, x$n_cells)
  )[seq_len(x$n_cells)],
  dominant_lineage = dominant_lineage,
  dominant_lineage_weight = dominant_weight
)

pair_map <- as.data.table(x$pair_gene_map)
temporal <- as.data.table(x$temporal_lr$scores)
wnt_pair_map <- pair_map[grepl("(^|[+_])Wnt[0-9A-Za-z]*($|[+_])", ligand_genes)]
wnt_output <- merge(temporal, wnt_pair_map, by = "pair_idx", all = FALSE)
root_output <- as.data.table(x$root_candidates)
root_output[, `:=`(
  root_role = "other candidate",
  start_cell = NA_integer_
)]
root_output[candidate == x$selected_root_potency_mode, `:=`(
  root_role = "selected",
  start_cell = as.integer(x$start_cell)
)]
alternative_candidate <- c(
  as.character(x$root_alternative$candidate),
  NA_character_
)[1L]
root_output[
  !is.na(alternative_candidate) & candidate == alternative_candidate,
  `:=`(
    root_role = "highest-supported alternative",
    start_cell = as.integer(x$root_alternative$start_cell)
  )
]

dir.create(dirname(EXPORT_CONFIG$paths$cell_output), recursive = TRUE, showWarnings = FALSE)
fwrite(cell_output, EXPORT_CONFIG$paths$cell_output)
fwrite(wnt_output, EXPORT_CONFIG$paths$wnt_output)
fwrite(root_output, EXPORT_CONFIG$paths$root_output)

cat(sprintf("cells=%d lineages=%d wnt_rows=%d wnt_pairs=%d\n",
            nrow(cell_output), ncol(lineage_weights), nrow(wnt_output), uniqueN(wnt_output$pair_idx)))
