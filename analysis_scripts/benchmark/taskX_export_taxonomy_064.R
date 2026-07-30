#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

.llr_cmd <- commandArgs(trailingOnly = FALSE)
.llr_file <- sub("^--file=", "", grep("^--file=", .llr_cmd, value = TRUE)[1])
.llr_here <- dirname(normalizePath(.llr_file, winslash = "/", mustWork = TRUE))
source(file.path(.llr_here, "../config/llr_bootstrap.R"), local = FALSE)
llr_source_config()
PROJECT_CONFIG <- CONFIG

CONFIG <- list(
  paths = list(
    rlib = file.path(PROJECT_CONFIG$dir$assets, "rlib"),
    rds = file.path(
      PROJECT_CONFIG$dir$linealinkr,
      "pancreas",
      "linealinkr_complete.rds"
    ),
    metadata = file.path(
      PROJECT_CONFIG$dir$linealinkr,
      "pancreas",
      "extract",
      "pancreas_meta.csv"
    ),
    output = file.path(PROJECT_CONFIG$dir$assets, "results", "biology", "taxonomy")
  ),
  threads = parallel::detectCores(),
  seed = 12345L,
  h_fork = 0.5,
  min_n_eff_cells = 5,
  roles = c("instructive", "reinforcing", "feedback", "lateral", "bystander")
)

set.seed(CONFIG$seed)
Sys.setenv(OMP_NUM_THREADS = as.character(CONFIG$threads))
data.table::setDTthreads(CONFIG$threads)
.libPaths(c(CONFIG$paths$rlib, .libPaths()))
dir.create(CONFIG$paths$output, recursive = TRUE, showWarnings = FALSE)

result <- readRDS(CONFIG$paths$rds)
metadata <- fread(CONFIG$paths$metadata, encoding = "UTF-8")
stopifnot(
  identical(as.character(result$package_version), PROJECT_CONFIG$llr$package_version),
  identical(as.character(result$benchmark_protocol), PROJECT_CONFIG$llr$benchmark_protocol),
  identical(as.character(result$result_revision), PROJECT_CONFIG$llr$result_revision),
  identical(as.character(metadata[[1L]]), as.character(result$cells))
)

celltype_labels <- as.character(metadata[[result$celltype_column]])
celltype_labels[is.na(celltype_labels) | celltype_labels == ""] <- "unassigned"
celltypes <- as.integer(as.factor(celltype_labels))

role_features <- LineaLinkR::ComputeRoleFeatures(
  temporal_scores = result$temporal_lr$scores,
  lineage_weights = result$lineage_weights,
  pseudotime = result$pseudotime,
  celltypes = celltypes,
  intermediate = result$temporal_lr$intermediate,
  h_fork = CONFIG$h_fork,
  compute_lag = FALSE,
  timing_source = "cell_state_tau"
)
calibration_mask <-
  role_features$structural_valid %in% TRUE &
  role_features$expression_supported %in% TRUE &
  is.finite(role_features$spec_self) &
  is.finite(role_features$spec_max_other) &
  is.finite(role_features$delta_tau) &
  is.finite(role_features$n_eff_cells) &
  role_features$n_eff_cells >= CONFIG$min_n_eff_cells
scoring_features <- role_features
scoring_features$fully_supported <- calibration_mask
role_inference <- LineaLinkR::ComputeRoleInference(
  role_features = scoring_features,
  h_fork = CONFIG$h_fork,
  min_n_eff_cells = CONFIG$min_n_eff_cells,
  timing_source = "cell_state_tau"
)
taxonomy <- as.data.table(role_inference)
taxonomy[, fully_supported := role_features$fully_supported]
taxonomy[, `:=`(
  instructive = instructive_like,
  reinforcing = reinforcing_like,
  feedback = feedback_like,
  lateral = lateral_like,
  bystander = bystander_like
)]

role_matrix <- as.matrix(taxonomy[, CONFIG$roles, with = FALSE])
finite_role_rows <- rowSums(is.finite(role_matrix)) == length(CONFIG$roles)
taxonomy[, calibration_eligible := calibration_mask]
taxonomy[, dominant_role := NA_character_]
taxonomy[
  calibration_eligible == TRUE,
  dominant_role := CONFIG$roles[max.col(
    role_matrix[calibration_mask, , drop = FALSE],
    ties.method = "first"
  )]
]
taxonomy[, `:=`(
  taxonomy_calibration_status = fifelse(
    calibration_eligible,
    "calibration_eligible",
    "not_calibrated"
  ),
  role_inference_status = fifelse(
    calibration_eligible,
    "exploratory_calibration",
    "not_calibrated"
  ),
  calibration_boundary = paste(
    "coarse pathway-prior alignment;",
    "not trajectory ground-truth accuracy"
  ),
  package_version = as.character(result$package_version),
  benchmark_protocol = as.character(result$benchmark_protocol),
  result_revision = as.character(result$result_revision)
)]

taxonomy <- merge(
  taxonomy,
  as.data.table(result$pair_gene_map),
  by = "pair_idx",
  all.x = TRUE,
  sort = FALSE
)
setorder(taxonomy, pair_idx, source, target, lineage, bin)

stopifnot(
  nrow(taxonomy) == nrow(result$temporal_lr$scores),
  sum(calibration_mask) > 0L,
  identical(calibration_mask, finite_role_rows),
  all(is.finite(taxonomy$spec_self[taxonomy$calibration_eligible])),
  all(taxonomy$spec_self[taxonomy$calibration_eligible] >= 0),
  all(taxonomy$spec_self[taxonomy$calibration_eligible] <= 1),
  all(!is.na(taxonomy$dominant_role[taxonomy$calibration_eligible])),
  all(is.na(taxonomy$dominant_role[!taxonomy$calibration_eligible]))
)

fwrite(taxonomy, file.path(CONFIG$paths$output, "pancreas_taxonomy_full.csv"))
fwrite(taxonomy, file.path(CONFIG$paths$output, "pancreas_taxonomy_rows_064.csv"))
cat(sprintf(
  "TAXONOMY_064_EXPORT_OK rows=%d eligible=%d lineages=%d\n",
  nrow(taxonomy),
  sum(taxonomy$calibration_eligible),
  result$n_lineages
))
