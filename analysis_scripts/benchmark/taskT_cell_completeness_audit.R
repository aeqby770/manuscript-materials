script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg), winslash = "/", mustWork = TRUE))
Sys.setenv(LLR_ASSETS = normalizePath(file.path(script_dir, "../.."), winslash = "/", mustWork = TRUE))
source(file.path(script_dir, "../config/llr_bootstrap.R"), local = TRUE)
assets <- llr_assets_dir()
.llr_config_env <- new.env(parent = globalenv())
sys.source(file.path(assets, "pipeline", "config", "00_config.R"), envir = .llr_config_env)
PROJECT_CONFIG <- .llr_config_env$CONFIG
CORE_ALGORITHMS <- .llr_config_env$BENCHMARK_PSEUDOTIME_ALGORITHMS
CORE_DATASETS <- .llr_config_env$BENCHMARK_COMPLETE_CORE_DATASETS
EXTENSION_DATASETS <- .llr_config_env$BENCHMARK_ORDERING_DATASETS
TABLE_PROTOCOL_TAG <- gsub("[^A-Za-z0-9_]+", "_", PROJECT_CONFIG$llr$benchmark_protocol)
rm(.llr_config_env)

library(data.table)

set.seed(12345)
Sys.setenv(OMP_NUM_THREADS = parallel::detectCores())
data.table::setDTthreads(parallel::detectCores())

CONFIG <- list(
  paths = list(
    input = file.path(
      PROJECT_CONFIG$dir$table1_main,
      sprintf("table1_pseudotime_%s.csv", TABLE_PROTOCOL_TAG)
    ),
    output = file.path(
      assets,
      "results/paper_tables/supplement/tables/cell_completeness_audit"
    )
  ),
  threads = parallel::detectCores(),
  seed = 12345L,
  benchmark_protocol = PROJECT_CONFIG$llr$benchmark_protocol,
  result_revision = PROJECT_CONFIG$llr$result_revision,
  plot = list(),
  thresholds = c(0.50, 0.75, 0.90, 1.00)
)

dir.create(CONFIG$paths$output, recursive = TRUE, showWarnings = FALSE)

pt <- fread(CONFIG$paths$input, encoding = "UTF-8")
stopifnot(
  identical(CORE_ALGORITHMS, c("LineaLinkR", PROJECT_CONFIG$algorithm_groups$ti_trajectory)),
  !anyDuplicated(CORE_ALGORITHMS),
  all(CORE_DATASETS %in% EXTENSION_DATASETS),
  all(c(
    "dataset", "algorithm", "Spearman", "Spearman_signed", "orientation_flipped",
    "n_aligned", "n_reference", "cell_coverage", "cell_alignment_mode",
    "direction_policy", "is_duplicate"
  ) %in% names(pt))
)
pt <- pt[
  is_duplicate == FALSE &
    dataset %in% EXTENSION_DATASETS &
    algorithm %in% CORE_ALGORITHMS
]
pt_key_counts <- pt[, .N, by = .(dataset, algorithm)]
stopifnot(
  "Pseudotime result table contains duplicate dataset-algorithm rows" =
    all(pt_key_counts$N == 1L)
)

expected_grid <- CJ(
  dataset = EXTENSION_DATASETS,
  algorithm = CORE_ALGORITHMS,
  unique = TRUE,
  sorted = FALSE
)
per_output <- merge(
  expected_grid,
  pt[, .(
    dataset, algorithm, Spearman, Spearman_signed, orientation_flipped,
    n_aligned, n_reference, cell_coverage, cell_alignment_mode,
    direction_policy
  )],
  by = c("dataset", "algorithm"),
  all.x = TRUE,
  sort = FALSE
)
per_output[, n_reference_cells := n_reference]
per_output[, cell_completeness := n_aligned / n_reference_cells]
per_output[, output_status := fcase(
  !is.finite(Spearman) | !is.finite(cell_completeness) | !is.finite(cell_coverage),
  "missing_or_invalid",
  n_aligned != n_reference_cells | cell_completeness != 1 | cell_coverage != 1,
  "partial_output",
  default = "complete_output"
)]
per_output[, benchmark_scope := fifelse(
  dataset %in% CORE_DATASETS, "complete_core", "scale_extension"
)]
setcolorder(per_output, c(
  "dataset", "algorithm", "benchmark_scope", "output_status", "Spearman",
  "Spearman_signed", "orientation_flipped", "n_aligned", "n_reference_cells",
  "cell_completeness", "cell_coverage", "cell_alignment_mode", "direction_policy"
))
per_output[, algorithm_order := match(algorithm, CORE_ALGORITHMS)]
per_output[, dataset_order := match(dataset, EXTENSION_DATASETS)]
setorder(per_output, algorithm_order, dataset_order)
per_output[, c("algorithm_order", "dataset_order") := NULL]

core_rows <- per_output[dataset %in% CORE_DATASETS]
stopifnot(
  "The fixed complete-core matrix is missing or contains partial outputs" =
    nrow(core_rows) == length(CORE_DATASETS) * length(CORE_ALGORITHMS) &&
    setequal(core_rows$dataset, CORE_DATASETS) &&
    setequal(core_rows$algorithm, CORE_ALGORITHMS) &&
    all(core_rows$output_status == "complete_output")
)

algorithm_summary <- per_output[
  , .(
    n_planned_datasets = uniqueN(dataset),
    n_complete_outputs = sum(output_status == "complete_output"),
    n_partial_outputs = sum(output_status == "partial_output"),
    n_missing_or_invalid_outputs = sum(output_status == "missing_or_invalid"),
    min_cell_completeness = min(cell_completeness, na.rm = TRUE),
    median_cell_completeness = median(cell_completeness, na.rm = TRUE),
    mean_cell_completeness = mean(cell_completeness, na.rm = TRUE)
  ),
  by = algorithm
][order(-n_complete_outputs, algorithm)]

thresholds <- data.table(cell_completeness_threshold = CONFIG$thresholds)
threshold_rows <- cbind(
  per_output[rep(seq_len(.N), each = nrow(thresholds))],
  thresholds[rep(seq_len(.N), times = nrow(per_output))]
)[
  cell_completeness >= cell_completeness_threshold
]
threshold_summary <- threshold_rows[
  , .(
    n_planned_datasets = length(EXTENSION_DATASETS),
    n_qualifying_datasets = uniqueN(dataset)
  ),
  by = .(cell_completeness_threshold, algorithm)
][order(cell_completeness_threshold, -n_qualifying_datasets, algorithm)]

core_rows[, dataset_rank := frank(-Spearman, ties.method = "average"), by = dataset]
full_cell_core_summary <- core_rows[
  , .(
    n_common_datasets = uniqueN(dataset),
    mean_spearman = mean(Spearman),
    median_spearman = median(Spearman),
    mean_rank = mean(dataset_rank),
    min_cell_completeness = min(cell_completeness)
  ),
  by = algorithm
][order(mean_rank)]

fwrite(
  per_output,
  file.path(CONFIG$paths$output, "pseudotime_cell_completeness_by_output.csv")
)
fwrite(
  algorithm_summary,
  file.path(CONFIG$paths$output, "pseudotime_cell_completeness_by_algorithm.csv")
)
fwrite(
  threshold_summary,
  file.path(CONFIG$paths$output, "pseudotime_cell_completeness_threshold_sensitivity.csv")
)
fwrite(
  full_cell_core_summary,
  file.path(CONFIG$paths$output, "pseudotime_full_cell_core_summary.csv")
)

print(algorithm_summary)
print(full_cell_core_summary)
