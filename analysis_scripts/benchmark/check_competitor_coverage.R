#!/usr/bin/env Rscript
# Derive manuscript coverage tables from the current fair-comparison dataset set.
suppressPackageStartupMessages(library(data.table))

SCRIPT_PATH <- normalizePath(
  sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1L]),
  winslash = "/", mustWork = TRUE
)
PIPELINE_DIR <- normalizePath(file.path(dirname(SCRIPT_PATH), ".."), winslash = "/", mustWork = TRUE)
source(file.path(PIPELINE_DIR, "config", "00_config.R"))

set.seed(CONFIG$seed)
Sys.setenv(OMP_NUM_THREADS = CONFIG$n_threads)
data.table::setDTthreads(CONFIG$n_threads)

REG <- file.path(CONFIG$dir$paper_tables, "_registry")
FAIR_MATRIX_PATH <- file.path(REG, "fair_comparison_matrix.csv")
RUN_REGISTRY_PATH <- file.path(
  REG,
  sprintf("benchmark_results_%s.csv", CONFIG$llr$benchmark_protocol)
)
COMPETITORS <- setdiff(BENCHMARK_PSEUDOTIME_ALGORITHMS, "LineaLinkR")
FAIR_REQUIRED_COLUMNS <- c(
  "algorithm", "dataset", "dataset_role", "output_status",
  "current_run_status", "benchmark_protocol", "result_revision"
)
RUN_REQUIRED_COLUMNS <- c(
  "algorithm", "dataset", "status", "benchmark_protocol", "result_revision"
)

stopifnot(file.exists(FAIR_MATRIX_PATH), file.exists(RUN_REGISTRY_PATH))
fair <- fread(FAIR_MATRIX_PATH, showProgress = FALSE)
registry <- fread(RUN_REGISTRY_PATH, showProgress = FALSE)
stopifnot(
  all(FAIR_REQUIRED_COLUMNS %in% names(fair)),
  all(RUN_REQUIRED_COLUMNS %in% names(registry)),
  nrow(fair) == length(BENCHMARK_PSEUDOTIME_ALGORITHMS) * length(BENCHMARK_ORDERING_DATASETS),
  uniqueN(fair, by = c("algorithm", "dataset")) == nrow(fair),
  all(fair$benchmark_protocol == CONFIG$llr$benchmark_protocol),
  all(fair$result_revision == CONFIG$llr$result_revision)
)

coverage <- fair[algorithm %chin% COMPETITORS]
stopifnot(
  nrow(coverage) == length(COMPETITORS) * length(BENCHMARK_ORDERING_DATASETS),
  setequal(coverage$dataset, BENCHMARK_ORDERING_DATASETS),
  setequal(coverage$algorithm, COMPETITORS),
  all(coverage$dataset_role %chin% c("complete_core", "scale_extension"))
)

registry <- registry[
  benchmark_protocol == CONFIG$llr$benchmark_protocol &
    result_revision == CONFIG$llr$result_revision &
    algorithm %chin% COMPETITORS &
    dataset %chin% BENCHMARK_ORDERING_DATASETS,
  .(
    algorithm = as.character(algorithm),
    dataset = as.character(dataset),
    registered_run_status = as.character(status)
  )
]
stopifnot(!anyDuplicated(registry, by = c("algorithm", "dataset")))
coverage <- merge(
  coverage,
  registry,
  by = c("algorithm", "dataset"),
  all.x = TRUE,
  sort = FALSE
)
stopifnot(all(
  (is.na(coverage$current_run_status) & is.na(coverage$registered_run_status)) |
    (!is.na(coverage$current_run_status) & !is.na(coverage$registered_run_status) &
       coverage$current_run_status == coverage$registered_run_status)
))

core <- coverage[dataset_role == "complete_core"]
stopifnot(
  nrow(core) == length(COMPETITORS) * length(BENCHMARK_COMPLETE_CORE_DATASETS),
  all(core$output_status == "complete_output"),
  all(core$current_run_status == "OK"),
  all(core$registered_run_status == "OK")
)

coverage[, coverage_status := fifelse(
  output_status == "complete_output" & current_run_status == "OK",
  "complete_current_run",
  output_status
)]
setorder(coverage, dataset_role, algorithm, dataset)
summary <- coverage[, .(
  n_datasets = .N,
  n_complete_current_run = sum(coverage_status == "complete_current_run"),
  n_noncomplete = sum(coverage_status != "complete_current_run")
), by = .(dataset_role, algorithm)]
setorder(summary, dataset_role, algorithm)

fwrite(coverage, file.path(REG, "competitor_coverage_matrix.csv"))
fwrite(summary, file.path(REG, "competitor_coverage_summary.csv"))
cat("CURRENT_COMPETITOR_COVERAGE_WRITTEN\n")
