#!/usr/bin/env Rscript

CONFIG <- list(
  paths = list(
    result = "F:/biotrade/0.64LLR/cache/linealinkr_holdout/GSE216606_celltag_lsk/linealinkr_complete.rds",
    cells = "F:/biotrade/0.64LLR/results/biology/gse216606/gse216606_linealinkr_cells.csv",
    summary = "F:/biotrade/0.64LLR/results/biology/gse216606/gse216606_linealinkr_run_summary.csv",
    root_candidates = "F:/biotrade/0.64LLR/results/biology/gse216606/gse216606_root_candidates.csv"
  ),
  threads = parallel::detectCores(),
  seed = 12345L,
  plot = list(),
  expected = list(
    dataset = "GSE216606_celltag_lsk",
    package_version = "0.6.4",
    benchmark_protocol = "linealinkr064_submission",
    result_revision = "linealinkr064_terminal_leaf_endpoint"
  )
)

set.seed(CONFIG$seed)
Sys.setenv(OMP_NUM_THREADS = CONFIG$threads)
data.table::setDTthreads(CONFIG$threads)

result <- readRDS(CONFIG$paths$result)
stopifnot(
  identical(result$dataset, CONFIG$expected$dataset),
  identical(as.character(result$package_version), CONFIG$expected$package_version),
  identical(as.character(result$benchmark_protocol), CONFIG$expected$benchmark_protocol),
  identical(as.character(result$result_revision), CONFIG$expected$result_revision),
  length(result$cells) == result$n_cells,
  length(result$pseudotime) == result$n_cells,
  nrow(result$lineage_weights) == result$n_cells,
  all(is.finite(result$lineage_weights))
)

inferred_lineage <- max.col(result$lineage_weights, ties.method = "first")
predictions <- data.table::data.table(
  cell_id = result$cells,
  pseudotime = as.numeric(result$pseudotime),
  inferred_lineage = inferred_lineage,
  max_lineage_weight = result$lineage_weights[
    cbind(seq_len(nrow(result$lineage_weights)), inferred_lineage)
  ]
)

run_summary <- data.table::data.table(
  dataset = result$dataset,
  package_version = as.character(result$package_version),
  benchmark_protocol = result$benchmark_protocol,
  result_revision = result$result_revision,
  config_hash = result$config_hash,
  n_cells = result$n_cells,
  n_genes = result$n_genes,
  n_types = result$n_types,
  n_lineages = result$n_lineages,
  start_cell = result$start_cell,
  start_cell_id = result$cells[result$start_cell],
  n_root_candidates = nrow(result$root_candidates),
  has_alternative_root = !is.null(result$root_alternative),
  root_support = result$root_support,
  direction_margin = result$direction_margin,
  direction_status = result$direction_status,
  root_evidence_mode = result$root_evidence_mode,
  root_input_source = result$root_input_source,
  celltype_column = result$celltype_column,
  lineage_cluster_source = result$lineage_cluster_source,
  input_content_schema = result$input_manifest$content_schema,
  input_matrix_source = result$input_manifest$matrix_source,
  input_matrix_nnz = result$input_manifest$matrix_nnz,
  input_matrix_sum = result$input_manifest$matrix_sum
)

dir.create(dirname(CONFIG$paths$cells), recursive = TRUE, showWarnings = FALSE)
data.table::fwrite(predictions, CONFIG$paths$cells)
data.table::fwrite(run_summary, CONFIG$paths$summary)
data.table::fwrite(data.table::as.data.table(result$root_candidates), CONFIG$paths$root_candidates)

cat(sprintf(
  "GSE216606_EXPORT_OK cells=%d lineages=%d direction=%s\n",
  nrow(predictions),
  result$n_lineages,
  result$direction_status
))
