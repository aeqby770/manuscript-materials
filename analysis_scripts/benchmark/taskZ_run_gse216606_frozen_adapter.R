#!/usr/bin/env Rscript

CONFIG_ADAPTER <- list(
  runner = "F:/biotrade/0.64LLR/pipeline/manuscript/02_run_linealinkr.R",
  dataset = "GSE216606_celltag_lsk",
  canonical = "F:/biotrade/0.64LLR/cache/biology/GSE216606_celltag_lsk_canonical.h5ad",
  output = "F:/biotrade/0.64LLR/cache/linealinkr_holdout",
  species = "mouse",
  threads = parallel::detectCores(),
  seed = 12345L
)

set.seed(CONFIG_ADAPTER$seed)
Sys.setenv(OMP_NUM_THREADS = CONFIG_ADAPTER$threads)
data.table::setDTthreads(CONFIG_ADAPTER$threads)

stopifnot(file.exists(CONFIG_ADAPTER$canonical))

runner_lines <- readLines(
  CONFIG_ADAPTER$runner,
  encoding = "UTF-8",
  warn = FALSE
)
source_marker <- which(trimws(runner_lines) == "llr_source_config()")
stopifnot(length(source_marker) == 1L)

adapter_lines <- c(
  "llr_source_config()",
  sprintf(
    "datasets[[\"%s\"]] <- \"%s\"",
    CONFIG_ADAPTER$dataset,
    CONFIG_ADAPTER$canonical
  ),
  sprintf(
    "CONFIG$dataset_narrative_registry[[\"%s\"]] <- normalizePath(\"%s\", winslash = \"/\", mustWork = TRUE)",
    CONFIG_ADAPTER$dataset,
    CONFIG_ADAPTER$canonical
  ),
  sprintf(
    "CONFIG$dataset_species[[\"%s\"]] <- \"%s\"",
    CONFIG_ADAPTER$dataset,
    CONFIG_ADAPTER$species
  ),
  sprintf(
    "CONFIG$dir$linealinkr <- \"%s\"",
    CONFIG_ADAPTER$output
  )
)
runner_lines <- append(
  runner_lines[-source_marker],
  adapter_lines,
  after = source_marker - 1L
)
run_marker <- which(trimws(runner_lines) == "RunLineaLinkRBenchmark()")
stopifnot(length(run_marker) == 1L)
runner_lines <- append(
  runner_lines,
  "filter_datasets <- CONFIG_ADAPTER$dataset",
  after = run_marker - 1L
)

eval(parse(text = runner_lines), envir = globalenv())
