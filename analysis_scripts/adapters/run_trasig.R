#!/usr/bin/env Rscript

# path bootstrap (relative; set LLR_ASSETS or LLR_PIPELINE to override)
.llr_cmd <- commandArgs(trailingOnly = FALSE)
.llr_file <- sub("^--file=", "", grep("^--file=", .llr_cmd, value = TRUE)[1])
if (!length(.llr_file) || is.na(.llr_file) || !nzchar(.llr_file)) {
  stop("Run this script via Rscript so --file= is set")
}
.llr_here <- dirname(normalizePath(.llr_file, winslash = "/", mustWork = TRUE))
source(file.path(.llr_here, "../config/llr_bootstrap.R"), local = FALSE)

# TraSig: slingshot prep (Python, no rpy2) + optional TraSig main.py
llr_source_config()
args <- commandArgs(trailingOnly = TRUE)
ds_name <- args[1]
if (is.na(ds_name)) stop("Usage: Rscript run_trasig.R <dataset_name> [output_dir]")
out_dir <- if (length(args) >= 2 && !is.na(args[2])) args[2] else file.path(CONFIG$dir$algorithm_out, "TraSig", ds_name)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

h5ad <- DatasetPath(ds_name)
status <- sys::exec_wait(
  CONFIG$exe$python,
  args = c(
    file.path(llr_pipeline_dir(), "adapters/run_trasig_prep.py"),
    h5ad, out_dir, ds_name
  ),
  std_out = file.path(CONFIG$dir$log, sprintf("trasig_prep_%s.log", ds_name)),
  std_err = file.path(CONFIG$dir$log, sprintf("trasig_prep_%s.err", ds_name))
)
prep_ok <- file.exists(file.path(out_dir, "hid_var.pickle"))
res_files <- list.files(file.path(out_dir, "results"), pattern = "permutation_results\\.pickle$", full.names = TRUE)
res_ok <- length(res_files) > 0L
if (status == 0 && res_ok) {
  writeLines(as.character(Sys.time()), file.path(out_dir, "trasig_done.txt"))
  cat("OK (full)\n")
  quit(status = 0)
}
cat("FAIL\n")
quit(status = 1)
