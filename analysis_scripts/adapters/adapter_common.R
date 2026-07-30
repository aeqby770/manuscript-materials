
# path bootstrap (relative; set LLR_ASSETS or LLR_PIPELINE to override)
.llr_cmd <- commandArgs(trailingOnly = FALSE)
.llr_file <- sub("^--file=", "", grep("^--file=", .llr_cmd, value = TRUE)[1])
if (!length(.llr_file) || is.na(.llr_file) || !nzchar(.llr_file)) {
  stop("Run this script via Rscript so --file= is set")
}
.llr_here <- dirname(normalizePath(.llr_file, winslash = "/", mustWork = TRUE))
source(file.path(.llr_here, "../config/llr_bootstrap.R"), local = FALSE)

# Shared paths for benchmark adapters (source once from each run_*.R)
if (!exists("CONFIG") || is.null(CONFIG)) {
  llr_source_config()
}
EXPORT_R_SCRIPT <- file.path(CONFIG$dir$pipeline, "adapters", "00_export_R.R")
RSCRIPT <- CONFIG$exe$rscript
PYTHON <- CONFIG$exe$python

EnsureRExport <- function(ds) {
  rd <- file.path(CONFIG$dir$r_export, ds)
  if (!file.exists(file.path(rd, "meta.csv"))) {
    if (!file.exists(EXPORT_R_SCRIPT)) stop("Missing export script: ", EXPORT_R_SCRIPT)
    status <- system2(RSCRIPT, args = c(EXPORT_R_SCRIPT, ds))
    if (status != 0) stop("00_export_R failed for ", ds)
  }
  rd
}

ClusterCol <- function(meta) {
  for (col in c("clusters", "cluster", "cell_type", "celltype", "CellType")) {
    if (col %in% names(meta)) return(col)
  }
  names(meta)[1]
}
