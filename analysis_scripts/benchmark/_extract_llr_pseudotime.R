# Minimal extract: write LLR pseudotime CSVs (avoid printing large objects)
# path bootstrap (relative; set LLR_ASSETS to override)
.llr_cmd <- commandArgs(trailingOnly = FALSE)
.llr_file <- sub("^--file=", "", grep("^--file=", .llr_cmd, value = TRUE)[1])
if (!length(.llr_file) || is.na(.llr_file) || !nzchar(.llr_file)) {
  stop("Run via Rscript so --file= is set")
}
.llr_here <- dirname(normalizePath(.llr_file, winslash = "/", mustWork = TRUE))
source(file.path(.llr_here, "../config/llr_bootstrap.R"), local = FALSE)
.llr_assets <- llr_assets_dir()

set.seed(12345)
Sys.setenv(OMP_NUM_THREADS = as.character(parallel::detectCores()))
data.table::setDTthreads(parallel::detectCores())

CONFIG <- list(
  paths = list(
    workspace = file.path(.llr_assets, "cache/linealinkr_workspace"),
    out_dir = file.path(.llr_assets, "results/paper_tables/supplement/tables/temporal_only_orthogonal_validation")
  ),
  datasets = c("pancreas", "dentategyrus")
)

dir.create(CONFIG$paths$out_dir, recursive = TRUE, showWarnings = FALSE)

for (ds in CONFIG$datasets) {
  rds_path <- file.path(CONFIG$paths$workspace, ds, "linealinkr_complete.rds")
  rds <- readRDS(rds_path)
  cells <- as.character(rds$cells)
  pt <- as.numeric(rds$pseudotime)
  out <- data.frame(cell = cells, pseudotime = pt, stringsAsFactors = FALSE)
  out_path <- file.path(CONFIG$paths$out_dir, paste0(ds, "_llr_pseudotime.csv"))
  utils::write.csv(out, out_path, row.names = FALSE)
  message("wrote ", out_path, " n=", nrow(out))
  rm(rds)
  gc(FALSE)
}
message("done")
