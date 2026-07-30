script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg), winslash = "/", mustWork = TRUE))
Sys.setenv(LLR_ASSETS = normalizePath(file.path(script_dir, "../.."), winslash = "/", mustWork = TRUE))
source(file.path(script_dir, "../config/llr_bootstrap.R"), local = TRUE)
assets <- llr_assets_dir()

library(data.table)

CONFIG <- list(
  paths = list(
    workspace = file.path(assets, "cache/linealinkr_workspace"),
    output = file.path(assets, "results/paper_tables/supplement/tables/uncertainty")
  ),
  dataset = "pancreas",
  threads = parallel::detectCores(),
  seed = 12345L
)

set.seed(CONFIG$seed)
Sys.setenv(OMP_NUM_THREADS = CONFIG$threads)
data.table::setDTthreads(CONFIG$threads)

ds <- CONFIG$dataset
rds <- readRDS(file.path(CONFIG$paths$workspace, ds, "linealinkr_complete.rds"))
lw <- as.data.frame(rds$lineage_weights)
colnames(lw) <- paste0("L", seq_len(ncol(lw)))
lw$cell <- as.character(rds$cells)
lw$pseudotime <- as.numeric(rds$pseudotime)
utils::write.csv(lw, file.path(CONFIG$paths$output, paste0(ds, "_lineage_weights.csv")), row.names = FALSE)
message("wrote lineage weights n=", nrow(lw), " ncol=", ncol(rds$lineage_weights))
