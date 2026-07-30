#!/usr/bin/env Rscript
# Build retry queue: batch01 jobs still missing expected artifacts (skip successes).
# path bootstrap (relative; set LLR_ASSETS to override)
.llr_cmd <- commandArgs(trailingOnly = FALSE)
.llr_file <- sub("^--file=", "", grep("^--file=", .llr_cmd, value = TRUE)[1])
if (!length(.llr_file) || is.na(.llr_file) || !nzchar(.llr_file)) {
  stop("Run via Rscript so --file= is set")
}
.llr_here <- dirname(normalizePath(.llr_file, winslash = "/", mustWork = TRUE))
source(file.path(.llr_here, "../config/llr_bootstrap.R"), local = FALSE)
.llr_assets <- llr_assets_dir()
llr_source_config()

suppressPackageStartupMessages(library(data.table))

queue_csv <- file.path(CONFIG$dir$paper_tables, "_registry", "batch01_missing_output_queue.csv")
out_csv <- file.path(CONFIG$dir$paper_tables, "_registry", "batch01_retry_queue.csv")
if (!file.exists(queue_csv)) stop("Missing queue: ", queue_csv)

OUT_ALG <- CONFIG$dir$algorithm_out
OUT_MISC <- CONFIG$dir$algorithm_misc

artifact_path <- function(algo, ds) {
  misc <- function(f) file.path(OUT_MISC, f)
  alg <- function(...) file.path(OUT_ALG, algo, ds, ...)
  switch(algo,
    DPT = misc(sprintf("dpt_%s_pseudotime.csv", ds)),
    slingshot = alg("slingshot_pseudotime.csv"),
    SCORPIUS = alg("scorpius_pseudotime.csv"),
    TSCAN = alg("tscan_pseudotime.csv"),
    FateID = alg("fateid_pseudotime.csv"),
    monocle3 = alg("monocle3_cds.rds"),
    Palantir = misc(sprintf("palantir_%s_result.h5ad", ds)),
    MultistageOT = misc(sprintf("multistageot_%s_pseudotime.csv", ds)),
    scFates = misc(sprintf("scfates_%s_pseudotime.csv", ds)),
    pCreode = misc(sprintf("pcreode_%s_pseudotime.csv", ds)),
    TSCAN_official = alg("tscan_official_pseudotime.csv"),
    NA_character_
  )
}

artifact_ok <- function(algo, ds) {
  p <- artifact_path(algo, ds)
  if (is.na(p) || !nzchar(p)) return(FALSE)
  file.exists(p) && !is.na(file.info(p)$size) && file.info(p)$size > 100
}

q <- fread(queue_csv)
q <- q[algorithm != "LineaLinkR"]
q[, has_output := mapply(artifact_ok, algorithm, dataset)]
retry <- q[has_output == FALSE, .(algorithm, dataset)]
fwrite(retry, out_csv)

cat(sprintf("Batch01 retry queue: %d jobs (skipped %d with valid output)\n",
            nrow(retry), sum(q$has_output)))
print(retry[, .N, by = algorithm][order(-N)])
