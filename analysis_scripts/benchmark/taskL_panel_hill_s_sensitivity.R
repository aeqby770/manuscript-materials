# Full-panel post-hoc Hill-S multiplier sensitivity from audited LineaLinkR RDS outputs.
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
Sys.setenv(OMP_NUM_THREADS = parallel::detectCores())
data.table::setDTthreads(parallel::detectCores())

library(data.table)

CONFIG <- list(
  paths = list(
    workspace = file.path(.llr_assets, "cache/linealinkr_workspace"),
    panel = file.path(.llr_assets, "audit/panel_completion.tsv"),
    out = file.path(.llr_assets, "results/paper_tables/supplement/tables/sensitivity_sweep")
  ),
  multipliers = c(0.5, 1, 2, 5),
  top_n = 30L,
  seed = 12345L
)

dir.create(CONFIG$paths$out, recursive = TRUE, showWarnings = FALSE)

panel <- fread(CONFIG$paths$panel)
stopifnot(nrow(panel) > 0L, all(panel$status == "complete_current"))
datasets <- panel$dataset
rds_paths <- file.path(CONFIG$paths$workspace, datasets, "linealinkr_complete.rds")
stopifnot(all(file.exists(rds_paths)))
available_datasets <- datasets

score_schema <- data.table(
  pair_idx = integer(),
  source = character(),
  target = character(),
  lineage = integer(),
  bin = integer(),
  score = numeric()
)
pair_schema <- data.table(
  pair_idx = integer(),
  ligand_genes = character(),
  receptor_genes = character()
)

results <- lapply(available_datasets, \(dataset_id) {
  result_path <- file.path(CONFIG$paths$workspace, dataset_id, "linealinkr_complete.rds")
  result <- readRDS(result_path)
  expected <- panel[dataset == dataset_id]
  stopifnot(
    nrow(expected) == 1L,
    identical(as.character(result$dataset), dataset_id),
    identical(as.character(result$package_version), expected$expected_package_version),
    identical(as.character(result$benchmark_protocol), expected$expected_benchmark_protocol),
    identical(as.character(result$result_revision), expected$expected_result_revision)
  )
  scores <- rbindlist(list(score_schema, as.data.table(result$temporal_lr$scores)), fill = TRUE)
  pair_map <- rbindlist(list(pair_schema, as.data.table(result$pair_gene_map)), fill = TRUE)
  scores <- scores[is.finite(score) & score > 0 & score < 1]
  scores <- merge(scores, pair_map, by = "pair_idx", all.x = TRUE)
  multipliers <- data.table(hill_multiplier = CONFIG$multipliers, join_key = 1L)
  scores[, join_key := 1L]
  expanded <- merge(scores, multipliers, by = "join_key", allow.cartesian = TRUE)
  expanded[, score_k := score / (hill_multiplier - (hill_multiplier - 1) * score)]
  pair_scores <- expanded[, .(
    mean_score = mean(score_k),
    median_score = median(score_k),
    active_fraction = mean(score_k > 0.2)
  ), by = .(hill_multiplier, ligand_genes, receptor_genes)]
  pair_scores[, pair_id := paste(ligand_genes, receptor_genes, sep = "___")]
  setorder(pair_scores, hill_multiplier, -mean_score, pair_id)
  pair_scores[, pair_rank := seq_len(.N), by = hill_multiplier]
  top_pairs <- pair_scores[pair_rank <= CONFIG$top_n]
  base_pairs <- unique(top_pairs[hill_multiplier == 1, pair_id])
  summary <- top_pairs[, .(
    n_top = .N,
    n_overlap = sum(pair_id %chin% base_pairs),
    n_union = length(unique(c(pair_id, base_pairs))),
    mean_top_score = mean(mean_score),
    median_top_score = median(median_score),
    mean_active_fraction = mean(active_fraction)
  ), by = hill_multiplier]
  summary[, `:=`(
    dataset = dataset_id,
    jaccard_top30_vs_k1 = n_overlap / n_union
  )]
  coverage <- data.table(
    dataset = dataset_id,
    n_score_rows = nrow(scores),
    n_pairs = uniqueN(scores[, .(ligand_genes, receptor_genes)]),
    n_lineages = uniqueN(scores$lineage),
    n_cells = result$n_cells,
    evaluable = nrow(scores) > 0,
    availability = "output_present"
  )
  list(summary = summary, coverage = coverage)
})

summary <- rbindlist(lapply(results, `[[`, "summary"), fill = TRUE)
coverage <- rbindlist(lapply(results, `[[`, "coverage"), fill = TRUE)
coverage <- rbindlist(list(
  coverage,
  data.table(
    dataset = setdiff(datasets, coverage$dataset),
    n_score_rows = NA_integer_,
    n_pairs = NA_integer_,
    n_lineages = NA_integer_,
    n_cells = NA_integer_,
    evaluable = FALSE,
    availability = "missing_llr_output"
  )
), fill = TRUE)
grid <- CJ(dataset = datasets, hill_multiplier = CONFIG$multipliers)
summary <- merge(grid, summary, by = c("dataset", "hill_multiplier"), all.x = TRUE)
summary <- merge(summary, coverage, by = "dataset", all.x = TRUE)
setorder(summary, dataset, hill_multiplier)

panel_summary <- summary[evaluable == TRUE, .(
  n_datasets = uniqueN(dataset),
  median_jaccard = median(jaccard_top30_vs_k1, na.rm = TRUE),
  q1_jaccard = quantile(jaccard_top30_vs_k1, 0.25, na.rm = TRUE),
  q3_jaccard = quantile(jaccard_top30_vs_k1, 0.75, na.rm = TRUE),
  min_jaccard = min(jaccard_top30_vs_k1, na.rm = TRUE),
  median_mean_top_score = median(mean_top_score, na.rm = TRUE)
), by = hill_multiplier]

fwrite(summary, file.path(CONFIG$paths$out, "panel_hill_s_multiplier_sensitivity.csv"))
fwrite(coverage, file.path(CONFIG$paths$out, "panel_hill_s_evaluable_coverage.csv"))
fwrite(panel_summary, file.path(CONFIG$paths$out, "panel_hill_s_multiplier_summary.csv"))
print(panel_summary)
print(coverage[evaluable == FALSE])
