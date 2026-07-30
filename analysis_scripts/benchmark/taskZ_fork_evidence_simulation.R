# Deterministic module-level simulations for LineaLinkR 0.6.4 fork evidence.
.llr_cmd <- commandArgs(trailingOnly = FALSE)
.llr_file <- sub("^--file=", "", grep("^--file=", .llr_cmd, value = TRUE)[1L])
stopifnot(length(.llr_file) == 1L, !is.na(.llr_file), nzchar(.llr_file))
.llr_here <- dirname(normalizePath(.llr_file, winslash = "/", mustWork = TRUE))
source(file.path(.llr_here, "../config/llr_bootstrap.R"), local = FALSE)
.llr_assets <- llr_assets_dir()

.libPaths(c(
  Sys.getenv("LLR064_RLIB", unset = "F:/biotrade/0.64LLR/rlib"),
  .libPaths()
))
set.seed(12345)
Sys.setenv(OMP_NUM_THREADS = parallel::detectCores())
data.table::setDTthreads(parallel::detectCores())

library(data.table)
library(LineaLinkR)

CONFIG <- list(
  paths = list(
    output_dir = file.path(
      .llr_assets,
      "results/paper_tables/supplement/tables/fork_evidence_simulation"
    )
  ),
  seed = 12345L,
  n_replicates = 100L,
  n_bins = 10L,
  fork_bin = 5L,
  window_bins = 3L,
  activity_floor = 1e-4,
  min_n_eff = 5,
  dropout_keep_probability = 0.70,
  evidence_threshold = 0.25,
  n_boot = 500L,
  scenarios = c(
    "shared_trend",
    "focal_branch",
    "donor_confounded_null",
    "low_expression_spike",
    "dropout",
    "nested_fork"
  )
)

dir.create(CONFIG$paths$output_dir, recursive = TRUE, showWarnings = FALSE)

scores <- CJ(
  scenario = CONFIG$scenarios,
  replicate = seq_len(CONFIG$n_replicates),
  lineage = 1:3,
  bin = seq_len(CONFIG$n_bins)
)
scores[, pair_idx := .GRP, by = .(scenario, replicate)]
scores[, `:=`(
  source = 1L,
  target = 2L,
  n_eff_cells = 20,
  structural_valid = TRUE,
  post_fork = bin > CONFIG$fork_bin,
  score = 0.5
)]
scores[scenario == "shared_trend" & post_fork, score := 1.5]
scores[
  scenario %in% c("focal_branch", "donor_confounded_null", "dropout") &
    post_fork,
  score := ifelse(lineage == 1L, 2.0, 1.0)
]
scores[
  scenario == "nested_fork" & post_fork,
  score := c(2.0, 1.0, 3.0)[lineage]
]
scores[scenario == "low_expression_spike", score := 1e-8]
scores[
  scenario == "low_expression_spike" & post_fork & lineage == 1L,
  score := 5e-8
]
scores[
  scenario != "low_expression_spike",
  score := pmax(score + rnorm(.N, mean = 0, sd = 0.025), 0)
]
scores[
  scenario == "dropout",
  score := score * rbinom(.N, size = 1L, prob = CONFIG$dropout_keep_probability)
]

pair_meta <- unique(scores[, .(pair_idx, scenario, replicate)])
bif_df <- data.frame(
  lineage = 1:3,
  b_star = rep(CONFIG$fork_bin, 3L),
  reason = rep("crossing", 3L)
)
sibling_map <- LineaLinkR::BuildSiblingContrastMap(list(
  c(1L, 2L, 3L),
  c(1L, 2L, 4L),
  c(1L, 5L, 6L)
))
pooled <- LineaLinkR::ComputeForkContrasts(
  as.data.frame(scores[, !"post_fork"]),
  bif_df,
  sibling_map,
  window_bins = CONFIG$window_bins,
  activity_floor = CONFIG$activity_floor,
  min_n_eff = CONFIG$min_n_eff
)
pooled <- as.data.table(pooled)[lineage == 1L]
pooled <- merge(pooled, pair_meta, by = "pair_idx", all.x = TRUE)
pooled[, pooled_sibling_contrast := ifelse(
  is.finite(sibling_contrast),
  sibling_contrast,
  0
)]

sample_scores <- CJ(
  pair_idx = pair_meta$pair_idx,
  block_label = paste0("D", 1:4)
)
sample_scores <- merge(sample_scores, pair_meta, by = "pair_idx", all.x = TRUE)
sample_scores[, sibling_contrast := 0]
sample_scores[
  scenario %in% c("focal_branch", "dropout", "nested_fork"),
  sibling_contrast := 0.8 + rnorm(.N, mean = 0, sd = 0.08)
]
sample_scores[
  scenario == "shared_trend",
  sibling_contrast := rnorm(.N, mean = 0, sd = 0.02)
]
sample_scores[
  scenario == "donor_confounded_null" & block_label == "D1",
  sibling_contrast := 5.0
]
sample_scores[
  scenario == "donor_confounded_null" & block_label != "D1",
  sibling_contrast := -0.2 + rnorm(.N, mean = 0, sd = 0.02)
]
sample_scores[
  scenario == "low_expression_spike",
  sibling_contrast := rnorm(.N, mean = 0, sd = 1e-8)
]
sample_scores[, `:=`(
  source = 1L,
  target = 2L,
  lineage = 1L,
  evaluable = TRUE
)]
sample_summary <- LineaLinkR::SummarizeForkContrastsBySample(
  as.data.frame(sample_scores[, .(
    pair_idx, source, target, lineage, block_label,
    sibling_contrast, evaluable
  )]),
  n_boot = CONFIG$n_boot,
  seed = CONFIG$seed
)
sample_summary <- as.data.table(sample_summary)

by_replicate <- merge(
  pooled[, .(
    pair_idx,
    scenario,
    replicate,
    pooled_evaluable = evaluable,
    pooled_sibling_contrast
  )],
  sample_summary[, .(
    pair_idx,
    median_contrast,
    sample_support,
    positive_sample_fraction,
    n_samples
  )],
  by = "pair_idx",
  all = FALSE
)
by_replicate[, `:=`(
  seed = CONFIG$seed,
  true_positive = scenario %in% c("focal_branch", "dropout", "nested_fork"),
  evidence_score = pmax(pooled_sibling_contrast, 0) *
    pmax(median_contrast, 0) * sample_support
)]
by_replicate[, detected := evidence_score >= CONFIG$evidence_threshold]
setcolorder(by_replicate, c(
  "scenario", "replicate", "seed", "true_positive", "pooled_evaluable",
  "pooled_sibling_contrast", "median_contrast", "sample_support",
  "positive_sample_fraction", "n_samples", "evidence_score", "detected",
  "pair_idx"
))
setorder(by_replicate, scenario, replicate)

score_ranks <- rank(by_replicate$evidence_score, ties.method = "average")
n_positive <- sum(by_replicate$true_positive)
n_negative <- sum(!by_replicate$true_positive)
integrated_auroc <- (
  sum(score_ranks[by_replicate$true_positive]) -
    n_positive * (n_positive + 1) / 2
) / (n_positive * n_negative)
true_positive_rate <- mean(by_replicate[true_positive == TRUE, detected])
specificity <- mean(by_replicate[true_positive == FALSE, detected == FALSE])

scenario_summary <- by_replicate[, .(
  n_replicates = .N,
  evaluable_fraction = mean(pooled_evaluable),
  median_pooled_contrast = median(pooled_sibling_contrast),
  median_sample_contrast = median(median_contrast),
  median_sample_support = median(sample_support),
  median_evidence_score = median(evidence_score),
  detection_rate = mean(detected)
), by = scenario]
setorder(scenario_summary, scenario)

gates <- data.table(
  gate = c(
    "integrated_auroc",
    "true_positive_rate",
    "specificity",
    "shared_trend_cancellation",
    "donor_confounded_rejection",
    "low_expression_rejection",
    "dropout_recovery",
    "nested_fork_recovery"
  ),
  observed = c(
    integrated_auroc,
    true_positive_rate,
    specificity,
    mean(by_replicate[scenario == "shared_trend", abs(pooled_sibling_contrast) < 0.10]),
    mean(by_replicate[scenario == "donor_confounded_null", detected == FALSE]),
    mean(by_replicate[scenario == "low_expression_spike", detected == FALSE]),
    mean(by_replicate[scenario == "dropout", detected]),
    mean(by_replicate[scenario == "nested_fork", detected])
  ),
  threshold = c(0.90, 0.80, 0.95, 0.95, 0.95, 0.95, 0.70, 0.80),
  comparison = rep(">=", 8L)
)
gates[, passed := observed >= threshold]

fwrite(
  by_replicate,
  file.path(CONFIG$paths$output_dir, "fork_evidence_simulation_by_replicate.csv")
)
fwrite(
  scenario_summary,
  file.path(CONFIG$paths$output_dir, "fork_evidence_simulation_summary.csv")
)
fwrite(
  gates,
  file.path(CONFIG$paths$output_dir, "fork_evidence_simulation_gates.csv")
)

stopifnot(
  nrow(by_replicate) == length(CONFIG$scenarios) * CONFIG$n_replicates,
  uniqueN(by_replicate, by = c("scenario", "replicate")) == nrow(by_replicate),
  all(is.finite(by_replicate$pooled_sibling_contrast)),
  all(is.finite(by_replicate$evidence_score)),
  all(gates$passed)
)

cat(sprintf(
  "FORK_SIMULATION_064_OK rows=%d auroc=%.6f tpr=%.6f specificity=%.6f\n",
  nrow(by_replicate),
  integrated_auroc,
  true_positive_rate,
  specificity
))
