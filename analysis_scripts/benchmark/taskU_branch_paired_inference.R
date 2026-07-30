library(data.table)

set.seed(12345)
Sys.setenv(OMP_NUM_THREADS = parallel::detectCores())
data.table::setDTthreads(parallel::detectCores())

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg), winslash = "/", mustWork = TRUE))
Sys.setenv(LLR_ASSETS = normalizePath(file.path(script_dir, "../.."), winslash = "/", mustWork = TRUE))
source(file.path(script_dir, "../config/llr_bootstrap.R"), local = TRUE)
assets <- llr_assets_dir()
.llr_config_env <- new.env(parent = globalenv())
sys.source(file.path(assets, "pipeline", "config", "00_config.R"), envir = .llr_config_env)
PROJECT_CONFIG <- .llr_config_env$CONFIG
BRANCH_DATASETS <- .llr_config_env$BENCHMARK_BRANCH_DATASETS
BRANCH_TERMINAL_EXPECTED_N <- .llr_config_env$BENCHMARK_BRANCH_TERMINAL_EXPECTED_N
TABLE_PROTOCOL_TAG <- gsub("[^A-Za-z0-9_]+", "_", PROJECT_CONFIG$llr$benchmark_protocol)
rm(.llr_config_env)

CONFIG <- list(
  paths = list(
    branch_f1 = file.path(
      PROJECT_CONFIG$dir$table1_main,
      sprintf("table1_branch_f1_%s.csv", TABLE_PROTOCOL_TAG)
    ),
    lineage = file.path(
      PROJECT_CONFIG$dir$table1_main,
      sprintf("table1_lineage_%s.csv", TABLE_PROTOCOL_TAG)
    ),
    output = file.path(
      assets,
      "results/paper_tables/supplement/tables/branch_paired_inference"
    )
  ),
  threads = parallel::detectCores(),
  seed = 12345L,
  benchmark_protocol = PROJECT_CONFIG$llr$benchmark_protocol,
  result_revision = PROJECT_CONFIG$llr$result_revision,
  n_boot = 10000L,
  numeric_tolerance = 1e-12,
  plot = list()
)

dir.create(CONFIG$paths$output, recursive = TRUE, showWarnings = FALSE)

branch <- fread(CONFIG$paths$branch_f1, encoding = "UTF-8")
lineage <- fread(CONFIG$paths$lineage, encoding = "UTF-8")
stopifnot(
  setequal(unique(branch$dataset), BRANCH_DATASETS),
  setequal(unique(lineage$dataset), BRANCH_DATASETS),
  all(branch$n_aligned == unname(BRANCH_TERMINAL_EXPECTED_N[branch$dataset])),
  all(lineage$n_aligned == unname(BRANCH_TERMINAL_EXPECTED_N[lineage$dataset]))
)

branch_wide <- dcast(
  branch,
  dataset ~ algorithm,
  value.var = "Branch_F1"
)[, `:=`(
  metric = "one_to_one_branch_segment_macro_F1",
  LineaLinkR_minus_slingshot = LineaLinkR - slingshot
)]

lineage_wide <- dcast(
  lineage,
  dataset ~ algorithm,
  value.var = "ARI"
)[, `:=`(
  metric = "lineage_ARI",
  LineaLinkR_minus_slingshot = LineaLinkR - slingshot
)]

paired <- rbindlist(
  list(
    branch_wide[, .(dataset, metric, LineaLinkR, slingshot, LineaLinkR_minus_slingshot)],
    lineage_wide[, .(dataset, metric, LineaLinkR, slingshot, LineaLinkR_minus_slingshot)]
  )
)[order(metric, dataset)]

n_datasets <- uniqueN(branch_wide$dataset)
bootstrap_indices <- matrix(
  sample.int(n_datasets, n_datasets * CONFIG$n_boot, replace = TRUE),
  nrow = n_datasets,
  ncol = CONFIG$n_boot
)
sign_matrix <- as.matrix(expand.grid(rep(list(c(-1, 1)), n_datasets)))

branch_diff <- branch_wide$LineaLinkR_minus_slingshot
branch_boot <- colMeans(matrix(branch_diff[bootstrap_indices], nrow = n_datasets))
branch_null <- as.numeric(sign_matrix %*% branch_diff) / n_datasets
branch_wilcox <- wilcox.test(
  branch_wide$LineaLinkR,
  branch_wide$slingshot,
  paired = TRUE,
  alternative = "two.sided",
  exact = TRUE,
  conf.int = FALSE
)

lineage_diff <- lineage_wide$LineaLinkR_minus_slingshot
lineage_boot <- colMeans(matrix(lineage_diff[bootstrap_indices], nrow = n_datasets))
lineage_null <- as.numeric(sign_matrix %*% lineage_diff) / n_datasets
lineage_wilcox <- wilcox.test(
  lineage_wide$LineaLinkR,
  lineage_wide$slingshot,
  paired = TRUE,
  alternative = "two.sided",
  exact = FALSE,
  conf.int = FALSE
)

summary <- data.table(
  metric = c("one_to_one_branch_segment_macro_F1", "lineage_ARI"),
  n_paired_datasets = n_datasets,
  LineaLinkR_mean = c(mean(branch_wide$LineaLinkR), mean(lineage_wide$LineaLinkR)),
  slingshot_mean = c(mean(branch_wide$slingshot), mean(lineage_wide$slingshot)),
  mean_paired_difference = c(mean(branch_diff), mean(lineage_diff)),
  bootstrap_ci95_lo = c(
    unname(quantile(branch_boot, 0.025)),
    unname(quantile(lineage_boot, 0.025))
  ),
  bootstrap_ci95_hi = c(
    unname(quantile(branch_boot, 0.975)),
    unname(quantile(lineage_boot, 0.975))
  ),
  exact_sign_flip_p_two_sided = c(
    mean(abs(branch_null) >= abs(mean(branch_diff)) - CONFIG$numeric_tolerance),
    mean(abs(lineage_null) >= abs(mean(lineage_diff)) - CONFIG$numeric_tolerance)
  ),
  paired_wilcoxon_p_two_sided = c(branch_wilcox$p.value, lineage_wilcox$p.value),
  n_positive_differences = c(sum(branch_diff > 0), sum(lineage_diff > 0)),
  n_zero_differences = c(sum(branch_diff == 0), sum(lineage_diff == 0)),
  n_negative_differences = c(sum(branch_diff < 0), sum(lineage_diff < 0)),
  bootstrap_seed = CONFIG$seed,
  bootstrap_replicates = CONFIG$n_boot,
  inference_unit = "dataset"
)

fwrite(
  paired,
  file.path(CONFIG$paths$output, "branch_metrics_paired_by_dataset.csv")
)
fwrite(
  summary,
  file.path(CONFIG$paths$output, "branch_metrics_paired_inference_summary.csv")
)

print(summary)
