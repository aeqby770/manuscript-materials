# Calibration denominators, evidence gates, and observed permutation counts.

test_that("lag null p uses observed permutation count", {
  side <- data.frame(
    pair_idx = rep(1L, 10L),
    source = rep(1L, 10L),
    target = rep(2L, 10L),
    lineage = rep(1L, 10L),
    bin = rep(1:5, each = 2L),
    score = rep(0.5, 10L),
    ligand_score = c(0.1, 0.3, 1.0, 0.4, 0.2, 0.1, 0.2, 0.5, 0.9, 0.3),
    receptor_score = c(0.2, 0.4, 0.5, 0.8, 1.0, 0.1, 0.3, 0.4, 0.6, 1.0),
    stringsAsFactors = FALSE
  )
  pt <- seq(0, 1, length.out = 40L)
  ct <- rep(1:2, each = 20L)
  lw <- matrix(1, nrow = 40L, ncol = 1L)
  lag <- ComputeLigandTargetLag(
    side, lw, pt, ct, method = "peak", n_bins = 5L,
    compute_lag_null = TRUE, n_lag_perm = 49L, lag_null_seed = 2L
  )
  expect_true("n_lag_perm_observed" %in% colnames(lag))
  expect_true(lag$n_lag_perm_observed[1L] >= 1L)
  expect_true(lag$n_lag_perm_observed[1L] <= 49L)
})

test_that("temporal_prioritized is not aliased to variation-only column", {
  df <- data.frame(
    pair_idx = 1:2,
    source = 1L,
    target = 2L,
    lineage = 1L,
    n_eff_cells = 10,
    lineage_entropy = 0.2,
    stringsAsFactors = FALSE
  )
  fake_cal <- df[, c("pair_idx", "source", "target", "lineage"), drop = FALSE]
  fake_cal$p_fake <- c(0.01, 0.5)
  fake_cal$q_fake <- c(0.02, 0.6)
  fake_cal$trajectory_evaluable <- TRUE
  temp_cal <- fake_cal
  temp_cal$p_temporal <- c(0.01, 0.01)
  temp_cal$q_temporal <- c(0.02, 0.02)
  temp_cal$trajectory_stat_name <- c("mean_score", "mean_score")
  temp_cal$calibration_hypothesis <- "trajectory_strength"
  out <- MergeEvidence(
    df,
    fake_calibration = fake_cal,
    temporal_calibration = temp_cal
  )
  expect_true(out$temporal_prioritized[1])
  expect_true(out$trajectory_strength_prioritized[1])
  expect_false(isTRUE(out$temporal_variation_prioritized[1]))
})

test_that("dual_layer_evaluable is exposed when both calibrations requested", {
  df <- data.frame(
    pair_idx = 1L, source = 1L, target = 2L, lineage = 1L,
    n_eff_cells = 10, lineage_entropy = 0.2, stringsAsFactors = FALSE
  )
  cal <- df[, c("pair_idx", "source", "target", "lineage"), drop = FALSE]
  cal$p_fake <- 0.01
  cal$q_fake <- 0.02
  cal$trajectory_evaluable <- TRUE
  cal$p_temporal <- 0.02
  cal$q_temporal <- 0.04
  out <- MergeEvidence(df, fake_calibration = cal, temporal_calibration = cal)
  expect_true("dual_layer_evaluable" %in% colnames(out))
  expect_true(out$dual_layer_evaluable[1])
})
