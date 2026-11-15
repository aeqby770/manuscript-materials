# Adaptive two-stage temporal permutation and FakeLR wrappers.

make_adaptive_fixture <- function() {
  skip_if_not_installed("Matrix")
  set.seed(42)
  n_genes <- 30L
  n_cells <- 200L
  n_types <- 2L
  x <- Matrix::rsparsematrix(n_genes, n_cells, density = 0.35)
  x@x <- abs(x@x)
  ct <- sample.int(n_types, n_cells, replace = TRUE)
  pt <- seq(0, 1, length.out = n_cells)
  lw <- matrix(1, nrow = n_cells, ncol = 1L)
  lr_db <- list(
    ligand_rows = list(1L, 5L, 10L),
    receptor_rows = list(2L, 6L, 11L)
  )
  list(x = x, ct = ct, pt = pt, lw = lw, lr_db = lr_db, n_bins = 4L)
}

test_that("ComputeAdaptiveTemporalPermutation returns adaptive audit columns", {
  f <- make_adaptive_fixture()
  out <- ComputeAdaptiveTemporalPermutation(
    f$x, f$ct, f$lr_db, f$pt, f$lw,
    n_bins = f$n_bins,
    min_pct = 0.0,
    min_n_eff = 3.0,
    min_valid_bins = 2L,
    min_valid_fraction = 0.1,
    n_perms_screen = 19L,
    n_perms_confirm = 29L,
    confirm_top_k = 2L,
    q_confirm_threshold = 1.0,
    seed = 7L,
    min_perm_observed = 1L,
    min_perm_coverage = 0
  )
  traj <- out$trajectories
  expect_true(is.data.frame(traj))
  expect_true(all(c(
    "perm_stage", "n_perm_used", "min_achievable_p", "monte_carlo_se",
    "adaptive_perm_status", "p_temporal_screen", "adaptive_perm_config_hash"
  ) %in% colnames(traj)))
  expect_true(is.list(out$adaptive_perm_audit))
  expect_equal(out$adaptive_perm_audit$n_perms_screen, 19L)
  expect_equal(out$adaptive_perm_audit$n_perms_confirm, 29L)
  expect_true(out$adaptive_perm_audit$n_candidates >= 1L)
  expect_equal(
    out$adaptive_perm_audit$n_confirmed + out$adaptive_perm_audit$n_screened_out,
    out$adaptive_perm_audit$n_trajectories
  )
})

test_that("adaptive temporal permutation assigns confirm vs screened_out stages", {
  f <- make_adaptive_fixture()
  out <- ComputeAdaptiveTemporalPermutation(
    f$x, f$ct, f$lr_db, f$pt, f$lw,
    n_bins = f$n_bins,
    min_pct = 0.0,
    min_n_eff = 3.0,
    min_valid_bins = 2L,
    min_valid_fraction = 0.1,
    n_perms_screen = 17L,
    n_perms_confirm = 31L,
    confirm_top_k = 1L,
    q_confirm_threshold = 0.0,
    seed = 11L,
    min_perm_observed = 1L,
    min_perm_coverage = 0
  )
  traj <- out$trajectories
  confirmed <- traj[traj$adaptive_perm_status == "confirmed", , drop = FALSE]
  screened <- traj[traj$adaptive_perm_status == "screened_out", , drop = FALSE]
  expect_true(nrow(confirmed) >= 1L)
  expect_true(nrow(screened) >= 1L)
  expect_true(all(confirmed$perm_stage == "confirm"))
  expect_equal(confirmed$n_perm_used, confirmed$n_perm_observed)
  expect_equal(screened$n_perm_used, screened$n_perm_observed)
  expect_equal(confirmed$min_achievable_p, 1 / (confirmed$n_perm_observed + 1))
  expect_equal(screened$min_achievable_p, 1 / (screened$n_perm_observed + 1))
  if (any(is.finite(confirmed$p_temporal))) {
    se <- sqrt(confirmed$p_temporal * (1 - confirmed$p_temporal) / confirmed$n_perm_used)
    expect_equal(confirmed$monte_carlo_se, se, tolerance = 1e-10)
  }
})

test_that("ComputeAdaptiveTrajectoryFakeLRNull runs two-stage FakeLR calibration", {
  f <- make_adaptive_fixture()
  out <- NULL
  out <- ComputeAdaptiveTrajectoryFakeLRNull(
      f$x, f$ct, f$lr_db,
      candidate_fake_genes = 15:25,
      pseudotime = f$pt,
      lineage_weights = f$lw,
      n_bins = f$n_bins,
      n_fake_screen = 100L,
      n_fake_confirm = 100L,
      confirm_top_k = 1L,
      q_confirm_threshold = 1.0,
      n_strata = 3L,
      caliper_k = 5L,
      min_pct = 0,
      min_n_eff = 3.0,
      min_valid_bins = 2L,
      min_valid_fraction = 0.1,
      seed = 7L,
      min_fake_observed = 1L,
      min_fake_coverage = 0
  )
  traj <- out$trajectories
  expect_true(all(c(
    "perm_stage", "n_perm_used", "min_achievable_p", "monte_carlo_se",
    "adaptive_perm_status", "p_fake_screen"
  ) %in% colnames(traj)))
  expect_true(out$adaptive_perm_audit$n_candidates >= 1L)
  expect_true(any(traj$adaptive_perm_status == "confirmed"))
  expect_equal(traj$n_perm_used, traj$n_fake_observed)
  expect_true(all(c("screen", "confirm") %in% out$matching_audit$adaptive_stage))
  expect_true(length(out$stage_fake_null_config_hashes) >= 2L)
})
