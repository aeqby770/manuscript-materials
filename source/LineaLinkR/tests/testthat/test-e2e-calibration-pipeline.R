# End-to-end calibration pipeline without synthetic p-values.

make_e2e_fixture <- function(seed = 42L, n_cells = 200L) {
  skip_if_not_installed("Matrix")
  set.seed(seed)
  n_genes <- 30L
  n_types <- 2L
  n_bins <- 4L
  x <- Matrix::rsparsematrix(n_genes, n_cells, density = 0.35)
  x@x <- abs(x@x)
  ct <- sample.int(n_types, n_cells, replace = TRUE)
  pt <- seq(0, 1, length.out = n_cells)
  lw <- matrix(1, nrow = n_cells, ncol = 1L)
  lr_db <- list(
    ligand_rows = list(1L, 5L, 10L),
    receptor_rows = list(2L, 6L, 11L)
  )
  list(x = x, ct = ct, pt = pt, lw = lw, lr_db = lr_db, n_bins = n_bins)
}

test_that("E2E: TemporalLR -> FakeLR -> TemporalPerm -> MergeEvidence", {
  f <- make_e2e_fixture()
  scores <- ComputeTemporalLR(
    f$x, f$ct, f$lr_db, f$pt, f$lw,
    n_bins = f$n_bins, min_pct = 0.0
  )
  skip_if(nrow(scores) == 0L, "no LR scores in fixture")

  fake <- ComputeTrajectoryFakeLRNull(
    f$x, f$ct, f$lr_db,
    candidate_fake_genes = 15:25,
    pseudotime = f$pt,
    lineage_weights = f$lw,
    n_bins = f$n_bins,
    n_fake_per_real = 100L,
    n_strata = 3L,
    caliper_k = 5L,
    min_pct = 0.0,
    min_n_eff = 3.0,
    min_valid_bins = 2L,
    min_valid_fraction = 0.1,
    seed = 7L,
    min_fake_observed = 1L,
    min_fake_coverage = 0
  )

  temporal <- ComputeTemporalPermutation(
    f$x, f$ct, f$lr_db, f$pt, f$lw,
    n_bins = f$n_bins,
    min_pct = 0.0,
    min_n_eff = 3.0,
    min_valid_bins = 2L,
    min_valid_fraction = 0.1,
    n_perms = 49L,
    seed = 11L,
    min_perm_observed = 20L
  )

  merged <- MergeEvidence(
    scores,
    fake_calibration = fake$trajectories,
    temporal_calibration = temporal$trajectories
  )

  expect_equal(nrow(merged), nrow(scores))
  expect_gt(sum(is.finite(merged$p_fake)), 0L)
  expect_gt(sum(is.finite(merged$p_temporal)), 0L)
  expect_true("identity_prioritized" %in% colnames(merged))
  expect_true("n_perm_observed" %in% colnames(temporal$trajectories))
})

test_that("E2E: no temporal structure — few significant trajectories", {
  f <- make_e2e_fixture(seed = 101L, n_cells = 250L)
  temporal <- ComputeTemporalPermutation(
    f$x, f$ct, f$lr_db, f$pt, f$lw,
    n_bins = f$n_bins, min_pct = 0.0,
    min_n_eff = 3.0, min_valid_bins = 2L, min_valid_fraction = 0.1,
    n_perms = 99L, seed = 3L, min_perm_observed = 50L
  )
  traj <- temporal$trajectories
  traj <- traj[is.finite(traj$p_temporal), , drop = FALSE]
  expect_gte(nrow(traj), 5L)
  frac_sig <- mean(traj$q_temporal <= 0.05, na.rm = TRUE)
  expect_lt(frac_sig, 0.15)
})

test_that("E2E: injected temporal structure is detectable", {
  skip_if_not_installed("Matrix")
  set.seed(202L)
  n_genes <- 20L
  n_cells <- 120L
  n_types <- 2L
  n_bins <- 6L
  x <- Matrix::Matrix(0, nrow = n_genes, ncol = n_cells, sparse = TRUE)
  ct <- rep(1:2, each = n_cells / 2)
  pt <- rep(seq(0, 1, length.out = n_cells / 2), 2)
  lw <- matrix(1, nrow = n_cells, ncol = 1L)
  dynamic <- pmax(0, pt - 0.3) * 5
  lig_idx <- 1L
  rec_idx <- 2L
  for (j in seq_len(n_cells)) {
    if (ct[j] == 1L) {
      x[lig_idx, j] <- dynamic[j]
    } else {
      x[rec_idx, j] <- dynamic[j]
    }
  }
  lr_db <- list(ligand_rows = list(lig_idx), receptor_rows = list(rec_idx))

  temporal <- ComputeTemporalPermutation(
    x, ct, lr_db, pt, lw,
    n_bins = n_bins, min_pct = 0.0,
    min_n_eff = 3.0, min_valid_bins = 2L, min_valid_fraction = 0.1,
    n_perms = 199L, seed = 5L
  )
  target <- temporal$trajectories
  target <- target[
    target$pair_idx == 1L & target$source == 1L &
      target$target == 2L & target$lineage == 1L,
  , drop = FALSE]
  expect_equal(nrow(target), 1L)
  expect_true(as.logical(target$trajectory_evaluable[1]))
  expect_lt(target$p_temporal[1], 0.05)
})
