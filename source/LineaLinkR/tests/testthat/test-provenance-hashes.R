make_provenance_fixture <- function() {
  set.seed(11)
  n_genes <- 30L
  n_cells <- 40L
  x <- Matrix::Matrix(
    abs(rnorm(n_genes * n_cells)),
    nrow = n_genes,
    sparse = TRUE
  )
  ct <- rep(1:4, each = 10L)
  pt <- rep(seq(0, 1, length.out = 10L), 4L)
  lw <- matrix(runif(n_cells * 2L), nrow = n_cells, ncol = 2L)
  lw <- lw / rowSums(lw)
  lr_db <- list(
    ligand_rows = list(1L, 2L),
    receptor_rows = list(3L, 4L)
  )
  list(x = x, ct = ct, pt = pt, lw = lw, lr_db = lr_db)
}

test_that("ComputeTemporalPermutation exposes split provenance hashes", {
  f <- make_provenance_fixture()
  out <- ComputeTemporalPermutation(
    f$x, f$ct, f$lr_db, f$pt, f$lw,
    n_bins = 4L, min_pct = 0, n_perms = 2L, seed = 7L,
    min_perm_observed = 1L, min_perm_coverage = 0
  )
  expect_true(grepl("^llr_", out$analysis_data_hash))
  expect_true(grepl("^llr_", out$scoring_config_hash))
  expect_true(grepl("^llr_", out$temporal_null_config_hash))
  expect_equal(nchar(out$analysis_data_hash), 36L)
  expect_identical(out$calibration_config_hash, out$scoring_config_hash)
  expect_equal(
    unique(out$trajectories$analysis_data_hash),
    out$analysis_data_hash
  )
  expect_equal(
    unique(out$trajectories$scoring_config_hash),
    out$scoring_config_hash
  )
})

test_that("ComputeFakeLRNull exposes split provenance hashes", {
  f <- make_provenance_fixture()
  suppressWarnings(
    out <- ComputeFakeLRNull(
      f$x, f$ct, f$lr_db,
      candidate_fake_genes = 5:20,
      pseudotime = f$pt,
      lineage_weights = f$lw,
      n_bins = 4L,
      n_fake_per_real = 100L,
      n_strata = 2L,
      caliper_k = 2L,
      min_pct = 0,
      seed = 7L,
      min_fake_observed = 1L,
      min_fake_coverage = 0
    )
  )
  expect_true(grepl("^llr_", out$analysis_data_hash))
  expect_true(grepl("^llr_", out$scoring_config_hash))
  expect_true(grepl("^llr_", out$fake_null_config_hash))
  expect_identical(out$calibration_config_hash, out$scoring_config_hash)
  expect_equal(
    unique(out$trajectory_calibration$analysis_data_hash),
    out$analysis_data_hash
  )
})

test_that("matching pipelines share analysis and scoring hashes", {
  f <- make_provenance_fixture()
  temporal <- ComputeTemporalPermutation(
    f$x, f$ct, f$lr_db, f$pt, f$lw,
    n_bins = 4L, min_pct = 0, n_perms = 2L, seed = 7L,
    min_perm_observed = 1L, min_perm_coverage = 0
  )
  suppressWarnings(
    fake <- ComputeFakeLRNull(
      f$x, f$ct, f$lr_db,
      candidate_fake_genes = 5:20,
      pseudotime = f$pt,
      lineage_weights = f$lw,
      n_bins = 4L,
      n_fake_per_real = 100L,
      n_strata = 2L,
      caliper_k = 2L,
      min_pct = 0,
      seed = 7L,
      min_fake_observed = 1L,
      min_fake_coverage = 0
    )
  )
  expect_identical(
    temporal$analysis_data_hash,
    fake$analysis_data_hash
  )
  expect_identical(
    temporal$scoring_config_hash,
    fake$scoring_config_hash
  )
  expect_false(identical(
    temporal$temporal_null_config_hash,
    fake$fake_null_config_hash
  ))
})

test_that("ComputeRoleScores exposes role_config_hash", {
  roles <- ComputeRoleScores(
    scores_df = data.frame(
      bin = 1L,
      lineage = 1L,
      n_eff_cells = 10
    ),
    spec_self = 0.8,
    spec_max_other = 0.1,
    delta_tau = 0.2,
    bif_df = data.frame(
      lineage = 1L,
      b_star = 2L,
      reason = "crossing"
    ),
    tau_0 = 0.1,
    h_fork = 0.5,
    min_n_eff_cells = 5
  )
  expect_true(grepl("^llr_", attr(roles, "role_config_hash")))
})

test_that("analysis_data_hash changes with sparse values and positions", {
  x0 <- Matrix::sparseMatrix(
    i = c(1L, 2L), j = c(1L, 2L), x = c(1, 2), dims = c(3L, 3L)
  )
  x_value <- x0
  x_value@x[1L] <- x_value@x[1L] + 0.25
  x_position <- Matrix::sparseMatrix(
    i = c(1L, 3L), j = c(1L, 2L), x = c(1, 2), dims = c(3L, 3L)
  )
  ct <- c(1L, 1L, 2L)
  pt <- c(0, 0.5, 1)
  lw <- matrix(1, nrow = 3L, ncol = 1L)
  lr_db <- list(ligand_rows = list(1L), receptor_rows = list(2L))
  hash0 <- LineaLinkR:::.compute_analysis_data_hash(x0, ct, lr_db, pt, lw)
  hash_value <- LineaLinkR:::.compute_analysis_data_hash(
    x_value, ct, lr_db, pt, lw
  )
  hash_position <- LineaLinkR:::.compute_analysis_data_hash(
    x_position, ct, lr_db, pt, lw
  )
  expect_false(identical(hash0, hash_value))
  expect_false(identical(hash0, hash_position))
})

test_that("role_config_hash includes every lineage tau", {
  bif <- data.frame(
    lineage = c(1L, 2L), b_star = c(3L, 4L), reason = "crossing"
  )
  hash_a <- LineaLinkR:::.compute_role_config_hash(
    c(`1` = 0.1, `2` = 0.2), 0.5, 5, timing_source = "cell_state_tau",
    bif_df = bif
  )
  hash_b <- LineaLinkR:::.compute_role_config_hash(
    c(`1` = 0.3, `2` = 0.2), 0.5, 5, timing_source = "cell_state_tau",
    bif_df = bif
  )
  expect_false(identical(hash_a, hash_b))
})
