test_that("cardinality_matched complex FakeLR yields p_fake_complex", {
  skip_if_not_installed("Matrix")
  set.seed(7)
  n_genes <- 24
  n_cells <- 48
  x <- Matrix::rsparsematrix(n_genes, n_cells, density = 0.35)
  x@x <- abs(x@x) + 0.1
  ct <- sample.int(3L, n_cells, replace = TRUE)
  pt <- seq(0, 1, length.out = n_cells)
  lw <- matrix(c(0.7, 0.3), nrow = n_cells, ncol = 2, byrow = TRUE)
  lr_db <- list(
    ligand_rows = list(c(1L, 3L), 5L),
    receptor_rows = list(c(2L, 4L), 6L)
  )
  fakes <- 13:22

  skip_res <- ComputeFakeLRNull(
    x, ct, lr_db,
    candidate_fake_genes = fakes,
    pseudotime = pt,
    lineage_weights = lw,
    n_bins = 4L,
    n_fake_per_real = 100L,
    n_strata = 3L,
    caliper_k = 8L,
    max_match_distance = 2,
    min_pct = 0,
    min_n_eff = 0.1,
    min_valid_bins = 2L,
    min_valid_fraction = 0.1,
    min_fake_observed = 0L,
    min_fake_coverage = 0,
    min_fake_observed_generic = 0L,
    min_fake_coverage_generic = 0,
    min_unique_fake_draws = 0L,
    min_unique_fake_fraction = 0,
    seed = 11L,
    complex_null_mode = "skip"
  )
  complex_rows <- skip_res$trajectory_calibration$complex_pair
  expect_true(any(complex_rows))
  expect_true(all(is.na(skip_res$trajectory_calibration$p_fake[complex_rows])))

  match_res <- ComputeFakeLRNull(
    x, ct, lr_db,
    candidate_fake_genes = fakes,
    pseudotime = pt,
    lineage_weights = lw,
    n_bins = 4L,
    n_fake_per_real = 100L,
    n_strata = 3L,
    caliper_k = 8L,
    max_match_distance = 2,
    min_pct = 0,
    min_n_eff = 0.1,
    min_valid_bins = 2L,
    min_valid_fraction = 0.1,
    min_fake_observed = 0L,
    min_fake_coverage = 0,
    min_fake_observed_generic = 0L,
    min_fake_coverage_generic = 0,
    min_unique_fake_draws = 0L,
    min_unique_fake_fraction = 0,
    seed = 11L,
    complex_null_mode = "cardinality_matched",
    fake_null_mode = "context_generic"
  )
  cal <- match_res$trajectory_calibration
  expect_true("p_fake_complex" %in% names(cal))
  expect_true("complex_null_mode" %in% names(cal))
  complex_cal <- cal[cal$complex_pair, , drop = FALSE]
  expect_gt(nrow(complex_cal), 0L)
  expect_true(any(is.finite(complex_cal$p_fake_complex)))
  expect_true(any(complex_cal$identity_calibration_evaluable))
  simple_cal <- cal[!cal$complex_pair, , drop = FALSE]
  if (nrow(simple_cal) > 0L) {
    expect_true(any(is.finite(simple_cal$p_fake_generic) | is.finite(simple_cal$p_fake)))
  }
})

test_that("cardinality_matched supports up to 8 subunits per side", {
  skip_if_not_installed("Matrix")
  set.seed(19)
  n_genes <- 30
  n_cells <- 40
  x <- Matrix::rsparsematrix(n_genes, n_cells, density = 0.3)
  x@x <- abs(x@x) + 0.1
  ct <- sample.int(2L, n_cells, replace = TRUE)
  pt <- seq(0, 1, length.out = n_cells)
  lw <- matrix(c(0.6, 0.4), nrow = n_cells, ncol = 2, byrow = TRUE)
  lr_db <- list(
    ligand_rows = list(c(1L, 3L, 5L)),
    receptor_rows = list(c(2L, 4L, 6L, 7L))
  )
  res <- ComputeFakeLRNull(
    x, ct, lr_db,
    candidate_fake_genes = 10:25,
    pseudotime = pt,
    lineage_weights = lw,
    n_bins = 4L,
    n_fake_per_real = 100L,
    n_strata = 2L,
    caliper_k = 6L,
    max_match_distance = 3,
    min_pct = 0,
    min_n_eff = 0.1,
    min_valid_bins = 2L,
    min_valid_fraction = 0.1,
    min_fake_observed = 0L,
    min_fake_coverage = 0,
    min_fake_observed_generic = 0L,
    min_fake_coverage_generic = 0,
    min_unique_fake_draws = 0L,
    min_unique_fake_fraction = 0,
    seed = 3L,
    complex_null_mode = "cardinality_matched",
    fake_null_mode = "context_generic"
  )
  cal <- res$trajectory_calibration
  expect_true(all(cal$complex_pair))
  expect_true(any(is.finite(cal$p_fake_complex)))
})

test_that("cardinality_matched complex supports edge and all null modes", {
  skip_if_not_installed("Matrix")
  set.seed(21)
  n_genes <- 28
  n_cells <- 56
  x <- Matrix::rsparsematrix(n_genes, n_cells, density = 0.35)
  x@x <- abs(x@x) + 0.1
  ct <- sample.int(3L, n_cells, replace = TRUE)
  pt <- seq(0, 1, length.out = n_cells)
  lw <- matrix(c(0.7, 0.3), nrow = n_cells, ncol = 2, byrow = TRUE)
  lr_db <- list(
    ligand_rows = list(c(1L, 3L)),
    receptor_rows = list(c(2L, 4L))
  )
  both_res <- ComputeFakeLRNull(
    x, ct, lr_db,
    candidate_fake_genes = 12:24,
    pseudotime = pt,
    lineage_weights = lw,
    n_bins = 4L,
    n_fake_per_real = 100L,
    n_strata = 3L,
    caliper_k = 8L,
    max_match_distance = 2,
    min_pct = 0,
    min_n_eff = 0.1,
    min_valid_bins = 2L,
    min_valid_fraction = 0.1,
    min_fake_observed = 0L,
    min_fake_coverage = 0,
    min_fake_observed_generic = 0L,
    min_fake_coverage_generic = 0,
    min_unique_fake_draws = 0L,
    min_unique_fake_fraction = 0,
    seed = 9L,
    complex_null_mode = "cardinality_matched",
    fake_null_mode = "both"
  )
  cal <- both_res$trajectory_calibration
  expect_true(all(cal$complex_pair))
  expect_true(any(is.finite(cal$p_fake_complex)))
  expect_true(any(is.finite(cal$p_edge_iu) | is.finite(cal$p_fake_edge)))

  all_res <- ComputeFakeLRNull(
    x, ct, lr_db,
    candidate_fake_genes = 12:24,
    pseudotime = pt,
    lineage_weights = lw,
    n_bins = 4L,
    n_fake_per_real = 100L,
    n_strata = 3L,
    caliper_k = 8L,
    max_match_distance = 2,
    min_pct = 0,
    min_n_eff = 0.1,
    min_valid_bins = 2L,
    min_valid_fraction = 0.1,
    min_fake_observed = 0L,
    min_fake_coverage = 0,
    min_fake_observed_generic = 0L,
    min_fake_coverage_generic = 0,
    min_unique_fake_draws = 0L,
    min_unique_fake_fraction = 0,
    seed = 9L,
    complex_null_mode = "cardinality_matched",
    fake_null_mode = "all"
  )
  cal_all <- all_res$trajectory_calibration
  expect_true(any(is.finite(cal_all$p_network_iu) | is.finite(cal_all$p_fake_network)))
})

test_that("complex_null_mode skip leaves complex identity NA without cardinality warning", {
  skip_if_not_installed("Matrix")
  set.seed(3)
  n_genes <- 20
  n_cells <- 40
  x <- Matrix::rsparsematrix(n_genes, n_cells, density = 0.3)
  x@x <- abs(x@x) + 0.1
  ct <- sample.int(2L, n_cells, replace = TRUE)
  pt <- seq(0, 1, length.out = n_cells)
  lw <- matrix(c(0.6, 0.4), nrow = n_cells, ncol = 2, byrow = TRUE)
  lr_db <- list(
    ligand_rows = list(c(1L, 3L)),
    receptor_rows = list(2L)
  )
  res <- suppressWarnings(
    ComputeFakeLRNull(
      x, ct, lr_db,
      candidate_fake_genes = 10:18,
      pseudotime = pt,
      lineage_weights = lw,
      n_bins = 4L,
      n_fake_per_real = 100L,
      min_pct = 0,
      min_n_eff = 0.1,
      min_valid_bins = 2L,
      min_valid_fraction = 0.1,
      min_fake_observed = 0L,
      min_fake_coverage = 0,
      seed = 1L,
      complex_null_mode = "skip"
    )
  )
  complex_rows <- res$trajectory_calibration$complex_pair
  expect_true(any(complex_rows))
  expect_true(all(is.na(res$trajectory_calibration$p_fake[complex_rows])))
})
