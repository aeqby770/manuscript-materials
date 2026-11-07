# Formal API contracts for null calibration and adaptive methods.

test_that("ComputeFakeLRNull rejects n_fake_per_real < 100", {
  skip_if_not_installed("Matrix")
  set.seed(1)
  n_cells <- 50L
  x <- Matrix::rsparsematrix(15L, n_cells, density = 0.3)
  x@x <- abs(x@x) + 0.05
  ct <- rep(1:2, each = 25L)
  pt <- seq(0, 1, length.out = n_cells)
  lw <- matrix(1, n_cells, 1L)
  lr_db <- list(ligand_rows = list(1L), receptor_rows = list(2L))
  expect_error(
    ComputeFakeLRNull(
      x, ct, lr_db, pt, lw,
      candidate_fake_genes = 3:8,
      n_fake_per_real = 10L,
      n_bins = 4L, min_pct = 0, min_n_eff = 0.1,
      min_valid_bins = 2L, min_valid_fraction = 0,
      seed = 1L
    ),
    "n_fake_per_real >= 100"
  )
})

test_that("ComputeTrajectoryFakeLRNull rejects n_fake_per_real < 100", {
  skip_if_not_installed("Matrix")
  set.seed(2)
  n_cells <- 50L
  x <- Matrix::rsparsematrix(15L, n_cells, density = 0.3)
  x@x <- abs(x@x) + 0.05
  ct <- rep(1:2, each = 25L)
  pt <- seq(0, 1, length.out = n_cells)
  lw <- matrix(1, n_cells, 1L)
  lr_db <- list(ligand_rows = list(1L), receptor_rows = list(2L))
  expect_error(
    ComputeTrajectoryFakeLRNull(
      x, ct, lr_db, 3:8, pt, lw,
      n_fake_per_real = 50L,
      n_bins = 4L, min_pct = 0, min_n_eff = 0.1,
      min_valid_bins = 2L, min_valid_fraction = 0,
      min_fake_observed = 0L, min_fake_coverage = 0,
      seed = 1L
    ),
    "n_fake_per_real >= 100"
  )
})

test_that("adaptive APIs run without experimental warning", {
  f <- local({
    skip_if_not_installed("Matrix")
    set.seed(3)
    n_cells <- 120L
    x <- Matrix::rsparsematrix(25L, n_cells, density = 0.3)
    x@x <- abs(x@x)
    ct <- sample.int(2L, n_cells, replace = TRUE)
    pt <- seq(0, 1, length.out = n_cells)
    lw <- matrix(1, n_cells, 1L)
    lr_db <- list(ligand_rows = list(1L), receptor_rows = list(2L))
    list(x = x, ct = ct, pt = pt, lw = lw, lr_db = lr_db)
  })
  expect_no_warning(
    ComputeAdaptiveTemporalPermutation(
      f$x, f$ct, f$lr_db, f$pt, f$lw,
      n_bins = 4L, min_pct = 0, min_n_eff = 0.1,
      min_valid_bins = 2L, min_valid_fraction = 0,
      n_perms_screen = 19L, n_perms_confirm = 29L,
      confirm_top_k = 2L, q_confirm_threshold = 1,
      min_perm_observed = 1L, min_perm_coverage = 0, seed = 4L
    )
  )
})
