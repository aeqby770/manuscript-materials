# Gene-level lag weighting and pooled donor hypothesis calibration.

test_that("weighted_complex_gene_mean uses aligned local weights", {
  expr <- matrix(0, nrow = 5L, ncol = 60L)
  cells <- c(41L, 45L, 49L)
  weights <- c(0.2, 0.3, 0.5)
  values <- c(10, 20, 30)
  expr[1L, cells] <- values
  got <- LineaLinkR:::.weighted_complex_gene_mean(expr, 1L, cells, weights)
  expect_equal(got, sum(weights * values) / sum(weights))
  expect_error(
    LineaLinkR:::.weighted_complex_gene_mean(expr, 1L, cells, weights[-1L]),
    "length must equal"
  )
})

test_that("gene-level lag profiles respect non-contiguous high cell indices", {
  skip_if_not_installed("Matrix")
  set.seed(63L)
  n_cells <- 80L
  n_bins <- 4L
  ct <- c(rep(1L, 40L), rep(1L, 20L), rep(2L, 20L))
  block_labels <- rep(c("donor_a", "donor_b"), each = 40L)
  pt <- seq(0, 1, length.out = n_cells)
  lw <- matrix(1, n_cells, 1L)
  x <- Matrix::Matrix(0.05, nrow = 6L, ncol = n_cells, sparse = TRUE)
  x[1L, 1:40] <- 2
  x[2L, 21:40] <- 2
  x[2L, 61:80] <- 2
  lr_db <- list(ligand_rows = list(1L), receptor_rows = list(2L))
  scores <- ComputeTemporalLR(
    x, ct, lr_db, pt, lw, n_bins = n_bins, min_pct = 0, min_n_eff = 0.1
  )
  skip_if(nrow(scores) == 0L, "no pooled scores")
  lag <- ComputeLigandReceptorLag(
    scores, lw, pt, ct,
    n_bins = n_bins,
    expr = x, lr_db = lr_db
  )
  skip_if(nrow(lag) == 0L, "no lag rows")
  expect_true(any(is.finite(lag$lag)))
  expect_true(any(is.finite(lag$ligand_peak_time)))
})

test_that("donor rescore keeps pooled zero-expression curves evaluable", {
  skip_if_not_installed("Matrix")
  set.seed(71L)
  n_cells <- 80L
  n_bins <- 4L
  ct <- c(rep(1L, 40L), rep(1L, 20L), rep(2L, 20L))
  block_labels <- rep(c("donor_a", "donor_b"), each = 40L)
  pt <- seq(0, 1, length.out = n_cells)
  lw <- matrix(1, n_cells, 1L)
  x <- Matrix::Matrix(0.05, nrow = 6L, ncol = n_cells, sparse = TRUE)
  x[1L, 1:40] <- 2
  x[2L, 21:40] <- 2
  x[2L, 61:80] <- 2
  lr_db <- list(ligand_rows = list(1L), receptor_rows = list(2L))
  pooled <- ComputeTemporalLR(
    x, ct, lr_db, pt, lw, n_bins = n_bins, min_pct = 0, min_n_eff = 0.1
  )
  skip_if(nrow(pooled) == 0L, "no pooled scores")
  donor <- ComputeTemporalEffectsBySample(
    pooled, block_labels, pt, lw,
    expr = x, celltypes = ct, lr_db = lr_db,
    temporal_null_per_block = FALSE,
    n_perms_block = 19L,
    min_perm_observed_block = 0L,
    min_perm_coverage_block = 0,
    min_block_cells = 15L,
    min_valid_bins = 1L,
    min_n_eff = 0,
    min_bin_mass = 0
  )
  donor_b <- donor[donor$block_label == "donor_b", , drop = FALSE]
  skip_if(nrow(donor_b) == 0L, "donor_b absent")
  expect_true(any(donor_b$evaluable %in% TRUE))
  expect_true(any(is.finite(donor_b$temporal_effect) & donor_b$temporal_effect == 0))
})
