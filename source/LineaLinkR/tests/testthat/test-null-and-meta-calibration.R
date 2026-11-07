# Null calibration, spatial block permutations, role block bootstrap, and lag meta-analysis.

test_that("ComputeFakeLRNull exposes fake_null_median", {
  skip_if_not_installed("Matrix")
  set.seed(5)
  n_cells <- 70L
  x <- Matrix::rsparsematrix(18L, n_cells, density = 0.35)
  x@x <- abs(x@x) + 0.05
  ct <- sample.int(2L, n_cells, replace = TRUE)
  pt <- seq(0, 1, length.out = n_cells)
  lw <- matrix(0.9, n_cells, 1L)
  lr_db <- list(ligand_rows = list(1L), receptor_rows = list(2L))
  fake_genes <- 3:10
  out <- ComputeFakeLRNull(
    x, ct, lr_db, pt, lw,
    candidate_fake_genes = fake_genes,
    n_bins = 4L, n_fake_per_real = 100L,
    min_pct = 0, min_n_eff = 0.1, min_valid_bins = 2L,
    min_valid_fraction = 0, min_fake_observed = 0L,
    seed = 9L
  )
  cal <- out$trajectory_calibration
  expect_true("fake_null_median" %in% colnames(cal))
  expect_true(any(is.finite(cal$fake_null_median)))
})

test_that("donor identity null_centered uses fake_null_median when available", {
  skip_if_not_installed("Matrix")
  set.seed(13L)
  n_cells <- 90L
  x <- Matrix::rsparsematrix(20L, n_cells, density = 0.3)
  x@x <- abs(x@x) + 0.02
  ct <- sample.int(2L, n_cells, replace = TRUE)
  pt <- seq(0, 1, length.out = n_cells)
  lw <- matrix(0.9, n_cells, 1L)
  block_labels <- rep(c("A", "B", "C"), each = 30L)
  lr_db <- list(ligand_rows = list(1L), receptor_rows = list(2L))
  fake_genes <- 3:12
  scores <- ComputeTemporalLR(x, ct, lr_db, pt, lw, n_bins = 4L, min_pct = 0)
  skip_if(nrow(scores) == 0L, "no scores")
  out <- ComputeTemporalEffectsBySample(
    scores, block_labels, pt, lw,
    expr = x, celltypes = ct, lr_db = lr_db,
    candidate_fake_genes = fake_genes,
    fake_null_per_block = TRUE,
    n_fake_per_real_block = 100L,
    min_fake_observed_block = 1L,
    min_fake_coverage_block = 0,
    min_unique_fake_draws_block = 1L,
    min_unique_fake_fraction_block = 0,
    min_pct_block = 0,
    n_perms_block = 19L, min_perm_observed_block = 0L, min_perm_coverage_block = 0,
    min_block_cells = 10L,
    min_valid_bins = 1L, min_n_eff = 0, min_bin_mass = 0
  )
  expect_true("fake_null_median" %in% colnames(out))
  expect_true(any(is.finite(out$identity_effect_null_centered)))
})

test_that("spatial block_label_permutation null runs", {
  set.seed(2)
  n_cells <- 60L
  ct <- sample.int(3L, n_cells, replace = TRUE)
  coords <- cbind(runif(n_cells), runif(n_cells))
  blocks <- rep(c("D1", "D2"), each = 30L)
  null <- ComputeSpatialAdjacencyNull(
    ct, coords, method = "knn", k = 4L,
    null_mode = "block_label_permutation",
    block_labels = blocks, n_perm = 29L, seed = 1L
  )
  expect_equal(null$null_mode, "block_label_permutation")
  expect_true(all(dim(null$null_mean) == c(3L, 3L)))
})

test_that("BootstrapRoleInference block mode requires rebuild_context", {
  rf <- data.frame(
    pair_idx = 1L, source = 1L, target = 2L, lineage = 1L, bin = 1:4,
    spec_self = 0.5, spec_max_other = 0.2, delta_tau = 0.1,
    n_eff_cells = 10, stringsAsFactors = FALSE
  )
  expect_error(
    BootstrapRoleInference(
      rf, n_boot = 5L, bootstrap_mode = "block",
      block_labels = rep(1:2, each = 2L)
    ),
    "rebuild_context"
  )
})

test_that("SummarizeLagByDonor passes compute_lag_null", {
  skip_if_not_installed("Matrix")
  set.seed(17)
  n_cells <- 80L
  x <- Matrix::rsparsematrix(15L, n_cells, density = 0.3)
  x@x <- abs(x@x)
  ct <- sample.int(2L, n_cells, replace = TRUE)
  pt <- seq(0, 1, length.out = n_cells)
  lw <- matrix(1, n_cells, 1L)
  blocks <- rep(c("A", "B"), each = 40L)
  lr_db <- list(ligand_rows = list(1L), receptor_rows = list(2L))
  scores <- ComputeTemporalLR(x, ct, lr_db, pt, lw, n_bins = 5L, min_pct = 0)
  skip_if(nrow(scores) == 0L, "no scores")
  lag_by_donor <- SummarizeLagByDonor(
    scores, blocks, lw, pt, ct,
    n_bins = 5L, min_block_cells = 15L,
    compute_lag_null = TRUE, n_lag_perm = 39L
  )
  expect_true("p_lag" %in% colnames(lag_by_donor))
})

test_that("MetaAnalyzeLagEffects defaults to stouffer on p_lag", {
  df <- data.frame(
    pair_idx = rep(1L, 3L), source = 1L, target = 2L, lineage = 1L,
    block_id = 1:3, block_label = c("A", "B", "C"),
    lag = c(0.2, 0.3, 0.1),
    lag_null_centered = c(0.4, 0.5, 0.2),
    p_lag = c(0.04, 0.02, 0.06),
    stringsAsFactors = FALSE
  )
  meta <- MetaAnalyzeLagEffects(df, include_lodo = FALSE)
  expect_equal(meta$meta_method, "stouffer")
  expect_true(is.finite(meta$meta_p))
  expect_true(is.finite(meta$q_meta))
})
