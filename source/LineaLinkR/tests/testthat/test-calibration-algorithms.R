# Calibration algorithms: permutation null medians, lag nulls, and identity meta-analysis.

test_that("ComputeTemporalPermutation exposes temporal_null_median", {
  skip_if_not_installed("Matrix")
  set.seed(7)
  n_cells <- 80L
  x <- Matrix::rsparsematrix(15L, n_cells, density = 0.3)
  x@x <- abs(x@x) + 0.05
  ct <- sample.int(2L, n_cells, replace = TRUE)
  pt <- seq(0, 1, length.out = n_cells)
  lw <- matrix(0.85, n_cells, 1L)
  lr_db <- list(ligand_rows = list(1L), receptor_rows = list(2L))
  res <- ComputeTemporalPermutation(
    x, ct, lr_db, pt, lw,
    n_bins = 4L, n_perms = 29L, min_perm_observed = 0L,
    min_perm_coverage = 0, min_pct = 0, seed = 3L
  )
  expect_true("temporal_null_median" %in% colnames(res$trajectories))
  expect_true(any(is.finite(res$trajectories$temporal_null_median)))
})

test_that("formal donor-rescore auto-enables temporal_null_per_block", {
  skip_if_not_installed("Matrix")
  set.seed(11L)
  n_cells <- 90L
  x <- Matrix::rsparsematrix(20L, n_cells, density = 0.3)
  x@x <- abs(x@x)
  ct <- sample.int(2L, n_cells, replace = TRUE)
  pt <- seq(0, 1, length.out = n_cells)
  lw <- matrix(0.9, n_cells, 1L)
  block_labels <- rep(c("A", "B", "C"), each = 30L)
  lr_db <- list(ligand_rows = list(1L), receptor_rows = list(2L))
  scores <- ComputeTemporalLR(x, ct, lr_db, pt, lw, n_bins = 4L, min_pct = 0)
  skip_if(nrow(scores) == 0L, "no scores")
  out <- ComputeTemporalEffectsBySample(
    scores, block_labels, pt, lw,
    expr = x, celltypes = ct, lr_db = lr_db,
    n_perms_block = 19L, min_block_cells = 10L,
    min_valid_bins = 1L, min_n_eff = 0, min_bin_mass = 0
  )
  expect_equal(unique(out$effect_mode), "donor_rescore")
  expect_true(any(out$temporal_null_source == "permutation", na.rm = TRUE))
  expect_true(any(is.finite(out$temporal_effect_null_centered)))
})

test_that("ComputeLigandTargetLag profile permutation yields p_lag", {
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
  expect_true(all(c("p_lag", "q_lag", "lag_null_mode") %in% colnames(lag)))
  expect_true(is.finite(lag$p_lag[1L]))
  expect_equal(lag$lag_null_mode[1L], "circular_shift")
  expect_true(is.finite(lag$lag_null_centered[1L]))
})

test_that("MetaAnalyzeIdentityEffects defaults to stouffer on p_identity", {
  df <- data.frame(
    pair_idx = rep(1L, 3L), source = 1L, target = 2L, lineage = 1L,
    block_id = 1:3, block_label = c("A", "B", "C"),
    identity_effect = c(1, 1.2, 0.8),
    identity_effect_null_centered = c(1.5, 1.8, 1.1),
    effect_se = 0.1,
    p_identity = c(0.04, 0.03, 0.05),
    identity_source = "block_fake_null",
    evaluable_identity = TRUE,
    evaluable = TRUE, stringsAsFactors = FALSE
  )
  meta <- MetaAnalyzeIdentityEffects(df, include_lodo = FALSE)
  expect_equal(meta$meta_method, "stouffer")
  expect_true(is.finite(meta$meta_p))
  expect_true(is.finite(meta$q_meta))
})
