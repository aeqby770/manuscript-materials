# Ligand-target timing lag model.

make_lag_fixture <- function() {
  skip_if_not_installed("Matrix")
  set.seed(42L)
  n_genes <- 20L
  n_cells <- 80L
  n_types <- 2L
  n_bins <- 5L
  x <- Matrix::rsparsematrix(n_genes, n_cells, density = 0.35)
  x@x <- abs(x@x)
  ct <- rep(1:2, each = n_cells / 2)
  pt <- seq(0, 1, length.out = n_cells)
  lw <- matrix(c(rep(c(0.85, 0.15), n_cells / 2)), nrow = n_cells, ncol = 1L)
  lr_db <- list(
    ligand_rows = list(1L, 3L),
    receptor_rows = list(2L, 4L)
  )
  scores <- ComputeTemporalLR(
    x, ct, lr_db, pt, lw, n_bins = n_bins, min_pct = 0.0
  )
  skip_if(nrow(scores) == 0L, "no scores in fixture")
  block_labels <- rep(c("donor_A", "donor_B"), each = n_cells / 2)
  list(
    x = x, ct = ct, pt = pt, lw = lw, lr_db = lr_db,
    scores = scores, n_bins = n_bins, block_labels = block_labels
  )
}

make_side_profile_scores <- function(n_bins = 5L) {
  bins <- rep(seq_len(n_bins), each = 2L)
  lig <- c(0.1, 0.3, 1.0, 0.4, 0.2, 0.1, 0.2, 0.5, 0.9, 0.3)
  rec <- c(0.2, 0.4, 0.5, 0.8, 1.0, 0.1, 0.3, 0.4, 0.6, 1.0)
  data.frame(
    pair_idx = rep(1L, length(bins)),
    source = rep(1L, length(bins)),
    target = rep(2L, length(bins)),
    lineage = rep(1L, length(bins)),
    bin = bins,
    score = pmax(lig, rec) * 0.5,
    ligand_score = lig,
    receptor_score = rec,
    stringsAsFactors = FALSE
  )
}

test_that("ComputeLigandTargetLag peak method returns trajectory lag columns", {
  f <- make_lag_fixture()
  lag <- ComputeLigandTargetLag(
    f$scores, f$lw, f$pt, f$ct, method = "peak"
  )
  expect_true(is.data.frame(lag))
  expect_equal(nrow(lag), length(unique(do.call(
    paste, c(f$scores[, c("pair_idx", "source", "target", "lineage")], sep = "\001")
  ))))
  expect_true(all(c(
    "ligand_peak_time", "receptor_peak_time", "lag_peak", "lag",
    "lag_se", "lag_significant", "lag_threshold", "bin_width", "lag_method"
  ) %in% colnames(lag)))
  expect_equal(lag$lag, lag$lag_peak)
  expect_equal(lag$lag_peak, lag$receptor_peak_time - lag$ligand_peak_time,
               tolerance = 1e-8, ignore_attr = TRUE)
  expect_true(all(is.finite(lag$lag_se) | is.na(lag$lag)))
})

test_that("ComputeLigandTargetLag respects explicit side profiles", {
  f <- make_lag_fixture()
  side <- make_side_profile_scores(f$n_bins)
  pt <- seq(0, 1, length.out = 40L)
  ct <- rep(1:2, each = 20L)
  lw <- matrix(1, nrow = 40L, ncol = 1L)
  lag <- ComputeLigandTargetLag(
    side, lw, pt, ct, method = "peak", n_bins = f$n_bins
  )
  expect_equal(nrow(lag), 1L)
  expect_true(is.finite(lag$lag_peak[1L]))
  expect_equal(lag$lag_method[1L], "peak")
})

test_that("ComputeLigandTargetLag centroid and cross_cor methods run", {
  f <- make_lag_fixture()
  side <- make_side_profile_scores(f$n_bins)
  pt <- seq(0, 1, length.out = 40L)
  ct <- rep(1:2, each = 20L)
  lw <- matrix(1, nrow = 40L, ncol = 1L)
  lag_c <- ComputeLigandTargetLag(
    side, lw, pt, ct, method = "centroid", n_bins = f$n_bins
  )
  lag_x <- ComputeLigandTargetLag(
    side, lw, pt, ct, method = "cross_cor", n_bins = f$n_bins
  )
  expect_true(is.finite(lag_c$lag_centroid[1L]))
  expect_equal(lag_c$lag, lag_c$lag_centroid[1L])
  expect_true(is.finite(lag_x$lag_cross_cor[1L]))
  expect_true(is.finite(lag_x$cross_cor_max[1L]))
})

test_that("cross correlation preserves original missing-bin positions", {
  timing <- seq(0, 1, length.out = 6L)
  est <- LineaLinkR:::.lag_from_cross_cor(
    c(1, 3, 2, NA, NA, NA),
    c(NA, NA, 1, 3, 2, NA),
    timing, timing, bin_width = 1, max_lag_bins = 3L
  )
  expect_equal(est$lag_cross_cor, 2)
  expect_equal(est$n_valid_bins, 3L)
})

test_that("cross-correlation ties do not prefer negative lag", {
  timing <- seq(0, 1, length.out = 5L)
  est <- LineaLinkR:::.lag_from_cross_cor(
    c(0, 0, 0, 1, 0),
    c(1, 0, 1, 0, 1),
    timing, timing, bin_width = 1, max_lag_bins = 2L
  )

  expect_equal(est$lag_cross_cor, 0)
  expect_true(est$lag_tie)
  expect_equal(est$n_equivalent_lags, 2L)
})

test_that("small circular-shift null enumerates unique states", {
  timing <- seq(0, 1, length.out = 5L)
  null <- LineaLinkR:::.lag_permutation_null(
    c(1, 2, 3, 2, 1), c(2, 1, 2, 3, 2), timing, timing,
    method = "peak", bin_width = 0.25, max_lag_bins = 2L,
    observed_lag = 0.25, n_perm = 199L, seed = 3L
  )

  expect_equal(null$n_lag_perm_observed, 4L)
})

test_that("ligand-receptor lag name is the canonical API", {
  f <- make_lag_fixture()
  side <- make_side_profile_scores(f$n_bins)
  pt <- seq(0, 1, length.out = 40L)
  ct <- rep(1:2, each = 20L)
  lw <- matrix(1, nrow = 40L, ncol = 1L)
  canonical <- ComputeLigandReceptorLag(
    side, lw, pt, ct, method = "peak", n_bins = f$n_bins
  )
  legacy <- ComputeLigandTargetLag(
    side, lw, pt, ct, method = "peak", n_bins = f$n_bins
  )
  expect_equal(canonical, legacy)
})

test_that("lag_significant uses default threshold of one bin width", {
  side <- make_side_profile_scores(5L)
  pt <- seq(0, 1, length.out = 40L)
  ct <- rep(1:2, each = 20L)
  lw <- matrix(1, nrow = 40L, ncol = 1L)
  lag <- ComputeLigandTargetLag(
    side, lw, pt, ct, method = "peak", n_bins = 5L
  )
  expect_equal(lag$lag_threshold[1L], lag$bin_width[1L])
  if (is.finite(lag$lag[1L])) {
    expect_equal(lag$lag_significant[1L], abs(lag$lag[1L]) > lag$lag_threshold[1L])
  }
})

test_that("MergeLagIntoRoleFeatures and compute_lag integrate with role features", {
  f <- make_lag_fixture()
  res <- ComputeTemporalLRWithIntermediate(
    f$x, f$ct, f$lr_db, f$pt, f$lw, n_bins = f$n_bins, min_pct = 0.0
  )
  scores <- res$scores
  skip_if(nrow(scores) == 0L, "no scores")
  feats <- ComputeRoleFeatures(
    scores, f$lw, f$pt, f$ct,
    intermediate = res$intermediate,
    compute_lag = TRUE,
    lag_method = "peak",
    expr = f$x,
    lr_db = f$lr_db
  )
  expect_true(all(c("lag", "lag_se", "lag_significant", "lag_peak") %in% colnames(feats)))
  expect_equal(nrow(feats), nrow(scores))
  lag_df <- attr(feats, "lag_df")
  expect_true(is.data.frame(lag_df))
  merged <- MergeLagIntoRoleFeatures(
    feats[, setdiff(colnames(feats), c("lag", "lag_se", "lag_significant"))],
    lag_df
  )
  expect_equal(merged$lag, feats$lag)
})

test_that("SummarizeLagByDonor returns per-donor lag rows", {
  f <- make_lag_fixture()
  donor_lag <- SummarizeLagByDonor(
    f$scores,
    block_labels = f$block_labels,
    lineage_weights = f$lw,
    pseudotime = f$pt,
    celltypes = f$ct,
    method = "peak",
    min_block_cells = 10L,
    expr = f$x,
    lr_db = f$lr_db
  )
  expect_true(all(c("block_label", "block_id", "n_cells", "lag", "lag_se") %in%
                    colnames(donor_lag)))
  expect_gt(length(unique(donor_lag$block_label)), 1L)
  expect_true(any(is.finite(donor_lag$lag)))
})

test_that("SummarizeLagByDonor requires block_labels", {
  f <- make_lag_fixture()
  expect_error(
    SummarizeLagByDonor(
      f$scores,
      block_labels = NULL,
      lineage_weights = f$lw,
      pseudotime = f$pt,
      celltypes = f$ct
    ),
    "block_labels"
  )
})

test_that("ComputeLigandTargetLag reports lag_profile_source", {
  f <- make_lag_fixture()
  lag_expr <- ComputeLigandTargetLag(
    f$scores, f$lw, f$pt, f$ct, method = "peak",
    expr = f$x, lr_db = f$lr_db
  )
  expect_true("lag_profile_source" %in% colnames(lag_expr))
  expect_equal(unique(lag_expr$lag_profile_source), "expr_lr_db")

  side <- make_side_profile_scores(f$n_bins)
  lag_side <- ComputeLigandTargetLag(
    side, matrix(1, 40L, 1L), seq(0, 1, length.out = 40L), rep(1:2, each = 20L),
    method = "peak", n_bins = f$n_bins
  )
  expect_equal(unique(lag_side$lag_profile_source), "explicit_columns")
})

test_that("ComputeRoleFeatures auto-enriches lag from expr and lr_db", {
  f <- make_lag_fixture()
  inter <- ComputeTemporalLRWithIntermediate(
    f$x, f$ct, f$lr_db, f$pt, f$lw, n_bins = f$n_bins, min_pct = 0.0
  )
  feats <- ComputeRoleFeatures(
    inter$scores, f$lw, f$pt, f$ct,
    intermediate = inter$intermediate,
    compute_lag = TRUE,
    expr = f$x,
    lr_db = f$lr_db
  )
  expect_true(all(c("lag", "lag_profile_source") %in% colnames(feats)))
  expect_equal(unique(feats$lag_profile_source), "expr_lr_db")
})
