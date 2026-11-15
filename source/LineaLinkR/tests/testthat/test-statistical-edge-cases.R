test_that("intersection-union p is NA unless both sides are evaluable", {
  p_lig <- c(0.01, 0.02, NA_real_, 0.05)
  p_rec <- c(0.03, NA_real_, 0.04, NA_real_)
  ok <- c(TRUE, TRUE, TRUE, FALSE)
  out <- LineaLinkR:::.intersection_union_p(p_lig, p_rec, ok)
  expect_equal(out[1], 0.03)
  expect_true(is.na(out[2]))
  expect_true(is.na(out[3]))
  expect_true(is.na(out[4]))
})

test_that("dominant role prefers lateral/bystander when main roles are zero", {
  mat <- data.frame(
    instructive_like = 0,
    transition_like = 0,
    reinforcing_like = 0,
    feedback_like = 0,
    lateral_like = 0.2,
    bystander_like = 0.8,
    stringsAsFactors = FALSE
  )
  expect_equal(LineaLinkR:::.dominant_role_like(mat), "bystander_like")
})

test_that("radius spatial adjacency excludes self on diagonal", {
  skip_if_not_installed("Matrix")
  ct <- rep(1L, 5L)
  coords <- cbind(seq(0, 4), rep(0, 5))
  adj <- ComputeSpatialAdjacency(
    ct, coords, method = "radius", radius = 0.6, binary = FALSE
  )
  expect_equal(adj[1, 1], 0)
})

test_that("adaptive temporal merge recomputes q on full family", {
  f <- make_temporal_fixture()
  res <- suppressWarnings(
    ComputeAdaptiveTemporalPermutation(
      f$x, f$ct, f$lr_db,
      pseudotime = f$pt,
      lineage_weights = f$lw,
      n_bins = 4L,
      n_perms_screen = 19L,
      n_perms_confirm = 49L,
      q_confirm_threshold = 1,
      confirm_top_k = 5L,
      min_perm_observed = 0L,
      min_perm_coverage = 0,
      seed = 1L
    )
  )
  traj <- res$trajectories
  expect_equal(traj$q_temporal, stats::p.adjust(traj$p_temporal, method = "BH"))
})

test_that("BootstrapSamples duplicates resampled blocks", {
  sample_effects <- data.frame(
    pair_idx = 1L,
    source = 1L,
    target = 2L,
    lineage = 1L,
    block_id = rep(1:2, each = 2L),
    block_label = rep(c("A", "B"), each = 2L),
    temporal_effect = c(1, 1.1, -0.5, -0.4),
    effect_se = rep(0.2, 4L),
    evaluable = TRUE,
    stringsAsFactors = FALSE
  )
  set.seed(42)
  boot <- BootstrapSamples(sample_effects, n_boot = 200L, seed = 42L)
  expect_false(is.null(boot$summary))
  expect_true(all(is.finite(boot$summary$meta_effect_q025)))
})

test_that("formal evaluable filter drops FALSE rows", {
  df <- data.frame(
    evaluable = c(TRUE, FALSE, NA),
    x = 1:3,
    stringsAsFactors = FALSE
  )
  out <- LineaLinkR:::.filter_formal_evaluable_rows(df)
  expect_equal(nrow(out), 1L)
  expect_true(out$evaluable)
})
