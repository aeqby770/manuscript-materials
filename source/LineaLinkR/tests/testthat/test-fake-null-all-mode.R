test_that("all mode conserves fake budget across three null layers", {
  f <- make_temporal_fixture()
  n_fake <- 100L
  res <- NULL
  res <- ComputeFakeLRNull(
      f$x, f$ct, f$lr_db,
      candidate_fake_genes = 5:20,
      pseudotime = f$pt,
      lineage_weights = f$lw,
      n_bins = 4L,
      n_fake_per_real = n_fake,
      n_strata = 1L,
      caliper_k = 10L,
      max_match_distance = 2,
      fake_null_mode = "all",
      min_pct = 0,
      seed = 44L,
      min_fake_observed_generic = 0L,
      min_fake_coverage_generic = 0,
      min_fake_observed_edge = 0L,
      min_fake_coverage_edge = 0,
      min_unique_fake_draws = 0L,
      min_unique_fake_fraction = 0
  )
  cal <- res$trajectory_calibration
  expect_true(all(
    cal$n_fake_expected_generic + cal$n_fake_expected_edge + cal$n_fake_expected_network ==
      n_fake
  ))
  expect_true(all(cal$n_fake_expected_generic == 33L))
  expect_true(all(cal$n_fake_expected_edge == 33L))
  expect_true(all(cal$n_fake_expected_network == 34L))
  expect_true(all(cal$fake_null_mode == "all"))
  kinds <- unique(res$fake_pairs$null_kind)
  expect_true(any(kinds %in% 0L))
  expect_true(any(kinds %in% c(1L, 2L)))
  expect_true(any(kinds %in% c(3L, 4L)))
})

test_that("triple is an alias for all fake null mode", {
  f <- make_temporal_fixture()
  res_all <- NULL
  res_triple <- NULL
  res_all <- ComputeFakeLRNull(
      f$x, f$ct, f$lr_db,
      candidate_fake_genes = 5:20,
      pseudotime = f$pt,
      lineage_weights = f$lw,
      n_bins = 4L,
      n_fake_per_real = 100L,
      n_strata = 1L,
      caliper_k = 10L,
      max_match_distance = 2,
      fake_null_mode = "all",
      min_pct = 0,
      seed = 55L,
      min_fake_observed_generic = 0L,
      min_fake_coverage_generic = 0,
      min_fake_observed_edge = 0L,
      min_fake_coverage_edge = 0,
      min_unique_fake_draws = 0L,
      min_unique_fake_fraction = 0
  )
  res_triple <- ComputeFakeLRNull(
      f$x, f$ct, f$lr_db,
      candidate_fake_genes = 5:20,
      pseudotime = f$pt,
      lineage_weights = f$lw,
      n_bins = 4L,
      n_fake_per_real = 100L,
      n_strata = 1L,
      caliper_k = 10L,
      max_match_distance = 2,
      fake_null_mode = "triple",
      min_pct = 0,
      seed = 55L,
      min_fake_observed_generic = 0L,
      min_fake_coverage_generic = 0,
      min_fake_observed_edge = 0L,
      min_fake_coverage_edge = 0,
      min_unique_fake_draws = 0L,
      min_unique_fake_fraction = 0
  )
  expect_true(all(res_triple$trajectory_calibration$fake_null_mode == "all"))
  expect_equal(
    res_triple$trajectory_calibration[, c(
      "n_fake_expected_generic", "n_fake_expected_edge", "n_fake_expected_network"
    )],
    res_all$trajectory_calibration[, c(
      "n_fake_expected_generic", "n_fake_expected_edge", "n_fake_expected_network"
    )]
  )
  expect_equal(res_triple$fake_pairs, res_all$fake_pairs)
})

test_that("all fake null modes produce fake draws", {
  f <- make_temporal_fixture()
  modes <- c(
    "context_generic", "edge_rewire", "both",
    "lr_network_rewire", "all"
  )
  for (mode in modes) {
    res <- NULL
    res <- ComputeFakeLRNull(
        f$x, f$ct, f$lr_db,
        candidate_fake_genes = 5:20,
        pseudotime = f$pt,
        lineage_weights = f$lw,
        n_bins = 4L,
        n_fake_per_real = 100L,
        n_strata = 1L,
        caliper_k = 10L,
        max_match_distance = 2,
        fake_null_mode = mode,
        min_pct = 0,
        seed = 66L,
        min_fake_observed_generic = 0L,
        min_fake_coverage_generic = 0,
        min_fake_observed_edge = 0L,
        min_fake_coverage_edge = 0,
        min_unique_fake_draws = 0L,
        min_unique_fake_fraction = 0
    )
    expect_gt(nrow(res$fake_pairs), 0L)
    expect_gt(nrow(res$trajectory_calibration), 0L)
  }
})
