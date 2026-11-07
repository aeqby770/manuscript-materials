# Donor-level meta-analysis: null-centered effects, p-value combination, and q_meta.

test_that(".stouffer_combine_p and .fisher_combine_p require >= 2 donors", {
  expect_true(is.na(LineaLinkR:::.stouffer_combine_p(0.01)))
  expect_true(is.na(LineaLinkR:::.fisher_combine_p(0.01)))
  p <- c(0.01, 0.04, 0.2)
  st <- LineaLinkR:::.stouffer_combine_p(p)
  fi <- LineaLinkR:::.fisher_combine_p(p)
  expect_true(is.finite(st))
  expect_true(is.finite(fi))
  expect_lt(st, min(p))
  expect_lt(fi, min(p))
})

test_that("one-sided Stouffer uses upper-tail p-values", {
  p <- c(0.05, 0.05)
  expected_z <- sum(stats::qnorm(1 - p)) / sqrt(length(p))
  expected_p <- stats::pnorm(expected_z, lower.tail = FALSE)
  expect_equal(
    LineaLinkR:::.stouffer_combine_p(p, sided = "one"),
    expected_p,
    tolerance = 1e-12
  )
})

test_that(".block_trajectory_null_median yields finite median for bin scores", {
  y <- c(0.1, 0.5, 0.9, 1.2)
  med <- LineaLinkR:::.block_trajectory_null_median(
    y, stat = "amplitude", top_k = 3L, n_draws = 50L, seed = 3L
  )
  expect_true(is.finite(med))
  obs <- LineaLinkR:::.trajectory_stat_from_bins(y, "amplitude", 3L)
  expect_lte(med, obs)
})

test_that("ComputeTemporalEffectsBySample requires formal rescoring inputs", {
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

  expect_error(
    ComputeTemporalEffectsBySample(
      scores, block_labels, pt, lw,
      min_block_cells = 10L, min_valid_bins = 1L, min_n_eff = 0, min_bin_mass = 0
    ),
    "requires"
  )

  rescore <- ComputeTemporalEffectsBySample(
    scores, block_labels, pt, lw,
    expr = x, celltypes = ct, lr_db = lr_db,
    temporal_null_per_block = FALSE,
    min_block_cells = 10L, min_valid_bins = 1L, min_n_eff = 0, min_bin_mass = 0
  )
  expect_equal(unique(rescore$effect_mode), "donor_rescore")
  expect_true(any(is.finite(rescore$temporal_effect)))
})

test_that("MetaAnalyzeTemporalEffects adds q_meta and requires >= 2 donors", {
  df_one <- data.frame(
    pair_idx = 1L, source = 1L, target = 2L, lineage = 1L,
    block_id = 1L, block_label = "A",
    temporal_effect = 0.5, temporal_effect_null_centered = 0.2,
    effect_se = 0.1, p_temporal_block = 0.03,
    evaluable = TRUE, stringsAsFactors = FALSE
  )
  meta_one <- MetaAnalyzeTemporalEffects(df_one, effect_col = "temporal_effect")
  expect_true(is.na(meta_one$meta_p))
  expect_true(is.na(meta_one$q_meta))

  df_two <- rbind(
    df_one,
    transform(df_one, block_id = 2L, block_label = "B",
              temporal_effect = 0.6, temporal_effect_null_centered = 0.3,
              p_temporal_block = 0.02)
  )
  meta_two <- MetaAnalyzeTemporalEffects(
    df_two, effect_col = "temporal_effect_null_centered", include_lodo = FALSE
  )
  expect_true(is.finite(meta_two$meta_p))
  expect_true(is.finite(meta_two$q_meta))
  expect_equal(meta_two$meta_method, "fisher")
  expect_equal(meta_two$q_meta, stats::p.adjust(meta_two$meta_p, method = "BH"))
})

test_that("MetaAnalyzeTemporalEffects stouffer combines block p-values", {
  df <- data.frame(
    pair_idx = rep(1L, 3L), source = 1L, target = 2L, lineage = 1L,
    block_id = 1:3, block_label = c("A", "B", "C"),
    temporal_effect = c(0.4, 0.5, 0.6),
    temporal_effect_null_centered = c(0.1, 0.2, 0.15),
    effect_se = 0.1, p_temporal_block = c(0.04, 0.03, 0.05),
    evaluable = TRUE, stringsAsFactors = FALSE
  )
  meta <- MetaAnalyzeTemporalEffects(
    df, method = "stouffer", include_lodo = FALSE
  )
  expect_equal(meta$meta_method, "stouffer")
  expect_true(is.finite(meta$meta_p))
  expect_lt(meta$meta_p, min(df$p_temporal_block))
  expect_true(is.finite(meta$q_meta))
})

test_that("MetaAnalyzeTemporalEffects warns for random effects with < 3 donors", {
  df <- data.frame(
    pair_idx = rep(1L, 2L), source = 1L, target = 2L, lineage = 1L,
    block_id = 1:2, block_label = c("A", "B"),
    temporal_effect = c(0.4, 0.6),
    temporal_effect_null_centered = c(0.1, 0.2),
    effect_se = c(0.1, 0.12), p_temporal_block = c(0.04, 0.03),
    evaluable = TRUE, stringsAsFactors = FALSE
  )
  expect_warning(
    MetaAnalyzeTemporalEffects(
      df, method = "random", include_lodo = FALSE
    ),
    "n_samples >= 3"
  )
})

test_that("temporal_null_per_block populates p_temporal_block", {
  skip_if_not_installed("Matrix")
  set.seed(13L)
  n_cells <- 90L
  x <- Matrix::rsparsematrix(20L, n_cells, density = 0.35)
  x@x <- abs(x@x) + 0.05
  ct <- sample.int(2L, n_cells, replace = TRUE)
  pt <- seq(0, 1, length.out = n_cells)
  lw <- matrix(0.85, n_cells, 1L)
  block_labels <- rep(c("A", "B", "C"), each = 30L)
  lr_db <- list(ligand_rows = list(1L), receptor_rows = list(2L))
  scores <- ComputeTemporalLR(x, ct, lr_db, pt, lw, n_bins = 4L, min_pct = 0)
  skip_if(nrow(scores) == 0L, "no scores")

  out <- ComputeTemporalEffectsBySample(
    scores, block_labels, pt, lw,
    expr = x, celltypes = ct, lr_db = lr_db,
    temporal_null_per_block = TRUE,
    n_perms_block = 19L,
    min_perm_observed_block = 0L,
    min_perm_coverage_block = 0,
    n_null_median_draws = 30L,
    min_block_cells = 10L, min_valid_bins = 1L, min_n_eff = 0, min_bin_mass = 0
  )
  expect_true("p_temporal_block" %in% colnames(out))
  expect_true(any(is.finite(out$p_temporal_block)))
})
