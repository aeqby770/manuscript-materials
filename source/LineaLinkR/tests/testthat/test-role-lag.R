# Role inference and ligand-target lag.

test_that("lineage specificity is NA when fewer than two evaluable lineages", {
  scores <- data.frame(
    pair_idx = rep(1L, 3L),
    source = rep(1L, 3L),
    target = rep(2L, 3L),
    lineage = c(1L, 1L, 2L),
    bin = c(1L, 2L, 2L),
    score = c(1.0, 0.8, 0.2),
    stringsAsFactors = FALSE
  )
  spec <- LineaLinkR:::.compute_lineage_specificity(scores, n_total_lineages = 3L)
  single_bin_idx <- 1L
  dual_bin_idx <- 2L
  expect_true(is.na(spec$spec_self[single_bin_idx]))
  expect_true(is.na(spec$spec_max_other[single_bin_idx]))
  expect_equal(spec$n_evaluable_lineages[single_bin_idx], 1L)
  expect_equal(spec$n_total_lineages[single_bin_idx], 3L)
  expect_false(is.na(spec$spec_self[dual_bin_idx]))
  expect_equal(spec$spec_self[dual_bin_idx], 0.8, tolerance = 1e-8)
})

test_that("lineage specificity retains valid zero-score lineages", {
  scores <- data.frame(
    pair_idx = 1L, source = 1L, target = 2L,
    lineage = c(1L, 2L), bin = 1L, score = c(0.8, 0),
    stringsAsFactors = FALSE
  )
  spec <- LineaLinkR:::.compute_lineage_specificity(scores, n_total_lineages = 2L)

  expect_equal(spec$spec_self, c(1, 0))
  expect_equal(spec$n_evaluable_lineages, rep(2L, 2L))
  expect_equal(spec$n_supported_lineages, rep(1L, 2L))
})

test_that("ComputeRoleFeatures outputs lineage count columns", {
  skip_if_not_installed("Matrix")
  set.seed(88)
  n_cells <- 40L
  x <- Matrix::rsparsematrix(10L, n_cells, density = 0.3)
  x@x <- abs(x@x)
  ct <- rep(1:2, each = n_cells / 2)
  pt <- seq(0, 1, length.out = n_cells)
  lw <- matrix(c(rep(c(0.8, 0.2), n_cells / 2)), nrow = n_cells, ncol = 2,
               byrow = TRUE)
  lr_db <- list(ligand_rows = list(1L), receptor_rows = list(2L))
  scores <- ComputeTemporalLR(x, ct, lr_db, pt, lw, n_bins = 4L, min_pct = 0.0)
  skip_if(nrow(scores) == 0L, "no scores")
  feats <- ComputeRoleFeatures(
    scores, lw, pt, ct,
    h_tilde = matrix(runif(8), nrow = 2, ncol = 4)
  )
  expect_true(all(c(
    "n_total_lineages", "n_evaluable_lineages", "n_supported_lineages"
  ) %in% colnames(feats)))
  expect_equal(unique(feats$n_total_lineages), 2L)
})

test_that("ComputeRoleInference requires n_evaluable_lineages >= 2 for formal roles", {
  scores_df <- data.frame(
    pair_idx = 1L, source = 1L, target = 2L,
    bin = 1L, lineage = 1L, score = 1.0, n_eff_cells = 10,
    stringsAsFactors = FALSE
  )
  bif_df <- data.frame(
    lineage = 1L, b_star = 2L, reason = "crossing",
    stringsAsFactors = FALSE
  )
  feats <- data.frame(
    scores_df,
    spec_self = 0.9,
    spec_max_other = 0.1,
    n_total_lineages = 2L,
    n_evaluable_lineages = 1L,
    n_supported_lineages = 1L,
    source_timing = 0.2,
    target_timing = 0.5,
    signed_lag = 0.3,
    delta_tau = 0.3,
    b_star = 2L,
    transition_relative_position = -1,
    bif_position = "pre_bifurcation",
    fully_supported = TRUE,
    stringsAsFactors = FALSE
  )
  attr(feats, "bif_df") <- bif_df
  infer <- ComputeRoleInference(feats, tau_0 = 0.05, min_n_eff_cells = 0)
  expect_equal(infer$role_inference_status, "exploratory")
  expect_true(is.na(infer$dominant_role_like[1L]))
})

test_that("timing_source ligand_target_lag wires lag into delta_tau", {
  feats <- data.frame(
    pair_idx = 1L, source = 1L, target = 2L,
    bin = 1L, lineage = 1L, score = 1.0, n_eff_cells = 10,
    spec_self = 0.9, spec_max_other = 0.1,
    n_total_lineages = 2L, n_evaluable_lineages = 2L, n_supported_lineages = 2L,
    source_timing = 0.1, target_timing = 0.2,
    signed_lag = 0.1, delta_tau = 0.1,
    lag = 0.45,
    ligand_time = 0.2, receptor_time = 0.65,
    b_star = 3L, transition_relative_position = -2,
    bif_position = "pre_bifurcation",
    fully_supported = TRUE,
    stringsAsFactors = FALSE
  )
  attr(feats, "bif_df") <- data.frame(
    lineage = 1L, b_star = 3L, reason = "crossing",
    stringsAsFactors = FALSE
  )
  infer_tau <- ComputeRoleInference(
    feats, timing_source = "cell_state_tau",
    tau_0 = 0.05, min_n_eff_cells = 0
  )
  infer_lag <- ComputeRoleInference(
    feats, timing_source = "ligand_target_lag",
    tau_0 = 0.05, min_n_eff_cells = 0
  )
  expect_gt(infer_lag$instructive_like, infer_tau$instructive_like)
  expect_equal(attr(infer_lag, "timing_source"), "ligand_receptor_lag")
})

test_that("ComputeLigandReceptorLag exposes lag_exceeds_threshold", {
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
  expect_no_warning(
    lag <- ComputeLigandReceptorLag(
      side, lw, pt, ct, method = "peak", n_bins = 5L
    )
  )
  expect_true("lag_exceeds_threshold" %in% colnames(lag))
  expect_equal(lag$lag_exceeds_threshold, lag$lag_significant)
})

test_that("lag profiles keep NA bins instead of zero-fill for peak method", {
  side <- data.frame(
    pair_idx = rep(1L, 6L),
    source = rep(1L, 6L),
    target = rep(2L, 6L),
    lineage = rep(1L, 6L),
    bin = rep(1:3, each = 2L),
    score = rep(0.5, 6L),
    ligand_score = c(NA, NA, 0.1, 0.2, 1.0, 0.5),
    receptor_score = c(NA, NA, 0.2, 0.3, 0.4, 1.0),
    stringsAsFactors = FALSE
  )
  pt <- seq(0, 1, length.out = 20L)
  ct <- rep(1:2, each = 10L)
  lw <- matrix(1, nrow = 20L, ncol = 1L)
  lag <- suppressWarnings(ComputeLigandTargetLag(
    side, lw, pt, ct, method = "peak", n_bins = 3L
  ))
  expect_equal(nrow(lag), 1L)
  expect_true(is.finite(lag$lag_peak[1L]))
})

test_that("BootstrapRoleInference reports formal stability mode and Wilson CI", {
  skip_if_not_installed("Matrix")
  set.seed(99)
  n_cells <- 60L
  x <- Matrix::rsparsematrix(20L, n_cells, density = 0.3)
  x@x <- abs(x@x)
  ct <- rep(1:2, each = n_cells / 2)
  pt <- seq(0, 1, length.out = n_cells)
  lw <- matrix(c(rep(c(0.8, 0.2), n_cells / 2)), nrow = n_cells, ncol = 2,
               byrow = TRUE)
  lr_db <- list(ligand_rows = list(1L), receptor_rows = list(2L))
  res <- ComputeTemporalLRWithIntermediate(
    x, ct, lr_db, pt, lw, n_bins = 6L, min_pct = 0.0
  )
  skip_if(nrow(res$scores) == 0L, "no scores")
  merged <- MergeEvidence(res$scores)
  feats <- ComputeRoleFeatures(
    res$scores, lw, pt, ct,
    intermediate = res$intermediate,
    evidence = merged
  )
  boot <- BootstrapRoleInference(
    feats, merged_evidence = merged,
    n_boot = 50L, seed = 7L, min_n_eff_cells = 0.0
  )
  expect_equal(attr(boot, "role_stability_mode"), "sensitivity")
  stab <- boot$stability
  expect_true(all(c(
    "modal_role", "role_win_frequency", "role_dominant_wilson_ci"
  ) %in% colnames(stab)))
  formal <- stab$role_bootstrap_status == "sensitivity"
  if (any(formal)) {
    wilson <- stab$role_dominant_wilson_ci[formal][[1L]]
    expect_named(wilson, c("lower", "upper"))
    expect_true(wilson["lower"] <= wilson["upper"])
    ci <- stab$role_ci[formal][[1L]]
    expect_type(ci, "list")
    expect_true(all(vapply(ci, function(x) x["lower"] <= x["upper"], logical(1L))))
  }
})
