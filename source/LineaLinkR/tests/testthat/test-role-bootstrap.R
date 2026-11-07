# Role bootstrap stability.

make_role_bootstrap_fixture <- function() {
  skip_if_not_installed("Matrix")
  set.seed(77)
  n_genes <- 20L
  n_cells <- 60L
  x <- Matrix::rsparsematrix(n_genes, n_cells, density = 0.3)
  x@x <- abs(x@x)
  ct <- rep(1:2, each = n_cells / 2)
  pt <- seq(0, 1, length.out = n_cells)
  lw <- matrix(c(rep(c(0.8, 0.2), n_cells / 2)), nrow = n_cells, ncol = 2,
               byrow = TRUE)
  lr_db <- list(
    ligand_rows = list(1L, 3L),
    receptor_rows = list(2L, 4L)
  )
  list(x = x, ct = ct, pt = pt, lw = lw, lr_db = lr_db)
}

test_that("ComputeRoleInference adds in_transition_window", {
  scores_df <- data.frame(
    pair_idx = 1L, source = 1L, target = 2L,
    bin = c(1L, 2L, 3L), lineage = 1L, score = 1.0, n_eff_cells = 10,
    stringsAsFactors = FALSE
  )
  bif_df <- data.frame(
    lineage = 1L, b_star = 2L, reason = "crossing",
    stringsAsFactors = FALSE
  )
  feats <- data.frame(
    scores_df,
    spec_self = c(0.7, 0.9, 0.5),
    spec_max_other = c(0.2, 0.1, 0.3),
    source_timing = 0.2,
    target_timing = 0.5,
    signed_lag = 0.3,
    delta_tau = 0.3,
    b_star = 2L,
    transition_relative_position = c(-1, 0, 1),
    bif_position = c("pre_bifurcation", "transition", "post_bifurcation"),
    fully_supported = TRUE,
    stringsAsFactors = FALSE
  )
  attr(feats, "bif_df") <- bif_df

  infer0 <- ComputeRoleInference(feats, tau_0 = 0.05, min_n_eff_cells = 0,
                                 transition_window = 0L)
  infer1 <- ComputeRoleInference(feats, tau_0 = 0.05, min_n_eff_cells = 0,
                                 transition_window = 1L)
  expect_true("in_transition_window" %in% colnames(infer0))
  expect_equal(infer0$in_transition_window, c(FALSE, TRUE, FALSE))
  expect_equal(infer1$in_transition_window, c(TRUE, TRUE, TRUE))
  expect_gt(infer1$transition_like[1L], infer0$transition_like[1L])
  expect_equal(infer0$transition_like[1L], 0)
})

test_that("ComputeRoleScores transition_window affects Rust transition_like", {
  scores_df <- data.frame(
    bin = c(1L, 2L, 3L), lineage = 1L, n_eff_cells = 10,
    stringsAsFactors = FALSE
  )
  bif_df <- data.frame(
    lineage = 1L, b_star = 2L, reason = "crossing",
    stringsAsFactors = FALSE
  )
  role0 <- ComputeRoleScores(
    scores_df,
    spec_self = c(0.7, 0.9, 0.5),
    spec_max_other = c(0.2, 0.1, 0.3),
    delta_tau = rep(0.3, 3L),
    bif_df = bif_df,
    tau_0 = 0.05,
    min_n_eff_cells = 0,
    transition_window = 0L
  )
  role1 <- ComputeRoleScores(
    scores_df,
    spec_self = c(0.7, 0.9, 0.5),
    spec_max_other = c(0.2, 0.1, 0.3),
    delta_tau = rep(0.3, 3L),
    bif_df = bif_df,
    tau_0 = 0.05,
    min_n_eff_cells = 0,
    transition_window = 1L
  )
  expect_equal(role0$transition_like[1L], 0)
  expect_gt(role1$transition_like[1L], 0)
  expect_equal(role0$transition_like[2L], role1$transition_like[2L])
})

test_that("BootstrapRoleInference returns stability for formal trajectories", {
  f <- make_role_bootstrap_fixture()
  res <- ComputeTemporalLRWithIntermediate(
    f$x, f$ct, f$lr_db, f$pt, f$lw, n_bins = 6L, min_pct = 0.0
  )
  scores <- res$scores
  skip_if(nrow(scores) == 0L, "no scores in fixture")

  merged <- MergeEvidence(scores)
  feats <- ComputeRoleFeatures(
    scores, f$lw, f$pt, f$ct,
    intermediate = res$intermediate,
    evidence = merged
  )

  boot <- BootstrapRoleInference(
    feats, merged_evidence = merged,
    n_boot = 40L, seed = 42L,
    min_n_eff_cells = 0.0
  )

  expect_type(boot, "list")
  expect_true(all(c("inference", "stability") %in% names(boot)))
  expect_equal(attr(boot, "n_boot"), 40L)
  expect_equal(attr(boot, "bootstrap_mode"), "row")

  stab <- boot$stability
  expect_true(all(c(
    "dominant_role_like", "role_bootstrap_fraction", "role_stability",
    "role_ci", "role_bootstrap_status", "n_boot"
  ) %in% colnames(stab)))

  formal <- stab$role_bootstrap_status == "sensitivity"
  exploratory <- stab$role_bootstrap_status == "exploratory"
  if (any(formal)) {
    expect_true(all(is.finite(stab$role_stability[formal])))
    expect_true(all(stab$role_stability[formal] >= 0 &
                      stab$role_stability[formal] <= 1))
    fr <- stab$role_bootstrap_fraction[formal][[1L]]
    expect_type(fr, "double")
    expect_named(fr)
    expect_equal(sum(fr), 1, tolerance = 1e-6)
    ci <- stab$role_ci[formal][[1L]]
    expect_type(ci, "list")
    expect_true(all(c("instructive_like", "transition_like") %in% names(ci)))
  }
  if (any(exploratory)) {
    expect_true(all(is.na(stab$role_stability[exploratory])))
    expect_true(all(is.na(stab$dominant_role_like[exploratory])))
  }
})

test_that("BootstrapRoleInference cell mode requires pseudotime and lw", {
  f <- make_role_bootstrap_fixture()
  scores <- ComputeTemporalLR(
    f$x, f$ct, f$lr_db, f$pt, f$lw, n_bins = 4L, min_pct = 0.0
  )
  skip_if(nrow(scores) == 0L, "no scores in fixture")
  feats <- ComputeRoleFeatures(
    scores, f$lw, f$pt, f$ct,
    h_tilde = matrix(runif(8), nrow = 2, ncol = 4)
  )

  expect_error(
    BootstrapRoleInference(feats, bootstrap_mode = "cell", n_boot = 5L),
    "pseudotime"
  )

  block_labels <- rep(c("A", "B"), each = length(f$pt) / 2)
  boot <- BootstrapRoleInference(
    feats,
    bootstrap_mode = "cell",
    pseudotime = f$pt,
    lineage_weights = f$lw,
    block_labels = block_labels,
    n_boot = 20L,
    seed = 1L,
    min_n_eff_cells = 0.0
  )
  expect_equal(attr(boot, "bootstrap_mode"), "cell")
  expect_equal(attr(boot, "role_stability_mode"), "timing_sensitivity")
  expect_gt(nrow(boot$stability), 0L)
})

test_that("BootstrapRoleInference is backward compatible with ComputeRoleInference", {
  f <- make_role_bootstrap_fixture()
  scores <- ComputeTemporalLR(
    f$x, f$ct, f$lr_db, f$pt, f$lw, n_bins = 4L, min_pct = 0.0
  )
  skip_if(nrow(scores) == 0L, "no scores in fixture")
  feats <- ComputeRoleFeatures(
    scores, f$lw, f$pt, f$ct,
    h_tilde = matrix(runif(8), nrow = 2, ncol = 4)
  )

  infer <- ComputeRoleInference(feats, min_n_eff_cells = 0.0)
  boot <- BootstrapRoleInference(feats, n_boot = 10L, min_n_eff_cells = 0.0)

  expect_equal(
    colnames(infer),
    colnames(boot$inference)
  )
  expect_equal(infer$instructive_like, boot$inference$instructive_like)
})

test_that("BootstrapRoleInference block mode runs with rebuild_context", {
  f <- make_role_bootstrap_fixture()
  res <- ComputeTemporalLRWithIntermediate(
    f$x, f$ct, f$lr_db, f$pt, f$lw, n_bins = 6L, min_pct = 0.0
  )
  scores <- res$scores
  skip_if(nrow(scores) == 0L, "no scores in fixture")

  merged <- MergeEvidence(scores)
  feats <- ComputeRoleFeatures(
    scores, f$lw, f$pt, f$ct,
    intermediate = res$intermediate,
    evidence = merged
  )
  block_labels <- rep(c("A", "B", "C"), length.out = ncol(f$x))

  boot <- BootstrapRoleInference(
    feats,
    merged_evidence = merged,
    bootstrap_mode = "block",
    block_labels = block_labels,
    rebuild_context = list(
      expr = f$x,
      celltypes = f$ct,
      lr_db = f$lr_db,
      pseudotime = f$pt,
      lineage_weights = f$lw,
      n_bins = 6L,
      use_intermediate = TRUE,
      evidence = merged
    ),
    n_boot = 15L,
    seed = 9L,
    min_n_eff_cells = 0.0
  )

  expect_equal(attr(boot, "bootstrap_mode"), "block")
  expect_equal(attr(boot, "role_stability_mode"), "conditional_role_stability")
  expect_true("n_boot_valid" %in% colnames(boot$stability))
  formal <- boot$stability$role_bootstrap_status == "conditional_role_stability"
  if (any(formal)) {
    expect_true(all(boot$stability$n_boot_valid[formal] >= 2L))
    expect_true(all(is.finite(boot$stability$role_stability[formal])))
  }
})

test_that("role bootstrap includes lateral and bystander dominance", {
  main_cols <- LineaLinkR:::.role_main_like_cols()
  expect_equal(
    main_cols,
    c("instructive_like", "transition_like", "reinforcing_like",
      "feedback_like", "lateral_like", "bystander_like")
  )
  maxima <- matrix(
    c(0, 0, 0, 0, 0, 1), nrow = 1L,
    dimnames = list(NULL, main_cols)
  )
  expect_equal(LineaLinkR:::.trajectory_dominant_role(maxima), "bystander_like")
})

test_that("block-stratified cell resampling returns requested draw count", {
  set.seed(3L)
  idx <- LineaLinkR:::.resample_cell_indices(
    12L, block_labels = rep(c("A", "B", "C"), each = 4L), n_draw = 25L
  )
  expect_length(idx, 25L)
  expect_true(all(idx >= 1L & idx <= 12L))
})
