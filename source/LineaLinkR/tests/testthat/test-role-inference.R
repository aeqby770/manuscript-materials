# Role feature engineering and end-to-end inference.

make_role_fixture <- function() {
  skip_if_not_installed("Matrix")
  set.seed(77)
  n_genes <- 20L
  n_cells <- 60L
  n_types <- 2L
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

test_that("ComputeRoleFeatures derives timing and specificity", {
  f <- make_role_fixture()
  res <- ComputeTemporalLRWithIntermediate(
    f$x, f$ct, f$lr_db, f$pt, f$lw, n_bins = 6L, min_pct = 0.0
  )
  scores <- res$scores
  skip_if(nrow(scores) == 0L, "no scores in fixture")

  feats <- ComputeRoleFeatures(
    scores, f$lw, f$pt, f$ct, intermediate = res$intermediate
  )
  expect_equal(nrow(feats), nrow(scores))
  expect_true(all(c(
    "spec_self", "spec_max_other", "source_timing", "target_timing",
    "signed_lag", "delta_tau", "transition_relative_position", "bif_position"
  ) %in% colnames(feats)))
  expect_equal(feats$signed_lag, feats$delta_tau)
  specificity_evaluable <- is.finite(feats$spec_self)
  expect_true(any(specificity_evaluable))
  expect_true(all(
    feats$spec_self[specificity_evaluable] >= 0 &
      feats$spec_self[specificity_evaluable] <= 1
  ))
  expect_true(all(is.finite(feats$delta_tau)))
  expect_false(is.null(attr(feats, "bif_df")))
})

test_that("ComputeRoleFeatures specificity sums to one per (pair, s, t, bin)", {
  f <- make_role_fixture()
  scores <- ComputeTemporalLR(
    f$x, f$ct, f$lr_db, f$pt, f$lw, n_bins = 4L, min_pct = 0.0
  )
  skip_if(nrow(scores) == 0L, "no scores in fixture")
  feats <- ComputeRoleFeatures(
    scores, f$lw, f$pt, f$ct,
    h_tilde = matrix(runif(8), nrow = 2, ncol = 4)
  )
  key <- paste(feats$pair_idx, feats$source, feats$target, feats$bin, sep = ".")
  for (k in unique(key)) {
    idx <- which(key == k)
    finite_idx <- idx[is.finite(feats$spec_self[idx])]
    if (length(finite_idx) > 1L) {
      expect_equal(sum(feats$spec_self[finite_idx]), 1, tolerance = 1e-6)
    } else if (!length(finite_idx)) {
      expect_true(all(is.na(feats$spec_self[idx])))
    }
  }
})

test_that("ComputeRoleInference gates formal roles on fully_supported", {
  f <- make_role_fixture()
  res <- ComputeTemporalLRWithIntermediate(
    f$x, f$ct, f$lr_db, f$pt, f$lw, n_bins = 6L, min_pct = 0.0
  )
  scores <- res$scores
  skip_if(nrow(scores) == 0L, "no scores in fixture")

  feats <- ComputeRoleFeatures(
    scores, f$lw, f$pt, f$ct, intermediate = res$intermediate
  )
  merged <- MergeEvidence(scores)
  feats <- ComputeRoleFeatures(
    scores, f$lw, f$pt, f$ct,
    intermediate = res$intermediate,
    evidence = merged
  )

  infer <- ComputeRoleInference(feats, merged_evidence = merged,
                                min_n_eff_cells = 0.0)
  expect_true(all(c(
    "instructive_like", "transition_like", "reinforcing_like",
    "feedback_like", "bystander_like", "role_inference_status",
    "dominant_role_like", "fully_supported"
  ) %in% colnames(infer)))

  formal <- infer$role_inference_status == "formal"
  exploratory <- infer$role_inference_status == "exploratory"
  expect_true(all(is.na(infer$instructive_like[exploratory]) |
                    !is.finite(infer$fully_supported[exploratory])))
  if (any(formal)) {
    expect_true(all(is.finite(infer$instructive_like[formal]) |
                      infer$instructive_like[formal] == 0))
    expect_true(all(!is.na(infer$dominant_role_like[formal])))
  }
  if (any(exploratory)) {
    expect_true(all(is.na(infer$dominant_role_like[exploratory])))
  }
})

test_that("ComputeRoleScores remains backward compatible", {
  scores_df <- data.frame(
    bin = 1L, lineage = 1L, n_eff_cells = 10,
    stringsAsFactors = FALSE
  )
  bif_df <- data.frame(
    lineage = 1L, b_star = 2L, reason = "crossing",
    stringsAsFactors = FALSE
  )
  role <- ComputeRoleScores(
    scores_df,
    spec_self = 0.8,
    spec_max_other = 0.1,
    delta_tau = 0.3,
    bif_df = bif_df,
    tau_0 = 0.1,
    min_n_eff_cells = 0
  )
  expect_equal(
    colnames(role),
    c("instructive", "transition_like", "reinforcing", "feedback",
      "lateral", "bystander")
  )
  expect_gt(role$instructive, role$feedback)
})

test_that("ComputeRoleInference sender leads increases instructive_like", {
  scores_df <- data.frame(
    pair_idx = 1L, source = 1L, target = 2L,
    bin = 1L, lineage = 1L, score = 1.0, n_eff_cells = 10,
    stringsAsFactors = FALSE
  )
  bif_df <- data.frame(
    lineage = 1L, b_star = 3L, reason = "crossing",
    stringsAsFactors = FALSE
  )
  feats <- data.frame(
    scores_df,
    spec_self = 0.9,
    spec_max_other = 0.1,
    source_timing = 0.2,
    target_timing = 0.5,
    signed_lag = 0.3,
    delta_tau = 0.3,
    b_star = 3L,
    transition_relative_position = -2,
    bif_position = "pre_bifurcation",
    fully_supported = TRUE,
    stringsAsFactors = FALSE
  )
  attr(feats, "bif_df") <- bif_df
  sender <- ComputeRoleInference(feats, tau_0 = 0.05, min_n_eff_cells = 0)
  feats$delta_tau <- -0.3
  feats$signed_lag <- -0.3
  receiver <- ComputeRoleInference(feats, tau_0 = 0.05, min_n_eff_cells = 0)
  expect_gt(sender$instructive_like, receiver$instructive_like)
})

test_that("ComputeRoleInference attaches missing calibration columns", {
  feats <- data.frame(
    pair_idx = 1L, source = 1L, target = 2L, lineage = 1L,
    bin = 1L, n_eff_cells = 10, spec_self = 0.8,
    spec_max_other = 0.1, delta_tau = 0.2,
    fully_supported = TRUE, stringsAsFactors = FALSE
  )
  attr(feats, "bif_df") <- data.frame(
    lineage = 1L, b_star = 3L, reason = "crossing"
  )
  merged <- feats[c("pair_idx", "source", "target", "lineage")]
  merged$fully_supported <- TRUE
  merged$p_joint <- 0.01
  merged$q_joint <- 0.02
  merged$p_identity <- 0.03
  merged$q_identity <- 0.04
  merged$p_temporal <- 0.05
  merged$q_temporal <- 0.06
  out <- ComputeRoleInference(
    feats, merged_evidence = merged, tau_0 = 0.1, min_n_eff_cells = 0
  )
  expect_equal(out$q_joint, 0.02)
  expect_equal(out$q_identity, 0.04)
  expect_equal(out$q_temporal, 0.06)
})
