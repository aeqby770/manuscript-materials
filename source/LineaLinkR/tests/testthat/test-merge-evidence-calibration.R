# Key-based MergeEvidence and calibration semantics.

make_scores_df <- function(n = 4L) {
  data.frame(
    pair_idx = rep(1:2, each = 2L),
    source = rep(1:2, times = 2L),
    target = rep(1:2, times = 2L),
    lineage = rep(1:2, each = 2L),
    bin = rep(1:2, times = 2L),
    score = runif(n),
    n_eff_cells = rep(10, n),
    lineage_entropy = rep(0.2, n),
    stringsAsFactors = FALSE
  )
}

make_fake_cal <- function(scores_df, p = NULL, config_hash = "fake-h") {
  keys <- unique(scores_df[, c("pair_idx", "source", "target", "lineage")])
  if (is.null(p)) {
    p <- runif(nrow(keys), min = 0.001, max = 0.5)
  } else if (length(p) != nrow(keys)) {
    p <- rep(p, length.out = nrow(keys))
  }
  keys$p_fake <- p
  keys$q_fake <- stats::p.adjust(keys$p_fake, method = "BH")
  keys$trajectory_evaluable <- TRUE
  keys$calibration_hypothesis <- "lr_identity"
  keys$calibration_config_hash <- config_hash
  keys
}

make_temp_cal <- function(keys, p = NULL, config_hash = "temp-h") {
  if (is.null(p)) {
    p <- runif(nrow(keys), min = 0.001, max = 0.5)
  } else if (length(p) != nrow(keys)) {
    p <- rep(p, length.out = nrow(keys))
  }
  keys$p_temporal <- p
  keys$q_temporal <- stats::p.adjust(keys$p_temporal, method = "BH")
  keys$trajectory_evaluable <- TRUE
  keys$calibration_hypothesis <- "temporal_variation"
  keys$trajectory_stat_name <- "amplitude"
  keys$calibration_config_hash <- config_hash
  keys
}

test_that("MergeEvidence key order is invariant", {
  scores <- make_scores_df()
  cal <- make_fake_cal(scores)
  cal_shuf <- cal[sample(nrow(cal)), , drop = FALSE]
  out1 <- MergeEvidence(scores, fake_calibration = cal)
  out2 <- MergeEvidence(scores, fake_calibration = cal_shuf)
  expect_equal(out1$q_fake, out2$q_fake)
  expect_equal(out1$prioritized_candidate, out2$prioritized_candidate)
})

test_that("MergeEvidence preserves input row order", {
  scores <- make_scores_df()
  scores$bin <- c(2L, 1L, 2L, 1L)
  scores$score <- c(0.11, 0.22, 0.33, 0.44)
  out <- MergeEvidence(scores, fake_calibration = make_fake_cal(scores))
  expect_identical(out$bin, scores$bin)
  expect_identical(out$score, scores$score)
})

test_that("MergeEvidence rejects duplicate calibration keys", {
  scores <- make_scores_df()
  cal <- make_fake_cal(scores)
  cal <- rbind(cal, cal[1, , drop = FALSE])
  expect_error(
    MergeEvidence(scores, fake_calibration = cal),
    "duplicate trajectory keys"
  )
})

test_that("MergeEvidence missing calibration keys become NA", {
  scores <- make_scores_df()
  cal <- make_fake_cal(scores)[1, , drop = FALSE]
  out <- MergeEvidence(scores, fake_calibration = cal)
  expect_true(is.na(out$q_fake[3]))
  expect_equal(out$calibration_level[3], "missing_calibration_key")
  expect_false(out$prioritized_candidate[3])
})

test_that("trajectory q-values are copied to bins without re-BH", {
  scores <- make_scores_df()
  cal <- make_fake_cal(scores, p = c(0.001, 0.04, 0.2, 0.3))
  expected_q <- stats::p.adjust(cal$p_fake, method = "BH")
  out <- MergeEvidence(scores, fake_calibration = cal)
  for (i in seq_len(nrow(cal))) {
    rows <- out$pair_idx == cal$pair_idx[i] &
      out$source == cal$source[i] &
      out$target == cal$target[i] &
      out$lineage == cal$lineage[i]
    expect_equal(unique(out$q_fake[rows]), expected_q[i])
  }
})

test_that("calibration failure does not fall back to structural candidate", {
  df <- data.frame(
    pair_idx = 1:2,
    source = 1:2,
    target = 1:2,
    lineage = 1:2,
    n_eff_cells = 10,
    lineage_entropy = 0.2,
    stringsAsFactors = FALSE
  )
  fake_cal <- data.frame(
    pair_idx = 1:2, source = 1:2, target = 1:2, lineage = 1:2,
    p_fake = c(0.001, NA_real_),
    q_fake = c(0.002, NA_real_),
    trajectory_evaluable = c(TRUE, FALSE)
  )
  out <- MergeEvidence(df, fake_calibration = fake_cal, alpha = 0.05)
  expect_true(out$prioritized_candidate[1])
  expect_false(out$prioritized_candidate[2])
  expect_equal(out$calibration_level[2], "attempted_but_unevaluable")
})

test_that("prioritized_candidate uses OR when both calibrations present", {
  df <- data.frame(
    pair_idx = 1:2,
    source = 1:2,
    target = 1:2,
    lineage = 1:2,
    n_eff_cells = 10,
    lineage_entropy = 0.2,
    stringsAsFactors = FALSE
  )
  fake_cal <- data.frame(
    pair_idx = 1:2, source = 1:2, target = 1:2, lineage = 1:2,
    p_fake = c(0.001, 0.9),
    q_fake = c(0.002, 0.9),
    trajectory_evaluable = TRUE,
    calibration_hypothesis = "lr_identity",
    calibration_config_hash = "shared"
  )
  temp_cal <- data.frame(
    pair_idx = 1:2, source = 1:2, target = 1:2, lineage = 1:2,
    p_temporal = c(0.001, 0.001),
    q_temporal = c(0.002, 0.002),
    trajectory_evaluable = TRUE,
    calibration_hypothesis = "temporal_variation",
    trajectory_stat_name = "amplitude",
    calibration_config_hash = "shared"
  )
  out <- MergeEvidence(df, fake_calibration = fake_cal,
                       temporal_calibration = temp_cal, alpha = 0.05)
  expect_true(out$prioritized_candidate[1])
  expect_true(out$prioritized_candidate[2])
  expect_true(out$high_confidence[1])
  expect_false(out$high_confidence[2])
  expect_true(out$prioritized_candidate[2])
  expect_false(out$identity_prioritized[2])
  expect_true(out$temporal_prioritized[2])
})

test_that("unevaluable FakeLR does not suppress temporal OR support", {
  df <- data.frame(
    pair_idx = 1L,
    source = 1L,
    target = 1L,
    lineage = 1L,
    n_eff_cells = 10,
    lineage_entropy = 0.2,
    stringsAsFactors = FALSE
  )
  fake_cal <- data.frame(
    pair_idx = 1L, source = 1L, target = 1L, lineage = 1L,
    p_fake = NA_real_,
    q_fake = NA_real_,
    trajectory_evaluable = FALSE,
    calibration_hypothesis = "lr_identity",
    calibration_config_hash = "shared"
  )
  temp_cal <- data.frame(
    pair_idx = 1L, source = 1L, target = 1L, lineage = 1L,
    p_temporal = 0.001,
    q_temporal = 0.001,
    trajectory_evaluable = TRUE,
    calibration_hypothesis = "temporal_variation",
    trajectory_stat_name = "amplitude",
    calibration_config_hash = "shared"
  )
  out <- MergeEvidence(df, fake_calibration = fake_cal,
                       temporal_calibration = temp_cal, alpha = 0.05)
  expect_true(out$prioritized_candidate[1])
  expect_true(out$temporal_prioritized[1])
  expect_equal(out$calibration_level[1], "partially_calibrated")
})

test_that("q-only fake calibration is evaluable", {
  df <- data.frame(
    pair_idx = 1L,
    source = 1L,
    target = 1L,
    lineage = 1L,
    n_eff_cells = 10,
    lineage_entropy = 0.2,
    stringsAsFactors = FALSE
  )
  fake_cal <- data.frame(
    pair_idx = 1L, source = 1L, target = 1L, lineage = 1L,
    q_fake = 0.01,
    trajectory_evaluable = TRUE
  )
  out <- MergeEvidence(df, fake_calibration = fake_cal, alpha = 0.05)
  expect_true(out$identity_evaluable[1])
  expect_true(out$prioritized_candidate[1])
})

test_that("deprecated vector interface still works with trajectory keys", {
  df <- data.frame(
    pair_idx = 1:2,
    source = 1:2,
    target = 1:2,
    lineage = 1:2,
    n_eff_cells = 10,
    lineage_entropy = 0.2,
    stringsAsFactors = FALSE
  )
  suppressWarnings(
    out <- MergeEvidence(df, p_fake = c(0.001, 0.5), p_temporal = c(0.001, 0.5))
  )
  expect_equal(out$calibration_level[1], "fully_calibrated")
  expect_true(out$high_confidence[1])
})
test_that("MergeEvidence prefers edge FakeLR for joint identity inference", {
  df <- data.frame(
    pair_idx = 1:2,
    source = 1:2,
    target = 1:2,
    lineage = 1:2,
    n_eff_cells = 10,
    lineage_entropy = 0.2,
    stringsAsFactors = FALSE
  )
  fake_cal <- data.frame(
    pair_idx = 1:2, source = 1:2, target = 1:2, lineage = 1:2,
    p_fake = c(0.02, 0.02),
    p_fake_edge = c(0.001, 0.9),
    p_fake_generic = c(0.9, 0.001),
    q_fake = c(0.002, 0.002),
    q_fake_edge = c(0.002, 0.9),
    q_fake_generic = c(0.9, 0.002),
    fake_null_mode = "both",
    trajectory_evaluable = TRUE,
    calibration_hypothesis = "lr_identity",
    calibration_config_hash = "shared"
  )
  temp_cal <- data.frame(
    pair_idx = 1:2, source = 1:2, target = 1:2, lineage = 1:2,
    p_temporal = c(0.001, 0.001),
    q_temporal = c(0.002, 0.002),
    trajectory_evaluable = TRUE,
    calibration_hypothesis = "temporal_variation",
    trajectory_stat_name = "amplitude",
    calibration_config_hash = "shared"
  )
  out <- MergeEvidence(df, fake_calibration = fake_cal,
                       temporal_calibration = temp_cal, alpha = 0.05)
  expect_equal(out$identity_inference_mode, c("edge_rewire", "edge_rewire"))
  expect_true(out$lr_identity_edge_supported[1])
  expect_false(out$lr_identity_edge_supported[2])
  expect_false(out$lr_identity_generic_supported[1])
  expect_true(out$lr_identity_generic_supported[2])
  expect_true(out$lr_identity_supported[1])
  expect_false(out$lr_identity_supported[2])
  expect_equal(out$p_joint[1], max(0.001, 0.001))
  expect_equal(out$p_joint[2], max(0.9, 0.001))
  expect_true(out$fully_supported[1])
  expect_false(out$fully_supported[2])
})

test_that("MergeEvidence respects per-row identity_calibration_evaluable", {
  df <- data.frame(
    pair_idx = 1:3,
    source = 1L,
    target = 1L,
    lineage = 1L,
    n_eff_cells = 10,
    lineage_entropy = 0.2,
    stringsAsFactors = FALSE
  )
  fake_cal <- data.frame(
    pair_idx = 1:3,
    source = 1L,
    target = 1L,
    lineage = 1L,
    p_fake = c(0.01, 0.02, 0.03),
    q_fake = c(0.01, 0.02, 0.03),
    identity_calibration_evaluable = c(TRUE, FALSE, TRUE),
    calibration_config_hash = "h"
  )
  out <- MergeEvidence(df, fake_calibration = fake_cal)
  expect_equal(out$p_fake, c(0.01, NA_real_, 0.03))
  expect_equal(out$q_fake, c(0.01, NA_real_, 0.03))
})

test_that("MergeEvidence separates temporal variation from trajectory strength", {
  scores <- data.frame(
    pair_idx = 1L,
    source = 1L,
    target = 2L,
    lineage = 1L,
    n_eff_cells = 10,
    lineage_entropy = 0.2
  )
  fake <- scores[c("pair_idx", "source", "target", "lineage")]
  fake$p_fake <- 0.01
  fake$calibration_config_hash <- "shared"
  fake$calibration_hypothesis <- "lr_identity"
  temporal <- scores[c("pair_idx", "source", "target", "lineage")]
  temporal$p_temporal <- 0.01
  temporal$calibration_config_hash <- "shared"

  temporal$calibration_hypothesis <- "temporal_variation"
  variation <- MergeEvidence(
    scores,
    fake_calibration = fake,
    temporal_calibration = temporal
  )
  expect_true(variation$fully_supported)
  expect_true(is.na(variation$identity_and_strength_supported))

  temporal$calibration_hypothesis <- "trajectory_strength"
  strength <- MergeEvidence(
    scores,
    fake_calibration = fake,
    temporal_calibration = temporal
  )
  expect_true(is.na(strength$fully_supported))
  expect_true(strength$identity_and_strength_supported)
})

test_that("q_identity applies BH on mixed edge/generic primary family", {
  fake <- expand.grid(
    pair_idx = 1:100,
    source = 1L,
    target = 2L,
    lineage = 1L,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  fake$p_fake_generic <- 0.9
  fake$p_fake_edge <- NA_real_
  fake$p_fake_edge[1L] <- 0.04
  fake$p_fake_generic[1L] <- NA_real_
  fake$identity_calibration_evaluable <- TRUE

  df <- data.frame(
    pair_idx = rep(1:100, each = 2L),
    source = 1L,
    target = 2L,
    lineage = 1L,
    bin = rep(1:2, 100L),
    n_eff_cells = 10,
    lineage_entropy = 0.2,
    stringsAsFactors = FALSE
  )
  out <- MergeEvidence(df, fake_calibration = fake, alpha = 0.05)
  traj1 <- out[out$pair_idx == 1L, ][1L, ]
  expect_equal(traj1$p_identity, 0.04)
  expect_gt(traj1$q_identity, 0.05)
  expect_false(traj1$lr_identity_supported)
  expect_equal(traj1$q_fake_edge, 0.04)
})

test_that("MergeEvidence merges p_fake_complex into p_identity", {
  fake <- data.frame(
    pair_idx = 1L,
    source = 1L,
    target = 2L,
    lineage = 1L,
    p_fake = NA_real_,
    p_fake_complex = 0.03,
    q_fake_complex = 0.03,
    identity_calibration_evaluable = TRUE,
    complex_pair = TRUE,
    stringsAsFactors = FALSE
  )
  df <- data.frame(
    pair_idx = 1L,
    source = 1L,
    target = 2L,
    lineage = 1L,
    bin = 1L,
    n_eff_cells = 10,
    lineage_entropy = 0.2,
    stringsAsFactors = FALSE
  )
  out <- MergeEvidence(df, fake_calibration = fake, alpha = 0.05)
  expect_true("p_fake_complex" %in% colnames(out))
  expect_equal(out$p_identity[1L], 0.03)
  expect_equal(out$q_identity[1L], 0.03)
})

test_that("legacy unspecified temporal calibration is not variation", {
  df <- data.frame(
    pair_idx = 1L,
    source = 1L,
    target = 2L,
    lineage = 1L,
    n_eff_cells = 10,
    lineage_entropy = 0.2,
    stringsAsFactors = FALSE
  )
  temporal <- df[c("pair_idx", "source", "target", "lineage")]
  temporal$p_temporal <- 0.01
  temporal$q_temporal <- 0.02
  temporal$calibration_config_hash <- "shared"
  out <- MergeEvidence(df, temporal_calibration = temporal, alpha = 0.05)
  expect_equal(out$temporal_hypothesis, "legacy_unspecified")
  expect_false(out$temporal_variation_evaluable)
  expect_true(is.na(out$fully_supported))
})

test_that("deprecated p_temporal defaults to temporal_variation hypothesis", {
  df <- data.frame(
    pair_idx = 1:2,
    source = 1:2,
    target = 1:2,
    lineage = 1:2,
    n_eff_cells = 10,
    lineage_entropy = 0.2,
    stringsAsFactors = FALSE
  )
  suppressWarnings(
    out <- MergeEvidence(df, p_temporal = c(0.01, 0.5), alpha = 0.05)
  )
  expect_equal(out$temporal_hypothesis, rep("temporal_variation", 2))
  expect_true(out$temporal_variation_evaluable[1])
})
