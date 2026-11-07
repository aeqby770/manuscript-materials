test_that("MergeEvidence splits assignment confidence from fate entropy", {
  df <- data.frame(
    pair_idx = 1:2,
    source = 1:2,
    target = 1:2,
    lineage = 1:2,
    n_eff_cells = 10,
    lineage_entropy = c(0.9, 0.1),
    assignment_confidence = c(0.95, 0.05),
    stringsAsFactors = FALSE
  )
  out <- MergeEvidence(df, max_entropy = 0.8)
  expect_equal(out$fate_entropy, df$lineage_entropy)
  expect_equal(out$fate_ambiguity, df$lineage_entropy)
  expect_equal(out$assignment_confidence, df$assignment_confidence)
  expect_equal(out$pass_assignment, c(TRUE, FALSE))
  expect_equal(out$pass_entropy, c(FALSE, TRUE))
  expect_equal(out$structurally_eligible, out$pass_n_eff & out$pass_assignment)
})

test_that("missing assignment_confidence does not reuse fate entropy as a gate", {
  df <- data.frame(
    pair_idx = 1:2,
    source = 1:2,
    target = 1:2,
    lineage = 1:2,
    n_eff_cells = 10,
    lineage_entropy = c(0.2, 0.9),
    stringsAsFactors = FALSE
  )
  out <- MergeEvidence(df, max_entropy = 0.8)
  expect_true(all(is.na(out$assignment_confidence)))
  expect_equal(out$assignment_confidence_source, rep("missing", 2L))
  expect_true(all(is.na(out$pass_assignment)))
  expect_equal(out$pass_entropy, c(TRUE, FALSE))
  expect_true(all(out$structurally_eligible))
  expect_equal(out$fate_entropy, df$lineage_entropy)
})

test_that("MergeEvidence copies provenance hashes from calibration tables", {
  df <- data.frame(
    pair_idx = 1L,
    source = 1L,
    target = 2L,
    lineage = 1L,
    n_eff_cells = 10,
    lineage_entropy = 0.2,
    stringsAsFactors = FALSE
  )
  fake <- df[c("pair_idx", "source", "target", "lineage")]
  fake$p_fake <- 0.01
  fake$analysis_data_hash <- "llr_analysis"
  fake$scoring_config_hash <- "llr_scoring"
  fake$calibration_config_hash <- "llr_scoring"
  fake$fake_null_config_hash <- "llr_fake"
  temporal <- fake[c("pair_idx", "source", "target", "lineage",
                     "analysis_data_hash", "scoring_config_hash",
                     "calibration_config_hash")]
  temporal$p_temporal <- 0.01
  temporal$temporal_null_config_hash <- "llr_temporal"
  out <- MergeEvidence(
    df,
    fake_calibration = fake,
    temporal_calibration = temporal
  )
  expect_equal(out$analysis_data_hash, "llr_analysis")
  expect_equal(out$scoring_config_hash, "llr_scoring")
  expect_equal(unique(out$fake_null_config_hash), "llr_fake")
  expect_equal(unique(out$temporal_null_config_hash), "llr_temporal")
})

test_that("MergeEvidence rejects mismatched scoring hashes", {
  df <- data.frame(
    pair_idx = 1L,
    source = 1L,
    target = 2L,
    lineage = 1L,
    n_eff_cells = 10,
    lineage_entropy = 0.2,
    stringsAsFactors = FALSE
  )
  fake <- df[c("pair_idx", "source", "target", "lineage")]
  fake$p_fake <- 0.01
  fake$scoring_config_hash <- "a"
  temporal <- fake
  temporal$p_temporal <- 0.01
  temporal$scoring_config_hash <- "b"
  expect_error(
    MergeEvidence(df, fake_calibration = fake, temporal_calibration = temporal),
    "scoring configuration hashes do not match"
  )
})
