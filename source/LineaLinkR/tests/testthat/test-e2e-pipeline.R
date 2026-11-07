# End-to-end tests using the pancreas benchmark RDS.
# Covers ComputeLineageEntropyRaw, SmoothMA3, FindBifurcationBins,
# ComputeRoleScores, and MergeEvidence.

PANCREAS_RDS <- Sys.getenv("LINEALINKR_PANCREAS_RDS", unset = "")

skip_if_no_pancreas <- function() {
  if (!nzchar(PANCREAS_RDS) || !file.exists(PANCREAS_RDS)) {
    skip("Set LINEALINKR_PANCREAS_RDS to run pancreas e2e tests")
  }
}

# ── Load fixture once ───────────────────────────────────────────────────
load_pancreas <- local({
 cache <- NULL
 function() {
 skip_if_no_pancreas()
 if (is.null(cache)) cache <<- readRDS(PANCREAS_RDS)
 cache
 }
})

# ── E2E: entropy → smooth → bifurcation reproduces stored bif_df ───────
test_that("E2E: entropy pipeline reproduces stored bif_df on pancreas", {
 px <- load_pancreas()
 inter <- px$temporal_lr$intermediate

 h_raw <- ComputeLineageEntropyRaw(inter)
 expect_true(is.matrix(h_raw))
 expect_equal(nrow(h_raw), inter$n_lineages)
 expect_equal(ncol(h_raw), inter$n_bins)

 # Compare with stored entropy_raw
 expect_equal(dim(h_raw), dim(px$entropy_raw))
 expect_true(all(abs(h_raw - px$entropy_raw) < 1e-5 | (is.nan(h_raw) & is.nan(px$entropy_raw))),
 info = "Recomputed h_raw should match stored entropy_raw")

 h_smooth <- SmoothMA3(h_raw)
 expect_equal(dim(h_smooth), dim(px$h_tilde))
 expect_true(all(abs(h_smooth - px$h_tilde) < 1e-5 | (is.nan(h_smooth) & is.nan(px$h_tilde))),
 info = "Recomputed h_tilde should match stored h_tilde")

 bif <- FindBifurcationBins(h_smooth, h_fork = 0.5)
 expect_equal(nrow(bif), inter$n_lineages)
 expect_equal(bif$b_star, px$bif_df$b_star)
 expect_equal(bif$reason, px$bif_df$reason)
})

# ── E2E: ComputeRoleScores on pancreas temporal_lr ─────────────────────
test_that("E2E: ComputeRoleScores produces valid output on pancreas", {
 px <- load_pancreas()
 scores_df <- px$temporal_lr$scores
 bif_df <- px$bif_df
 n <- nrow(scores_df)

 # Minimal spec/delta_tau for smoke testing (uniform values)
 role <- ComputeRoleScores(
 scores_df,
 spec_self = rep(0.5, n),
 spec_max_other = rep(0.3, n),
 delta_tau = rep(0.05, n),
 bif_df = bif_df,
 tau_0 = 1.0,
 h_fork = 0.5,
 min_n_eff_cells = 5.0
 )
 expect_true(is.data.frame(role))
 expect_equal(nrow(role), n)
 expect_true(all(c("instructive", "reinforcing", "feedback",
 "lateral", "bystander") %in% colnames(role)))

 # bystander must be in [0, 1] for finite rows
 finite_rows <- is.finite(role$bystander)
 expect_true(all(role$bystander[finite_rows] >= -1e-6))
 expect_true(all(role$bystander[finite_rows] <= 1.0 + 1e-6))

 # Rust evaluates the gate in f32; exclude values within one rounding unit of 5.
 low_neff <- is.finite(scores_df$n_eff_cells) & scores_df$n_eff_cells < 5.0 - 1e-6
 if (any(low_neff)) {
 expect_true(all(role$bystander[low_neff] == 1.0),
 info = " gate: low n_eff rows must be bystander=1")
 expect_true(all(role$instructive[low_neff] == 0.0))
 }
})

# ── E2E: MergeEvidence on pancreas temporal_lr (structural only) ───────
test_that("E2E: MergeEvidence structural_only on pancreas", {
 px <- load_pancreas()
 scores_df <- px$temporal_lr$scores
 out <- MergeEvidence(scores_df)
 expect_true(is.data.frame(out))
 expect_equal(nrow(out), nrow(scores_df))
 expect_equal(out$calibration_level[1], "structural_only")

 # n_eff gate
 expect_true(all(out$pass_n_eff == (scores_df$n_eff_cells >= 5.0)))
 # entropy gate (biological; separate from assignment confidence)
 expect_true(all(out$pass_entropy == (scores_df$lineage_entropy <= 0.8)))
 expect_true(all(is.na(out$pass_assignment)))

 # Missing assignment confidence does not block the n_eff structural gate.
 expect_equal(out$structurally_eligible, out$pass_n_eff)
 expect_true(all(is.na(out$high_confidence)))
 expect_false(any(out$prioritized_candidate %in% TRUE))

 cat(" MergeEvidence: ", sum(out$prioritized_candidate), "/", nrow(out),
 " rows prioritized (structural_only)\n")
})

# ── E2E: MergeEvidence with synthetic p_fake + p_temporal ─────────────
test_that("E2E: MergeEvidence full calibration with synthetic p-values", {
 px <- load_pancreas()
 scores_df <- px$temporal_lr$scores
 key_cols <- c("pair_idx", "source", "target", "lineage")
 cal <- unique(scores_df[key_cols])

 # Synthetic trajectory-level FakeLR and temporal-null calibrations.
 cal$p_fake <- runif(nrow(cal), min = 0.001, max = 0.5)
 cal$q_fake <- stats::p.adjust(cal$p_fake, method = "BH")
 cal$trajectory_evaluable <- TRUE
 fake_cal <- cal
 temporal_cal <- cal[key_cols]
 temporal_cal$p_temporal <- runif(nrow(cal), min = 0.001, max = 0.5)
 temporal_cal$q_temporal <- stats::p.adjust(temporal_cal$p_temporal, method = "BH")
 temporal_cal$trajectory_evaluable <- TRUE
 temporal_cal$calibration_hypothesis <- "temporal_variation"
 temporal_cal$trajectory_stat_name <- "amplitude"

 out <- MergeEvidence(scores_df,
                      fake_calibration = fake_cal,
                      temporal_calibration = temporal_cal)
 expect_equal(out$calibration_level[1], "fully_calibrated")
 expect_true("q_fake" %in% colnames(out))
 expect_true("q_temporal" %in% colnames(out))
 expect_true("pass_fake" %in% colnames(out))
 expect_true("pass_temporal" %in% colnames(out))

 # high_confidence requires ALL gates
 hc <- out$pass_n_eff & out$pass_entropy &
       out$pass_fake & out$pass_temporal
 expect_true(all(out$high_confidence == hc))
})

# ── E2E: ComputeFakePValues pure-R logic ───────────────────────────────
test_that("E2E: ComputeFakePValues Laplace smoothing on pancreas-shaped data", {
 px <- load_pancreas()
 scores_df <- px$temporal_lr$scores
 # Take 20 finite rows with unique full keys for a one-to-one null match.
 key_cols <- c("pair_idx", "source", "target", "lineage", "bin")
 keep <- is.finite(scores_df$score) & !duplicated(scores_df[key_cols])
 real <- head(scores_df[keep, , drop = FALSE], 20)

 # Construct matching fake scores (same key columns, different scores)
 set.seed(42)
 fake_rows <- do.call(rbind, lapply(seq_len(nrow(real)), function(i) {
 data.frame(pair_idx = real$pair_idx[i], source = real$source[i],
 target = real$target[i], lineage = real$lineage[i],
 bin = real$bin[i],
 score = real$score[i] * runif(5, 0.5, 1.5))
 }))

 pf <- ComputeFakePValues(real, fake_rows, n_fake_total = nrow(fake_rows))
 expect_length(pf, nrow(real))
 expect_true(all(pf > 0 & pf <= 1))

 # Five matched fake rows per real row: minimum p is 1/(5+1).
 expect_true(all(pf >= 1 / 6))
})
