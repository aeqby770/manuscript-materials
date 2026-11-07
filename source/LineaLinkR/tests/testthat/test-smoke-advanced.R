# Smoke tests for advanced LineaLinkR functions.
# Tests: ComputeLineageWeightsTrunkAware, ComputeTemporalLRWithIntermediate,
# ComputeLineageEntropyRaw, SmoothMA3, FindBifurcationBins,
# ComputeRoleScores, ComputeFakePValues, MergeEvidence, TrajectoryDE.

# ── Shared fixtures ──────────────────────────────────────────────────────
make_temporal_fixture <- function() {
 skip_if_not_installed("Matrix")
 n_genes <- 20; n_cells <- 40; n_types <- 3
 set.seed(42)
 x <- Matrix::rsparsematrix(n_genes, n_cells, density = 0.3)
 x@x <- abs(x@x)
 ct <- sample.int(n_types, n_cells, replace = TRUE)
 pt <- seq(0, 1, length.out = n_cells)
 lw <- matrix(c(rep(c(0.7, 0.3), n_cells / 2)), nrow = n_cells, ncol = 2,
 byrow = TRUE)
 lr_db <- list(
 ligand_rows = list(1L, 3L),
 receptor_rows = list(2L, 4L)
 )
 list(x = x, ct = ct, pt = pt, lw = lw, lr_db = lr_db,
 n_genes = n_genes, n_cells = n_cells, n_types = n_types)
}

# ── ComputeLineageWeightsTrunkAware ─────────────────────────────────────
test_that("ComputeLineageWeightsTrunkAware returns matrix with rows summing to 1", {
 set.seed(42)
 n_cells <- 50; n_dims <- 5; n_clusters <- 4
 emb <- matrix(rnorm(n_cells * n_dims), nrow = n_cells)
 cl <- sample.int(n_clusters, n_cells, replace = TRUE)
 lw <- ComputeLineageWeightsTrunkAware(emb, cl, start_cluster = 1L, sigma = 1.0)
 expect_true(is.matrix(lw))
 expect_equal(nrow(lw), n_cells)
 row_sums <- rowSums(lw)
 expect_true(all(abs(row_sums - 1.0) < 1e-4))
})

# ── ComputeTemporalLRWithIntermediate ───────────────────────────────────
test_that("ComputeTemporalLRWithIntermediate returns scores + intermediate", {
 f <- make_temporal_fixture()
 res <- ComputeTemporalLRWithIntermediate(
 f$x, f$ct, f$lr_db, f$pt, f$lw, n_bins = 4L, min_pct = 0.0)
 expect_type(res, "list")
 expect_true("scores" %in% names(res))
 expect_true("intermediate" %in% names(res))
 expect_true(is.data.frame(res$scores))
 inter <- res$intermediate
 expect_true(all(c("hwsums", "wsums", "n_types", "n_lineages", "n_bins")
 %in% names(inter)))
 expect_equal(inter$n_bins, 4L)
})

# ── ComputeLineageEntropyRaw + SmoothMA3 ───────────────────────────────
test_that("ComputeLineageEntropyRaw and SmoothMA3 produce correct shape", {
 f <- make_temporal_fixture()
 res <- ComputeTemporalLRWithIntermediate(
 f$x, f$ct, f$lr_db, f$pt, f$lw, n_bins = 8L, min_pct = 0.0)
 h_raw <- ComputeLineageEntropyRaw(res$intermediate)
 expect_true(is.matrix(h_raw))
 expect_equal(nrow(h_raw), res$intermediate$n_lineages)
 expect_equal(ncol(h_raw), 8L)

 h_smooth <- SmoothMA3(h_raw)
 expect_true(is.matrix(h_smooth))
 expect_equal(dim(h_smooth), dim(h_raw))
})

test_that("SmoothMA3 bypasses when n_bins < 6", {
 h_raw <- matrix(c(0.1, 0.2, 0.3, 0.4, 0.5), nrow = 1)
 h_smooth <- SmoothMA3(h_raw)
  expect_equal(h_smooth, h_raw, tolerance = 1e-6)
})

# ── FindBifurcationBins ─────────────────────────────────────────────────
test_that("FindBifurcationBins returns data.frame with lineage/b_star/reason", {
 h <- matrix(c(0.9, 0.85, 0.7, 0.55, 0.45, 0.40, 0.35, 0.30), nrow = 1)
 bif <- FindBifurcationBins(h, h_fork = 0.5)
 expect_true(is.data.frame(bif))
 expect_true(all(c("lineage", "b_star", "reason") %in% colnames(bif)))
 expect_equal(nrow(bif), 1L)
 expect_true(bif$reason[1] %in% c("crossing", "c2_not_informative",
 "c3_no_low_entropy", "c4_boundary_onset"))
})

test_that("FindBifurcationBins returns NA b_star for flat curve (C2)", {
 h <- matrix(rep(0.48, 8), nrow = 1)
 bif <- FindBifurcationBins(h, h_fork = 0.5)
 expect_true(is.na(bif$b_star[1]))
 expect_equal(bif$reason[1], "c2_not_informative")
})

# ── ComputeRoleScores ───────────────────────────────────────────────────
test_that("ComputeRoleScores returns 6-column data.frame", {
 f <- make_temporal_fixture()
 res <- ComputeTemporalLRWithIntermediate(
 f$x, f$ct, f$lr_db, f$pt, f$lw, n_bins = 8L, min_pct = 0.0)
 scores_df <- res$scores
 if (nrow(scores_df) == 0) skip("No LR scores produced (min_pct too high)")

 h_raw <- ComputeLineageEntropyRaw(res$intermediate)
 h_smooth <- SmoothMA3(h_raw)
 bif_df <- FindBifurcationBins(h_smooth, h_fork = 0.5)

 n <- nrow(scores_df)
 role <- ComputeRoleScores(
 scores_df,
 spec_self = rep(0.5, n),
 spec_max_other = rep(0.3, n),
 delta_tau = rep(0.1, n),
 bif_df = bif_df,
 tau_0 = 1.0,
 h_fork = 0.5,
 min_n_eff_cells = 0.0
 )
 expect_true(is.data.frame(role))
 expect_equal(nrow(role), n)
 expect_true(all(c("instructive", "transition_like", "reinforcing", "feedback",
 "lateral", "bystander") %in% colnames(role)))
  expect_true(all(role$bystander >= 0, na.rm = TRUE))
})

test_that("ComputeDeltaTau and role scores use target minus source", {
  expect_equal(
    ComputeDeltaTau(source_pseudotime = 0.3, target_pseudotime = 0.6),
    0.3
  )
  scores_df <- data.frame(bin = 1L, lineage = 1L, n_eff_cells = 10)
  bif_df <- data.frame(lineage = 1L, b_star = 2L, reason = "crossing")
  sender_leads <- ComputeRoleScores(
    scores_df, spec_self = 0.8, spec_max_other = 0.1, delta_tau = 0.3,
    bif_df = bif_df, tau_0 = 0.1, min_n_eff_cells = 5
  )
  receiver_leads <- ComputeRoleScores(
    scores_df, spec_self = 0.8, spec_max_other = 0.1, delta_tau = -0.3,
    bif_df = bif_df, tau_0 = 0.1, min_n_eff_cells = 5
  )

  expect_gt(sender_leads$instructive, receiver_leads$instructive)
})

# ── ComputeFakePValues ──────────────────────────────────────────────────
test_that("ComputeFakePValues returns p-values in (0, 1]", {
 real <- data.frame(pair_idx = c(1, 1), source = c(1, 2),
 target = c(2, 1), lineage = c(1, 1),
 bin = c(1, 2), score = c(0.5, 0.3))
  fake <- data.frame(pair_idx = c(1, 1, 1), source = c(1, 1, 2),
  target = c(2, 2, 1), lineage = c(1, 1, 1),
  bin = c(1, 1, 2), score = c(0.4, 0.6, 0.1))
  pf <- ComputeFakePValues(real, fake, n_fake_total = 3L)
  expect_length(pf, nrow(real))
  expect_true(all(pf > 0 & pf <= 1))
  expect_equal(as.numeric(pf), c(2 / 3, 1 / 2), tolerance = 1e-8)
  expect_equal(attr(pf, "n_fake_observed"), c(2L, 1L))
})

test_that("ComputeFakePValues returns NA without matched fake rows", {
  real <- data.frame(pair_idx = 1, source = 1, target = 2,
                     lineage = 1, bin = 1, score = 0.5)
  fake <- data.frame(pair_idx = 1, source = 2, target = 1,
                     lineage = 1, bin = 1, score = 0.9)
  expect_true(is.na(ComputeFakePValues(real, fake)[1]))
})

test_that("ComputeTrajectoryFakePValues is deprecated for formal inference", {
  real <- data.frame(
    pair_idx = rep(1, 3), source = rep(1, 3), target = rep(2, 3),
    lineage = rep(1, 3), bin = 1:3, score = c(0.0, 1.0, 0.0)
  )
  fake <- data.frame(
    fake_id = rep(1:3, each = 3),
    pair_idx = rep(1, 9), source = rep(1, 9), target = rep(2, 9),
    lineage = rep(1, 9), bin = rep(1:3, 3),
    score = c(0.2, 0.3, 0.2, 0.1, 0.4, 0.2, 0.2, 0.5, 0.3)
  )
  expect_error(
    suppressWarnings(LineaLinkR:::ComputeTrajectoryFakePValues(
      real, fake, n_bins = 3L
    )),
    "must not be used for formal inference"
  )
})

test_that("ComputeFakeLRNull returns fake scores and pair metadata", {
  f <- make_temporal_fixture()
  res <- ComputeFakeLRNull(
    f$x, f$ct, f$lr_db,
    candidate_fake_genes = 5:20,
    pseudotime = f$pt,
    lineage_weights = f$lw,
    n_bins = 4L,
    n_fake_per_real = 100L,
    n_strata = 5L,
    caliper_k = 5L,
    min_pct = 0.0,
    seed = 1L
  )
  expect_true(is.data.frame(res$fake_scores))
  expect_true(is.data.frame(res$fake_pairs))
  expect_true(is.data.frame(res$trajectory_calibration))
  expect_true(is.data.frame(res$matching_audit))
  expect_true(is.data.frame(res$matching_balance))
  expect_true(all(c("ligand_gene", "receptor_gene", "shadowing_pair", "null_kind")
                  %in% colnames(res$fake_pairs)))
  expect_true(all(c("p_fake", "q_fake", "n_fake_observed", "fake_coverage",
                    "n_fake_observed_generic", "n_fake_observed_edge",
                    "p_replace_ligand", "p_replace_receptor", "p_edge_iu",
                    "generic_null_evaluable", "edge_null_evaluable",
                    "replace_ligand_null_evaluable", "replace_receptor_null_evaluable",
                    "n_draws", "edge_fraction")
                  %in% colnames(res$trajectory_calibration)))
  expect_equal(res$n_fake_total, nrow(res$fake_pairs))
  expect_equal(res$n_real_pairs, length(f$lr_db$ligand_rows))
  expect_true(all(res$fake_pairs$ligand_match_distance <= 0.5))
  expect_true(all(res$fake_pairs$receptor_match_distance <= 0.5))
})

test_that("temporal scores distinguish Kish ESS from weight mass", {
  f <- make_temporal_fixture()
  out <- ComputeTemporalLR(
    f$x, f$ct, f$lr_db, f$pt, f$lw, n_bins = 4L, min_pct = 0
  )
  expect_true(all(c("n_eff_cells", "weight_mass_cells") %in% names(out)))
  expect_true(all(out$n_eff_cells >= out$weight_mass_cells - 1e-6))
})

test_that("temporal permutation preserves supplied blocks", {
  f <- make_temporal_fixture()
  blocks <- rep(c("donor_a", "donor_b"), length.out = ncol(f$x))
  out <- ComputeTemporalPermutation(
    f$x, f$ct, f$lr_db, f$pt, f$lw,
    n_bins = 4L,
    min_pct = 0,
    n_perms = 2L,
    min_perm_observed = 0L,
    block_labels = blocks
  )
  expect_true(all(out$trajectories$permutation_block_mode ==
                    "block_x_celltype_x_dominant_lineage"))
  expect_true(all(out$trajectories$n_permutation_blocks == 2L))
  expect_true(is.data.frame(out$permutation_block_audit))
  expect_equal(out$permutation_block_audit$n_blocks, 2L)
  expect_true(out$permutation_block_audit$n_shuffle_strata >= 2L)
  expect_true(out$permutation_block_audit$non_permutable_cell_fraction >= 0)
})

test_that("temporal permutation accepts block_celltype stratum", {
  f <- make_temporal_fixture()
  blocks <- rep(c("donor_a", "donor_b"), length.out = ncol(f$x))
  out <- ComputeTemporalPermutation(
    f$x, f$ct, f$lr_db, f$pt, f$lw,
    n_bins = 4L,
    min_pct = 0,
    n_perms = 2L,
    min_perm_observed = 0L,
    block_labels = blocks,
    permutation_stratum = "block_celltype"
  )
  expect_true(all(out$trajectories$permutation_stratum == "block_celltype"))
  expect_true(all(out$trajectories$permutation_block_mode == "block_x_celltype"))
  expect_true(all(out$trajectories$n_permutation_blocks == 2L))
})

test_that("FakeLR context matching respects block_labels", {
  f <- make_temporal_fixture()
  blocks <- rep(c("donor_a", "donor_b"), length.out = ncol(f$x))
  res <- NULL
  res <- ComputeFakeLRNull(
      f$x, f$ct, f$lr_db,
      candidate_fake_genes = c(3L, 4L, 5L, 6L, 7L, 8L),
      pseudotime = f$pt,
      lineage_weights = f$lw,
      n_bins = 4L,
      n_fake_per_real = 100L,
      n_strata = 2L,
      caliper_k = 2L,
      min_pct = 0,
      seed = 3L,
      min_fake_observed = 1L,
      min_fake_coverage = 0,
      block_labels = blocks
  )
  expect_equal(res$matching_audit$n_matching_blocks, 2L)
  expect_equal(res$matching_audit$matching_block_mode, "block_x_context")
})

test_that("both FakeLR mode splits odd budget before block allocation", {
  f <- make_temporal_fixture()
  blocks <- rep(c("donor_a", "donor_b", "donor_c"), length.out = ncol(f$x))
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
      fake_null_mode = "both",
      min_pct = 0,
      seed = 11L,
      min_fake_observed_generic = 0L,
      min_fake_coverage_generic = 0,
      min_fake_observed_edge = 0L,
      min_fake_coverage_edge = 0,
      min_unique_fake_draws = 0L,
      min_unique_fake_fraction = 0,
      block_labels = blocks
  )
  expect_true(all(res$trajectory_calibration$n_fake_expected_generic == 50L))
  expect_true(all(res$trajectory_calibration$n_fake_expected_edge == 50L))
  expect_true(all(res$trajectory_calibration$n_fake_observed_generic <= 50L))
  expect_true(all(res$trajectory_calibration$n_fake_observed_edge <= 50L))
  expect_true(all(res$trajectory_calibration$fake_coverage_generic <= 1,
                  na.rm = TRUE))
  expect_true(all(res$trajectory_calibration$fake_coverage_edge <= 1,
                  na.rm = TRUE))
})

test_that("FakeLR endpoint null splits ligand and receptor replacement", {
  f <- make_temporal_fixture()
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
      fake_null_mode = "both",
      min_pct = 0,
      seed = 21L,
      min_fake_observed_generic = 0L,
      min_fake_coverage_generic = 0,
      min_fake_observed_edge = 0L,
      min_fake_coverage_edge = 0,
      min_unique_fake_draws = 0L,
      min_unique_fake_fraction = 0
  )
  cal <- res$trajectory_calibration
  expect_true(all(cal$n_fake_expected_replace_ligand == 25L))
  expect_true(all(cal$n_fake_expected_replace_receptor == 25L))
  expect_true(all(
    cal$n_fake_observed_replace_ligand + cal$n_fake_observed_replace_receptor <=
      cal$n_fake_observed_edge
  ))
  both_finite <- is.finite(cal$p_replace_ligand) & is.finite(cal$p_replace_receptor)
  if (any(both_finite)) {
    expect_equal(
      cal$p_edge_iu[both_finite],
      pmax(cal$p_replace_ligand[both_finite], cal$p_replace_receptor[both_finite])
    )
  }
  edge_kinds <- unique(res$fake_pairs$null_kind)
  expect_true(all(edge_kinds %in% c(0L, 1L, 2L)))
})

test_that("lr_network_rewire null draws from lr_db pools and calibrates directionally", {
  f <- make_temporal_fixture()
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
      fake_null_mode = "lr_network_rewire",
      min_pct = 0,
      seed = 33L,
      min_fake_observed_generic = 0L,
      min_fake_coverage_generic = 0,
      min_fake_observed_edge = 0L,
      min_fake_coverage_edge = 0,
      min_unique_fake_draws = 0L,
      min_unique_fake_fraction = 0
  )
  cal <- res$trajectory_calibration
  expect_true(all(cal$n_fake_expected_network == 100L))
  expect_true(all(cal$n_fake_expected_generic == 0L))
  expect_true(all(cal$n_fake_expected_edge == 0L))
  network_kinds <- unique(res$fake_pairs$null_kind)
  expect_true(all(network_kinds %in% c(3L, 4L)))
  lig_pool <- unique(unlist(f$lr_db$ligand_rows))
  rec_pool <- unique(unlist(f$lr_db$receptor_rows))
  net_lig <- res$fake_pairs$ligand_gene[res$fake_pairs$null_kind == 4L]
  net_rec <- res$fake_pairs$receptor_gene[res$fake_pairs$null_kind == 3L]
  if (length(net_lig) > 0L) {
    expect_true(all(net_lig %in% lig_pool))
  }
  if (length(net_rec) > 0L) {
    expect_true(all(net_rec %in% rec_pool))
  }
  both_finite <- is.finite(cal$p_network_rewire_ligand) &
    is.finite(cal$p_network_rewire_receptor)
  if (any(both_finite)) {
    expect_equal(
      cal$p_network_iu[both_finite],
      pmax(cal$p_network_rewire_ligand[both_finite],
           cal$p_network_rewire_receptor[both_finite])
    )
  }
})

test_that("edge null cannot replace generic evidence when edge coverage fails", {
  f <- make_temporal_fixture()
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
      fake_null_mode = "both",
      min_pct = 0,
      min_n_eff = 0.1,
      min_valid_bins = 2L,
      min_valid_fraction = 0,
      seed = 12L,
      min_fake_observed_generic = 1L,
      min_fake_coverage_generic = 0,
      min_fake_observed_edge = 999L,
      min_fake_coverage_edge = 1,
      min_unique_fake_draws = 0L,
      min_unique_fake_fraction = 0
  )
  cal <- res$trajectory_calibration
  expect_true(all(is.na(cal$p_fake_edge)))
  keep <- is.finite(cal$p_fake_generic)
  expect_true(any(keep))
  expect_equal(cal$p_fake[keep], cal$p_fake_generic[keep])
  expect_true(all(cal$identity_inference_mode[keep] == "context_generic"))
})

test_that("FakeLR unique-draw gate rejects repeated null support", {
  f <- make_temporal_fixture()
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
      fake_null_mode = "both",
      min_pct = 0,
      min_n_eff = 0.1,
      min_valid_bins = 2L,
      min_valid_fraction = 0,
      seed = 13L,
      min_fake_observed_generic = 0L,
      min_fake_coverage_generic = 0,
      min_fake_observed_edge = 0L,
      min_fake_coverage_edge = 0,
      min_unique_fake_draws = 100L,
      min_unique_fake_fraction = 0
  )
  cal <- res$trajectory_calibration
  expect_true(all(cal$n_unique_fake_draws_generic <= 50L))
  expect_true(all(cal$n_unique_fake_draws_edge <= 50L))
  expect_true(all(is.na(cal$p_fake_generic)))
  expect_true(all(is.na(cal$p_fake_edge)))
  expect_true(all(is.na(cal$p_fake)))
})

test_that("ComputeFakeLRNull cleans and audits candidate genes", {
  f <- make_temporal_fixture()
  res <- NULL
  res <- ComputeFakeLRNull(
      f$x, f$ct, f$lr_db,
      candidate_fake_genes = c(1L, 5L, 5L, 6L, 2L, 7L),
      pseudotime = f$pt,
      lineage_weights = f$lw,
      n_bins = 4L,
      n_fake_per_real = 100L,
      n_strata = 2L,
      caliper_k = 2L,
      min_pct = 0,
      seed = 1L,
      min_fake_observed = 1L,
      min_fake_coverage = 0
  )
  expect_equal(res$n_candidate_input, 6L)
  expect_equal(res$n_candidate_unique, 5L)
  expect_equal(res$n_duplicates_removed, 1L)
  expect_equal(res$n_real_lr_removed, 2L)
  expect_equal(res$candidate_pool_size, 3L)
  expect_true(all(res$fake_pairs$ligand_gene %in% 5:7))
  expect_true(all(res$fake_pairs$receptor_gene %in% 5:7))

  expect_error(
    ComputeFakeLRNull(
      f$x, f$ct, f$lr_db,
      candidate_fake_genes = c(1L, 2L, 3L, 4L),
      pseudotime = f$pt,
      lineage_weights = f$lw,
      n_bins = 4L
    ),
    "at least two non-LR genes"
  )
})

test_that("ComputeTrajectoryFakeLRNull streams trajectory calibration", {
  f <- make_temporal_fixture()
  res <- ComputeTrajectoryFakeLRNull(
    f$x, f$ct, f$lr_db,
    candidate_fake_genes = 5:20,
    pseudotime = f$pt,
    lineage_weights = f$lw,
    n_bins = 4L,
    n_fake_per_real = 100L,
    n_strata = 5L,
    caliper_k = 5L,
    min_pct = 0.0,
    seed = 1L,
    min_fake_observed = 1L,
    min_fake_coverage = 0
  )
  expect_true(is.data.frame(res$trajectories))
  expect_true(all(res$trajectories$n_fake_observed <= 100L))
  expect_true(all(res$trajectories$fake_coverage <= 1))
  eval_p <- res$trajectories$p_fake[is.finite(res$trajectories$p_fake)]
  if (length(eval_p) > 0L) {
    expect_true(all(eval_p >= 1 / 4))
  }
})

test_that("ComputeTrajectoryFakeLRNull accepts mean_score", {
  f <- make_temporal_fixture()
  res <- NULL
  res <- ComputeTrajectoryFakeLRNull(
      f$x, f$ct, f$lr_db,
      candidate_fake_genes = 5:20,
      pseudotime = f$pt,
      lineage_weights = f$lw,
      n_bins = 4L,
      n_fake_per_real = 100L,
      n_strata = 2L,
      caliper_k = 2L,
      min_pct = 0,
      trajectory_stat = "mean_score",
      min_fake_observed = 1L,
      min_fake_coverage = 0
  )
  expect_true(is.data.frame(res$trajectories))
})

# ── EvidenceDecomp ───────────────────────────────────────────────────────
test_that("EvidenceDecomp returns aligned_r rather than partial_r", {
 f <- make_temporal_fixture()
 res <- ComputeTemporalLRWithIntermediate(
 f$x, f$ct, f$lr_db, f$pt, f$lw, n_bins = 4L, min_pct = 0.0)
 if (nrow(res$scores) == 0L) skip("No LR scores produced")
 ev <- EvidenceDecomp(
 res$scores, f$x, f$ct, f$lr_db, f$pt, f$lw, n_bins = 4L, min_pct = 0.0)
  expect_true("aligned_r" %in% colnames(ev))
  expect_true("profile_concordance" %in% colnames(ev))
  expect_true("n_microbins" %in% colnames(ev))
  expect_false("partial_r" %in% colnames(ev))
  expect_true(all(is.na(ev$z_score)))
  expect_equal(nrow(ev), nrow(res$scores))
})

test_that("EvidenceDecomp shares normalized input semantics with TemporalLR", {
  f <- make_temporal_fixture()
  res <- ComputeTemporalLRWithIntermediate(
    f$x, f$ct, f$lr_db, f$pt, f$lw, n_bins = 4L, min_pct = 0
  )
  if (nrow(res$scores) == 0L) skip("No LR scores produced")

  normalized <- EvidenceDecomp(
    res$scores, f$x, f$ct, f$lr_db, f$pt, f$lw,
    n_bins = 4L, min_pct = 0
  )
  rescaled <- EvidenceDecomp(
    res$scores, f$x, f$ct, f$lr_db, f$pt, 10 * f$lw,
    n_bins = 4L, min_pct = 0
  )
  expect_equal(rescaled, normalized, tolerance = 1e-7)
  expect_error(
    EvidenceDecomp(res$scores, f$x, f$ct, f$lr_db, f$pt, f$lw, n_bins = 1L),
    "n_bins"
  )
})

test_that("public direction and evidence helpers are exported", {
  exports <- getNamespaceExports("LineaLinkR")
  expect_true("ComputeDeltaTau" %in% exports)
  expect_true("EvidenceDecomp" %in% exports)
})

# ── MergeEvidence ───────────────────────────────────────────────────────
test_that("MergeEvidence adds expected columns", {
 df <- data.frame(
   pair_idx = 1:3,
   source = 1:3,
   target = 1:3,
   lineage = 1:3,
   n_eff_cells = c(10, 2, 15),
   lineage_entropy = c(0.3, 0.9, 0.5),
   pvalue = c(0.01, 0.5, 0.02),
   score = c(0.5, 0.6, 0.7))
  out <- MergeEvidence(df)
  expect_true(is.data.frame(out))
  expect_true(all(c("pass_n_eff", "pass_entropy", "structurally_eligible",
  "lr_identity_calibrated", "temporal_specificity_calibrated",
  "prioritized_candidate", "high_confidence",
  "calibration_level") %in% colnames(out)))
  expect_equal(out$calibration_level[1], "structural_only")
  expect_true(out$pass_n_eff[1])
  expect_false(out$pass_n_eff[2])
  expect_true(all(is.na(out$pass_assignment)))
  expect_true(out$pass_entropy[1])
  expect_false(out$pass_entropy[2])
  expect_true(out$structurally_eligible[1])
  expect_false(out$structurally_eligible[2])
  expect_true(all(is.na(out$high_confidence)))
  expect_false(any(out$prioritized_candidate %in% TRUE))
})

test_that("MergeEvidence requires fake and temporal support for high confidence", {
  df <- data.frame(
    pair_idx = 1:2,
    source = 1:2,
    target = 1:2,
    lineage = 1:2,
    n_eff_cells = c(10, 10),
    lineage_entropy = c(0.3, 0.3),
    score = c(0.8, 0.2)
  )
  pf <- c(0.001, 0.5)
  pt <- c(0.001, 0.5)
  suppressWarnings(out <- MergeEvidence(df, p_fake = pf, p_temporal = pt))
  expect_equal(out$calibration_level[1], "fully_calibrated")
  expect_true(out$high_confidence[1])
  expect_false(out$high_confidence[2])
  expect_equal(out$p_joint, pmax(pf, pt))
  expect_equal(out$q_joint, stats::p.adjust(pmax(pf, pt), method = "BH"))
  expect_equal(out$joint_supported, c(TRUE, FALSE))
  expect_true(out$fully_supported[1])
  expect_false(out$fully_supported[2])
  expect_equal(out$high_confidence, out$fully_supported)
})

test_that("MergeEvidence checks calibration hashes before merging", {
  scores <- data.frame(
    pair_idx = 1:2,
    source = 1:2,
    target = 1:2,
    lineage = 1:2,
    n_eff_cells = c(10, 10),
    lineage_entropy = c(0.2, 0.2)
  )
  fake <- scores[, c("pair_idx", "source", "target", "lineage")]
  fake$p_fake <- c(0.01, 0.2)
  fake$calibration_config_hash <- "shared"
  temporal <- scores[, c("pair_idx", "source", "target", "lineage")]
  temporal$p_temporal <- c(0.02, 0.1)
  temporal$calibration_config_hash <- "shared"

  out <- MergeEvidence(
    scores,
    fake_calibration = fake,
    temporal_calibration = temporal
  )
  expect_equal(out$p_joint, c(0.02, 0.2))
  temporal$calibration_config_hash <- "different"
  expect_error(
    MergeEvidence(
      scores,
      fake_calibration = fake,
      temporal_calibration = temporal
    ),
    "scoring configuration hashes do not match"
  )
})

test_that("MergeEvidence calibration level is assigned per row", {
  df <- data.frame(
    pair_idx = 1:2,
    source = 1:2,
    target = 1:2,
    lineage = 1:2,
    n_eff_cells = c(10, 10),
    lineage_entropy = c(0.3, 0.3)
  )
  suppressWarnings(out <- MergeEvidence(
    df,
    p_fake = c(0.001, NA_real_),
    p_temporal = c(0.001, 0.001)
  ))
  expect_equal(
    out$calibration_level,
    c("fully_calibrated", "temporal_calibrated")
  )
  expect_true(out$lr_identity_calibrated[1])
  expect_false(out$lr_identity_calibrated[2])
  expect_true(is.na(out$high_confidence[2]))
})

test_that("ComputeTemporalPermutation returns trajectory and maxT bin tables", {
  f <- make_temporal_fixture()
  out <- ComputeTemporalPermutation(
    f$x, f$ct, f$lr_db, f$pt, f$lw,
    n_bins = 4L, min_pct = 0.0,
    n_perms = 3L, seed = 7L
  )
  expect_true(is.data.frame(out$trajectories))
  expect_true(is.data.frame(out$bins))
  expect_true(all(c("p_temporal", "q_temporal", "temporal_stat",
                    "p_temporal_variation", "q_temporal_variation",
                    "trajectory_stat_name", "calibration_hypothesis",
                    "analysis_data_hash", "scoring_config_hash",
                    "calibration_config_hash", "temporal_null_config_hash",
                    "hill_scale",
                    "trajectory_evaluable", "n_valid_bins",
                    "structural_valid")
                  %in% c(colnames(out$trajectories), colnames(out$bins))))
  expect_true(all(c("p_temporal_bin", "p_temporal_bin_maxT", "temporal_deviation",
                    "structural_valid", "expression_supported")
                  %in% colnames(out$bins)))
  expect_true(all(out$bins$localization_method == "within_trajectory_maxT"))
  eval_traj <- out$trajectories[is.finite(out$trajectories$p_temporal), , drop = FALSE]
  if (nrow(eval_traj) > 0L) {
    expect_true(all(eval_traj$p_temporal >= 1 / 4))
  }
  eval_bins <- out$bins[is.finite(out$bins$p_temporal_bin), , drop = FALSE]
  if (nrow(eval_bins) > 0L) {
    expect_true(all(eval_bins$p_temporal_bin >= 1 / 4))
  }
})

test_that("FakeLR and temporal permutation expose the same shared configuration", {
  f <- make_temporal_fixture()
  temporal <- ComputeTemporalPermutation(
    f$x, f$ct, f$lr_db, f$pt, f$lw,
    n_bins = 4L, min_pct = 0,
    n_perms = 2L, seed = 7L,
    min_perm_observed = 1L,
    min_perm_coverage = 0
  )
  fake <- NULL
  fake <- ComputeFakeLRNull(
      f$x, f$ct, f$lr_db,
      candidate_fake_genes = 5:20,
      pseudotime = f$pt,
      lineage_weights = f$lw,
      n_bins = 4L,
      n_fake_per_real = 100L,
      n_strata = 2L,
      caliper_k = 2L,
      min_pct = 0,
      seed = 7L,
      min_fake_observed = 1L,
      min_fake_coverage = 0
  )
  expect_identical(
    temporal$trajectories$scoring_config_hash[1],
    fake$trajectory_calibration$scoring_config_hash[1]
  )
  expect_identical(
    temporal$trajectories$analysis_data_hash[1],
    fake$trajectory_calibration$analysis_data_hash[1]
  )
  expect_equal(
    temporal$trajectories$hill_scale[1],
    fake$trajectory_calibration$hill_scale[1]
  )
})

# ── TrajectoryDE (column name verification) ─────────────────────────────
test_that("TrajectoryDE column names are rho/pval/qval", {
 n_genes <- 5; n_cells <- 30
 set.seed(42)
 expr <- matrix(abs(rnorm(n_genes * n_cells)), nrow = n_genes)
 pt <- seq(0, 1, length.out = n_cells)
 res <- TrajectoryDE(expr, pt, min_cells = 5L)
 expect_equal(colnames(res), c("rho", "pval", "qval"))
 expect_equal(nrow(res), n_genes)
})
