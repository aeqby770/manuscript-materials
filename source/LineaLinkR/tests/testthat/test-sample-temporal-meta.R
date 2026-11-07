# Sample-level temporal effects and donor meta-analysis.

make_block_fixture <- function(n_cells = 120L, n_blocks = 3L) {
  skip_if_not_installed("Matrix")
  set.seed(7L)
  n_genes <- 20L
  n_types <- 2L
  n_bins <- 4L
  x <- Matrix::rsparsematrix(n_genes, n_cells, density = 0.3)
  x@x <- abs(x@x)
  ct <- sample.int(n_types, n_cells, replace = TRUE)
  pt <- seq(0, 1, length.out = n_cells)
  lw <- matrix(c(rep(c(0.8, 0.2), length.out = n_cells)), ncol = 1L)
  block_labels <- rep(paste0("donor_", seq_len(n_blocks)),
                      each = ceiling(n_cells / n_blocks))[seq_len(n_cells)]
  lr_db <- list(
    ligand_rows = list(1L, 3L),
    receptor_rows = list(2L, 4L)
  )
  scores <- ComputeTemporalLR(
    x, ct, lr_db, pt, lw, n_bins = n_bins, min_pct = 0.0
  )
  skip_if(nrow(scores) == 0L, "no LR scores in fixture")
  fake_cal <- unique(scores[, c("pair_idx", "source", "target", "lineage")])
  fake_cal$p_fake <- runif(nrow(fake_cal), 0.001, 0.2)
  fake_cal$q_fake <- stats::p.adjust(fake_cal$p_fake, method = "BH")
  list(
    x = x, ct = ct, pt = pt, lw = lw, lr_db = lr_db,
    block_labels = block_labels, scores = scores,
    fake_cal = fake_cal, n_bins = n_bins
  )
}

donor_effects <- function(f, temporal_null_per_block = FALSE, ...) {
  ComputeTemporalEffectsBySample(
    f$scores,
    block_labels = f$block_labels,
    pseudotime = f$pt,
    lineage_weights = f$lw,
    expr = f$x,
    celltypes = f$ct,
    lr_db = f$lr_db,
    temporal_null_per_block = temporal_null_per_block,
    n_perms_block = 19L,
    min_perm_observed_block = 0L,
    min_perm_coverage_block = 0,
    ...
  )
}

test_that("ComputeTemporalEffectsBySample returns per-block trajectory rows", {
  f <- make_block_fixture()
  out <- donor_effects(
    f,
    fake_calibration = f$fake_cal,
    min_block_cells = 10L,
    min_valid_bins = 1L,
    min_n_eff = 0,
    min_bin_mass = 0
  )
  expect_true(is.data.frame(out))
  expect_true(all(c(
    "block_label", "temporal_effect", "identity_effect", "peak_time",
    "evaluable", "effect_se", "p_identity"
  ) %in% colnames(out)))
  expect_gt(length(unique(out$block_label)), 1L)
  expect_true(any(is.finite(out$temporal_effect)))
})

test_that("MetaAnalyzeTemporalEffects aggregates donors", {
  f <- make_block_fixture()
  sample_fx <- donor_effects(
    f,
    fake_calibration = f$fake_cal,
    min_block_cells = 10L,
    min_valid_bins = 1L,
    min_n_eff = 0,
    min_bin_mass = 0
  )
  meta <- MetaAnalyzeTemporalEffects(
    sample_fx, method = "fixed", include_lodo = TRUE
  )
  expect_true(all(c(
    "meta_effect", "meta_p", "q_meta", "heterogeneity_I2", "sign_consistency",
    "lodo_stability", "meta_method", "tau2", "meta_p_random", "meta_se_random"
  ) %in% colnames(meta)))
  expect_equal(meta$meta_method[1L], "fixed")
  expect_true(is.na(meta$meta_p_random[1L]))
  expect_equal(nrow(meta), length(unique(do.call(
    paste,
    c(sample_fx[, c("pair_idx", "source", "target", "lineage")], sep = "\001")
  ))))
  expect_true(any(is.finite(meta$meta_p)))
})

test_that("MetaAnalyzeTemporalEffects random method returns DL tau2", {
  f <- make_block_fixture()
  sample_fx <- donor_effects(
    f,
    fake_calibration = f$fake_cal,
    min_block_cells = 10L,
    min_valid_bins = 1L,
    min_n_eff = 0,
    min_bin_mass = 0
  )
  meta <- MetaAnalyzeTemporalEffects(sample_fx, method = "random", include_lodo = FALSE)
  expect_equal(meta$meta_method[1L], "random")
  expect_true(is.finite(meta$tau2[1L]) || identical(meta$tau2[1L], 0))
  expect_true(is.finite(meta$meta_p_random[1L]))
  expect_true(is.finite(meta$meta_se_random[1L]))
  expect_true(meta$meta_se_random[1L] >= meta$meta_se[1L])
})

test_that("pooled identity calibration is annotation-only for donors", {
  f <- make_block_fixture()
  sample_fx <- donor_effects(
    f,
    fake_calibration = f$fake_cal,
    min_block_cells = 10L,
    min_valid_bins = 1L,
    min_n_eff = 0,
    min_bin_mass = 0
  )
  expect_true(all(is.na(sample_fx$p_identity)))
  expect_true(any(is.finite(sample_fx$p_identity_pooled)))
  expect_equal(unique(sample_fx$identity_source), "pooled_annotation")
  expect_error(
    MetaAnalyzeIdentityEffects(sample_fx),
    "No evaluable block-specific"
  )
})

test_that("effect_se uses bin variance when >= 2 valid bins", {
  f <- make_block_fixture()
  out <- donor_effects(
    f,
    min_block_cells = 10L,
    min_valid_bins = 2L,
    min_n_eff = 0,
    min_bin_mass = 0
  )
  multi_bin <- out[out$n_valid_bins >= 2L & is.finite(out$temporal_effect), , drop = FALSE]
  skip_if(nrow(multi_bin) < 1L, "no multi-bin rows in fixture")
  expect_true(all(multi_bin$effect_se >= 1e-3))
  expect_true(all(is.finite(multi_bin$effect_se)))
})

test_that("donor bin mass and peak reuse weighted Rust-equivalent bounds", {
  pt <- seq(0, 1, length.out = 6L)
  lw <- cbind(
    c(0.01, 0.01, 0.01, 0.01, 0.95, 0.01),
    c(0.99, 0.99, 0.99, 0.99, 0.05, 0.99)
  )
  mass <- LineaLinkR:::.block_lineage_bin_mass(
    seq_along(pt), pt, lw, n_bins = 3L
  )
  peak <- LineaLinkR:::.block_peak_time(
    seq_along(pt), pt, lw, lineage_id = 1L, peak_bin = 2L, n_bins = 3L
  )

  expect_equal(mass[1L, ], c(0.04, 0.95, 0.01), tolerance = 1e-12)
  expect_equal(peak, pt[5L], tolerance = 1e-12)
})

test_that("LeaveOneSampleOut reports stability rows", {
  f <- make_block_fixture()
  sample_fx <- donor_effects(
    f,
    min_block_cells = 10L,
    min_valid_bins = 1L,
    min_n_eff = 0,
    min_bin_mass = 0
  )
  lodo <- LeaveOneSampleOut(sample_fx)
  expect_true(is.data.frame(lodo))
  expect_true(all(c("left_out_block_id", "meta_effect", "sign_stable") %in% colnames(lodo)))
  expect_gt(nrow(lodo), 0L)
})

test_that("BootstrapSamples returns summary quantiles", {
  f <- make_block_fixture()
  sample_fx <- donor_effects(
    f,
    min_block_cells = 10L,
    min_valid_bins = 1L,
    min_n_eff = 0,
    min_bin_mass = 0
  )
  boot <- BootstrapSamples(sample_fx, n_boot = 50L, seed = 1L)
  expect_type(boot, "list")
  expect_true(!is.null(boot$summary))
  expect_true(all(c("meta_effect_mean", "meta_effect_q025", "meta_effect_q975")
                    %in% colnames(boot$summary)))
  expect_true(!is.null(boot$replicates))
})

test_that("ComputeTemporalEffectsBySample rejects missing block_labels length", {
  f <- make_block_fixture()
  expect_error(
    ComputeTemporalEffectsBySample(
      f$scores,
      block_labels = f$block_labels[-1L],
      pseudotime = f$pt,
      lineage_weights = f$lw,
      expr = f$x,
      celltypes = f$ct,
      lr_db = f$lr_db
    ),
    "`block_labels`"
  )
})

test_that("optional per-block rescoring with expr runs", {
  f <- make_block_fixture()
  out <- donor_effects(
    f,
    min_block_cells = 10L,
    min_valid_bins = 1L,
    min_n_eff = 0,
    min_bin_mass = 0
  )
  expect_gt(sum(is.finite(out$temporal_effect)), 0L)
})

test_that("fake_null_per_block computes block-stratified p_identity", {
  f <- make_block_fixture()
  out <- donor_effects(
    f,
    temporal_null_per_block = FALSE,
    candidate_fake_genes = 10:18,
    fake_null_per_block = TRUE,
    n_fake_per_real_block = 100L,
    min_fake_observed_block = 1L,
    min_fake_coverage_block = 0,
    min_unique_fake_draws_block = 1L,
    min_unique_fake_fraction_block = 0,
    min_pct_block = 0,
    min_block_cells = 10L,
    min_valid_bins = 1L,
    min_n_eff = 0,
    min_bin_mass = 0
  )
  expect_true("identity_source" %in% colnames(out))
  expect_equal(unique(out$identity_source), "block_fake_null")
  expect_true(any(is.finite(out$p_identity)))
  expect_true(any(out$evaluable_identity))
})

test_that("fake_null_per_block requires expr lr_db and candidate genes", {
  f <- make_block_fixture()
  expect_error(
    donor_effects(
      f,
      fake_null_per_block = TRUE,
      min_block_cells = 10L
    ),
    "candidate_fake_genes"
  )
})
