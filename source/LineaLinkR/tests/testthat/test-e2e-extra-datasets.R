DENTATEGYRUS_RDS <- Sys.getenv("LINEALINKR_DENTATEGYRUS_RDS", unset = "")
STEMNET_RDS <- Sys.getenv("LINEALINKR_STEMNET_RDS", unset = "")

load_dataset_fixture <- function(path) {
  if (!nzchar(path) || !file.exists(path)) {
    skip("External benchmark RDS path not configured or unavailable")
  }
  readRDS(path)
}

check_entropy_pipeline <- function(obj) {
  inter <- obj$temporal_lr$intermediate
  h_raw <- ComputeLineageEntropyRaw(inter)
  expect_true(is.matrix(h_raw))
  expect_equal(nrow(h_raw), inter$n_lineages)
  expect_equal(ncol(h_raw), inter$n_bins)
  expect_equal(dim(h_raw), dim(obj$entropy_raw))
  expect_true(all(abs(h_raw - obj$entropy_raw) < 1e-5 | (is.nan(h_raw) & is.nan(obj$entropy_raw))))

  h_smooth <- SmoothMA3(h_raw)
  expect_equal(dim(h_smooth), dim(obj$h_tilde))
  expect_true(all(abs(h_smooth - obj$h_tilde) < 1e-5 | (is.nan(h_smooth) & is.nan(obj$h_tilde))))

  bif <- FindBifurcationBins(h_smooth, h_fork = 0.5)
  expect_equal(nrow(bif), inter$n_lineages)
  expect_equal(bif$b_star, obj$bif_df$b_star)
  expect_equal(bif$reason, obj$bif_df$reason)
}

test_that("E2E: entropy pipeline reproduces stored outputs on dentategyrus", {
  dg <- load_dataset_fixture(DENTATEGYRUS_RDS)
  check_entropy_pipeline(dg)
})

test_that("E2E: entropy pipeline reproduces stored outputs on cellrank_stemnet", {
  st <- load_dataset_fixture(STEMNET_RDS)
  check_entropy_pipeline(st)
})

test_that("E2E: MergeEvidence structural_only remains consistent on dentategyrus and stemnet", {
  dg <- load_dataset_fixture(DENTATEGYRUS_RDS)
  st <- load_dataset_fixture(STEMNET_RDS)
  for (obj in list(dg, st)) {
    scores_df <- obj$temporal_lr$scores
    out <- MergeEvidence(scores_df)
    expect_true(is.data.frame(out))
    expect_equal(nrow(out), nrow(scores_df))
    expect_true(all(out$pass_n_eff == (scores_df$n_eff_cells >= 5.0)))
    expect_true(all(out$pass_entropy == (scores_df$lineage_entropy <= 0.8)))
    expect_true(all(is.na(out$pass_assignment)))
    expect_true(all(is.na(out$high_confidence)))
    expect_false(any(out$prioritized_candidate %in% TRUE))
  }
})
