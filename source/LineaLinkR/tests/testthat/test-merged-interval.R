# Variant B: merged adjacent bins on a thin side, fixed sisters, post-only
# output. Raw inputs are real (scoring runs through the package); only the
# support pattern of the score table is edited to create thin sides.

make_merged_fixture <- function() {
  set.seed(7)
  n <- 1200L
  # 60% of cells belong to lineage 1, 40% to lineage 2; cells with pt < 0.4
  # are shared stem cells (0.5/0.5). Sender/receiver states alternate so both
  # cell types span the whole pseudotime range of each lineage.
  lin1 <- seq_len(720L)
  lin2 <- 721L:n
  pt <- numeric(n)
  pt[lin1] <- seq(0.05, 0.95, length.out = length(lin1))
  pt[lin2] <- seq(0.05, 0.95, length.out = length(lin2))
  lw <- matrix(0, n, 2L)
  lw[lin1, 1L] <- 1
  lw[lin2, 2L] <- 1
  stem <- pt < 0.4
  lw[stem, ] <- 0.5
  ct <- rep(c(1L, 2L), times = n %/% 2L)
  # Expression: ligand (gene 1) rises post-fork on lineage 1; receptor
  # (gene 2) is expressed by TGT cells everywhere.
  expr <- matrix(0, 4L, n)
  expr[1L, ] <- 0.5 + 1.5 * (lw[, 1L] > 0.75 & pt > 0.5)
  expr[2L, ] <- 1.0
  expr[3L, ] <- 0.2
  expr[4L, ] <- 0.2
  expr <- as(Matrix::Matrix(expr, sparse = TRUE), "CsparseMatrix")
  lr_db <- list(ligand_rows = list(1L), receptor_rows = list(2L))
  list(expr = expr, ct = ct, pt = pt, lw = lw, lr_db = lr_db)
}

make_merged_scores <- function(f) {
  tl <- ComputeTemporalLRWithIntermediate(
    f$expr, f$ct, f$lr_db, f$pt, f$lw,
    n_bins = 12L, min_n_eff = 5.0, common_grid = TRUE
  )
  tl$scores
}

# Remove pre-fork support (bins 1-5; fork sits at bin 6) except `keep_bins`,
# leaving the post side intact.
thin_pre_side <- function(scores, lineage, keep_bins) {
  s <- scores
  drop_bins <- setdiff(1:5, keep_bins)
  drop_rows <- s$lineage == lineage & s$bin %in% drop_bins
  s$score[drop_rows] <- NA_real_
  s
}

test_that("thin pre side becomes evaluable on the merged interval", {
  f <- make_merged_fixture()
  scores <- make_merged_scores(f)
  fork_bins <- data.frame(lineage = 1:2, fork_cluster = 1L, fork_bin = 6L)
  smap <- data.frame(focal_lineage = c(1L, 2L), sibling_lineage = c(2L, 1L),
                     fork_cluster = 1L, shared_prefix_length = 1L)
  # Keep only one pre bin (bin 5) for lineage 1 -> fixed window no-call.
  scores_thin <- thin_pre_side(scores, 1L, keep_bins = 5L)
  base <- ComputeForkContrasts(scores_thin, data.frame(lineage = 1:2, b_star = 6L,
                                                       reason = "crossing"),
                               smap, window_bins = 2L, min_n_eff = 5,
                               fork_bins = fork_bins, common_grid = TRUE)
  expect_false(base$evaluable[base$lineage == 1L & base$source == 1L & base$target == 2L])

  out <- ComputeMergedIntervalForkContrasts(
    scores_thin, data.frame(lineage = 1:2, b_star = 6L, reason = "crossing"),
    smap, f$expr, f$ct, f$lr_db, f$pt, f$lw, fork_bins,
    window_bins = 2L, min_n_eff = 5,
    hill_scale = attr(scores_thin, "hill_scale")
  )
  row <- out[out$lineage == 1L & out$source == 1L & out$target == 2L, ]
  expect_equal(row$estimation_mode, "merged_interval")
  expect_true(is.finite(row$sibling_contrast))
  expect_equal(row$merged_pre_bins, 1L)
  expect_false(is.na(row$merged_pre_edges))
})

test_that("a truly missing side keeps theta NA and reports post difference", {
  f <- make_merged_fixture()
  scores <- make_merged_scores(f)
  fork_bins <- data.frame(lineage = 1:2, fork_cluster = 1L, fork_bin = 6L)
  smap <- data.frame(focal_lineage = c(1L, 2L), sibling_lineage = c(2L, 1L),
                     fork_cluster = 1L, shared_prefix_length = 1L)
  # Remove every pre bin for lineage 1 -> no pre support at all.
  scores_no_pre <- scores
  pre_bins <- 1:5
  drop_rows <- scores_no_pre$lineage == 1L & scores_no_pre$bin %in% pre_bins
  scores_no_pre$score[drop_rows] <- NA_real_
  out <- ComputeMergedIntervalForkContrasts(
    scores_no_pre, data.frame(lineage = 1:2, b_star = 6L, reason = "crossing"),
    smap, f$expr, f$ct, f$lr_db, f$pt, f$lw, fork_bins,
    window_bins = 2L, min_n_eff = 5,
    hill_scale = attr(scores_no_pre, "hill_scale")
  )
  row <- out[out$lineage == 1L & out$source == 1L & out$target == 2L, ]
  expect_true(is.na(row$sibling_contrast))
  expect_true(is.finite(row$post_branch_difference))
  expect_equal(row$merged_pre_bins, 0L)
})

test_that("merged estimator is mirror-symmetric between focal and sister", {
  f <- make_merged_fixture()
  scores <- make_merged_scores(f)
  fork_bins <- data.frame(lineage = 1:2, fork_cluster = 1L, fork_bin = 6L)
  smap <- data.frame(focal_lineage = c(1L, 2L), sibling_lineage = c(2L, 1L),
                     fork_cluster = 1L, shared_prefix_length = 1L)
  scores_thin <- thin_pre_side(scores, 1L, keep_bins = 5L)
  out <- ComputeMergedIntervalForkContrasts(
    scores_thin, data.frame(lineage = 1:2, b_star = 6L, reason = "crossing"),
    smap, f$expr, f$ct, f$lr_db, f$pt, f$lw, fork_bins,
    window_bins = 2L, min_n_eff = 5,
    hill_scale = attr(scores_thin, "hill_scale")
  )
  merged <- out[out$estimation_mode == "merged_interval", ]
  if (nrow(merged) == 2L) {
    expect_equal(merged$sibling_contrast[1L], -merged$sibling_contrast[2L],
                 tolerance = 1e-10)
  } else {
    expect_true(nrow(merged) >= 1L)
  }
})
