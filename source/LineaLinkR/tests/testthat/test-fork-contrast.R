test_that("BuildSiblingContrastMap uses the nearest shared fork", {
  paths <- list(c(1L, 2L, 3L), c(1L, 2L, 4L), c(1L, 5L, 6L))
  out <- BuildSiblingContrastMap(paths)

  expect_true(all(c(
    "focal_lineage", "sibling_lineage", "fork_cluster",
    "shared_prefix_length"
  ) %in% colnames(out)))
  expect_equal(out$sibling_lineage[out$focal_lineage == 1L], 2L)
  expect_equal(out$sibling_lineage[out$focal_lineage == 2L], 1L)
  expect_setequal(out$sibling_lineage[out$focal_lineage == 3L], c(1L, 2L))
  expect_equal(unique(out$fork_cluster[out$focal_lineage == 1L]), 2L)
  expect_equal(unique(out$fork_cluster[out$focal_lineage == 3L]), 1L)
})

test_that("BuildSiblingContrastMap separates roots and rejects duplicates", {
  expect_equal(nrow(BuildSiblingContrastMap(list(c(1L, 2L), c(9L, 3L)))), 0L)
  forest <- BuildSiblingContrastMap(list(
    c(1L, 2L), c(1L, 3L), c(9L, 4L), c(9L, 5L)
  ))
  expect_equal(forest$focal_lineage, 1:4)
  expect_equal(forest$sibling_lineage, c(2L, 1L, 4L, 3L))
  expect_error(
    BuildSiblingContrastMap(list(c(1L, 2L), c(1L, 2L))),
    "duplicate lineage path"
  )
})

test_that("shared temporal trends cancel in sibling contrast", {
  f <- make_fork_score_fixture()
  out <- ComputeForkContrasts(
    f$scores, f$bif, f$sibling,
    window_bins = 2L, min_n_eff = 5
  )
  row <- out[out$pair_idx == 1L & out$lineage == 1L, ]
  expect_equal(row$sibling_contrast, 0, tolerance = 1e-12)
  expect_gt(row$post_fork_change, 0)
})

test_that("branch-specific activation remains after sibling subtraction", {
  f <- make_fork_score_fixture()
  out <- ComputeForkContrasts(
    f$scores, f$bif, f$sibling,
    window_bins = 2L, min_n_eff = 5
  )
  focal <- out[out$pair_idx == 2L & out$lineage == 1L, ]
  sibling <- out[out$pair_idx == 2L & out$lineage == 2L, ]
  expect_gt(focal$sibling_contrast, 0)
  expect_lt(sibling$sibling_contrast, 0)
  expect_gt(abs(focal$sibling_contrast), 1)
})

test_that("data.table scores retain all tied nearest siblings", {
  f <- make_fork_score_fixture()
  third <- f$scores[f$scores$lineage == 1L, , drop = FALSE]
  third$lineage <- 3L
  scores <- data.table::as.data.table(rbind(f$scores, third))
  bif <- rbind(
    f$bif,
    transform(f$bif[f$bif$lineage == 1L, , drop = FALSE], lineage = 3L)
  )
  sibling <- BuildSiblingContrastMap(list(
    c(1L, 2L), c(1L, 3L), c(1L, 4L)
  ))

  out <- ComputeForkContrasts(
    scores, bif, sibling,
    window_bins = 2L, min_n_eff = 5
  )

  expect_equal(nrow(out), 9L)
  expect_setequal(out$lineage, 1:3)
})

test_that("complete fork grids avoid per-trajectory R splitting", {
  scores <- expand.grid(
    pair_idx = seq_len(1000L),
    source = 1L,
    target = 2L,
    lineage = 1:2,
    bin = 1:12,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  scores$score <- ((scores$pair_idx + scores$lineage + scores$bin) %% 17) / 17
  scores$n_eff_cells <- 20
  scores$structural_valid <- TRUE
  bif <- data.frame(
    lineage = 1:2,
    b_star = 6L,
    reason = "crossing"
  )
  sibling <- data.frame(
    focal_lineage = 1:2,
    sibling_lineage = 2:1,
    fork_cluster = 1L,
    shared_prefix_length = 1L
  )

  elapsed <- system.time({
    out <- ComputeForkContrasts(
      scores, bif, sibling,
      window_bins = 2L, min_n_eff = 5
    )
  })[["elapsed"]]

  expect_equal(nrow(out), 2000L)
  expect_lt(elapsed, 3)
})

test_that("low activity and incomplete windows are not high-confidence evidence", {
  f <- make_fork_score_fixture()
  f$scores$score[f$scores$pair_idx == 3L] <- 1e-8
  out <- ComputeForkContrasts(
    f$scores, f$bif, f$sibling,
    window_bins = 2L, activity_floor = 1e-4,
    min_n_eff = 5
  )
  low <- out[out$pair_idx == 3L, ]
  expect_true(all(!low$evaluable))
  expect_true(all(low$reason == "below_activity_floor"))

  f$bif$b_star[1L] <- 1L
  bad <- ComputeForkContrasts(
    f$scores, f$bif, f$sibling,
    window_bins = 2L, min_n_eff = 5
  )
  expect_true(any(bad$reason == "insufficient_pre_bins"))
})

test_that("fork contrast respects the production structural-valid gate", {
  f <- make_fork_score_fixture()
  colnames(f$scores)[colnames(f$scores) == "structurally_valid"] <-
    "structural_valid"
  f$scores$structural_valid[
    f$scores$pair_idx == 1L & f$scores$lineage == 1L &
      f$scores$bin %in% 1:2
  ] <- FALSE
  out <- ComputeForkContrasts(
    f$scores, f$bif, f$sibling,
    window_bins = 2L, min_n_eff = 5
  )
  row <- out[out$pair_idx == 1L & out$lineage == 1L, ]
  expect_false(row$evaluable)
  expect_equal(row$reason, "insufficient_pre_bins")
})

test_that("fixed pooled bounds reproduce pooled temporal scores", {
  f <- make_temporal_fixture()
  pooled <- ComputeTemporalLRWithIntermediate(
    f$x, f$ct, f$lr_db, f$pt, f$lw,
    n_bins = 4L, min_pct = 0, min_n_eff = 0.1,
    return_complete_grid = TRUE
  )
  bounds <- attr(pooled$scores, "bin_bounds")
  fixed <- ComputeTemporalLRFixedBounds(
    f$x, f$ct, f$lr_db, f$pt, f$lw,
    bin_bounds = bounds,
    hill_scale = attr(pooled$scores, "hill_scale"),
    n_bins = 4L, min_pct = 0, min_n_eff = 0.1,
    return_complete_grid = TRUE
  )
  key <- c("pair_idx", "source", "target", "lineage", "bin")
  a <- pooled$scores[do.call(order, pooled$scores[key]), ]
  b <- fixed[do.call(order, fixed[key]), ]
  expect_equal(a[key], b[key])
  expect_equal(a$score, b$score, tolerance = 1e-7)
})

test_that("requested contexts accept data.table input", {
  f <- make_temporal_fixture()
  requested <- data.table::data.table(
    pair_idx = 1L,
    source = 1L,
    target = 2L
  )
  out <- ComputeTemporalLRWithIntermediate(
    f$x, f$ct, f$lr_db, f$pt, f$lw,
    n_bins = 4L, min_pct = 0, min_n_eff = 0.1,
    requested_contexts = requested
  )
  expect_true(nrow(out$scores) > 0L)
})

test_that("R feature engineering uses production Rust bin boundaries", {
  f <- make_temporal_fixture()
  pooled <- ComputeTemporalLRWithIntermediate(
    f$x, f$ct, f$lr_db, f$pt, f$lw,
    n_bins = 4L, min_pct = 0, min_n_eff = 0.1
  )
  expect_equal(
    LineaLinkR:::.compute_lineage_bin_bounds(f$pt, f$lw, 4L),
    attr(pooled$scores, "bin_bounds"),
    tolerance = 1e-7
  )
})

test_that("sample subsets retain pooled boundaries", {
  f <- make_temporal_fixture()
  pooled <- ComputeTemporalLRWithIntermediate(
    f$x, f$ct, f$lr_db, f$pt, f$lw,
    n_bins = 4L, min_pct = 0, min_n_eff = 0.1,
    return_complete_grid = TRUE
  )
  bounds <- attr(pooled$scores, "bin_bounds")
  donor_a <- c(seq_len(10L), 31:40)
  donor_b <- 11:30
  a <- ComputeTemporalLRFixedBounds(
    f$x[, donor_a], f$ct[donor_a], f$lr_db,
    f$pt[donor_a], f$lw[donor_a, , drop = FALSE],
    bounds, attr(pooled$scores, "hill_scale"),
    n_bins = 4L, min_pct = 0, min_n_eff = 0.1
  )
  b <- ComputeTemporalLRFixedBounds(
    f$x[, donor_b], f$ct[donor_b], f$lr_db,
    f$pt[donor_b], f$lw[donor_b, , drop = FALSE],
    bounds, attr(pooled$scores, "hill_scale"),
    n_bins = 4L, min_pct = 0, min_n_eff = 0.1
  )
  local_a <- LineaLinkR:::.compute_lineage_bin_bounds(
    f$pt[donor_a], f$lw[donor_a, , drop = FALSE], 4L
  )
  local_b <- LineaLinkR:::.compute_lineage_bin_bounds(
    f$pt[donor_b], f$lw[donor_b, , drop = FALSE], 4L
  )
  expect_false(identical(local_a, local_b))
  expect_equal(attr(a, "bin_bounds"), bounds)
  expect_equal(attr(b, "bin_bounds"), bounds)
  expect_equal(
    LineaLinkR:::.lineage_cell_bin_index(
      0.25, attr(a, "bin_bounds")[[1L]], 4L
    ),
    LineaLinkR:::.lineage_cell_bin_index(
      0.25, attr(b, "bin_bounds")[[1L]], 4L
    )
  )
})

test_that("sample support exposes a single-donor positive contrast", {
  x <- data.frame(
    pair_idx = rep(1:2, each = 4),
    source = 1L,
    target = 2L,
    lineage = 1L,
    block_label = rep(paste0("D", 1:4), 2),
    sibling_contrast = c(1.0, 0.8, 1.1, 0.9, 5.0, -0.2, -0.3, -0.1),
    evaluable = TRUE
  )
  out <- SummarizeForkContrastsBySample(x, n_boot = 200L, seed = 12345L)
  expect_gt(out$sample_support[out$pair_idx == 1L], 0.9)
  expect_equal(out$sample_support[out$pair_idx == 2L], 0.75)
  expect_lt(out$positive_sample_fraction[out$pair_idx == 2L], 0.5)
  expect_gt(out$median_contrast[out$pair_idx == 1L], 0)
  expect_lt(out$median_contrast[out$pair_idx == 2L], 0)
})

test_that("sample support accepts data.table input from the production runner", {
  x <- data.table::data.table(
    pair_idx = rep(1:2, each = 4),
    source = 1L,
    target = 2L,
    lineage = 1L,
    block_label = rep(paste0("D", 1:4), 2),
    sibling_contrast = c(1.0, 0.8, 1.1, 0.9, 5.0, -0.2, -0.3, -0.1),
    evaluable = TRUE
  )
  out <- SummarizeForkContrastsBySample(x, n_boot = 200L, seed = 12345L)
  expect_s3_class(out, "data.frame")
  expect_equal(out$sample_support[out$pair_idx == 1L], 1)
  expect_equal(out$sample_support[out$pair_idx == 2L], 0.75)
})

test_that("single-sample contrast is descriptive only", {
  x <- data.frame(
    pair_idx = 1L,
    source = 1L,
    target = 2L,
    lineage = 1L,
    block_label = "D1",
    sibling_contrast = 0.7,
    evaluable = TRUE
  )
  out <- SummarizeForkContrastsBySample(x, n_boot = 200L, seed = 12345L)
  expect_equal(out$inference_unit, "single_sample_descriptive")
  expect_true(is.na(out$ci_low) && is.na(out$ci_high))
  expect_false("p_value" %in% colnames(out))
})

test_that("evidence composition preserves missing support layers", {
  fork <- data.frame(
    pair_idx = 1:2,
    source = 1L,
    target = 2L,
    lineage = 1L,
    activity_strength = c(2, 1),
    pre_fork_enrichment = c(0.2, 0.7),
    fork_locality = c(0.5, 0.1),
    post_fork_change = c(1, -0.2),
    sibling_contrast = c(0.8, -0.1),
    post_fork_persistence = c(0.9, 0.4),
    evaluable = TRUE
  )
  response <- data.frame(
    pair_idx = 1L,
    source = 1L,
    target = 2L,
    lineage = 1L,
    receiver_response_support = 0.75
  )
  out <- ComposeForkEvidence(fork, response_support = response)
  expect_equal(out$evidence_completeness, c("fork+response", "fork_only"))
  expect_true(is.na(out$receiver_response_support[2L]))
  expect_false(any(grepl(
    "instructive|reinforcing|feedback|lateral",
    colnames(out), ignore.case = TRUE
  )))
})

test_that("evidence composition preserves pooled and sample inference units", {
  fork <- data.frame(
    pair_idx = 1L,
    source = 1L,
    target = 2L,
    lineage = 1L,
    activity_strength = 2,
    pre_fork_enrichment = 0.2,
    fork_locality = 0.5,
    post_fork_change = 1,
    sibling_contrast = 0.8,
    post_fork_persistence = 0.9,
    evaluable = TRUE,
    inference_unit = "pooled_cells_descriptive"
  )
  sample <- data.frame(
    pair_idx = 1L,
    source = 1L,
    target = 2L,
    lineage = 1L,
    sample_support = 0.8,
    inference_unit = "biological_sample_bootstrap"
  )
  out <- ComposeForkEvidence(fork, sample_support = sample)
  expect_equal(out$inference_unit, "pooled_cells_descriptive")
  expect_equal(out$sample_inference_unit, "biological_sample_bootstrap")
  expect_equal(out$evidence_completeness, "fork+sample")
})

test_that("shared-window contrast cancels a common smooth trend exactly", {
  f <- make_common_grid_fork_fixture()
  out <- ComputeForkContrasts(
    f$scores, f$bif, f$sibling,
    window_bins = 2L, min_n_eff = 5,
    fork_bins = data.frame(lineage = 1:2, fork_bin = 4L)
  )
  focal <- out[out$lineage == 1L, ]
  expect_equal(focal$sibling_contrast, 0, tolerance = 1e-12)
  expect_equal(focal$fork_time_source, "topology")
  expect_equal(focal$temporal_coordinate, "common_pseudotime_grid")
  expect_equal(focal$n_pre_bins_shared, 2L)
  expect_equal(focal$n_post_bins_shared, 2L)
})

test_that("per-lineage windows misreport the audit counterexample", {
  f <- make_common_grid_fork_fixture()
  legacy <- ComputeForkContrasts(
    f$scores, f$bif, f$sibling,
    window_bins = 2L, min_n_eff = 5, common_grid = FALSE
  )
  shared <- ComputeForkContrasts(
    f$scores, f$bif, f$sibling,
    window_bins = 2L, min_n_eff = 5,
    fork_bins = data.frame(lineage = 1:2, fork_bin = 4L)
  )
  spurious <- legacy$sibling_contrast[legacy$lineage == 1L]
  corrected <- shared$sibling_contrast[shared$lineage == 1L]
  expect_true(is.finite(spurious) && abs(spurious) > 1e-3)
  expect_equal(corrected, 0, tolerance = 1e-12)
  expect_equal(
    legacy$fork_time_source[legacy$lineage == 1L], "per_lineage"
  )
})

test_that("entropy fallback anchors sisters at the focal fork bin", {
  f <- make_common_grid_fork_fixture()
  out <- ComputeForkContrasts(
    f$scores, f$bif, f$sibling,
    window_bins = 2L, min_n_eff = 5, common_grid = TRUE
  )
  focal <- out[out$lineage == 1L, ]
  expect_equal(focal$fork_time_source, "entropy_shared")
  expect_equal(focal$b_star, 4L)
  expect_equal(focal$sibling_contrast, 0, tolerance = 1e-12)
})

test_that("all-invalid entropy anchors return no-call rows instead of aborting", {
  f <- make_common_grid_fork_fixture()
  bif_invalid <- f$bif
  bif_invalid$b_star <- NA_integer_
  bif_invalid$reason <- "c4_boundary_onset"
  out <- ComputeForkContrasts(
    f$scores, bif_invalid, f$sibling,
    window_bins = 2L, min_n_eff = 5, common_grid = TRUE
  )
  expect_equal(nrow(out), 2L)
  expect_true(all(!out$evaluable))
  expect_true(all(out$reason == "invalid_bifurcation"))
  expect_true(all(is.na(out$sibling_contrast)))
})

test_that("positive and negative shared-window effects have the right sign", {
  f <- make_common_grid_fork_fixture()
  bump <- f$scores
  bump$score[bump$lineage == 1L & bump$bin >= 5L] <-
    bump$score[bump$lineage == 1L & bump$bin >= 5L] + 2
  out <- ComputeForkContrasts(
    bump, f$bif, f$sibling,
    window_bins = 2L, min_n_eff = 5,
    fork_bins = data.frame(lineage = 1:2, fork_bin = 4L)
  )
  expect_equal(
    out$sibling_contrast[out$lineage == 1L], 2, tolerance = 1e-12
  )
  expect_lt(out$sibling_contrast[out$lineage == 2L], 0)

  drop <- f$scores
  drop$score[drop$lineage == 1L & drop$bin >= 5L] <-
    drop$score[drop$lineage == 1L & drop$bin >= 5L] - 2
  out2 <- ComputeForkContrasts(
    drop, f$bif, f$sibling,
    window_bins = 2L, min_n_eff = 5,
    fork_bins = data.frame(lineage = 1:2, fork_bin = 4L)
  )
  expect_lt(out2$sibling_contrast[out2$lineage == 1L], 0)
})

test_that("missing focal and sister support get distinct no-call reasons", {
  f <- make_common_grid_fork_fixture()
  f$scores$score[f$scores$lineage == 1L & f$scores$bin %in% 2:3] <- NA_real_
  out <- ComputeForkContrasts(
    f$scores, f$bif, f$sibling,
    window_bins = 2L, min_n_eff = 5,
    fork_bins = data.frame(lineage = 1:2, fork_bin = 4L)
  )
  expect_equal(out$reason[out$lineage == 1L], "insufficient_pre_bins")

  g <- make_common_grid_fork_fixture()
  g$scores$score[g$scores$lineage == 2L & g$scores$bin %in% 2:3] <- NA_real_
  out2 <- ComputeForkContrasts(
    g$scores, g$bif, g$sibling,
    window_bins = 2L, min_n_eff = 5,
    fork_bins = data.frame(lineage = 1:2, fork_bin = 4L)
  )
  expect_equal(
    out2$reason[out2$lineage == 1L], "insufficient_shared_pre_bins"
  )
  expect_equal(out2$n_pre_bins_shared[out2$lineage == 1L], 0L)
})

test_that("nested forks anchor sisters at each focal's own event", {
  f <- make_common_grid_fork_fixture()
  third <- f$scores[f$scores$lineage == 1L, , drop = FALSE]
  third$lineage <- 3L
  attr(third, "temporal_coordinate") <- NULL
  scores <- rbind(f$scores, third)
  attr(scores, "temporal_coordinate") <- "common_pseudotime_grid"
  scores$score[scores$lineage == 3L & scores$bin >= 5L] <-
    scores$score[scores$lineage == 3L & scores$bin >= 5L] + 1.5
  bif <- data.frame(
    lineage = 1:3, b_star = c(5L, 5L, 3L), reason = "crossing"
  )
  sibling <- BuildSiblingContrastMap(list(
    c(1L, 2L, 3L), c(1L, 2L, 4L), c(1L, 5L, 6L)
  ))
  out <- ComputeForkContrasts(
    scores, bif, sibling,
    window_bins = 2L, min_n_eff = 5, common_grid = TRUE
  )
  focal3 <- out[out$lineage == 3L, ]
  expect_equal(focal3$b_star, 3L)
  expect_equal(focal3$n_sisters_mapped, 2L)
  expect_equal(focal3$fork_cluster, 1L)
  expect_gt(focal3$sibling_contrast, 0)
})

test_that("shared-window contrast is invariant to row order", {
  f <- make_common_grid_fork_fixture()
  fb <- data.frame(lineage = 1:2, fork_bin = 4L)
  a <- ComputeForkContrasts(
    f$scores, f$bif, f$sibling,
    window_bins = 2L, min_n_eff = 5, fork_bins = fb
  )
  b <- ComputeForkContrasts(
    f$scores[sample.int(nrow(f$scores)), ], f$bif, f$sibling,
    window_bins = 2L, min_n_eff = 5, fork_bins = fb
  )
  key <- c("pair_idx", "source", "target", "lineage")
  cols <- c(
    "sibling_contrast", "post_fork_change", "sibling_change",
    "n_pre_bins_shared", "n_post_bins_shared", "reason"
  )
  expect_equal(
    a[do.call(order, a[key]), cols],
    b[do.call(order, b[key]), cols]
  )
})

test_that("a missing topology anchor is a no-call, not a zero", {
  f <- make_common_grid_fork_fixture()
  out <- ComputeForkContrasts(
    f$scores, f$bif, f$sibling,
    window_bins = 2L, min_n_eff = 5,
    fork_bins = data.frame(lineage = 2L, fork_bin = 4L)
  )
  row <- out[out$lineage == 1L, ]
  expect_false(row$evaluable)
  expect_equal(row$reason, "fork_bin_unknown")
  expect_true(is.na(row$sibling_contrast))
})

test_that("LocateForkBinsOnGrid maps fork clusters to shared bins", {
  sibling <- data.frame(
    focal_lineage = c(1L, 2L),
    sibling_lineage = c(2L, 1L),
    fork_cluster = c(2L, 2L),
    shared_prefix_length = 1L
  )
  cluster_labels <- c(1L, 1L, 2L, 2L, 2L, 3L)
  pseudotime <- c(0.10, 0.20, 0.50, 0.60, 0.70, 0.90)
  bounds <- seq(0, 1, by = 0.25)
  out <- LocateForkBinsOnGrid(
    sibling, cluster_labels, pseudotime, bounds
  )
  expect_equal(out$fork_cluster, c(2L, 2L))
  expect_equal(out$fork_pseudotime, c(0.6, 0.6), tolerance = 1e-12)
  expect_equal(out$fork_bin, c(3L, 3L))
})

test_that("fork lookup retains scorer bin indices at repeated boundaries", {
  sibling <- data.frame(focal_lineage = 1L, fork_cluster = 2L)
  bounds <- c(0, 0.25, 0.25, 0.5, 0.75, 1)
  out <- LocateForkBinsOnGrid(sibling, c(2L, 2L), c(0.25, 0.25), bounds)
  expect_equal(out$fork_bin, 3L)
  expect_equal(out$fork_bin, sum(bounds <= 0.25))
  matrix_out <- LocateForkBinsOnGrid(sibling, c(2L, 2L), c(0.25, 0.25),
    rbind(bounds, bounds))
  expect_equal(matrix_out, out)
  expect_error(LocateForkBinsOnGrid(sibling, 2L, 0.25,
    c(0, 0.5, 0.25, 1)), "non-decreasing")
})

test_that("late persistence retains observed expression-gated zeros", {
  events <- data.frame(pair_idx = 1L, source = 1L, target = 2L,
    event_node = 1L, child = 2L, fork_time = 0.25,
    sibling_children = "3", evaluable = TRUE, sibling_contrast = 0.4,
    focal_pre = 0.1, other_pre = 0.1, post_fork_change = 0.4)
  attr(events, "event_units") <- list(sibling_selection = "all")
  scores <- data.frame(pair_idx = 1L, source = 1L, target = 2L,
    lineage = rep(1:2, each = 4L), bin = rep(1:4, 2L),
    score = 0, score_positive = 0, weight_mass_cells = 10,
    structural_valid = TRUE, expression_supported = FALSE)
  temporal <- list(scores = scores,
    intermediate = list(bin_bounds = rep(list(seq(0, 1, by = 0.25)), 2L)))
  paths <- list(c(1L, 2L), c(1L, 3L))
  out <- AddLateForkPersistence(events, temporal, paths)
  positive <- AddLateForkPersistence(events, temporal, paths,
    expression_mean_mode = "positive_cells")
  expect_equal(out$late_contrast, 0)
  expect_equal(out$n_late_bins, 2L)
  expect_equal(out$persistent_positive_score, 0)
  expect_equal(positive$persistent_positive_score, 0)
  expect_equal(out$sibling_contrast, 0.4)
  expect_equal(out$late_change_contrast, 0)
  expect_equal(out$phase_pattern, "early_positive_only")
})

test_that("late expression loss does not shorten the observed trajectory", {
  events <- data.frame(pair_idx = 1L, source = 1L, target = 2L,
    event_node = 1L, child = 2L, fork_time = 0.25,
    sibling_children = "3", evaluable = TRUE, sibling_contrast = 0.4,
    focal_pre = 0.1, other_pre = 0.1, post_fork_change = 0.4)
  attr(events, "event_units") <- list(sibling_selection = "all")
  scores <- data.frame(pair_idx = 1L, source = 1L, target = 2L,
    lineage = rep(1:2, each = 4L), bin = rep(1:4, 2L),
    score = c(0.6, 0.6, 0.6, 0, rep(0.1, 4L)),
    weight_mass_cells = 10, structural_valid = TRUE,
    expression_supported = c(TRUE, TRUE, TRUE, FALSE, rep(TRUE, 4L)))
  temporal <- list(scores = scores,
    intermediate = list(bin_bounds = rep(list(seq(0, 1, by = 0.25)), 2L)))
  out <- AddLateForkPersistence(events, temporal,
    list(c(1L, 2L), c(1L, 3L)))
  expect_equal(out$late_contrast, 0.2)
  expect_equal(out$n_late_bins, 2L)
  expect_equal(out$persistent_positive_score, 0.2)
})

test_that("late persistence keeps an unobserved sister distinct from zero", {
  events <- data.frame(pair_idx = 1L, source = 1L, target = 2L,
    event_node = 1L, child = 2L, fork_time = 0.25,
    sibling_children = "3", evaluable = TRUE, sibling_contrast = 0.4,
    focal_pre = 0.1, other_pre = 0.1, post_fork_change = 0.4)
  attr(events, "event_units") <- list(sibling_selection = "all")
  scores <- data.frame(pair_idx = 1L, source = 1L, target = 2L,
    lineage = rep(1:2, each = 4L), bin = rep(1:4, 2L),
    score = c(rep(0.6, 4L), rep(NA_real_, 4L)),
    weight_mass_cells = c(rep(10, 4L), rep(0, 4L)),
    structural_valid = rep(c(TRUE, FALSE), each = 4L),
    expression_supported = rep(c(TRUE, FALSE), each = 4L))
  temporal <- list(scores = scores,
    intermediate = list(bin_bounds = rep(list(seq(0, 1, by = 0.25)), 2L)))
  out <- AddLateForkPersistence(events, temporal,
    list(c(1L, 2L), c(1L, 3L)))
  expect_true(is.na(out$late_contrast))
  expect_true(is.na(out$n_late_bins))
  expect_equal(out$persistent_positive_score, 0.4)
  expect_equal(out$phase_pattern, "late_unobserved")
  expect_true(is.na(out$late_change_contrast))
})

test_that("phase summaries recover late-only changes on the native sister panel", {
  events <- data.frame(pair_idx = 1L, source = 1L, target = 2L,
    event_node = 1L, child = 2L, fork_time = 0.25,
    sibling_children = "3,4", evaluable = TRUE, sibling_contrast = -0.2,
    focal_pre = 0.2, other_pre = 0.1, post_fork_change = -0.2)
  attr(events, "event_units") <- list(sibling_selection = "all")
  scores <- data.frame(pair_idx = 1L, source = 1L, target = 2L,
    lineage = rep(1:4, each = 4L), bin = rep(1:4, 4L),
    score = rep(c(0.8, 0.1, 0.3, 10), each = 4L),
    weight_mass_cells = rep(c(10, 10, 100, 1000), each = 4L),
    structural_valid = TRUE, expression_supported = TRUE)
  temporal <- list(scores = scores,
    intermediate = list(bin_bounds = rep(list(seq(0, 1, by = 0.25)), 4L)))
  out <- AddLateForkPersistence(events, temporal,
    list(c(1L, 2L), c(1L, 3L), c(1L, 4L), c(1L, 5L)))
  expect_equal(out$late_focal_activity_native, 0.8)
  expect_equal(out$late_other_activity_native, 0.2)
  expect_equal(out$late_change_contrast, 0.5)
  expect_equal(out$phase_pattern, "late_positive_only")
  expect_equal(out$n_late_phase_bins, 2L)
  expect_equal(out$persistent_positive_score, 0)
  positive_legacy <- AddLateForkPersistence(events, temporal,
    list(c(1L, 2L), c(1L, 3L), c(1L, 4L), c(1L, 5L)),
    persistence_mode = "positive")
  expect_equal(positive_legacy$persistent_positive_score, -0.2)
  expect_equal(out$late_contrast,0.8-(0.1+0.3)/2)
  mass_control <- AddLateForkPersistence(events,temporal,
    list(c(1L,2L),c(1L,3L),c(1L,4L),c(1L,5L)),sibling_scope="pooled_native")
  expect_equal(mass_control$late_contrast,0.8-(0.1*10+0.3*100)/110)
})

test_that("branch dynamics retains a pulse beyond an early-ending sister", {
  events <- data.frame(pair_idx = 1:3, source = 1L, target = 2L,
    event_node = 1L, child = 2L, sibling_children = "3,4",
    evaluable = TRUE, sibling_contrast = 0.2)
  scores <- data.table::CJ(pair_idx = 1:3, lineage = 1:3, bin = 1:9)
  scores[, `:=`(source = 1L, target = 2L, time = (bin - 1) / 8)]
  scores[, score := 0.2 + 0.2 * sin(2 * pi * time) + 0.03 * lineage +
    (lineage == 1L) * data.table::fcase(pair_idx == 1L,
      pmin(time, 0.5) - pmax(time - 0.5, 0),
      pair_idx == 2L, 0.4 * pmin(time, 0.5), default = 0)]
  scores[, `:=`(weight_mass_cells = 10, structural_valid = lineage != 3L | bin <= 3L)]
  temporal <- list(scores = scores,
    intermediate = list(bin_bounds = rep(list(seq(0, 1, length.out = 10L)), 3L)))
  paths <- list(c(1L, 2L), c(1L, 3L), c(1L, 4L))
  out <- SummarizeForkDynamics(events, temporal, paths)
  expect_equal(out$dynamic_pattern, c("pulse", "activation_plateau", "constant"))
  expect_equal(out$dynamic_n_bins, rep(9, 3L))
  expect_gt(out$pulse_bic_gain[1L], 0)
  expect_equal(out$pulse_bic_gain[2:3], c(0, 0))
  expect_equal(out$dynamic_change_time[1:2], c(0.5, 0.5))
  expect_equal(out$sibling_contrast, events$sibling_contrast)
  single <- SummarizeForkDynamics(events[1L, ], temporal, paths)
  expect_equal(single$dynamic_pattern, out$dynamic_pattern[1L])
  expect_equal(single$pulse_bic_gain, out$pulse_bic_gain[1L], tolerance = 1e-8)
  temporal$scores <- scores[bin <= 4L]
  short <- SummarizeForkDynamics(events, temporal, paths)
  expect_equal(short$pulse_bic_gain, rep(0, 3L))
  expect_equal(short$dynamic_status, rep("observed", 3L))
  expect_equal(short$dynamic_pattern, c("increasing", "increasing", "constant"))
  expect_equal(short$dynamic_n_bins, rep(4, 3L))
  expect_equal(short$dynamic_time_end, rep(7 / 18, 3L))
  temporal$scores <- scores[bin == 1L]
  one_time <- SummarizeForkDynamics(events, temporal, paths)
  expect_equal(one_time$dynamic_status, rep("insufficient_observed_overlap", 3L))
  expect_true(all(is.na(one_time$pulse_bic_gain)))
})

test_that("focal activity timing remains distinct from a later relative peak", {
  events <- data.frame(pair_idx = 1L, source = 1L, target = 2L,
    event_node = 1L, child = 2L, sibling_children = "3",
    evaluable = TRUE, sibling_contrast = 0.2, event_ranking_score = 0.1)
  scores <- data.table::CJ(lineage = 1:2, bin = 1:10)
  scores[, `:=`(pair_idx = 1L, source = 1L, target = 2L, time = (bin - 0.5) / 10)]
  scores[, score := 0.5 + 0.3 * pmin(time, 0.35) - 0.3 * pmax(time - 0.35, 0) -
    (lineage == 2L) * (0.04 + 0.25 * pmin(time, 0.65) - 0.3 * pmax(time - 0.65, 0))]
  scores[, `:=`(weight_mass_cells = 10, structural_valid = TRUE)]
  temporal <- list(scores = scores,
    intermediate = list(bin_bounds = rep(list(seq(0, 1, by = 0.1)), 2L)))
  paths <- list(c(1L, 2L), c(1L, 3L))
  out <- SummarizeForkDynamics(events, temporal, paths)
  expect_equal(out$dynamic_pattern, "pulse")
  expect_equal(out$dynamic_change_time, 0.65)
  expect_equal(out$focal_dynamic_pattern, "pulse")
  expect_equal(out$focal_peak_time, 0.35)
  expect_equal(out$focal_dynamic_n_bins, 10)
  expect_equal(out$event_ranking_score, events$event_ranking_score)
  tabular <- SummarizeForkDynamics(data.table::as.data.table(events), temporal, paths)
  expect_equal(tabular$focal_peak_time, out$focal_peak_time)
  expect_equal(tabular$dynamic_change_time, out$dynamic_change_time)
  temporal$scores[lineage == 2L & bin > 4L, structural_valid := FALSE]
  short <- SummarizeForkDynamics(events, temporal, paths)
  expect_equal(short$dynamic_status, "observed")
  expect_equal(short$dynamic_n_bins, 4)
  expect_equal(short$dynamic_time_end, 0.35)
  expect_equal(short$pulse_bic_gain, 0)
  expect_equal(short$focal_dynamic_status, "observed")
  expect_equal(short$focal_peak_time, 0.35)
})

test_that("fixed sisters give the same persistent contrast under unequal cell masses", {
  fc <- data.frame(evaluable=TRUE,sibling_contrast=.1,post_fork_change=.3,
    pair_idx=1L,source=1L,target=2L,event_node=1L,child=2L,fork_time=.1,
    sibling_children="3,4",unit_id=1L,focal_pre=0,other_pre=0)
  attr(fc,"event_units") <- list(activity_transform="log_exposure",sibling_selection="all")
  paths <- list(c(1L,2L),c(1L,3L),c(1L,4L))
  timing <- list(scores=data.frame(pair_idx=1L,source=1L,target=2L,
    lineage=rep(1:3,each=2),bin=rep(1:2,3),score=rep(-expm1(-c(.3,0,.4)),each=2),
    weight_mass_cells=rep(c(5,1,9),each=2),structural_valid=TRUE,expression_supported=TRUE),
    intermediate=list(bin_bounds=rep(list(c(0,.5,1)),3L)))
  expect_equal(AddLateForkPersistence(fc,timing,paths)$persistent_positive_score,.1,tolerance=1e-9)
  timing$scores$weight_mass_cells <- rep(c(5,9,1),each=2)
  expect_equal(AddLateForkPersistence(fc,timing,paths)$persistent_positive_score,.1,tolerance=1e-9)
})

test_that("within-sample event amplitudes remove a pure sample-mixture signal", {
  groups <- expand.grid(state=1:3,sample=c("A","B"),type=1:2,stringsAsFactors=FALSE)
  counts <- c(20,70,30,20,30,70)
  cells <- groups[rep(seq_len(nrow(groups)),rep(counts,2)),]
  expression <- Matrix::Matrix(rbind(rep(1,nrow(cells)),ifelse(cells$sample=="A",9,1)),sparse=TRUE)
  celltype <- as.integer(cells$type)
  paths <- list(c(1L,2L),c(1L,3L))
  weights <- cbind(ifelse(cells$state==1,.5,as.numeric(cells$state==2)),
    ifelse(cells$state==1,.5,as.numeric(cells$state==3)))
  time <- ifelse(cells$state==1,.4,.8)
  grid <- rep(list(c(0,.333,.666,1)),2)
  db <- list(ligand_rows=list(1L),receptor_rows=list(2L))
  pooled <- ComputeForkContrastsEvent(expression,celltype,db,time,weights,lineage_paths=paths,
    lineage_clusters=cells$state,bin_bounds=grid,event_scope="all",boundary="tree_segment_interval",
    scoring_mode="interval",hill_scale=1,activity_transform="log_exposure",activity_mode="continuous")
  within <- lapply(c("A","B"),function(id) {
    i <- which(cells$sample==id)
    local <- ComputeForkContrastsEvent(expression[,i],celltype[i],db,time[i],weights[i,],
      lineage_paths=paths,lineage_clusters=cells$state[i],units_override=attr(pooled,"event_units"),hill_scale=1)
    curve <- ComputeTemporalLRFixedBounds(expression[,i],celltype[i],db,time[i],weights[i,],
      bin_bounds=grid,hill_scale=1,return_complete_grid=TRUE)
    local <- AddLateForkPersistence(local,list(scores=curve,intermediate=list(bin_bounds=grid)),paths)
    local$block_label <- id
    local
  })
  observations <- data.table::rbindlist(within)
  score <- LineaLinkR:::.paired_event_amplitudes(pooled,observations,c("A","B"),measurement_error=FALSE)
  k <- pooled$source==1L & pooled$target==2L & pooled$child==2L
  expect_equal(pooled$sibling_contrast[k],log1p(sqrt(6.6))-log1p(sqrt(3.4)),tolerance=1e-6)
  expect_equal(score$sibling_contrast[k],0,tolerance=1e-9)
  expect_equal(score$persistent_positive_score[k],0,tolerance=1e-9)
  expect_equal(score$n_paired_samples[k],2L)
  planned <- pooled
  planned$sibling_children[k] <- paste0(planned$sibling_children[k],",99")
  partial <- LineaLinkR:::.paired_event_amplitudes(planned,observations,c("A","B"),measurement_error=FALSE)
  expect_equal(partial$n_paired_samples[k],2L)
  expect_equal(partial$n_full_panel_samples[k],0L)
  expect_equal(partial$paired_sibling_children[k],pooled$sibling_children[k])
  observations$evaluable <- FALSE
  empty <- LineaLinkR:::.paired_event_amplitudes(pooled,observations,c("A","B"),measurement_error=FALSE)
  expect_true(all(empty$n_paired_samples==0L & !empty$evaluable & is.na(empty$sibling_contrast)))
})
