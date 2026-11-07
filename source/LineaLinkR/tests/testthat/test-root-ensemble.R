test_that("attractor diagnostics distinguish nearest and mean distance", {
  f <- make_start_cell_fixture()
  out <- InferStartCellAttractor(
    f$embedding, f$clusters, f$nn_idx, f$nn_dists, f$potency,
    complete_terminal_coverage = TRUE
  )

  expect_true(all(c(
    "nearest_attractor_distance", "mean_attractor_distance"
  ) %in% colnames(out$scores)))
  expect_false(isTRUE(all.equal(
    out$scores$nearest_attractor_distance,
    out$scores$mean_attractor_distance
  )))
  expect_false("potential" %in% colnames(out$scores))
})

test_that("attractor score has five independent evidence weights", {
  f <- make_start_cell_fixture()
  expect_error(
    InferStartCellAttractor(
      f$embedding, f$clusters, f$nn_idx, f$nn_dists,
      f$potency, weights = rep(1, 6), complete_terminal_coverage = TRUE
    ),
    "length 5"
  )
})

test_that("root ensemble reports ambiguity for a symmetric path", {
  f <- make_symmetric_root_fixture()
  out <- InferRootEnsemble(
    f$embedding, f$clusters, f$nn_idx, f$nn_dists,
    fate_entropy = f$fate_entropy,
    commitment = f$commitment,
    potency = f$potency,
    k_grid = c(4L, 6L, 8L), pc_grid = c(1L, 2L),
    n_resamples = 12L, subsample_fraction = 0.8,
    min_direction_margin = 0.05, seed = 12345L
  )

  expect_equal(out$direction_status, "ambiguous")
  expect_lte(out$direction_margin, 0.05)
  expect_true(all(c(
    "candidate_region", "selection_frequency", "coverage",
    "rank_stability", "entropy_gradient", "commitment_gradient",
    "potency_gradient"
  ) %in% colnames(out$candidates)))
  expect_equal(out$direction_score_components,
               c("attractor_direction", "potency_gradient"))
  expect_equal(out$confidence_components, c("coverage", "rank_stability"))
  expect_equal(out$effective_direction_components, 2L)
})

test_that("root direction score excludes confidence and merges attractor axes", {
  perturbations <- data.frame(
    perturbation = rep(1L, 2L),
    candidate_region = 1:2,
    coverage = c(0.91, 1.00),
    rank_stability = c(0.10, 1.00),
    entropy_gradient = c(0.80, 0.20),
    commitment_gradient = c(0.80, 0.20),
    potency_gradient = c(NA_real_, NA_real_),
    successful = TRUE
  )

  scored <- LineaLinkR:::.root_score_perturbations(perturbations)

  expect_equal(scored$attractor_direction, c(0.90, 0.60))
  expect_equal(scored$direction_component_count, c(1L, 1L))
  expect_equal(scored$winner_share, c(1, 0))
  expect_equal(scored$aggregate_score, scored$attractor_direction)
  expect_equal(scored$confidence_score, c(0.505, 1.000))
})

test_that("ambiguous root fallback can regularize direction with confidence", {
  perturbations <- data.frame(
    perturbation = rep(1L, 2L),
    candidate_region = 1:2,
    coverage = c(0.91, 1.00),
    rank_stability = c(0.10, 1.00),
    entropy_gradient = c(0.80, 0.20),
    commitment_gradient = c(0.80, 0.20),
    potency_gradient = c(NA_real_, NA_real_),
    successful = TRUE
  )

  scored <- LineaLinkR:::.root_score_perturbations(
    perturbations,
    selection_score_mode = "confidence_regularized"
  )

  expect_equal(scored$winner_share, c(0, 1))
  expect_equal(scored$selection_score_mode, rep(
    "confidence_regularized", 2L
  ))
  expect_equal(scored$aggregate_score, c(0.7025, 0.8000))
})

test_that("cached maximum-k neighbours reproduce direct smaller-k graphs", {
  f <- make_directional_y_fixture()
  cached <- LineaLinkR:::.root_build_knn_grid(
    f$embedding[, 1:2, drop = FALSE], c(4L, 7L, 9L)
  )

  expect_equal(cached[["4"]], LineaLinkR:::.root_build_knn(
    f$embedding[, 1:2, drop = FALSE], 4L
  ))
  expect_equal(cached[["7"]], LineaLinkR:::.root_build_knn(
    f$embedding[, 1:2, drop = FALSE], 7L
  ))
  expect_equal(cached[["9"]], LineaLinkR:::.root_build_knn(
    f$embedding[, 1:2, drop = FALSE], 9L
  ))
})

test_that("batched Rust pseudotime equals repeated single-root calls", {
  f <- make_directional_y_fixture()
  starts <- c(1L, 20L, 40L)
  batched <- LineaLinkR:::ComputePseudotimeFuzzyMulti(
    f$nn_idx, f$nn_dists, starts, rank_normalize = TRUE
  )
  repeated <- vapply(starts, function(start_cell) {
    ComputePseudotimeFuzzy(
      f$nn_idx, f$nn_dists, start_cell, rank_normalize = TRUE
    )
  }, numeric(nrow(f$embedding)))
  single_batch <- LineaLinkR:::ComputePseudotimeFuzzyMulti(
    f$nn_idx, f$nn_dists, starts[[1L]], rank_normalize = TRUE
  )

  expect_equal(batched, repeated, tolerance = 0)
  expect_equal(dim(single_batch), c(nrow(f$embedding), 1L))
  expect_equal(single_batch[, 1L], repeated[, 1L], tolerance = 0)
})

test_that("native kNN component count detects connected and split graphs", {
  connected <- matrix(c(2L, 1L, 4L, 3L, 3L, 4L, 1L, 2L), nrow = 4L)
  split <- matrix(c(2L, 1L, 4L, 3L), nrow = 4L)

  expect_equal(LineaLinkR:::KnnComponentCount(connected), 1L)
  expect_equal(LineaLinkR:::KnnComponentCount(split), 2L)
  expect_equal(
    LineaLinkR:::KnnComponentLabels(split),
    c(1L, 1L, 2L, 2L)
  )
})

test_that("root ensemble selects the upstream trunk when evidence is directional", {
  f <- make_directional_y_fixture()
  out <- InferRootEnsemble(
    f$embedding, f$clusters, f$nn_idx, f$nn_dists,
    fate_entropy = f$fate_entropy,
    commitment = f$commitment,
    potency = f$potency,
    k_grid = c(5L, 7L, 9L), pc_grid = c(1L, 2L),
    n_resamples = 12L, subsample_fraction = 0.8,
    min_direction_margin = 0.05, seed = 12345L
  )

  expect_equal(out$direction_status, "supported")
  expect_true(out$start_cell %in% f$expected_root_cells)
  expect_gte(out$root_support, 0.7)
})

test_that("one blind attractor axis cannot support a signed direction claim", {
  f <- make_directional_y_fixture()
  out <- InferRootEnsemble(
    f$embedding, f$clusters, f$nn_idx, f$nn_dists,
    fate_entropy = f$fate_entropy,
    commitment = f$commitment,
    potency = NULL,
    k_grid = c(5L, 7L), pc_grid = c(1L, 2L),
    n_resamples = 6L, min_resamples = 3L,
    adaptive_resampling = TRUE,
    subsample_fraction = 0.8,
    min_direction_margin = 0.05, seed = 12345L
  )

  expect_equal(out$effective_direction_components, 1L)
  expect_equal(out$direction_status, "ambiguous")
})

test_that("adaptive root ensemble stops after stable minimum resamples", {
  f <- make_directional_y_fixture()
  out <- InferRootEnsemble(
    f$embedding, f$clusters, f$nn_idx, f$nn_dists,
    fate_entropy = f$fate_entropy,
    commitment = f$commitment,
    potency = f$potency,
    k_grid = c(5L, 7L), pc_grid = c(1L, 2L),
    n_resamples = 12L, min_resamples = 3L,
    adaptive_resampling = TRUE,
    subsample_fraction = 0.8,
    min_direction_margin = 0.05, seed = 12345L
  )

  expect_equal(max(out$perturbations$resample), 3L)
  expect_true(out$adaptive_stopped)
  expect_equal(out$resamples_used, 3L)
})

test_that("adaptive root ensemble retains full audit for ambiguous direction", {
  f <- make_symmetric_root_fixture()
  out <- InferRootEnsemble(
    f$embedding, f$clusters, f$nn_idx, f$nn_dists,
    fate_entropy = f$fate_entropy,
    commitment = f$commitment,
    potency = f$potency,
    k_grid = c(4L, 6L), pc_grid = c(1L, 2L),
    n_resamples = 6L, min_resamples = 3L,
    adaptive_resampling = TRUE,
    subsample_fraction = 0.8,
    min_direction_margin = 0.05, seed = 12345L
  )

  expect_equal(max(out$perturbations$resample), 6L)
  expect_false(out$adaptive_stopped)
  expect_equal(out$resamples_used, 6L)
})

test_that("root ensemble is invariant to PC sign and cell-row permutation", {
  f <- make_directional_y_fixture()
  args <- list(
    fate_entropy = f$fate_entropy,
    commitment = f$commitment,
    potency = f$potency,
    k_grid = c(5L, 7L), pc_grid = c(1L, 2L),
    n_resamples = 8L, subsample_fraction = 0.8,
    min_direction_margin = 0.05, seed = 12345L
  )
  base <- do.call(InferRootEnsemble, c(list(
    f$embedding, f$clusters, f$nn_idx, f$nn_dists
  ), args))

  flipped_embedding <- f$embedding
  flipped_embedding[, 1L] <- -flipped_embedding[, 1L]
  flipped <- do.call(InferRootEnsemble, c(list(
    flipped_embedding, f$clusters, f$nn_idx, f$nn_dists
  ), args))

  permutation <- c(seq(2L, nrow(f$embedding), by = 2L),
                   seq(1L, nrow(f$embedding), by = 2L))
  permuted_embedding <- f$embedding[permutation, , drop = FALSE]
  permuted_knn <- make_knn_fixture(permuted_embedding, 8L)
  permuted_args <- args
  permuted_args$fate_entropy <- f$fate_entropy[permutation]
  permuted_args$commitment <- f$commitment[permutation]
  permuted_args$potency <- f$potency[permutation]
  permuted <- do.call(InferRootEnsemble, c(list(
    permuted_embedding, f$clusters[permutation],
    permuted_knn$nn_idx, permuted_knn$nn_dists
  ), permuted_args))

  expect_equal(f$clusters[base$start_cell], f$clusters[flipped$start_cell])
  expect_equal(f$clusters[base$start_cell],
               f$clusters[permutation[permuted$start_cell]])
  expect_equal(base$direction_status, flipped$direction_status)
  expect_equal(base$direction_status, permuted$direction_status)
})

test_that("root ensemble remains evaluable when subsample kNN components split", {
  set.seed(12345)
  embedding <- rbind(
    cbind(rnorm(20, 0, 0.05), rnorm(20, 0, 0.05)),
    cbind(rnorm(20, 10, 0.05), rnorm(20, 0, 0.05)),
    cbind(rnorm(20, 20, 0.05), rnorm(20, 0, 0.05))
  )
  clusters <- rep(1:3, each = 20)
  raw_knn <- make_knn_fixture(embedding, 5L)
  connected <- ConnectKnnComponents(
    raw_knn$nn_idx, raw_knn$nn_dists, embedding
  )

  out <- InferRootEnsemble(
    embedding, clusters, connected$nn_idx, connected$nn_dists,
    fate_entropy = seq(1, 0, length.out = nrow(embedding)),
    commitment = seq(0, 1, length.out = nrow(embedding)),
    k_grid = c(3L, 4L), pc_grid = c(1L, 2L),
    n_resamples = 4L, subsample_fraction = 0.8,
    min_direction_margin = 0.05, seed = 12345L
  )

  expect_true(is.finite(out$root_support))
  expect_true(is.finite(out$direction_margin))
  expect_true(out$direction_status %in% c("ambiguous", "weak", "supported"))
})

test_that("induced kNN components ignore zero-padded neighbours", {
  nn_idx <- matrix(
    c(
      2L, 0L, 4L,
      1L, 0L, 0L,
      4L, 0L, 0L,
      3L, 0L, 1L
    ),
    nrow = 4L,
    byrow = TRUE
  )

  expect_warning(
    components <- LineaLinkR:::.knn_components_induced(
      nn_idx, subset = 1:4, n_cells = 4L
    ),
    NA
  )
  expect_equal(components, list(c(1L, 2L, 4L, 3L)))
})
