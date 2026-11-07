# Spatial permutation null validation.

make_spatial_null_coords <- function() {
  rbind(
    c(0, 0),
    c(1, 0),
    c(10, 0),
    c(11, 0),
    c(20, 0)
  )
}

test_that("ComputeSpatialAdjacencyNull returns null summary fields", {
  coords <- make_spatial_null_coords()
  types <- c(1L, 1L, 2L, 2L, 3L)
  null <- ComputeSpatialAdjacencyNull(
    types,
    coords,
    method = "knn",
    k = 2L,
    null_mode = "label_permutation",
    n_perm = 49L,
    seed = 42L
  )
  expect_equal(null$observed, ComputeSpatialAdjacency(types, coords, k = 2L))
  expect_equal(dim(null$null_mean), dim(null$observed))
  expect_equal(dim(null$null_sd), dim(null$observed))
  expect_equal(dim(null$p_spatial), dim(null$observed))
  expect_equal(dim(null$q_spatial), dim(null$observed))
  expect_equal(null$null_percentile, null$p_spatial)
  expect_equal(null$null_mode, "label_permutation")
  expect_equal(null$n_perm, 49L)
  expect_true(all(null$null_percentile >= 0 & null$null_percentile <= 1))
})

test_that("NormalizeSpatialAdjacency computes z-scores from null_result", {
  coords <- make_spatial_null_coords()
  types <- c(1L, 1L, 2L, 2L, 3L)
  null <- ComputeSpatialAdjacencyNull(
    types, coords, method = "knn", k = 2L,
    null_mode = "label_permutation", n_perm = 49L, seed = 7L
  )
  norm <- NormalizeSpatialAdjacency(null_result = null)
  expect_equal(norm$adjacency, null$observed)
  expect_equal(dim(norm$z_score), dim(null$observed))
  expect_equal(norm$null_percentile, null$null_percentile)
  expect_equal(norm$p_spatial, null$p_spatial)
  expect_equal(norm$q_spatial, null$q_spatial)
  expect_gt(norm$z_score["1", "2"], norm$z_score["1", "3"])
})

test_that("FilterLRBySpatialSupport adds null columns in label_permutation mode", {
  coords <- make_spatial_null_coords()
  types <- c(1L, 1L, 2L, 2L, 3L)
  adj <- ComputeSpatialAdjacency(types, coords, method = "knn", k = 2L)
  df <- data.frame(
    source = c(1L, 1L),
    target = c(2L, 3L),
    score = c(1, 2),
    stringsAsFactors = FALSE
  )
  out <- FilterLRBySpatialSupport(
    df,
    adj,
    null_mode = "label_permutation",
    celltypes = types,
    spatial_coords = coords,
    n_perm = 49L,
    seed = 11L,
    max_spatial_q = 1
  )
  expect_true(all(c(
    "spatial_support_z", "p_spatial", "q_spatial",
    "spatial_support_null_pctl", "spatial_null_mode"
  ) %in% colnames(out)))
  expect_equal(out$spatial_null_mode, rep("label_permutation", 2L))
  expect_true(is.finite(out$spatial_support_z[1]))
  expect_lt(out$spatial_support_null_pctl[1], out$spatial_support_null_pctl[2])
  expect_true(all(out$q_spatial >= out$p_spatial))
})

test_that("FilterLRBySpatialSupport raw mode is backward compatible", {
  coords <- make_spatial_null_coords()
  types <- c(1L, 1L, 2L, 2L, 3L)
  adj <- ComputeSpatialAdjacency(types, coords, method = "knn", k = 2L)
  df <- data.frame(source = 1L, target = 2L, stringsAsFactors = FALSE)
  out <- FilterLRBySpatialSupport(df, adj, null_mode = "raw")
  expect_false(any(c(
    "spatial_support_z", "spatial_support_null_pctl", "spatial_null_mode"
  ) %in% colnames(out)))
})

test_that("FilterLRBySpatialSupport requires coords for null modes", {
  coords <- make_spatial_null_coords()
  types <- c(1L, 1L, 2L, 2L, 3L)
  adj <- ComputeSpatialAdjacency(types, coords, method = "knn", k = 2L)
  df <- data.frame(source = 1L, target = 2L, stringsAsFactors = FALSE)
  expect_error(
    FilterLRBySpatialSupport(
      df, adj, null_mode = "size_matched", celltypes = types
    ),
    "spatial_coords"
  )
})

test_that("block_label_permutation preserves labels within blocks", {
  coords <- rbind(
    c(0, 0), c(1, 0), c(2, 0),
    c(10, 0), c(11, 0), c(12, 0)
  )
  types <- c(1L, 2L, 3L, 1L, 2L, 3L)
  blocks <- rep(c("A", "B"), each = 3L)
  set.seed(99)
  perm <- LineaLinkR:::.permute_celltypes_within_blocks(types, blocks)
  expect_equal(sort(perm[blocks == "A"]), sort(types[blocks == "A"]))
  expect_equal(sort(perm[blocks == "B"]), sort(types[blocks == "B"]))
})

test_that("ComputeSpatialAdjacencyNull block mode differs from global permutation", {
  coords <- rbind(
    c(0, 0), c(1, 0), c(2, 0),
    c(10, 0), c(11, 0), c(12, 0)
  )
  types <- c(1L, 1L, 2L, 2L, 3L, 3L)
  blocks <- rep(c("D1", "D2"), each = 3L)
  global <- ComputeSpatialAdjacencyNull(
    types, coords, method = "knn", k = 2L,
    null_mode = "label_permutation", n_perm = 79L, seed = 3L
  )
  block <- ComputeSpatialAdjacencyNull(
    types, coords, method = "knn", k = 2L,
    null_mode = "block_label_permutation",
    block_labels = blocks, n_perm = 79L, seed = 3L
  )
  expect_equal(block$null_mode, "block_label_permutation")
  expect_false(isTRUE(all.equal(global$null_mean, block$null_mean)))
})

test_that("FilterLRBySpatialSupport block mode requires block_labels", {
  coords <- make_spatial_null_coords()
  types <- c(1L, 1L, 2L, 2L, 3L)
  adj <- ComputeSpatialAdjacency(types, coords, method = "knn", k = 2L)
  df <- data.frame(source = 1L, target = 2L, stringsAsFactors = FALSE)
  expect_error(
    FilterLRBySpatialSupport(
      df, adj,
      null_mode = "block_label_permutation",
      celltypes = types, spatial_coords = coords
    ),
    "block_labels"
  )
})
