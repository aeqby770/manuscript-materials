# Spatial adjacency constraints for ligand-receptor pairs.

make_spatial_coords <- function() {
  # Type 1 left, type 2 right; type 3 far right (isolated from type 1).
  rbind(
    c(0, 0),
    c(1, 0),
    c(10, 0),
    c(11, 0),
    c(20, 0)
  )
}

test_that("ComputeSpatialAdjacency kNN detects nearby types", {
  coords <- make_spatial_coords()
  types <- c(1L, 1L, 2L, 2L, 3L)
  adj <- ComputeSpatialAdjacency(types, coords, method = "knn", k = 2L)

  expect_equal(dim(adj), c(3L, 3L))
  expect_equal(rownames(adj), colnames(adj))
  expect_gt(adj["1", "2"], 0)
  expect_equal(adj["1", "3"], 0)
  expect_equal(attr(adj, "method"), "knn")
  expect_equal(attr(adj, "k"), 2L)
})

test_that("ComputeSpatialAdjacency radius method requires radius", {
  coords <- make_spatial_coords()
  types <- c(1L, 1L, 2L, 2L, 3L)
  expect_error(
    ComputeSpatialAdjacency(types, coords, method = "radius"),
    "radius"
  )
  adj <- ComputeSpatialAdjacency(
    types, coords, method = "radius", radius = 10
  )
  expect_equal(attr(adj, "method"), "radius")
  expect_gt(adj["1", "2"], 0)
  expect_equal(adj["1", "3"], 0)
})

test_that("ComputeSpatialAdjacency binary mode is 0/1", {
  coords <- make_spatial_coords()
  types <- c(1L, 1L, 2L, 2L, 3L)
  adj <- ComputeSpatialAdjacency(
    types, coords, method = "knn", k = 2L, binary = TRUE
  )
  expect_true(all(adj %in% c(0, 1)))
  expect_equal(adj["1", "2"], 1)
  expect_equal(adj["1", "3"], 0)
})

test_that("FilterLRBySpatialSupport flags and filters rows", {
  coords <- make_spatial_coords()
  types <- c(1L, 1L, 2L, 2L, 3L)
  adj <- ComputeSpatialAdjacency(types, coords, method = "knn", k = 2L)
  df <- data.frame(
    source = c(1L, 1L, 1L),
    target = c(2L, 3L, 2L),
    score = c(1, 2, 3),
    stringsAsFactors = FALSE
  )
  flagged <- FilterLRBySpatialSupport(df, adj, min_spatial_support = 0.05)
  expect_equal(nrow(flagged), 3L)
  expect_true(flagged$spatially_supported[1])
  expect_false(flagged$spatially_supported[2])
  expect_true(flagged$spatially_supported[3])
  expect_equal(flagged$spatial_support[1], adj["1", "2"])

  filtered <- FilterLRBySpatialSupport(
    df, adj, min_spatial_support = 0.05, filter = TRUE
  )
  expect_equal(nrow(filtered), 2L)
  expect_true(all(filtered$spatially_supported))
})

test_that("MergeEvidence keeps spatial adjacency as orthogonal annotation", {
  coords <- make_spatial_coords()
  types <- c(1L, 1L, 2L, 2L, 3L)
  adj <- ComputeSpatialAdjacency(types, coords, method = "knn", k = 2L)
  df <- data.frame(
    pair_idx = 1:2,
    source = c(1L, 1L),
    target = c(2L, 3L),
    lineage = c(1L, 1L),
    n_eff_cells = c(10, 10),
    lineage_entropy = c(0.2, 0.2),
    stringsAsFactors = FALSE
  )

  base <- MergeEvidence(df)
  expect_false(any(base$prioritized_candidate %in% TRUE))
  expect_false("spatially_eligible" %in% colnames(base))

  fake_cal <- df[, c("pair_idx", "source", "target", "lineage"), drop = FALSE]
  fake_cal$p_fake <- c(0.01, 0.01)
  fake_cal$q_fake <- c(0.02, 0.02)
  fake_cal$trajectory_evaluable <- TRUE
  spatial <- MergeEvidence(
    df,
    fake_calibration = fake_cal,
    spatial_adjacency = adj,
    min_spatial_support = 0.05
  )
  expect_true(spatial$spatially_eligible[1])
  expect_false(spatial$spatially_eligible[2])
  expect_true(spatial$prioritized_candidate[1])
  expect_true(spatial$prioritized_candidate[2])
  expect_true(spatial$identity_prioritized[1])
  expect_true(spatial$identity_prioritized[2])
  expect_equal(spatial$spatial_support[1], adj["1", "2"])
})

test_that("ComputeTemporalLR attaches spatial columns when adjacency given", {
  skip_if_not_installed("Matrix")
  set.seed(1)
  n_genes <- 12L
  n_cells <- 30L
  x <- Matrix::rsparsematrix(n_genes, n_cells, density = 0.3)
  x@x <- abs(x@x)
  coords <- cbind(runif(n_cells), runif(n_cells))
  ct <- sample.int(2L, n_cells, replace = TRUE)
  pt <- seq(0, 1, length.out = n_cells)
  lw <- matrix(1, n_cells, 1L)
  lr_db <- list(ligand_rows = list(1L), receptor_rows = list(2L))
  adj <- ComputeSpatialAdjacency(ct, coords, method = "knn", k = 3L)
  scores <- ComputeTemporalLR(
    x, ct, lr_db, pt, lw, n_bins = 4L, min_pct = 0,
    spatial_adjacency = adj, min_spatial_support = 0.05
  )
  skip_if(nrow(scores) == 0L, "no scores")
  expect_true(all(c("spatial_support", "spatially_eligible") %in% colnames(scores)))
  expect_true(any(is.finite(scores$spatial_support)))
})
