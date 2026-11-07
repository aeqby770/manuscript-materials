test_that("unified tree projection follows the rooted segment formula", {
  embedding <- matrix(c(0, 0.2, 0.5, 1), ncol = 1L)
  result <- ComputeUnifiedTreeProjectionTime(
    embedding = embedding,
    cluster_labels = c(1L, 1L, 2L, 2L),
    lineage_paths = list(c(1L, 2L)),
    start_cell = 1L
  )

  expect_equal(result$pseudotime, c(0, 0.2 / 0.75, 0.5 / 0.75, 1))
  expect_equal(result$root_cluster, 1L)
  expect_equal(result$longest_root_to_leaf, 0.75)
})

test_that("cells cannot jump to a non-incident sister edge", {
  embedding <- rbind(
    c(0, 0), c(1, 0), c(2, 1), c(2, -1), c(1.9, -0.9)
  )
  result <- ComputeUnifiedTreeProjectionTime(
    embedding = embedding,
    cluster_labels = c(1L, 2L, 3L, 4L, 3L),
    lineage_paths = list(c(1L, 2L, 3L), c(1L, 2L, 4L)),
    start_cell = 1L
  )

  expect_true(all(is.finite(result$pseudotime)))
  expect_true(all(result$pseudotime >= 0 & result$pseudotime <= 1))
  expect_equal(result$projection$child[5L], 3L)
  expect_equal(nrow(result$tree_edges), 3L)
})

test_that("dominant-lineage constraint changes a mixed cluster assignment", {
  embedding <- rbind(
    c(0, 0), c(1, 0), c(2, 1), c(2, -1), c(0.2, 0)
  )
  incident <- ComputeUnifiedTreeProjectionTime(
    embedding = embedding,
    cluster_labels = c(1L, 2L, 3L, 4L, 3L),
    lineage_paths = list(c(1L, 2L, 3L), c(1L, 2L, 4L)),
    start_cell = 1L
  )
  dominant <- ComputeUnifiedTreeProjectionTime(
    embedding = embedding,
    cluster_labels = c(1L, 2L, 3L, 4L, 3L),
    lineage_paths = list(c(1L, 2L, 3L), c(1L, 2L, 4L)),
    start_cell = 1L,
    lineage_weights = cbind(c(1, 1, 1, 0, 1), c(1, 1, 0, 1, 0)),
    edge_constraint = "dominant_lineage_path"
  )

  expect_equal(incident$projection$parent[5L], 2L)
  expect_equal(dominant$projection$parent[5L], 1L)
  expect_lt(dominant$pseudotime[5L], incident$pseudotime[5L])
})
