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


test_that("event segments preserve terminals and use observed graph boundaries", {
  paths <- list(c(1L, 2L, 3L, 4L, 5L), c(1L, 2L, 3L, 6L, 7L))
  graph <- Matrix::sparseMatrix(i = c(3L, 4L, 3L, 6L),
    j = c(4L, 3L, 6L, 3L), x = 1, dims = c(7L, 7L))
  value <- LineaLinkR:::.build_branch_segment_tree(paths, seq_len(7L),
    c(0, 0.1, 0.3, 0.4, 0.7, 0.5, 0.8), graph)
  expect_identical(value$lineage_paths, list(c(3L, 5L), c(3L, 7L)))
  expect_identical(value$lineage_clusters, c(3L, 3L, 3L, 5L, 5L, 7L, 7L))
  expect_equal(value$anchors$fork_time, 0.375)
  expect_equal(value$anchors$n_observed_children, 2L)
})


test_that("a single RNA path contracts without introducing a fork", {
  graph <- Matrix::sparseMatrix(i=c(1L,2L,2L,3L,3L,4L),
    j=c(2L,1L,3L,2L,4L,3L),x=1,dims=c(4L,4L))
  value <- LineaLinkR:::.build_branch_segment_tree(list(1:4),1:4,
    c(0,.3,.6,1),graph)
  expect_identical(value$lineage_paths,list(4L))
  expect_identical(value$lineage_clusters,rep(4L,4L))
  expect_equal(nrow(value$anchors),0L)
})
