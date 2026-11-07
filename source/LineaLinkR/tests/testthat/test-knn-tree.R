# Variant A: cluster tree from the native cell kNN graph, and consistency of
# paths / weights / fork placement derived from it.

make_knn_tree_fixture <- function() {
  # Three clusters in a chain: 1-2-3 plus an isolated cluster 4.
  # Cluster 1: cells 1-6, cluster 2: cells 7-12, cluster 3: cells 13-18,
  # cluster 4: cells 19-24 (no observed edge to the others).
  set.seed(42)
  n <- 24
  emb <- cbind(
    c(rnorm(6, 0), rnorm(6, 1), rnorm(6, 2), rnorm(6, 5)),
    c(rnorm(6, 0), rnorm(6, 0.2), rnorm(6, 0.4), rnorm(6, 5))
  )
  cl <- rep(1:4, each = 6)
  # kNN connectivity: strong 1-2 and 2-3 mass, dense inside clusters,
  # nothing connecting cluster 4.
  conn <- Matrix::Matrix(0, n, n, sparse = TRUE)
  add_edge <- function(i, j, w) {
    conn[i, j] <<- w
    conn[j, i] <<- w
  }
  within <- list(1:6, 7:12, 13:18, 19:24)
  for (g in within) {
    for (i in g) for (j in g) if (i < j) add_edge(i, j, 1)
  }
  for (i in 1:6) for (j in 7:12) add_edge(i, j, 2)
  for (i in 7:12) for (j in 13:18) add_edge(i, j, 2)
  list(emb = emb, cl = cl, conn = conn)
}

test_that("cluster tree uses observed support and keeps components", {
  f <- make_knn_tree_fixture()
  tree <- BuildKnnClusterTree(f$conn, f$cl, 1L)
  # The fixture observes 1-2 and 2-3 only: the tree keeps those two edges.
  expect_equal(nrow(tree$edges), 2L)
  expect_equal(sort(tree$edges$from), c(1L, 2L))
  expect_equal(sort(tree$edges$to), c(2L, 3L))
  expect_true(all(tree$edges$support > 0))
  expect_equal(tree$n_components, 2L)         # cluster 4 stays isolated
  expect_equal(tree$root_component, tree$components[1L])
  # No invented edge to the isolated cluster.
  expect_false(any(tree$edges$from == 4L | tree$edges$to == 4L))
})

test_that("paths and weights come from the same tree", {
  f <- make_knn_tree_fixture()
  paths <- ExtractLineagesFromKnnGraph(f$conn, f$cl, 1L)
  # Chain 1-2-3 with root 1 -> exactly one root-to-leaf path.
  expect_equal(length(paths$paths), 1L)
  expect_equal(paths$paths[[1]], c(1L, 2L, 3L))

  w <- ComputeLineageWeightsFromKnnGraph(
    f$conn, f$emb, f$cl, 1L, sigma = 1.0, tree = paths$tree
  )
  expect_equal(ncol(w), 1L)
  expect_equal(nrow(w), length(f$cl))
  expect_true(all(abs(rowSums(w) - 1) < 1e-8))
  # Cluster 4 has no path: geometric fall-back keeps finite normalised rows.
  expect_true(all(is.finite(w)))
})

test_that("trunk clusters receive geometric weights, own-path clusters 1.0", {
  f <- make_knn_tree_fixture()
  # Two leaves from a shared root: 2 is a trunk cluster on both paths.
  cl <- c(rep(1, 6), rep(2, 6), rep(3, 6))
  conn <- Matrix::Matrix(0, 18, 18, sparse = TRUE)
  for (i in 1:6) for (j in 7:12) {
    conn[i, j] <- 2; conn[j, i] <- 2
  }
  for (i in 1:6) for (j in 13:18) {
    conn[i, j] <- 2; conn[j, i] <- 2
  }
  # Cluster 3 sits far from the other lineage's polyline so its off-path
  # geometric fall-back is negligible.
  emb <- cbind(c(rnorm(6, 0), rnorm(6, 1), rnorm(6, -3)),
               c(rnorm(6, 0), rnorm(6, 1), rnorm(6, -3)))
  paths <- ExtractLineagesFromKnnGraph(conn, cl, 1L)
  expect_equal(length(paths$paths), 2L)
  w <- ComputeLineageWeightsFromKnnGraph(conn, emb, cl, 1L, sigma = 1.0,
                                         tree = paths$tree)
  expect_equal(ncol(w), 2L)
  # Leaf cells on exactly one path dominate that lineage's weight.
  leaf3 <- which(cl == 3)
  expect_true(all(w[leaf3, 1] < 0.05))
  expect_true(all(w[leaf3, 2] > 0.95))
  # Root cluster cells (on both paths) are geometric and strictly between 0/1.
  root <- which(cl == 1)
  expect_true(all(w[root, ] > 0))
  expect_true(all(w[root, ] < 1))
})
