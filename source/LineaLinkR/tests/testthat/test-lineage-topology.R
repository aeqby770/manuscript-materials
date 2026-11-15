test_that("ClusterLineageTopology is deterministic and compact", {
  set.seed(12345)
  embedding <- rbind(
    matrix(rnorm(80, -3, 0.15), ncol = 2),
    matrix(rnorm(80, -1, 0.15), ncol = 2),
    matrix(rnorm(80, 1, 0.15), ncol = 2),
    matrix(rnorm(80, 3, 0.15), ncol = 2)
  )

  first <- ClusterLineageTopology(
    embedding,
    k = 4L,
    nstart = 1L,
    iter_max = 500L,
    seed = 2L
  )
  second <- ClusterLineageTopology(
    embedding,
    k = 4L,
    nstart = 1L,
    iter_max = 500L,
    seed = 2L
  )

  expect_identical(first$labels, second$labels)
  expect_identical(sort(unique(first$labels)), seq_len(first$occupied_k))
  expect_identical(first$requested_k, 4L)
  expect_identical(first$occupied_k, 4L)
  expect_length(first$labels, nrow(embedding))
  expect_equal(dim(first$centers), c(4L, 2L))
  expect_identical(first$empty_start_warnings, 0L)
})

test_that("ClusterLineageTopology does not accept truth inputs", {
  expect_false(any(grepl(
    "truth|time|stage|day|branch",
    names(formals(ClusterLineageTopology)),
    ignore.case = TRUE
  )))
})

test_that("state geometry preserves directed observations, components and explicit roots", {
  embedding <- cbind(position = rep(0:5, 2L),
    compartment = rep(c(0, 10), each = 6L))
  connectivity <- Matrix::sparseMatrix(
    i = c(1:5, 2:6, 7:11, 8:12),
    j = c(2:6, 1:5, 8:12, 7:11),
    x = rep(c(2, 1, 2, 1), each = 5L),
    dims = c(12L, 12L)) + Matrix::Diagonal(12L)
  result <- RefineStateGraph(embedding, connectivity, start_cell = 3L,
    k = 4L, n_neighbors = 3L, project_root = FALSE)
  expected <- as.matrix(Matrix::Diagonal(
    x = 1 / Matrix::rowSums(connectivity)) %*% connectivity %*% embedding)
  component <- rep(1:2, each = 6L)
  observed <- which(result$nn_idx > 0L, arr.ind = TRUE)
  expect_equal(result$embedding, expected, tolerance = 1e-12)
  expect_equal(dim(result$embedding), dim(embedding))
  expect_identical(result$start_cell, 3L)
  expect_equal(component[result$nn_idx[observed]], component[observed[, 1L]])
  expect_true(all(result$nn_idx[observed] != observed[, 1L]))
  expect_equal(rowSums(result$lineage_weights), rep(1, nrow(embedding)),
    tolerance = 1e-6)
})


test_that("conditional composition retains cell differences within a lineage", {
  outcomes <- diag(2)
  colnames(outcomes) <- c("A","B")
  query <- matrix(.5,nrow=2,ncol=2)
  neighbors <- rbind(c(1L,2L),c(2L,1L))
  weights <- rbind(c(.9,.1),c(.9,.1))
  result <- PredictLineageComposition(query,diag(2),outcomes,neighbors,weights,TRUE)
  expect_equal(rowSums(result$prediction),c(1,1),tolerance=1e-12)
  expect_true(all(result$prediction>0))
  expect_gt(result$prediction[1,"A"],result$prediction[2,"A"])
  permuted <- PredictLineageComposition(query,diag(2)[2:1,,drop=FALSE],
    outcomes[2:1,,drop=FALSE],3L-neighbors,weights)
  expect_equal(permuted,result$prediction,tolerance=1e-12)
})
