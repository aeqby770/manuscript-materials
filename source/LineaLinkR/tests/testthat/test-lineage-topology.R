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
