make_toy_knn_per_lineage <- function(n_cells = 50, k = 5) {
  set.seed(42)
  idx <- matrix(sample.int(n_cells, n_cells * k, replace = TRUE), nrow = n_cells, ncol = k)
  dist <- matrix(abs(rnorm(n_cells * k, mean = 1, sd = 0.3)), nrow = n_cells, ncol = k)
  for (i in seq_len(n_cells)) {
    idx[i, idx[i, ] == i] <- ((i %% n_cells) + 1L)
  }
  list(idx = idx, dist = dist, n_cells = n_cells, k = k)
}

test_that("ComputePseudotimePerLineage returns matrix with root zero and inactive lineage NaN mask", {
  skip("ComputePseudotimePerLineage is not exported in this release")
  knn <- make_toy_knn_per_lineage(50, 5)
  lw <- cbind(
    rep(1, knn$n_cells),
    c(rep(1, knn$n_cells / 2), rep(0, knn$n_cells / 2))
  )
  pt_lin <- ComputePseudotimePerLineage(
    knn$idx,
    knn$dist,
    lineage_weights = lw,
    start_cell = 1L,
    threshold = 0.05
  )
  expect_true(is.matrix(pt_lin))
  expect_equal(nrow(pt_lin), knn$n_cells)
  expect_equal(ncol(pt_lin), 2L)
  expect_equal(pt_lin[1, 1], 0.0)
  expect_equal(pt_lin[1, 2], 0.0)
  expect_true(all(is.finite(pt_lin[1:(knn$n_cells / 2), 2]) | is.nan(pt_lin[1:(knn$n_cells / 2), 2])))
  expect_true(all(is.nan(pt_lin[(knn$n_cells / 2 + 1):knn$n_cells, 2])))
})
