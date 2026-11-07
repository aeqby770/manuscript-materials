test_that("InferStartCell does not use palantir pseudotime column", {
  set.seed(1)
  n <- 80
  emb <- cbind(
    seq(0, 1, length.out = n) + rnorm(n, sd = 0.05),
    rnorm(n)
  )
  cl <- rep(1:4, each = n / 4)
  # The adversarial Palantir label marks the last cell as earliest.
  pal_pt <- seq(1, 0, length.out = n)
  oracle_root <- which.min(pal_pt)

  blind <- InferStartCell(emb, cl, method = "mst_extremum", seed = 1L)
  expect_type(blind, "integer")
  expect_true(blind >= 1L && blind <= n)
  expect_false(blind == oracle_root)
})

test_that("ComputePseudotimeFuzzy rank_normalize=FALSE changes scale not order", {
  knn <- local({
    set.seed(42)
    n_cells <- 40
    k <- 5
    idx  <- matrix(sample.int(n_cells, n_cells * k, replace = TRUE),
                   nrow = n_cells, ncol = k)
    dist <- matrix(abs(rnorm(n_cells * k, 0.5, 0.2)), nrow = n_cells, ncol = k)
    list(idx = idx, dist = dist, n_cells = n_cells, k = k)
  })
  pt_rank <- ComputePseudotimeFuzzy(knn$idx, knn$dist, start_cell = 1L, rank_normalize = TRUE)
  pt_lin  <- ComputePseudotimeFuzzy(knn$idx, knn$dist, start_cell = 1L, rank_normalize = FALSE)
  v <- is.finite(pt_rank) & is.finite(pt_lin)
  expect_true(sum(v) > 10)
  expect_equal(cor(pt_rank[v], pt_lin[v], method = "spearman"), 1, tolerance = 1e-6)
  expect_false(isTRUE(all.equal(pt_rank[v], pt_lin[v])))
})
