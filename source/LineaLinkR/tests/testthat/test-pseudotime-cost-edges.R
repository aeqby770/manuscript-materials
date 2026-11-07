test_that("realized pseudotime costs expose one symmetrized edge per pair", {
  idx <- matrix(c(
    2L, 3L,
    1L, 3L,
    1L, 2L
  ), nrow = 3L, byrow = TRUE)
  dists <- matrix(c(
    1, 3,
    1, 2,
    3, 2
  ), nrow = 3L, byrow = TRUE)

  edges <- LineaLinkR:::ComputePseudotimeCostEdges(idx, dists)

  expect_equal(nrow(edges), 3L)
  expect_equal(paste(edges$from, edges$to, sep = ">"), c("1>2", "1>3", "2>3"))
  expect_equal(edges$raw_distance, c(1, 3, 2))
  expect_true(all(edges$fuzzy_cost >= edges$raw_distance))
})
