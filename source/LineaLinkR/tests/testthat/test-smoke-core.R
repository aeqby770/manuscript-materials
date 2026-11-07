# Smoke tests for exported core API shapes and input validation.
# Numerical behavior is covered by the Rust unit tests.

# ── Helpers ──────────────────────────────────────────────────────────────
make_toy_knn <- function(n_cells = 50, k = 5) {
  set.seed(42)
  idx  <- matrix(sample.int(n_cells, n_cells * k, replace = TRUE),
                 nrow = n_cells, ncol = k)
  dist <- matrix(abs(rnorm(n_cells * k, mean = 1, sd = 0.3)),
                 nrow = n_cells, ncol = k)
  # Remove self-loops from the fixture.

  for (i in seq_len(n_cells)) {
    idx[i, idx[i, ] == i] <- ((i %% n_cells) + 1L)
  }
  list(idx = idx, dist = dist, n_cells = n_cells, k = k)
}

make_toy_embedding <- function(n_cells = 50, n_dims = 5, n_clusters = 3) {
  set.seed(42)
  emb <- matrix(rnorm(n_cells * n_dims), nrow = n_cells, ncol = n_dims)
  cl  <- sample.int(n_clusters, n_cells, replace = TRUE)
  list(emb = emb, cl = cl, n_cells = n_cells, n_dims = n_dims,
       n_clusters = n_clusters)
}

# ── ComputePseudotimeFuzzy ───────────────────────────────────────────────
test_that("ComputePseudotimeFuzzy returns numeric vector of correct length", {
  knn <- make_toy_knn(50, 5)
  pt <- ComputePseudotimeFuzzy(knn$idx, knn$dist, start_cell = 1L)
  expect_type(pt, "double")
  expect_length(pt, knn$n_cells)
  expect_true(all(is.finite(pt) | is.nan(pt)))
  # The root cell has pseudotime zero.
  expect_equal(pt[1], 0.0)
})

test_that("ComputePseudotimeFuzzy ignores 0 and NA padded neighbours", {
  idx_valid <- matrix(c(2L, 1L, 2L), nrow = 3L, ncol = 1L, byrow = TRUE)
  dist_valid <- matrix(1, nrow = 3L, ncol = 1L)
  idx_padded <- matrix(c(2L, 0L, 1L, NA_integer_, 2L, 0L),
                       nrow = 3L, ncol = 2L, byrow = TRUE)
  dist_padded <- matrix(c(1, 0, 1, 0, 1, 0), nrow = 3L, ncol = 2L, byrow = TRUE)

  pt_valid <- ComputePseudotimeFuzzy(idx_valid, dist_valid, start_cell = 1L,
                                     rank_normalize = FALSE)
  pt_padded <- ComputePseudotimeFuzzy(idx_padded, dist_padded, start_cell = 1L,
                                      rank_normalize = FALSE)

  expect_equal(pt_padded, pt_valid, tolerance = 1e-6)
})

test_that("R wrappers reject out-of-range cell indices before Rust", {
  knn <- make_toy_knn(10, 2)
  expect_error(
    ComputePseudotimeFuzzy(knn$idx, knn$dist, start_cell = 0L),
    "start_cell"
  )
  bad_idx <- knn$idx
  bad_idx[1, 1] <- 99L
  expect_error(
    ComputePseudotimeFuzzy(bad_idx, knn$dist, start_cell = 1L),
    "nn_idx"
  )
})

test_that("BuildKnnFromSparse skips self neighbours", {
  skip("BuildKnnFromSparse is not exported in this release")
  skip_if_not_installed("Matrix")
  dist <- Matrix::sparseMatrix(
    i = c(1L, 2L, 3L), j = c(2L, 3L, 1L),
    x = c(0.2, 0.3, 0.4), dims = c(3L, 3L)
  )
  conn <- Matrix::sparseMatrix(
    i = c(1L, 2L, 3L, 1L), j = c(1L, 2L, 3L, 3L),
    x = c(1, 1, 1, 1), dims = c(3L, 3L)
  )
  knn <- BuildKnnFromSparse(dist, conn, n_cells = 3L, k = 2L)

  expect_false(any(knn$nn_idx > 0L & knn$nn_idx == row(knn$nn_idx)))
  expect_equal(knn$nn_idx[, 1L], c(2L, 3L, 1L))
})

# ── ExtractLineages ──────────────────────────────────────────────────────
test_that("ExtractLineages returns list of integer vectors", {
  toy <- make_toy_embedding(50, 5, 4)
  lineages <- ExtractLineages(toy$emb, toy$cl, start_cluster = 1L)
  expect_type(lineages, "list")
  expect_true(length(lineages) >= 1)
  for (lin in lineages) {
    expect_type(lin, "integer")
    expect_true(lin[1] == 1L)  # starts at root cluster
  }
})

test_that("ExtractLineages preserves non-contiguous input cluster labels", {
  emb <- matrix(c(0, 0, 1, 0, 2, 0), nrow = 3L, ncol = 2L, byrow = TRUE)
  cl <- c(10L, 20L, 30L)
  lineages <- ExtractLineages(emb, cl, start_cluster = 10L)

  expect_true(length(lineages) >= 1L)
  expect_equal(lineages[[1L]][1L], 10L)
  expect_true(all(unlist(lineages, use.names = FALSE) %in% cl))
})

# ── ComputeLineageWeights ───────────────────────────────────────────────
test_that("ComputeLineageWeights returns matrix with rows summing to 1", {
  toy <- make_toy_embedding(50, 5, 4)
  lw <- ComputeLineageWeights(toy$emb, toy$cl, start_cluster = 1L, sigma = 1.0)
  expect_true(is.matrix(lw))
  expect_equal(nrow(lw), toy$n_cells)
  expect_true(ncol(lw) >= 1)
  row_sums <- rowSums(lw)
  expect_true(all(abs(row_sums - 1.0) < 1e-4))
})

test_that("ComputeLineageWeightsTrunkAware accepts non-contiguous cluster labels", {
  emb <- matrix(c(0, 0, 1, 0, 2, 0, 3, 0), nrow = 4L, ncol = 2L, byrow = TRUE)
  cl <- c(10L, 10L, 30L, 50L)
  lw <- ComputeLineageWeightsTrunkAware(emb, cl, start_cluster = 10L, sigma = 1.0)

  expect_true(is.matrix(lw))
  expect_equal(nrow(lw), nrow(emb))
  expect_true(all(is.finite(lw)))
  expect_true(all(abs(rowSums(lw) - 1.0) < 1e-4))
})

# ── ComputeFateProbabilities ─────────────────────────────────────────────
test_that("ComputeFateProbabilities returns matrix with valid probabilities", {
  knn <- make_toy_knn(50, 5)
  pt <- ComputePseudotimeFuzzy(knn$idx, knn$dist, start_cell = 1L)
  terminals <- c(10L, 20L)
  fate <- ComputeFateProbabilities(knn$idx, knn$dist, pt,
                                   terminal_cells = terminals,
                                   bandwidth = 0.1, max_iter = 100L)
  expect_true(is.matrix(fate))
  expect_equal(nrow(fate), knn$n_cells)
  expect_equal(ncol(fate), length(terminals))
  expect_true(all(fate >= -1e-6))  # non-negative (allow float tolerance)
})

test_that("ComputeFateProbabilities rejects invalid terminal cells", {
  knn <- make_toy_knn(10, 2)
  pt <- ComputePseudotimeFuzzy(knn$idx, knn$dist, start_cell = 1L)
  expect_error(
    ComputeFateProbabilities(knn$idx, knn$dist, pt, terminal_cells = 11L),
    "terminal_cells"
  )
})

# ── ComputeLRScores ──────────────────────────────────────────────────────
test_that("ComputeLRScores returns data.frame with expected columns", {
  skip_if_not_installed("Matrix")
  n_genes <- 20
  n_cells <- 30
  n_types <- 3
  set.seed(42)
  x <- Matrix::rsparsematrix(n_genes, n_cells, density = 0.3)
  x@x <- abs(x@x)
  ct <- sample.int(n_types, n_cells, replace = TRUE)
  lr_db <- list(
    ligand_rows   = list(1L, 3L),
    receptor_rows = list(2L, 4L)
  )
  res <- ComputeLRScores(x, ct, lr_db, n_perms = 10L, seed = 1, min_pct = 0.0)
  expect_true(is.data.frame(res))
  expect_true(all(c("source", "target", "pair_idx", "score", "pvalue") %in% colnames(res)))
  if (nrow(res) > 0) {
    expect_true(all(res$score >= 0))
    expect_true(all(res$pvalue >= 0 & res$pvalue <= 1))
  }
})

test_that("ComputeLRScores rejects out-of-range LR gene indices", {
  skip_if_not_installed("Matrix")
  x <- Matrix::rsparsematrix(5, 10, density = 0.3)
  x@x <- abs(x@x)
  ct <- sample.int(2L, 10L, replace = TRUE)
  bad_lr <- list(ligand_rows = list(1L), receptor_rows = list(6L))
  expect_error(ComputeLRScores(x, ct, bad_lr, n_perms = 1L), "lr_db")
})

test_that("ComputeLRScores validates permutation parameters before Rust", {
  skip_if_not_installed("Matrix")
  x <- Matrix::Matrix(matrix(seq_len(40), nrow = 4L), sparse = TRUE)
  ct <- rep(1:2, each = 5L)
  lr_db <- list(ligand_rows = list(1L), receptor_rows = list(2L))

  expect_error(ComputeLRScores(x, ct, lr_db, n_perms = -1L), "n_perms")
  expect_error(ComputeLRScores(x, ct, lr_db, min_pct = 1.1), "min_pct")
  expect_error(ComputeLRScores(x, ct, lr_db, seed = Inf), "seed")
})

test_that("LR wrappers reject invalid expression and drop explicit zeros", {
  skip_if_not_installed("Matrix")
  x <- Matrix::sparseMatrix(
    i = c(1L, 2L, 1L, 2L),
    j = c(1L, 1L, 2L, 2L),
    x = c(1, 1, 2, 2),
    dims = c(2L, 4L)
  )
  x@x[1L] <- 0
  ct <- rep(1L, 4L)
  lr_db <- list(ligand_rows = list(1L), receptor_rows = list(2L))

  observed <- ComputeLRScores(x, ct, lr_db, n_perms = 1L, min_pct = 0)
  expected <- ComputeLRScores(Matrix::drop0(x), ct, lr_db, n_perms = 1L, min_pct = 0)
  expect_equal(observed, expected, tolerance = 1e-7)

  negative <- x
  negative@x[2L] <- -1
  nonfinite <- x
  nonfinite@x[2L] <- Inf
  expect_error(ComputeLRScores(negative, ct, lr_db, n_perms = 1L), "non-negative")
  expect_error(ComputeLRScores(nonfinite, ct, lr_db, n_perms = 1L), "finite")
  expect_error(ComputeLRScores(as.matrix(x), ct, lr_db, n_perms = 1L), "dgCMatrix")
})

test_that("LR database rejects empty complexes and deduplicates subunits", {
  skip_if_not_installed("Matrix")
  x <- Matrix::Matrix(matrix(seq_len(40), nrow = 4L), sparse = TRUE)
  ct <- rep(1:2, each = 5L)
  empty_lr <- list(ligand_rows = list(integer()), receptor_rows = list(2L))
  expect_error(ComputeLRScores(x, ct, empty_lr, n_perms = 1L), "empty complexes")

  duplicate_lr <- list(ligand_rows = list(c(1L, 1L)), receptor_rows = list(c(2L, 2L)))
  unique_lr <- list(ligand_rows = list(1L), receptor_rows = list(2L))
  expect_equal(
    ComputeLRScores(x, ct, duplicate_lr, n_perms = 1L),
    ComputeLRScores(x, ct, unique_lr, n_perms = 1L),
    tolerance = 1e-7
  )
})

test_that("CCC wrappers compact cell types and restore original labels", {
  skip_if_not_installed("Matrix")
  x <- Matrix::Matrix(matrix(seq_len(72), nrow = 6L), sparse = TRUE)
  ct <- rep(c(2L, 10L, 25L), 4L)
  lr_db <- list(ligand_rows = list(1L), receptor_rows = list(2L))
  cluster <- ComputeLRScores(x, ct, lr_db, n_perms = 1L, min_pct = 0)
  expect_setequal(unique(c(cluster$source, cluster$target)), c(2L, 10L, 25L))
  expect_equal(attr(cluster, "celltype_levels"), c(2L, 10L, 25L))

  temporal <- ComputeTemporalLR(
    x, ct, lr_db, seq(0, 1, length.out = ncol(x)),
    matrix(1, nrow = ncol(x), ncol = 1L),
    n_bins = 2L, min_pct = 0
  )
  expect_true(all(temporal$source %in% c(2L, 10L, 25L)))
  expect_true(all(temporal$target %in% c(2L, 10L, 25L)))
})

# ── ComputeTemporalLR ───────────────────────────────────────────────────
test_that("ComputeTemporalLR returns data.frame with expected columns", {
  skip_if_not_installed("Matrix")
  n_genes <- 20
  n_cells <- 30
  n_types <- 3
  set.seed(42)
  x <- Matrix::rsparsematrix(n_genes, n_cells, density = 0.3)
  x@x <- abs(x@x)
  ct <- sample.int(n_types, n_cells, replace = TRUE)
  pt <- seq(0, 1, length.out = n_cells)
  lw <- matrix(1.0, nrow = n_cells, ncol = 1)  # single lineage
  lr_db <- list(
    ligand_rows   = list(1L, 3L),
    receptor_rows = list(2L, 4L)
  )
  res <- ComputeTemporalLR(x, ct, lr_db, pt, lw, n_bins = 4L, min_pct = 0.0)
  expect_true(is.data.frame(res))
  expect_true(all(c(
    "pair_idx", "source", "target", "lineage", "bin", "score",
    "structural_valid", "expression_supported"
  ) %in% colnames(res)))
})

test_that("complete temporal grid preserves structural zeros used by calibration", {
  skip_if_not_installed("Matrix")
  expr <- Matrix::Matrix(
    rbind(
      c(1, 1, 0, 0, 2, 2),
      c(1, 1, 1, 1, 2, 2)
    ),
    sparse = TRUE
  )
  celltypes <- rep(1L, 6L)
  pseudotime <- seq(0, 1, length.out = 6L)
  lineage_weights <- matrix(1, nrow = 6L, ncol = 1L)
  lr_db <- list(ligand_rows = list(1L), receptor_rows = list(2L))

  complete <- ComputeTemporalLR(
    expr, celltypes, lr_db, pseudotime, lineage_weights,
    n_bins = 3L, min_pct = 0.5, min_n_eff = 1,
    return_complete_grid = TRUE
  )
  sparse <- ComputeTemporalLR(
    expr, celltypes, lr_db, pseudotime, lineage_weights,
    n_bins = 3L, min_pct = 0.5, min_n_eff = 1,
    return_complete_grid = FALSE
  )
  calibrated <- ComputeTemporalPermutation(
    expr, celltypes, lr_db, pseudotime, lineage_weights,
    n_bins = 3L, min_pct = 0.5, min_n_eff = 1,
    min_valid_bins = 2L, min_valid_fraction = 0,
    n_perms = 2L, min_perm_observed = 0L, min_perm_coverage = 0
  )

  expect_equal(complete$bin, 1:3)
  expect_equal(complete$structural_valid, rep(TRUE, 3L))
  expect_equal(complete$expression_supported, c(TRUE, FALSE, TRUE))
  expect_equal(complete$score[2L], 0)
  expect_equal(nrow(sparse), 2L)
  expect_equal(
    calibrated$trajectories$temporal_stat[1L],
    max(complete$score) - min(complete$score),
    tolerance = 1e-6
  )
})

test_that("ComputeTemporalLR requires one-dimensional pseudotime", {
  skip_if_not_installed("Matrix")
  n_genes <- 20
  n_cells <- 30
  n_types <- 3
  set.seed(42)
  x <- Matrix::rsparsematrix(n_genes, n_cells, density = 0.3)
  x@x <- abs(x@x)
  ct <- sample.int(n_types, n_cells, replace = TRUE)
  pt <- cbind(seq(0, 1, length.out = n_cells), seq(1, 0, length.out = n_cells))
  lw <- matrix(0.5, nrow = n_cells, ncol = 2)
  lr_db <- list(
    ligand_rows   = list(1L, 3L),
    receptor_rows = list(2L, 4L)
  )

  expect_error(
    ComputeTemporalLR(x, ct, lr_db, pt, lw, n_bins = 4L, min_pct = 0.0),
    regexp = "pseudotime"
  )
})

test_that("ComputeTemporalLR skips disconnected pseudotime cells", {
  skip_if_not_installed("Matrix")
  nn_idx <- matrix(c(2L, 1L, 2L, 3L, 6L, 5L), ncol = 1L)
  nn_dists <- matrix(1, nrow = 6L, ncol = 1L)
  pt <- ComputePseudotimeFuzzy(nn_idx, nn_dists, start_cell = 1L)
  expect_true(all(is.nan(pt[5:6])))

  expr <- Matrix::Matrix(
    matrix(c(
      1, 2, 3, 4, 100, 200,
      2, 3, 4, 5, 200, 300,
      3, 4, 5, 6, 300, 400,
      4, 5, 6, 7, 400, 500
    ), nrow = 4L, byrow = TRUE),
    sparse = TRUE
  )
  celltypes <- c(1L, 2L, 1L, 2L, 1L, 2L)
  lineage_weights <- matrix(1, nrow = 6L, ncol = 1L)
  lr_db <- list(ligand_rows = list(1L), receptor_rows = list(2L))

  observed <- ComputeTemporalLR(
    expr, celltypes, lr_db, pt, lineage_weights,
    n_bins = 2L, min_pct = 0
  )
  expected <- ComputeTemporalLR(
    expr[, 1:4], celltypes[1:4], lr_db, pt[1:4], lineage_weights[1:4, , drop = FALSE],
    n_bins = 2L, min_pct = 0
  )
  expect_equal(observed, expected, tolerance = 1e-7)
})

test_that("ComputeTemporalLR rejects infinite and out-of-range pseudotime", {
  skip_if_not_installed("Matrix")
  x <- Matrix::Matrix(matrix(seq_len(40), nrow = 4L), sparse = TRUE)
  ct <- rep(1:2, each = 5L)
  lw <- matrix(1, nrow = 10L, ncol = 1L)
  lr_db <- list(ligand_rows = list(1L), receptor_rows = list(2L))
  pt <- seq(0, 1, length.out = 10L)

  expect_error(
    ComputeTemporalLR(x, ct, lr_db, replace(pt, 1L, Inf), lw, n_bins = 2L),
    "Inf"
  )
  expect_error(
    ComputeTemporalLR(x, ct, lr_db, replace(pt, 1L, 1.1), lw, n_bins = 2L),
    "\\[0, 1\\]"
  )
})

# ── TrajectoryDE ─────────────────────────────────────────────────────────
test_that("TrajectoryDE returns data.frame with rho/pval/qval", {
  n_genes <- 10
  n_cells <- 50
  set.seed(42)
  expr <- matrix(abs(rnorm(n_genes * n_cells)), nrow = n_genes, ncol = n_cells)
  pt <- seq(0, 1, length.out = n_cells)
  res <- TrajectoryDE(expr, pt, min_cells = 5L)
  expect_true(is.data.frame(res))
  expect_true(all(c("rho", "pval", "qval") %in% colnames(res)))
  expect_equal(nrow(res), n_genes)
  expect_true(all(res$pval >= 0 & res$pval <= 1))
  expect_true(all(res$qval >= 0 & res$qval <= 1))
})
test_that("public integer inputs are rejected before coercion", {
  idx <- matrix(c(2, 1, 2, 1), nrow = 2)
  d <- matrix(c(1, 1, 2, 2), nrow = 2)
  expect_error(
    ComputePseudotimeFuzzy(idx, d, start_cell = 1.5),
    "integer-valued"
  )
  expect_error(
    ComputePseudotimeFuzzy(idx, d, start_cell = 1L, rank_normalize = NA),
    "TRUE or FALSE"
  )
})

test_that("kNN validation rejects negative distances and bad indices", {
  idx <- matrix(c(2, 1), nrow = 2)
  expect_error(
    ComputePseudotimeFuzzy(idx, matrix(c(-1, 1), nrow = 2), 1L),
    "non-negative"
  )
  expect_error(
    ComputePseudotimeFuzzy(matrix(c(3, 1), nrow = 2),
                           matrix(c(1, 1), nrow = 2), 1L),
    "outside"
  )
})
