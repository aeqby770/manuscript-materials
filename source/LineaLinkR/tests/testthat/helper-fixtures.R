make_temporal_fixture <- function() {
  skip_if_not_installed("Matrix")
  n_genes <- 20
  n_cells <- 40
  n_types <- 3
  set.seed(42)
  x <- Matrix::rsparsematrix(n_genes, n_cells, density = 0.3)
  x@x <- abs(x@x)
  ct <- sample.int(n_types, n_cells, replace = TRUE)
  pt <- seq(0, 1, length.out = n_cells)
  lw <- matrix(c(rep(c(0.7, 0.3), n_cells / 2)), nrow = n_cells, ncol = 2, byrow = TRUE)
  lr_db <- list(
    ligand_rows = list(1L, 3L),
    receptor_rows = list(2L, 4L)
  )
  list(
    x = x, ct = ct, pt = pt, lw = lw, lr_db = lr_db,
    n_genes = n_genes, n_cells = n_cells, n_types = n_types
  )
}

make_start_cell_fixture <- function() {
  set.seed(12345)
  embedding <- rbind(
    cbind(seq(0, 1, length.out = 30), 0),
    cbind(seq(1.1, 3, length.out = 20), 1),
    cbind(seq(1.1, 2, length.out = 20), -1)
  )
  cluster_labels <- c(rep(1:3, each = 10), rep(4:5, each = 10),
                      rep(6:7, each = 10))
  distance_matrix <- as.matrix(stats::dist(embedding))
  diag(distance_matrix) <- Inf
  nn_idx <- t(vapply(
    seq_len(nrow(embedding)),
    function(i) order(distance_matrix[i, ])[1:6],
    integer(6)
  ))
  nn_dists <- matrix(
    distance_matrix[cbind(rep(seq_len(nrow(embedding)), each = 6),
                          as.vector(t(nn_idx)))],
    nrow = nrow(embedding), byrow = TRUE
  )
  list(
    embedding = embedding,
    clusters = cluster_labels,
    nn_idx = nn_idx,
    nn_dists = nn_dists,
    potency = rep(0, nrow(embedding))
  )
}

make_knn_fixture <- function(embedding, k = 8L) {
  distance_matrix <- as.matrix(stats::dist(embedding))
  diag(distance_matrix) <- Inf
  nn_idx <- t(vapply(
    seq_len(nrow(embedding)),
    function(i) order(distance_matrix[i, ])[seq_len(k)],
    integer(k)
  ))
  nn_dists <- matrix(
    distance_matrix[cbind(rep(seq_len(nrow(embedding)), each = k),
                          as.vector(t(nn_idx)))],
    nrow = nrow(embedding), byrow = TRUE
  )
  list(nn_idx = nn_idx, nn_dists = nn_dists)
}

make_symmetric_root_fixture <- function() {
  x <- seq(-1, 1, length.out = 80)
  embedding <- cbind(x, 0.02 * cos(seq(0, 4 * pi, length.out = 80)))
  knn <- make_knn_fixture(embedding, 8L)
  list(
    embedding = embedding,
    clusters = rep(1:4, each = 20),
    nn_idx = knn$nn_idx,
    nn_dists = knn$nn_dists,
    fate_entropy = 1 - abs(x),
    commitment = abs(x),
    potency = 1 - abs(x)
  )
}

make_directional_y_fixture <- function() {
  trunk_x <- seq(0, 1, length.out = 40)
  branch_x <- seq(1.03, 2, length.out = 20)
  embedding <- rbind(
    cbind(trunk_x, 0.02 * sin(seq(0, 2 * pi, length.out = 40))),
    cbind(branch_x, seq(0.03, 0.9, length.out = 20)),
    cbind(branch_x, -seq(0.03, 0.9, length.out = 20))
  )
  clusters <- c(rep(1:4, each = 10), rep(5:6, each = 10),
                rep(7:8, each = 10))
  development <- c(seq(0, 0.55, length.out = 40),
                   seq(0.60, 1, length.out = 20),
                   seq(0.60, 1, length.out = 20))
  knn <- make_knn_fixture(embedding, 8L)
  list(
    embedding = embedding,
    clusters = clusters,
    nn_idx = knn$nn_idx,
    nn_dists = knn$nn_dists,
    fate_entropy = 1 - development,
    commitment = development,
    potency = 1 - development,
    expected_root_cells = which(clusters == 1L)
  )
}
