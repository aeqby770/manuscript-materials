test_that("OrientPseudotimeBlind flips anti-PC1 direction", {
  set.seed(2)
  n <- 120
  emb <- cbind(seq(0, 1, length.out = n) + rnorm(n, sd = 0.03), rnorm(n))
  pt <- seq(1, 0, length.out = n)
  expect_lt(BlindEmbeddingScore(pt, emb), 0)
  pt2 <- OrientPseudotimeBlind(pt, emb)
  expect_gt(BlindEmbeddingScore(pt2, emb), 0)
})

test_that("InferStartCellAuto picks PC1-min on linear trajectory", {
  set.seed(3)
  n <- 60
  k <- 5
  emb <- cbind(seq(0, 1, length.out = n), rnorm(n, sd = 0.1))
  cl <- rep(1:3, each = n / 3)
  d <- as.matrix(stats::dist(emb))
  diag(d) <- Inf
  ord <- t(apply(d, 1L, order))
  idx <- ord[, seq_len(k)]
  dist <- matrix(d[cbind(rep(seq_len(n), each = k), as.vector(t(idx)))],
                 nrow = n, byrow = TRUE)
  auto <- InferStartCellAuto(emb, cl, idx, dist, rank_normalize = FALSE, seed = 1L)
  expect_equal(auto$method, "cluster_pc1_min")
  expect_true(auto$blind_score > 0.5)
  expect_true(auto$ordering_score > 0.5)
  expect_true(!is.null(auto$root_direction_status))
})

test_that("InferStartCellAuto is stable under embedding sign flips", {
  set.seed(31)
  n <- 80
  k <- 6
  x <- seq(0, 1, length.out = n)
  emb <- cbind(x, sin(pi * x) * 0.02)
  cl <- rep(1:4, each = n / 4)
  d <- as.matrix(stats::dist(emb))
  diag(d) <- Inf
  ord <- t(apply(d, 1L, order))
  idx <- ord[, seq_len(k)]
  dist <- matrix(d[cbind(rep(seq_len(n), each = k), as.vector(t(idx)))],
                 nrow = n, byrow = TRUE)
  auto_a <- InferStartCellAuto(emb, cl, idx, dist, rank_normalize = FALSE, seed = 1L)
  emb_flip <- emb
  emb_flip[, 1L] <- -emb_flip[, 1L]
  auto_b <- InferStartCellAuto(emb_flip, cl, idx, dist, rank_normalize = FALSE, seed = 1L)
  expect_equal(auto_a$start_cell, auto_b$start_cell)
  expect_true(isTRUE(auto_a$direction_ambiguous))
  expect_true(isTRUE(auto_b$direction_ambiguous))
  expect_equal(auto_a$root_direction_status, "ambiguous_pc_sign")
  expect_equal(auto_b$root_direction_status, "ambiguous_pc_sign")
  expect_true(auto_a$ordering_score > 0.45)
  expect_true(auto_b$ordering_score > 0.45)
})

test_that("InferStartCellAuto fails when no candidate root covers the graph", {
  n <- 12
  emb <- cbind(seq(0, 1, length.out = n), rep(0, n))
  cl <- rep(1:3, each = 4L)
  idx <- matrix(0L, nrow = n, ncol = 2L)
  dist <- matrix(0, nrow = n, ncol = 2L)

  expect_error(
    InferStartCellAuto(emb, cl, idx, dist, rank_normalize = FALSE, seed = 1L),
    "each cell must have at least one valid non-self neighbour"
  )
})

test_that("AuditPseudotimeDirection reports potency and PC conflicts", {
  pt <- seq(0, 1, length.out = 100)
  emb <- cbind(pt, 0)

  current <- AuditPseudotimeDirection(pt, potency = 1 - pt, embedding = emb)
  expect_equal(current$potency_status, "supports_current_root")
  expect_equal(current$pc_convention_status, "supports_current_root")
  expect_equal(current$direction_status, "supports_current_root")

  conflict <- AuditPseudotimeDirection(pt, potency = pt, embedding = emb)
  expect_equal(conflict$potency_status, "supports_flipped_root")
  expect_equal(conflict$pc_convention_status, "supports_current_root")
  expect_equal(conflict$direction_status, "pc_convention_conflict")
})

test_that("InferStartCellAttractor selects high-potency upstream region", {
  set.seed(4)
  trunk <- cbind(seq(0, 1, length.out = 40), rnorm(40, sd = 0.03))
  branch_a <- cbind(seq(1.05, 2, length.out = 30), seq(0.05, 0.9, length.out = 30))
  branch_b <- cbind(seq(1.05, 2, length.out = 30), -seq(0.05, 0.9, length.out = 30))
  emb <- rbind(trunk, branch_a, branch_b)
  n <- nrow(emb)
  k <- 8
  cl <- c(rep(1:4, each = 10), rep(5:7, each = 10), rep(8:10, each = 10))
  potency <- 1 - ((emb[, 1] - min(emb[, 1])) / diff(range(emb[, 1])))
  idx <- matrix(0L, n, k)
  dist <- matrix(0, n, k)
  for (i in seq_len(n)) {
    d <- sqrt(rowSums((emb - matrix(emb[i, ], nrow = n, ncol = ncol(emb), byrow = TRUE))^2))
    o <- order(d)[2:(k + 1L)]
    idx[i, ] <- o
    dist[i, ] <- d[o]
  }
  root <- InferStartCellAttractor(emb, cl, idx, dist, potency, top_fraction = 0.1, seed = 1L)
  expect_equal(root$method, "attractor_energy")
  expect_true(length(root$root_region) > 1L)
  expect_true(emb[root$start_cell, 1] < 0.5)
  expect_true(all(root$attractor_clusters %in% cl))
  expect_true(root$cellwise_start_cell %in% seq_len(n))
  expect_true(length(root$cellwise_root_region) >= 1L)
  expected_score <- with(
    root$scores,
    0.25 * potency + 0.20 * fate_entropy +
      0.30 * nearest_attractor_distance + 0.15 * source_flow +
      0.10 * (1 - commitment)
  )
  expect_equal(root$scores$root_score, expected_score)
})

test_that("conservative root fallback requires strong positive blind evidence", {
  candidates <- data.frame(
    candidate = c(
      "detected_genes_high_state",
      "detected_genes_low_state",
      "detected_genes_high_cellwise"
    ),
    blind_score = c(-0.01, 0.23, -0.17),
    component_score = c(0.60, 0.58, 0.55)
  )

  retained <- SelectRootCandidateConservative(
    candidates,
    default_candidate = "detected_genes_high_state",
    min_blind_score = 0,
    min_gain = 0.30
  )
  expect_equal(retained$selected_candidate, "detected_genes_high_state")
  expect_equal(retained$selection_status, "default_retained")

  candidates$blind_score[[2L]] <- 0.50
  switched_low <- SelectRootCandidateConservative(
    candidates,
    default_candidate = "detected_genes_high_state",
    min_blind_score = 0,
    min_gain = 0.30
  )
  expect_equal(switched_low$selected_candidate, "detected_genes_low_state")
  expect_equal(switched_low$selection_status, "strong_blind_fallback")
  expect_equal(switched_low$blind_gain, 0.51)

  candidates$blind_score <- c(-0.37, -0.28, 0.25)
  switched_cellwise <- SelectRootCandidateConservative(
    candidates,
    default_candidate = "detected_genes_high_state",
    min_blind_score = 0,
    min_gain = 0.30
  )
  expect_equal(
    switched_cellwise$selected_candidate,
    "detected_genes_high_cellwise"
  )
  expect_equal(switched_cellwise$blind_gain, 0.62)
})

test_that("ordering score is invariant to embedding and trajectory signs", {
  set.seed(41)
  n <- 160L
  latent <- seq(0, 1, length.out = n)
  embedding <- cbind(
    latent + rnorm(n, sd = 0.01),
    latent^2 + rnorm(n, sd = 0.01),
    sin(pi * latent) + rnorm(n, sd = 0.01)
  )
  shuffled <- sample(latent)
  flipped_embedding <- embedding
  flipped_embedding[, 1L] <- -flipped_embedding[, 1L]

  score <- BlindOrderingScore(latent, embedding)
  expect_equal(
    score,
    BlindOrderingScore(latent, flipped_embedding),
    tolerance = 1e-10
  )
  expect_equal(
    score,
    BlindOrderingScore(1 - latent, embedding),
    tolerance = 1e-10
  )
  expect_gt(score, BlindOrderingScore(shuffled, embedding) + 0.30)
})

test_that("ordering selector cannot use low-potency direction as an automatic fallback", {
  candidates <- data.frame(
    candidate = c(
      "detected_genes_high_state",
      "detected_genes_high_cellwise",
      "detected_genes_low_state"
    ),
    ordering_score = c(0.60, 0.82, 0.95),
    component_score = c(0.70, 0.80, 0.90)
  )

  selected <- SelectRootCandidateByOrdering(
    candidates,
    default_candidate = "detected_genes_high_state",
    allowed_candidates = c(
      "detected_genes_high_state",
      "detected_genes_high_cellwise"
    ),
    min_gain = 0.10
  )
  expect_equal(
    selected$selected_candidate,
    "detected_genes_high_cellwise"
  )
  expect_equal(selected$selection_status, "ordering_fallback")
  expect_equal(selected$ordering_gain, 0.22)

  candidates$ordering_score[[2L]] <- 0.66
  retained <- SelectRootCandidateByOrdering(
    candidates,
    default_candidate = "detected_genes_high_state",
    allowed_candidates = c(
      "detected_genes_high_state",
      "detected_genes_high_cellwise"
    ),
    min_gain = 0.10
  )
  expect_equal(retained$selected_candidate, "detected_genes_high_state")
  expect_equal(retained$selection_status, "default_retained")
})

test_that("endpoint selector switches only for strong potency-alignment gain", {
  candidates <- data.frame(
    candidate = c(
      "detected_genes_high_state",
      "detected_genes_high_cellwise",
      "detected_genes_low_state",
      "detected_genes_low_cellwise"
    ),
    potency_mode = rep(c("detected_genes_high", "detected_genes_low"), each = 2L),
    root_representation = rep(c("state_aggregate", "cellwise_component"), 2L),
    potency_spearman = c(-0.31, -0.30, -0.63, -0.42),
    ordering_score = c(0.565, 0.566, 0.522, 0.345),
    component_score = c(0.50, 0.84, 0.61, 0.88)
  )

  selected <- SelectRootEndpointByPotencyAlignment(
    candidates,
    min_alignment_gain = 0.20,
    min_ordering_gain = 0.10
  )
  expect_equal(selected$selected_candidate, "detected_genes_low_state")
  expect_equal(selected$selected_potency_mode, "detected_genes_low")
  expect_equal(selected$endpoint_alignment_gain, 0.32, tolerance = 1e-12)
  expect_equal(selected$endpoint_selection_status, "strong_alignment_switch")
  expect_equal(selected$representation_selection_status, "state_retained")

  candidates$potency_spearman[3L] <- -0.50
  retained <- SelectRootEndpointByPotencyAlignment(
    candidates,
    min_alignment_gain = 0.20,
    min_ordering_gain = 0.10
  )
  expect_equal(retained$selected_candidate, "detected_genes_high_state")
  expect_equal(retained$selected_potency_mode, "detected_genes_high")
  expect_equal(retained$endpoint_selection_status, "default_mode_retained")
})

test_that("endpoint selector blocks an internal-trunk terminal candidate", {
  candidates <- data.frame(
    candidate = c(
      "detected_genes_high_state",
      "detected_genes_high_cellwise",
      "detected_genes_low_state",
      "detected_genes_low_cellwise"
    ),
    potency_mode = rep(c("detected_genes_high", "detected_genes_low"), each = 2L),
    root_representation = rep(c("state_aggregate", "cellwise_component"), 2L),
    potency_spearman = c(-0.31, -0.30, -0.63, -0.42),
    ordering_score = c(0.565, 0.566, 0.522, 0.345),
    component_score = c(0.50, 0.84, 0.61, 0.88),
    terminal_leaf_eligible = c(FALSE, FALSE, FALSE, FALSE)
  )

  blocked <- SelectRootEndpointByPotencyAlignment(
    candidates,
    min_alignment_gain = 0.20,
    min_ordering_gain = 0.10
  )
  expect_equal(blocked$selected_candidate, "detected_genes_high_state")
  expect_equal(blocked$unconstrained_best_alignment_mode, "detected_genes_low")
  expect_equal(blocked$endpoint_topology_gate_status, "nonterminal_switch_blocked")

  candidates$terminal_leaf_eligible[3L] <- TRUE
  allowed <- SelectRootEndpointByPotencyAlignment(
    candidates,
    min_alignment_gain = 0.20,
    min_ordering_gain = 0.10
  )
  expect_equal(allowed$selected_candidate, "detected_genes_low_state")
  expect_equal(allowed$endpoint_topology_gate_status, "terminal_switch_allowed")
})

test_that("endpoint selector is row-order invariant and keeps ordering gate", {
  candidates <- data.frame(
    candidate = c(
      "detected_genes_high_state",
      "detected_genes_high_cellwise",
      "detected_genes_low_state",
      "detected_genes_low_cellwise"
    ),
    potency_mode = rep(c("detected_genes_high", "detected_genes_low"), each = 2L),
    root_representation = rep(c("state_aggregate", "cellwise_component"), 2L),
    potency_spearman = c(-0.40, -0.39, -0.65, -0.70),
    ordering_score = c(0.45, 0.56, 0.52, 0.57),
    component_score = c(0.70, 0.80, 0.75, 0.85)
  )

  selected <- SelectRootEndpointByPotencyAlignment(
    candidates,
    min_alignment_gain = 0.20,
    min_ordering_gain = 0.10
  )
  permuted <- SelectRootEndpointByPotencyAlignment(
    candidates[c(4L, 2L, 1L, 3L), ],
    min_alignment_gain = 0.20,
    min_ordering_gain = 0.10
  )
  expect_equal(selected$selected_candidate, "detected_genes_low_state")
  expect_equal(permuted$selected_candidate, selected$selected_candidate)

  candidates$ordering_score[4L] <- 0.63
  cellwise <- SelectRootEndpointByPotencyAlignment(
    candidates,
    min_alignment_gain = 0.20,
    min_ordering_gain = 0.10,
    min_representation_alignment_gain = 0.00
  )
  expect_equal(cellwise$selected_candidate, "detected_genes_low_cellwise")
  expect_equal(
    cellwise$representation_selection_status,
    "ordering_alignment_switch"
  )
  expect_equal(
    cellwise$representation_alignment_gain,
    0.05,
    tolerance = 1e-12
  )
})

test_that("endpoint selector blocks cellwise geometry gain when potency alignment degrades", {
  candidates <- data.frame(
    candidate = c(
      "detected_genes_high_state",
      "detected_genes_high_cellwise"
    ),
    potency_mode = rep("detected_genes_high", 2L),
    root_representation = c("state_aggregate", "cellwise_component"),
    potency_spearman = c(-0.58, -0.12),
    ordering_score = c(0.30, 0.51),
    component_score = c(0.68, 0.85)
  )

  selected <- SelectRootEndpointByPotencyAlignment(
    candidates,
    min_alignment_gain = 0.20,
    min_ordering_gain = 0.10,
    min_representation_alignment_gain = 0.00
  )
  expect_equal(selected$selected_candidate, "detected_genes_high_state")
  expect_equal(selected$representation_selection_status, "state_retained")
  expect_equal(selected$representation_alignment_gain, 0)
})

test_that("attractor root selection aggregates evidence at state level", {
  root_score <- c(0.99, 0.10, 0.10, 0.10, rep(0.70, 4L))
  clusters <- rep(1:2, each = 4L)
  neighbors <- matrix(c(
    2L, 3L,
    1L, 3L,
    1L, 2L,
    1L, 2L,
    6L, 7L,
    5L, 7L,
    5L, 6L,
    5L, 6L
  ), nrow = 8L, byrow = TRUE)
  potency <- c(rep(0.5, 5L), 0.9, 0.5, 0.5)

  selected <- LineaLinkR:::.select_cluster_aggregate_root(
    root_score,
    clusters,
    neighbors,
    potency
  )

  expect_equal(selected$root_cluster, 2L)
  expect_equal(selected$start_cell, 6L)
  expect_equal(selected$root_region, 5:8)
})

test_that("attractor absorption uses a row-stochastic sparse transition", {
  idx <- matrix(c(
    2L, 3L,
    1L, 3L,
    1L, 2L,
    0L, 0L
  ), nrow = 4L, byrow = TRUE)
  dist <- matrix(c(
    1, 2,
    1, 1,
    2, 1,
    Inf, Inf
  ), nrow = 4L, byrow = TRUE)

  transition <- LineaLinkR:::.knn_sparse_transition(idx, dist, 4L)

  expect_s4_class(transition, "dgCMatrix")
  expect_equal(as.numeric(Matrix::rowSums(transition)), rep(1, 4L))
  expect_equal(as.numeric(transition[4L, ]), c(0, 0, 0, 1))

  expected <- LineaLinkR:::.markov_fate_probs(
    3L, transition[1:3, 1:3], c(1L, 3L), 300L, 1e-5
  )
  native <- matrix(
    LineaLinkR:::rs_compute_attractor_fate_f64(
      as.integer(t(idx[1:3, , drop = FALSE])),
      as.numeric(t(dist[1:3, , drop = FALSE])),
      as.integer(c(1L, 3L)),
      3L,
      2L,
      300L,
      1e-5
    ),
    nrow = 3L,
    byrow = TRUE
  )
  expect_equal(native, expected, tolerance = 1e-10)
  observed <- ComputeFateProbabilities(
    idx[1:3, , drop = FALSE],
    dist[1:3, , drop = FALSE],
    pseudotime = rep(0, 3L),
    terminal_cells = c(1L, 3L),
    bandwidth = 1,
    max_iter = 300L,
    tol = 1e-5
  )
  expect_equal(observed, expected, tolerance = 2e-5)
})

test_that("absorption failures reach R with their numerical cause", {
  expect_error(
    LineaLinkR:::rs_compute_attractor_fate_f64(
      as.integer(c(2, 3, 1, 3, 1, 2)), c(1, 1, 1, 2, 1, 1),
      as.integer(c(1, 3)), 3L, 2L, 0L, 1e-5
    ),
    "did not converge"
  )
  expect_error(
    LineaLinkR:::rs_compute_attractor_fate_f64(
      1:3, rep(1, 3), as.integer(c(1, 3)), 3L, 1L, 300L, 1e-5
    ),
    "cannot reach any terminal"
  )
})

test_that("undefined root candidates retain missing scores and unreachable cells", {
  embedding <- matrix(c(-1, 0, 1), ncol = 1L)
  neighbors <- matrix(1:3, ncol = 1L)
  distances <- matrix(1, nrow = 3L, ncol = 1L)
  expect_error(
    InferStartCellAttractor(embedding, 1:3, neighbors, distances, rep(0, 3)),
    "1 cells cannot reach any terminal"
  )
  unavailable <- InferStartCellAttractor(
    embedding, 1:3, neighbors, distances, rep(0, 3), on_unreachable = "return"
  )
  expect_identical(unavailable$status, "unreachable_terminals")
  expect_identical(unavailable$unreachable_cells, 2L)
  expect_identical(unavailable$start_cell, NA_integer_)
  expect_identical(unavailable$cellwise_start_cell, NA_integer_)
  expect_null(unavailable$scores)
  expect_identical(unavailable$attractor_clusters, c(1L, 3L))
  expect_identical(
    LineaLinkR:::rs_attractor_unreachable_cells(c(1L, 2L, 2L, 3L), rep(1, 4),
      c(1L, 4L), 4L, 1L), c(2L, 3L)
  )
  connected <- matrix(c(2L, 3L, 1L, 3L, 1L, 2L), nrow = 3L, byrow = TRUE)
  recovered <- InferStartCellAttractor(embedding, 1:3, connected,
    matrix(c(1, 1, 1, 2, 1, 1), 3L, 2L, byrow = TRUE), rep(0, 3),
    max_iter = 0L, on_unreachable = "return")
  expect_identical(recovered$status, "available")
  expect_identical(recovered$fate_solver, "sparse_lu")
  expect_lte(recovered$fate_residual_max, 1e-10)
})

test_that("unavailable root modes cannot win or bypass the terminal criterion", {
  candidates <- data.frame(
    candidate = c("high_state", "high_cell", "low_state", "low_cell"),
    potency_mode = rep(c("detected_genes_high", "detected_genes_low"), each = 2L),
    root_representation = rep(c("state_aggregate", "cellwise_component"), 2L),
    potency_spearman = c(-0.4, -0.4, NA, NA),
    ordering_score = c(0.5, 0.5, NA, NA),
    component_score = c(0.7, 0.7, NA, NA),
    terminal_leaf_eligible = c(FALSE, FALSE, FALSE, FALSE),
    available = c(TRUE, TRUE, FALSE, FALSE)
  )
  retained <- SelectRootEndpointByPotencyAlignment(candidates)
  expect_identical(retained$selected_candidate, "high_state")
  expect_identical(retained$unavailable_candidates, c("low_state", "low_cell"))
  expect_true(retained$default_available)

  candidates$available <- c(FALSE, FALSE, TRUE, TRUE)
  candidates$potency_spearman <- c(NA, NA, -0.6, -0.6)
  candidates$ordering_score <- c(NA, NA, 0.6, 0.6)
  candidates$component_score <- c(NA, NA, 0.8, 0.7)
  expect_error(SelectRootEndpointByPotencyAlignment(candidates), "terminal-leaf criterion")
  candidates$terminal_leaf_eligible[3L] <- TRUE
  selected <- SelectRootEndpointByPotencyAlignment(candidates[c(4L, 2L, 3L, 1L), ])
  expect_identical(selected$selected_candidate, "low_state")
  expect_identical(selected$endpoint_selection_status, "default_unavailable")
  expect_true(is.na(selected$endpoint_alignment_gain))
  expect_false(selected$default_available)
  candidates$available <- rep(FALSE, 4L)
  expect_error(SelectRootEndpointByPotencyAlignment(candidates), "No available root candidate")
})

test_that("retention terminals reflect state transitions rather than self edges", {
  embedding <- matrix(seq(-3, 4), ncol = 1L)
  clusters <- c(10L, 10L, 10L, 20L, 20L, 30L, 30L, 40L)
  neighbors <- matrix(c(2L, 3L, 1L, 3L, 2L, 4L, 3L, 5L,
                        4L, 6L, 5L, 7L, 6L, 8L, 7L, 6L), ncol = 2L, byrow = TRUE)
  distances <- matrix(c(1, 2, 1, 1, 1, 1, 1, 1,
                        1, 1, 2, 1, 1, 2, 1, 2), ncol = 2L, byrow = TRUE)
  selected <- InferStartCellAttractor(embedding, clusters, neighbors, distances,
    rep(0, 8L), n_attractors = 2L, attractor_selection = "retention")
  self_edges <- InferStartCellAttractor(embedding, clusters,
    cbind(seq_len(8L), neighbors), cbind(0, distances), rep(0, 8L),
    n_attractors = 2L, attractor_selection = "retention")
  expect_identical(selected$attractor_clusters, c(10L, 30L))
  expect_identical(selected$attractor_selection, "retention")
  expect_equal(selected$attractor_states$retention[c(1L, 2L, 4L)], c(5 / 6, 0.5, 0))
  expect_equal(selected$attractor_states, self_edges$attractor_states)
  expect_identical(selected$attractor_clusters, self_edges$attractor_clusters)
  direct <- InferStartCellAttractor(embedding, clusters,
    cbind(seq_len(8L), neighbors), cbind(0, distances), rep(0, 8L),
    n_attractors = 2L, attractor_selection = "retention", fate_solver = "sparse_lu")
  expect_identical(direct$fate_solver, "sparse_lu")
  expect_equal(direct$scores$fate_entropy, selected$scores$fate_entropy, tolerance = 1e-8)
  expect_equal(direct$scores$commitment, selected$scores$commitment, tolerance = 1e-8)
  expect_identical(direct$start_cell, selected$start_cell)
  expect_lt(direct$fate_residual_max, 1e-12)
})
