make_fork_score_fixture <- function() {
  grid <- expand.grid(
    pair_idx = 1:3,
    source = 1L,
    target = 2L,
    lineage = 1:2,
    bin = 1:6,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  shared <- c(1, 1, 1, 3, 4, 5)
  focal <- c(1, 1, 1, 5, 6, 7)
  flat <- rep(0.05, 6)
  profiles <- list(
    `1.1` = shared,
    `1.2` = shared,
    `2.1` = focal,
    `2.2` = shared,
    `3.1` = flat,
    `3.2` = flat
  )
  key <- paste(grid$pair_idx, grid$lineage, sep = ".")
  grid$score <- unlist(
    Map(function(k, b) profiles[[k]][b], key, grid$bin),
    use.names = FALSE
  )
  grid$n_eff_cells <- 20
  grid$structurally_valid <- TRUE
  list(
    scores = grid,
    bif = data.frame(
      lineage = 1:2,
      b_star = c(3L, 3L),
      reason = "crossing"
    ),
    sibling = data.frame(
      focal_lineage = c(1L, 2L),
      sibling_lineage = c(2L, 1L),
      fork_cluster = 1L,
      shared_prefix_length = 1L
    )
  )
}

# Two-lineage score table on one shared grid.  Both branches evaluate the
# same smooth activity curve `g(t) = 0.1 + 0.6 * t^2` at identical bins,
# mirroring the audit counterexample where per-lineage windows produce a
# spurious contrast.
make_common_grid_fork_fixture <- function(n_bins = 8L) {
  bounds <- seq(0, 1, length.out = n_bins + 1L)
  mid <- (bounds[-1L] + bounds[-length(bounds)]) / 2
  g <- 0.1 + 0.6 * mid^2
  grid <- expand.grid(
    pair_idx = 1L,
    source = 1L,
    target = 2L,
    lineage = 1:2,
    bin = seq_len(n_bins),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  grid$score <- g[grid$bin]
  grid$n_eff_cells <- 20
  grid$structural_valid <- TRUE
  attr(grid, "temporal_coordinate") <- "common_pseudotime_grid"
  list(
    scores = grid,
    mid = mid,
    bif = data.frame(
      lineage = 1:2,
      b_star = c(4L, 6L),
      reason = "crossing"
    ),
    sibling = data.frame(
      focal_lineage = c(1L, 2L),
      sibling_lineage = c(2L, 1L),
      fork_cluster = 2L,
      shared_prefix_length = 1L
    )
  )
}
