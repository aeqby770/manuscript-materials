# kNN graph helpers for trajectory inference (benchmark + production).

#' Add cluster-MST bridge edges between lineage-adjacent clusters.
#'
#' When Scanpy kNN splits a bifurcated dataset into disconnected components,
#' Dijkstra pseudotime leaves many cells as NaN. For each consecutive cluster
#' pair on an MST lineage path, connect the closest cell pair in embedding space.
#'
#' @param nn_idx integer matrix (`n_cells` x `k`), 1-based neighbor ids (0 = pad).
#' @param nn_dists numeric matrix (`n_cells` x `k`), L2 distances.
#' @param embedding numeric matrix (`n_cells` x `n_dims`).
#' @param cluster_labels integer vector, 1-based cluster ids.
#' @param start_cluster integer, 1-based root cluster for MST lineages.
#' @return list with `nn_idx` and `nn_dists`.
#' @export
AugmentKnnWithMstBridges <- function(nn_idx, nn_dists, embedding,
                                     cluster_labels, start_cluster) {
  stopifnot(is.matrix(nn_idx), is.matrix(nn_dists),
            nrow(nn_idx) == nrow(nn_dists))
  emb <- as.matrix(embedding)
  cl  <- as.integer(cluster_labels)
  n_cells <- nrow(emb)
  lineages <- ExtractLineages(emb, cl, start_cluster)
  if (!length(lineages)) {
    return(list(nn_idx = nn_idx, nn_dists = nn_dists))
  }

  cluster_cells <- split(seq_len(n_cells), cl)
  for (path in lineages) {
    if (length(path) < 2L) next
    for (t in seq_len(length(path) - 1L)) {
      ca <- path[t]; cb <- path[t + 1L]
      ia <- cluster_cells[[as.character(ca)]]
      ib <- cluster_cells[[as.character(cb)]]
      if (is.null(ia) || is.null(ib) || !length(ia) || !length(ib)) next
      best <- .nearest_cell_pair(emb, ia, ib)
      br <- .insert_neighbor(nn_idx, nn_dists, best$i, best$j, best$d)
      nn_idx <- br$nn_idx; nn_dists <- br$nn_dists
      br <- .insert_neighbor(nn_idx, nn_dists, best$j, best$i, best$d)
      nn_idx <- br$nn_idx; nn_dists <- br$nn_dists
    }
  }
  list(nn_idx = nn_idx, nn_dists = nn_dists)
}

#' Connect disconnected kNN components with nearest-embedding bridges.
#'
#' @param nn_idx integer matrix (`n_cells` x `k`), 1-based neighbor ids (0 = pad).
#' @param nn_dists numeric matrix (`n_cells` x `k`), L2 distances.
#' @param embedding numeric matrix (`n_cells` x `n_dims`).
#' @return list with `nn_idx`, `nn_dists`, and `n_components_before`.
#' @export
ConnectKnnComponents <- function(nn_idx, nn_dists, embedding,
                                 fixed_main = FALSE, bridge_candidates = NULL) {
  stopifnot(is.matrix(nn_idx), is.matrix(nn_dists),
            nrow(nn_idx) == nrow(nn_dists))
  emb <- as.matrix(embedding)
  component_labels <- KnnComponentLabels(nn_idx)
  n_before <- max(component_labels)
  if (n_before <= 1L) {
    return(list(nn_idx = nn_idx, nn_dists = nn_dists, n_components_before = n_before))
  }
  comps <- split(seq_len(nrow(nn_idx)), component_labels)

  # The fixed-main option anchors every bridge to the initial largest component.
  sizes <- vapply(comps, length, integer(1))
  main <- comps[[which.max(sizes)]]
  for (ii in seq_along(comps)) {
    if (identical(comps[[ii]], main)) next
    main_cells <- if (is.null(bridge_candidates)) main else
      main[bridge_candidates[main]]
    target_cells <- if (is.null(bridge_candidates)) comps[[ii]] else
      comps[[ii]][bridge_candidates[comps[[ii]]]]
    best <- .nearest_cell_pair(emb, main_cells, target_cells)
    br <- .insert_neighbor(nn_idx, nn_dists, best$i, best$j, best$d)
    nn_idx <- br$nn_idx; nn_dists <- br$nn_dists
    br <- .insert_neighbor(nn_idx, nn_dists, best$j, best$i, best$d)
    nn_idx <- br$nn_idx; nn_dists <- br$nn_dists
    if (!fixed_main) main <- c(main, comps[[ii]])
  }
  list(nn_idx = nn_idx, nn_dists = nn_dists, n_components_before = n_before)
}

KnnComponentCount <- function(nn_idx) {
  idx <- .as_integerish(nn_idx, "nn_idx", allow_na = TRUE)
  n_cells <- nrow(idx)
  k <- ncol(idx)
  usable <- !is.na(idx) & idx != 0L & idx != row(idx) &
    idx >= 1L & idx <= n_cells
  idx[!usable] <- 0L
  as.integer(rs_knn_component_count(
    as.integer(t(idx)), as.integer(n_cells), as.integer(k)
  ))
}

KnnComponentLabels <- function(nn_idx) {
  idx <- .as_integerish(nn_idx, "nn_idx", allow_na = TRUE)
  n_cells <- nrow(idx)
  k <- ncol(idx)
  usable <- !is.na(idx) & idx != 0L & idx != row(idx) &
    idx >= 1L & idx <= n_cells
  idx[!usable] <- 0L
  as.integer(rs_knn_component_labels(
    as.integer(t(idx)), as.integer(n_cells), as.integer(k)
  ))
}

#' @keywords internal
.nearest_cell_pair <- function(emb, idx_a, idx_b) {
  sub_a <- emb[idx_a, , drop = FALSE]
  sub_b <- emb[idx_b, , drop = FALSE]
  nearest <- RANN::nn2(sub_a, sub_b, k = 1L)
  query_id <- which.min(nearest$nn.dists[, 1L])
  list(
    i = idx_a[nearest$nn.idx[query_id, 1L]],
    j = idx_b[query_id],
    d = nearest$nn.dists[query_id, 1L]
  )
}

#' @keywords internal
.insert_neighbor <- function(nn_idx, nn_dists, i, j, d) {
  k <- ncol(nn_idx)
  row_i <- nn_idx[i, ]
  if (j %in% row_i) {
    return(list(nn_idx = nn_idx, nn_dists = nn_dists))
  }
  empty <- which(row_i == 0L)
  slot <- if (length(empty)) empty[1L] else which.max(nn_dists[i, ])
  nn_idx[i, slot] <- j
  nn_dists[i, slot] <- d
  list(nn_idx = nn_idx, nn_dists = nn_dists)
}

#' @keywords internal
.knn_components <- function(nn_idx) {
  n <- nrow(nn_idx)
  adj <- vector("list", n)
  for (i in seq_len(n)) {
    nz <- nn_idx[i, ]
    nz <- nz[nz > 0L & nz <= n]
    if (length(nz)) {
      adj[[i]] <- unique(c(adj[[i]], nz))
      for (j in nz) adj[[j]] <- unique(c(adj[[j]], i))
    }
  }
  seen <- rep(FALSE, n)
  comps <- list()
  for (s in seq_len(n)) {
    if (seen[s]) next
    q <- s
    seen[s] <- TRUE
    comp <- s
    while (length(q)) {
      u <- q[1]
      q <- q[-1]
      for (v in adj[[u]]) {
        if (!seen[v]) {
          seen[v] <- TRUE
          comp <- c(comp, v)
          q <- c(q, v)
        }
      }
    }
    comps[[length(comps) + 1L]] <- comp
  }
  comps
}
