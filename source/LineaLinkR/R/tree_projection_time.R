#' Unified pseudotime from projection onto an existing rooted lineage tree
#'
#' Keeps the supplied root, cluster labels and root-to-leaf topology fixed.
#' Cluster nodes are placed at their cell centroids in `embedding`, except that
#' the root node is placed at the supplied root cell. Each cell is projected
#' only onto tree edges incident to its assigned cluster. Distances are measured
#' in the same embedding used to define the nodes, and the whole tree shares one
#' longest-root-to-leaf normalization.
#'
#' @param embedding Numeric matrix with cells in rows.
#' @param cluster_labels Cluster identifier for each cell.
#' @param lineage_paths List of root-to-leaf cluster paths.
#' @param start_cell One-based index of the fixed root cell.
#' @param lineage_weights Optional cell-by-lineage weight matrix. Required when
#'   `edge_constraint = "dominant_lineage_path"`.
#' @param edge_constraint Allowed edge set for each cell. The original
#'   candidate uses edges incident to the assigned cluster. The representation
#'   correction uses all edges on the cell's fixed maximum-weight lineage.
#' @return A list containing `pseudotime`, per-cell `projection`, rooted
#'   `tree_edges`, `node_table`, and `longest_root_to_leaf`.
#' @export
ComputeUnifiedTreeProjectionTime <- function(
    embedding, cluster_labels, lineage_paths, start_cell,
    lineage_weights = NULL,
    edge_constraint = c("incident_cluster", "dominant_lineage_path")) {
  embedding <- as.matrix(embedding)
  storage.mode(embedding) <- "double"
  cluster_labels <- as.integer(cluster_labels)
  start_cell <- as.integer(start_cell)
  edge_constraint <- match.arg(edge_constraint)
  if (nrow(embedding) != length(cluster_labels)) {
    stop("`cluster_labels` length must equal nrow(embedding)", call. = FALSE)
  }
  if (length(start_cell) != 1L || is.na(start_cell) ||
      start_cell < 1L || start_cell > nrow(embedding)) {
    stop("`start_cell` must index one row of `embedding`", call. = FALSE)
  }
  if (!is.list(lineage_paths) || !length(lineage_paths)) {
    stop("`lineage_paths` must be a non-empty list", call. = FALSE)
  }

  root_cluster <- cluster_labels[start_cell]
  paths <- lapply(lineage_paths, as.integer)
  roots <- vapply(paths, `[`, integer(1L), 1L)
  if (any(lengths(paths) < 2L) || any(roots != root_cluster)) {
    stop("Each lineage path must start at the fixed root and contain an edge",
         call. = FALSE)
  }

  node_sums <- rowsum(embedding, cluster_labels, reorder = TRUE)
  node_ids <- as.integer(rownames(node_sums))
  node_counts <- tabulate(
    match(cluster_labels, node_ids), nbins = length(node_ids)
  )
  node_centers <- node_sums / node_counts
  node_centers[match(root_cluster, node_ids), ] <- embedding[start_cell, ]

  tree_edges <- unique(data.table::rbindlist(lapply(paths, function(path) {
    data.table::data.table(parent = head(path, -1L), child = tail(path, -1L))
  })))
  tree_edges[, `:=`(
    tie_node_low = pmin(parent, child),
    tie_node_high = pmax(parent, child)
  )]
  data.table::setorder(
    tree_edges, tie_node_low, tie_node_high, parent, child
  )
  tree_edges[, edge_id := seq_len(.N)]
  parent_coordinates <- node_centers[
    match(tree_edges$parent, node_ids), , drop = FALSE
  ]
  child_coordinates <- node_centers[
    match(tree_edges$child, node_ids), , drop = FALSE
  ]
  tree_edges[, edge_length := sqrt(rowSums(
    (child_coordinates - parent_coordinates)^2
  ))]
  tree_edges[, edge_key := paste(parent, child, sep = ">")]

  path_node_distances <- data.table::rbindlist(lapply(paths, function(path) {
    path_keys <- paste(head(path, -1L), tail(path, -1L), sep = ">")
    path_lengths <- tree_edges$edge_length[
      match(path_keys, tree_edges$edge_key)
    ]
    data.table::data.table(
      node = path,
      root_distance = c(0, cumsum(path_lengths))
    )
  }))
  node_table <- path_node_distances[, .(
    root_distance = mean(root_distance),
    path_distance_spread = diff(range(root_distance))
  ), by = node]
  child_parent_counts <- tree_edges[, data.table::uniqueN(parent), by = child]
  if (!setequal(node_ids, node_table$node) ||
      max(node_table$path_distance_spread) > 1e-10 ||
      root_cluster %in% tree_edges$child ||
      any(tree_edges$child == tree_edges$parent) ||
      any(child_parent_counts$V1 != 1L)) {
    stop("`lineage_paths` do not define one consistent rooted tree",
         call. = FALSE)
  }

  leaf_nodes <- vapply(paths, tail, n = 1L, FUN.VALUE = integer(1L))
  longest_root_to_leaf <- max(
    node_table[node %in% leaf_nodes, root_distance]
  )
  if (!is.finite(longest_root_to_leaf) || longest_root_to_leaf <= 0) {
    stop("The rooted tree must have positive total length", call. = FALSE)
  }
  tree_edges[
    node_table, on = .(parent = node),
    parent_root_distance := i.root_distance
  ]
  tree_edges[
    node_table, on = .(child = node),
    child_root_distance := i.root_distance
  ]

  cells <- data.table::data.table(
    cell_index = seq_along(cluster_labels),
    cluster = cluster_labels
  )
  if (edge_constraint == "incident_cluster") {
    outgoing <- merge(
      cells, tree_edges,
      by.x = "cluster", by.y = "parent", allow.cartesian = TRUE
    )
    data.table::setnames(outgoing, "cluster", "parent")
    incoming <- merge(
      cells, tree_edges,
      by.x = "cluster", by.y = "child", allow.cartesian = TRUE
    )
    data.table::setnames(incoming, "cluster", "child")
    candidates <- unique(data.table::rbindlist(list(
      outgoing[, .(
        cell_index, parent, child, edge_id, edge_length,
        parent_root_distance, tie_node_low, tie_node_high
      )],
      incoming[, .(
        cell_index, parent, child, edge_id, edge_length,
        parent_root_distance, tie_node_low, tie_node_high
      )]
    )), by = c("cell_index", "edge_id"))
  } else {
    lineage_weights <- as.matrix(lineage_weights)
    if (nrow(lineage_weights) != nrow(embedding) ||
        ncol(lineage_weights) != length(paths) ||
        any(!is.finite(lineage_weights))) {
      stop("`lineage_weights` must match cells and lineage paths",
           call. = FALSE)
    }
    cells[, dominant_lineage := max.col(
      lineage_weights, ties.method = "first"
    )]
    path_edges <- unique(data.table::rbindlist(lapply(
      seq_along(paths),
      function(lineage) {
        path <- paths[[lineage]]
        data.table::data.table(
          dominant_lineage = lineage,
          parent = head(path, -1L),
          child = tail(path, -1L)
        )
      }
    )))
    path_edges <- merge(
      path_edges, tree_edges,
      by = c("parent", "child"), all.x = TRUE, sort = FALSE
    )
    candidates <- merge(
      cells[, .(cell_index, dominant_lineage)],
      path_edges,
      by = "dominant_lineage", allow.cartesian = TRUE
    )[, .(
      cell_index, parent, child, edge_id, edge_length,
      parent_root_distance, tie_node_low, tie_node_high
    )]
  }
  if (data.table::uniqueN(candidates$cell_index) != nrow(embedding)) {
    stop("At least one cluster has no incident lineage-tree edge",
         call. = FALSE)
  }

  cell_coordinates <- embedding[candidates$cell_index, , drop = FALSE]
  parent_coordinates <- node_centers[
    match(candidates$parent, node_ids), , drop = FALSE
  ]
  child_coordinates <- node_centers[
    match(candidates$child, node_ids), , drop = FALSE
  ]
  edge_vectors <- child_coordinates - parent_coordinates
  alpha_raw <- rowSums(
    (cell_coordinates - parent_coordinates) * edge_vectors
  ) / rowSums(edge_vectors^2)
  alpha_raw[!is.finite(alpha_raw)] <- 0
  alpha <- pmin(1, pmax(0, alpha_raw))
  projected_coordinates <- parent_coordinates + edge_vectors * alpha
  candidates[, `:=`(
    alpha = alpha,
    projection_distance_sq = rowSums(
      (cell_coordinates - projected_coordinates)^2
    )
  )]
  data.table::setorder(
    candidates, cell_index, projection_distance_sq,
    tie_node_low, tie_node_high, edge_id
  )
  projection <- candidates[!duplicated(cell_index)]
  data.table::setorder(projection, cell_index)
  pseudotime <- pmin(1, pmax(0, (
    projection$parent_root_distance +
      projection$alpha * projection$edge_length
  ) / longest_root_to_leaf))

  projection[, `:=`(
    cluster = cluster_labels,
    projection_distance = sqrt(projection_distance_sq),
    pseudotime = pseudotime
  )]
  node_table[, `:=`(
    is_root = node == root_cluster,
    is_leaf = node %in% leaf_nodes
  )]
  list(
    pseudotime = pseudotime,
    projection = projection[, .(
      cell_index, cluster, parent, child, edge_id, alpha,
      projection_distance, pseudotime
    )],
    tree_edges = tree_edges[, .(
      edge_id, parent, child, edge_length,
      parent_root_distance, child_root_distance,
      tie_node_low, tie_node_high
    )],
    node_table = node_table,
    root_cluster = root_cluster,
    longest_root_to_leaf = longest_root_to_leaf,
    edge_constraint = edge_constraint
  )
}
