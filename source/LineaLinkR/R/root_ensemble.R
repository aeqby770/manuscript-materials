# Root-direction perturbation ensemble.

#' Infer a root with perturbation support and explicit ambiguity.
#'
#' Candidate root regions are generated without using pseudotime ground truth.
#' Their induced pseudotime orderings are compared across neighbourhood size,
#' representation dimension, and deterministic cluster-stratified subsamples.
#' Selection frequency is a perturbation diagnostic, not a calibrated
#' probability. Symmetric evidence is returned as `"ambiguous"` rather than
#' silently oriented by a coordinate sign.
#'
#' @param embedding finite numeric cell-by-dimension matrix.
#' @param cluster_labels positive integer cluster label per cell.
#' @param nn_idx,nn_dists full-data kNN index and distance matrices.
#' @param fate_entropy optional vector expected to decrease from an early root.
#' @param commitment optional vector expected to increase from an early root.
#' @param potency optional vector expected to decrease from an early root.
#' @param k_grid neighbourhood sizes used for perturbations.
#' @param pc_grid leading embedding dimensions used for perturbations.
#' @param n_resamples number of deterministic stratified subsamples.
#' @param min_resamples minimum subsamples evaluated before adaptive stopping.
#' @param adaptive_resampling whether to stop after `min_resamples` when the
#'   same direction has at least 0.9 perturbation support and clears the
#'   registered direction margin.
#' @param selection_score_mode candidate-ranking score. `"direction_only"`
#'   keeps confidence separate from signed direction evidence;
#'   `"confidence_regularized"` adds coverage and rank stability for an
#'   explicitly ambiguous, annotation-free fallback.
#' @param subsample_fraction fraction retained within every cluster.
#' @param min_direction_margin evidence margin below which direction is
#'   reported as ambiguous.
#' @param seed non-negative integer seed.
#' @return A list containing the selected cell and region, support diagnostics,
#'   pseudotime, candidate table, perturbation table, and evidence mode.
#' @export
InferRootEnsemble <- function(embedding,
                              cluster_labels,
                              nn_idx,
                              nn_dists,
                              fate_entropy = NULL,
                              commitment = NULL,
                              potency = NULL,
                              k_grid = NULL,
                              pc_grid = NULL,
                              n_resamples = 20L,
                              min_resamples = 3L,
                              adaptive_resampling = TRUE,
                              selection_score_mode = c(
                                "direction_only",
                                "confidence_regularized"
                              ),
                              subsample_fraction = 0.8,
                              min_direction_margin = 0.05,
                              seed = 12345L) {
  emb <- .prepare_embedding(embedding)
  cl_map <- .compact_cluster_labels(cluster_labels)
  cl <- cl_map$labels
  knn <- .prepare_knn_inputs(nn_idx, nn_dists)
  n_cells <- nrow(emb)
  if (length(cl) != n_cells || nrow(knn$indices) != n_cells) {
    stop("embedding, cluster labels, and kNN rows must describe the same cells",
         call. = FALSE)
  }

  fate_entropy <- .root_optional_axis(fate_entropy, n_cells, "fate_entropy")
  commitment <- .root_optional_axis(commitment, n_cells, "commitment")
  potency <- .root_optional_axis(potency, n_cells, "potency")
  n_resamples <- .validate_positive_int(n_resamples, "n_resamples")
  min_resamples <- .validate_positive_int(min_resamples, "min_resamples")
  if (min_resamples > n_resamples) {
    stop("`min_resamples` cannot exceed `n_resamples`", call. = FALSE)
  }
  if (length(adaptive_resampling) != 1L || is.na(adaptive_resampling)) {
    stop("`adaptive_resampling` must be TRUE or FALSE", call. = FALSE)
  }
  adaptive_resampling <- isTRUE(adaptive_resampling)
  selection_score_mode <- match.arg(selection_score_mode)
  subsample_fraction <- .validate_unit_interval(
    subsample_fraction, "subsample_fraction"
  )
  if (subsample_fraction <= 0) {
    stop("`subsample_fraction` must be greater than zero", call. = FALSE)
  }
  min_direction_margin <- as.numeric(min_direction_margin)
  if (length(min_direction_margin) != 1L ||
      !is.finite(min_direction_margin) || min_direction_margin < 0) {
    stop("`min_direction_margin` must be a non-negative numeric scalar",
         call. = FALSE)
  }
  seed <- .validate_seed(seed)

  max_k <- min(n_cells - 1L, ncol(knn$indices))
  k_grid <- if (is.null(k_grid)) {
    unique(as.integer(pmax(2L, pmin(max_k, c(10L, 20L, max_k)))))
  } else {
    unique(.as_integerish(k_grid, "k_grid"))
  }
  if (!length(k_grid) || any(k_grid < 2L) || any(k_grid >= n_cells)) {
    stop("`k_grid` must contain integers in [2, n_cells - 1]", call. = FALSE)
  }
  pc_grid <- if (is.null(pc_grid)) {
    unique(as.integer(pmin(ncol(emb), c(2L, 5L, ncol(emb)))))
  } else {
    unique(.as_integerish(pc_grid, "pc_grid"))
  }
  if (!length(pc_grid) || any(pc_grid < 1L) || any(pc_grid > ncol(emb))) {
    stop("`pc_grid` must contain valid leading embedding dimensions",
         call. = FALSE)
  }

  candidate_regions <- .root_candidate_regions(emb, cl)
  evaluate_combinations <- function(combinations) {
    resample_ids <- unique(combinations$resample)
    selected_cache <- setNames(lapply(resample_ids, function(resample_id) {
      .root_stratified_subsample(
        emb, cl, subsample_fraction, resample_id, seed
      )
    }), as.character(resample_ids))
    graph_keys <- unique(combinations[c("pc", "resample")])
    graph_cache <- setNames(lapply(seq_len(nrow(graph_keys)), function(key_id) {
      key <- graph_keys[key_id, ]
      selected <- selected_cache[[as.character(key$resample)]]
      sub_emb <- emb[selected, seq_len(key$pc), drop = FALSE]
      key_k <- combinations$k[
        combinations$pc == key$pc & combinations$resample == key$resample
      ]
      list(
        selected = selected,
        embedding = sub_emb,
        knn_grid = .root_build_knn_grid(
          sub_emb, pmin(key_k, length(selected) - 1L)
        )
      )
    }), paste(graph_keys$pc, graph_keys$resample, sep = ":"))
    perturbations <- lapply(seq_len(nrow(combinations)), function(combo_id) {
      combo <- combinations[combo_id, ]
      graph <- graph_cache[[paste(combo$pc, combo$resample, sep = ":")]]
      selected <- graph$selected
      sub_emb <- graph$embedding
      sub_knn <- graph$knn_grid[[as.character(
        min(combo$k, length(selected) - 1L)
      )]]
      sub_knn <- ConnectKnnComponents(
        sub_knn$indices, sub_knn$distances, sub_emb
      )
      n_components_before <- sub_knn$n_components_before
      sub_knn <- list(
        indices = sub_knn$nn_idx,
        distances = sub_knn$nn_dists
      )

      start_subs <- vapply(candidate_regions, function(region) {
        .root_region_medoid(sub_emb, which(cl[selected] == region))
      }, integer(1L))
      pt_matrix <- ComputePseudotimeFuzzyMulti(
        sub_knn$indices, sub_knn$distances, start_subs,
        rank_normalize = TRUE
      )
      lapply(seq_along(candidate_regions), function(region_id) {
        region <- candidate_regions[[region_id]]
        pt <- pt_matrix[, region_id]
        full_pt <- rep(NA_real_, n_cells)
        full_pt[selected] <- pt
        finite <- is.finite(pt)
        coverage <- mean(finite)
        successful <- coverage >= 0.90 && sum(finite) >= 10L
        data.frame(
          perturbation = combo_id,
          k = combo$k,
          pc = combo$pc,
          resample = combo$resample,
          candidate_region = region,
          n_components_before = n_components_before,
          coverage = coverage,
          entropy_gradient = .root_gradient(
            pt, fate_entropy[selected], -1
          ),
          commitment_gradient = .root_gradient(
            pt, commitment[selected], 1
          ),
          potency_gradient = .root_gradient(pt, potency[selected], -1),
          successful = successful,
          stringsAsFactors = FALSE,
          pt_vector = I(list(full_pt))
        )
      })
    })
    do.call(rbind, unlist(perturbations, recursive = FALSE))
  }

  initial_resamples <- if (adaptive_resampling) {
    min_resamples
  } else {
    n_resamples
  }
  combinations <- expand.grid(
    k = k_grid,
    pc = pc_grid,
    resample = seq_len(initial_resamples),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  perturbation_df <- evaluate_combinations(combinations)
  adaptive_stopped <- FALSE
  if (adaptive_resampling && min_resamples < n_resamples) {
    preview <- perturbation_df
    preview$rank_stability <- 0
    preview <- .root_score_perturbations(
      preview,
      selection_score_mode = selection_score_mode
    )
    preview_valid <- preview[preview$successful, , drop = FALSE]
    preview_scores <- tapply(
      preview_valid$aggregate_score,
      preview_valid$candidate_region,
      mean
    )
    preview_support <- tapply(
      preview_valid$winner_share,
      preview_valid$candidate_region,
      sum
    ) / length(unique(preview_valid$perturbation))
    preview_order <- order(preview_scores, decreasing = TRUE)
    preview_best <- preview_scores[[preview_order[[1L]]]]
    preview_second <- preview_scores[[preview_order[[2L]]]]
    preview_margin <- max(
      0,
      (preview_best - preview_second) /
        max(abs(preview_best), .Machine$double.eps)
    )
    adaptive_stopped <- max(preview_support) >= 0.9 &&
      preview_margin > min_direction_margin
    if (!adaptive_stopped) {
      combinations <- expand.grid(
        k = k_grid,
        pc = pc_grid,
        resample = seq_len(n_resamples),
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
      )
      perturbation_df <- evaluate_combinations(combinations)
    }
  }
  rownames(perturbation_df) <- NULL

  stability <- vapply(candidate_regions, function(region) {
    rows <- perturbation_df$candidate_region == region &
      perturbation_df$successful
    pt_list <- perturbation_df$pt_vector[rows]
    if (length(pt_list) < 2L) return(0)
    pt_matrix <- do.call(cbind, pt_list)
    correlation <- suppressWarnings(stats::cor(
      pt_matrix, method = "spearman", use = "pairwise.complete.obs"
    ))
    values <- abs(correlation[upper.tri(correlation)])
    values <- values[is.finite(values)]
    if (length(values)) stats::median(values) else 0
  }, numeric(1L))
  names(stability) <- as.character(candidate_regions)
  perturbation_df$rank_stability <- unname(
    stability[as.character(perturbation_df$candidate_region)]
  )

  scored <- .root_score_perturbations(
    perturbation_df,
    selection_score_mode = selection_score_mode
  )
  candidate_table <- do.call(rbind, lapply(candidate_regions, function(region) {
    rows <- scored$candidate_region == region & scored$successful
    valid <- scored[rows, , drop = FALSE]
    data.frame(
      candidate_region = as.integer(cl_map$levels[region]),
      selection_frequency = sum(valid$winner_share) /
        length(unique(scored$perturbation[scored$successful])),
      evidence_score = mean(valid$aggregate_score),
      coverage = mean(valid$coverage),
      rank_stability = stability[[as.character(region)]],
      entropy_gradient = .root_mean_finite(valid$entropy_gradient),
      commitment_gradient = .root_mean_finite(valid$commitment_gradient),
      potency_gradient = .root_mean_finite(valid$potency_gradient),
      attractor_direction = .root_mean_finite(valid$attractor_direction),
      confidence_score = .root_mean_finite(valid$confidence_score),
      direction_component_count = max(valid$direction_component_count),
      n_successful = nrow(valid),
      stringsAsFactors = FALSE
    )
  }))
  compact_lookup <- match(candidate_table$candidate_region, cl_map$levels)
  candidate_table$compact_region <- compact_lookup
  ordering <- order(
    candidate_table$evidence_score,
    candidate_table$selection_frequency,
    -candidate_table$candidate_region,
    decreasing = TRUE
  )
  candidate_table <- candidate_table[ordering, , drop = FALSE]
  best <- candidate_table[1L, , drop = FALSE]
  second_score <- if (nrow(candidate_table) > 1L) {
    candidate_table$evidence_score[[2L]]
  } else {
    best$evidence_score[[1L]]
  }
  direction_margin <- max(
    0,
    (best$evidence_score[[1L]] - second_score) /
      max(abs(best$evidence_score[[1L]]), .Machine$double.eps)
  )
  root_support <- best$selection_frequency[[1L]]
  effective_direction_components <- as.integer(1L + !is.null(potency))
  direction_status <- if (effective_direction_components < 2L ||
                          direction_margin <= min_direction_margin ||
                          root_support < 0.5) {
    "ambiguous"
  } else if (root_support < 0.7) {
    "weak"
  } else {
    "supported"
  }

  best_compact <- best$compact_region[[1L]]
  full_region_cells <- which(cl == best_compact)
  start_cell <- .root_region_medoid(emb, full_region_cells)
  pseudotime <- ComputePseudotimeFuzzy(
    knn$indices, knn$distances, start_cell, rank_normalize = TRUE
  )
  scored$pt_vector <- NULL
  candidate_table$compact_region <- NULL

  list(
    start_cell = as.integer(start_cell),
    root_region = as.integer(best$candidate_region[[1L]]),
    root_support = root_support,
    direction_margin = direction_margin,
    direction_status = direction_status,
    pseudotime = pseudotime,
    candidates = candidate_table,
    perturbations = scored,
    direction_score_components = c(
      "attractor_direction",
      if (!is.null(potency)) "potency_gradient"
    ),
    confidence_components = c("coverage", "rank_stability"),
    selection_score_mode = selection_score_mode,
    effective_direction_components = effective_direction_components,
    adaptive_stopped = adaptive_stopped,
    resamples_used = max(perturbation_df$resample),
    evidence_mode = if (is.null(potency)) "blind" else "anchored"
  )
}

.root_optional_axis <- function(x, n_cells, name) {
  if (is.null(x)) return(NULL)
  value <- as.numeric(x)
  if (length(value) != n_cells || any(is.infinite(value))) {
    stop(sprintf("`%s` must have one finite-or-NA value per cell", name),
         call. = FALSE)
  }
  value
}

.root_candidate_regions <- function(emb, cl) {
  pc1 <- .pc1_scores_signed(emb)
  cluster_ids <- sort(unique(cl))
  cluster_pc1 <- vapply(cluster_ids, function(region) {
    mean(pc1[cl == region])
  }, numeric(1L))
  arbitrary_root <- cluster_ids[[1L]]
  lineages <- ExtractLineages(emb, cl, start_cluster = arbitrary_root)
  edges <- do.call(rbind, lapply(lineages, function(path) {
    if (length(path) < 2L) return(NULL)
    cbind(path[-length(path)], path[-1L])
  }))
  mst_leaves <- if (is.null(edges)) {
    integer(0)
  } else {
    degree <- table(c(edges[, 1L], edges[, 2L]))
    as.integer(names(degree)[degree == 1L])
  }
  sort(unique(c(
    cluster_ids[[which.min(cluster_pc1)]],
    cluster_ids[[which.max(cluster_pc1)]],
    mst_leaves
  )))
}

.root_stratified_subsample <- function(emb, cl, fraction, resample, seed) {
  global <- colMeans(emb)
  selected <- unlist(lapply(sort(unique(cl)), function(region) {
    cells <- which(cl == region)
    centroid <- colMeans(emb[cells, , drop = FALSE])
    radial <- sqrt(rowSums((emb[cells, , drop = FALSE] - matrix(
      global, nrow = length(cells), ncol = ncol(emb), byrow = TRUE
    ))^2))
    local <- sqrt(rowSums((emb[cells, , drop = FALSE] - matrix(
      centroid, nrow = length(cells), ncol = ncol(emb), byrow = TRUE
    ))^2))
    canonical <- cells[order(radial, local)]
    retain <- min(length(cells), max(2L, floor(length(cells) * fraction)))
    set.seed(as.integer((seed + 104729 * resample) %%
                          .Machine$integer.max))
    sort(canonical[sample.int(length(canonical), retain, replace = FALSE)])
  }), use.names = FALSE)
  sort(as.integer(selected))
}

.root_build_knn <- function(embedding, k) {
  k <- min(as.integer(k), nrow(embedding) - 1L)
  raw <- RANN::nn2(embedding, embedding, k = k + 1L)
  indices <- t(vapply(seq_len(nrow(embedding)), function(i) {
    keep <- which(raw$nn.idx[i, ] != i)[seq_len(k)]
    as.integer(raw$nn.idx[i, keep])
  }, integer(k)))
  distances <- t(vapply(seq_len(nrow(embedding)), function(i) {
    keep <- which(raw$nn.idx[i, ] != i)[seq_len(k)]
    as.numeric(raw$nn.dists[i, keep])
  }, numeric(k)))
  list(indices = indices, distances = distances)
}

.root_build_knn_grid <- function(embedding, k_grid) {
  k_grid <- sort(unique(as.integer(k_grid)))
  max_k <- min(max(k_grid), nrow(embedding) - 1L)
  full <- .root_build_knn(embedding, max_k)
  setNames(lapply(k_grid, function(k) {
    k <- min(k, max_k)
    list(
      indices = full$indices[, seq_len(k), drop = FALSE],
      distances = full$distances[, seq_len(k), drop = FALSE]
    )
  }), as.character(k_grid))
}

.root_region_medoid <- function(embedding, region_cells) {
  centroid <- colMeans(embedding[region_cells, , drop = FALSE])
  distance <- rowSums((embedding[region_cells, , drop = FALSE] - matrix(
    centroid, nrow = length(region_cells), ncol = ncol(embedding), byrow = TRUE
  ))^2)
  as.integer(region_cells[[which.min(distance)]])
}

.root_gradient <- function(pseudotime, axis, expected_sign) {
  if (is.null(axis)) return(NA_real_)
  valid <- is.finite(pseudotime) & is.finite(axis)
  if (sum(valid) < 10L) return(NA_real_)
  value <- suppressWarnings(stats::cor(
    pseudotime[valid], axis[valid], method = "spearman"
  ))
  if (is.finite(value)) expected_sign * value else NA_real_
}

.root_component_rank <- function(x) {
  finite <- is.finite(x)
  out <- rep(NA_real_, length(x))
  if (!any(finite)) return(out)
  if (diff(range(x[finite])) <= .Machine$double.eps) {
    out[finite] <- 0.5
  } else {
    out[finite] <- (rank(x[finite], ties.method = "average") - 1) /
      (sum(finite) - 1)
  }
  out
}

.root_score_perturbations <- function(
    perturbation_df,
    selection_score_mode = c("direction_only", "confidence_regularized")) {
  selection_score_mode <- match.arg(selection_score_mode)
  split_rows <- split(seq_len(nrow(perturbation_df)), perturbation_df$perturbation)
  scored <- lapply(split_rows, function(rows) {
    block <- perturbation_df[rows, , drop = FALSE]
    attractor_matrix <- cbind(
      entropy_gradient = .root_gradient_evidence(block$entropy_gradient),
      commitment_gradient = .root_gradient_evidence(block$commitment_gradient)
    )
    attractor_count <- rowSums(is.finite(attractor_matrix))
    block$attractor_direction <- rowMeans(attractor_matrix, na.rm = TRUE)
    block$attractor_direction[attractor_count == 0L] <- NA_real_
    direction_matrix <- cbind(
      attractor_direction = block$attractor_direction,
      potency_gradient = .root_gradient_evidence(block$potency_gradient)
    )
    block$direction_component_count <- rowSums(is.finite(direction_matrix))
    block$aggregate_score <- rowMeans(direction_matrix, na.rm = TRUE)
    block$aggregate_score[block$direction_component_count == 0L] <- NA_real_
    block$confidence_score <- rowMeans(cbind(
      coverage = block$coverage,
      rank_stability = block$rank_stability
    ))
    confidence_regularized_matrix <- cbind(
      coverage = block$coverage,
      rank_stability = block$rank_stability,
      entropy_gradient = attractor_matrix[, "entropy_gradient"],
      commitment_gradient = attractor_matrix[, "commitment_gradient"],
      potency_gradient = direction_matrix[, "potency_gradient"]
    )
    block$confidence_regularized_score <- rowMeans(
      confidence_regularized_matrix,
      na.rm = TRUE
    )
    block$aggregate_score <- switch(
      selection_score_mode,
      direction_only = block$aggregate_score,
      confidence_regularized = block$confidence_regularized_score
    )
    block$selection_score_mode <- selection_score_mode
    block$aggregate_score[!block$successful] <- NA_real_
    if (!any(is.finite(block$aggregate_score))) {
      block$winner_share <- 0
      return(block)
    }
    best <- max(block$aggregate_score, na.rm = TRUE)
    tied <- is.finite(block$aggregate_score) &
      abs(block$aggregate_score - best) <= 1e-12
    block$winner_share <- ifelse(tied, 1 / sum(tied), 0)
    block
  })
  do.call(rbind, scored)
}

.root_gradient_evidence <- function(x) {
  finite <- is.finite(x)
  out <- rep(NA_real_, length(x))
  out[finite] <- (pmax(-1, pmin(1, x[finite])) + 1) / 2
  out
}

.root_mean_finite <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) mean(x) else NA_real_
}
