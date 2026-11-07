# Blind root-cell inference for trajectory pipelines (no pseudotime GT root protocol).
#
# Does NOT read pseudotime ground-truth columns from metadata — callers must
# never pass palantir_pseudotime / true_pseudotime into this function.

#' Infer a root cell without using pseudotime ground truth.
#'
#' @param embedding numeric matrix \[n_cells x n_dims\] (UMAP or PCA).
#' @param cluster_labels integer vector of length n_cells, 1-based cluster ids.
#' @param method character; `"mst_extremum"`, `"cluster_pc1_min"`, `"cluster_pc1_max"`, or `"auto"`.
#' @param seed integer seed for tie-breaking.
#' @return integer scalar, 1-based cell index.
#' @export
InferStartCell <- function(embedding,
                           cluster_labels,
                           method = c("mst_extremum", "cluster_pc1_min", "cluster_pc1_max", "auto"),
                           seed = 1L) {
  method <- match.arg(method)
  emb <- .prepare_embedding(embedding)
  cl <- .as_integerish(cluster_labels, "cluster_labels")
  n_cells <- nrow(emb)
  if (n_cells < 2L || length(cl) != n_cells || any(cl < 1L)) {
    stop("`cluster_labels` must contain one positive integer per embedding row",
         call. = FALSE)
  }
  seed <- .validate_non_negative_int(seed, "seed")

  if (method == "cluster_pc1_min") {
    return(.infer_start_cluster_pc1_min(emb, cl, seed))
  }
  if (method == "cluster_pc1_max") {
    return(.infer_start_cluster_pc1_max(emb, cl, seed))
  }
  if (method == "auto") {
    stop("InferStartCell(..., method='auto') requires nn graph; use InferStartCellAuto()")
  }

  .infer_start_mst_extremum(emb, cl, seed)
}

#' Blind score: mean Spearman(pseudotime, sign-anchored PC1..PCk). Higher is better.
#' @param pseudotime numeric vector.
#' @param embedding numeric matrix \[n_cells x n_dims\].
#' @param n_pc number of PCs to average (default 3).
#' @return numeric scalar; -Inf if undefined.
#' @export
BlindEmbeddingScore <- function(pseudotime, embedding, n_pc = 3L) {
  emb <- .prepare_embedding(embedding)
  if (!is.numeric(pseudotime) || length(pseudotime) != nrow(emb) ||
      any(is.infinite(pseudotime))) {
    stop("`pseudotime` length must equal nrow(embedding)", call. = FALSE)
  }
  n_pc <- .validate_positive_int(n_pc, "n_pc")
  vv <- is.finite(pseudotime)
  if (sum(vv) < 50L || ncol(emb) < 1L) return(-Inf)
  nr <- min(n_pc, ncol(emb), sum(vv) - 1L)
  if (nr < 1L) return(-Inf)
  pcs_obj <- tryCatch(
    stats::prcomp(emb[vv, , drop = FALSE], center = TRUE, rank. = nr),
    error = function(e) NULL
  )
  if (is.null(pcs_obj) || !ncol(pcs_obj$x)) return(-Inf)
  pcs <- pcs_obj$x
  rot <- pcs_obj$rotation[, seq_len(ncol(pcs)), drop = FALSE]
  anchors <- max.col(t(abs(rot)), ties.method = "first")
  signs <- sign(rot[cbind(anchors, seq_along(anchors))])
  signs[signs == 0] <- 1
  pcs <- sweep(pcs, 2L, signs, `*`)
  cors <- vapply(seq_len(ncol(pcs)), function(j) {
    suppressWarnings(stats::cor(pseudotime[vv], pcs[, j], method = "spearman"))
  }, numeric(1))
  cors <- cors[is.finite(cors)]
  if (!length(cors)) return(-Inf)
  mean(cors)
}

#' Direction-invariant ordering score in expression-derived embedding space.
#'
#' Unlike `BlindEmbeddingScore()`, this score does not assign a biological
#' direction from an arbitrary principal-component sign. It measures the root
#' mean square of absolute Spearman correlations between pseudotime and the
#' leading embedding PCs. Reversing pseudotime or any embedding axis therefore
#' leaves the score unchanged.
#'
#' @param pseudotime numeric vector.
#' @param embedding numeric matrix \[n_cells x n_dims\].
#' @param n_pc number of PCs used (default 3).
#' @return numeric scalar in \[0, 1\]; `-Inf` if undefined.
#' @export
BlindOrderingScore <- function(pseudotime, embedding, n_pc = 3L) {
  emb <- .prepare_embedding(embedding)
  if (!is.numeric(pseudotime) || length(pseudotime) != nrow(emb) ||
      any(is.infinite(pseudotime))) {
    stop("`pseudotime` length must equal nrow(embedding)", call. = FALSE)
  }
  n_pc <- .validate_positive_int(n_pc, "n_pc")
  valid <- is.finite(pseudotime)
  if (sum(valid) < 50L || ncol(emb) < 1L) return(-Inf)
  rank_limit <- min(n_pc, ncol(emb), sum(valid) - 1L)
  if (rank_limit < 1L) return(-Inf)

  pcs <- stats::prcomp(
    emb[valid, , drop = FALSE],
    center = TRUE,
    rank. = rank_limit
  )$x
  if (!ncol(pcs)) return(-Inf)

  pseudotime_rank <- rank(pseudotime[valid], ties.method = "average")
  pseudotime_centered <- pseudotime_rank - mean(pseudotime_rank)
  pc_ranks <- apply(pcs, 2L, rank, ties.method = "average")
  pc_ranks <- matrix(pc_ranks, nrow = nrow(pcs))
  pc_centered <- sweep(pc_ranks, 2L, colMeans(pc_ranks), `-`)
  numerator <- colSums(pc_centered * pseudotime_centered)
  denominator <- sqrt(
    sum(pseudotime_centered^2) * colSums(pc_centered^2)
  )
  correlations <- abs(numerator / denominator)
  correlations <- correlations[is.finite(correlations)]
  if (!length(correlations)) return(-Inf)
  sqrt(mean(correlations^2))
}

#' Select a root candidate only when blind evidence strongly beats the default.
#'
#' The high-potency state-aggregate candidate remains the default. An
#' alternative is accepted only when its blind embedding score is positive and
#' exceeds the default by at least `min_gain`. This conservative gate avoids
#' replacing a biological potency prior for small changes in a sign-convention
#' diagnostic. Ground-truth time, lineage, clone and fate labels are not inputs.
#'
#' @param candidates data frame with `candidate`, `blind_score`, and
#'   `component_score` columns.
#' @param default_candidate character scalar naming the default row.
#' @param min_blind_score minimum blind score required for a fallback.
#' @param min_gain minimum improvement over the default blind score.
#' @return list describing the selected candidate and blind-score gain.
#' @export
SelectRootCandidateConservative <- function(candidates,
                                            default_candidate =
                                              "detected_genes_high_state",
                                            min_blind_score = 0,
                                            min_gain = 0.30) {
  cand <- as.data.frame(candidates, stringsAsFactors = FALSE)
  required <- c("candidate", "blind_score", "component_score")
  stopifnot(
    all(required %in% names(cand)),
    nrow(cand) >= 1L,
    length(default_candidate) == 1L,
    is.finite(min_blind_score),
    is.finite(min_gain),
    min_gain >= 0
  )
  cand$candidate <- as.character(cand$candidate)
  cand$blind_score <- as.numeric(cand$blind_score)
  cand$component_score <- as.numeric(cand$component_score)
  default_index <- match(default_candidate, cand$candidate)
  stopifnot(length(default_index) == 1L, is.finite(default_index))
  default_score <- cand$blind_score[[default_index]]
  stopifnot(is.finite(default_score))

  cand$blind_gain <- cand$blind_score - default_score
  cand$eligible <- cand$candidate != default_candidate &
    is.finite(cand$blind_score) &
    cand$blind_score > min_blind_score &
    cand$blind_gain >= min_gain
  eligible_order <- order(
    !cand$eligible,
    -cand$blind_score,
    -cand$component_score,
    seq_len(nrow(cand))
  )
  selected_index <- c(
    eligible_order[cand$eligible[eligible_order]],
    default_index
  )[[1L]]
  selected_is_default <- selected_index == default_index

  list(
    selected_candidate = cand$candidate[[selected_index]],
    default_candidate = default_candidate,
    selected_blind_score = cand$blind_score[[selected_index]],
    default_blind_score = default_score,
    blind_gain = cand$blind_gain[[selected_index]],
    selection_status = ifelse(
      selected_is_default,
      "default_retained",
      "strong_blind_fallback"
    ),
    min_blind_score = as.numeric(min_blind_score),
    min_gain = as.numeric(min_gain)
  )
}

#' Select a high-potency root representation by ordering quality.
#'
#' Candidate ordering is compared with `BlindOrderingScore()`, which is
#' invariant to trajectory direction and embedding signs. The caller supplies
#' the allowed candidate names; production annotation-free inference restricts
#' this set to high-potency state and high-potency cellwise roots, so a
#' low-potency direction cannot be introduced automatically.
#'
#' @param candidates data frame with `candidate`, `ordering_score`, and
#'   `component_score` columns.
#' @param default_candidate character scalar naming the default row.
#' @param allowed_candidates character vector of candidates eligible for
#'   automatic selection.
#' @param min_gain minimum ordering-score improvement over the default.
#' @return list describing the selected candidate and score gain.
#' @export
SelectRootCandidateByOrdering <- function(candidates,
                                          default_candidate =
                                            "detected_genes_high_state",
                                          allowed_candidates = c(
                                            "detected_genes_high_state",
                                            "detected_genes_high_cellwise"
                                          ),
                                          min_gain = 0.10) {
  cand <- as.data.frame(candidates, stringsAsFactors = FALSE)
  required <- c("candidate", "ordering_score", "component_score")
  stopifnot(
    all(required %in% names(cand)),
    nrow(cand) >= 1L,
    length(default_candidate) == 1L,
    length(allowed_candidates) >= 1L,
    is.finite(min_gain),
    min_gain >= 0
  )
  cand$candidate <- as.character(cand$candidate)
  cand$ordering_score <- as.numeric(cand$ordering_score)
  cand$component_score <- as.numeric(cand$component_score)
  eligible_rows <- cand$candidate %in% unique(allowed_candidates)
  stopifnot(
    default_candidate %in% cand$candidate,
    default_candidate %in% allowed_candidates,
    any(eligible_rows)
  )
  default_index <- match(default_candidate, cand$candidate)
  default_score <- cand$ordering_score[[default_index]]
  stopifnot(is.finite(default_score))

  cand$ordering_gain <- cand$ordering_score - default_score
  cand$eligible <- eligible_rows &
    cand$candidate != default_candidate &
    is.finite(cand$ordering_score) &
    cand$ordering_gain >= min_gain
  eligible_order <- order(
    !cand$eligible,
    -cand$ordering_score,
    -cand$component_score,
    seq_len(nrow(cand))
  )
  selected_index <- c(
    eligible_order[cand$eligible[eligible_order]],
    default_index
  )[[1L]]
  selected_is_default <- selected_index == default_index

  list(
    selected_candidate = cand$candidate[[selected_index]],
    default_candidate = default_candidate,
    selected_ordering_score = cand$ordering_score[[selected_index]],
    default_ordering_score = default_score,
    ordering_gain = cand$ordering_gain[[selected_index]],
    selection_status = ifelse(
      selected_is_default,
      "default_retained",
      "ordering_fallback"
    ),
    allowed_candidates = unique(as.character(allowed_candidates)),
    min_gain = as.numeric(min_gain)
  )
}

#' Select a bidirectional potency endpoint without ground-truth time.
#'
#' The high-potency state root remains the default. A different potency mode is
#' accepted only when its state-aggregate pseudotime has a Spearman correlation
#' with potency that is at least `min_alignment_gain` more negative than the
#' default. Within the selected potency mode, the cellwise representation is
#' accepted only when its direction-invariant ordering score improves by at
#' least `min_ordering_gain` and its negative potency alignment improves by at
#' least `min_representation_alignment_gain`. Ground-truth time, lineage, clone
#' and fate labels are not inputs.
#'
#' @param candidates data frame with `candidate`, `potency_mode`,
#'   `root_representation`, `potency_spearman`, `ordering_score`, and
#'   `component_score` columns. An optional logical
#'   `terminal_leaf_eligible` column marks candidates whose oriented anchor is
#'   a terminal leaf in the expression-derived lineage topology. Each potency
#'   mode must contain one
#'   `state_aggregate` row; an optional `cellwise_component` row may accompany
#'   it. An optional `available` column excludes undefined candidates before
#'   score comparisons. If the default mode is unavailable, only an available
#'   terminal-leaf mode can be selected; no alignment gain is reported.
#' @param default_mode character scalar naming the default potency mode.
#' @param min_alignment_gain minimum improvement in negative potency
#'   correlation required to switch potency modes.
#' @param min_ordering_gain minimum ordering-score improvement required to use
#'   the cellwise representation within the selected mode.
#' @param min_representation_alignment_gain minimum improvement in negative
#'   potency correlation required to use the cellwise representation. Zero
#'   prevents a geometry-only switch when potency alignment becomes weaker.
#' @return list describing the selected endpoint mode and representation.
#' @export
SelectRootEndpointByPotencyAlignment <- function(
    candidates,
    default_mode = "detected_genes_high",
    min_alignment_gain = 0.20,
    min_ordering_gain = 0.10,
    min_representation_alignment_gain = 0.00) {
  cand <- as.data.frame(candidates, stringsAsFactors = FALSE)
  unavailable_candidates <- character()
  if ("available" %in% names(cand)) {
    stopifnot(is.logical(cand$available), !anyNA(cand$available))
    unavailable_candidates <- as.character(cand$candidate[!cand$available])
    cand <- cand[cand$available, , drop = FALSE]
  }
  if (!nrow(cand)) stop("No available root candidate", call. = FALSE)
  required <- c(
    "candidate",
    "potency_mode",
    "root_representation",
    "potency_spearman",
    "ordering_score",
    "component_score"
  )
  stopifnot(
    all(required %in% names(cand)),
    nrow(cand) >= 1L,
    length(default_mode) == 1L,
    is.finite(min_alignment_gain),
    min_alignment_gain >= 0,
    is.finite(min_ordering_gain),
    min_ordering_gain >= 0,
    is.finite(min_representation_alignment_gain),
    min_representation_alignment_gain >= 0
  )
  cand$candidate <- as.character(cand$candidate)
  cand$potency_mode <- as.character(cand$potency_mode)
  cand$root_representation <- as.character(cand$root_representation)
  cand$potency_spearman <- as.numeric(cand$potency_spearman)
  cand$ordering_score <- as.numeric(cand$ordering_score)
  cand$component_score <- as.numeric(cand$component_score)
  if (!"terminal_leaf_eligible" %in% names(cand)) {
    cand$terminal_leaf_eligible <- rep(TRUE, nrow(cand))
  }
  cand$terminal_leaf_eligible <- as.logical(cand$terminal_leaf_eligible)
  stopifnot(
    !anyDuplicated(cand$candidate),
    all(is.finite(cand$potency_spearman)),
    all(is.finite(cand$ordering_score)),
    all(is.finite(cand$component_score)),
    !anyNA(cand$terminal_leaf_eligible)
  )

  state_rows <- cand$root_representation == "state_aggregate"
  state_candidates <- cand[state_rows, , drop = FALSE]
  stopifnot(
    nrow(state_candidates) >= 1L,
    !anyDuplicated(state_candidates$potency_mode)
  )
  default_state_index <- match(default_mode, state_candidates$potency_mode)
  default_available <- !is.na(default_state_index)
  default_alignment <- state_candidates$potency_spearman[default_state_index]
  unconstrained_endpoint_order <- order(
    state_candidates$potency_spearman,
    -state_candidates$component_score,
    state_candidates$candidate
  )
  unconstrained_best_state_index <- unconstrained_endpoint_order[[1L]]
  unconstrained_best_mode <- state_candidates$potency_mode[[
    unconstrained_best_state_index
  ]]
  endpoint_eligible <-
    state_candidates$potency_mode == default_mode |
    state_candidates$terminal_leaf_eligible
  if (!any(endpoint_eligible)) {
    stop("No available root candidate satisfies the terminal-leaf criterion", call. = FALSE)
  }
  endpoint_order <- order(
    !endpoint_eligible,
    state_candidates$potency_spearman,
    -state_candidates$component_score,
    state_candidates$candidate
  )
  best_state_index <- endpoint_order[[1L]]
  best_mode <- state_candidates$potency_mode[[best_state_index]]
  endpoint_alignment_gain <- default_alignment -
    state_candidates$potency_spearman[[best_state_index]]
  endpoint_switch <- best_mode != default_mode &&
    (!default_available || endpoint_alignment_gain >= min_alignment_gain)
  selected_mode <- ifelse(endpoint_switch, best_mode, default_mode)

  mode_candidates <- cand[
    cand$potency_mode == selected_mode,
    ,
    drop = FALSE
  ]
  selected_state_index <- which(
    mode_candidates$root_representation == "state_aggregate"
  )[[1L]]
  state_ordering_score <- mode_candidates$ordering_score[[
    selected_state_index
  ]]
  representation_gain <- mode_candidates$ordering_score -
    state_ordering_score
  state_potency_spearman <- mode_candidates$potency_spearman[[
    selected_state_index
  ]]
  representation_alignment_gain <- state_potency_spearman -
    mode_candidates$potency_spearman
  representation_eligible <-
    mode_candidates$root_representation == "state_aggregate" |
    (
      mode_candidates$root_representation == "cellwise_component" &
        (
          selected_mode == default_mode |
          mode_candidates$terminal_leaf_eligible
        ) &
        representation_gain >= min_ordering_gain &
        representation_alignment_gain >=
          min_representation_alignment_gain
    )
  representation_order <- order(
    !representation_eligible,
    -mode_candidates$ordering_score,
    -mode_candidates$component_score,
    mode_candidates$candidate
  )
  selected_index <- representation_order[[1L]]
  selected_candidate <- mode_candidates$candidate[[selected_index]]
  selected_representation <- mode_candidates$root_representation[[
    selected_index
  ]]

  list(
    selected_candidate = selected_candidate,
    selected_potency_mode = selected_mode,
    selected_root_representation = selected_representation,
    default_potency_mode = default_mode,
    default_available = default_available,
    unavailable_candidates = unavailable_candidates,
    best_alignment_mode = best_mode,
    unconstrained_best_alignment_mode = unconstrained_best_mode,
    endpoint_alignment_gain = as.numeric(endpoint_alignment_gain),
    endpoint_selection_status = ifelse(
      !default_available,
      "default_unavailable",
      ifelse(endpoint_switch, "strong_alignment_switch", "default_mode_retained")
    ),
    endpoint_topology_gate_status = ifelse(
      endpoint_switch,
      "terminal_switch_allowed",
      ifelse(
        unconstrained_best_mode != default_mode &&
          !state_candidates$terminal_leaf_eligible[[
            unconstrained_best_state_index
          ]],
        "nonterminal_switch_blocked",
        "default_mode_retained"
      )
    ),
    selected_terminal_leaf_eligible = as.logical(
      mode_candidates$terminal_leaf_eligible[[selected_index]]
    ),
    representation_ordering_gain = as.numeric(
      representation_gain[[selected_index]]
    ),
    representation_alignment_gain = as.numeric(
      representation_alignment_gain[[selected_index]]
    ),
    representation_selection_status = ifelse(
      selected_representation == "cellwise_component",
      "ordering_alignment_switch",
      "state_retained"
    ),
    min_alignment_gain = as.numeric(min_alignment_gain),
    min_ordering_gain = as.numeric(min_ordering_gain),
    min_representation_alignment_gain = as.numeric(
      min_representation_alignment_gain
    )
  )
}

#' Flip pseudotime when reversed direction scores higher on embedding PCs (no GT).
#'
#' @param pseudotime numeric vector of per-cell pseudotime values.
#' @param embedding numeric matrix (`n_cells` x `n_dims`).
#' @param n_pc number of principal components averaged in the blind score.
#' @param min_delta minimum blind-score improvement required before flipping.
#' @return numeric vector, same length as `pseudotime`.
#' @export
OrientPseudotimeBlind <- function(pseudotime, embedding, n_pc = 3L, min_delta = 0.05) {
  emb <- as.matrix(embedding)
  pt <- as.numeric(pseudotime)
  vv <- is.finite(pt)
  if (sum(vv) < 50L) return(pt)
  s0 <- BlindEmbeddingScore(pt, emb, n_pc)
  pt_flip <- pt
  pt_flip[vv] <- 1 - pt[vv]
  s1 <- BlindEmbeddingScore(pt_flip, emb, n_pc)
  if (is.finite(s0) && s0 >= 0) return(pt)
  if (is.finite(s1) && s1 > 0 && (!is.finite(s0) || s1 - s0 >= min_delta)) pt_flip else pt
}

#' Orient pseudotime using an explicitly requested count-depth anchor.
#'
#' This manual helper orients pseudotime to be positively correlated with total
#' counts. Library size is a technical covariate rather than a universal
#' maturity marker, so this helper is excluded from automatic root inference.
#' Biological interpretation requires an independent root prior.
#'
#' @param pseudotime numeric vector.
#' @param total_counts numeric vector of per-cell total counts, same length.
#' @param min_cor minimum |Spearman(pt, total_counts)| to act on.
#' @return numeric vector (possibly globally flipped; ranking preserved).
#' @export
OrientPseudotimeByDepth <- function(pseudotime, total_counts, min_cor = 0.05) {
  pt <- as.numeric(pseudotime); tc <- as.numeric(total_counts)
  vv <- is.finite(pt) & is.finite(tc)
  if (sum(vv) < 50L) return(pt)
  rho <- suppressWarnings(stats::cor(pt[vv], tc[vv], method = "spearman"))
  if (is.finite(rho) && rho < -min_cor) { m <- max(pt[vv]); pt[vv] <- m - pt[vv] }
  pt
}

#' Orient pseudotime using a biological maturity (potency) prior.
#'
#' The root (pseudotime minimum) is assumed to sit at the high-potency end.
#' This assumption fails in some reverse-depth datasets, so the function is a
#' manual prior-driven orientation helper and is not used automatically.
#' It changes only the global sign, not the fuzzy-Dijkstra ordering.
#'
#' @param pseudotime numeric vector.
#' @param potency numeric vector (higher = more plastic/early), same length.
#' @param min_cor minimum Spearman(pt, potency) above which the direction is flipped.
#' @return numeric vector (possibly globally flipped; ranking preserved).
#' @export
OrientPseudotimeByPotency <- function(pseudotime, potency, min_cor = 0.10) {
  pt  <- as.numeric(pseudotime)
  pot <- as.numeric(potency)
  vv  <- is.finite(pt) & is.finite(pot)
  if (sum(vv) < 50L) return(pt)
  rho <- suppressWarnings(stats::cor(pt[vv], pot[vv], method = "spearman"))
  if (is.finite(rho) && rho > min_cor) {
    m <- max(pt[vv]); pt[vv] <- m - pt[vv]
  }
  pt
}

#' Audit whether pseudotime direction is supported by blind diagnostics.
#'
#' This diagnostic does not flip pseudotime. It separates two pieces of
#' evidence that are easy to conflate: a potency prior supplied by the caller
#' and the sign-anchored PC convention used by `BlindEmbeddingScore()`. The PC
#' convention is useful for detecting internally conflicting orientations, but
#' it is not a biological ground truth.
#'
#' @param pseudotime numeric vector of per-cell pseudotime values.
#' @param potency optional numeric vector; higher values are interpreted as
#'   more root-like for this audit.
#' @param embedding optional numeric matrix used for the blind PC score.
#' @param n_pc number of PCs passed to `BlindEmbeddingScore()`.
#' @param endpoint_fraction lower/upper pseudotime fraction used as root and
#'   terminal endpoint summaries.
#' @param min_potency_delta minimum scaled endpoint potency difference.
#' @param min_blind_delta minimum blind-score margin for a PC-convention call.
#' @return named list with potency, PC-convention and combined status fields.
#' @export
AuditPseudotimeDirection <- function(pseudotime,
                                     potency = NULL,
                                     embedding = NULL,
                                     n_pc = 3L,
                                     endpoint_fraction = 0.05,
                                     min_potency_delta = 0.10,
                                     min_blind_delta = 0.05) {
  pt <- as.numeric(pseudotime)
  vv <- is.finite(pt)
  endpoint_fraction <- as.numeric(endpoint_fraction)
  min_potency_delta <- as.numeric(min_potency_delta)
  min_blind_delta <- as.numeric(min_blind_delta)
  if (sum(vv) < 50L ||
      length(endpoint_fraction) != 1L || !is.finite(endpoint_fraction) ||
      endpoint_fraction <= 0 || endpoint_fraction >= 0.5) {
    return(list(
      direction_status = "insufficient_cells",
      potency_status = "not_tested",
      pc_convention_status = "not_tested",
      potency_spearman = NA_real_,
      endpoint_potency_delta = NA_real_,
      blind_score_current = NA_real_,
      blind_score_flipped = NA_real_,
      blind_score_margin = NA_real_
    ))
  }

  potency_status <- "not_tested"
  potency_spearman <- NA_real_
  endpoint_potency_delta <- NA_real_
  root_potency_median <- NA_real_
  terminal_potency_median <- NA_real_
  if (!is.null(potency)) {
    pot <- as.numeric(potency)
    pp <- is.finite(pt) & is.finite(pot)
    if (sum(pp) >= 50L) {
      pot_n <- .scale01(pot)
      q_root <- stats::quantile(pt[pp], endpoint_fraction, names = FALSE,
                                type = 8, na.rm = TRUE)
      q_term <- stats::quantile(pt[pp], 1 - endpoint_fraction, names = FALSE,
                                type = 8, na.rm = TRUE)
      root_idx <- pp & pt <= q_root
      term_idx <- pp & pt >= q_term
      root_potency_median <- stats::median(pot_n[root_idx], na.rm = TRUE)
      terminal_potency_median <- stats::median(pot_n[term_idx], na.rm = TRUE)
      endpoint_potency_delta <- root_potency_median - terminal_potency_median
      potency_spearman <- suppressWarnings(stats::cor(
        pt[pp], pot_n[pp], method = "spearman"
      ))
      potency_status <- if (!is.finite(endpoint_potency_delta) ||
                            abs(endpoint_potency_delta) < min_potency_delta) {
        "uninformative"
      } else if (endpoint_potency_delta > 0) {
        "supports_current_root"
      } else {
        "supports_flipped_root"
      }
    } else {
      potency_status <- "insufficient_cells"
    }
  }

  pc_convention_status <- "not_tested"
  blind_score_current <- NA_real_
  blind_score_flipped <- NA_real_
  blind_score_margin <- NA_real_
  if (!is.null(embedding)) {
    emb <- .prepare_embedding(embedding)
    if (nrow(emb) == length(pt)) {
      blind_score_current <- BlindEmbeddingScore(pt, emb, n_pc = n_pc)
      pt_flip <- pt
      pt_flip[vv] <- max(pt[vv], na.rm = TRUE) - pt[vv]
      blind_score_flipped <- BlindEmbeddingScore(pt_flip, emb, n_pc = n_pc)
      blind_score_margin <- blind_score_current - blind_score_flipped
      pc_convention_status <- if (!is.finite(blind_score_margin) ||
                                  abs(blind_score_margin) < min_blind_delta) {
        "ambiguous"
      } else if (blind_score_margin > 0) {
        "supports_current_root"
      } else {
        "supports_flipped_root"
      }
    } else {
      pc_convention_status <- "length_mismatch"
    }
  }

  direction_status <- if (potency_status %in% c("supports_current_root", "supports_flipped_root") &&
                          pc_convention_status %in% c("supports_current_root", "supports_flipped_root") &&
                          !identical(potency_status, pc_convention_status)) {
    "pc_convention_conflict"
  } else if (potency_status %in% c("supports_current_root", "supports_flipped_root")) {
    potency_status
  } else if (pc_convention_status %in% c("supports_current_root", "supports_flipped_root")) {
    paste0("pc_convention_", pc_convention_status)
  } else if (potency_status == "uninformative" || pc_convention_status == "ambiguous") {
    "ambiguous"
  } else {
    "not_tested"
  }

  list(
    direction_status = direction_status,
    potency_status = potency_status,
    pc_convention_status = pc_convention_status,
    potency_spearman = potency_spearman,
    endpoint_potency_delta = endpoint_potency_delta,
    root_potency_median = root_potency_median,
    terminal_potency_median = terminal_potency_median,
    blind_score_current = blind_score_current,
    blind_score_flipped = blind_score_flipped,
    blind_score_margin = blind_score_margin,
    endpoint_fraction = endpoint_fraction,
    min_potency_delta = min_potency_delta,
    min_blind_delta = min_blind_delta
  )
}

#' Pick a root cell without labels or pseudotime ground truth.
#'
#' Candidate roots are scored by how well their Dijkstra ordering follows the
#' embedding manifold. The signed PC score is retained as a diagnostic only:
#' PC signs are arbitrary, so opposite endpoints with similar absolute scores are
#' reported as direction-ambiguous instead of being silently resolved by a
#' coordinate-sign convention.
#'
#' @param embedding numeric matrix (`n_cells` x `n_dims`).
#' @param cluster_labels integer vector of length `n_cells`, 1-based cluster ids.
#' @param nn_idx integer matrix (`n_cells` x `k`), 1-based neighbour ids; 0/NA = pad.
#' @param nn_dists numeric matrix (`n_cells` x `k`), L2 distances.
#' @param rank_normalize logical; rank-normalise pseudotime before blind scoring.
#' @param seed integer seed for tie-breaking.
#' @param min_delta absolute-score margin below which opposite signed candidates
#'   are treated as indistinguishable by the blind embedding convention.
#' @return list with `start_cell`, `method`, signed `scores`, absolute
#'   `ordering_scores`, `blind_score`, `ordering_score`, and direction diagnostics.
#' @export
InferStartCellAuto <- function(embedding,
                                 cluster_labels,
                                 nn_idx,
                                 nn_dists,
                                 rank_normalize = FALSE,
                                 seed = 1L,
                                 min_delta = 0.05) {
  emb <- .prepare_embedding(embedding)
  cl <- .as_integerish(cluster_labels, "cluster_labels")
  if (length(cl) != nrow(emb) || any(cl < 1L)) {
    stop("`cluster_labels` must contain one positive integer per embedding row",
         call. = FALSE)
  }
  knn <- .prepare_knn_inputs(nn_idx, nn_dists)
  if (nrow(knn$indices) != nrow(emb)) {
    stop("kNN rows must equal nrow(embedding)", call. = FALSE)
  }
  if (!is.logical(rank_normalize) || length(rank_normalize) != 1L ||
      is.na(rank_normalize)) {
    stop("`rank_normalize` must be TRUE or FALSE", call. = FALSE)
  }
  seed <- .validate_non_negative_int(seed, "seed")
  min_delta <- as.numeric(min_delta)
  if (length(min_delta) != 1L || !is.finite(min_delta) || min_delta < 0) {
    stop("`min_delta` must be a non-negative numeric scalar", call. = FALSE)
  }
  methods <- c("cluster_pc1_min", "cluster_pc1_max", "mst_extremum")
  cells <- vapply(methods, function(method) {
    InferStartCell(emb, cl, method, seed)
  }, integer(1L))
  best_rank <- isTRUE(rank_normalize)
  scores <- vapply(seq_along(methods), function(i) {
    pt <- ComputePseudotimeFuzzy(knn$indices, knn$distances, cells[[i]])
    if (best_rank) {
      ok <- is.finite(pt)
      if (sum(ok) >= 2L) {
        r <- rank(pt[ok], ties.method = "average")
        pt[ok] <- (r - 1) / (length(r) - 1)
      }
    }
    ff <- mean(is.finite(pt))
    if (ff < 0.99) -Inf else BlindEmbeddingScore(pt, emb)
  }, numeric(1L))
  names(scores) <- methods
  ordering_scores <- abs(scores)
  ordering_scores[!is.finite(scores)] <- -Inf
  names(ordering_scores) <- methods
  lineage_counts <- vapply(cells, function(cell) {
    start_cluster <- as.integer(cl[[cell]])
    length(ExtractLineages(emb, cl, start_cluster = start_cluster))
  }, integer(1L))
  names(lineage_counts) <- methods
  if (!any(is.finite(scores))) {
    stop("InferStartCellAuto failed: no candidate root produced finite pseudotime coverage")
  }
  best_order <- max(ordering_scores[is.finite(ordering_scores)])
  near_best <- which(is.finite(ordering_scores) &
                       (best_order - ordering_scores) <= min_delta)
  near_signs <- sign(scores[near_best])
  near_signs <- near_signs[near_signs != 0]
  direction_ambiguous <- length(near_best) > 1L &&
    length(unique(near_signs)) > 1L
  best_i <- if (isTRUE(direction_ambiguous)) {
    near_best[which.min(cells[near_best])]
  } else {
    which.max(ordering_scores)
  }
  best_method <- methods[[best_i]]
  best_cell <- cells[[best_i]]
  best_score <- scores[[best_i]]
  best_ordering_score <- ordering_scores[[best_i]]
  root_direction_status <- if (isTRUE(direction_ambiguous)) {
    "ambiguous_pc_sign"
  } else {
    "pc_convention_selected"
  }
  list(
    start_cell = best_cell,
    method = best_method,
    rank_normalize = best_rank,
    scores = as.list(scores),
    ordering_scores = as.list(ordering_scores),
    lineage_counts = as.list(lineage_counts),
    blind_score = best_score,
    ordering_score = best_ordering_score,
    root_direction_status = root_direction_status,
    direction_ambiguous = direction_ambiguous,
    direction_tie_margin = min_delta
  )
}

#' Attractor-guided energy root selection without ground-truth pseudotime.
#'
#' By default, terminal attractors are approximated as peripheral, low-potency
#' cluster basins. The experimental retention option replaces peripherality
#' with mean within-state transition retention after excluding self edges.
#' Cell-level evidence combines high potency, high fate entropy,
#' source-like expected potency drop under the local transition distribution,
#' distance from the nearest attractor
#' (a maximin separation criterion), and weakly committed to any one basin.
#' The root state maximizes the median cell score. Its representative cell
#' maximizes within-state directed kNN centrality, with potency and cell index as
#' deterministic tie-breakers. This prevents an isolated high-score cell from
#' changing the inferred root state.
#'
#' Fate absorption is evaluated by a native double-precision sparse harmonic
#' solve. Failure to converge is reported instead of returning the initial prior;
#' state aggregation and deterministic representative selection are performed
#' in R. Suitable for heuristic root selection when no biological root prior is
#' available.
#'
#' @param embedding numeric matrix (`n_cells` x `n_dims`).
#' @param cluster_labels integer vector of length `n_cells`, 1-based cluster ids.
#' @param nn_idx integer matrix (`n_cells` x `k`), 1-based neighbour ids; 0/NA = pad.
#' @param nn_dists numeric matrix (`n_cells` x `k`), L2 distances.
#' @param potency numeric vector of length `n_cells`; higher means more plastic.
#' @param top_fraction fraction defining the high-score enrichment diagnostic.
#' @param n_attractors number of terminal attractor clusters (`0` = auto).
#' @param max_iter maximum GMRES steps for each independent fate probability.
#' @param tol probability tolerance; the backward residual tolerance is tightened
#'   by the square root of machine precision for rare absorbing events.
#' @param weights numeric length-5 vector weighting potency, fate entropy,
#'   nearest-attractor distance, source flow, and (1 - commitment).
#' @param seed integer seed for deterministic tie-breaking.
#' @param on_unreachable `"error"` (default) or `"return"`. The latter returns
#'   `status = "unreachable_terminals"`, missing root scores, and the indices of
#'   all unreachable cells when the candidate has no defined absorption model.
#'   It does not remove cells or supply replacement probabilities. Numerical
#'   nonconvergence still raises an error.
#' @param attractor_selection `"peripheral"` (default) or experimental
#'   `"retention"`. Both subtract mean state potency from their state evidence
#'   and retain the same number of terminal candidates. Retention averages
#'   each cell's probability of staying within its state conditional on a
#'   non-self transition; a cell with no outgoing non-self edge retains mass 1.
#' @param fate_solver `"gmres"` (default) or experimental `"sparse_lu"`.
#'   Sparse LU solves the same absorption equations after analytically removing
#'   self transitions and normalizing the remaining rows. It uses more memory
#'   on large graphs and does not repair unreachable terminal states.
#' @param complete_terminal_coverage Add a deterministic minimum-potency cell
#'   from each uncovered sink strongly connected component. Defaults to FALSE.
#' @return list with start cell, root state, attractor clusters, cell and cluster
#'   scores, graph components, availability status, and unreachable cell indices.
#' @export
InferStartCellAttractor <- function(embedding,
                                   cluster_labels,
                                   nn_idx,
                                   nn_dists,
                                   potency,
                                   top_fraction = 0.05,
                                   n_attractors = 0L,
                                   max_iter = 3000L,
                                   tol = 1e-5,
                                   weights = c(0.25, 0.20, 0.30, 0.15, 0.10),
                                   seed = 1L,
                                   on_unreachable = c("error", "return"),
                                   attractor_selection = c("peripheral", "retention"),
                                   fate_solver = c("gmres", "sparse_lu"),
                                   complete_terminal_coverage = FALSE) {
  on_unreachable <- match.arg(on_unreachable)
  attractor_selection <- match.arg(attractor_selection)
  fate_solver <- match.arg(fate_solver)
  set.seed(seed)
  emb <- as.matrix(embedding)
  cl_map <- .compact_cluster_labels(cluster_labels)
  cl <- cl_map$labels
  n_cells <- nrow(emb)
  n_dims <- ncol(emb)
  k <- ncol(nn_idx)
  n_clusters <- length(cl_map$levels)
  stopifnot(length(cl) == n_cells, length(potency) == n_cells,
            is.matrix(nn_idx), is.matrix(nn_dists),
            nrow(nn_idx) == n_cells, nrow(nn_dists) == n_cells,
            ncol(nn_idx) == ncol(nn_dists))
  .assert_indices_in_range(nn_idx, 1L, n_cells, "nn_idx", allow_zero = TRUE, allow_na = TRUE)
  out <- .infer_start_cell_attractor_impl(
    emb, cl, cl_map$levels, nn_idx, nn_dists, potency,
    top_fraction = top_fraction,
    n_attractors = n_attractors,
    max_iter = max_iter,
    tol = tol,
    weights = weights,
    on_unreachable = on_unreachable,
    attractor_selection = attractor_selection,
    fate_solver = fate_solver,
    complete_terminal_coverage = complete_terminal_coverage
  )
  out$attractor_clusters <- as.integer(cl_map$levels[out$attractor_clusters])
  out
}

#' @keywords internal
.pc1_scores_signed <- function(emb) {
  pcs <- stats::prcomp(emb, center = TRUE, rank. = min(5L, ncol(emb), nrow(emb) - 1L))
  pc1 <- pcs$x[, 1L]
  loading <- pcs$rotation[, 1L]
  anchor <- which.max(abs(loading))
  sgn <- sign(loading[[anchor]])
  if (sgn == 0) sgn <- 1
  pc1 * sgn
}

#' @keywords internal
.infer_start_cluster_pc1_min <- function(emb, cl, seed) {
  set.seed(seed)
  pc1 <- .pc1_scores_signed(emb)
  cl_means <- tapply(pc1, cl, mean)
  start_cluster <- as.integer(names(cl_means)[which.min(cl_means)])
  .cell_extreme_pc1_in_cluster(cl, start_cluster, pc1, which.min)
}

#' @keywords internal
.infer_start_cluster_pc1_max <- function(emb, cl, seed) {
  set.seed(seed)
  pc1 <- .pc1_scores_signed(emb)
  cl_means <- tapply(pc1, cl, mean)
  start_cluster <- as.integer(names(cl_means)[which.max(cl_means)])
  .cell_extreme_pc1_in_cluster(cl, start_cluster, pc1, which.max)
}

#' @keywords internal
.infer_start_mst_extremum <- function(emb, cl, seed) {
  set.seed(seed)
  pc1 <- .pc1_scores_signed(emb)
  cl_means <- tapply(pc1, cl, mean)
  start_cluster <- as.integer(names(cl_means)[which.min(cl_means)])

  lineages <- ExtractLineages(emb, cl, start_cluster = start_cluster)
  if (length(lineages) == 0L) {
    return(.infer_start_cluster_pc1_min(emb, cl, seed))
  }

  leaf_clusters <- unique(vapply(
    lineages,
    function(path) as.integer(path[length(path)]),
    integer(1L)
  ))

  global <- colMeans(emb)
  cluster_ids <- sort(unique(cl))
  centroid_mat <- do.call(
    rbind,
    lapply(cluster_ids, function(cid) colMeans(emb[cl == cid, , drop = FALSE]))
  )
  rownames(centroid_mat) <- as.character(cluster_ids)

  dists <- vapply(
    leaf_clusters,
    function(lc) {
      cvec <- centroid_mat[as.character(lc), , drop = TRUE]
      sqrt(sum((cvec - global)^2))
    },
    numeric(1)
  )
  root_cluster <- leaf_clusters[which.max(dists)]
  .cell_farthest_from_global_in_cluster(emb, cl, root_cluster, global, pc1)
}

#' @keywords internal
.cell_extreme_pc1_in_cluster <- function(cl, cluster_id, pc1_scores, picker) {
  idx <- which(cl == cluster_id)
  if (length(idx) == 0L) {
    return(picker(pc1_scores))
  }
  idx[picker(pc1_scores[idx])]
}

#' @keywords internal
.cell_farthest_from_global_in_cluster <- function(emb, cl, cluster_id, global, pc1_scores) {
  idx <- which(cl == cluster_id)
  if (length(idx) == 0L) {
    return(which.max(abs(pc1_scores)))
  }
  sub <- emb[idx, , drop = FALSE]
  d2 <- rowSums((sub - matrix(global, nrow = length(idx), ncol = ncol(emb), byrow = TRUE))^2)
  idx[which.max(d2)]
}

#' @keywords internal
.scale01 <- function(x) {
  r <- range(x, finite = TRUE)
  if (!is.finite(r[1L]) || diff(r) <= 0) return(rep(0.5, length(x)))
  (x - r[1L]) / diff(r)
}

#' @keywords internal
.select_cluster_aggregate_root <- function(root_score, cl, nn_idx, potency,
                                           top_fraction = 0.05) {
  cluster_scores <- tapply(root_score, cl, stats::median)
  n_top <- max(2L, ceiling(length(root_score) * top_fraction))
  high_score_cutoff <- sort(root_score, decreasing = TRUE)[[
    min(n_top, length(root_score))
  ]]
  high_score_fraction <- tapply(
    root_score >= high_score_cutoff,
    cl,
    mean
  )
  root_cluster <- as.integer(names(cluster_scores)[which.max(cluster_scores)])
  root_region <- which(cl == root_cluster)

  n_cells <- length(cl)
  valid_neighbor <- nn_idx > 0L & nn_idx <= n_cells
  safe_neighbor <- nn_idx
  safe_neighbor[!valid_neighbor] <- 1L
  root_source <- matrix(
    cl == root_cluster,
    nrow = n_cells,
    ncol = ncol(nn_idx)
  )
  root_target <- matrix(
    cl[as.vector(safe_neighbor)] == root_cluster,
    nrow = n_cells,
    ncol = ncol(nn_idx)
  )
  root_edge <- valid_neighbor & root_source & root_target
  root_centrality <- rowSums(root_edge) + tabulate(
    safe_neighbor[root_edge],
    nbins = n_cells
  )
  start_cell <- root_region[[order(
    -root_centrality[root_region],
    -potency[root_region],
    root_region
  )[1L]]]

  list(
    start_cell = as.integer(start_cell),
    root_region = as.integer(root_region),
    root_cluster = as.integer(root_cluster),
    cluster_scores = as.numeric(cluster_scores),
    cluster_ids = as.integer(names(cluster_scores)),
    high_score_fraction = as.numeric(high_score_fraction),
    root_centrality = as.numeric(root_centrality)
  )
}

#' @keywords internal
.select_cellwise_component_root <- function(root_score, nn_idx,
                                            top_fraction = 0.05) {
  n_top <- max(2L, ceiling(length(root_score) * top_fraction))
  top_cells <- order(root_score, decreasing = TRUE)[seq_len(n_top)]
  components <- .knn_components_induced(
    nn_idx,
    top_cells,
    length(root_score)
  )
  component_medians <- vapply(
    components,
    function(component) stats::median(root_score[component]),
    numeric(1L)
  )
  root_region <- components[[which.max(component_medians)]]
  start_cell <- root_region[[which.max(root_score[root_region])]]
  list(
    start_cell = as.integer(start_cell),
    root_region = as.integer(root_region),
    component_median_score = as.numeric(max(component_medians))
  )
}

#' @keywords internal
.infer_start_cell_attractor_impl <- function(emb, cl, cl_levels, nn_idx, nn_dists, potency,
                                             top_fraction, n_attractors, max_iter, tol, weights,
                                             on_unreachable = "error",
                                             attractor_selection = "peripheral",
                                             fate_solver = "gmres",
                                             complete_terminal_coverage = FALSE) {
  n_cells <- nrow(emb)
  pot_n <- .scale01(as.numeric(potency))
  global <- colMeans(emb)
  dist_global <- sqrt(rowSums((emb - matrix(global, n_cells, ncol(emb), byrow = TRUE))^2))
  peripheral <- .scale01(dist_global)

  cl_ids <- sort(unique(cl))
  cl_peri <- vapply(cl_ids, function(cid) mean(peripheral[cl == cid]), numeric(1))
  cl_pot <- vapply(cl_ids, function(cid) mean(pot_n[cl == cid]), numeric(1))
  cl_retention <- rep(NA_real_, length(cl_ids))
  cl_evidence <- cl_peri
  if (identical(attractor_selection, "retention")) {
    nonself <- nn_idx > 0L & nn_idx <= n_cells & is.finite(nn_dists) &
      nn_idx != row(nn_idx)
    neighbor <- nn_idx
    neighbor[!nonself] <- 1L
    transition_weight <- 1 / (nn_dists + 1e-6)
    transition_weight[!nonself] <- 0
    outgoing_weight <- rowSums(transition_weight)
    within_state <- matrix(cl[as.vector(neighbor)], nrow = n_cells) == cl
    retention <- rowSums(transition_weight * within_state) / outgoing_weight
    retention[outgoing_weight == 0] <- 1
    cl_retention <- as.numeric(rowsum(retention, cl, reorder = TRUE)) / tabulate(cl)
    cl_evidence <- cl_retention
  }
  attractor_score <- cl_evidence - cl_pot
  n_attr <- if (as.integer(n_attractors) > 0L) {
    min(as.integer(n_attractors), length(cl_ids))
  } else {
    max(2L, ceiling(length(cl_ids) * 0.2))
  }
  attractor_clusters <- cl_ids[order(attractor_score, decreasing = TRUE)[seq_len(n_attr)]]
  attractor_states <- data.frame(
    cluster = as.integer(cl_levels[cl_ids]), n_cells = tabulate(cl),
    potency = cl_pot, peripheral = cl_peri, retention = cl_retention,
    selection_score = attractor_score, selected = cl_ids %in% attractor_clusters
  )

  terminal_cells <- vapply(attractor_clusters, function(cid) {
    idx <- which(cl == cid)
    idx[which.min(pot_n[idx])]
  }, integer(1))

  # Every sink strongly connected component needs an absorbing representative.
  # Weak graph connectivity alone does not make the selected terminals reachable.
  terminal_coverage_added <- integer()
  uncovered_cells <- rs_attractor_unreachable_cells(
    as.integer(t(nn_idx)), as.numeric(t(nn_dists)),
    as.integer(terminal_cells), as.integer(n_cells), as.integer(ncol(nn_idx)))
  switch(as.character(isTRUE(complete_terminal_coverage) &&
                      length(uncovered_cells) > 0L), "TRUE" = {
    valid_edge <- nn_idx > 0L & nn_idx <= n_cells & is.finite(nn_dists)
    edge_from <- row(nn_idx)[valid_edge]
    edge_to <- nn_idx[valid_edge]
    adjacency <- Matrix::sparseMatrix(i = edge_from, j = edge_to, x = 1,
      dims = c(n_cells, n_cells))
    directed_graph <- igraph::graph_from_adjacency_matrix(adjacency,
      mode = "directed", weighted = TRUE)
    components <- igraph::components(directed_graph, mode = "strong")
    component <- components$membership
    outgoing <- unique(component[edge_from][
      component[edge_from] != component[edge_to]])
    uncovered_sinks <- setdiff(setdiff(seq_len(components$no), outgoing),
      component[terminal_cells])
    candidate_order <- order(component, pot_n, seq_len(n_cells))
    representatives <- candidate_order[!duplicated(component[candidate_order])]
    terminal_coverage_added <- as.integer(representatives[
      component[representatives] %in% uncovered_sinks])
    terminal_cells <- c(terminal_cells, terminal_coverage_added)
    attractor_clusters <- unique(c(attractor_clusters, cl[terminal_coverage_added]))
    attractor_states$selected <- cl_ids %in% attractor_clusters
  })

  if (identical(on_unreachable, "return")) {
    unreachable <- rs_attractor_unreachable_cells(
      as.integer(t(nn_idx)), as.numeric(t(nn_dists)),
      as.integer(terminal_cells), as.integer(n_cells), as.integer(ncol(nn_idx))
    )
    if (length(unreachable)) {
      return(list(
        status = "unreachable_terminals",
        unreachable_cells = unreachable,
        start_cell = NA_integer_,
        cellwise_start_cell = NA_integer_,
        cellwise_component_score = NA_real_,
        root_region = integer(),
        root_cluster = NA_integer_,
        attractor_clusters = as.integer(attractor_clusters),
        attractor_selection = attractor_selection,
        fate_solver = fate_solver,
        fate_residual_max = NA_real_,
        attractor_states = attractor_states,
        terminal_cells = as.integer(terminal_cells),
        scores = NULL,
        cluster_scores = NULL,
        components = data.frame(median_score = NA_real_, selected = FALSE)
      ))
    }
  }

  fate_residual_max <- NA_real_
  gmres_result <- NULL
  if (fate_solver == "gmres") {
    gmres_result <- try(rs_compute_attractor_fate_f64(
      as.integer(t(nn_idx)),
      as.numeric(t(nn_dists)),
      as.integer(terminal_cells),
      as.integer(n_cells),
      as.integer(ncol(nn_idx)),
      as.integer(max_iter),
      as.numeric(tol)
    ), silent = TRUE)
    if (inherits(gmres_result, "try-error")) {
      if (!grepl("Absorption solve|Absorption preconditioner",
                 as.character(gmres_result))) stop(gmres_result)
      fate_solver <- "sparse_lu"
    }
  }
  fate_mat <- switch(fate_solver, gmres = matrix(
    gmres_result,
    nrow = n_cells,
    ncol = length(terminal_cells),
    byrow = TRUE
  ), sparse_lu = {
    unreachable <- rs_attractor_unreachable_cells(
      as.integer(t(nn_idx)), as.numeric(t(nn_dists)),
      as.integer(terminal_cells), as.integer(n_cells), as.integer(ncol(nn_idx)))
    if (length(unreachable)) {
      stop("Absorption is undefined: cells cannot reach the selected terminal states", call. = FALSE)
    }
    nonself <- is.finite(nn_idx) & nn_idx > 0L & nn_idx <= n_cells &
      is.finite(nn_dists) & nn_idx != row(nn_idx)
    transition <- Matrix::sparseMatrix(i = row(nn_idx)[nonself], j = nn_idx[nonself],
      x = 1 / (nn_dists[nonself] + 1e-6), dims = c(n_cells, n_cells))
    transition <- Matrix::Diagonal(x = 1 / pmax(Matrix::rowSums(transition),
      .Machine$double.xmin)) %*% transition
    transient <- setdiff(seq_len(n_cells), terminal_cells)
    independent <- seq_len(length(terminal_cells) - 1L)
    probability <- matrix(0, n_cells, length(terminal_cells))
    probability[cbind(terminal_cells, seq_along(terminal_cells))] <- 1
    if (length(transient) > 0L && length(independent) > 0L) {
      system <- Matrix::Diagonal(length(transient)) - transition[transient, transient, drop = FALSE]
      rhs <- as.matrix(transition[transient, terminal_cells[independent], drop = FALSE])
      probability[transient, independent] <- rs_sparse_absorption_lu64(system, rhs)
    }
    probability[transient, length(terminal_cells)] <-
      1 - rowSums(probability[transient, independent, drop = FALSE])
    residual <- probability[transient, , drop = FALSE] -
      transition[transient, , drop = FALSE] %*% probability
    fate_residual_max <- max(c(0, abs(as.matrix(residual))))
    residual_limit <- max(8 * .Machine$double.eps, tol * sqrt(.Machine$double.eps))
    if (!all(is.finite(probability)) || min(probability) < -tol || max(probability) > 1 + tol ||
        !is.finite(fate_residual_max) || fate_residual_max > residual_limit) {
      stop("Sparse LU absorption solve did not satisfy probability or residual tolerance", call. = FALSE)
    }
    probability
  })
  fate_entropy <- apply(fate_mat, 1L, function(p) {
    p <- p[p > 1e-12]
    if (!length(p)) return(0)
    -sum(p * log(p))
  })
  fate_entropy <- .scale01(fate_entropy)
  commitment <- apply(fate_mat, 1L, max)

  att_centroids <- do.call(
    rbind,
    lapply(attractor_clusters, function(cid) colMeans(emb[cl == cid, , drop = FALSE]))
  )
  attractor_dist_sq <- outer(
    rowSums(att_centroids^2), rowSums(emb^2), "+"
  ) - 2 * tcrossprod(att_centroids, emb)
  attractor_dist <- sqrt(pmax(attractor_dist_sq, 0))
  if (is.null(dim(attractor_dist))) {
    attractor_dist <- matrix(attractor_dist, nrow = 1L)
  }
  nearest_att_dist_n <- .scale01(apply(attractor_dist, 2L, min))
  mean_att_dist_n <- .scale01(colMeans(attractor_dist))

  valid_neighbor <- nn_idx > 0L & nn_idx <= n_cells & is.finite(nn_dists)
  safe_neighbor <- nn_idx
  safe_neighbor[!valid_neighbor] <- 1L
  source_weight <- 1 / (nn_dists + 1e-6)
  source_weight[!valid_neighbor] <- 0
  source_weight_total <- rowSums(source_weight)
  source_delta <- matrix(pot_n, nrow = n_cells, ncol = ncol(nn_idx)) -
    matrix(
      pot_n[as.vector(safe_neighbor)],
      nrow = n_cells,
      ncol = ncol(nn_idx)
    )
  source_flow <- rowSums(source_weight * source_delta) /
    pmax(source_weight_total, .Machine$double.eps)
  source_flow <- .scale01(source_flow)

  w <- as.numeric(weights)
  if (length(w) != 5L || !all(is.finite(w)) || sum(w) <= 0) {
    stop("`weights` must be a finite numeric vector of length 5 with positive sum",
         call. = FALSE)
  }
  w <- w / sum(w)
  root_score <- w[1L] * pot_n + w[2L] * fate_entropy +
    w[3L] * nearest_att_dist_n + w[4L] * source_flow +
    w[5L] * (1 - commitment)

  root_selection <- .select_cluster_aggregate_root(
    root_score,
    cl,
    nn_idx,
    pot_n,
    top_fraction
  )
  root_region <- root_selection$root_region
  start_cell <- root_selection$start_cell
  cellwise_selection <- .select_cellwise_component_root(
    root_score,
    nn_idx,
    top_fraction
  )

  all_comps <- split(seq_len(n_cells), KnnComponentLabels(nn_idx))
  selected_comp <- which.max(vapply(all_comps, function(comp) {
    if (any(comp %in% root_region)) stats::median(root_score[comp]) else -Inf
  }, numeric(1)))

  list(
    status = "available",
    unreachable_cells = integer(),
    start_cell = as.integer(start_cell),
    cellwise_start_cell = cellwise_selection$start_cell,
    cellwise_root_region = cellwise_selection$root_region,
    cellwise_component_score = cellwise_selection$component_median_score,
    method = "attractor_energy",
    root_region = as.integer(root_region),
    root_cluster = as.integer(root_selection$root_cluster),
    attractor_clusters = as.integer(attractor_clusters),
    attractor_selection = attractor_selection,
    fate_solver = fate_solver,
    fate_residual_max = fate_residual_max,
    attractor_states = attractor_states,
    terminal_cells = as.integer(terminal_cells),
    terminal_coverage_added = terminal_coverage_added,
    scores = data.frame(
      potency = pot_n,
      fate_entropy = fate_entropy,
      source_flow = source_flow,
      nearest_attractor_distance = nearest_att_dist_n,
      mean_attractor_distance = mean_att_dist_n,
      commitment = commitment,
      root_score = root_score
    ),
    cluster_scores = data.frame(
      cluster = root_selection$cluster_ids,
      median_score = root_selection$cluster_scores,
      high_score_fraction = root_selection$high_score_fraction,
      selected = root_selection$cluster_ids == root_selection$root_cluster
    ),
    components = data.frame(
      component = seq_along(all_comps),
      size = vapply(all_comps, length, integer(1)),
      median_score = vapply(all_comps, function(comp) stats::median(root_score[comp]), numeric(1)),
      selected = seq_along(all_comps) == selected_comp
    )
  )
}

#' @keywords internal
.knn_row_stochastic <- function(nn_idx, nn_dists, n_cells) {
  k <- ncol(nn_idx)
  lapply(seq_len(n_cells), function(i) {
    nbrs <- nn_idx[i, ]
    d <- nn_dists[i, ]
    ok <- nbrs > 0L & nbrs <= n_cells & is.finite(d)
    if (!any(ok)) return(list(idx = integer(0), w = numeric(0)))
    nbrs <- as.integer(nbrs[ok])
    w <- 1 / (d[ok] + 1e-6)
    w <- w / sum(w)
    list(idx = nbrs, w = w)
  })
}

#' @keywords internal
.knn_sparse_transition <- function(nn_idx, nn_dists, n_cells) {
  valid <- nn_idx > 0L & nn_idx <= n_cells & is.finite(nn_dists)
  source <- row(nn_idx)[valid]
  target <- as.integer(nn_idx[valid])
  inverse_distance <- 1 / (as.numeric(nn_dists[valid]) + 1e-6)
  compact_total <- rowsum(inverse_distance, source, reorder = FALSE)
  row_total <- numeric(n_cells)
  row_total[as.integer(rownames(compact_total))] <- compact_total[, 1L]
  isolated <- which(row_total <= 0)
  Matrix::sparseMatrix(
    i = c(source, isolated),
    j = c(target, isolated),
    x = c(inverse_distance / row_total[source], rep(1, length(isolated))),
    dims = c(n_cells, n_cells),
    giveCsparse = TRUE
  )
}

#' @keywords internal
.markov_fate_probs <- function(n_cells, trans, terminal_cells, max_iter, tol) {
  n_term <- length(terminal_cells)
  if (n_term == 0L) {
    return(matrix(0, nrow = n_cells, ncol = 0L))
  }
  prob <- matrix(1 / n_term, nrow = n_cells, ncol = n_term)
  for (k in seq_len(n_term)) {
    prob[terminal_cells[[k]], ] <- 0
    prob[terminal_cells[[k]], k] <- 1
  }
  for (iter in seq_len(max_iter)) {
    new_prob <- as.matrix(trans %*% prob)
    new_prob[terminal_cells, ] <- 0
    new_prob[cbind(terminal_cells, seq_len(n_term))] <- 1
    max_delta <- max(abs(new_prob - prob))
    prob <- new_prob
    if (max_delta < tol) break
  }
  prob
}

#' @keywords internal
.knn_components_induced <- function(nn_idx, subset, n_cells) {
  in_set <- rep(FALSE, n_cells)
  in_set[subset] <- TRUE
  seen <- rep(FALSE, n_cells)
  comps <- list()
  for (s in subset) {
    if (seen[s]) next
    q <- s
    seen[s] <- TRUE
    comp <- s
    while (length(q)) {
      u <- q[1L]
      q <- q[-1L]
      nbrs <- nn_idx[u, ]
      valid_nbrs <- !is.na(nbrs) & nbrs > 0L & nbrs <= n_cells
      nbrs <- nbrs[valid_nbrs]
      nbrs <- nbrs[in_set[nbrs] & !seen[nbrs]]
      if (length(nbrs)) {
        seen[nbrs] <- TRUE
        comp <- c(comp, nbrs)
        q <- c(q, nbrs)
      }
    }
    comps[[length(comps) + 1L]] <- comp
  }
  comps
}

#' @keywords internal
.cell_closest_to_cluster_centroid <- function(emb, cl, cluster_id, pc1_scores) {
  idx <- which(cl == cluster_id)
  if (length(idx) == 0L) {
    return(which.min(pc1_scores))
  }
  sub <- emb[idx, , drop = FALSE]
  centroid <- colMeans(sub)
  d2 <- rowSums((sub - matrix(centroid, nrow = length(idx), ncol = ncol(emb), byrow = TRUE))^2)
  idx[which.min(d2)]
}
