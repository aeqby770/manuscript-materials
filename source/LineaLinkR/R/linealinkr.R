# LineaLinkR public R API.
#
# Each function validates R inputs, converts 1-based indices and column-major
# matrices to the Rust boundary contract, and calls the corresponding `rs_*`
# shim in R/wrap__LineaLinkR.R.
#
# ── Index / matrix conventions at this boundary ──────────────────────────
# - R supplies 1-based indices (cells, clusters, genes, terminal cells);
#   the −1 translation happens inside src/extendr_bindings.rs.
# - R matrices are column-major; Rust expects row-major flat slices.
#   `as.numeric(t(m))` materialises row-major inputs, and
#   `matrix(flat, ..., byrow = TRUE)` reshapes Rust outputs.

.as_integerish <- function(x, name, allow_na = FALSE) {
  if (!is.numeric(x)) {
    stop(sprintf("`%s` must be numeric integer-valued input", name), call. = FALSE)
  }
  raw <- as.numeric(x)
  finite_or_na <- is.finite(raw) | (allow_na & is.na(raw))
  integer_or_na <- is.na(raw) | raw == round(raw)
  in_integer_range <- is.na(raw) |
    (raw >= -.Machine$integer.max - 1 & raw <= .Machine$integer.max)
  if (!all(finite_or_na & integer_or_na & in_integer_range)) {
    stop(sprintf("`%s` must contain only integer-valued entries", name), call. = FALSE)
  }
  out <- as.integer(raw)
  dim(out) <- dim(x)
  dimnames(out) <- dimnames(x)
  out
}

.compact_cluster_labels <- function(cluster_labels) {
  raw <- .as_integerish(cluster_labels, "cluster_labels")
  if (any(raw < 1L)) {
    stop("`cluster_labels` must be positive 1-based integers")
  }
  levels <- sort(unique(raw))
  list(labels = as.integer(match(raw, levels)), levels = as.integer(levels))
}

.compact_cluster_id <- function(cluster_id, levels) {
  cluster_id <- .as_integerish(cluster_id, "cluster_id")
  if (length(cluster_id) != 1L) {
    stop("`cluster_id` must be a scalar", call. = FALSE)
  }
  as.integer(match(cluster_id, levels))
}

.assert_indices_in_range <- function(x, lower, upper, name,
                                     allow_zero = FALSE, allow_na = FALSE) {
  ix <- .as_integerish(x, name, allow_na = allow_na)
  valid <- rep(FALSE, length(ix))
  if (allow_na) valid <- valid | is.na(ix)
  if (allow_zero) valid <- valid | (!is.na(ix) & ix == 0L)
  valid <- valid | (!is.na(ix) & ix >= lower & ix <= upper)
  if (!all(valid)) {
    stop(sprintf("`%s` contains indices outside [%d, %d]", name, lower, upper))
  }
  invisible(ix)
}

.assert_lr_db_indices <- function(lr_db, n_genes) {
  for (field in c("ligand_rows", "receptor_rows")) {
    for (i in seq_along(lr_db[[field]])) {
      .assert_indices_in_range(lr_db[[field]][[i]], 1L, n_genes,
                               sprintf("lr_db$%s[[%d]]", field, i))
    }
  }
  invisible(lr_db)
}

.assert_celltypes <- function(celltypes, n_cells) {
  ct <- .as_integerish(celltypes, "celltypes")
  if (length(ct) != n_cells || any(ct < 1L)) {
    stop("`celltypes` must be a length-n_cells vector of positive 1-based integers")
  }
  levels <- sort(unique(ct))
  compact <- as.integer(match(ct, levels))
  attr(compact, "levels") <- levels
  compact
}

.prepare_embedding <- function(embedding) {
  if (!is.numeric(embedding)) {
    stop("`embedding` must be numeric", call. = FALSE)
  }
  emb <- as.matrix(embedding)
  if (nrow(emb) < 1L || ncol(emb) < 1L || !all(is.finite(emb))) {
    stop("`embedding` must be a non-empty finite numeric matrix", call. = FALSE)
  }
  emb
}

.prepare_knn_inputs <- function(nn_idx, nn_dists) {
  if (!is.matrix(nn_idx) || !is.matrix(nn_dists) ||
      nrow(nn_idx) != nrow(nn_dists) || ncol(nn_idx) != ncol(nn_dists)) {
    stop("`nn_idx` and `nn_dists` must be matrices with identical dimensions",
         call. = FALSE)
  }
  if (!is.numeric(nn_idx) || !is.numeric(nn_dists)) {
    stop("`nn_idx` and `nn_dists` must be numeric", call. = FALSE)
  }
  n_cells <- nrow(nn_idx)
  if (n_cells < 2L || ncol(nn_idx) < 1L) {
    stop("kNN input must contain at least two cells and one neighbour column",
         call. = FALSE)
  }
  idx <- .as_integerish(nn_idx, "nn_idx", allow_na = TRUE)
  valid_index <- is.na(idx) | idx == 0L | (idx >= 1L & idx <= n_cells)
  if (!all(valid_index)) {
    stop(sprintf("`nn_idx` contains indices outside [1, %d]", n_cells),
         call. = FALSE)
  }
  dists <- matrix(as.numeric(nn_dists), nrow = n_cells, ncol = ncol(nn_idx))
  supplied <- !is.na(idx) & idx != 0L
  if (any(supplied & (!is.finite(dists) | dists < 0))) {
    stop("distances paired with valid neighbours must be finite and non-negative",
         call. = FALSE)
  }
  usable <- supplied & idx != row(idx)
  if (any(rowSums(usable) < 1L)) {
    stop("each cell must have at least one valid non-self neighbour", call. = FALSE)
  }
  idx[!usable] <- 0L
  dists[!usable] <- 0
  list(indices = idx, distances = dists)
}

.prepare_lr_expr <- function(expr) {
  if (!inherits(expr, "dgCMatrix")) {
    stop("`expr` must inherit from `dgCMatrix`", call. = FALSE)
  }
  if (!all(is.finite(expr@x))) {
    stop("`expr` stored values must be finite", call. = FALSE)
  }
  if (any(expr@x < 0)) {
    stop("`expr` must contain non-negative expression values", call. = FALSE)
  }
  Matrix::drop0(expr)
}

.assert_pseudotime_1d <- function(pseudotime, n_cells) {
  if (is.matrix(pseudotime) && ncol(pseudotime) > 1L) {
    stop("`pseudotime` must be a one-dimensional vector; ",
         "per-lineage pseudotime matrices are not supported",
         call. = FALSE)
  }
  pt <- as.numeric(pseudotime)
  if (length(pt) != n_cells) {
    stop("`pseudotime` length must equal ncol(expr) (", n_cells, ")",
         call. = FALSE)
  }
  if (any(is.infinite(pt))) {
    stop("`pseudotime` cannot contain Inf or -Inf", call. = FALSE)
  }
  bad_range <- is.finite(pt) & (pt < 0 | pt > 1)
  if (any(bad_range)) {
    stop("finite pseudotime values must be in [0, 1]", call. = FALSE)
  }
  pt
}

.validate_positive_int <- function(x, name, min = 1L) {
  x <- .as_integerish(x, name)
  if (length(x) != 1L || x < min) {
    stop(sprintf("`%s` must be a finite integer >= %d", name, min), call. = FALSE)
  }
  x
}

.validate_non_negative_int <- function(x, name) {
  x <- .as_integerish(x, name)
  if (length(x) != 1L || x < 0L) {
    stop(sprintf("`%s` must be a finite integer >= 0", name), call. = FALSE)
  }
  x
}

.validate_unit_interval <- function(x, name) {
  x <- as.numeric(x)
  if (length(x) != 1L || anyNA(x) || x < 0 || x > 1) {
    stop(sprintf("`%s` must be a finite numeric in [0, 1]", name), call. = FALSE)
  }
  x
}

.validate_positive_num <- function(x, name) {
  x <- as.numeric(x)
  if (length(x) != 1L || !is.finite(x) || x <= 0) {
    stop(sprintf("`%s` must be a finite numeric > 0", name), call. = FALSE)
  }
  x
}

.validate_seed <- function(seed) {
  seed <- as.numeric(seed)
  if (length(seed) != 1L || !is.finite(seed) || seed < 0 ||
      seed != floor(seed) || seed > 2^53 - 1) {
    stop("`seed` must be an exactly representable non-negative integer",
         call. = FALSE)
  }
  seed
}

.normalize_trajectory_stat <- function(stat) {
  stat <- match.arg(
    stat,
    c("amplitude", "max", "mean_score", "auc", "topk_mean")
  )
  if (identical(stat, "auc")) {
    warning(
      "`stat = 'auc'` is deprecated; use `stat = 'mean_score'` ",
      "(equal-weight bin average, not a pseudotime integral).",
      call. = FALSE
    )
    stat <- "mean_score"
  }
  list(
    stat = stat,
    code = match(stat, c("amplitude", "max", "mean_score", "topk_mean")) - 1L
  )
}

.assert_lineage_weights <- function(lineage_weights, n_cells) {
  lw <- as.matrix(lineage_weights)
  if (nrow(lw) != n_cells) {
    stop("`lineage_weights` nrow must equal length of `pseudotime` / `n_cells` (",
         n_cells, ")", call. = FALSE)
  }
  if (ncol(lw) < 1L) {
    stop("`lineage_weights` must have at least one column", call. = FALSE)
  }
  if (!all(is.finite(lw))) {
    stop("`lineage_weights` must contain only finite values", call. = FALSE)
  }
  if (any(lw < 0)) {
    stop("`lineage_weights` must be non-negative", call. = FALSE)
  }
  rs <- rowSums(lw)
  if (any(rs <= 0)) {
    stop("each cell must have positive total lineage weight", call. = FALSE)
  }
  lw / rs
}

.prepare_block_labels <- function(block_labels, n_cells) {
  if (is.null(block_labels)) {
    out <- rep.int(1L, n_cells)
    attr(out, "levels") <- "all_cells"
    return(out)
  }
  if (!is.atomic(block_labels) || length(block_labels) != n_cells || anyNA(block_labels)) {
    stop("`block_labels` must be NULL or a complete length-n_cells atomic vector",
         call. = FALSE)
  }
  levels <- unique(as.character(block_labels))
  out <- as.integer(match(as.character(block_labels), levels))
  attr(out, "levels") <- levels
  out
}

.calibration_evaluable <- function(p, q) {
  is.finite(p) | is.finite(q)
}

.trajectory_pmax_finite <- function(a, b) {
  ifelse(
    is.finite(a) & is.finite(b),
    pmax(a, b),
    ifelse(is.finite(a), a, ifelse(is.finite(b), b, NA_real_))
  )
}

.intersection_union_p <- function(p_a, p_b, both_evaluable) {
  out <- rep(NA_real_, length(p_a))
  ok <- both_evaluable %in% TRUE
  if (any(ok)) {
    out[ok] <- pmax(p_a[ok], p_b[ok])
  }
  out
}

.filter_formal_evaluable_rows <- function(df) {
  if (!"evaluable" %in% colnames(df)) {
    return(df)
  }
  df[df$evaluable %in% TRUE, , drop = FALSE]
}

.recompute_all_bh_q <- function(df) {
  p_cols <- grep("^p_", colnames(df), value = TRUE)
  for (p_col in p_cols) {
    q_col <- sub("^p_", "q_", p_col)
    df[[q_col]] <- stats::p.adjust(df[[p_col]], method = "BH")
  }
  df
}

.calibration_supported <- function(evaluable, p, q, alpha) {
  ifelse(
    evaluable,
    ifelse(is.finite(q), q <= alpha,
           ifelse(is.finite(p), p <= alpha, NA)),
    NA
  )
}

.layer_prioritized <- function(structurally_eligible, evaluable, supported) {
  ifelse(evaluable, structurally_eligible & supported, NA)
}

.fake_identity_detail_cols <- function() {
  c(
    "p_fake_generic", "p_fake_edge", "p_fake_complex",
    "p_replace_ligand", "p_replace_receptor",
    "p_edge_iu", "p_fake_network", "p_network_rewire_ligand", "p_network_rewire_receptor",
    "p_network_iu", "fake_null_mode", "complex_null_mode", "complex_pair",
    "mean_ligand_match_distance", "mean_receptor_match_distance",
    "n_unique_fake_draws"
  )
}

.normalize_fake_null_mode <- function(mode) {
  mode <- match.arg(
    mode,
    c("context_generic", "edge_rewire", "both", "lr_network_rewire", "all", "triple")
  )
  if (identical(mode, "triple")) "all" else mode
}

.fake_null_mode_code <- function(mode) {
  match(mode, c("context_generic", "edge_rewire", "both", "lr_network_rewire", "all")) - 1L
}

.primary_p_fake <- function(mode, p_network_iu, p_edge_iu, p_fake_edge, p_fake_generic,
                            network_ok, edge_ok, generic_ok) {
  use_network <- mode %in% c("all", "lr_network_rewire") & is.finite(p_network_iu)
  ifelse(
    use_network,
    p_network_iu,
    ifelse(
      is.finite(p_edge_iu),
      p_edge_iu,
      ifelse(edge_ok, p_fake_edge, ifelse(generic_ok, p_fake_generic, NA_real_))
    )
  )
}

.demote_suffixed_cols <- function(df, cols, layer) {
  for (col in cols) {
    col2 <- paste0(col, ".", layer)
    if (col2 %in% colnames(df) && !col %in% colnames(df)) {
      df[[col]] <- df[[col2]]
    }
    if (col2 %in% colnames(df)) {
      df[[col2]] <- NULL
    }
  }
  df
}

.assign_trajectory_q <- function(df, by, p_col, q_col) {
  if (!p_col %in% colnames(df)) {
    return(df)
  }
  if (!q_col %in% colnames(df)) {
    df[[q_col]] <- NA_real_
  }
  keys <- unique(df[, c(by, p_col), drop = FALSE])
  keys <- .validate_trajectory_keys(keys, paste0(q_col, " assignment"), by)
  keys[[q_col]] <- stats::p.adjust(keys[[p_col]], method = "BH")
  df[[q_col]] <- keys[[q_col]][match(
    do.call(paste, c(df[, by, drop = FALSE], sep = "\001")),
    do.call(paste, c(keys[, by, drop = FALSE], sep = "\001"))
  )]
  df
}

.resolve_fake_identity_layer <- function(out, by) {
  out <- .demote_suffixed_cols(out, .fake_identity_detail_cols(), "fake")
  if (!"p_fake" %in% colnames(out)) {
    out$p_fake <- NA_real_
  }
  if (!"q_fake" %in% colnames(out)) {
    out$q_fake <- NA_real_
  }
  for (col in c("p_fake_edge", "p_fake_generic", "p_fake_complex", "q_fake_edge",
                "q_fake_generic", "q_fake_complex",
                "p_replace_ligand", "p_replace_receptor", "q_replace_ligand",
                "q_replace_receptor", "p_edge_iu", "q_edge_iu",
                "p_fake_network", "p_network_rewire_ligand", "p_network_rewire_receptor",
                "q_fake_network", "q_network_rewire_ligand", "q_network_rewire_receptor",
                "p_network_iu", "q_network_iu")) {
    if (!col %in% colnames(out)) {
      out[[col]] <- NA_real_
    }
  }
  for (pair in list(
    c("p_fake_edge", "q_fake_edge"),
    c("p_fake_generic", "q_fake_generic"),
    c("p_fake_complex", "q_fake_complex"),
    c("p_replace_ligand", "q_replace_ligand"),
    c("p_replace_receptor", "q_replace_receptor"),
    c("p_edge_iu", "q_edge_iu"),
    c("p_fake_network", "q_fake_network"),
    c("p_network_rewire_ligand", "q_network_rewire_ligand"),
    c("p_network_rewire_receptor", "q_network_rewire_receptor"),
    c("p_network_iu", "q_network_iu")
  )) {
    out <- .assign_trajectory_q(out, by, pair[[1L]], pair[[2L]])
  }
  out$p_identity <- ifelse(
    is.finite(out$p_network_iu), out$p_network_iu,
    ifelse(is.finite(out$p_fake_network), out$p_fake_network,
           ifelse(is.finite(out$p_edge_iu), out$p_edge_iu,
                  ifelse(is.finite(out$p_fake_edge), out$p_fake_edge,
                         ifelse(is.finite(out$p_fake_complex), out$p_fake_complex,
                                ifelse(is.finite(out$p_fake_generic), out$p_fake_generic,
                                       out$p_fake)))))
  )
  out$identity_inference_mode <- ifelse(
    is.finite(out$p_network_iu) | is.finite(out$p_fake_network), "lr_network_rewire",
    ifelse(is.finite(out$p_edge_iu) | is.finite(out$p_fake_edge), "edge_rewire",
           ifelse(is.finite(out$p_fake_generic), "context_generic",
                  ifelse(is.finite(out$p_fake), "legacy_primary", NA_character_)))
  )
  out <- .assign_trajectory_q(out, by, "p_identity", "q_identity")
  missing_q <- !is.finite(out$q_identity) & is.finite(out$q_fake)
  out$q_identity[missing_q] <- out$q_fake[missing_q]
  out
}

.fake_p_for_joint <- function(out) {
  if ("p_identity" %in% colnames(out)) {
    out$p_identity
  } else if ("p_network_iu" %in% colnames(out)) {
    ifelse(
      is.finite(out$p_network_iu), out$p_network_iu,
      ifelse(is.finite(out$p_fake_network), out$p_fake_network,
             ifelse(is.finite(out$p_edge_iu), out$p_edge_iu,
                    ifelse(is.finite(out$p_fake_edge), out$p_fake_edge,
                           ifelse(is.finite(out$p_fake_complex), out$p_fake_complex,
                                  out$p_fake))))
    )
  } else if ("p_edge_iu" %in% colnames(out)) {
    ifelse(is.finite(out$p_edge_iu), out$p_edge_iu,
           ifelse(is.finite(out$p_fake_edge), out$p_fake_edge,
                  ifelse(is.finite(out$p_fake_complex), out$p_fake_complex,
                         out$p_fake)))
  } else if ("p_fake_edge" %in% colnames(out)) {
    ifelse(is.finite(out$p_fake_edge), out$p_fake_edge,
           ifelse(is.finite(out$p_fake_complex), out$p_fake_complex, out$p_fake))
  } else if ("p_fake_complex" %in% colnames(out)) {
    ifelse(is.finite(out$p_fake_complex), out$p_fake_complex, out$p_fake)
  } else {
    out$p_fake
  }
}

.audit_permutation_blocks <- function(blocks, celltypes, lineage_weights,
                                      permutation_stratum) {
  block_tab <- as.integer(table(blocks))
  dominant_lineage <- max.col(lineage_weights, ties.method = "first")
  strata <- if (identical(permutation_stratum, "block_celltype")) {
    interaction(blocks, celltypes, drop = TRUE, lex.order = TRUE)
  } else {
    interaction(blocks, celltypes, dominant_lineage,
                drop = TRUE, lex.order = TRUE)
  }
  stratum_tab <- table(strata)
  cell_stratum_size <- as.integer(stratum_tab[as.character(strata)])
  data.frame(
    n_blocks = length(block_tab),
    min_block_size = min(block_tab),
    median_block_size = stats::median(block_tab),
    max_block_size = max(block_tab),
    singleton_block_fraction = mean(block_tab == 1L),
    n_shuffle_strata = length(stratum_tab),
    min_shuffle_stratum_size = min(stratum_tab),
    median_shuffle_stratum_size = stats::median(as.integer(stratum_tab)),
    max_shuffle_stratum_size = max(stratum_tab),
    singleton_shuffle_stratum_fraction = mean(stratum_tab == 1L),
    non_permutable_cell_fraction = mean(cell_stratum_size <= 1L),
    stringsAsFactors = FALSE
  )
}

.provenance_hash <- function(parts) {
  tmp <- tempfile()
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(parts, tmp)
  paste0("llr_", unname(tools::md5sum(tmp)))
}

.lr_db_hash_parts <- function(lr_db) {
  list(
    ligand_rows = lapply(lr_db$ligand_rows, as.integer),
    receptor_rows = lapply(lr_db$receptor_rows, as.integer)
  )
}

.compute_analysis_data_hash <- function(expr, celltypes, lr_db, pseudotime,
                                        lineage_weights, block_labels = NULL) {
  ct <- .assert_celltypes(celltypes, ncol(expr))
  lw <- .assert_lineage_weights(lineage_weights, ncol(expr))
  pt <- .assert_pseudotime_1d(pseudotime, ncol(expr))
  blocks <- .prepare_block_labels(block_labels, ncol(expr))
  lr_db <- .validate_lr_db(lr_db)
  .provenance_hash(list(
    tag = "analysis_data_v2",
    expr_dim = as.integer(dim(expr)),
    expr_p = as.integer(expr@p),
    expr_i = as.integer(expr@i),
    expr_x = as.numeric(expr@x),
    celltypes = ct,
    celltype_levels = attr(ct, "levels"),
    lr_db = .lr_db_hash_parts(lr_db),
    pseudotime = pt,
    lineage_weights = as.numeric(lw),
    block_labels = blocks,
    block_levels = attr(blocks, "levels")
  ))
}

.compute_scoring_config_hash <- function(n_bins, min_pct, min_n_eff,
                                         min_valid_bins, min_valid_fraction,
                                         stat_code, top_k) {
  .provenance_hash(list(
    tag = "scoring_config_v1",
    n_bins = n_bins,
    min_pct = min_pct,
    min_n_eff = min_n_eff,
    min_valid_bins = min_valid_bins,
    min_valid_fraction = min_valid_fraction,
    stat_code = stat_code,
    top_k = top_k
  ))
}

.compute_fake_null_config_hash <- function(candidate_fake_genes, fake_null_mode,
                                           n_fake_per_real, seed,
                                           max_match_distance, block_labels,
                                           n_cells, complex_null_mode = "skip",
                                           n_strata, caliper_k,
                                           min_fake_observed_generic,
                                           min_fake_coverage_generic,
                                           min_fake_observed_edge,
                                           min_fake_coverage_edge,
                                           min_unique_fake_draws,
                                           min_unique_fake_fraction) {
  blocks <- .prepare_block_labels(block_labels, n_cells)
  .provenance_hash(list(
    tag = "fake_null_config_v3",
    candidate_fake_genes = as.integer(candidate_fake_genes),
    fake_null_mode = fake_null_mode,
    complex_null_mode = complex_null_mode,
    n_fake_per_real = n_fake_per_real,
    n_strata = n_strata,
    caliper_k = caliper_k,
    seed = seed,
    max_match_distance = max_match_distance,
    min_fake_observed_generic = min_fake_observed_generic,
    min_fake_coverage_generic = min_fake_coverage_generic,
    min_fake_observed_edge = min_fake_observed_edge,
    min_fake_coverage_edge = min_fake_coverage_edge,
    min_unique_fake_draws = min_unique_fake_draws,
    min_unique_fake_fraction = min_unique_fake_fraction,
    block_labels = blocks,
    block_levels = attr(blocks, "levels")
  ))
}

.compute_temporal_null_config_hash <- function(n_perms, stat, block_labels,
                                               permutation_stratum, seed,
                                               n_cells) {
  blocks <- .prepare_block_labels(block_labels, n_cells)
  .provenance_hash(list(
    tag = "temporal_null_config_v1",
    n_perms = n_perms,
    stat = stat,
    block_labels = blocks,
    block_levels = attr(blocks, "levels"),
    permutation_stratum = permutation_stratum,
    seed = seed
  ))
}

.compute_role_config_hash <- function(tau_0, h_fork, min_n_eff_cells,
                                      transition_window = 0L,
                                      timing_source = NA_character_,
                                      bif_df = NULL,
                                      evidence_hashes = NULL) {
  .provenance_hash(list(
    tag = "role_config_v3",
    tau_0 = tau_0,
    h_fork = h_fork,
    min_n_eff_cells = min_n_eff_cells,
    transition_window = as.integer(transition_window),
    timing_source = timing_source,
    bifurcation = bif_df,
    evidence_hashes = evidence_hashes
  ))
}

.attach_provenance_columns <- function(df, provenance) {
  for (nm in names(provenance)) {
    df[[nm]] <- provenance[[nm]]
  }
  df
}

.scoring_hash_from_df <- function(df) {
  if ("scoring_config_hash" %in% colnames(df)) {
    unique(as.character(df$scoring_config_hash))
  } else if ("calibration_config_hash" %in% colnames(df)) {
    unique(as.character(df$calibration_config_hash))
  } else {
    character(0)
  }
}

.analysis_hash_from_df <- function(df) {
  if ("analysis_data_hash" %in% colnames(df)) {
    unique(as.character(df$analysis_data_hash))
  } else {
    character(0)
  }
}

.derive_entropy_semantics <- function(scores_df) {
  if ("fate_entropy" %in% colnames(scores_df)) {
    fate <- scores_df$fate_entropy
  } else if ("lineage_entropy" %in% colnames(scores_df)) {
    fate <- scores_df$lineage_entropy
  } else {
    stop("`scores_df` must contain `lineage_entropy` or `fate_entropy`",
         call. = FALSE)
  }
  if ("assignment_confidence" %in% colnames(scores_df)) {
    conf <- scores_df$assignment_confidence
    conf_source <- rep("explicit", nrow(scores_df))
  } else if ("lineage_weight_max" %in% colnames(scores_df)) {
    conf <- pmax(0, pmin(1, scores_df$lineage_weight_max))
    conf_source <- rep("lineage_weight_max", nrow(scores_df))
  } else {
    conf <- rep(NA_real_, nrow(scores_df))
    conf_source <- rep("missing", nrow(scores_df))
  }
  list(
    assignment_confidence = conf,
    assignment_confidence_source = conf_source,
    fate_entropy = fate,
    fate_ambiguity = fate
  )
}

.compute_fake_matching_balance <- function(fake_pairs) {
  empty <- data.frame(
    pair_idx = integer(),
    source = character(),
    target = character(),
    lineage = integer(),
    n_draws = integer(),
    n_generic = integer(),
    n_edge = integer(),
      n_rewire_receptor = integer(),
      n_rewire_ligand = integer(),
      n_network_rewire_receptor = integer(),
      n_network_rewire_ligand = integer(),
      edge_fraction = double(),
    mean_ligand_match_distance = double(),
    mean_receptor_match_distance = double(),
    max_ligand_match_distance = double(),
    max_receptor_match_distance = double(),
    stringsAsFactors = FALSE
  )
  if (is.null(fake_pairs) || nrow(fake_pairs) == 0L) {
    return(empty)
  }
  fp <- fake_pairs
  if (!"null_kind" %in% colnames(fp)) {
    fp$null_kind <- 0L
  }
  fp$null_kind <- as.integer(fp$null_kind)
  fp$pair_idx <- as.integer(fp$shadowing_pair)
  keys <- unique(fp[, c("pair_idx", "source", "target", "lineage"), drop = FALSE])
  rows <- lapply(seq_len(nrow(keys)), function(i) {
    k <- keys[i, , drop = FALSE]
    sel <- fp$pair_idx == k$pair_idx &
      fp$source == k$source &
      fp$target == k$target &
      fp$lineage == k$lineage
    sub <- fp[sel, , drop = FALSE]
    n_draws <- nrow(sub)
    data.frame(
      pair_idx = k$pair_idx,
      source = k$source,
      target = k$target,
      lineage = k$lineage,
      n_draws = n_draws,
      n_generic = sum(sub$null_kind == 0L),
      n_edge = sum(sub$null_kind != 0L),
      n_rewire_receptor = sum(sub$null_kind == 1L),
      n_rewire_ligand = sum(sub$null_kind == 2L),
      n_network_rewire_receptor = sum(sub$null_kind == 3L),
      n_network_rewire_ligand = sum(sub$null_kind == 4L),
      edge_fraction = if (n_draws > 0L) sum(sub$null_kind != 0L) / n_draws else NA_real_,
      mean_ligand_match_distance = mean(sub$ligand_match_distance),
      mean_receptor_match_distance = mean(sub$receptor_match_distance),
      max_ligand_match_distance = max(sub$ligand_match_distance),
      max_receptor_match_distance = max(sub$receptor_match_distance),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

#' Compute pseudotime from HNSW kNN with UMAP-fuzzy edge weights.
#'
#' @param nn_idx integer matrix \[n_cells × k\] of 1-based neighbour ids.
#' @param nn_dists numeric matrix \[n_cells × k\] of L2 distances.
#' @param start_cell integer scalar (1-based root cell index).
#' @param rank_normalize logical; if TRUE, rank-normalise finite pseudotime to the 0--1 range.
#' @return numeric vector of length `n_cells`; pseudotime in the range 0 to 1.
#'   Cells disconnected from `start_cell` receive `NaN`.
#' @examples
#' set.seed(1)
#' n <- 30L; k <- 4L
#' emb <- cbind(seq(0, 1, length.out = n), rnorm(n))
#' nn_idx <- matrix(0L, n, k); nn_dists <- matrix(0, n, k)
#' for (i in seq_len(n)) {
#'   d <- sqrt(rowSums((emb - matrix(emb[i, ], n, 2, byrow = TRUE))^2))
#'   o <- order(d)[2:(k + 1L)]
#'   nn_idx[i, ] <- o; nn_dists[i, ] <- d[o]
#' }
#' pt <- ComputePseudotimeFuzzy(nn_idx, nn_dists, start_cell = 1L)
#' stopifnot(length(pt) == n, pt[1L] == 0)
#' @export
ComputePseudotimeFuzzy <- function(nn_idx, nn_dists, start_cell, rank_normalize = FALSE) {
  knn <- .prepare_knn_inputs(nn_idx, nn_dists)
  n_cells <- nrow(knn$indices)
  k       <- ncol(knn$indices)
  .assert_indices_in_range(start_cell, 1L, n_cells, "start_cell")
  start_cell <- .as_integerish(start_cell, "start_cell")
  if (!is.logical(rank_normalize) || length(rank_normalize) != 1L ||
      is.na(rank_normalize)) {
    stop("`rank_normalize` must be TRUE or FALSE", call. = FALSE)
  }
  pt <- rs_compute_pseudotime_fuzzy(
    as.integer(t(knn$indices)),
    as.numeric(t(knn$distances)),
    as.integer(n_cells),
    as.integer(k),
    as.integer(start_cell)
  )
  if (isTRUE(rank_normalize)) {
    ok <- is.finite(pt)
    if (sum(ok) >= 2L) {
      r <- rank(pt[ok], ties.method = "average")
      pt[ok] <- (r - 1) / (length(r) - 1)
    }
  }
  pt
}

ComputePseudotimeFuzzyMulti <- function(nn_idx, nn_dists, start_cells,
                                        rank_normalize = FALSE) {
  knn <- .prepare_knn_inputs(nn_idx, nn_dists)
  n_cells <- nrow(knn$indices)
  k <- ncol(knn$indices)
  .assert_indices_in_range(start_cells, 1L, n_cells, "start_cells")
  start_cells <- .as_integerish(start_cells, "start_cells")
  if (!is.logical(rank_normalize) || length(rank_normalize) != 1L ||
      is.na(rank_normalize)) {
    stop("`rank_normalize` must be TRUE or FALSE", call. = FALSE)
  }
  pt <- rs_compute_pseudotime_fuzzy_multi(
    as.integer(t(knn$indices)),
    as.numeric(t(knn$distances)),
    as.integer(n_cells),
    as.integer(k),
    as.integer(start_cells)
  )
  pt <- matrix(pt, nrow = n_cells, ncol = length(start_cells))
  if (isTRUE(rank_normalize)) {
    pt <- apply(pt, 2L, function(values) {
      ok <- is.finite(values)
      if (sum(ok) >= 2L) {
        r <- rank(values[ok], ties.method = "average")
        values[ok] <- (r - 1) / (length(r) - 1)
      }
      values
    })
  }
  pt
}

#' Return the realized fuzzy and raw costs on one connected kNN graph.
#'
#' This internal diagnostic uses the same Rust graph construction as
#' `ComputePseudotimeFuzzy()`. It exposes each symmetrized undirected edge once
#' so that double-precision shortest paths can be compared without rebuilding
#' the fuzzy membership calculation in R.
#' @keywords internal
ComputePseudotimeCostEdges <- function(nn_idx, nn_dists) {
  knn <- .prepare_knn_inputs(nn_idx, nn_dists)
  n_cells <- nrow(knn$indices)
  k <- ncol(knn$indices)
  out <- rs_hnsw_cost_edges(
    as.integer(t(knn$indices)),
    as.numeric(t(knn$distances)),
    as.integer(n_cells),
    as.integer(k)
  )
  data.frame(
    from = as.integer(out$from),
    to = as.integer(out$to),
    fuzzy_cost = as.numeric(out$fuzzy_cost),
    raw_distance = as.numeric(out$raw_distance)
  )
}

#' Build a deterministic coarse clustering for lineage topology.
#'
#' Root selection and pseudotime can require a finer clustering than lineage
#' topology. This helper constructs the latter independently so that an MST is
#' not forced to treat small local partitions as separate biological states.
#' It uses only the supplied embedding and never accepts ordering or branch
#' annotations.
#'
#' @param embedding numeric matrix \[n_cells × n_dims\].
#' @param k requested number of topology clusters.
#' @param nstart number of random k-means initializations.
#' @param iter_max maximum Lloyd iterations per initialization.
#' @param seed integer random seed.
#' @return list with compact 1-based `labels`, `requested_k`, `occupied_k`,
#'   occupied-cluster `centers`, `tot_withinss`, and `algorithm`.
#' @export
ClusterLineageTopology <- function(embedding, k = 15L, nstart = 10L,
                                   iter_max = 500L, seed = 12345L) {
  emb <- .prepare_embedding(embedding)
  k <- .validate_positive_int(k, "k", min = 2L)
  nstart <- .validate_positive_int(nstart, "nstart")
  iter_max <- .validate_positive_int(iter_max, "iter_max")
  seed <- .as_integerish(seed, "seed")
  if (length(seed) != 1L) {
    stop("`seed` must be a finite integer scalar", call. = FALSE)
  }
  if (k > nrow(emb)) {
    stop("`k` cannot exceed nrow(embedding)", call. = FALSE)
  }

  set.seed(seed)
  empty_start_warnings <- 0L
  fit <- withCallingHandlers(
    stats::kmeans(
      emb,
      centers = k,
      nstart = nstart,
      iter.max = iter_max,
      algorithm = "Lloyd"
    ),
    warning = function(condition) {
      is_empty_start <- identical(
        conditionMessage(condition),
        "empty cluster: try a better set of initial centers"
      )
      if (is_empty_start) {
        empty_start_warnings <<- empty_start_warnings + 1L
        invokeRestart("muffleWarning")
      }
    }
  )
  occupied_levels <- sort(unique(fit$cluster))
  labels <- as.integer(match(fit$cluster, occupied_levels))

  list(
    labels = labels,
    requested_k = as.integer(k),
    occupied_k = as.integer(length(occupied_levels)),
    centers = fit$centers[occupied_levels, , drop = FALSE],
    tot_withinss = as.numeric(fit$tot.withinss),
    algorithm = "Lloyd",
    empty_start_warnings = as.integer(empty_start_warnings)
  )
}

#' Extract root→leaf lineage paths from cluster centroids via MST.
#'
#' @param embedding numeric matrix \[n_cells × n_dims\].
#' @param cluster_labels integer vector of length n_cells, 1-based.
#' @param start_cluster integer scalar (1-based root cluster).
#' @return list of integer vectors; each element is one lineage path
#'   (1-based cluster ids from root to leaf).
#' @export
ExtractLineages <- function(embedding, cluster_labels, start_cluster) {
  emb <- .prepare_embedding(embedding)
  if (length(cluster_labels) != nrow(emb)) {
    stop("`cluster_labels` length must equal nrow(embedding)", call. = FALSE)
  }
  cl_map <- .compact_cluster_labels(cluster_labels)
  cl <- cl_map$labels
  n_cells    <- nrow(emb)
  n_dims     <- ncol(emb)
  n_clusters <- length(cl_map$levels)
  start_cluster <- .compact_cluster_id(start_cluster, cl_map$levels)
  if (is.na(start_cluster)) stop("`start_cluster` is not present in `cluster_labels`")
  out <- rs_extract_lineages(
    as.numeric(t(emb)),
    cl,
    as.integer(n_cells),
    as.integer(n_dims),
    as.integer(n_clusters),
    as.integer(start_cluster)
  )
  lapply(out, function(path) as.integer(cl_map$levels[as.integer(path)]))
}

#' Soft per-cell lineage weights along each lineage's polyline.
#'
#' Internally re-runs cluster centroids, MST extraction, and lineage paths.
#' Cluster membership on a lineage path receives a hard weight of 1.0.
#' For geometry-based weighting on shared trunk clusters, use
#' `ComputeLineageWeightsTrunkAware`.
#'
#' @param embedding numeric matrix \[n_cells × n_dims\].
#' @param cluster_labels integer vector of length n_cells, 1-based.
#' @param start_cluster integer scalar (1-based root cluster).
#' @param sigma numeric scalar, Gaussian bandwidth on lineage-path distance.
#' @return numeric matrix \[n_cells × n_lineages\], rows sum to 1.
#' @export
ComputeLineageWeights <- function(embedding, cluster_labels, start_cluster, sigma = 1.0) {
  emb <- .prepare_embedding(embedding)
  sigma <- .validate_positive_num(sigma, "sigma")
  if (length(cluster_labels) != nrow(emb)) {
    stop("`cluster_labels` length must equal nrow(embedding)", call. = FALSE)
  }
  cl_map <- .compact_cluster_labels(cluster_labels)
  cl <- cl_map$labels
  n_cells    <- nrow(emb)
  n_dims     <- ncol(emb)
  n_clusters <- length(cl_map$levels)
  start_cluster <- .compact_cluster_id(start_cluster, cl_map$levels)
  if (is.na(start_cluster)) stop("`start_cluster` is not present in `cluster_labels`")
  flat <- rs_compute_lineage_weights(
    as.numeric(t(emb)),
    cl,
    as.integer(n_cells),
    as.integer(n_dims),
    as.integer(n_clusters),
    as.integer(start_cluster),
    as.numeric(sigma)
  )
  n_lineages <- length(flat) %/% n_cells
  matrix(flat, nrow = n_cells, ncol = n_lineages, byrow = TRUE)
}

#' Fate probabilities via absorbing-state Markov power iteration.
#'
#' @param nn_idx integer matrix \[n_cells × k\] of 1-based neighbour ids.
#' @param nn_dists numeric matrix \[n_cells × k\].
#' @param pseudotime numeric vector of length n_cells; `NaN`/`NA` cells are
#'   excluded from diffusion (matches Palantir).
#' @param terminal_cells integer vector of 1-based terminal cell indices.
#' @param bandwidth numeric scalar, sigmoid soft-gate bandwidth on
#'   pseudotime difference.
#' @param max_iter,tol convergence parameters.
#' @return numeric matrix \[n_cells × length(terminal_cells)\] of fate
#'   probabilities (rows sum to 1 over reachable terminals).
#' @export
ComputeFateProbabilities <- function(nn_idx, nn_dists, pseudotime, terminal_cells,
                                     bandwidth = 0.05, max_iter = 500L, tol = 1e-4) {
  knn <- .prepare_knn_inputs(nn_idx, nn_dists)
  n_cells     <- nrow(knn$indices)
  k           <- ncol(knn$indices)
  pseudotime <- .assert_pseudotime_1d(pseudotime, n_cells)
  bandwidth <- .validate_positive_num(bandwidth, "bandwidth")
  max_iter <- .validate_positive_int(max_iter, "max_iter")
  tol <- .validate_positive_num(tol, "tol")
  terminal_cells <- unique(.as_integerish(terminal_cells, "terminal_cells"))
  if (length(terminal_cells) < 1L) {
    stop("`terminal_cells` must contain at least one cell", call. = FALSE)
  }
  .assert_indices_in_range(terminal_cells, 1L, n_cells, "terminal_cells")
  n_terminals <- length(terminal_cells)
  flat <- rs_compute_fate_probabilities(
    as.integer(t(knn$indices)),
    as.numeric(t(knn$distances)),
    pseudotime,
    as.integer(terminal_cells),
    as.integer(n_cells),
    as.integer(k),
    as.numeric(bandwidth),
    as.integer(max_iter),
    as.numeric(tol)
  )
  matrix(flat, nrow = n_cells, ncol = n_terminals, byrow = TRUE)
}

# Internal: validate the lr_db structure expected by ComputeLRScores /
# ComputeTemporalLR.  keeps the contract minimal — gene-symbol → row
# resolution is intentionally not done here; it belongs in a full
# user-facing LR-database parser.
.validate_lr_db <- function(lr_db) {
  if (!is.list(lr_db) ||
      is.null(lr_db$ligand_rows) || is.null(lr_db$receptor_rows) ||
      !is.list(lr_db$ligand_rows) || !is.list(lr_db$receptor_rows) ||
      length(lr_db$ligand_rows) != length(lr_db$receptor_rows)) {
    stop("`lr_db` must be a list with `ligand_rows` and `receptor_rows`, ",
         "each a list of integer vectors (1-based gene row indices), ",
         "and the two sub-lists must have equal length.")
  }
  ligand_rows <- lapply(
    lr_db$ligand_rows,
    function(rows) unique(.as_integerish(rows, "lr_db$ligand_rows"))
  )
  receptor_rows <- lapply(
    lr_db$receptor_rows,
    function(rows) unique(.as_integerish(rows, "lr_db$receptor_rows"))
  )
  if (length(ligand_rows) < 1L ||
      any(lengths(ligand_rows) < 1L) ||
      any(lengths(receptor_rows) < 1L)) {
    stop("`lr_db` must contain at least one LR pair and no empty complexes",
         call. = FALSE)
  }
  lr_db$ligand_rows <- ligand_rows
  lr_db$receptor_rows <- receptor_rows
  lr_db
}

#' Cluster-level LR communication scores with permutation p-values.
#'
#' @param expr `dgCMatrix` (genes × cells), log-normalised.
#' @param celltypes integer vector of length ncol(expr), 1-based cluster id.
#' @param lr_db list with `ligand_rows` (list of int) and `receptor_rows`
#'   (list of int), one element per LR pair, each holding 1-based row
#'   indices into `expr` for the ligand / receptor subunits (length 1 for
#'   simple pairs; > 1 for multi-subunit complexes).
#' @param n_perms integer.
#' @param seed numeric (R 32-bit integers cannot hold a full u64 seed).
#' @param min_pct numeric, minimum per-cluster expressing fraction.
#' @return data.frame with columns source/target/pair_idx (1-based),
#'   score, pvalue, pct_ligand, pct_receptor.
#' @export
ComputeLRScores <- function(expr, celltypes, lr_db,
                            n_perms = 500L, seed = 12345L, min_pct = 0.05) {
  n_perms <- .validate_positive_int(n_perms, "n_perms")
  seed <- .validate_seed(seed)
  min_pct <- .validate_unit_interval(min_pct, "min_pct")
  expr <- .prepare_lr_expr(expr)
  lr_db <- .validate_lr_db(lr_db)
  .assert_lr_db_indices(lr_db, nrow(expr))
  ct      <- .assert_celltypes(celltypes, ncol(expr))
  ct_levels <- attr(ct, "levels")
  n_types <- max(ct)
  out <- rs_compute_lr_scores(
    expr,
    ct,
    as.integer(n_types),
    lapply(lr_db$ligand_rows,   as.integer),
    lapply(lr_db$receptor_rows, as.integer),
    n_perms,
    seed,
    min_pct
  )
  out <- data.frame(out, stringsAsFactors = FALSE)
  out$source <- ct_levels[out$source]
  out$target <- ct_levels[out$target]
  attr(out, "celltype_levels") <- ct_levels
  out
}

#' Pseudotime-binned LR scores per (lineage, bin).
#'
#' Returns `lineage_entropy = max(mean H(source), mean H(target))` per row.
#' For per-lineage entropy curves, use `ComputeTemporalLRWithIntermediate`,
#' then `ComputeLineageEntropyRaw` and `SmoothMA3`.
#'
#' @param expr `dgCMatrix` (genes × cells), log-normalised.
#' @param celltypes integer vector of length ncol(expr), 1-based.
#' @param lr_db see `ComputeLRScores`.
#' @param pseudotime numeric vector of length ncol(expr); `NA`/`NaN` cells are
#'   excluded from binning.
#' @param lineage_weights numeric matrix \[n_cells × n_lineages\] (output
#'   of `ComputeLineageWeights`).
#' @param n_bins integer, number of pseudotime bins per lineage.
#' @param min_pct numeric, minimum *weighted* per-bin expressing fraction.
#' @param min_n_eff minimum Kish effective sample size required for a
#'   structurally valid bin. Use the same value as downstream FakeLR and
#'   temporal permutation calibration (default 5).
#' @param return_complete_grid when `TRUE` (default), return every bin for each
#'   prescreened trajectory. Structural-invalid bins have `score = NA`; valid
#'   bins below the expression gate have `score = 0`.
#' @param requested_contexts optional data.frame with columns `pair_idx`,
#'   `source`, and `target` (1-based cell-type indices). When supplied with
#'   `return_complete_grid = TRUE`, the complete grid is expanded for this
#'   pooled hypothesis family instead of donor-local prescreen contexts.
#' @param spatial_adjacency optional `n_types` by `n_types` matrix from
#'   `ComputeSpatialAdjacency()`; when supplied, adds `spatial_support` and
#'   `spatially_eligible` columns to the score table.
#' @param min_spatial_support minimum adjacency value for `spatially_eligible`
#'   when `spatial_adjacency` is provided (default 0.05).
#' @param common_grid when `TRUE`, score every lineage on one shared
#'   boundary vector computed from the pooled on-trajectory mass (a cell's
#'   mass is its maximum lineage weight).  `FALSE` (default) reproduces the
#'   legacy per-lineage equal-frequency grid, in which a bin index is a
#'   lineage-relative quantile and is not time-comparable across branches.
#'   Fork-local contrasts require the shared coordinate.
#' @return data.frame with columns pair_idx/source/target/lineage/bin
#'   (1-based), score, `n_eff_cells` (Kish ESS), `weight_mass_cells`,
#'   lineage_entropy, `structural_valid`, and `expression_supported`.
#'   Attribute `temporal_coordinate` records the grid mode
#'   (`"common_pseudotime_grid"` or `"lineage_local_equal_frequency_bins"`).
#' @export
ComputeTemporalLR <- function(expr, celltypes, lr_db, pseudotime, lineage_weights,
                              n_bins = 12L, min_pct = 0.05,
                              min_n_eff = 5.0,
                              return_complete_grid = TRUE,
                              requested_contexts = NULL,
                              spatial_adjacency = NULL,
                              min_spatial_support = 0.05,
                              common_grid = FALSE) {
  n_bins <- .validate_positive_int(n_bins, "n_bins", min = 2L)
  min_pct <- .validate_unit_interval(min_pct, "min_pct")
  min_n_eff <- .validate_positive_num(min_n_eff, "min_n_eff")
  if (!is.logical(common_grid) || length(common_grid) != 1L ||
      is.na(common_grid)) {
    stop("`common_grid` must be TRUE or FALSE", call. = FALSE)
  }
  expr <- .prepare_lr_expr(expr)
  lr_db <- .validate_lr_db(lr_db)
  .assert_lr_db_indices(lr_db, nrow(expr))
  lw <- .assert_lineage_weights(lineage_weights, ncol(expr))
  n_lineages <- ncol(lw)
  ct      <- .assert_celltypes(celltypes, ncol(expr))
  ct_levels <- attr(ct, "levels")
  n_types <- max(ct)
  pt_vec  <- .assert_pseudotime_1d(pseudotime, ncol(expr))
  req_ctx <- if (is.null(requested_contexts)) {
    list(pair_idx = integer(), source = integer(), target = integer())
  } else {
  .pooled_lr_context_vectors(requested_contexts, ct_levels)
  }
  out <- rs_compute_temporal_lr(
    expr,
    ct,
    as.integer(n_types),
    lapply(lr_db$ligand_rows,   as.integer),
    lapply(lr_db$receptor_rows, as.integer),
    pt_vec,
    as.numeric(t(lw)),
    as.integer(n_lineages),
    n_bins,
    min_pct,
    min_n_eff,
    isTRUE(return_complete_grid),
    req_ctx$pair_idx,
    req_ctx$source,
    req_ctx$target,
    common_grid
  )
  out <- data.frame(out, stringsAsFactors = FALSE)
  out$structural_valid <- as.logical(out$structural_valid)
  out$expression_supported <- as.logical(out$expression_supported)
  out$source <- ct_levels[out$source]
  out$target <- ct_levels[out$target]
  attr(out, "celltype_levels") <- ct_levels
  attr(out, "temporal_coordinate") <- if (common_grid) {
    "common_pseudotime_grid"
  } else {
    "lineage_local_equal_frequency_bins"
  }
  .attach_spatial_columns(out, spatial_adjacency, min_spatial_support)
}

#' Trunk-aware soft per-cell lineage weights.
#'
#' Geometry-based lineage weights for shared trunk clusters.  Clusters that
#' appear on at least two lineage paths use the Gaussian distance weight
#' `exp(-d^2/2σ^2)`, while single-lineage clusters on their own path retain a
#' hard weight of 1.0 and off-path cells use the Gaussian weight.
#'
#' @param embedding,cluster_labels,start_cluster,sigma see `ComputeLineageWeights`.
#' @return numeric matrix \[n_cells × n_lineages\], rows sum to 1.
#' @export
ComputeLineageWeightsTrunkAware <- function(embedding, cluster_labels, start_cluster, sigma = 1.0) {
  emb <- .prepare_embedding(embedding)
  sigma <- .validate_positive_num(sigma, "sigma")
  if (length(cluster_labels) != nrow(emb)) {
    stop("`cluster_labels` length must equal nrow(embedding)", call. = FALSE)
  }
  cl_map <- .compact_cluster_labels(cluster_labels)
  cl <- cl_map$labels
  n_cells    <- nrow(emb)
  n_dims     <- ncol(emb)
  n_clusters <- length(cl_map$levels)
  start_cluster <- .compact_cluster_id(start_cluster, cl_map$levels)
  if (is.na(start_cluster)) stop("`start_cluster` is not present in `cluster_labels`")
  flat <- rs_compute_lineage_weights_trunk_aware(
    as.numeric(t(emb)),
    cl,
    as.integer(n_cells),
    as.integer(n_dims),
    as.integer(n_clusters),
    as.integer(start_cluster),
    as.numeric(sigma)
  )
  # Rust returns row-major flat [n_cells * n_lineages]; reshape to col-major
  # R matrix with rows = cells, cols = lineages (same layout as
  # ComputeLineageWeights).
  n_lineages <- length(flat) / n_cells
  matrix(flat, nrow = n_cells, ncol = n_lineages, byrow = TRUE)
}

#' Branch-local lineage weights from expression residuals.
#'
#' At each fork, expression relative to the cell's own cluster center is
#' projected onto outgoing branch directions. Cells already downstream of a
#' branch retain that assignment. Paths must be rooted and ordered.
#' With priming enabled, fork-centered expression gives continuous weights
#' within and between downstream states.
#'
#' @param embedding numeric matrix of cells by embedding dimensions.
#' @param cluster_labels integer cluster label per cell.
#' @param lineage_paths rooted cluster paths in output lineage order.
#' @param priming use continuous fork-centered weights for fate tendency.
#' @return numeric matrix of cells by lineages, with rows summing to one.
#' @export
ComputeLineageWeightsLocalBranchResidual <- function(embedding, cluster_labels,
                                                     lineage_paths, priming = FALSE) {
  emb <- .prepare_embedding(embedding)
  stopifnot(length(cluster_labels) == nrow(emb))
  cl_map <- .compact_cluster_labels(cluster_labels)
  paths <- lapply(lineage_paths, match, table = cl_map$levels)
  stopifnot(length(paths) > 0L, all(lengths(paths) > 0L),
            all(!is.na(unlist(paths))))
  flat <- rs_compute_lineage_weights_local_branch_residual(
    as.numeric(t(emb)), cl_map$labels,
    as.integer(nrow(emb)), as.integer(ncol(emb)),
    as.integer(length(cl_map$levels)), paths, as.logical(priming)
  )
  matrix(flat, nrow = nrow(emb), ncol = length(paths), byrow = TRUE)
}

ComputeIntermediateLineagePriming <- function(embedding, cluster_labels,
                                               lineage_paths, lineage_weights) {
  leaves <- vapply(lineage_paths, tail, integer(1L), n = 1L)
  nodes <- data.table(state = unlist(lineage_paths),
    lineage = rep(seq_along(lineage_paths), lengths(lineage_paths)))[,
      .(descendants = uniqueN(lineage)), by = state]
  intermediate <- nodes[descendants == 1L & !state %in% leaves, state]
  rows <- which(cluster_labels %in% intermediate)
  continuous <- ComputeLineageWeightsLocalBranchResidual(embedding, cluster_labels,
    lineage_paths, priming = TRUE)
  weights <- as.matrix(lineage_weights)
  weights[rows, ] <- continuous[rows, , drop = FALSE]
  weights
}

ComputeLineageProgrammeEvidence <- function(embedding, programme_labels, root_cells,
                                               lineage_weights, reference_weights,
                                               supported_terminals,
                                               condition = rep("all", nrow(embedding)),
                                               condition_support = matrix(TRUE,
                                                 nrow = length(unique(condition)),
                                                 ncol = length(supported_terminals),
                                                 dimnames = list(unique(condition),
                                                   as.character(supported_terminals))),
                                               query_cells = root_cells,
                                               calibration_context = NULL) {
  embedding <- as.matrix(embedding)
  anchor_rows <- which(!is.na(programme_labels))
  terminal_levels <- sort(unique(programme_labels[anchor_rows]))
  anchor_code <- match(programme_labels[anchor_rows], terminal_levels)
  centers <- rowsum(embedding[anchor_rows, , drop = FALSE], anchor_code) /
    tabulate(anchor_code, length(terminal_levels))
  residual <- embedding[anchor_rows, , drop = FALSE] - centers[anchor_code, , drop = FALSE]
  covariance <- crossprod(residual) / (length(anchor_rows) - length(terminal_levels))
  root_rows <- which(root_cells)
  origin <- colMeans(embedding[root_rows, , drop = FALSE])
  selected <- match(supported_terminals, terminal_levels)
  delta <- sweep(centers[selected, , drop = FALSE], 2L, origin, "-")
  directions <- solve(covariance, t(delta))
  distance <- rowSums(delta * t(directions))
  score <- sweep(embedding[root_rows, , drop = FALSE], 2L, origin, "-") %*% directions
  context_levels <- rownames(condition_support)
  context_values <- condition
  switch(as.character(!is.null(calibration_context)), "TRUE" = {
    context_values <- calibration_context
    context_levels <- sort(unique(context_values))
  })
  context <- match(context_values, context_levels)
  root_context <- context[root_rows]
  context_size <- tabulate(root_context, length(context_levels))
  context_mean <- rowsum(score, root_context) / context_size
  centered <- score - context_mean[root_context, , drop = FALSE]
  context_variance <- rowsum(centered^2, root_context) / context_size
  query_rows <- which(query_cells)
  query_context <- context[query_rows]
  query_score <- sweep(embedding[query_rows, , drop = FALSE], 2L, origin, "-") %*% directions
  query_score <- query_score - context_mean[query_context, , drop = FALSE]
  context_distance <- matrix(distance, nrow = length(context_levels),
    ncol = length(supported_terminals), byrow = TRUE)
  switch(as.character(!is.null(calibration_context)), "TRUE" = {
    context_distance <- context_distance - context_mean
  })
  noise_sd <- sqrt(context_variance[query_context, , drop = FALSE])
  signed_distance <- context_distance[query_context, , drop = FALSE]
  standardized <- query_score * sign(signed_distance) / noise_sd
  segment_length <- abs(signed_distance) / noise_sd
  lower <- -standardized
  upper <- segment_length - standardized
  reflect <- lower > 0
  interval_upper <- upper
  interval_lower <- lower
  interval_upper[reflect] <- -lower[reflect]
  interval_lower[reflect] <- -upper[reflect]
  log_upper <- pnorm(interval_upper, log.p = TRUE)
  log_interval <- log_upper + log(-expm1(pnorm(interval_lower, log.p = TRUE) - log_upper))
  log_evidence <- standardized^2 / 2 + log(2 * pi) / 2 -
    log(segment_length) + log_interval
  normalizer <- pmax(0, matrixStats::rowMaxs(log_evidence))
  ancestral <- exp(-normalizer)
  excess <- pmax(exp(log_evidence - normalizer) - ancestral, 0)
  support <- condition_support[, match(as.character(supported_terminals),
    colnames(condition_support)), drop = FALSE]
  support_context <- match(condition, rownames(support))
  excess <- excess * support[support_context[query_rows], , drop = FALSE]
  total <- ancestral + rowSums(excess)
  ancestral <- ancestral / total
  activated <- excess / total
  weights <- as.matrix(lineage_weights)
  weights[query_rows, ] <- weights[query_rows, , drop = FALSE] * ancestral
  programme <- matrix(0, nrow(embedding), length(supported_terminals))
  programme[query_rows, ] <- activated
  programme_reference <- programme * 0
  reference_rows <- which(programme_labels %in% supported_terminals)
  programme_reference[cbind(reference_rows,
    match(programme_labels[reference_rows], supported_terminals))] <- 1
  programme_reference <- programme_reference * support[support_context, , drop = FALSE]
  list(weights = cbind(weights, programme),
    reference_weights = cbind(as.matrix(reference_weights), programme_reference),
    programme_terminals = supported_terminals,
    root_rows = root_rows, query_rows = query_rows,
    activation = cbind(ancestral = ancestral, activated),
    model = list(centers = centers, terminal_levels = terminal_levels,
      covariance = covariance, origin = origin, directions = directions,
      distance = distance, context_levels = context_levels, context_mean = context_mean,
      context_variance = context_variance, context_support = support))
}

ComputeLineageStateProgrammeReadout <- function(embedding, programme_labels, cluster_labels,
                                                  lineage_paths, lineage_weights,
                                                  reference_weights, supported_terminals,
                                                  program_profiles, adjacency,
                                                  context = rep("all", nrow(embedding)),
                                                  condition = rep("all", nrow(embedding)),
                                                  condition_support = matrix(TRUE,
                                                    nrow = length(unique(condition)),
                                                    ncol = length(supported_terminals),
                                                    dimnames = list(unique(condition),
                                                      as.character(supported_terminals))),
                                                  rna_identity = NULL) {
  terminals <- vapply(lineage_paths, tail, integer(1L), n = 1L)
  evidence <- ComputeLineageProgrammeEvidence(embedding, programme_labels,
    cluster_labels == lineage_paths[[1L]][1L], lineage_weights, reference_weights,
    supported_terminals, condition, condition_support,
    query_cells = !cluster_labels %in% terminals, calibration_context = context)
  evidence$programme_labels <- programme_labels
  evidence$rna_identity <- rna_identity
  evidence$lineage_paths <- lineage_paths
  lineage_terminals <- c(terminals, supported_terminals)
  evidence <- ComputeLineageRootProgrammeEntry(evidence, embedding, cluster_labels,
    lineage_paths, lineage_terminals, lineage_weights, adjacency, context)
  ComputeLineageProgrammePhase(evidence, embedding, cluster_labels, lineage_paths,
    lineage_terminals, program_profiles, adjacency, context, rna_identity)
}

ComputeLineageRootProgrammeEntry <- function(evidence, embedding, cluster_labels,
                                               lineage_paths, lineage_terminals,
                                               base_weights, adjacency,
                                               context = rep("all", nrow(embedding))) {
  embedding <- as.matrix(embedding)
  base_weights <- as.matrix(base_weights)
  weights <- as.matrix(evidence$weights)
  paths <- data.table(node = unlist(lineage_paths, use.names = FALSE),
    path = rep(seq_along(lineage_paths), lengths(lineage_paths)),
    step = sequence(lengths(lineage_paths)))
  families <- paths[, .(terminal = tail(node, 1L), child = node[2L]), by = path]
  children <- sort(unique(families$child))
  family <- match(families$child[match(lineage_terminals, families$terminal)], children)
  nbase <- ncol(base_weights)
  base_family <- family[seq_len(nbase)]
  programme_family <- family[seq.int(nbase + 1L, ncol(weights))]
  rows <- evidence$root_rows
  context_levels <- sort(unique(context[rows]))
  context <- match(context, context_levels)
  root_context <- context[rows]
  root_centers <- rowsum(embedding[rows, , drop = FALSE], root_context) /
    tabulate(root_context, length(context_levels))
  grid <- CJ(context = seq_along(context_levels), child = children)
  grid[, grid_id := .I]
  branch_nodes <- merge(paths, paths[step == 2L, .(path, child = node)], by = "path")
  branch_nodes <- branch_nodes[step >= 2L, .(step = min(step)), by = .(child, node)]
  observed <- data.table(context = context, node = cluster_labels)[,
    .(cells = .N), by = .(context, node)]
  entrance <- merge(branch_nodes, observed, by = "node", allow.cartesian = TRUE)
  entrance <- entrance[, .SD[step == min(step)], by = .(context, child)]
  entrance <- merge(entrance, grid, by = c("context", "child"))
  graph <- (adjacency + Matrix::t(adjacency)) / 2
  edges <- as.data.table(Matrix::summary(graph[rows, , drop = FALSE]))
  edges[, `:=`(source_context = context[rows[i]], context = context[j],
    node = cluster_labels[j])]
  edges <- merge(edges[source_context == context],
    entrance[, .(context, node, grid_id)], by = c("context", "node"))
  target_membership <- Matrix::sparseMatrix(i = edges$j, j = edges$grid_id,
    x = edges$x, dims = c(nrow(embedding), nrow(grid)))
  boundary_mass <- Matrix::colSums(target_membership)
  boundary_end <- as.matrix(Matrix::crossprod(target_membership, embedding)) / boundary_mass
  offsets <- boundary_end - root_centers[grid$context, , drop = FALSE]
  common <- rowsum(offsets, grid$context, na.rm = TRUE) /
    as.numeric(rowsum(as.integer(is.finite(rowSums(offsets))), grid$context))
  precision_common <- t(solve(evidence$model$covariance, t(common)))
  shared_projection <- rowSums(offsets * precision_common[grid$context, , drop = FALSE]) /
    rowSums(common * precision_common)[grid$context]
  delta <- offsets - common[grid$context, , drop = FALSE] * shared_projection
  directions <- solve(evidence$model$covariance, t(delta))
  origin <- rowSums(root_centers[grid$context, , drop = FALSE] * t(directions))
  all_scores <- sweep(embedding[rows, , drop = FALSE] %*% directions, 2L, origin, "-")
  selected_columns <- outer((root_context - 1L) * length(children),
    seq_along(children), "+")
  score <- matrix(all_scores[cbind(rep(seq_along(rows), length(children)),
    as.vector(selected_columns))], nrow = length(rows))
  context_mean <- rowsum(score, root_context) / tabulate(root_context, length(context_levels))
  score <- score - context_mean[root_context, , drop = FALSE]
  variance <- rowsum(score^2, root_context) / tabulate(root_context, length(context_levels))
  noise <- sqrt(variance[root_context, , drop = FALSE])
  z <- score / noise
  distance <- matrix(rowSums(delta * t(directions)), nrow = length(context_levels), byrow = TRUE)
  segment <- distance[root_context, , drop = FALSE] / noise
  represented <- !is.na(segment)
  lower <- -z[represented]
  upper <- segment[represented] - z[represented]
  reflect <- lower > 0
  interval_upper <- upper
  interval_lower <- lower
  interval_upper[reflect] <- -lower[reflect]
  interval_lower[reflect] <- -upper[reflect]
  log_upper <- pnorm(interval_upper, log.p = TRUE)
  log_interval <- log_upper + log(-expm1(pnorm(interval_lower, log.p = TRUE) - log_upper))
  log_evidence <- matrix(-Inf, nrow = length(rows), ncol = length(children))
  log_evidence[represented] <- z[represented]^2 / 2 + log(2 * pi) / 2 +
    log_interval - log(segment[represented])
  normalizer <- pmax(0, matrixStats::rowMaxs(log_evidence))
  ancestral <- exp(-normalizer)
  excess <- pmax(exp(log_evidence - normalizer) - ancestral, 0)
  total <- ancestral + rowSums(excess)
  ancestral <- ancestral / total
  activated <- excess / total
  complete <- all(lineage_terminals[seq_len(nbase)] %in% evidence$programme_terminals)
  switch(as.character(complete),
    "TRUE" = {
      family_mass <- base_weights[rows, , drop = FALSE] %*%
        (outer(base_family, seq_along(children), "==") * 1)
      programme_columns <- match(evidence$programme_terminals, lineage_terminals[seq_len(nbase)])
      conditional <- base_weights[rows, programme_columns, drop = FALSE] /
        family_mass[, programme_family, drop = FALSE]
      weights[rows, ] <- cbind(base_weights[rows, , drop = FALSE] * ancestral,
        activated[, programme_family, drop = FALSE] * conditional)
    },
    "FALSE" = {
      original <- weights[rows, , drop = FALSE]
      family_mass <- original %*% (outer(family, seq_along(children), "==") * 1)
      conditional <- original / family_mass[, family, drop = FALSE]
      weights[rows, ] <- original * ancestral + activated[, family, drop = FALSE] * conditional
    })
  expected_mass <- switch(as.character(complete),
    "TRUE" = rowSums(base_weights[rows, , drop = FALSE]),
    "FALSE" = rowSums(evidence$weights[rows, , drop = FALSE])) * ancestral + rowSums(activated)
  stopifnot(all(is.finite(weights)), all(weights >= 0),
    max(abs(rowSums(weights[rows, , drop = FALSE]) - expected_mass)) < 1e-12)
  evidence$weights <- weights
  evidence$entry <- list(children = children, activation = cbind(ancestral, activated),
    standardized = z, segment_length = segment, complete_programme_basis = complete,
    context_levels = context_levels, entrance = entrance)
  evidence
}

ComputeLineageProgrammePhase <- function(evidence, embedding, cluster_labels,
                                           lineage_paths, lineage_terminals,
                                           program_profiles, adjacency,
                                           context = rep("all", nrow(embedding)),
                                           rna_identity = NULL) {
  embedding <- as.matrix(embedding)
  query <- as.matrix(evidence$weights)
  reference <- as.matrix(evidence$reference_weights)
  levels <- unique(lineage_terminals)
  selected <- evidence$programme_terminals
  nbase <- length(levels)
  programme_columns <- seq.int(nbase + 1L, ncol(query))
  stopifnot(identical(as.character(lineage_terminals[programme_columns]),
    as.character(selected)), nrow(query) == nrow(embedding),
    identical(dim(query), dim(reference)))
  paths <- data.table(node = unlist(lineage_paths, use.names = FALSE),
    path = rep(seq_along(lineage_paths), lengths(lineage_paths)),
    step = sequence(lengths(lineage_paths)))
  paths[, terminal := tail(node, 1L), by = path]
  shared <- paths[, .(n_paths = uniqueN(path)), by = node][n_paths > 1L, node]
  forks <- paths[node %in% shared, .(fork = tail(node, 1L),
    fork_step = max(step)), by = path]
  stage <- merge(paths[terminal %in% selected], forks, by = "path")
  setorder(stage, path, step)
  stage[, parent := shift(node, fill = node[1L]), by = path]
  edge_delta <- program_profiles[as.character(stage$node), , drop = FALSE] -
    program_profiles[as.character(stage$parent), , drop = FALSE]
  stage[, net_change := rowSums(edge_delta * abs(edge_delta))]
  stage[, score := cumsum(net_change), by = path]
  stage[, score := score - score[step == fork_step][1L], by = path]
  anchor_mass <- rowsum(reference[, programme_columns, drop = FALSE], cluster_labels)
  stage[, anchor_mass := anchor_mass[cbind(match(as.character(node),
    rownames(anchor_mass)), match(terminal, selected))]]
  chosen <- stage[step > fork_step & anchor_mass > 0][
    order(-score, step, node), .SD[1L], by = terminal]
  chosen <- chosen[match(selected, terminal)]
  stopifnot(!anyNA(chosen$node))
  memberships <- outer(cluster_labels, chosen$node, "==") * 1
  graph <- (adjacency + Matrix::t(adjacency)) / 2
  support <- as.matrix(graph %*% memberships) / Matrix::rowSums(graph)
  reference[, programme_columns] <- reference[, programme_columns, drop = FALSE] * support
  query_mass <- colSums(query)
  reference_mass <- colSums(reference)
  stopifnot(all(reference_mass > 0))
  scale <- query_mass / reference_mass
  membership <- outer(lineage_terminals, levels, "==") * 1
  pooled_query <- query %*% membership
  pooled_reference <- sweep(reference, 2L, scale, "*") %*% membership
  columns <- match(selected, levels)
  model <- evidence$model
  score <- sweep(embedding, 2L, model$origin, "-") %*% model$directions
  context <- match(context, sort(unique(context)))
  root_context <- context[evidence$root_rows]
  context_mean <- rowsum(score[evidence$root_rows, , drop = FALSE], root_context) /
    tabulate(root_context)
  context_distance <- matrix(model$distance, nrow = nrow(context_mean),
    ncol = length(selected), byrow = TRUE) - context_mean
  relative <- score - context_mean[context, , drop = FALSE]
  rank_data <- data.table(row = rep(seq_len(nrow(query)), length(selected)),
    terminal = rep(selected, each = nrow(query)),
    score = as.vector(relative / context_distance[context, , drop = FALSE]),
    mass = as.vector(pooled_query[, columns, drop = FALSE]))
  ties <- rank_data[, .(mass = sum(mass)), by = .(terminal, score)]
  setorder(ties, terminal, score)
  ties[, phase := (cumsum(mass) - mass / 2) / sum(mass), by = terminal]
  rank_data[ties, on = .(terminal, score), phase := i.phase]
  phase <- matrix(rank_data$phase, nrow = nrow(query))
  query_programme <- pooled_query[, columns, drop = FALSE] * phase
  reference_programme <- pooled_reference[, columns, drop = FALSE] * phase
  pooled_query[, columns] <- pooled_query[, columns, drop = FALSE] * (1 - phase)
  pooled_reference[, columns] <- pooled_reference[, columns, drop = FALSE] * (1 - phase)
  weights <- cbind(pooled_query, query_programme)
  reference <- cbind(pooled_reference, reference_programme)
  stopifnot(all(is.finite(weights)), all(is.finite(reference)),
    all(weights >= 0), all(reference >= 0),
    max(abs(rowSums(weights) - rowSums(query))) < 1e-12)
  evidence$weights <- weights
  evidence$reference_weights <- reference
  evidence$lineage_terminals <- c(levels, selected)
  evidence$phase <- phase
  evidence$reference_states <- chosen[, .(terminal, state = node, score)]
  evidence$channel_mass <- data.table(terminal = lineage_terminals,
    query_mass = query_mass, reference_mass = reference_mass, scale = scale)
  switch(as.character(!is.null(rna_identity)), "TRUE" = {
    evidence$identity_readout <- PoolLineageProgrammeIdentities(weights, reference,
      c(levels, selected), selected, rna_identity$mapping)
  })
  evidence
}

PoolLineageProgrammeIdentities <- function(lineage_weights, reference_weights,
                                             lineage_terminals, programme_terminals,
                                             identity_mapping) {
  mapping <- as.data.table(identity_mapping)
  channels <- seq_along(lineage_terminals)
  programme <- channels > length(channels) - length(programme_terminals)
  identities <- mapping$identity[match(lineage_terminals, mapping$state)]
  stopifnot(!anyNA(identities),
    identical(as.character(lineage_terminals[programme]), as.character(programme_terminals)))
  channel_mapping <- data.table(channel = channels, terminal = lineage_terminals,
    programme = programme, identity = identities,
    group_id = fifelse(programme, paste0("programme_", identities), paste0("base_", channels)))
  channel_mapping[, pool := .GRP, by = group_id]
  membership <- Matrix::sparseMatrix(i = channels, j = channel_mapping$pool, x = 1,
    dims = c(length(channels), max(channel_mapping$pool)))
  list(weights = as.matrix(lineage_weights %*% membership),
    reference_weights = as.matrix(reference_weights %*% membership),
    channel_mapping = channel_mapping)
}

RefineRootEntryMixture <- function(embedding, cluster_labels, lineage_paths,
                                   lineage_weights, lineage_terminals,
                                   adjacency, covariance,
                                   context = rep("all", nrow(embedding)),
                                   aggregation = "log",
                                   geometry = c("covariance_contrast", "euclidean_raw")) {
  geometry <- match.arg(geometry)
  embedding <- as.matrix(embedding)
  weights <- as.matrix(lineage_weights)
  paths <- data.table(node = unlist(lineage_paths, use.names = FALSE),
    path = rep(seq_along(lineage_paths), lengths(lineage_paths)),
    step = sequence(lengths(lineage_paths)))
  root_rows <- which(cluster_labels == lineage_paths[[1L]][1L])
  children <- sort(unique(paths[step == 2L, node]))
  context <- as.character(context)
  precision <- switch(geometry,
    covariance_contrast = solve(covariance),
    euclidean_raw = diag(ncol(embedding)))
  branch_nodes <- merge(paths, paths[step == 2L, .(path, child = node)], by = "path")
  branch_nodes <- branch_nodes[step >= 2L, .(step = min(step)), by = .(child, node)]
  observed <- data.table(context = context, node = cluster_labels)[, .(cells = .N),
    by = .(context, node)]
  entrance <- merge(branch_nodes, observed, by = "node", allow.cartesian = TRUE)
  entrance <- entrance[, .SD[step == min(step)], by = .(context, child)]
  graph <- (adjacency + Matrix::t(adjacency)) / 2
  edges <- as.data.table(Matrix::summary(graph[root_rows, , drop = FALSE]))
  edges[, `:=`(source_context = context[root_rows[i]], context = context[j],
    node = cluster_labels[j])]
  edges <- merge(edges[source_context == context], entrance[, .(context, child, node)],
    by = c("context", "node"))
  targets <- edges[, .(mass = sum(x)), by = .(context, child, target = j)]
  roots <- data.table(root_row = root_rows, context = context[root_rows])
  evidence <- roots[, {
    query_rows <- root_row
    current_context <- .BY$context
    origin <- colMeans(embedding[query_rows, , drop = FALSE])
    centered <- sweep(embedding[query_rows, , drop = FALSE], 2L, origin, "-")
    rays <- targets[context == current_context]
    ray_child <- match(rays$child, children)
    rays[, total := sum(mass), by = child]
    membership <- Matrix::sparseMatrix(i = seq_len(nrow(rays)), j = ray_child,
      x = rays$mass / rays$total, dims = c(nrow(rays), length(children)))
    delta <- sweep(embedding[rays$target, , drop = FALSE], 2L, origin, "-")
    branch_delta <- as.matrix(Matrix::crossprod(membership, delta))
    represented <- which(Matrix::colSums(membership) > 0)
    switch(geometry, covariance_contrast = {
      common <- colMeans(branch_delta[represented, , drop = FALSE])
      precision_common <- as.vector(precision %*% common)
      delta <- delta - (as.vector(delta %*% precision_common) /
        sum(common * precision_common)) %o% common
    })
    directions <- precision %*% t(delta)
    score <- centered %*% directions
    noise <- sqrt(colMeans(score^2))
    z <- sweep(score, 2L, noise, "/")
    lengths <- rowSums(delta * t(directions)) / noise
    segment <- matrix(lengths, nrow = length(query_rows), ncol = length(lengths), byrow = TRUE)
    lower <- -z
    upper <- segment - z
    reflect <- lower > 0
    interval_upper <- upper
    interval_lower <- lower
    interval_upper[reflect] <- -lower[reflect]
    interval_lower[reflect] <- -upper[reflect]
    log_upper <- pnorm(interval_upper, log.p = TRUE)
    log_interval <- log_upper + log(-expm1(pnorm(interval_lower, log.p = TRUE) - log_upper))
    log_evidence <- z^2 / 2 + log(2 * pi) / 2 + log_interval - log(segment)
    result <- matrix(0, nrow = length(query_rows), ncol = length(children))
    switch(aggregation,
      log = { result <- as.matrix(log_evidence %*% membership) },
      likelihood = {
        row_max <- matrixStats::rowMaxs(log_evidence)
        mixture <- as.matrix(exp(log_evidence - row_max) %*% membership)
        result[, represented] <- log(mixture[, represented, drop = FALSE]) + row_max
      })
    support <- matrix(rep(seq_along(children) %in% represented, each = length(query_rows)),
      nrow = length(query_rows), ncol = length(children))
    cbind(data.table(row = query_rows),
      setNames(as.data.table(result), paste0("evidence_", children)),
      setNames(as.data.table(support), paste0("observed_", children)))
  }, by = context]
  setorder(evidence, row)
  families <- merge(paths[, .(terminal = tail(node, 1L)), by = path],
    paths[step == 2L, .(path, child = node)], by = "path")
  family <- match(families$child[match(lineage_terminals, families$terminal)], children)
  values <- as.matrix(evidence[, paste0("evidence_", children), with = FALSE])
  support <- as.matrix(evidence[, paste0("observed_", children), with = FALSE])
  factors <- exp(values - matrixStats::rowMaxs(values))
  original <- weights[evidence$row, , drop = FALSE]
  active <- original * support[, family, drop = FALSE]
  revised <- active * factors[, family, drop = FALSE]
  weights[evidence$row, ] <- original * (!support[, family, drop = FALSE]) +
    revised / rowSums(revised) * rowSums(active)
  list(weights = weights, root_rows = evidence$row, children = children,
    log_evidence = values, observed = support, context = evidence$context)
}

ReconnectConservedLineageTips <- function(program_profiles, induced_program,
                                         cluster_labels, pseudotime,
                                         lineage_paths, lineage_weights,
                                         direction = "forward",
                                         prediction = max.col(lineage_weights,
                                           ties.method = "first"),
                                         identity_terminals = integer()) {
  cl_map <- .compact_cluster_labels(cluster_labels)
  cluster <- cl_map$labels
  paths <- lapply(lineage_paths, match, table = cl_map$levels)
  root <- paths[[1L]][1L]
  phase <- as.numeric(rowsum(pseudotime, cluster) / tabulate(cluster))
  profile <- as.matrix(program_profiles)
  program <- as.matrix(induced_program)
  raw_program <- pmax(sweep(profile, 2L, profile[root, ], "-"), 0)
  switch(direction, identity = {
    centered_phase <- phase - mean(phase)
    centered_profile <- sweep(profile, 2L, colMeans(profile), "-")
    slope <- as.numeric(crossprod(centered_phase, profile) / sum(centered_phase^2))
    residual <- centered_profile - outer(centered_phase, slope)
    program <- cbind(pmax(residual, 0), pmax(-residual, 0))
    raw_program <- cbind(pmax(centered_profile, 0), pmax(-centered_profile, 0))
  })
  nodes <- as.data.table(stack(setNames(paths, seq_along(paths))))
  setnames(nodes, c("values", "ind"), c("node", "lineage"))
  nodes[, from := shift(node), by = lineage]
  edges <- unique(nodes[!is.na(from), .(from, to = node)])
  source_leaves <- vapply(paths, tail, integer(1L), n = 1L)
  tips <- edges[to %in% source_leaves, .(parent = from, state = to)]
  pairs <- switch(direction,
    forward = , plateau = , identity = {
      ancestry <- as.data.table(stack(setNames(paths, source_leaves)))
      setnames(ancestry, c("values", "ind"), c("parent", "state_b"))
      ancestry[, state_b := as.integer(as.character(state_b))]
      merge(tips[, .(parent, state_a = state)], ancestry,
        by = "parent", allow.cartesian = TRUE)
    },
    backward = merge(tips, tips, by = "parent", suffixes = c("_a", "_b"),
      allow.cartesian = TRUE))
  pairs <- pairs[state_a != state_b]
  setnames(pairs, c("state_a", "state_b"), c("from", "to"))
  pairs[, `:=`(phase_from = phase[from], phase_parent = phase[parent],
    phase_to = phase[to],
    program_loss = rowSums(pmax(program[from, , drop = FALSE] -
      program[to, , drop = FALSE], 0)^2),
    reverse_loss = rowSums(pmax(program[to, , drop = FALSE] -
      program[from, , drop = FALSE], 0)^2),
    retained_program = rowSums(pmin(program[from, , drop = FALSE],
      program[to, , drop = FALSE])^2),
    program_distance2 = rowSums((program[from, , drop = FALSE] -
      program[to, , drop = FALSE])^2),
    raw_loss = rowSums(pmax(raw_program[from, , drop = FALSE] -
      raw_program[to, , drop = FALSE], 0)^2),
    raw_reverse_loss = rowSums(pmax(raw_program[to, , drop = FALSE] -
      raw_program[from, , drop = FALSE], 0)^2),
    raw_retained = rowSums(pmin(raw_program[from, , drop = FALSE],
      raw_program[to, , drop = FALSE])^2))]
  pairs[, phase_eligible := switch(direction,
    forward = phase_parent < phase_from & phase_from < phase_to,
    identity = to %in% match(identity_terminals, cl_map$levels) &
      !from %in% match(identity_terminals, cl_map$levels),
    plateau = phase_from <= phase_parent & phase_to <= phase_parent &
      (phase_from > phase_to | (phase_from == phase_to & from > to)),
    backward = phase_from <= phase_parent & phase_to > phase_parent)]
  pairs[, distinctive_identity_lost := FALSE]
  switch(direction, forward = {
    acquired <- pmax(program[pairs$from, , drop = FALSE] -
      program[pairs$parent, , drop = FALSE], 0)
    target_acquired <- pmax(program[pairs$to, , drop = FALSE] -
      program[pairs$parent, , drop = FALSE], 0)
    raw_acquired <- pmax(profile[pairs$from, , drop = FALSE] -
      profile[pairs$parent, , drop = FALSE], 0)
    raw_target_acquired <- pmax(profile[pairs$to, , drop = FALSE] -
      profile[pairs$parent, , drop = FALSE], 0)
    target_parent <- edges$from[match(pairs$to, edges$to)]
    target_persistent <- pmin(target_acquired,
      pmax(program[target_parent, , drop = FALSE] - program[pairs$parent, , drop = FALSE], 0))
    raw_target_persistent <- pmin(raw_target_acquired,
      pmax(profile[target_parent, , drop = FALSE] - profile[pairs$parent, , drop = FALSE], 0))
    pairs[, `:=`(acquired_energy = rowSums(acquired^2),
      inherited_energy = rowSums(pmin(program[from, , drop = FALSE],
        program[parent, , drop = FALSE])^2),
      acquired_loss = rowSums(pmax(acquired - target_acquired, 0)^2),
      acquired_retained = rowSums(pmin(acquired, target_acquired)^2),
      raw_acquired_loss = rowSums(pmax(raw_acquired - raw_target_acquired, 0)^2),
      raw_acquired_retained = rowSums(pmin(raw_acquired, raw_target_acquired)^2))]
    pairs[, recovered_identity_lost := to %in% match(identity_terminals, cl_map$levels) &
      rowSums(target_persistent^2) >= rowSums((target_acquired - target_persistent)^2) &
      rowSums(pmax(target_persistent - acquired, 0)^2) >
        rowSums(pmin(target_persistent, acquired)^2) &
      rowSums(pmax(raw_target_persistent - raw_acquired, 0)^2) >
        rowSums(pmin(raw_target_persistent, raw_acquired)^2)]
    pairs[, distinctive_identity_lost := (acquired_energy > inherited_energy | recovered_identity_lost) &
      acquired_loss > acquired_retained & raw_acquired_loss > raw_acquired_retained]
  })
  candidates <- pairs[phase_eligible & retained_program > 0 &
    !distinctive_identity_lost &
    program_loss <= retained_program & raw_loss <= raw_retained &
    (direction == "forward" | (reverse_loss <= retained_program &
      raw_reverse_loss <= raw_retained))][
      order(program_distance2, to), .SD[1L], by = from]
  insertion <- data.table()
  repaired_paths <- paths
  weights <- lineage_weights
  mapped_leaves <- source_leaves
  switch(direction, plateau = , identity = {
    target_nodes <- as.data.table(stack(setNames(paths, source_leaves)))
    setnames(target_nodes, c("values", "ind"), c("node", "target"))
    target_nodes[, target := as.integer(as.character(target))]
    target_nodes[, step := seq_len(.N), by = target]
    arm <- merge(candidates[, .(orphan = from, parent, target = to)], target_nodes, by = "target")
    arm[, parent_step := step[match(parent[1L], node)], by = orphan]
    arm <- arm[step > parent_step]
    child_counts <- edges[, .N, by = from]
    arm[, n_children := fcoalesce(child_counts$N[match(node, child_counts$from)], 0L)]
    admissible <- arm[, .(unforked = all(n_children <= 1L)), by = orphan][unforked == TRUE, orphan]
    arm <- arm[orphan %in% admissible]
    arm[, `:=`(
      program_loss = rowSums(pmax(program[node, , drop = FALSE] - program[target, , drop = FALSE], 0)^2),
      reverse_loss = rowSums(pmax(program[target, , drop = FALSE] - program[node, , drop = FALSE], 0)^2),
      retained_program = rowSums(pmin(program[node, , drop = FALSE], program[target, , drop = FALSE])^2),
      raw_loss = rowSums(pmax(raw_program[node, , drop = FALSE] - raw_program[target, , drop = FALSE], 0)^2),
      raw_reverse_loss = rowSums(pmax(raw_program[target, , drop = FALSE] - raw_program[node, , drop = FALSE], 0)^2),
      raw_retained = rowSums(pmin(raw_program[node, , drop = FALSE], raw_program[target, , drop = FALSE])^2))]
    admissible <- arm[, .(conserved = all(program_loss <= retained_program &
      reverse_loss <= retained_program & raw_loss <= raw_retained & raw_reverse_loss <= raw_retained)),
      by = orphan][conserved == TRUE, orphan]
    candidates <- candidates[from %in% admissible][order(parent, program_distance2, from),
      .SD[to == to[1L]], by = parent]
    insertion <- unique(rbind(arm[orphan %in% admissible & target %in% candidates$to,
      .(parent, target, node)], candidates[, .(parent, target = to, node = from)]))
    insertion[, node_phase := phase[node]]
    setorder(insertion, parent, node_phase, node)
    insertion[, from := shift(node, fill = parent[1L]), by = parent]
  })
  switch(as.character(nrow(candidates) > 0L & !direction %in% c("plateau", "identity")), "TRUE" = {
    target_paths <- paths[match(candidates$to, source_leaves)]
    target_nodes <- as.data.table(stack(setNames(target_paths, candidates$from)))
    setnames(target_nodes, c("values", "ind"), c("node", "orphan"))
    target_nodes[, orphan := as.integer(as.character(orphan))]
    target_nodes[, step := seq_len(.N), by = orphan]
    target_nodes[, from := shift(node), by = orphan]
    insertion <- switch(direction,
      forward = target_nodes[!is.na(from)],
      backward = target_nodes[!is.na(from) & phase[from] <= phase[orphan] &
        phase[node] > phase[orphan]][order(-step), .SD[1L], by = orphan])
    insertion[, `:=`(
      phase_gap = pmax(phase[from] - phase[orphan], phase[orphan] - phase[node], 0),
      incoming_loss = rowSums(pmax(raw_program[from, , drop = FALSE] -
        raw_program[orphan, , drop = FALSE], 0)^2),
      incoming_retained = rowSums(pmin(raw_program[from, , drop = FALSE],
        raw_program[orphan, , drop = FALSE])^2),
      outgoing_loss = rowSums(pmax(raw_program[orphan, , drop = FALSE] -
        raw_program[node, , drop = FALSE], 0)^2),
      outgoing_retained = rowSums(pmin(raw_program[orphan, , drop = FALSE],
        raw_program[node, , drop = FALSE])^2),
      incoming_program_loss = rowSums(pmax(program[from, , drop = FALSE] -
        program[orphan, , drop = FALSE], 0)^2),
      incoming_program_retained = rowSums(pmin(program[from, , drop = FALSE],
        program[orphan, , drop = FALSE])^2),
      outgoing_program_loss = rowSums(pmax(program[orphan, , drop = FALSE] -
        program[node, , drop = FALSE], 0)^2),
      outgoing_program_retained = rowSums(pmin(program[orphan, , drop = FALSE],
        program[node, , drop = FALSE])^2))]
    insertion <- insertion[incoming_loss <= incoming_retained &
      outgoing_loss <= outgoing_retained & (direction != "forward" |
        (incoming_program_loss <= incoming_program_retained &
          outgoing_program_loss <= outgoing_program_retained))][
            order(phase_gap, -step), .SD[1L], by = orphan]
    candidates <- candidates[from %in% insertion$orphan]
  })
  switch(as.character(nrow(insertion) > 0L), "TRUE" = {
    switch(direction, plateau = , identity = {
      edges[insertion, on = c(to = "node"), from := i.from]
    }, {
    setnames(insertion, c("from", "node"), c("insert_from", "insert_to"))
    insertion[, `:=`(source_insert_from = insert_from, source_insert_to = insert_to)]
    dependency <- matrix(0, nrow(insertion), nrow(insertion))
    linked <- which(insertion$insert_to %in% insertion$orphan)
    dependency[cbind(linked, match(insertion$insert_to[linked], insertion$orphan))] <- 1
    reach <- solve(diag(nrow(insertion)) - dependency)
    anchors <- which(!insertion$insert_to %in% insertion$orphan)
    anchor <- anchors[max.col(reach[, anchors, drop = FALSE], ties.method = "first")]
    insertion[, `:=`(insert_from = insert_from[anchor], insert_to = insert_to[anchor],
      orphan_phase = phase[orphan])]
    setorder(insertion, insert_from, insert_to, orphan_phase, orphan)
    insertion[, new_from := shift(orphan, fill = insert_from[1L]),
      by = .(insert_from, insert_to)]
    edges[insertion, on = c(to = "orphan"), from := i.new_from]
    following <- insertion[, .(from = last(orphan), to = insert_to[1L]),
      by = .(insert_from, insert_to)]
    edges[following, on = "to", from := i.from]
    })
    repaired_paths <- rs_extract_lineages_from_edges(as.integer(edges$from),
      as.integer(edges$to), as.integer(nrow(program)), as.integer(root))
    leaves <- vapply(repaired_paths, tail, integer(1L), n = 1L)
    dependency <- matrix(0, length(source_leaves), length(source_leaves))
    dependency[cbind(match(candidates$from, source_leaves),
      match(candidates$to, source_leaves))] <- 1
    reach <- solve(diag(length(source_leaves)) - dependency)
    surviving <- which(!source_leaves %in% candidates$from)
    mapped_leaves <- source_leaves[surviving[max.col(
      reach[, surviving, drop = FALSE], ties.method = "first")]]
    switch(direction, plateau = , identity = {
      plateau_leaves <- insertion[, .(source_leaf = target[1L], leaf = last(node)), by = parent]
      matched <- match(mapped_leaves, plateau_leaves$source_leaf)
      mapped_leaves[!is.na(matched)] <- plateau_leaves$leaf[matched[!is.na(matched)]]
    })
    wide <- as.data.table(lineage_weights)
    setnames(wide, as.character(seq_len(ncol(wide))))
    wide[, cell := seq_len(.N)]
    pooled <- melt(wide, id.vars = "cell", variable.name = "lineage", value.name = "weight")
    pooled[, lineage := match(mapped_leaves, leaves)[as.integer(as.character(lineage))]]
    pooled <- pooled[, .(weight = max(weight)), by = .(cell, lineage)]
    pooled[, weight := weight / sum(weight), by = cell]
    merged <- dcast(pooled, cell ~ lineage, value.var = "weight")
    setorder(merged, cell)
    weights <- as.matrix(merged[, as.character(seq_along(repaired_paths)), with = FALSE])
    prediction <- match(mapped_leaves[prediction], leaves)
    stopifnot(identical(sort(leaves), sort(unique(mapped_leaves))),
      identical(sort(unique(unlist(repaired_paths))), seq_len(nrow(program))))
  })
  dimnames(weights) <- list(rownames(lineage_weights), NULL)
  list(lineage_paths = relist(cl_map$levels[unlist(repaired_paths)], repaired_paths),
    lineage_weights = weights, prediction = prediction, candidates = candidates,
    insertion = insertion, repaired_edges = edges, cluster_levels = cl_map$levels,
    source_leaves = source_leaves, mapped_leaves = mapped_leaves, direction = direction)
}

RecoverPersistentLineageIdentity <- function(program_profiles, cluster_labels,
                                            pseudotime, lineage_paths,
                                            lineage_weights, embedding,
                                            ancestor_candidates, sigma = 1,
                                            prediction = max.col(lineage_weights,
                                              ties.method = "first")) {
  cl_map <- .compact_cluster_labels(cluster_labels)
  cluster <- cl_map$labels
  paths <- lapply(lineage_paths, match, table = cl_map$levels)
  root <- paths[[1L]][1L]
  profile <- as.matrix(program_profiles)
  phase <- as.numeric(rowsum(pseudotime, cluster) / tabulate(cluster))
  centered_phase <- phase - mean(phase)
  slope <- as.numeric(crossprod(centered_phase, profile) / sum(centered_phase^2))
  residual <- sweep(profile, 2L, colMeans(profile), "-") - outer(centered_phase, slope)
  program <- cbind(pmax(residual, 0), pmax(-residual, 0))
  nodes <- as.data.table(stack(setNames(paths, seq_along(paths))))
  setnames(nodes, c("values", "ind"), c("node", "lineage"))
  nodes[, step := seq_len(.N), by = lineage]
  nodes[, from := shift(node), by = lineage]
  edges <- unique(nodes[!is.na(from), .(from, to = node)])
  children <- edges[, .(n_children = .N), by = from]
  forks <- c(root, children[n_children > 1L, from])
  ancestry <- merge(nodes[, .(lineage, node, node_step = step)],
    nodes[, .(lineage, anchor = node, anchor_step = step)],
    by = "lineage", allow.cartesian = TRUE)
  ancestry <- unique(ancestry[anchor_step < node_step & anchor %in% forks,
    .(node, anchor, anchor_step)])
  anchors <- ancestry[order(node, -anchor_step), .SD[1L], by = node]
  baseline <- rep(root, nrow(program))
  baseline[anchors$node] <- anchors$anchor
  parent <- rep(root, nrow(program))
  parent[edges$to] <- edges$from
  boundaries <- merge(nodes[, .(lineage, node, node_step = step, anchor = baseline[node])],
    nodes[, .(lineage, anchor = node, anchor_step = step)], by = c("lineage", "anchor"))
  identity_ancestry <- merge(boundaries,
    nodes[, .(lineage, ancestor = node, ancestor_step = step)],
    by = "lineage", allow.cartesian = TRUE)
  identity_ancestry <- unique(identity_ancestry[
    ancestor_step >= anchor_step & ancestor_step <= node_step, .(node, ancestor)])
  grid <- data.table(node = rep(identity_ancestry$node, each = ncol(program)),
    coordinate = rep(seq_len(ncol(program)), times = nrow(identity_ancestry)),
    value = as.numeric(t(program[identity_ancestry$ancestor, , drop = FALSE])))
  grid <- grid[, .(value = min(value)), by = .(node, coordinate)]
  wide <- dcast(grid, node ~ coordinate, value.var = "value")
  setorder(wide, node)
  conserved <- as.matrix(wide[, as.character(seq_len(ncol(program))), with = FALSE])
  diagnostics <- copy(edges)
  diagnostics[, n_children := children$n_children[match(from, children$from)]]
  diagnostics[, `:=`(parent = parent[from], anchor = baseline[from],
    conserved_energy = rowSums(conserved[from, , drop = FALSE]^2),
    identity_change = rowSums((program[from, , drop = FALSE] -
      conserved[from, , drop = FALSE])^2),
    program_loss = rowSums(pmax(conserved[from, , drop = FALSE] -
      program[to, , drop = FALSE], 0)^2),
    retained_program = rowSums(pmin(conserved[from, , drop = FALSE],
      program[to, , drop = FALSE])^2))]
  leaves <- vapply(paths, tail, integer(1L), n = 1L)
  leaf_matches <- CJ(from = unique(edges$from), to = leaves)[
    from != to & phase[from] < phase[to]]
  leaf_matches[, `:=`(
    program_loss = rowSums(pmax(conserved[from, , drop = FALSE] -
      program[to, , drop = FALSE], 0)^2),
    retained_program = rowSums(pmin(conserved[from, , drop = FALSE],
      program[to, , drop = FALSE])^2))]
  matching <- unique(leaf_matches[program_loss <= retained_program, from])
  lost <- diagnostics[, .(all_lost = all(program_loss > retained_program)),
    by = from][all_lost == TRUE, from]
  candidates <- diagnostics[from != root & n_children > 1L & conserved_energy > 0 &
    identity_change <= conserved_energy & from %in% lost & !from %in% matching &
    phase[from] < phase[to]]
  candidates[, identity_kind := "persistent_identity"]
  acquired <- pmax(sweep(profile, 2L, profile[root, ], "-"), 0)
  branch_program <- pmax(acquired - acquired[baseline, , drop = FALSE], 0)
  branch_parent <- pmax(acquired[parent, , drop = FALSE] -
    acquired[baseline, , drop = FALSE], 0)
  branch_conserved <- pmin(branch_program, branch_parent)
  branch_diagnostics <- copy(edges)
  branch_diagnostics[, n_children := children$n_children[match(from, children$from)]]
  branch_diagnostics[, `:=`(parent = parent[from], anchor = baseline[from],
    conserved_energy = rowSums(branch_conserved[from, , drop = FALSE]^2),
    identity_change = rowSums((branch_program[from, , drop = FALSE] -
      branch_conserved[from, , drop = FALSE])^2),
    program_loss = rowSums(pmax(branch_conserved[from, , drop = FALSE] -
      pmax(acquired[to, , drop = FALSE] - acquired[baseline[from], , drop = FALSE], 0), 0)^2),
    retained_program = rowSums(pmin(branch_conserved[from, , drop = FALSE],
      pmax(acquired[to, , drop = FALSE] - acquired[baseline[from], , drop = FALSE], 0))^2))]
  branch_leaf_matches <- CJ(from = unique(edges$from), to = leaves)[
    from != to & phase[from] < phase[to]]
  branch_leaf_matches[, `:=`(
    program_loss = rowSums(pmax(branch_conserved[from, , drop = FALSE] -
      pmax(acquired[to, , drop = FALSE] - acquired[baseline[from], , drop = FALSE], 0), 0)^2),
    retained_program = rowSums(pmin(branch_conserved[from, , drop = FALSE],
      pmax(acquired[to, , drop = FALSE] - acquired[baseline[from], , drop = FALSE], 0))^2))]
  branch_matching <- unique(branch_leaf_matches[program_loss <= retained_program, from])
  branch_candidates <- branch_diagnostics[n_children == 1L & parent != anchor &
    conserved_energy > 0 & identity_change <= conserved_energy &
    program_loss > retained_program & !from %in% branch_matching]
  branch_candidates[, identity_kind := "acquired_branch"]
  candidates <- rbindlist(list(candidates, branch_candidates), fill = TRUE)
  ancestry <- merge(nodes[, .(lineage, to = node, to_step = step)],
    nodes[, .(lineage, from = node, from_step = step)],
    by = "lineage", allow.cartesian = TRUE)
  ancestry <- unique(ancestry[to_step - from_step > 1L,
    .(from, to, gap = to_step - from_step)])
  ancestor_pairs <- copy(as.data.table(ancestor_candidates))
  ancestor_pairs[, `:=`(from = match(from, cl_map$levels), to = match(to, cl_map$levels))]
  alternatives <- merge(ancestor_pairs, ancestry, by = c("from", "to"))
  alternatives <- merge(alternatives, candidates[, .(to, old_from = from)], by = "to")
  repairs <- alternatives[phase_eligible & program_loss <= retained_program][
    order(gap, distance2, from), .SD[1L], by = to]
  repaired_paths <- paths
  weights <- lineage_weights
  switch(as.character(nrow(repairs) > 0L), "TRUE" = {
    edges[repairs, on = "to", from := i.from]
    repaired_paths <- rs_extract_lineages_from_edges(as.integer(edges$from),
      as.integer(edges$to), as.integer(nrow(program)), as.integer(root))
    flat <- rs_compute_lineage_weights_trunk_aware_from_edges(
      as.numeric(t(embedding)), as.integer(cluster), as.integer(nrow(embedding)),
      as.integer(ncol(embedding)), as.integer(nrow(program)),
      as.integer(edges$from), as.integer(edges$to), as.integer(root), as.numeric(sigma))
    weights <- matrix(flat, nrow = nrow(embedding), ncol = length(repaired_paths), byrow = TRUE)
    prediction <- max.col(weights, ties.method = "first")
  })
  dimnames(weights) <- list(rownames(lineage_weights), NULL)
  repaired_leaves <- vapply(repaired_paths, tail, integer(1L), n = 1L)
  new_terminals <- setdiff(repaired_leaves, leaves)
  membership <- data.table()
  switch(as.character(length(new_terminals) > 0L), "TRUE" = {
    membership <- identity_ancestry[node %in% new_terminals,
      .(state = ancestor, terminal = node)]
    membership[, `:=`(
      program_loss = rowSums(pmax(program[state, , drop = FALSE] -
        conserved[terminal, , drop = FALSE], 0)^2),
      retained_program = rowSums(pmin(program[state, , drop = FALSE],
        conserved[terminal, , drop = FALSE])^2))]
    membership <- membership[program_loss <= retained_program][
      order(program_loss, terminal), .SD[1L], by = state]
    membership[, prediction := match(terminal, repaired_leaves)]
    state_index <- match(cluster, membership$state)
    prediction[!is.na(state_index)] <- membership$prediction[state_index[!is.na(state_index)]]
  })
  continuity <- list(candidates = data.table(), insertion = data.table())
  switch(as.character(length(new_terminals) > 0L), "TRUE" = {
    continuity <- ReconnectConservedLineageTips(profile, program, cluster_labels,
      pseudotime, relist(cl_map$levels[unlist(repaired_paths)], repaired_paths),
      weights, direction = "identity", prediction = prediction,
      identity_terminals = cl_map$levels[new_terminals])
    repaired_paths <- lapply(continuity$lineage_paths, match, table = cl_map$levels)
    edges <- continuity$repaired_edges
    weights <- continuity$lineage_weights
    prediction <- continuity$prediction
    repaired_leaves <- vapply(repaired_paths, tail, integer(1L), n = 1L)
    membership[, terminal := continuity$mapped_leaves[
      match(terminal, continuity$source_leaves)]]
    membership[, prediction := match(terminal, repaired_leaves)]
    new_terminals <- setdiff(repaired_leaves, leaves)
  })
  list(lineage_paths = relist(cl_map$levels[unlist(repaired_paths)], repaired_paths),
    lineage_weights = weights, prediction = prediction, candidates = candidates,
    repairs = repairs, repaired_edges = edges, cluster_levels = cl_map$levels,
    new_terminals = cl_map$levels[new_terminals], identity_membership = membership,
    edge_diagnostics = diagnostics, conserved_program = conserved,
    branch_diagnostics = branch_diagnostics,
    identity_continuation = continuity$candidates,
    identity_insertion = continuity$insertion)
}

ResolvePersistentLineageMembership <- function(program_profiles, cluster_labels,
                                               pseudotime, lineage_paths,
                                               lineage_weights,
                                               prediction = max.col(lineage_weights,
                                                 ties.method = "first"),
                                               cell_counts = NULL) {
  cl_map <- .compact_cluster_labels(cluster_labels)
  cluster <- cl_map$labels
  paths <- lapply(lineage_paths, match, table = cl_map$levels)
  root <- paths[[1L]][1L]
  profile <- as.matrix(program_profiles)
  inherited_prediction <- prediction
  phase <- as.numeric(rowsum(pseudotime, cluster) / tabulate(cluster))
  centered_phase <- phase - mean(phase)
  centered_profile <- sweep(profile, 2L, colMeans(profile), "-")
  slope <- as.numeric(crossprod(centered_phase, profile) / sum(centered_phase^2))
  residual <- centered_profile - outer(centered_phase, slope)
  program <- cbind(pmax(residual, 0), pmax(-residual, 0))
  raw_program <- cbind(pmax(centered_profile, 0), pmax(-centered_profile, 0))
  nodes <- as.data.table(stack(setNames(paths, seq_along(paths))))
  setnames(nodes, c("values", "ind"), c("node", "lineage"))
  nodes[, step := seq_len(.N), by = lineage]
  nodes[, from := shift(node), by = lineage]
  edges <- unique(nodes[!is.na(from), .(from, to = node)])
  children <- edges[, .(n_children = .N), by = from]
  forks <- c(root, children[n_children > 1L, from])
  leaves <- vapply(paths, tail, integer(1L), n = 1L)
  endpoints <- nodes[node %in% leaves, .(lineage, terminal = node, terminal_step = step)]
  anchors <- merge(endpoints,
    nodes[node %in% forks, .(lineage, anchor = node, anchor_step = step)], by = "lineage")
  anchors <- anchors[anchor_step < terminal_step][
    order(terminal, -anchor_step), .SD[1L], by = terminal]
  ancestry <- merge(anchors, nodes[, .(lineage, state = node, state_step = step)], by = "lineage")
  ancestry <- ancestry[state_step >= anchor_step & state_step <= terminal_step]
  grid <- data.table(terminal = rep(ancestry$terminal, each = ncol(program)),
    coordinate = rep(seq_len(ncol(program)), times = nrow(ancestry)),
    value = as.numeric(t(program[ancestry$state, , drop = FALSE])),
    raw_value = as.numeric(t(raw_program[ancestry$state, , drop = FALSE])))
  conserved <- grid[, .(
    value = ifelse(min(value) > 0, .N / sum(1 / value), 0),
    raw_value = ifelse(min(raw_value) > 0, .N / sum(1 / raw_value), 0)),
    by = .(terminal, coordinate)]
  wide <- dcast(conserved, terminal ~ coordinate, value.var = "value")
  raw_wide <- dcast(conserved, terminal ~ coordinate, value.var = "raw_value")
  template <- as.matrix(wide[, as.character(seq_len(ncol(program))), with = FALSE])
  raw_template <- as.matrix(raw_wide[, as.character(seq_len(ncol(program))), with = FALSE])
  terminal_index <- match(anchors$terminal, wide$terminal)
  diagnostics <- anchors[, .(terminal, anchor)]
  diagnostics[, `:=`(
    conserved_energy = rowSums(template[terminal_index, , drop = FALSE]^2),
    origin_change = rowSums((program[anchor, , drop = FALSE] -
      template[terminal_index, , drop = FALSE])^2),
    terminal_change = rowSums((program[terminal, , drop = FALSE] -
      template[terminal_index, , drop = FALSE])^2),
    raw_conserved_energy = rowSums(raw_template[terminal_index, , drop = FALSE]^2),
    raw_origin_change = rowSums((raw_program[anchor, , drop = FALSE] -
      raw_template[terminal_index, , drop = FALSE])^2),
    raw_terminal_change = rowSums((raw_program[terminal, , drop = FALSE] -
      raw_template[terminal_index, , drop = FALSE])^2))]
  eligible <- diagnostics[conserved_energy > 0 & origin_change <= conserved_energy &
    terminal_change <= conserved_energy & raw_conserved_energy > 0 &
    raw_origin_change <= raw_conserved_energy & raw_terminal_change <= raw_conserved_energy,
    terminal]
  descendants <- merge(nodes[, .(lineage, anchor = node, anchor_step = step)],
    nodes[, .(lineage, state = node, state_step = step)], by = "lineage", allow.cartesian = TRUE)
  descendants <- unique(descendants[state_step >= anchor_step, .(anchor, state)])
  membership <- merge(diagnostics[terminal %in% eligible, .(anchor, terminal)],
    descendants, by = "anchor")
  template_index <- match(membership$terminal, wide$terminal)
  membership[, `:=`(
    program_loss = rowSums(pmax(program[state, , drop = FALSE] -
      template[template_index, , drop = FALSE], 0)^2),
    reverse_loss = rowSums(pmax(template[template_index, , drop = FALSE] -
      program[state, , drop = FALSE], 0)^2),
    retained_program = rowSums(pmin(program[state, , drop = FALSE],
      template[template_index, , drop = FALSE])^2),
    raw_loss = rowSums(pmax(raw_program[state, , drop = FALSE] -
      raw_template[template_index, , drop = FALSE], 0)^2),
    raw_reverse_loss = rowSums(pmax(raw_template[template_index, , drop = FALSE] -
      raw_program[state, , drop = FALSE], 0)^2),
    raw_retained = rowSums(pmin(raw_program[state, , drop = FALSE],
      raw_template[template_index, , drop = FALSE])^2))]
  membership <- membership[program_loss <= retained_program & reverse_loss <= retained_program &
    raw_loss <= raw_retained & raw_reverse_loss <= raw_retained][
      order(program_loss, terminal), .SD[1L], by = state]
  membership[, prediction := match(terminal, leaves)]
  state_index <- match(cluster, membership$state)
  prediction[!is.na(state_index)] <- membership$prediction[state_index[!is.na(state_index)]]
  targets <- diagnostics[terminal %in% eligible][
    order(anchor, -phase[terminal]), .SD[1L], by = anchor]$terminal
  continuity <- ReconnectConservedLineageTips(profile, program, cluster_labels, pseudotime,
    lineage_paths, lineage_weights, direction = "identity", prediction = prediction,
    identity_terminals = cl_map$levels[targets])
  final_leaves <- match(vapply(continuity$lineage_paths, tail, integer(1L), n = 1L), cl_map$levels)
  membership[, terminal := continuity$mapped_leaves[match(terminal, continuity$source_leaves)]]
  membership[, prediction := match(terminal, final_leaves)]
  cell_identity <- data.table()
  terminal_overlap <- data.table()
  switch(as.character(!is.null(cell_counts)), "TRUE" = {
    template_lineage <- match(continuity$mapped_leaves[
      match(wide$terminal, continuity$source_leaves)], final_leaves)
    terminal_overlap <- CJ(template_terminal = wide$terminal, terminal = final_leaves)
    terminal_overlap[, lineage := template_lineage[match(template_terminal, wide$terminal)]]
    terminal_overlap <- terminal_overlap[terminal != final_leaves[lineage]]
    template_index <- match(terminal_overlap$template_terminal, wide$terminal)
    terminal_overlap[, `:=`(
      program_loss = rowSums(pmax(template[template_index, , drop = FALSE] -
        program[terminal, , drop = FALSE], 0)^2),
      retained_program = rowSums(pmin(template[template_index, , drop = FALSE],
        program[terminal, , drop = FALSE])^2),
      raw_loss = rowSums(pmax(raw_template[template_index, , drop = FALSE] -
        raw_program[terminal, , drop = FALSE], 0)^2),
      raw_retained = rowSums(pmin(raw_template[template_index, , drop = FALSE],
        raw_program[terminal, , drop = FALSE])^2))]
    shared <- unique(terminal_overlap[program_loss <= retained_program &
      raw_loss <= raw_retained, lineage])
    membership <- membership[!prediction %in% shared]
    counts <- as(cell_counts, "CsparseMatrix")
    stopifnot(nrow(counts) == ncol(profile), ncol(counts) == length(cluster))
    expression <- counts %*% Matrix::Diagonal(x = 1 / Matrix::colSums(counts))
    expression@x <- sqrt(expression@x)
    cell_template <- expm1(profile[final_leaves, , drop = FALSE])
    cell_template <- sqrt(cell_template / rowSums(cell_template))
    similarity <- as.matrix(Matrix::t(expression) %*% t(cell_template))
    closest <- max.col(similarity, ties.method = "first")
    inherited_prediction <- match(continuity$mapped_leaves[inherited_prediction], final_leaves)
    state_index <- match(cluster, membership$state)
    identity_prediction <- membership$prediction[state_index]
    assigned <- !is.na(state_index) & closest == identity_prediction
    prediction <- inherited_prediction
    prediction[assigned] <- identity_prediction[assigned]
    cell_identity <- data.table(cell = seq_along(cluster), cluster = cluster,
      inherited_prediction = inherited_prediction,
      state_identity_prediction = continuity$prediction,
      closest_terminal = final_leaves[closest], identity_prediction = identity_prediction,
      assigned = assigned, prediction = prediction)
    continuity$prediction <- prediction
  })
  list(lineage_paths = continuity$lineage_paths, lineage_weights = continuity$lineage_weights,
    prediction = continuity$prediction, candidates = diagnostics,
    identity_membership = membership, identity_targets = cl_map$levels[targets],
    continuity_candidates = continuity$candidates, insertion = continuity$insertion,
    repaired_edges = continuity$repaired_edges, cluster_levels = cl_map$levels,
    conserved_program = template, raw_conserved_program = raw_template,
    template_terminals = wide$terminal, cell_identity = cell_identity,
    terminal_overlap = terminal_overlap)
}

.repair_lineage_program_ancestry <- function(profile, program, persistent, phase,
                                             root, cluster, pseudotime, paths,
                                             pairs, cell_expression,
                                             shared_slope, mode, clock_mask = NULL) {
  nodes <- as.data.table(stack(setNames(paths, seq_along(paths))))
  setnames(nodes, c("values", "ind"), c("node", "lineage"))
  nodes[, step := seq_len(.N), by = lineage]
  nodes[, from := shift(node), by = lineage]
  edges <- unique(nodes[!is.na(from), .(from, to = node)])
  leaves <- vapply(paths, tail, integer(1L), n = 1L)
  pairs <- copy(pairs)
  pairs[, `:=`(
    program_loss = rowSums(pmax(persistent[from, , drop = FALSE] -
      program[to, , drop = FALSE], 0)^2),
    retained_program = rowSums(pmin(persistent[from, , drop = FALSE],
      program[to, , drop = FALSE])^2))]
  child_counts <- edges[, .N, by = from]
  current <- pairs[edges, on = .(from, to)]
  current[, parent_children := child_counts$N[match(from, child_counts$from)]]
  current[, retained_children := sum(program_loss <= retained_program), by = from]
  flagged <- switch(mode,
    fork = current[parent_children > 1L & retained_children > 0L &
      program_loss > retained_program],
    linear = current[parent_children == 1L & program_loss > retained_program],
    clock = current[program_loss > retained_program])
  switch(mode, clock = {
    switch(as.character(!is.null(clock_mask)), "TRUE" = {
      flagged <- flagged[clock_mask[cbind(from, to)]]
    })
  })
  flagged[, template := .I]
  leaf_matches <- CJ(from = unique(flagged$from), to = leaves)[from != to]
  leaf_matches[, `:=`(
    forward_loss = rowSums(pmax(persistent[from, , drop = FALSE] -
      program[to, , drop = FALSE], 0)^2),
    forward_retained = rowSums(pmin(persistent[from, , drop = FALSE],
      program[to, , drop = FALSE])^2),
    reverse_loss = rowSums(pmax(program[to, , drop = FALSE] -
      program[from, , drop = FALSE], 0)^2),
    reverse_retained = rowSums(pmin(program[to, , drop = FALSE],
      program[from, , drop = FALSE])^2),
    program_distance2 = rowSums((program[from, , drop = FALSE] -
      program[to, , drop = FALSE])^2))]
  leaf_matches <- leaf_matches[phase[from] < phase[to] &
    forward_loss <= forward_retained & reverse_loss <= reverse_retained][
      order(program_distance2, to), .SD[1L], by = from]
  ancestry <- merge(nodes[, .(lineage, from = node, from_step = step)],
    nodes[, .(lineage, to = node, to_step = step)], by = "lineage", allow.cartesian = TRUE)
  ancestry <- unique(ancestry[to_step - from_step > 1L,
    .(from, to, gap = to_step - from_step)])
  alternatives <- merge(pairs, ancestry, by = c("from", "to"))
  alternatives <- merge(alternatives,
    flagged[retained_children > 0L | from %in% leaf_matches$from,
      .(to, old_from = from)], by = "to")
  coherence <- flagged[, .(old_from = from, to,
    median_activity = numeric(.N), mean_activity = numeric(.N), n_cells = integer(.N))]
  switch(as.character(nrow(flagged) > 0L), "TRUE" = {
    expression <- cell_expression
    column <- rep(seq_len(ncol(expression)), diff(expression@p))
    gene <- expression@i + 1L
    expression@x <- pmin(pmax(log1p(expression@x) - profile[root, gene] -
      pmax(pseudotime[column] - phase[root], 0) * shared_slope[gene], 0),
      persistent[cbind(cluster[column], gene)])
    template <- persistent[flagged$from, , drop = FALSE] - 2 * pmin(
      persistent[flagged$from, , drop = FALSE], program[flagged$to, , drop = FALSE])
    activity <- as.data.table(as.matrix(Matrix::t(expression) %*% t(template)))
    setnames(activity, as.character(flagged$template))
    activity[, `:=`(cell = seq_len(.N), cluster = cluster)]
    activity <- melt(activity, id.vars = c("cell", "cluster"),
      variable.name = "template", value.name = "activity")
    activity[, template := as.integer(as.character(template))]
    activity <- merge(activity, flagged[, .(template, old_from = from, to)], by = "template")
    coherence <- activity[cluster == old_from, .(median_activity = median(activity),
      mean_activity = mean(activity), n_cells = .N), by = .(old_from, to)]
  })
  alternatives <- merge(alternatives, coherence, by = c("old_from", "to"), all.x = TRUE)
  alternatives[, acquired_compatible := TRUE]
  switch(mode, clock = {
    raw <- pmax(sweep(profile, 2L, profile[root, ], "-"), 0)
    alternatives[, acquired_compatible :=
      rowSums(pmax(raw[from, , drop = FALSE] - raw[to, , drop = FALSE], 0)^2) <=
      rowSums(pmin(raw[from, , drop = FALSE], raw[to, , drop = FALSE])^2)]
  })
  repairs <- alternatives[phase_eligible & program_loss <= retained_program &
    median_activity > 0 & acquired_compatible][order(gap, distance2, from), .SD[1L], by = to]
  edges[repairs, on = "to", from := i.from]
  orphans <- setdiff(unique(repairs$old_from), edges$from)
  continuation <- leaf_matches[from %in% orphans, .(orphan = from, leaf = to)]
  stopifnot(length(orphans) == nrow(continuation))
  switch(as.character(nrow(continuation) > 0L), "TRUE" = {
    repaired_paths <- rs_extract_lineages_from_edges(as.integer(edges$from),
      as.integer(edges$to), as.integer(nrow(profile)), as.integer(root))
    repaired_leaves <- vapply(repaired_paths, tail, integer(1L), n = 1L)
    stem_nodes <- as.data.table(stack(setNames(
      repaired_paths[match(continuation$orphan, repaired_leaves)], continuation$orphan)))
    setnames(stem_nodes, c("values", "ind"), c("node", "orphan"))
    stem_nodes[, orphan := as.integer(as.character(orphan))]
    stem_nodes[, step := seq_len(.N), by = orphan]
    degree <- tabulate(edges$from, nbins = nrow(profile))
    stem_nodes[, branch_step := max(step[node == root | degree[node] > 1L]), by = orphan]
    stems <- stem_nodes[step > branch_step, .(stem_head = first(node)), by = orphan]
    target_nodes <- as.data.table(stack(setNames(
      repaired_paths[match(continuation$leaf, repaired_leaves)], continuation$orphan)))
    setnames(target_nodes, c("values", "ind"), c("node", "orphan"))
    target_nodes[, orphan := as.integer(as.character(orphan))]
    target_nodes[, step := seq_len(.N), by = orphan]
    target_nodes[, from := shift(node), by = orphan]
    insertion <- target_nodes[!is.na(from) & phase[from] <= phase[orphan] &
      phase[node] > phase[orphan]][order(-step), .SD[1L], by = orphan]
    stopifnot(nrow(insertion) == nrow(continuation))
    setnames(insertion, c("from", "node"), c("insert_from", "insert_to"))
    insertion[, stem_head := stems$stem_head[match(orphan, stems$orphan)]]
    insertion[, orphan_phase := phase[orphan]]
    setorder(insertion, insert_from, insert_to, orphan_phase, orphan)
    insertion[, new_from := shift(orphan, fill = insert_from[1L]),
      by = .(insert_from, insert_to)]
    edges[insertion, on = c(to = "stem_head"), from := i.new_from]
    following <- insertion[, .(from = last(orphan), to = insert_to[1L]),
      by = .(insert_from, insert_to)]
    edges[following, on = "to", from := i.from]
  })
  paths <- rs_extract_lineages_from_edges(as.integer(edges$from), as.integer(edges$to),
    as.integer(nrow(profile)), as.integer(root))
  stopifnot(identical(sort(vapply(paths, tail, integer(1L), n = 1L)), sort(leaves)))
  list(lineage_paths = paths, repaired_edges = edges, repairs = repairs,
    commitment = coherence, edge_candidates = pairs, persistent_program = persistent)
}

.lineage_activity_identity_support <- function(profile, program, phase, root,
                                               pairs, activity_genes) {
  activity <- rowMeans(profile[, colnames(profile) %in% activity_genes, drop = FALSE])
  centered_phase <- phase - mean(phase)
  centered_activity <- activity - mean(activity)
  residual_activity <- centered_activity - centered_phase *
    sum(centered_phase * centered_activity) / sum(centered_phase^2)
  activity_slope <- pmax(as.numeric(crossprod(residual_activity, profile)) /
    sum(residual_activity^2), 0)
  identity <- pmax(program - outer(pmax(residual_activity - residual_activity[root], 0),
    activity_slope), 0)
  matched <- copy(pairs[phase_eligible == TRUE])
  matched[, coherent_retention := rowSums(pmin(identity[from, , drop = FALSE],
    identity[to, , drop = FALSE])^2) - rowSums(pmax(identity[from, , drop = FALSE] -
    identity[to, , drop = FALSE], 0)^2)]
  matched <- matched[order(-coherent_retention, distance2, to), .SD[1L], by = from]
  persistent <- identity * 0
  persistent[matched$from, ] <- pmin(identity[matched$from, , drop = FALSE],
    identity[matched$to, , drop = FALSE])
  loss <- rowSums(pmax(persistent[pairs$from, , drop = FALSE] -
    identity[pairs$to, , drop = FALSE], 0)^2)
  retained <- rowSums(pmin(persistent[pairs$from, , drop = FALSE],
    identity[pairs$to, , drop = FALSE])^2)
  mask <- matrix(FALSE, nrow(profile), nrow(profile))
  mask[cbind(pairs$from, pairs$to)] <- loss > retained
  list(mask = mask, activity = residual_activity, activity_slope = activity_slope,
    identity = identity, persistent = persistent, matched = matched,
    genes = colnames(profile)[colnames(profile) %in% activity_genes])
}

BuildRootRNAFeatures <- function(counts, root_cells, initial,
                                 n_components = 30L, expression_scale = 1e4) {
  expression <- as(counts, "CsparseMatrix")
  expression@x <- log1p(expression@x * rep(expression_scale /
    Matrix::colSums(expression), diff(expression@p)))
  root <- expression[, root_cells, drop = FALSE]
  gene_mean <- Matrix::rowMeans(root)
  mean_dot <- as.numeric(Matrix::crossprod(root, gene_mean))
  gram <- as.matrix(Matrix::crossprod(root))
  gram <- sweep(sweep(gram, 1L, mean_dot, "-"), 2L, mean_dot, "-")
  gram <- (gram + sum(gene_mean^2)) / (ncol(root) - 1L)
  start <- initial[root_cells] - mean(initial[root_cells])
  decomposition <- RSpectra::eigs_sym(gram, k = n_components, which = "LA",
    opts = list(initvec = start, tol = .Machine$double.eps))
  loadings <- as.matrix(root %*% decomposition$vectors) -
    tcrossprod(gene_mean, colSums(decomposition$vectors))
  loadings <- sweep(loadings, 2L,
    sqrt((ncol(root) - 1L) * decomposition$values), "/")
  features <- sweep(as.matrix(Matrix::crossprod(expression, loadings)), 2L,
    as.numeric(crossprod(gene_mean, loadings)), "-")
  colnames(features) <- colnames(loadings) <- paste0("pca_", seq_len(n_components))
  list(features = features, loadings = loadings, gene_mean = gene_mean,
    eigenvalues = decomposition$values, root_cells = sum(root_cells),
    total_root_variance = sum(diag(gram)))
}

BuildLineageProgramBasis <- function(counts, cluster_labels, pseudotime,
                                     lineage_paths, embedding,
                                     exclude_genes = character(), expression_scale = 1e4,
                                     activity_genes = NULL, reference_basis = NULL) {
  cl_map <- .compact_cluster_labels(cluster_labels)
  cluster <- cl_map$labels
  paths <- lapply(lineage_paths, match, table = cl_map$levels)
  root <- paths[[1L]][1L]
  size <- tabulate(cluster)
  counts <- as(counts, "CsparseMatrix")
  state_membership <- Matrix::sparseMatrix(i = cluster, j = seq_along(cluster), x = 1,
    dims = c(length(size), length(cluster)))
  bulk <- as.matrix(state_membership %*% Matrix::t(counts))
  intrinsic <- (bulk * (expression_scale / rowSums(bulk)))[,
    !colnames(bulk) %in% exclude_genes, drop = FALSE]
  counts <- counts[!rownames(counts) %in% exclude_genes, , drop = FALSE]
  log_intrinsic <- log(intrinsic)
  log_intrinsic[!is.finite(log_intrinsic)] <- NA_real_
  gene_reference <- switch(as.character(is.null(reference_basis)),
    "TRUE" = exp(colMeans(log_intrinsic, na.rm = TRUE)),
    "FALSE" = reference_basis$gene_reference[match(colnames(intrinsic),
      colnames(reference_basis$program_profiles))])
  state_ratios <- sweep(intrinsic, 2L, gene_reference, "/")
  ratios <- data.table(state = rep(seq_len(nrow(intrinsic)), each = ncol(intrinsic)),
    ratio = as.numeric(t(state_ratios)))
  state_factors <- ratios[is.finite(ratio) & ratio > 0, .(factor = median(ratio)), by = state]
  setorder(state_factors, state)
  profile <- log1p(intrinsic / state_factors$factor)
  dimnames(profile) <- list(as.character(cl_map$levels), colnames(intrinsic))
  cell_expression <- counts %*% Matrix::Diagonal(x = expression_scale / Matrix::colSums(counts))
  cell_ratios <- data.table(cell = rep(seq_len(ncol(cell_expression)), diff(cell_expression@p)),
    ratio = cell_expression@x / gene_reference[cell_expression@i + 1L])
  cell_factors <- cell_ratios[is.finite(ratio) & ratio > 0,
    .(factor = median(ratio)), by = cell]
  setorder(cell_factors, cell)
  cell_expression <- cell_expression %*% Matrix::Diagonal(x = 1 / cell_factors$factor)
  phase <- as.numeric(rowsum(pseudotime, cluster) / size)
  centers <- rowsum(embedding, cluster) / size
  phase_rank <- rank(phase, ties.method = "first")
  phase_rank[root] <- 0
  pairs <- CJ(from = seq_along(size), to = seq_along(size))[from != to]
  pairs[, `:=`(distance2 = rowSums((centers[to, , drop = FALSE] -
    centers[from, , drop = FALSE])^2), phase_eligible = phase_rank[from] < phase_rank[to])]
  raw_program <- pmax(sweep(profile, 2L, profile[root, ], "-"), 0)
  nodes <- as.data.table(stack(setNames(paths, seq_along(paths))))
  setnames(nodes, c("values", "ind"), c("node", "lineage"))
  nodes[, from := shift(node), by = lineage]
  edges <- unique(nodes[!is.na(from), .(from, to = node)])
  gene_columns <- as.character(seq_len(ncol(profile)))
  support <- as.data.table(raw_program[edges$to, , drop = FALSE])
  setnames(support, gene_columns)
  support[, from := edges$from]
  support <- dcast(melt(support, id.vars = "from"), from ~ variable,
    value.var = "value", fun.aggregate = max)
  persistent <- raw_program * 0
  persistent[support$from, ] <- pmin(raw_program[support$from, , drop = FALSE],
    as.matrix(support[, gene_columns, with = FALSE]))
  fork <- .repair_lineage_program_ancestry(profile, raw_program, persistent,
    phase, root, cluster, pseudotime, paths, pairs, cell_expression,
    numeric(ncol(profile)), "fork")
  support <- as.data.table(raw_program[pairs[phase_eligible == TRUE, to], , drop = FALSE])
  setnames(support, gene_columns)
  support[, from := pairs[phase_eligible == TRUE, from]]
  support <- dcast(melt(support, id.vars = "from"), from ~ variable,
    value.var = "value", fun.aggregate = max)
  persistent <- raw_program * 0
  persistent[support$from, ] <- pmin(raw_program[support$from, , drop = FALSE],
    as.matrix(support[, gene_columns, with = FALSE]))
  linear <- .repair_lineage_program_ancestry(profile, raw_program, persistent,
    phase, root, cluster, pseudotime, fork$lineage_paths, pairs, cell_expression,
    numeric(ncol(profile)), "linear")
  centered_phase <- phase - mean(phase)
  shared_slope <- switch(as.character(is.null(reference_basis)),
    "TRUE" = pmax(as.numeric(crossprod(centered_phase, profile) /
      sum(centered_phase^2)), 0),
    "FALSE" = reference_basis$shared_induction_slope[match(colnames(profile),
      colnames(reference_basis$program_profiles))])
  program <- pmax(raw_program - outer(pmax(phase - phase[root], 0), shared_slope), 0)
  matched <- copy(pairs[phase_eligible == TRUE])
  matched[, coherent_retention := rowSums(pmin(program[from, , drop = FALSE],
    program[to, , drop = FALSE])^2) - rowSums(pmax(program[from, , drop = FALSE] -
    program[to, , drop = FALSE], 0)^2)]
  matched <- matched[order(-coherent_retention, distance2, to), .SD[1L], by = from]
  persistent <- program * 0
  persistent[matched$from, ] <- pmin(program[matched$from, , drop = FALSE],
    program[matched$to, , drop = FALSE])
  activity_support <- NULL
  clock_mask <- NULL
  switch(as.character(is.null(activity_genes)), "TRUE" = {
    reference <- fread(system.file("extdata", "cell_cycle_activity.tsv", package = "LineaLinkR"))
    reference_ids <- unique(toupper(na.omit(unlist(reference, use.names = FALSE))))
    feature_ids <- toupper(sub("\\.[0-9]+$", "", colnames(profile)))
    activity_genes <- colnames(profile)[feature_ids %in% reference_ids]
  })
  switch(as.character(any(colnames(profile) %in% activity_genes)), "TRUE" = {
    activity_support <- .lineage_activity_identity_support(profile, program, phase,
      root, pairs, activity_genes)
    clock_mask <- activity_support$mask
  })
  clock <- .repair_lineage_program_ancestry(profile, program, persistent,
    phase, root, cluster, pseudotime, linear$lineage_paths, pairs, cell_expression,
    shared_slope, "clock", clock_mask)
  list(program_profiles = profile, induced_program = program,
    persistent_program = persistent, edge_candidates = clock$edge_candidates,
    repaired_edges = clock$repaired_edges, state_factors = state_factors$factor,
    gene_reference = gene_reference, shared_induction_slope = shared_slope,
    activity_identity_support = activity_support,
    cluster_levels = cl_map$levels, ancestor_stages = list(fork = fork, linear = linear,
      clock = clock))
}

.select_terminal_state_refinements <- function(program_profiles, cluster_labels,
                                                coarse_labels, lineage_paths, knn_conn) {
  cl_map <- .compact_cluster_labels(cluster_labels)
  coarse_map <- .compact_cluster_labels(coarse_labels)
  cluster <- cl_map$labels
  paths <- lapply(lineage_paths, match, table = cl_map$levels)
  profile <- as.matrix(program_profiles)
  states <- data.table(state = cluster, coarse = coarse_map$labels)[,
    .(coarse = coarse[1L], cells = .N), by = state][order(state)]
  adjacency <- as(knn_conn, "CsparseMatrix")
  diag(adjacency) <- 0
  fine_membership <- Matrix::sparseMatrix(i = cluster, j = seq_along(cluster), x = 1)
  coarse_membership <- Matrix::sparseMatrix(i = seq_along(cluster), j = coarse_map$labels, x = 1)
  contacts <- as.matrix(fine_membership %*% adjacency %*% coarse_membership)
  states[, local_edges := contacts[cbind(state, coarse)]]
  states[, external_edges := rowSums(contacts) - local_edges]
  states[, resident := local_edges > external_edges]
  leaves <- vapply(paths, tail, integer(1L), n = 1L)
  nodes <- data.table(state = unlist(paths), terminal = rep(leaves, lengths(paths)),
    step = sequence(lengths(paths)))
  shared <- merge(nodes, nodes, by = "state", allow.cartesian = TRUE)
  shared <- shared[terminal.x != terminal.y &
    states$coarse[terminal.x] == states$coarse[terminal.y]][
    order(-step.x), .SD[1L], by = .(terminal.x, terminal.y)]
  pairs <- shared[, .(from = terminal.x, to = terminal.y, anchor = state)]
  pairs[, `:=`(pair = seq_len(.N), supported_contrast = FALSE,
    both_resident = states$resident[from] & states$resident[to])]
  switch(as.character(nrow(pairs) > 0L), "TRUE" = {
    branch <- merge(pairs[, .(pair, terminal = from, anchor, to)],
      nodes, by = "terminal", allow.cartesian = TRUE)
    branch[, anchor_step := step[match(anchor[1L], state)], by = pair]
    branch <- branch[step > anchor_step]
    root <- paths[[1L]][1L]
    raw <- pmax(sweep(profile, 2L, profile[root, ], "-"), 0)
    acquired <- pmax(profile[pairs$from, , drop = FALSE] -
      profile[pairs$anchor, , drop = FALSE], 0)
    target <- pmax(profile[pairs$to, , drop = FALSE] -
      profile[pairs$anchor, , drop = FALSE], 0)
    branch_size <- branch[, .(states = .N), by = pair]
    pairs[, `:=`(private_states = branch_size$states[match(pair, branch_size$pair)],
      raw_loss = rowSums(pmax(raw[from, , drop = FALSE] - raw[to, , drop = FALSE], 0)^2),
      raw_retained = rowSums(pmin(raw[from, , drop = FALSE], raw[to, , drop = FALSE])^2),
      acquired_energy = rowSums(acquired^2),
      inherited_energy = rowSums(pmin(raw[from, , drop = FALSE], raw[anchor, , drop = FALSE])^2),
      acquired_loss = rowSums(pmax(acquired - target, 0)^2),
      acquired_retained = rowSums(pmin(acquired, target)^2))]
    pairs[, incompatible := raw_loss > raw_retained |
      (acquired_energy > inherited_energy & acquired_loss > acquired_retained)]
    reverse <- pairs[, .(from = to, to = from, reverse_incompatible = incompatible)]
    pairs[reverse, on = .(from, to),
      supported_contrast := both_resident & incompatible & i.reverse_incompatible]
  })
  terminals <- copy(states[match(leaves, state)])
  selection <- terminals[, .(candidate_terminals = .N, resident_terminals = sum(resident),
    supported_contrasts = 0L), by = coarse]
  contrast <- pairs[, .(supported_contrasts = sum(supported_contrast)),
    by = .(coarse = states$coarse[from])]
  selection[contrast, on = "coarse", supported_contrasts := i.supported_contrasts]
  selection[, selected := supported_contrasts > 0L]
  selection[, coarse := coarse_map$levels[coarse]]
  terminals[, `:=`(state = cl_map$levels[state], coarse = coarse_map$levels[coarse])]
  pairs[, `:=`(from = cl_map$levels[from], to = cl_map$levels[to], anchor = cl_map$levels[anchor])]
  list(selection = selection, terminals = terminals, pairs = pairs)
}

.partition_lineage_states <- function(cluster_labels, pseudotime, lineage_paths,
                                       adjacency, time, walk_steps = 4L) {
  root <- lineage_paths[[1L]][1L]
  states <- setdiff(sort(unique(cluster_labels)), root)
  rows <- split(which(cluster_labels %in% states),
    factor(cluster_labels[cluster_labels %in% states], levels = states))
  adjacency <- as(adjacency, "CsparseMatrix")
  diag(adjacency) <- 0
  symmetric <- (adjacency + Matrix::t(adjacency)) / 2
  local <- mapply("[", rep(list(symmetric), length(rows)), rows, rows,
    MoreArgs = list(drop = FALSE), SIMPLIFY = FALSE)
  graphs <- lapply(local, igraph::graph_from_adjacency_matrix,
    mode = "undirected", weighted = TRUE, diag = FALSE)
  partitions <- lapply(graphs, igraph::cluster_walktrap, steps = walk_steps)
  parts <- data.table(row = unlist(rows, use.names = FALSE),
    part = unlist(lapply(partitions, igraph::membership), use.names = FALSE))
  parts[, coarse := cluster_labels[row]]
  parts[, part_id := .GRP, by = .(coarse, part)]
  summary <- parts[, .(coarse = coarse[1L], cells = .N,
    phase = mean(pseudotime[row])), by = part_id][order(part_id)]
  nodes <- data.table(state = unlist(lineage_paths),
    lineage = rep(seq_along(lineage_paths), lengths(lineage_paths)))
  nodes[, from := shift(state), by = lineage]
  base_edges <- unique(nodes[!is.na(from), .(from, to = state)])
  fine <- Matrix::sparseMatrix(i = parts$part_id, j = parts$row, x = 1,
    dims = c(nrow(summary), length(cluster_labels)))
  coarse <- Matrix::sparseMatrix(i = seq_along(cluster_labels), j = cluster_labels, x = 1)
  contacts <- as.matrix(fine %*% adjacency %*% coarse)
  local_contacts <- as.matrix(fine %*% adjacency %*% Matrix::t(fine))
  local_contacts[outer(summary$coarse, summary$coarse, "!=")] <- 0
  summary[, parent := base_edges$from[match(coarse, base_edges$to)]]
  summary[, incoming := contacts[cbind(part_id, parent)]]
  early_contacts <- as.matrix(fine %*% Matrix::Diagonal(x = as.numeric(time == min(time))) %*%
    adjacency %*% coarse)
  summary[, entry_contact := early_contacts[cbind(part_id, parent)] / cells]
  summary[, unsupported := diag(local_contacts) == 0]
  summary[unsupported == TRUE, incoming := 0]
  entry <- summary[order(incoming == 0, -entry_contact, phase, part_id), .SD[1L], by = coarse]
  summary[, entry_part := entry$part_id[match(coarse, entry$coarse)]]
  components <- igraph::components(igraph::graph_from_adjacency_matrix(
    (local_contacts + t(local_contacts)) > 0, mode = "undirected", diag = FALSE))$membership
  summary[, unsupported := unsupported | components[part_id] != components[entry_part]]
  attachment_contacts <- local_contacts
  attachment_contacts[, summary$unsupported] <- 0
  attachment <- max.col(attachment_contacts, ties.method = "first")
  attachment[rowSums(attachment_contacts) == 0] <- summary$entry_part[rowSums(attachment_contacts) == 0]
  summary[, merged_part := fifelse(unsupported, attachment, part_id)]
  parts[, merged_part := summary$merged_part[part_id]]
  parts[, part := match(merged_part, sort(unique(merged_part))), by = coarse]
  mapping <- unique(parts[, .(coarse, merged_part, part)])[order(coarse, part)]
  mapping[, entry := merged_part == entry$part_id[match(coarse, entry$coarse)]]
  mapping[, state := coarse]
  mapping[entry == FALSE, state := max(cluster_labels) + seq_len(.N)]
  parts[mapping, on = .(coarse, part), state := i.state]
  local_rows <- split(parts$row, factor(parts$coarse, levels = states))
  local_parts <- split(parts$part, factor(parts$coarse, levels = states))
  local_graphs <- mapply("[", rep(list(adjacency), length(local_rows)), local_rows, local_rows,
    MoreArgs = list(drop = FALSE), SIMPLIFY = FALSE)
  local_roots <- mapping[entry == TRUE][match(states, coarse), part]
  trees <- mapply(BuildKnnClusterTree, local_graphs, local_parts, local_roots, SIMPLIFY = FALSE)
  names(trees) <- as.character(states)
  local_edges <- rbindlist(lapply(trees, "[[", "edges"), idcol = "coarse")
  local_edges[, coarse := as.integer(coarse)]
  local_edges[mapping, on = .(coarse, from = part), global_from := i.state]
  local_edges[mapping, on = .(coarse, to = part), global_to := i.state]
  mapping[, mapping_row := seq_len(.N)]
  parts[mapping, on = .(coarse, part), mapping_row := i.mapping_row]
  fine_ports <- Matrix::sparseMatrix(i = parts$mapping_row, j = parts$row, x = 1,
    dims = c(nrow(mapping), length(cluster_labels)))
  port_contacts <- as.matrix(fine_ports %*% adjacency %*% coarse)
  ports <- merge(base_edges, mapping, by.x = "from", by.y = "coarse", allow.cartesian = TRUE)
  ports[, contact := port_contacts[cbind(mapping_row, to)]]
  ports <- ports[order(-contact, part), .SD[1L], by = .(from, to)]
  edges <- copy(base_edges)
  edges[ports, on = .(from, to), from := i.state]
  edges <- rbind(edges, local_edges[, .(from = global_from, to = global_to)])
  cluster <- cluster_labels
  cluster[parts$row] <- parts$state
  paths <- rs_extract_lineages_from_edges(as.integer(edges$from), as.integer(edges$to),
    as.integer(max(cluster)), as.integer(root))
  list(cluster_labels = cluster, lineage_paths = paths, edges = edges,
    partition = mapping, child_ports = ports)
}

RefinePersistentLineageStatesFromCounts <- function(counts, cluster_labels, pseudotime,
                                                     coarse_model, embedding, knn_conn,
                                                     time, condition = rep("all", length(cluster_labels)),
                                                     expression_scale = 1e4, sigma = 1,
                                                     walk_steps = 4L) {
  proposal <- .partition_lineage_states(cluster_labels, pseudotime,
    coarse_model$lineage_paths, knn_conn, time, walk_steps)
  root <- proposal$lineage_paths[[1L]][1L]
  joint_weights <- ComputeLineageWeightsFromKnnGraph(NULL, embedding,
    proposal$cluster_labels, root, sigma = sigma, tree = list(edges = proposal$edges))
  exclude <- setdiff(rownames(counts), colnames(coarse_model$program_basis$program_profiles))
  joint <- RefineLineageProgramTopologyFromCounts(counts, proposal$cluster_labels, pseudotime,
    proposal$lineage_paths, joint_weights, embedding, exclude_genes = exclude,
    expression_scale = expression_scale, sigma = sigma, count_identity = TRUE,
    knn_conn = knn_conn, reference_basis = coarse_model$program_basis, coarse_labels = cluster_labels)
  selection <- SelectPersistentLineageRefinements(joint$program_basis$program_profiles,
    proposal$cluster_labels, cluster_labels, joint$lineage_paths, coarse_model$lineage_paths,
    knn_conn, time, condition)
  accepted <- selection$selection[selected == TRUE, coarse]
  mapping <- data.table(state = proposal$cluster_labels, coarse = cluster_labels)[,
    .(coarse = coarse[1L]), by = state][order(state)]
  mapping[, selected_state := fifelse(coarse %in% accepted, state, coarse)]
  mapping[, selected_state := match(selected_state, sort(unique(selected_state))) ]
  cluster <- mapping$selected_state[match(proposal$cluster_labels, mapping$state)]
  edges <- copy(proposal$edges)
  edges[, `:=`(from = mapping$selected_state[match(from, mapping$state)],
    to = mapping$selected_state[match(to, mapping$state)])]
  edges <- unique(edges[from != to])
  root <- mapping$selected_state[match(root, mapping$state)]
  paths <- rs_extract_lineages_from_edges(as.integer(edges$from), as.integer(edges$to),
    as.integer(max(cluster)), as.integer(root))
  weights <- ComputeLineageWeightsFromKnnGraph(NULL, embedding, cluster, root,
    sigma = sigma, tree = list(edges = edges))
  result <- RefineLineageProgramTopologyFromCounts(counts, cluster, pseudotime,
    paths, weights, embedding, exclude_genes = exclude, expression_scale = expression_scale,
    sigma = sigma, count_identity = TRUE, knn_conn = knn_conn,
    reference_basis = coarse_model$program_basis, coarse_labels = cluster_labels)
  result$cluster_labels <- cluster
  result$coarse_labels <- cluster_labels
  result$refinement_selection <- selection
  result$selected_coarse <- accepted
  result$state_mapping <- mapping
  result$graft_paths <- paths
  result$graft_weights <- weights
  result
}

SelectPersistentLineageRefinements <- function(program_profiles, cluster_labels,
                                                 coarse_labels, lineage_paths,
                                                 coarse_paths, adjacency, time,
                                                 condition = rep("all", length(cluster_labels))) {
  fine <- data.table(state = cluster_labels, coarse = coarse_labels)[,
    .(coarse = coarse[1L], cells = .N), by = state][order(state)]
  coarse_levels <- sort(unique(coarse_labels))
  profile <- as.matrix(program_profiles[as.character(fine$state), , drop = FALSE])
  adjacency <- as(adjacency, "CsparseMatrix")
  diag(adjacency) <- 0
  identity <- .select_terminal_state_refinements(profile, cluster_labels, coarse_labels,
    lineage_paths, adjacency)
  pairs <- copy(identity$pairs)
  pairs[, `:=`(from_coarse = fine$coarse[match(from, fine$state)],
    to_coarse = fine$coarse[match(to, fine$state)], net_supported = FALSE)]
  switch(as.character(nrow(pairs) > 0L), "TRUE" = {
    pairs[, signed_energy := rowSums((profile[as.character(from), , drop = FALSE] -
      profile[as.character(anchor), , drop = FALSE])^2)]
    pairs[, `:=`(contrast = acquired_loss > acquired_retained,
      net_acquired = acquired_energy > inherited_energy &
        acquired_energy > signed_energy - acquired_energy)]
    reverse <- pairs[, .(from = to, to = from, anchor,
      reverse_contrast = contrast, reverse_net_acquired = net_acquired)]
    pairs[reverse, on = .(from, to, anchor),
      net_supported := both_resident & contrast & i.reverse_contrast &
        (net_acquired | i.reverse_net_acquired)]
  })
  pairs <- pairs[from_coarse == to_coarse]
  fine_membership <- Matrix::sparseMatrix(i = seq_along(cluster_labels),
    j = match(cluster_labels, fine$state), x = 1)
  coarse_membership <- Matrix::sparseMatrix(i = seq_along(coarse_labels),
    j = match(coarse_labels, coarse_levels), x = 1)
  contacts <- as.matrix(Matrix::crossprod(fine_membership, adjacency %*% coarse_membership))
  coarse_nodes <- data.table(coarse = unlist(coarse_paths),
    path = rep(seq_along(coarse_paths), lengths(coarse_paths)))
  coarse_nodes[, parent := shift(coarse), by = path]
  edges <- unique(coarse_nodes[!is.na(parent), .(parent, child = coarse)])
  ports <- merge(fine, edges, by.x = "coarse", by.y = "parent", allow.cartesian = TRUE)
  ports[, contact := contacts[cbind(match(state, fine$state), match(child, coarse_levels))]]
  ports[, density := contact / cells]
  best <- ports[contact > 0][order(-density, state), .SD[1L], by = .(coarse, child)]
  resolved <- best[, .(exits = uniqueN(child), distinct_ports = uniqueN(state)),
    by = coarse][exits > 1L & distinct_ports > 1L, coarse]
  cells <- data.table(row = seq_along(cluster_labels), state = cluster_labels,
    coarse = coarse_labels, time = time, condition = condition)
  cells[, late := time == max(time), by = condition]
  cells[, context := .GRP, by = .(time, condition)]
  cells[, group := .GRP, by = .(state, context)]
  occupancy <- cells[, .(state = state[1L], coarse = coarse[1L], context = context[1L],
    late = late[1L], cells = .N), by = group][order(group)]
  terminals <- vapply(lineage_paths, tail, integer(1L), n = 1L)
  nodes <- data.table(state = unlist(lineage_paths),
    terminal = rep(terminals, lengths(lineage_paths)))[,
    .(paths = uniqueN(terminal), terminal = terminal[1L]), by = state]
  arm_groups <- merge(nodes[paths == 1L], occupancy[late == TRUE, .(state, group, context)],
    by = "state", allow.cartesian = TRUE)
  arm_groups[, arm_context := .GRP, by = .(terminal, context)]
  group_membership <- Matrix::sparseMatrix(i = cells$row, j = cells$group, x = 1,
    dims = c(nrow(cells), nrow(occupancy)))
  group_edges <- as.data.table(Matrix::summary(Matrix::crossprod(group_membership,
    adjacency %*% group_membership)))
  group_edges <- group_edges[occupancy$context[i] == occupancy$context[j] & occupancy$late[i]]
  contextual_graph <- Matrix::sparseMatrix(i = group_edges$i, j = group_edges$j,
    x = group_edges$x, dims = c(nrow(occupancy), nrow(occupancy)))
  arm_membership <- Matrix::sparseMatrix(i = arm_groups$group, j = arm_groups$arm_context,
    x = 1, dims = c(nrow(occupancy), uniqueN(arm_groups$arm_context)))
  support <- arm_groups[, .(terminal = terminal[1L], context = context[1L]),
    by = arm_context][order(arm_context)]
  support[, cells := as.numeric(Matrix::crossprod(arm_membership, occupancy$cells))]
  support[, inside := as.numeric(Matrix::colSums(arm_membership *
    (contextual_graph %*% arm_membership)))]
  support[, total := as.numeric(Matrix::crossprod(arm_membership, Matrix::rowSums(contextual_graph)))]
  support[, outside := total - inside]
  support[, resident := inside > outside]
  contexts <- unique(cells[, .(context, condition)])
  support[contexts, on = "context", condition := i.condition]
  supported <- support[, .(late_private_resident = any(resident)), by = terminal]
  terminal_table <- fine[state %in% terminals, .(terminal = state, coarse,
    late_private_resident = FALSE)]
  terminal_table[supported, on = "terminal", late_private_resident := i.late_private_resident]
  pairs[, temporal_support := net_supported &
    terminal_table$late_private_resident[match(from, terminal_table$terminal)] &
    terminal_table$late_private_resident[match(to, terminal_table$terminal)]]
  acquired <- pairs[, .(acquired_identity = any(temporal_support)), by = .(coarse = from_coarse)]
  local <- terminal_table[, .(local_terminals = .N,
    unsupported_terminals = sum(!late_private_resident)), by = coarse]
  selection <- data.table(coarse = coarse_levels, child_contact = coarse_levels %in% resolved,
    local_terminals = 0L, unsupported_terminals = 0L, acquired_identity = FALSE)
  selection[local, on = "coarse", `:=`(local_terminals = i.local_terminals,
    unsupported_terminals = i.unsupported_terminals)]
  selection[acquired, on = "coarse", acquired_identity := i.acquired_identity]
  selection[, structural_ports := child_contact & unsupported_terminals == 0L]
  selection[, `:=`(mechanism = "persistent_identity_or_structural_ports",
    selected = acquired_identity | structural_ports)]
  list(selection = selection, contexts = support, terminals = terminal_table,
    pairs = pairs, ports = ports)
}

AllocateRefinedLineageWeights <- function(embedding, cluster_labels, lineage_paths,
                                            lineage_weights, prediction, coarse_labels,
                                            coarse_paths, coarse_weights, coarse_prediction,
                                            method = "private_entry", time = NULL,
                                            refinement_roots = integer(),
                                            state_support = NULL, deferred_weights = NULL,
                                            deferred_paths = coarse_paths) {
  paths <- lineage_paths
  terminals <- vapply(paths, tail, integer(1L), n = 1L)
  old_terminals <- vapply(coarse_paths, tail, integer(1L), n = 1L)
  mapping <- data.table(state = cluster_labels, coarse = coarse_labels)[,
    .(coarse = coarse[1L]), by = state]
  nodes <- data.table(state = unlist(paths), lineage = rep(seq_along(paths), lengths(paths)),
    step = sequence(lengths(paths)))
  old_nodes <- data.table(state = unlist(coarse_paths),
    parent = rep(seq_along(coarse_paths), lengths(coarse_paths)),
    step = sequence(lengths(coarse_paths)))
  base <- as.matrix(coarse_weights)
  switch(as.character(!is.null(deferred_weights)), "TRUE" = {
    state_support <- copy(as.data.table(state_support))
    state_support[mapping, on = "state", coarse := i.coarse]
    deferred_nodes <- data.table(state = unlist(deferred_paths),
      parent = rep(seq_along(deferred_paths), lengths(deferred_paths)))
    deferred_terminals <- vapply(deferred_paths, tail, integer(1L), n = 1L)
    private <- deferred_nodes[, .(paths = uniqueN(parent)), by = state][
      paths == 1L & !state %in% deferred_terminals, state]
    state_support[, defer := private_interior & margin <= 0 & coarse %in% private]
    rows <- which(cluster_labels %in% state_support[defer == TRUE, state])
    base[rows, ] <- deferred_weights[rows, , drop = FALSE]
  })
  switch(method, private_entry = {
    children <- data.table(lineage = seq_along(paths), terminal = terminals,
      coarse = mapping$coarse[match(terminals, mapping$state)])
    children <- merge(children, old_nodes[, .(coarse = state, parent)],
      by = "coarse", allow.cartesian = TRUE)
    setorder(children, lineage)
    children[, family_size := .N, by = parent]
    refined <- children[family_size > 1L, lineage]
    nodes[children, on = "lineage", parent := i.parent]
    shared <- nodes[lineage %in% refined, .(lineages = uniqueN(lineage), step = max(step)),
      by = .(parent, state)]
    shared[, family_size := children$family_size[match(parent, children$parent)]]
    roots <- shared[lineages == family_size][order(-step), .SD[1L], by = parent]
    local_nodes <- merge(nodes[lineage %in% refined],
      roots[, .(parent, root = state, root_step = step)], by = "parent")
    local_nodes <- local_nodes[step >= root_step]
    local_nodes[, inherited_root := mapping$coarse[match(root, mapping$state)]]
    prefix <- merge(old_nodes[parent %in% roots$parent],
      unique(local_nodes[, .(parent, inherited_root)]), by = "parent")
    prefix[, entry_step := step[match(inherited_root[1L], state)], by = parent]
    prefix <- merge(prefix[step < entry_step, .(parent, state, step)],
      children[lineage %in% refined, .(parent, lineage)],
      by = "parent", allow.cartesian = TRUE)
    inherited <- merge(old_nodes, children[family_size == 1L, .(parent, lineage)], by = "parent")
    grafted <- rbind(inherited[, .(lineage, state, section = 1L, step)],
      prefix[, .(lineage, state, section = 1L, step)],
      local_nodes[, .(lineage, state, section = 2L, step)])
    setorder(grafted, lineage, section, step)
    paths <- unname(split(grafted$state, grafted$lineage))
    parent_scores <- as.matrix(lineage_weights)
    parent_scores[outer(coarse_prediction, children$parent, "!=")] <- -Inf
    constrained <- max.col(parent_scores, ties.method = "first")
    accepted <- children$parent[prediction] == coarse_prediction
    constrained[accepted] <- prediction[accepted]
    prediction <- constrained
    anchors <- data.table(row = seq_along(cluster_labels), state = cluster_labels,
      lineage = prediction, parent = children$parent[prediction])
    anchors <- merge(anchors, unique(local_nodes[, .(state, parent)]),
      by = c("state", "parent"))
    private <- unique(grafted[, .(state, lineage)])[,
      .(lineages = .N, lineage = lineage[1L]), by = state][lineages == 1L]
    training <- merge(anchors, private[, .(state, lineage)], by = c("state", "lineage"))
    training[, time := time[row]]
    training[, entry_time := min(time), by = lineage]
    training <- training[time == entry_time]
    compatible <- children[, .(lineage, parent)]
  }, native = {
    known <- data.table(lineage = seq_along(paths),
      parent = match(mapping$coarse[match(terminals, mapping$state)], old_terminals))[
      !is.na(parent)]
    support <- unique(merge(nodes[, .(state, lineage)], known, by = "lineage"))
    compatible <- merge(nodes, support[, .(state, parent)],
      by = "state", allow.cartesian = TRUE)
    compatible <- unique(compatible[, .SD[step == max(step)], by = lineage][, .(lineage, parent)])
    family <- compatible[, .(children = .N), by = parent]
    refined <- sort(unique(compatible[parent %in% family[children > 1L, parent], lineage]))
    private <- nodes[, .(paths = uniqueN(lineage), lineage = lineage[1L]), by = state][paths == 1L]
    labels <- private$lineage[match(cluster_labels, private$state)]
    local_paths <- merge(nodes,
      nodes[state %in% refinement_roots, .(lineage, root_step = step)],
      by = "lineage", allow.cartesian = TRUE)
    rows <- which(cluster_labels %in% local_paths[step >= root_step, state])
    labels[rows] <- prediction[rows]
    anchors <- data.table(row = seq_along(cluster_labels), lineage = labels)[!is.na(lineage)]
    training <- anchors[lineage %in% refined]
  })
  setorder(compatible, parent, lineage)
  compatible[, slot := .I]
  scores <- matrix(0, nrow = length(cluster_labels), ncol = length(paths))
  local <- anchors[lineage %in% refined]
  centers <- rowsum(embedding[local$row, , drop = FALSE], local$lineage) /
    as.numeric(table(factor(local$lineage, levels = refined)))
  switch(method, private_entry = {
    residual <- embedding[local$row, , drop = FALSE] -
      centers[match(local$lineage, refined), , drop = FALSE]
    covariance <- crossprod(residual) / nrow(residual)
    centers <- rowsum(embedding[training$row, , drop = FALSE], training$lineage) /
      as.numeric(table(factor(training$lineage, levels = refined)))
    directions <- solve(covariance, t(centers))
    scores[, refined] <- sweep(embedding %*% directions, 2L,
      rowSums(centers * t(directions)) / 2, "-")
  }, native = {
    origin <- colMeans(centers)
    directions <- sweep(centers, 2L, origin, "-")
    directions <- directions / sqrt(rowSums(directions^2))
    pairs <- CJ(left = seq_along(refined), right = seq_along(refined))[left < right]
    delta <- directions[pairs$left, , drop = FALSE] - directions[pairs$right, , drop = FALSE]
    projected <- embedding[local$row, , drop = FALSE] %*% t(delta)
    scale <- sqrt(mean(colMeans(projected^2) - colMeans(projected)^2))
    scores[, refined] <- sweep(embedding %*% t(directions), 2L,
      as.vector(directions %*% origin), "-") / scale
  })
  long <- data.table(row = rep(seq_along(cluster_labels), times = nrow(compatible)),
    parent = rep(compatible$parent, each = length(cluster_labels)),
    score = as.vector(scores[, compatible$lineage, drop = FALSE]))
  long[, weight := exp(score - max(score)), by = .(row, parent)]
  long[, weight := weight / sum(weight), by = .(row, parent)]
  conditional <- matrix(long$weight, nrow = length(cluster_labels))
  reference_conditional <- copy(conditional)
  owned <- merge(anchors[, .(row, lineage)], compatible, by = "lineage", allow.cartesian = TRUE)
  mask <- merge(unique(owned[, .(row, parent)]), compatible[, .(parent, slot)],
    by = "parent", allow.cartesian = TRUE)
  reference_conditional[cbind(mask$row, mask$slot)] <- 0
  reference_conditional[cbind(owned$row, owned$slot)] <- 1
  slot_weights <- conditional * base[, compatible$parent, drop = FALSE]
  slot_reference <- reference_conditional * base[, compatible$parent, drop = FALSE]
  membership <- outer(compatible$lineage, seq_along(paths), "==") * 1
  list(lineage_paths = paths, lineage_weights = slot_weights %*% membership,
    reference_weights = slot_reference %*% membership, prediction = prediction,
    coarse_weights = base, compatibility = compatible, anchors = anchors,
    training_anchors = training, state_support = state_support, method = method)
}

ResolveLineageRefinementByCondition <- function(lineage_weights, reference_weights,
                                                cluster_labels, coarse_labels, lineage_paths,
                                                condition, time,
                                                terminal_order = vapply(lineage_paths, tail, integer(1L), n = 1L)) {
  weights <- as.matrix(lineage_weights)
  reference <- as.matrix(reference_weights)
  states <- data.table(state = cluster_labels, coarse = coarse_labels)[,
    .(coarse = coarse[1L]), by = state]
  terminals <- data.table(terminal = terminal_order, column = seq_along(terminal_order),
    coarse = states$coarse[match(terminal_order, states$state)])
  terminals[, family_size := .N, by = coarse]
  selected <- terminals[family_size > 1L]
  support_table <- data.table()
  switch(as.character(nrow(selected) > 0L), "TRUE" = {
    families <- unique(selected$coarse)
    group <- match(selected$coarse, families)
    membership <- outer(group, seq_along(families), "==") * 1
    leaves <- vapply(lineage_paths, tail, integer(1L), n = 1L)
    nodes <- data.table(state = unlist(lineage_paths), terminal = rep(leaves, lengths(lineage_paths)))[,
      .(paths = uniqueN(terminal), terminal = terminal[1L]), by = state]
    private <- nodes[paths == 1L]
    cell_terminal <- private$terminal[match(cluster_labels, private$state)]
    conditions <- unique(condition)
    condition_id <- match(condition, conditions)
    latest <- data.table(condition = condition_id, time = time)[, .(time = max(time)), by = condition]
    late <- time == latest$time[match(condition_id, latest$condition)]
    support <- unclass(table(factor(condition_id[late], levels = seq_along(conditions)),
      factor(match(cell_terminal[late], selected$terminal), levels = seq_len(nrow(selected)))))
    supported_children <- (support > 0) %*% membership
    active <- support > 0 | supported_children[, group, drop = FALSE] == 0
    choices <- data.table(condition = rep(seq_along(conditions), each = nrow(selected)),
      child = rep(seq_len(nrow(selected)), times = length(conditions)),
      family = rep(group, times = length(conditions)), active = as.vector(t(active)))[
      active == TRUE, .(child = child[1L]), by = .(condition, family)]
    destination <- matrix(0, nrow = length(conditions), ncol = nrow(selected))
    destination[cbind(choices$condition, choices$child)] <- 1
    local <- weights[, selected$column, drop = FALSE]
    coarse_mass <- local %*% membership
    local <- local * active[condition_id, , drop = FALSE]
    denominator <- (local %*% membership)[, group, drop = FALSE]
    local <- local / replace(denominator, denominator == 0, 1) * coarse_mass[, group, drop = FALSE]
    collapse <- supported_children[condition_id, group, drop = FALSE] <= 1
    destinations <- destination[condition_id, , drop = FALSE]
    weights[, selected$column] <- ifelse(collapse,
      coarse_mass[, group, drop = FALSE] * destinations, local)
    reference_mass <- reference[, selected$column, drop = FALSE] %*% membership
    reference[, selected$column] <- ifelse(collapse,
      reference_mass[, group, drop = FALSE] * destinations,
      reference[, selected$column, drop = FALSE])
    support_table <- data.table(condition = rep(conditions, each = nrow(selected)),
      terminal = rep(selected$terminal, times = length(conditions)),
      coarse = rep(selected$coarse, times = length(conditions)),
      cells = as.vector(t(support)),
      resolved_children = as.vector(t(supported_children[, group, drop = FALSE])))
  })
  list(lineage_weights = weights, reference_weights = reference,
    terminal_order = terminal_order, branch_support = support_table)
}

.restrict_supported_lineage_families <- function(lineage_weights, cluster_labels,
                                                 lineage_paths, state_support,
                                                 terminal_order = vapply(lineage_paths,
                                                   tail, integer(1L), n = 1L)) {
  leaves <- vapply(lineage_paths, tail, integer(1L), n = 1L)
  nodes <- data.table(state = unlist(lineage_paths),
    lineage = rep(seq_along(leaves), lengths(lineage_paths)),
    terminal = rep(leaves, lengths(lineage_paths)))
  nodes[, `:=`(step = seq_len(.N), parent = shift(state)), by = lineage]
  private <- nodes[, .(n_paths = uniqueN(lineage), terminal = terminal[1L]),
    by = state][n_paths == 1L & terminal %in% terminal_order]
  private <- merge(private, as.data.table(state_support)[margin > 0, .(state)], by = "state")
  forks <- unique(nodes[!is.na(parent), .(parent, state)])[, .N, by = parent][N > 1L, parent]
  ancestors <- merge(nodes[state %in% private$state, .(lineage, state, step)],
    nodes[state %in% forks, .(lineage, fork = state, fork_step = step)],
    by = "lineage", allow.cartesian = TRUE)
  ancestors <- ancestors[fork_step < step][order(-fork_step), .SD[1L], by = state]
  family <- merge(ancestors[, .(state, fork)], nodes[, .(fork = state, terminal)],
    by = "fork", allow.cartesian = TRUE)
  supported_states <- unique(family$state)
  allowed <- matrix(0, length(supported_states), length(terminal_order))
  family <- family[terminal %in% terminal_order]
  allowed[cbind(match(family$state, supported_states),
    match(family$terminal, terminal_order))] <- 1
  rows <- which(cluster_labels %in% supported_states)
  columns <- seq_along(terminal_order)
  weights <- as.matrix(lineage_weights)
  local <- weights[rows, columns, drop = FALSE] *
    allowed[match(cluster_labels[rows], supported_states), , drop = FALSE]
  mass <- rowSums(local)
  weights[rows, columns] <- local / fifelse(mass > 0, mass, 1)
  weights
}

RefineLineageProgramTopologyFromCounts <- function(counts, cluster_labels, pseudotime,
                                                   lineage_paths, lineage_weights,
                                                   embedding, exclude_genes = character(),
                                                   expression_scale = 1e4, sigma = 1,
                                                   count_identity = FALSE,
                                                   knn_conn = NULL,
                                                   activity_genes = NULL,
                                                   reference_basis = NULL,
                                                   coarse_labels = NULL) {
  basis <- BuildLineageProgramBasis(counts, cluster_labels, pseudotime,
    lineage_paths, embedding, exclude_genes, expression_scale, activity_genes, reference_basis)
  cl_map <- .compact_cluster_labels(cluster_labels)
  paths <- lapply(lineage_paths, match, table = cl_map$levels)
  root <- paths[[1L]][1L]
  edges <- basis$repaired_edges
  switch(as.character(!is.null(coarse_labels)), "TRUE" = {
    coarse <- data.table(state = cl_map$labels, coarse = coarse_labels)[,
      .(coarse = coarse[1L]), by = state][match(seq_len(nrow(basis$program_profiles)), state), coarse]
    repairs <- as.data.table(basis$ancestor_stages$clock$repairs)
    local <- repairs[coarse[old_from] == coarse[to], .(from = old_from, to)]
    linear_edges <- basis$ancestor_stages$linear$repaired_edges
    local[, parent := linear_edges$from[match(from, linear_edges$to)]]
    profile <- basis$program_profiles
    program <- basis$induced_program
    raw <- pmax(sweep(profile, 2L, profile[root, ], "-"), 0)
    acquired <- pmax(program[local$from, , drop = FALSE] -
      program[local$parent, , drop = FALSE], 0)
    local[, `:=`(raw_loss = rowSums(pmax(raw[from, , drop = FALSE] - raw[to, , drop = FALSE], 0)^2),
      raw_retained = rowSums(pmin(raw[from, , drop = FALSE], raw[to, , drop = FALSE])^2),
      acquired_energy = rowSums(acquired^2),
      inherited_energy = rowSums(pmin(program[from, , drop = FALSE], program[parent, , drop = FALSE])^2))]
    restored <- local[raw_loss <= raw_retained & acquired_energy <= inherited_energy]
    edges[restored, on = "to", from := i.from]
    basis$repaired_edges <- edges
    basis$ancestor_stages$clock$repaired_edges <- edges
    basis$ancestor_stages$clock$repairs <- repairs[!to %in% restored$to]
    basis$refinement_continuity <- restored
  })
  repaired_paths <- rs_extract_lineages_from_edges(as.integer(edges$from),
    as.integer(edges$to), as.integer(nrow(basis$program_profiles)), as.integer(root))
  basis$ancestor_stages$clock$lineage_paths <- repaired_paths
  weights <- matrix(rs_compute_lineage_weights_trunk_aware_from_edges(
    as.numeric(t(embedding)), as.integer(cl_map$labels), as.integer(nrow(embedding)),
    as.integer(ncol(embedding)), as.integer(nrow(basis$program_profiles)),
    as.integer(edges$from), as.integer(edges$to), as.integer(root), as.numeric(sigma)),
    nrow = nrow(embedding), ncol = length(repaired_paths), byrow = TRUE)
  old_leaves <- vapply(paths, tail, integer(1L), n = 1L)
  new_leaves <- vapply(repaired_paths, tail, integer(1L), n = 1L)
  old_nodes <- as.data.table(stack(setNames(paths, old_leaves)))
  setnames(old_nodes, c("values", "ind"), c("state", "terminal"))
  old_nodes[, `:=`(terminal = as.integer(as.character(terminal)), old = TRUE)]
  new_nodes <- as.data.table(stack(setNames(repaired_paths, new_leaves)))
  setnames(new_nodes, c("values", "ind"), c("state", "terminal"))
  new_nodes[, `:=`(terminal = as.integer(as.character(terminal)), new = TRUE)]
  membership <- merge(unique(old_nodes), unique(new_nodes),
    by = c("state", "terminal"), all = TRUE)
  affected <- unique(membership[is.na(old) | is.na(new), state])
  unchanged <- !cl_map$labels %in% affected
  weights[unchanged, ] <- lineage_weights[unchanged, match(new_leaves, old_leaves), drop = FALSE]
  repaired_paths <- relist(cl_map$levels[unlist(repaired_paths)], repaired_paths)
  result <- RefineLineageProgramTopology(basis$program_profiles, basis$induced_program,
    cluster_labels, pseudotime, repaired_paths, weights, embedding, basis$edge_candidates,
    cell_counts = counts[match(colnames(basis$program_profiles), rownames(counts)), , drop = FALSE],
    sigma = sigma, knn_conn = knn_conn)
  result$program_basis <- basis
  switch(as.character(!is.null(coarse_labels) & !is.null(knn_conn)), "TRUE" = {
    result$refinement_identity <- .select_terminal_state_refinements(
      basis$program_profiles, cluster_labels, coarse_labels, result$lineage_paths, knn_conn)
  })
  switch(as.character(count_identity), "TRUE" = {
    result$count_identity <- .persistent_count_identity(
      counts[match(colnames(basis$program_profiles), rownames(counts)), , drop = FALSE],
      cluster_labels, result$lineage_paths, result$prediction)
    result$prediction <- result$count_identity$prediction
  })
  switch(as.character(!is.null(knn_conn)), "TRUE" = {
    result$rna_identity <- ResolveLineageRNAIdentities(
      counts[match(colnames(basis$program_profiles), rownames(counts)), , drop = FALSE],
      cluster_labels, vapply(result$lineage_paths, tail, integer(1L), n = 1L),
      knn_conn, result$prediction, program_profiles = basis$program_profiles,
      lineage_paths = result$lineage_paths)
  })
  result
}

ResolveLineageRNAIdentities <- function(counts, cluster_labels, terminal_states,
                                        knn_conn, prediction, template_counts = NULL,
                                        aggregate_pseudocount = 0.5,
                                        template_membership = NULL,
                                        program_profiles = NULL, lineage_paths = NULL) {
  cl_map <- .compact_cluster_labels(cluster_labels)
  cluster <- cl_map$labels
  terminals <- match(terminal_states, cl_map$levels)
  template_count <- length(terminals)
  template_index <- seq_len(template_count)
  anchor_code <- match(cluster, terminals)
  anchor_rows <- which(!is.na(anchor_code))
  anchor_membership <- switch(as.character(is.null(template_membership)),
    "TRUE" = Matrix::sparseMatrix(i = anchor_rows, j = anchor_code[anchor_rows], x = 1,
      dims = c(length(cluster), template_count)),
    "FALSE" = template_membership)
  aggregate <- switch(as.character(is.null(template_counts)),
    "TRUE" = as.matrix(counts %*% anchor_membership),
    "FALSE" = as.matrix(template_counts))
  smoothed <- aggregate + aggregate_pseudocount
  composition <- sweep(smoothed, 2L, colSums(smoothed), "/")
  expression_affinity <- crossprod(sqrt(composition))
  diag(expression_affinity) <- -Inf
  expression_neighbour <- max.col(expression_affinity, ties.method = "first")
  membership <- Matrix::sparseMatrix(i = cluster, j = seq_along(cluster), x = 1)
  quotient <- as.matrix(membership %*% knn_conn %*% Matrix::t(membership))
  quotient <- (quotient + t(quotient)) / 2
  diag(quotient) <- 0
  laplacian <- diag(rowSums(quotient), nrow = nrow(quotient)) - quotient
  interior <- setdiff(seq_len(nrow(laplacian)), terminals)
  reduced <- laplacian[terminals, terminals, drop = FALSE]
  switch(as.character(length(interior) > 0L), "TRUE" = {
    reduced <- reduced - laplacian[terminals, interior, drop = FALSE] %*%
      solve(laplacian[interior, interior, drop = FALSE], laplacian[interior, terminals, drop = FALSE])
  })
  graph_affinity <- pmax(-reduced, 0)
  diag(graph_affinity) <- 0
  graph_affinity <- graph_affinity / sqrt(outer(rowSums(graph_affinity), rowSums(graph_affinity)))
  diag(graph_affinity) <- -Inf
  graph_neighbour <- max.col(graph_affinity, ties.method = "first")
  direct_support <- quotient[cbind(terminals, terminals[expression_neighbour])] > 0
  joined <- expression_neighbour[expression_neighbour] == template_index &
    graph_neighbour == expression_neighbour & graph_neighbour[expression_neighbour] == template_index &
    direct_support
  acquisition_conflict <- rep(FALSE, template_count)
  switch(as.character(!is.null(program_profiles) & !is.null(lineage_paths)), "TRUE" = {
    profile <- as.matrix(program_profiles)[match(as.character(cl_map$levels),
      rownames(program_profiles)), , drop = FALSE]
    paths <- lapply(lineage_paths, match, table = cl_map$levels)
    nodes <- data.table(state = unlist(paths),
      terminal = rep(vapply(paths, tail, integer(1L), n = 1L), lengths(paths)),
      step = sequence(lengths(paths)))
    shared <- merge(nodes, nodes, by = "state", allow.cartesian = TRUE)
    shared <- shared[terminal.x != terminal.y][order(-step.x),
      .SD[1L], by = .(terminal.x, terminal.y)]
    pairs <- data.table(template = template_index[joined], from = terminals[joined],
      to = terminals[expression_neighbour[joined]])
    pairs[shared, on = .(from = terminal.x, to = terminal.y), anchor := i.state]
    raw <- pmax(sweep(profile, 2L, profile[paths[[1L]][1L], ], "-"), 0)
    acquired <- pmax(profile[pairs$from, , drop = FALSE] -
      profile[pairs$anchor, , drop = FALSE], 0)
    target <- pmax(profile[pairs$to, , drop = FALSE] -
      profile[pairs$anchor, , drop = FALSE], 0)
    acquisition_conflict[pairs$template] <- rowSums(acquired^2) >
      rowSums(pmin(raw[pairs$from, , drop = FALSE], raw[pairs$anchor, , drop = FALSE])^2) &
      rowSums(pmax(acquired - target, 0)^2) > rowSums(pmin(acquired, target)^2)
    joined <- joined & !acquisition_conflict & !acquisition_conflict[expression_neighbour]
  })
  representative <- fifelse(joined, pmin(template_index, expression_neighbour), template_index)
  identity <- match(representative, unique(representative))
  pooling <- Matrix::sparseMatrix(i = template_index, j = identity, x = 1,
    dims = c(template_count, max(identity)))
  pooled <- as.matrix(aggregate %*% pooling) + aggregate_pseudocount
  pooled_composition <- sweep(pooled, 2L, colSums(pooled), "/")
  score <- as.matrix(Matrix::t(counts) %*% log(pooled_composition))
  gene_weight <- rep(1, nrow(counts))
  added_cells <- numeric(max(identity))
  switch(as.character(is.null(template_counts) | !is.null(template_membership)), "TRUE" = {
    anchors <- anchor_membership %*% pooling
    initial_prediction <- max.col(score, ties.method = "first")
    initial_membership <- Matrix::sparseMatrix(i = seq_along(cluster),
      j = initial_prediction, x = 1, dims = c(length(cluster), max(identity)))
    adjacency <- (knn_conn + Matrix::t(knn_conn)) / 2
    diag(adjacency) <- 0
    neighbour_mass <- as.matrix(adjacency %*% initial_membership)
    own_mass <- neighbour_mass[cbind(seq_along(cluster), initial_prediction)]
    eligible <- own_mass > rowSums(neighbour_mass) / 2 & Matrix::rowSums(anchors) == 0
    supplement <- Matrix::Diagonal(x = as.numeric(eligible)) %*% initial_membership
    added_cells <- Matrix::colSums(supplement)
    supplement_counts <- as.matrix(counts %*% supplement)
    supplement_mass <- colSums(supplement_counts)
    supplement_composition <- sweep(supplement_counts, 2L,
      replace(supplement_mass, supplement_mass == 0, 1), "/")
    pooled_composition <- sweep(pooled_composition + supplement_composition, 2L,
      1 + as.numeric(supplement_mass > 0), "/")
    gene_weight <- log1p(ncol(counts) / (1 + Matrix::rowSums(counts > 0)))
    weighted_depth <- as.numeric(Matrix::t(counts) %*% gene_weight)
    score <- as.matrix(Matrix::t(counts) %*% (gene_weight * log(pooled_composition))) -
      outer(weighted_depth, log(colSums(gene_weight * pooled_composition)))
  })
  list(mapping = data.table(template = template_index, state = terminal_states, identity = identity),
    prediction = identity[prediction], rna_prediction = max.col(score, ties.method = "first"),
    composition = pooled_composition, score = score, gene_weight = gene_weight,
    added_cells = added_cells, acquisition_conflict = acquisition_conflict,
    expression_affinity = expression_affinity, graph_affinity = graph_affinity)
}

.persistent_count_identity <- function(counts, cluster_labels, lineage_paths,
                                       prediction, aggregate_pseudocount = 0.5) {
  leaves <- vapply(lineage_paths, tail, integer(1L), n = 1L)
  terminal_code <- match(cluster_labels, leaves)
  anchored <- which(!is.na(terminal_code))
  anchor_y <- Matrix::sparseMatrix(i = anchored, j = terminal_code[anchored], x = 1,
    dims = c(length(cluster_labels), length(leaves)))
  aggregate <- as.matrix(counts %*% anchor_y) + aggregate_pseudocount
  composition <- sweep(aggregate, 2L, colSums(aggregate), "/")
  score <- as.matrix(Matrix::t(counts) %*% log(composition))
  paths <- data.table(state = unlist(lineage_paths),
    lineage = rep(seq_along(leaves), lengths(lineage_paths)))
  membership <- dcast(unique(paths), state ~ lineage, value.var = "lineage",
    fun.aggregate = length)
  private <- rowSums(as.matrix(membership[, -1L])) == 1L &
    !membership$state %in% leaves
  owner <- max.col(as.matrix(membership[, -1L]), ties.method = "first")
  state_score <- data.table(state = cluster_labels, score)[,
    lapply(.SD, stats::median), by = state][match(membership$state, state)]
  state_matrix <- as.matrix(state_score[, -1L])
  other <- state_matrix
  other[cbind(seq_len(nrow(other)), owner)] <- -Inf
  margin <- state_matrix[cbind(seq_len(nrow(other)), owner)] -
    other[cbind(seq_len(nrow(other)), max.col(other, ties.method = "first"))]
  support <- data.table(state = membership$state, lineage = owner,
    private_interior = private, margin = margin)
  stable <- support[private_interior == TRUE,
    .(persistent = all(margin > 0)), by = lineage]
  persistent <- stable$persistent[match(seq_along(leaves), stable$lineage)]
  persistent[is.na(persistent)] <- FALSE
  rna_prediction <- max.col(score, ties.method = "first")
  use_rna <- persistent[rna_prediction] & persistent[prediction]
  list(prediction = fifelse(use_rna, rna_prediction, prediction),
    rna_prediction = rna_prediction, use_rna = use_rna,
    persistent = persistent, state_support = support,
    composition = composition, score = score)
}

.reattach_cross_lineage_tips <- function(program_profiles, induced_program,
                                         cluster_labels, pseudotime,
                                         lineage_paths, lineage_weights,
                                         embedding, knn_conn, ancestor_candidates, sigma = 1,
                                         prediction = max.col(lineage_weights,
                                           ties.method = "first")) {
  cl_map <- .compact_cluster_labels(cluster_labels)
  cluster <- cl_map$labels
  paths <- lapply(lineage_paths, match, table = cl_map$levels)
  root <- paths[[1L]][1L]
  raw <- pmax(sweep(program_profiles, 2L, program_profiles[root, ], "-"), 0)
  phase <- as.numeric(rowsum(pseudotime, cluster) / tabulate(cluster))
  nodes <- data.table(state = unlist(paths),
    lineage = rep(seq_along(paths), lengths(paths)))
  nodes[, parent := shift(state), by = lineage]
  leaves <- vapply(paths, tail, integer(1L), n = 1L)
  adjacency <- (knn_conn + Matrix::t(knn_conn) +
    abs(knn_conn - Matrix::t(knn_conn))) / 2
  membership <- Matrix::sparseMatrix(i = seq_along(cluster), j = cluster, x = 1)
  quotient <- as.matrix(Matrix::t(membership) %*% adjacency %*% membership)
  volume <- rowSums(quotient)
  affinity <- quotient / sqrt(outer(volume, volume))
  pairs <- CJ(from = leaves, to = leaves)[from != to]
  pairs[, parent := nodes$parent[match(from, nodes$state)]]
  pairs[, `:=`(incoming_affinity = affinity[cbind(from, parent)],
    cross_affinity = affinity[cbind(from, to)],
    raw_loss = rowSums(pmax(raw[from, , drop = FALSE] - raw[to, , drop = FALSE], 0)^2),
    raw_retained = rowSums(pmin(raw[from, , drop = FALSE], raw[to, , drop = FALSE])^2),
    program_loss = rowSums(pmax(induced_program[from, , drop = FALSE] -
      induced_program[to, , drop = FALSE], 0)^2),
    program_retained = rowSums(pmin(induced_program[from, , drop = FALSE],
      induced_program[to, , drop = FALSE])^2))]
  best <- pairs[order(from, -cross_affinity, to), .SD[1L], by = from]
  mutual <- merge(best, best[, .(from = to, to = from, target_parent = parent,
    target_incoming = incoming_affinity, reverse_raw_loss = raw_loss,
    reverse_program_loss = program_loss)], by = c("from", "to"))
  selected <- mutual[phase[from] < phase[to] &
    cross_affinity > incoming_affinity & cross_affinity > target_incoming &
    raw_loss <= raw_retained & reverse_raw_loss <= raw_retained &
    program_loss <= program_retained & reverse_program_loss <= program_retained]
  weights <- lineage_weights
  identity_diagnostics <- data.table()
  identity_repairs <- data.table()
  switch(as.character(nrow(selected) > 0L), "TRUE" = {
    edges <- unique(nodes[!is.na(parent), .(from = parent, to = state)])
    edges[selected, on = c(to = "from"), from := i.target_parent]
    edges[selected, on = "to", from := i.from]
    paths <- rs_extract_lineages_from_edges(as.integer(edges$from), as.integer(edges$to),
      as.integer(nrow(program_profiles)), as.integer(root))
    exposed <- intersect(selected$parent, vapply(paths, tail, integer(1L), n = 1L))
    parent <- rep(root, nrow(program_profiles))
    parent[edges$to] <- edges$from
    conserved <- pmin(raw, raw[parent, , drop = FALSE])
    identity_diagnostics <- copy(edges[to %in% exposed])
    identity_diagnostics[, `:=`(
      conserved_energy = rowSums(conserved[from, , drop = FALSE]^2),
      identity_change = rowSums((raw[from, , drop = FALSE] - conserved[from, , drop = FALSE])^2),
      program_loss = rowSums(pmax(conserved[from, , drop = FALSE] - raw[to, , drop = FALSE], 0)^2),
      retained_program = rowSums(pmin(conserved[from, , drop = FALSE], raw[to, , drop = FALSE])^2))]
    lost <- identity_diagnostics[conserved_energy > 0 & identity_change <= conserved_energy &
      program_loss > retained_program & phase[from] > phase[to], to]
    nodes <- data.table(state = unlist(paths),
      lineage = rep(seq_along(paths), lengths(paths)))
    nodes[, step := seq_len(.N), by = lineage]
    ancestry <- merge(nodes[, .(lineage, to = state, to_step = step)],
      nodes[, .(lineage, from = state, from_step = step)], by = "lineage", allow.cartesian = TRUE)
    ancestry <- unique(ancestry[to %in% lost & to_step - from_step > 1L,
      .(from, to, gap = to_step - from_step)])
    alternatives <- copy(as.data.table(ancestor_candidates))
    alternatives[, `:=`(from = match(from, cl_map$levels), to = match(to, cl_map$levels))]
    alternatives <- merge(alternatives, ancestry, by = c("from", "to"))
    identity_repairs <- alternatives[phase_eligible & program_loss <= retained_program][
      order(gap, distance2, from), .SD[1L], by = to]
    edges[identity_repairs, on = "to", from := i.from]
    paths <- rs_extract_lineages_from_edges(as.integer(edges$from), as.integer(edges$to),
      as.integer(nrow(program_profiles)), as.integer(root))
    weights <- ComputeLineageWeightsFromKnnGraph(NULL, embedding, cluster, root,
      sigma = sigma, tree = list(edges = edges))
    prediction <- max.col(weights, ties.method = "first")
  })
  list(lineage_paths = relist(cl_map$levels[unlist(paths)], paths),
    lineage_weights = weights, prediction = prediction, selected = selected,
    identity_diagnostics = identity_diagnostics, identity_repairs = identity_repairs,
    cluster_levels = cl_map$levels)
}

RefineLineageProgramTopology <- function(program_profiles, induced_program,
                                         cluster_labels, pseudotime,
                                         lineage_paths, lineage_weights,
                                         embedding, ancestor_candidates,
                                         cell_counts = NULL, sigma = 1,
                                         prediction = max.col(lineage_weights,
                                           ties.method = "first"), knn_conn = NULL) {
  backward <- ReconnectConservedLineageTips(program_profiles, induced_program,
    cluster_labels, pseudotime, lineage_paths, lineage_weights,
    direction = "backward", prediction = prediction)
  identity <- RecoverPersistentLineageIdentity(program_profiles, cluster_labels,
    pseudotime, backward$lineage_paths, backward$lineage_weights, embedding,
    ancestor_candidates, sigma = sigma, prediction = backward$prediction)
  forward <- ReconnectConservedLineageTips(program_profiles, induced_program,
    cluster_labels, pseudotime, identity$lineage_paths, identity$lineage_weights,
    direction = "forward", prediction = identity$prediction,
    identity_terminals = identity$new_terminals)
  plateau <- ReconnectConservedLineageTips(program_profiles, induced_program,
    cluster_labels, pseudotime, forward$lineage_paths, forward$lineage_weights,
    direction = "plateau", prediction = forward$prediction)
  cross_tip <- switch(as.character(is.null(knn_conn)),
    "TRUE" = plateau,
    "FALSE" = .reattach_cross_lineage_tips(program_profiles, induced_program,
      cluster_labels, pseudotime, plateau$lineage_paths, plateau$lineage_weights,
      embedding, knn_conn, ancestor_candidates, sigma = sigma, prediction = plateau$prediction))
  result <- ResolvePersistentLineageMembership(program_profiles, cluster_labels,
    pseudotime, cross_tip$lineage_paths, cross_tip$lineage_weights,
    prediction = cross_tip$prediction, cell_counts = cell_counts)
  result$program_stages <- list(backward = backward, identity = identity,
    forward = forward, plateau = plateau, cross_tip = cross_tip)
  result
}

ComputeForkContrastEmbedding <- function(embedding, cluster_labels, lineage_paths,
                                        include_state_changes = FALSE) {
  emb <- .prepare_embedding(embedding)
  stopifnot(length(cluster_labels) == nrow(emb))
  cl_map <- .compact_cluster_labels(cluster_labels)
  paths <- lapply(lineage_paths, match, table = cl_map$levels)
  stopifnot(length(paths) > 0L, all(lengths(paths) > 0L),
            all(!is.na(unlist(paths))))
  centers <- rowsum(emb, cl_map$labels) / tabulate(cl_map$labels)
  nodes <- data.table::rbindlist(lapply(paths, data.table::as.data.table), idcol = "path")
  nodes[, next_node := data.table::shift(V1, type = "lead"), by = path]
  edges <- unique(nodes[!is.na(next_node), .(parent = V1, child = next_node)])
  pairs <- merge(edges, edges, by = "parent", allow.cartesian = TRUE)[child.x < child.y]
  contrasts <- switch(as.character(include_state_changes),
    "FALSE" = centers[pairs$child.x, , drop = FALSE] -
      centers[pairs$child.y, , drop = FALSE],
    "TRUE" = centers[edges$child, , drop = FALSE] -
      centers[edges$parent, , drop = FALSE])
  basis <- matrix(numeric(), nrow = ncol(emb), ncol = 0L)
  contrast_rank <- 0L
  if (nrow(contrasts) > 0L) {
    decomposition <- svd(t(contrasts))
    contrast_rank <- sum(decomposition$d > max(dim(contrasts)) *
      max(decomposition$d) * .Machine$double.eps)
    basis <- decomposition$u[, seq_len(contrast_rank), drop = FALSE]
  }
  list(embedding = emb %*% basis, basis = basis, rank = contrast_rank,
    include_state_changes = include_state_changes,
    state_edges = data.frame(parent = cl_map$levels[edges$parent],
      child = cl_map$levels[edges$child]),
    fork_pairs = data.frame(parent = cl_map$levels[pairs$parent],
      child_left = cl_map$levels[pairs$child.x],
      child_right = cl_map$levels[pairs$child.y]),
    cluster_centers = centers, cluster_levels = cl_map$levels)
}

#' Cluster-level tree from the native cell kNN graph.
#'
#' Aggregates the observed cell-level kNN connectivity into cluster-to-cluster
#' connection mass `C_ab` and normalised support `S_ab = C_ab / sqrt(V_a V_b)`,
#' where `V_a` is the total connection volume of cluster `a`. A deterministic
#' maximum spanning tree is built on the positive-support edges of each
#' connected component; no edge is invented for clusters without observed
#' connections. No cell is dropped.
#'
#' @param knn_conn Cell-by-cell connectivity matrix (symmetric sparse, e.g.
#'   from `knn_conn.mtx`).
#' @param cluster_labels Integer vector of length n_cells, 1-based.
#' @param start_cluster Integer scalar (1-based root cluster); used to report
#'   whether the root component spans the tree.
#' @return list with `edges` (data.frame `from`, `to`, `support`), `components`
#'   (integer per cluster), `n_components`, `s_matrix`, `root_component`.
#' @export
BuildKnnClusterTree <- function(knn_conn, cluster_labels, start_cluster) {
  if (!inherits(knn_conn, "CsparseMatrix")) {
    knn_conn <- as(Matrix::Matrix(knn_conn, sparse = TRUE), "CsparseMatrix")
  }
  knn_conn_t <- Matrix::t(knn_conn)
  knn_conn <- (knn_conn + knn_conn_t + abs(knn_conn - knn_conn_t)) / 2
  n_cells <- nrow(knn_conn)
  if (length(cluster_labels) != n_cells) {
    stop("`cluster_labels` length must equal nrow(knn_conn)", call. = FALSE)
  }
  cl_map <- .compact_cluster_labels(cluster_labels)
  cl <- cl_map$labels
  n_clusters <- length(cl_map$levels)
  start_cluster <- .compact_cluster_id(start_cluster, cl_map$levels)
  if (is.na(start_cluster)) stop("`start_cluster` is not present in `cluster_labels`")

  # C = M^T %*% knn %*% M with sparse cluster indicator M (cells x clusters).
  m <- Matrix::sparseMatrix(
    i = seq_len(n_cells),
    j = cl,
    x = 1,
    dims = c(n_cells, n_clusters)
  )
  c_mat <- as.matrix(Matrix::t(m) %*% knn_conn %*% m)
  v <- rowSums(c_mat)
  denom <- sqrt(pmax(v, 0) %o% pmax(v, 0))
  s_mat <- c_mat / pmax(denom, .Machine$double.eps)
  diag(s_mat) <- 0

  edges <- which(s_mat > 0 & upper.tri(s_mat), arr.ind = TRUE)
  if (nrow(edges)) {
    edges <- data.table::data.table(
      from = as.integer(edges[, 1L]),
      to = as.integer(edges[, 2L]),
      support = s_mat[edges]
    )
    data.table::setorder(edges, -support, from, to)
  } else {
    edges <- data.table::data.table(
      from = integer(), to = integer(), support = numeric()
    )
  }

  # Kruskal maximum spanning tree with union-find (deterministic order).
  parent <- seq_len(n_clusters)
  find_root <- function(x) {
    while (parent[x] != x) {
      parent[x] <<- parent[parent[x]]
      x <- parent[x]
    }
    x
  }
  keep <- rep(FALSE, nrow(edges))
  if (nrow(edges)) {
    for (e in seq_len(nrow(edges))) {
      ra <- find_root(edges$from[e])
      rb <- find_root(edges$to[e])
      if (ra != rb) {
        parent[ra] <- rb
        keep[e] <- TRUE
      }
    }
  }
  tree <- edges[keep]
  comp <- vapply(seq_len(n_clusters), function(i) find_root(i), integer(1L))
  comp <- as.integer(match(comp, unique(comp)))

  list(
    edges = tree,
    components = comp,
    n_components = data.table::uniqueN(comp),
    s_matrix = s_mat,
    root_component = comp[start_cluster]
  )
}

#' Root-to-leaf lineage paths from the native cell kNN cluster tree.
#'
#' @inheritParams BuildKnnClusterTree
#' @return list of integer vectors (1-based cluster ids, root first), from the
#'   tree of the root cluster's connected component.
#' @export
ExtractLineagesFromKnnGraph <- function(knn_conn, cluster_labels, start_cluster) {
  cl_map <- .compact_cluster_labels(cluster_labels)
  start_cluster <- .compact_cluster_id(start_cluster, cl_map$levels)
  tree <- BuildKnnClusterTree(knn_conn, cl_map$labels, start_cluster)
  paths <- rs_extract_lineages_from_edges(
    as.integer(tree$edges$from),
    as.integer(tree$edges$to),
    as.integer(length(cl_map$levels)),
    as.integer(start_cluster)
  )
  list(
    paths = lapply(paths, function(path) as.integer(cl_map$levels[as.integer(path)])),
    tree = tree
  )
}

#' Trunk-aware lineage weights over the native cell kNN cluster tree.
#'
#' Same weighting arithmetic as `ComputeLineageWeightsTrunkAware`, but the
#' paths come from the kNN cluster tree so weights, paths and fork placement
#' share one observed-support structure.
#'
#' @inheritParams BuildKnnClusterTree
#' @param embedding numeric matrix \[n_cells × n_dims\].
#' @param sigma numeric scalar, Gaussian bandwidth on lineage-path distance.
#' @return numeric matrix \[n_cells × n_lineages\], rows sum to 1.
#' @export
ComputeLineageWeightsFromKnnGraph <- function(knn_conn, embedding, cluster_labels,
                                              start_cluster, sigma = 1.0,
                                              tree = NULL) {
  emb <- .prepare_embedding(embedding)
  sigma <- .validate_positive_num(sigma, "sigma")
  if (length(cluster_labels) != nrow(emb)) {
    stop("`cluster_labels` length must equal nrow(embedding)", call. = FALSE)
  }
  cl_map <- .compact_cluster_labels(cluster_labels)
  cl <- cl_map$labels
  n_clusters <- length(cl_map$levels)
  start_cluster <- .compact_cluster_id(start_cluster, cl_map$levels)
  if (is.na(start_cluster)) stop("`start_cluster` is not present in `cluster_labels`")
  if (is.null(tree)) tree <- BuildKnnClusterTree(knn_conn, cluster_labels, start_cluster)
  flat <- rs_compute_lineage_weights_trunk_aware_from_edges(
    as.numeric(t(emb)),
    cl,
    as.integer(nrow(emb)),
    as.integer(ncol(emb)),
    as.integer(n_clusters),
    as.integer(tree$edges$from),
    as.integer(tree$edges$to),
    as.integer(start_cluster),
    as.numeric(sigma)
  )
  n_cells <- nrow(emb)
  n_lineages <- length(flat) / n_cells
  weights <- matrix(flat, nrow = n_cells, ncol = n_lineages, byrow = TRUE)
  attr(weights, "knn_tree") <- tree
  weights
}

#' Binned LR scoring with intermediate tensors for entropy computation.
#'
#' Same scoring computation as `ComputeTemporalLR`, plus the
#' `LrBinIntermediate` (`hwsums` / `wsums` tensors) that per-lineage entropy needs to
#' marginalise `sum_{ct}` and build per-lineage entropy `H_l(b)`.  Call
#' `ComputeLineageEntropyRaw()` on the returned `intermediate` to get the
#' raw entropy matrix; pipe that through `SmoothMA3()` for `H̃_l(b)`.
#'
#' @param expr,celltypes,lr_db,pseudotime,lineage_weights,n_bins,min_pct,min_n_eff
#'   see `ComputeTemporalLR`.
#' @param return_complete_grid see `ComputeTemporalLR`.
#' @param requested_contexts see `ComputeTemporalLR`.
#' @param spatial_adjacency,min_spatial_support passed to `ComputeTemporalLR()`.
#' @param common_grid see `ComputeTemporalLR`.
#' @return list with two named elements:
#'   \itemize{
#'     \item `scores`: data.frame (same columns as `ComputeTemporalLR`)
#'     \item `intermediate`: list(`hwsums`, `wsums`, `n_types`,
#'       `n_lineages`, `n_bins`) — do NOT inspect directly; pass to
#'       `ComputeLineageEntropyRaw()`.
#'   }
#' @export
ComputeTemporalLRWithIntermediate <- function(expr, celltypes, lr_db, pseudotime,
                                              lineage_weights,
                                              n_bins = 12L, min_pct = 0.05,
                                              min_n_eff = 5.0,
                                              return_complete_grid = TRUE,
                                              requested_contexts = NULL,
                                              spatial_adjacency = NULL,
                                              min_spatial_support = 0.05,
                                              common_grid = FALSE) {
  n_bins <- .validate_positive_int(n_bins, "n_bins", min = 2L)
  min_pct <- .validate_unit_interval(min_pct, "min_pct")
  min_n_eff <- .validate_positive_num(min_n_eff, "min_n_eff")
  if (!is.logical(common_grid) || length(common_grid) != 1L ||
      is.na(common_grid)) {
    stop("`common_grid` must be TRUE or FALSE", call. = FALSE)
  }
  expr <- .prepare_lr_expr(expr)
  lr_db <- .validate_lr_db(lr_db)
  .assert_lr_db_indices(lr_db, nrow(expr))
  lw <- .assert_lineage_weights(lineage_weights, ncol(expr))
  n_lineages <- ncol(lw)
  ct      <- .assert_celltypes(celltypes, ncol(expr))
  ct_levels <- attr(ct, "levels")
  n_types <- max(ct)
  pt_vec  <- .assert_pseudotime_1d(pseudotime, ncol(expr))
  req_ctx <- if (is.null(requested_contexts)) {
    list(pair_idx = integer(), source = integer(), target = integer())
  } else {
    .pooled_lr_context_vectors(requested_contexts, ct_levels)
  }
  raw <- rs_compute_temporal_lr_with_inter(
    expr,
    ct,
    as.integer(n_types),
    lapply(lr_db$ligand_rows,   as.integer),
    lapply(lr_db$receptor_rows, as.integer),
    pt_vec,
    as.numeric(t(lw)),
    as.integer(n_lineages),
    n_bins,
    min_pct,
    min_n_eff,
    isTRUE(return_complete_grid),
    req_ctx$pair_idx,
    req_ctx$source,
    req_ctx$target,
    common_grid
  )
  scores <- data.frame(raw$scores, stringsAsFactors = FALSE)
  scores$structural_valid <- as.logical(scores$structural_valid)
  scores$expression_supported <- as.logical(scores$expression_supported)
  scores$source <- ct_levels[scores$source]
  scores$target <- ct_levels[scores$target]
  attr(scores, "celltype_levels") <- ct_levels
  attr(scores, "temporal_coordinate") <- if (common_grid) {
    "common_pseudotime_grid"
  } else {
    "lineage_local_equal_frequency_bins"
  }
  scores <- .attach_spatial_columns(
    scores, spatial_adjacency, min_spatial_support
  )
  intermediate <- raw$intermediate
  intermediate$celltype_levels <- ct_levels
  intermediate$hill_scale <- as.numeric(raw$hill_scale)
  attr(scores, "hill_scale") <- as.numeric(raw$hill_scale)
  bounds_matrix <- matrix(
    as.numeric(intermediate$bin_bounds),
    nrow = n_lineages, byrow = TRUE
  )
  pooled_bounds <- lapply(seq_len(n_lineages), function(i) {
    bounds_matrix[i, ]
  })
  intermediate$bin_bounds <- pooled_bounds
  attr(scores, "bin_bounds") <- pooled_bounds
  list(scores = scores, intermediate = intermediate)
}

#' Temporal LR scoring with pooled bin boundaries and Hill scale.
#'
#' This sample-rescoring API keeps all samples on boundaries and saturation
#' scale estimated from the pooled data. It does not estimate a new temporal
#' coordinate system within each sample.
#'
#' @param expr,celltypes,lr_db,pseudotime,lineage_weights,min_pct,min_n_eff
#'   see `ComputeTemporalLR()`.
#' @param bin_bounds List or matrix with one row per lineage and `n_bins + 1`
#'   finite, non-decreasing boundaries from 0 to 1.
#' @param hill_scale Positive pooled Hill saturation scale.
#' @param n_bins Number of temporal bins. Inferred from `bin_bounds` when NULL.
#' @param return_complete_grid,requested_contexts see `ComputeTemporalLR()`.
#' @param spatial_adjacency,min_spatial_support passed to the spatial evidence
#'   attachment layer.
#' @param shared_celltypes Explicit cell-type labels shared across lineage
#'   queries within this scoring context. Their own lineage probabilities do
#'   not restrict their availability to another population's branch.
#' @param normalize_lineage_weights Normalize each cell across lineages. Set FALSE
#'   for overlapping observed-subtree membership indicators.
#' @return Temporal LR score data frame on the supplied coordinate system.
#' @export
ComputeTemporalLRFixedBounds <- function(
    expr, celltypes, lr_db, pseudotime, lineage_weights,
    bin_bounds, hill_scale, n_bins = NULL, min_pct = 0.05,
    min_n_eff = 5.0, return_complete_grid = TRUE,
    requested_contexts = NULL, spatial_adjacency = NULL,
    min_spatial_support = 0.05, shared_celltypes = NULL,
    normalize_lineage_weights = TRUE) {
  expr <- .prepare_lr_expr(expr)
  lr_db <- .validate_lr_db(lr_db)
  .assert_lr_db_indices(lr_db, nrow(expr))
  lw <- .assert_lineage_weights(lineage_weights, ncol(expr))
  switch(as.character(normalize_lineage_weights), "FALSE" = {
    lw <- as.matrix(lineage_weights)
  })
  n_lineages <- ncol(lw)
  bounds <- if (is.matrix(bin_bounds)) {
    lapply(seq_len(nrow(bin_bounds)), function(i) {
      as.numeric(bin_bounds[i, ])
    })
  } else {
    bin_bounds
  }
  if (!is.list(bounds) || length(bounds) != n_lineages ||
      !all(vapply(bounds, is.numeric, logical(1L)))) {
    stop("`bin_bounds` must contain one numeric vector per lineage",
         call. = FALSE)
  }
  inferred_n_bins <- unique(vapply(bounds, length, integer(1L))) - 1L
  if (length(inferred_n_bins) != 1L || inferred_n_bins < 2L) {
    stop("all `bin_bounds` rows must share length >= 3", call. = FALSE)
  }
  n_bins <- if (is.null(n_bins)) inferred_n_bins else
    .validate_positive_int(n_bins, "n_bins", min = 2L)
  if (n_bins != inferred_n_bins) {
    stop("`n_bins` does not match `bin_bounds`", call. = FALSE)
  }
  valid_bounds <- vapply(bounds, function(x) {
    all(is.finite(x)) && all(diff(x) >= 0) &&
      abs(x[[1L]]) <= 1e-12 && abs(x[[length(x)]] - 1) <= 1e-12
  }, logical(1L))
  if (!all(valid_bounds)) {
    stop("each boundary row must be finite, non-decreasing, and span 0 to 1",
         call. = FALSE)
  }
  hill_scale <- as.numeric(hill_scale)
  if (length(hill_scale) != 1L || !is.finite(hill_scale) ||
      hill_scale <= 0) {
    stop("`hill_scale` must be one positive finite value", call. = FALSE)
  }
  min_pct <- .validate_unit_interval(min_pct, "min_pct")
  min_n_eff <- .validate_positive_num(min_n_eff, "min_n_eff")
  ct <- .assert_celltypes(celltypes, ncol(expr))
  ct_levels <- attr(ct, "levels")
  n_types <- max(ct)
  shared_codes <- match(shared_celltypes, ct_levels)
  stopifnot(!anyNA(shared_codes))
  # A declared context-wide population is observed in each lineage query.
  # Its cells are counted once per query, independently of recipient ancestry.
  lw[ct %in% shared_codes, ] <- 1
  pt_vec <- .assert_pseudotime_1d(pseudotime, ncol(expr))
  req_ctx <- if (is.null(requested_contexts)) {
    list(pair_idx = integer(), source = integer(), target = integer())
  } else {
    .pooled_lr_context_vectors(requested_contexts, ct_levels)
  }
  raw <- rs_compute_temporal_lr_fixed_bounds(
    expr,
    ct,
    as.integer(n_types),
    lapply(lr_db$ligand_rows, as.integer),
    lapply(lr_db$receptor_rows, as.integer),
    pt_vec,
    as.numeric(t(lw)),
    as.integer(n_lineages),
    as.integer(n_bins),
    as.numeric(t(do.call(rbind, bounds))),
    hill_scale,
    min_pct,
    min_n_eff,
    isTRUE(return_complete_grid),
    req_ctx$pair_idx,
    req_ctx$source,
    req_ctx$target
  )
  scores <- data.frame(raw, stringsAsFactors = FALSE)
  scores$structural_valid <- as.logical(scores$structural_valid)
  scores$expression_supported <- as.logical(scores$expression_supported)
  scores$source <- ct_levels[scores$source]
  scores$target <- ct_levels[scores$target]
  attr(scores, "celltype_levels") <- ct_levels
  attr(scores, "hill_scale") <- hill_scale
  attr(scores, "bin_bounds") <- bounds
  if (length(shared_codes)) attr(scores, "shared_celltypes") <- ct_levels[shared_codes]
  attr(scores, "temporal_coordinate") <- if (
    length(bounds) <= 1L ||
      all(vapply(bounds[-1L], function(row) {
        isTRUE(all.equal(row, bounds[[1L]], tolerance = 1e-12))
      }, logical(1L)))
  ) {
    "common_pseudotime_grid"
  } else {
    "lineage_local_equal_frequency_bins"
  }
  .attach_spatial_columns(scores, spatial_adjacency, min_spatial_support)
}

#' Temporal permutation null for LR trajectory dynamics.
#'
#' Permutes pseudotime within each `cell type x dominant lineage` stratum while
#' fixing expression, cell types, lineage weights, and the observed Hill
#' saturation scale.  The trajectory table is the inferential layer; the bin
#' table localises dynamic intervals using a max-statistic within each
#' trajectory.
#'
#' @param expr,celltypes,lr_db,pseudotime,lineage_weights,n_bins,min_pct
#'   see `ComputeTemporalLR`.
#' @param min_n_eff minimum Kish effective sample size required for a
#'   structurally valid bin.
#' @param min_valid_bins minimum number of structurally valid bins required
#'   for a trajectory.
#' @param min_valid_fraction minimum fraction of structurally valid bins
#'   required for a trajectory.
#' @param n_perms number of temporal permutations.
#' @param seed integer seed.
#' @param stat trajectory statistic: `"amplitude"` (temporal variation),
#'   `"max"` / `"mean_score"` / `"topk_mean"` (trajectory strength; not
#'   temporal variation), or deprecated alias `"auc"` (= `"mean_score"`).
#' @param top_k number of bins used by `"topk_mean"`.
#' @param min_perm_observed minimum observed valid permutations per trajectory.
#' @param min_perm_coverage minimum `n_perm_observed / n_perm_expected`.
#' @param block_labels optional donor, sample, or batch labels. Pseudotime is
#'   permuted only within the chosen `permutation_stratum` (see below).
#' @param requested_trajectories optional data.frame with columns `pair_idx`,
#'   `source`, `target`, and `lineage` (1-based indices). When supplied,
#'   permutation evaluates this pooled hypothesis family instead of locally
#'   prescreened trajectories.
#' @param permutation_stratum `"block_celltype_lineage"` (default) shuffles
#'   within `block × cell type × dominant lineage`; `"block_celltype"` shuffles
#'   within `block × cell type` only and is appropriate for soft lineage weights.
#' @return list with `trajectories` and `bins` data.frames. Both expose
#'   `analysis_data_hash`, `scoring_config_hash` (alias `calibration_config_hash`
#'   from the Rust scoring layer), `temporal_null_config_hash`, and dataset-level
#'   `hill_scale` used for compatibility checks; the trajectory table also names
#'   the tested hypothesis and provides statistic-specific p/q aliases.
#' @export
ComputeTemporalPermutation <- function(expr, celltypes, lr_db, pseudotime,
                                       lineage_weights,
                                       n_bins = 12L, min_pct = 0.05,
                                       min_n_eff = 5.0,
                                       min_valid_bins = 3L,
                                       min_valid_fraction = 0.25,
                                       n_perms = 999L, seed = 12345L,
                                       stat = c("amplitude", "max", "mean_score", "auc", "topk_mean"),
                                       top_k = 3L,
                                       min_perm_observed = 100L,
                                       min_perm_coverage = 0.8,
                                       block_labels = NULL,
                                       requested_trajectories = NULL,
                                       permutation_stratum = c(
                                         "block_celltype_lineage", "block_celltype"
                                       )) {
  stat_info <- .normalize_trajectory_stat(stat)
  stat_code <- stat_info$code
  permutation_stratum <- match.arg(permutation_stratum)
  stratum_code <- match(permutation_stratum, c(
    "block_celltype_lineage", "block_celltype"
  )) - 1L
  n_bins <- .validate_positive_int(n_bins, "n_bins", min = 2L)
  n_perms <- .validate_positive_int(n_perms, "n_perms", min = 1L)
  top_k <- .validate_positive_int(top_k, "top_k", min = 1L)
  min_valid_bins <- .validate_positive_int(min_valid_bins, "min_valid_bins", min = 2L)
  min_n_eff <- .validate_positive_num(min_n_eff, "min_n_eff")
  min_pct <- .validate_unit_interval(min_pct, "min_pct")
  min_valid_fraction <- .validate_unit_interval(min_valid_fraction, "min_valid_fraction")
  min_perm_observed <- .validate_non_negative_int(min_perm_observed, "min_perm_observed")
  min_perm_coverage <- .validate_unit_interval(min_perm_coverage, "min_perm_coverage")
  seed <- .validate_seed(seed)
  expr <- .prepare_lr_expr(expr)
  lr_db <- .validate_lr_db(lr_db)
  .assert_lr_db_indices(lr_db, nrow(expr))
  lw <- .assert_lineage_weights(lineage_weights, ncol(expr))
  n_lineages <- ncol(lw)
  ct <- .assert_celltypes(celltypes, ncol(expr))
  ct_levels <- attr(ct, "levels")
  n_types <- max(ct)
  pt_vec <- .assert_pseudotime_1d(pseudotime, ncol(expr))
  blocks <- .prepare_block_labels(block_labels, ncol(expr))
  block_levels <- attr(blocks, "levels")
  traj_by <- c("pair_idx", "source", "target", "lineage")
  req_traj <- if (is.null(requested_trajectories)) {
    stats::setNames(rep(list(integer()), length(traj_by)), traj_by)
  } else {
    .pooled_trajectory_key_vectors(requested_trajectories, traj_by, ct_levels)
  }
  raw <- rs_compute_temporal_permutation(
    expr,
    ct,
    as.integer(n_types),
    lapply(lr_db$ligand_rows, as.integer),
    lapply(lr_db$receptor_rows, as.integer),
    pt_vec,
    as.numeric(t(lw)),
    blocks,
    as.integer(n_lineages),
    n_bins,
    min_pct,
    min_n_eff,
    min_valid_bins,
    min_valid_fraction,
    n_perms,
    as.numeric(seed),
    as.integer(stat_code),
    top_k,
    as.integer(stratum_code),
    req_traj$pair_idx,
    req_traj$source,
    req_traj$target,
    req_traj$lineage
  )
  trajectories <- data.frame(raw$trajectories, stringsAsFactors = FALSE)
  trajectories$source <- ct_levels[trajectories$source]
  trajectories$target <- ct_levels[trajectories$target]
  insufficient <- trajectories$n_perm_observed < min_perm_observed |
    trajectories$perm_coverage < min_perm_coverage
  trajectories$p_temporal[insufficient] <- NA_real_
  trajectories$q_temporal <- stats::p.adjust(trajectories$p_temporal, method = "BH")
  temporal_hypothesis <- ifelse(
    identical(stat_info$stat, "amplitude"),
    "temporal_variation",
    "trajectory_strength"
  )
  p_alias <- ifelse(
    identical(stat_info$stat, "amplitude"),
    "p_temporal_variation",
    "p_trajectory_strength"
  )
  q_alias <- ifelse(
    identical(stat_info$stat, "amplitude"),
    "q_temporal_variation",
    "q_trajectory_strength"
  )
  trajectories[[p_alias]] <- trajectories$p_temporal
  trajectories[[q_alias]] <- trajectories$q_temporal
  trajectories$trajectory_stat_name <- stat_info$stat
  trajectories$calibration_hypothesis <- temporal_hypothesis
  provenance <- list(
    analysis_data_hash = .compute_analysis_data_hash(
      expr, ct, lr_db, pt_vec, lw, block_labels
    ),
    scoring_config_hash = .compute_scoring_config_hash(
      n_bins, min_pct, min_n_eff, min_valid_bins, min_valid_fraction,
      stat_code, top_k
    ),
    temporal_null_config_hash = .compute_temporal_null_config_hash(
      n_perms, stat_info$stat, block_labels, permutation_stratum, seed,
      ncol(expr)
    )
  )
  trajectories <- .attach_provenance_columns(trajectories, provenance)
  trajectories$calibration_config_hash <- trajectories$scoring_config_hash
  trajectories$hill_scale <- as.numeric(raw$hill_scale)
  trajectories$permutation_stratum <- permutation_stratum
  trajectories$permutation_block_mode <- if (permutation_stratum == "block_celltype") {
    if (is.null(block_labels)) "celltype" else "block_x_celltype"
  } else if (is.null(block_labels)) {
    "dominant_lineage_within_celltype"
  } else {
    "block_x_celltype_x_dominant_lineage"
  }
  trajectories$n_permutation_blocks <- length(block_levels)
  bins <- data.frame(raw$bins, stringsAsFactors = FALSE)
  bins$source <- ct_levels[bins$source]
  bins$target <- ct_levels[bins$target]
  bins <- .attach_provenance_columns(bins, provenance)
  bins$calibration_config_hash <- bins$scoring_config_hash
  bins$hill_scale <- as.numeric(raw$hill_scale)
  bins$permutation_stratum <- permutation_stratum
  bins$permutation_block_mode <- trajectories$permutation_block_mode[1]
  bins$n_permutation_blocks <- length(block_levels)
  bins$p_temporal_bin_maxT <- bins$p_temporal_bin
  bins$localization_method <- "within_trajectory_maxT"
  trajectories$hypothesis_family <- "prescreened_evaluable_trajectories"
  if (!identical(stat_info$stat, "amplitude")) {
    bins$p_temporal_bin <- NA_real_
    bins$p_temporal_bin_maxT <- NA_real_
    attr(bins, "localization_formal") <- FALSE
  } else {
    attr(bins, "localization_formal") <- TRUE
  }
  attr(trajectories, "celltype_levels") <- ct_levels
  attr(bins, "celltype_levels") <- ct_levels
  permutation_block_audit <- .audit_permutation_blocks(
    blocks, ct, lw, permutation_stratum
  )
  if (!is.null(block_labels) && permutation_block_audit$singleton_block_fraction > 0.5) {
    warning(
      "More than half of `block_labels` are singletons; ",
      "donor-level permutation within block x celltype x lineage may have low power.",
      call. = FALSE
    )
  }
  list(
    trajectories = trajectories,
    bins = bins,
    permutation_block_audit = permutation_block_audit,
    analysis_data_hash = provenance$analysis_data_hash,
    scoring_config_hash = provenance$scoring_config_hash,
    calibration_config_hash = provenance$scoring_config_hash,
    temporal_null_config_hash = provenance$temporal_null_config_hash,
    hill_scale = as.numeric(raw$hill_scale)
  )
}

.subset_lr_db <- function(lr_db, pair_idx) {
  pair_idx <- unique(as.integer(pair_idx))
  pair_idx <- pair_idx[pair_idx >= 1L & pair_idx <= length(lr_db$ligand_rows)]
  if (!length(pair_idx)) {
    stop("candidate subset contains no valid LR pair indices", call. = FALSE)
  }
  list(
    ligand_rows = lr_db$ligand_rows[pair_idx],
    receptor_rows = lr_db$receptor_rows[pair_idx]
  )
}

.remap_lr_db_pair_idx <- function(df, pair_idx_map) {
  if (!"pair_idx" %in% colnames(df) || !nrow(df)) {
    return(df)
  }
  local_idx <- seq_along(pair_idx_map)
  observed <- unique(as.integer(df$pair_idx))
  if (all(observed %in% local_idx)) {
    mapped <- pair_idx_map[as.character(df$pair_idx)]
    if (any(is.na(mapped))) {
      stop("confirm-stage LR subset returned unexpected pair_idx values", call. = FALSE)
    }
    df$pair_idx <- as.integer(mapped)
  }
  df
}

.adaptive_perm_select_candidates <- function(traj, p_col = "p_temporal",
                                              threshold = 0.1, top_k = 20L,
                                              by = c("pair_idx", "source", "target", "lineage")) {
  if (!p_col %in% colnames(traj)) {
    stop(sprintf("candidate table must contain `%s`", p_col), call. = FALSE)
  }
  p <- as.numeric(traj[[p_col]])
  ok <- is.finite(p)
  if (!any(ok)) {
    return(character(0))
  }
  keys <- .trajectory_key_col(traj, by)
  by_threshold <- ok & p <= threshold
  top_k <- .validate_non_negative_int(top_k, "confirm_top_k")
  by_top <- rep(FALSE, length(keys))
  if (top_k > 0L) {
    ord <- order(p[ok], na.last = NA)
    by_top[which(ok)[ord[seq_len(min(top_k, length(ord)))]]] <- TRUE
  }
  unique(keys[by_threshold | by_top])
}

.adaptive_perm_diagnostics <- function(p, n_perm_used) {
  n <- as.numeric(n_perm_used)
  list(
    min_achievable_p = ifelse(is.finite(n) & n >= 0, 1 / (n + 1), NA_real_),
    monte_carlo_se = ifelse(
      is.finite(p) & is.finite(n) & n > 0,
      sqrt(p * (1 - p) / n),
      NA_real_
    )
  )
}

.adaptive_perm_attach_diagnostics <- function(traj, p_col, n_perm_used, perm_stage) {
  p <- as.numeric(traj[[p_col]])
  n <- rep(as.numeric(n_perm_used), length.out = nrow(traj))
  diag <- .adaptive_perm_diagnostics(p, n)
  traj$perm_stage <- perm_stage
  traj$n_perm_used <- as.integer(n)
  traj$min_achievable_p <- diag$min_achievable_p
  traj$monte_carlo_se <- diag$monte_carlo_se
  traj
}

.merge_adaptive_temporal_perm <- function(screen_out, confirm_out,
                                          candidate_keys,
                                          n_perms_screen, n_perms_confirm,
                                          p_col = "p_temporal",
                                          by = c("pair_idx", "source", "target", "lineage")) {
  traj <- screen_out$trajectories
  bins <- screen_out$bins
  traj[[paste0(p_col, "_screen")]] <- traj[[p_col]]
  q_col <- sub("^p_", "q_", p_col)
  if (q_col %in% colnames(traj)) {
    traj[[paste0(q_col, "_screen")]] <- traj[[q_col]]
  }
  screen_observed <- if ("n_perm_observed" %in% colnames(traj)) {
    traj$n_perm_observed
  } else {
    n_perms_screen
  }
  traj <- .adaptive_perm_attach_diagnostics(
    traj, p_col, screen_observed, "screen"
  )
  keys <- .trajectory_key_col(traj, by)
  is_candidate <- keys %in% candidate_keys
  traj$adaptive_perm_status <- ifelse(is_candidate, "pending_confirm", "screened_out")

  if (length(candidate_keys) > 0L && !is.null(confirm_out)) {
    conf_traj <- confirm_out$trajectories
    conf_keys <- .trajectory_key_col(conf_traj, by)
    hit <- match(keys[is_candidate], conf_keys)
    idx <- which(is_candidate)
    ok <- !is.na(hit)
    if (any(ok)) {
      replace_cols <- setdiff(
        intersect(colnames(conf_traj), colnames(traj)),
        by
      )
      for (col in replace_cols) {
        traj[[col]][idx[ok]] <- conf_traj[[col]][hit[ok]]
      }
      traj$perm_stage[idx[ok]] <- "confirm"
      confirm_observed <- if ("n_perm_observed" %in% colnames(conf_traj)) {
        conf_traj$n_perm_observed[hit[ok]]
      } else {
        rep(n_perms_confirm, sum(ok))
      }
      traj$n_perm_used[idx[ok]] <- as.integer(confirm_observed)
      traj$adaptive_perm_status[idx[ok]] <- "confirmed"
      diag <- .adaptive_perm_diagnostics(
        traj[[p_col]][idx[ok]], confirm_observed
      )
      traj$min_achievable_p[idx[ok]] <- diag$min_achievable_p
      traj$monte_carlo_se[idx[ok]] <- diag$monte_carlo_se
    }

    conf_bins <- confirm_out$bins
    if (nrow(conf_bins) > 0L) {
      bin_by <- c(by, "bin")
      screen_bin_key <- do.call(
        paste,
        c(bins[, bin_by, drop = FALSE], sep = "\002")
      )
      conf_bin_key <- do.call(
        paste,
        c(conf_bins[, bin_by, drop = FALSE], sep = "\002")
      )
      bin_traj_keys <- .trajectory_key_col(bins, by)
      is_cand_bin <- bin_traj_keys %in% candidate_keys
      conf_hit <- match(screen_bin_key, conf_bin_key)
      update_bins <- is_cand_bin & !is.na(conf_hit)
      bin_replace_cols <- intersect(
        colnames(conf_bins),
        c("p_temporal_bin", "p_temporal_bin_maxT", "temporal_deviation",
          "score", "deviation", "structural_valid", "expression_supported")
      )
      for (col in bin_replace_cols) {
        bins[[col]][update_bins] <- conf_bins[[col]][conf_hit[update_bins]]
      }
    }
  }

  traj$adaptive_perm_status[traj$adaptive_perm_status == "pending_confirm"] <- "screened_out"
  if ("n_perm_observed" %in% colnames(traj)) {
    use_observed <- is.finite(traj$n_perm_observed) & traj$n_perm_observed >= 0L
    traj$n_perm_used[use_observed] <- as.integer(traj$n_perm_observed[use_observed])
  }
  diag_all <- .adaptive_perm_diagnostics(traj[[p_col]], traj$n_perm_used)
  traj$min_achievable_p <- diag_all$min_achievable_p
  traj$monte_carlo_se <- diag_all$monte_carlo_se
  traj <- .recompute_all_bh_q(traj)
  list(trajectories = traj, bins = bins)
}

.compute_adaptive_perm_config_hash <- function(n_perms_screen, n_perms_confirm,
                                               q_confirm_threshold, confirm_top_k,
                                               seed, kind = "temporal") {
  .provenance_hash(list(
    tag = paste0("adaptive_", kind, "_v2"),
    n_perms_screen = n_perms_screen,
    n_perms_confirm = n_perms_confirm,
    q_confirm_threshold = q_confirm_threshold,
    confirm_top_k = confirm_top_k,
    seed = seed
  ))
}

#' Adaptive two-stage temporal permutation null.
#'
#' Formal two-stage inference: screen with `n_perms_screen`, confirm candidates
#' with `n_perms_confirm`. Final `q_temporal` uses BH on the full trajectory
#' family after merging confirm-stage p-values. Non-candidates retain screen
#' statistics and `adaptive_perm_status = "screened_out"`.
#'
#' Stage 1 screens all trajectories with `n_perms_screen` permutations. Stage 2
#' re-runs permutation on the same full LR universe and substitutes confirm-stage
#' results only for candidate trajectories
#' (finite screen p-value at or below `q_confirm_threshold`, union the
#' `confirm_top_k` smallest screen p-values). Non-candidates retain screen-stage
#' p/q values and are tagged `adaptive_perm_status = "screened_out"`. Monte
#' Carlo diagnostics use each trajectory's observed permutation count, and all
#' confirm-stage null summaries and provenance columns are replaced together.
#' Confirmation currently recomputes the full LR universe; this API allocates
#' Monte Carlo precision adaptively but is not a computational acceleration.
#'
#' @inheritParams ComputeTemporalPermutation
#' @param n_perms_screen permutations for the screening stage (default 199).
#' @param n_perms_confirm permutations for the confirmation stage (default 4999).
#' @param q_confirm_threshold screen-stage p-value cutoff for confirmation
#'   (default 0.1; applied to raw screen p, not BH q).
#' @param confirm_top_k additionally confirm the smallest `confirm_top_k`
#'   finite screen p-values (default 20; set 0 to disable).
#' @return list as `ComputeTemporalPermutation()` plus `adaptive_perm_audit`
#'   with stage counts and `adaptive_perm_config_hash`.
#' @export
ComputeAdaptiveTemporalPermutation <- function(expr, celltypes, lr_db, pseudotime,
                                               lineage_weights,
                                               n_bins = 12L, min_pct = 0.05,
                                               min_n_eff = 5.0,
                                               min_valid_bins = 3L,
                                               min_valid_fraction = 0.25,
                                               n_perms_screen = 199L,
                                               n_perms_confirm = 4999L,
                                               q_confirm_threshold = 0.1,
                                               confirm_top_k = 20L,
                                               seed = 12345L,
                                               stat = c("amplitude", "max", "mean_score", "auc", "topk_mean"),
                                               top_k = 3L,
                                               min_perm_observed = 100L,
                                               min_perm_coverage = 0.8,
                                               block_labels = NULL,
                                               permutation_stratum = c(
                                                 "block_celltype_lineage", "block_celltype"
                                               )) {
  n_perms_screen <- .validate_positive_int(n_perms_screen, "n_perms_screen", min = 1L)
  n_perms_confirm <- .validate_positive_int(n_perms_confirm, "n_perms_confirm", min = 1L)
  q_confirm_threshold <- .validate_unit_interval(q_confirm_threshold, "q_confirm_threshold")
  confirm_top_k <- .validate_non_negative_int(confirm_top_k, "confirm_top_k")
  seed <- .validate_seed(seed)
  perm_args <- list(
    expr = expr,
    celltypes = celltypes,
    lr_db = lr_db,
    pseudotime = pseudotime,
    lineage_weights = lineage_weights,
    n_bins = n_bins,
    min_pct = min_pct,
    min_n_eff = min_n_eff,
    min_valid_bins = min_valid_bins,
    min_valid_fraction = min_valid_fraction,
    stat = stat,
    top_k = top_k,
    min_perm_observed = min_perm_observed,
    min_perm_coverage = min_perm_coverage,
    block_labels = block_labels,
    permutation_stratum = permutation_stratum
  )
  screen_out <- do.call(
    ComputeTemporalPermutation,
    c(perm_args, list(n_perms = n_perms_screen, seed = seed))
  )
  by <- c("pair_idx", "source", "target", "lineage")
  candidate_keys <- .adaptive_perm_select_candidates(
    screen_out$trajectories,
    p_col = "p_temporal",
    threshold = q_confirm_threshold,
    top_k = confirm_top_k,
    by = by
  )
  confirm_out <- NULL
  if (length(candidate_keys) > 0L) {
    confirm_out <- do.call(
      ComputeTemporalPermutation,
      utils::modifyList(perm_args, list(
        n_perms = n_perms_confirm,
        seed = seed + 1
      ))
    )
  }
  merged <- .merge_adaptive_temporal_perm(
    screen_out,
    confirm_out,
    candidate_keys,
    n_perms_screen,
    n_perms_confirm,
    p_col = "p_temporal",
    by = by
  )
  adaptive_hash <- .compute_adaptive_perm_config_hash(
    n_perms_screen, n_perms_confirm, q_confirm_threshold, confirm_top_k, seed,
    kind = "temporal"
  )
  merged$trajectories$adaptive_perm_config_hash <- adaptive_hash
  merged$bins$adaptive_perm_config_hash <- adaptive_hash
  n_confirmed <- sum(merged$trajectories$adaptive_perm_status == "confirmed")
  list(
    trajectories = merged$trajectories,
    bins = merged$bins,
    permutation_block_audit = screen_out$permutation_block_audit,
    analysis_data_hash = screen_out$analysis_data_hash,
    scoring_config_hash = screen_out$scoring_config_hash,
    calibration_config_hash = screen_out$calibration_config_hash,
    temporal_null_config_hash = adaptive_hash,
    stage_temporal_null_config_hashes = unique(
      merged$trajectories$temporal_null_config_hash
    ),
    adaptive_perm_config_hash = adaptive_hash,
    hill_scale = screen_out$hill_scale,
    adaptive_perm_audit = list(
      n_trajectories = nrow(merged$trajectories),
      n_candidates = length(candidate_keys),
      n_confirmed = n_confirmed,
      n_screened_out = sum(merged$trajectories$adaptive_perm_status == "screened_out"),
      n_perms_screen = n_perms_screen,
      n_perms_confirm = n_perms_confirm,
      q_confirm_threshold = q_confirm_threshold,
      confirm_top_k = confirm_top_k
    )
  )
}

#' Raw per-lineage entropy curves H_l(b).
#'
#' Marginalises over cell types to give
#' `H_l(b) = sum_ct hwsums[ct,l,b] / sum_ct wsums[ct,l,b]`. Bins with zero
#' denominator carry `NaN` (detect with `is.finite()`).
#'
#' @param intermediate the `intermediate` element returned by
#'   `ComputeTemporalLRWithIntermediate()`.
#' @return numeric matrix `[n_lineages × n_bins]` (rows = lineages,
#'   columns = pseudotime bins); pipe through `SmoothMA3()` for the
#'    smoothed curve `H̃_l(b)`.
#' @export
ComputeLineageEntropyRaw <- function(intermediate) {
  stopifnot(is.list(intermediate),
            all(c("hwsums", "wsums", "n_types", "n_lineages", "n_bins")
                %in% names(intermediate)))
  n_l <- as.integer(intermediate$n_lineages)
  n_b <- as.integer(intermediate$n_bins)
  flat <- rs_compute_lineage_entropy_raw(
    as.numeric(intermediate$hwsums),
    as.numeric(intermediate$wsums),
    as.integer(intermediate$n_types),
    n_l, n_b
  )
  # Rust returns row-major `l * n_bins + b`; reshape byrow=TRUE gives
  # rows = lineages, columns = bins.
  matrix(flat, nrow = n_l, ncol = n_b, byrow = TRUE)
}

#' MA(3) smoother for H_l(b) curves.
#'
#' Interior windows are full 3-bin `(H[b-1]+H[b]+H[b+1])/3`; boundary
#' windows at `b=1` and `b=n_bins` are 2-bin partial
#' `(H[0]+H[1])/2` / `(H[n-2]+H[n-1])/2`.  Any window containing `NaN`
#' propagates `NaN`.  When `n_bins < 6` the raw input is returned
#' unchanged.
#'
#' @param h_raw numeric matrix `[n_lineages × n_bins]` (typically from
#'   `ComputeLineageEntropyRaw()`).
#' @return numeric matrix same shape as `h_raw`, smoothed.
#' @export
SmoothMA3 <- function(h_raw) {
  m <- as.matrix(h_raw)
  n_l <- nrow(m)
  n_b <- ncol(m)
  # R col-major matrix → row-major flat (lineage outer, bin inner) via t().
  flat <- rs_smooth_ma3(
    as.numeric(t(m)),
    as.integer(n_l),
    as.integer(n_b)
  )
  out <- matrix(flat, nrow = n_l, ncol = n_b, byrow = TRUE)
  dimnames(out) <- dimnames(m)
  out
}

#' Locate bifurcation bins b*(l) per lineage.
#'
#' Applies a continuous two-bin rule to the smoothed entropy curve: the first
#' bin whose value and the next bin both fall strictly below `h_fork` wins.
#' The returned `reason` records edge cases:
#' \itemize{
#'   \item `"crossing"` — normal case, `b_star` is a 1-based bin index.
#'   \item `"c2_not_informative"` — curve has fewer than 4 finite bins or
#'     range below 0.05; `b_star = NA`.
#'   \item `"c3_no_low_entropy"` — informative curve but never satisfies
#'     the continuous-2-bin rule; `b_star = NA`.
#'   \item `"c4_boundary_onset"` — crossing at b <= 2 (1-based);
#'     `b_star` reported, but pre-fork taxonomy terms must be gated off
#'     downstream.
#' }
#'
#' @param h_tilde numeric matrix `[n_lineages × n_bins]`, typically
#'   `SmoothMA3(ComputeLineageEntropyRaw(intermediate))`.
#' @param h_fork numeric threshold (default 0.5).
#' @return data.frame with columns `lineage` (1-based integer),
#'   `b_star` (1-based integer or NA), `reason` (character).
#' @export
FindBifurcationBins <- function(h_tilde, h_fork = 0.5) {
  m   <- as.matrix(h_tilde)
  n_l <- nrow(m)
  n_b <- ncol(m)
  raw <- rs_find_bifurcation_bins(
    as.numeric(t(m)),  # col-major matrix → row-major flat (lineage outer)
    as.integer(n_l),
    as.integer(n_b),
    as.numeric(h_fork)
  )
  # Rust returns `b_star` as numeric with NaN encoding None; suppress
  # the coerce-to-integer NaN warning since NaN → NA_integer_ is exactly
  # the intent.
  b_int <- suppressWarnings(as.integer(raw$b_star))
  data.frame(
    lineage = seq_len(n_l),
    b_star  = b_int,
    reason  = raw$reason,
    stringsAsFactors = FALSE
  )
}

#' Compute sender-receiver pseudotime direction.
#'
#' Computes target pseudotime minus source pseudotime. Positive values mean
#' the receiver is later along pseudotime than the sender (sender leads).
#'
#' @param source_pseudotime numeric vector of sender pseudotime values.
#' @param target_pseudotime numeric vector of receiver pseudotime values.
#' @return numeric vector equal to `target_pseudotime - source_pseudotime`.
#' @export
ComputeDeltaTau <- function(source_pseudotime, target_pseudotime) {
  src <- as.numeric(source_pseudotime)
  tgt <- as.numeric(target_pseudotime)
  if (length(src) != length(tgt)) {
    stop("`source_pseudotime` and `target_pseudotime` must have equal length",
         call. = FALSE)
  }
  if (!all(is.finite(src)) || !all(is.finite(tgt))) {
    stop("`source_pseudotime` and `target_pseudotime` must be finite",
         call. = FALSE)
  }
  tgt - src
}

#' Compute six soft role scores per LR row.
#'
#' Returns `instructive`, `transition_like`, `reinforcing`, `feedback`,
#' `lateral`, and `bystander` scores with reliability and bifurcation gating. The
#' `bystander` score is QC-only and must not be used for LR ranking.
#'
#' @param scores_df data.frame produced by `ComputeTemporalLRWithIntermediate`
#'   (its `scores` element).  Must contain columns `bin` (1-based),
#'   `lineage` (1-based), and `n_eff_cells` (numeric Kish effective sample
#'   size for the row-level gate).
#' @param spec_self numeric vector: lineage specificity on the row's lineage.
#' @param spec_max_other numeric vector length `nrow(scores_df)`:
#'   `max_{l' != l} spec(l')`.  For single-lineage trajectories pass
#'   `rep(0, nrow(scores_df))`.
#' @param delta_tau numeric vector: **target pseudotime minus source
#'   pseudotime** per row. Positive values mean the receiver is later along
#'   pseudotime than the sender (sender leads). Use [ComputeDeltaTau()].
#' @param bif_df data.frame produced by `FindBifurcationBins`, with
#'   columns `lineage` (1-based), `b_star` (1-based integer or `NA`),
#'   `reason` (character: one of `"crossing"`, `"c2_not_informative"`,
#'   `"c3_no_low_entropy"`, `"c4_boundary_onset"`).
#' @param tau_0 numeric sigmoid scale (default 1.0; production callers often
#'   use `0.1 * sd(pseudotime)` per lineage).
#' @param h_fork numeric bifurcation threshold associated with `bif_df`
#'   (default 0.5); bifurcation bins are not recomputed here.
#' @param min_n_eff_cells numeric row-level reliability gate (default 5.0).
#' @param transition_window non-negative integer; bins with
#'   `abs(bin - b_star) <= transition_window` receive non-zero `transition_like`
#'   in the Rust core (default `0` = exact bifurcation bin only).
#' @return data.frame with 6 columns matching `scores_df` rows:
#'   `instructive`, `transition_like`, `reinforcing`, `feedback`, `lateral`,
#'   `bystander`. Attribute `role_config_hash` records the role parameters.
#' @export
ComputeRoleScores <- function(scores_df,
                              spec_self,
                              spec_max_other,
                              delta_tau,
                              bif_df,
                              tau_0           = 1.0,
                              h_fork          = 0.5,
                              min_n_eff_cells = 5.0,
                              transition_window = 0L) {
  stopifnot(is.data.frame(scores_df),
            all(c("bin", "lineage", "n_eff_cells") %in% colnames(scores_df)),
            is.data.frame(bif_df),
            all(c("lineage", "b_star", "reason") %in% colnames(bif_df)))

  n_rows     <- nrow(scores_df)
  n_lineages <- nrow(bif_df)
  stopifnot(length(spec_self)      == n_rows,
            length(spec_max_other) == n_rows,
            length(delta_tau)      == n_rows)
  tau_0 <- .validate_positive_num(tau_0, "tau_0")
  h_fork <- .validate_unit_interval(h_fork, "h_fork")
  min_n_eff_cells <- as.numeric(min_n_eff_cells)
  if (length(min_n_eff_cells) != 1L || !is.finite(min_n_eff_cells) ||
      min_n_eff_cells < 0) {
    stop("`min_n_eff_cells` must be a finite numeric value >= 0", call. = FALSE)
  }
  transition_window <- .validate_positive_int(transition_window, "transition_window",
                                                min = 0L)
  scores_df$bin <- .as_integerish(scores_df$bin, "scores_df$bin")
  scores_df$lineage <- .assert_indices_in_range(
    scores_df$lineage, 1L, n_lineages, "scores_df$lineage"
  )
  if (any(scores_df$bin < 1L)) {
    stop("`scores_df$bin` must contain positive 1-based integers", call. = FALSE)
  }
  if (!all(is.finite(spec_self)) || !all(spec_self >= 0 & spec_self <= 1)) {
    stop("`spec_self` must be finite values in [0, 1]", call. = FALSE)
  }
  if (!all(is.finite(spec_max_other)) ||
      !all(spec_max_other >= 0 & spec_max_other <= 1)) {
    stop("`spec_max_other` must be finite values in [0, 1]", call. = FALSE)
  }
  if (!all(is.finite(delta_tau))) {
    stop("`delta_tau` must be finite; use `ComputeDeltaTau()`", call. = FALSE)
  }
  if (!all(is.finite(scores_df$n_eff_cells))) {
    stop("`scores_df$n_eff_cells` must be finite", call. = FALSE)
  }

  # Ensure bif_df is ordered by lineage
  # aligns with the Rust side's 1-based `row_lineage` lookup.
  bif_df$lineage <- .as_integerish(bif_df$lineage, "bif_df$lineage")
  bif_df <- bif_df[order(bif_df$lineage), , drop = FALSE]
  if (!identical(bif_df$lineage, seq_len(n_lineages))) {
    stop("`bif_df$lineage` must contain each 1-based lineage exactly once",
         call. = FALSE)
  }

  # Encode reason as integer code (same mapping as Rust side).
  reason_levels <- c("crossing", "c2_not_informative",
                     "c3_no_low_entropy", "c4_boundary_onset")
  reason_code   <- match(as.character(bif_df$reason), reason_levels) - 1L
  if (any(is.na(reason_code))) {
    stop("ComputeRoleScores: bif_df$reason contains unknown levels: ",
         paste(unique(bif_df$reason[is.na(reason_code)]), collapse = ", "))
  }

  # NA b_star (C2/C3) → NaN on the Rust side, which decodes to None.
  b_star_num <- as.numeric(.as_integerish(
    bif_df$b_star, "bif_df$b_star", allow_na = TRUE
  ))
  if (any(b_star_num < 1, na.rm = TRUE)) {
    stop("finite `bif_df$b_star` values must be positive 1-based integers",
         call. = FALSE)
  }
  b_star_num[is.na(b_star_num)] <- NaN

  raw <- rs_compute_role_scores(
    as.integer(scores_df$bin),
    as.integer(scores_df$lineage),
    as.numeric(spec_self),
    as.numeric(spec_max_other),
    as.numeric(delta_tau),
    as.numeric(scores_df$n_eff_cells),  # row-level reliability gate
    b_star_num,
    as.integer(reason_code),
    as.numeric(tau_0),
    as.numeric(h_fork),
    as.numeric(min_n_eff_cells),
    as.integer(transition_window)
  )
  out <- data.frame(
    instructive = raw$instructive,
    transition_like = raw$transition_like,
    reinforcing = raw$reinforcing,
    feedback    = raw$feedback,
    lateral     = raw$lateral,
    bystander   = raw$bystander,
    stringsAsFactors = FALSE
  )
  attr(out, "role_config_hash") <- .compute_role_config_hash(
    tau_0, h_fork, min_n_eff_cells, transition_window
  )
  out
}

# ─── Role feature engineering ───────────────────────────────────────────

.compute_lineage_bin_bounds <- function(pseudotime, lineage_weights, n_bins) {
  n_bins <- .validate_positive_int(n_bins, "n_bins", min = 2L)
  pt <- .assert_pseudotime_1d(pseudotime, length(pseudotime))
  lw <- .assert_lineage_weights(lineage_weights, length(pt))
  raw <- rs_compute_bin_bounds(
    pt,
    as.numeric(t(lw)),
    as.integer(ncol(lw)),
    n_bins
  )
  lapply(raw, as.numeric)
}

#' Shared temporal grid boundaries for fork-local contrasts.
#'
#' Returns the single `n_bins + 1` boundary vector used when
#' `common_grid = TRUE`.  The grid is equal-frequency over the pooled
#' on-trajectory mass (a cell's mass is its maximum lineage weight), so every
#' lineage is scored on the same pseudotime intervals and a bin index denotes
#' the same relative graph-time interval on every branch.
#'
#' @param pseudotime numeric vector; `NA`/`NaN` cells are excluded.
#' @param lineage_weights numeric matrix \[n_cells × n_lineages\].
#' @param n_bins integer, number of shared temporal bins.
#' @return numeric vector of length `n_bins + 1` from 0 to 1.  To reuse it for
#'   per-sample rescoring, replicate it into one row per lineage and pass as
#'   `bin_bounds` to `ComputeTemporalLRFixedBounds()`.
#' @export
ComputeCommonBinBounds <- function(pseudotime, lineage_weights, n_bins = 12L) {
  n_bins <- .validate_positive_int(n_bins, "n_bins", min = 2L)
  pt <- .assert_pseudotime_1d(pseudotime, length(pseudotime))
  lw <- .assert_lineage_weights(lineage_weights, length(pt))
  as.numeric(rs_compute_common_bin_bounds(
    pt,
    as.numeric(t(lw)),
    as.integer(ncol(lw)),
    n_bins
  ))
}

.lineage_cell_bin_index <- function(pt, bounds, n_bins) {
  if (!is.finite(pt)) {
    return(NA_integer_)
  }
  bin <- sum(bounds <= pt) - 1L
  bin <- max(0L, min(bin, n_bins - 1L)) + 1L
  as.integer(bin)
}

.compute_celltype_lineage_bin_timing <- function(celltypes, pseudotime,
                                                 lineage_weights, n_bins,
                                                 ct_levels_full = NULL) {
  ct <- .assert_celltypes(celltypes, length(pseudotime))
  ct_levels_local <- attr(ct, "levels")
  if (!is.null(ct_levels_full)) {
    ct_levels_full <- as.character(ct_levels_full)
    ct_ix <- match(as.character(ct_levels_local[ct]), ct_levels_full)
    if (any(is.na(ct_ix))) {
      stop("Block `celltypes` levels must be a subset of `ct_levels_full`",
           call. = FALSE)
    }
    n_types <- length(ct_levels_full)
    ct_levels <- ct_levels_full
  } else {
    ct_ix <- ct
    n_types <- max(ct)
    ct_levels <- as.character(ct_levels_local)
  }
  n_lineages <- ncol(lineage_weights)
  bounds <- .compute_lineage_bin_bounds(pseudotime, lineage_weights, n_bins)
  wt_sum <- array(0, dim = c(n_types, n_lineages, n_bins))
  pt_sum <- array(0, dim = c(n_types, n_lineages, n_bins))
  n_cells <- length(pseudotime)
  for (cell in seq_len(n_cells)) {
    pt <- pseudotime[cell]
    if (!is.finite(pt)) {
      next
    }
    ct_i <- ct_ix[cell]
    for (lin in seq_len(n_lineages)) {
      w <- lineage_weights[cell, lin]
      if (w < 1e-8) {
        next
      }
      bin <- .lineage_cell_bin_index(pt, bounds[[lin]], n_bins)
      if (is.na(bin)) {
        next
      }
      wt_sum[ct_i, lin, bin] <- wt_sum[ct_i, lin, bin] + w
      pt_sum[ct_i, lin, bin] <- pt_sum[ct_i, lin, bin] + w * pt
    }
  }
  timing <- pt_sum / pmax(wt_sum, 1e-12)
  list(timing = timing, ct_levels = ct_levels, wt_sum = wt_sum)
}

.resolve_scores_celltype_index <- function(scores_col, ct_levels) {
  if (is.numeric(scores_col) &&
      all(scores_col %in% seq_along(ct_levels))) {
    return(as.integer(scores_col))
  }
  idx <- match(as.character(scores_col), as.character(ct_levels))
  if (any(is.na(idx))) {
    stop("Could not map `source`/`target` labels to `celltypes` levels",
         call. = FALSE)
  }
  idx
}

.compute_lineage_specificity <- function(scores_df, n_total_lineages = NULL) {
  if (!"score" %in% colnames(scores_df)) {
    stop("`temporal_scores` must contain a `score` column for lineage specificity",
         call. = FALSE)
  }
  n <- nrow(scores_df)
  spec_self <- rep(NA_real_, n)
  spec_max_other <- rep(NA_real_, n)
  n_total_out <- rep(NA_integer_, n)
  n_evaluable_out <- rep(NA_integer_, n)
  n_supported_out <- rep(NA_integer_, n)
  if (is.null(n_total_lineages)) {
    n_total_lineages <- length(unique(
      .as_integerish(scores_df$lineage, "temporal_scores$lineage")
    ))
  } else {
    n_total_lineages <- .validate_positive_int(
      n_total_lineages, "n_total_lineages", min = 1L
    )
  }
  key <- paste(
    scores_df$pair_idx, scores_df$source, scores_df$target, scores_df$bin,
    sep = "\001"
  )
  for (k in unique(key)) {
    idx <- which(key == k)
    sc_raw <- as.numeric(scores_df$score[idx])
    sc <- pmax(sc_raw, 0)
    evaluable <- is.finite(sc_raw)
    n_eval <- sum(evaluable)
    n_supp <- sum(sc > 0, na.rm = TRUE)
    n_total_out[idx] <- n_total_lineages
    n_evaluable_out[idx] <- n_eval
    n_supported_out[idx] <- n_supp
    if (n_eval < 2L) {
      next
    }
    denom <- sum(sc, na.rm = TRUE)
    if (denom <= 0 || !is.finite(denom)) {
      next
    }
    sp <- sc / denom
    for (j in seq_along(idx)) {
      spec_self[idx[j]] <- sp[j]
      other <- sp[-j]
      spec_max_other[idx[j]] <- if (length(other)) max(other, na.rm = TRUE) else 0
    }
  }
  list(
    spec_self = spec_self,
    spec_max_other = spec_max_other,
    n_total_lineages = n_total_out,
    n_evaluable_lineages = n_evaluable_out,
    n_supported_lineages = n_supported_out
  )
}

.default_tau_0_from_pseudotime <- function(pseudotime, lineage_weights) {
  n_lineages <- ncol(lineage_weights)
  vapply(seq_len(n_lineages), function(lin) {
    w <- lineage_weights[, lin]
    ok <- is.finite(pseudotime) & w > 0
    if (!any(ok)) {
      return(0.1)
    }
    pt <- pseudotime[ok]
    ww <- w[ok]
    mu <- stats::weighted.mean(pt, ww)
    var <- stats::weighted.mean((pt - mu)^2, ww)
    sd_pt <- sqrt(max(var, 0))
    if (sd_pt <= 0) {
      0.1
    } else {
      0.1 * sd_pt
    }
  }, numeric(1))
}

.resolve_bifurcation_df <- function(bifurcation, intermediate, h_tilde, h_fork) {
  if (!is.null(bifurcation)) {
    stopifnot(is.data.frame(bifurcation),
              all(c("lineage", "b_star", "reason") %in% colnames(bifurcation)))
    return(bifurcation)
  }
  if (!is.null(h_tilde)) {
    return(FindBifurcationBins(h_tilde, h_fork = h_fork))
  }
  if (!is.null(intermediate)) {
    h_raw <- ComputeLineageEntropyRaw(intermediate)
    h_smooth <- SmoothMA3(h_raw)
    return(FindBifurcationBins(h_smooth, h_fork = h_fork))
  }
  stop(
    "Provide `bifurcation`, `h_tilde`, or `intermediate` to locate bifurcation bins",
    call. = FALSE
  )
}

.attach_fully_supported <- function(df, evidence, by) {
  if (is.null(evidence)) {
    df$fully_supported <- NA
    return(df)
  }
  stopifnot(is.data.frame(evidence))
  if ("fully_supported" %in% colnames(evidence)) {
    keep <- unique(c(by, "fully_supported"))
    keep <- intersect(keep, colnames(evidence))
    cal <- evidence[, keep, drop = FALSE]
    cal <- cal[!duplicated(.trajectory_key_col(cal, by)), , drop = FALSE]
    key <- .trajectory_key_col(df, by)
    cal_key <- .trajectory_key_col(cal, by)
    df$fully_supported <- cal$fully_supported[match(key, cal_key)]
    return(df)
  }
  df$fully_supported <- NA
  df
}

#' Auto-compute role-inference features from temporal LR output.
#'
#' Derives lineage specificity (`spec_self`, `spec_max_other`), sender/receiver
#' bin timing, signed lag (`delta_tau`), and transition-relative position
#' (`bin - b*(l)`) per Methods §5.4. Specificity uses
#' `score(p,s,t,l,b) / sum_{l'} score(p,s,t,l',b)`.
#'
#' @param temporal_scores data.frame from `ComputeTemporalLR()` or
#'   `ComputeTemporalLRWithIntermediate()$scores`.
#' @param lineage_weights numeric matrix `[n_cells × n_lineages]`.
#' @param pseudotime numeric vector of length `n_cells`.
#' @param celltypes integer vector (1-based) of length `n_cells`.
#' @param bifurcation optional `FindBifurcationBins()` output; computed from
#'   `h_tilde` or `intermediate` when omitted.
#' @param intermediate optional `ComputeTemporalLRWithIntermediate()$intermediate`.
#' @param h_tilde optional smoothed lineage entropy matrix for bifurcation.
#' @param h_fork bifurcation threshold when inferring `bifurcation`.
#' @param n_bins number of pseudotime bins (inferred from `temporal_scores` when NULL).
#' @param evidence optional `MergeEvidence()` output; attaches trajectory-level
#'   `fully_supported` when present.
#' @param by trajectory key columns for evidence merge.
#' @param compute_lag when `TRUE`, attach ligand-receptor timing lag columns via
#'   [ComputeLigandReceptorLag()] and [MergeLagIntoRoleFeatures()].
#' @param lag_method lag estimator passed to `ComputeLigandReceptorLag()`.
#' @param lag_threshold minimum `|lag|` (pseudotime units) for `lag_exceeds_threshold`;
#'   defaults to one bin width when `NULL`. Deprecated alias: `lag_significant`.
#' @param timing_source timing signal for role scoring:
#'   `"cell_state_tau"` (cell-type pseudotime centroids),
#'   `"ligand_receptor_lag"` (`receptor_time - ligand_time` from lag side profiles),
#'   `"ligand_target_lag"` is retained as a legacy alias of
#'   `"ligand_receptor_lag"`; no downstream target programme is inferred.
#' @param expr optional expression matrix for side-profile rescoring when
#'   `ligand_score` / `receptor_score` are absent from `temporal_scores`.
#' @param lr_db optional LR database (required with `expr` for gene-level profiles).
#' @return data.frame parallel to `temporal_scores` with feature columns and
#'   attributes `bif_df`, `tau_0_by_lineage`, `n_bins`. When `compute_lag = TRUE`,
#'   trajectory-level lag columns are replicated on each bin-row.
#' @export
ComputeRoleFeatures <- function(temporal_scores,
                                lineage_weights,
                                pseudotime,
                                celltypes,
                                bifurcation = NULL,
                                intermediate = NULL,
                                h_tilde = NULL,
                                h_fork = 0.5,
                                n_bins = NULL,
                                evidence = NULL,
                                by = c("pair_idx", "source", "target", "lineage"),
                                compute_lag = FALSE,
                                lag_method = c("peak", "centroid", "cross_cor"),
                                lag_threshold = NULL,
                                timing_source = c("cell_state_tau",
                                                  "ligand_receptor_lag",
                                                  "ligand_target_lag"),
                                expr = NULL,
                                lr_db = NULL) {
  stopifnot(is.data.frame(temporal_scores))
  req <- c("pair_idx", "source", "target", "lineage", "bin")
  missing <- setdiff(req, colnames(temporal_scores))
  if (length(missing) > 0L) {
    stop("`temporal_scores` must contain columns: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  by <- as.character(by)
  timing_source <- match.arg(timing_source)
  n_cells <- length(pseudotime)
  lw <- .assert_lineage_weights(lineage_weights, n_cells)
  pt <- .assert_pseudotime_1d(pseudotime, n_cells)
  h_fork <- .validate_unit_interval(h_fork, "h_fork")
  n_bins <- .infer_n_bins(temporal_scores, n_bins)
  bif_df <- .resolve_bifurcation_df(bifurcation, intermediate, h_tilde, h_fork)

  spec <- .compute_lineage_specificity(temporal_scores, n_total_lineages = ncol(lw))
  timing <- .compute_celltype_lineage_bin_timing(celltypes, pt, lw, n_bins)
  ct_levels <- timing$ct_levels
  if (!is.null(attr(temporal_scores, "celltype_levels"))) {
    ct_levels <- attr(temporal_scores, "celltype_levels")
  }
  src_ix <- .resolve_scores_celltype_index(temporal_scores$source, ct_levels)
  tgt_ix <- .resolve_scores_celltype_index(temporal_scores$target, ct_levels)
  lin <- .as_integerish(temporal_scores$lineage, "temporal_scores$lineage")
  bin <- .as_integerish(temporal_scores$bin, "temporal_scores$bin")
  n <- nrow(temporal_scores)
  source_timing <- numeric(n)
  target_timing <- numeric(n)
  for (i in seq_len(n)) {
    l <- lin[i]
    b <- bin[i]
    source_timing[i] <- timing$timing[src_ix[i], l, b]
    target_timing[i] <- timing$timing[tgt_ix[i], l, b]
  }
  delta_tau <- ComputeDeltaTau(source_timing, target_timing)
  b_star <- bif_df$b_star[match(lin, bif_df$lineage)]
  transition_relative_position <- as.numeric(bin) - as.numeric(b_star)
  bif_position <- ifelse(
    is.na(b_star),
    "time_only",
    ifelse(bin < b_star, "pre_bifurcation",
           ifelse(bin == b_star, "transition",
                  "post_bifurcation"))
  )

  out <- data.frame(
    temporal_scores,
    spec_self = spec$spec_self,
    spec_max_other = spec$spec_max_other,
    n_total_lineages = spec$n_total_lineages,
    n_evaluable_lineages = spec$n_evaluable_lineages,
    n_supported_lineages = spec$n_supported_lineages,
    source_timing = source_timing,
    target_timing = target_timing,
    signed_lag = delta_tau,
    delta_tau = delta_tau,
    b_star = b_star,
    transition_relative_position = transition_relative_position,
    bif_position = bif_position,
    stringsAsFactors = FALSE
  )
  out <- .attach_fully_supported(out, evidence, by)
  attr(out, "bif_df") <- bif_df
  attr(out, "tau_0_by_lineage") <- .default_tau_0_from_pseudotime(pt, lw)
  attr(out, "n_bins") <- n_bins
  attr(out, "timing_source") <- timing_source
  if (isTRUE(compute_lag)) {
    lag_scores <- temporal_scores
    if (!is.null(expr) && !is.null(lr_db)) {
      lag_scores <- .enrich_temporal_scores_lr_sides(
        temporal_scores, expr, lr_db, celltypes, pt, lw, n_bins, ct_levels
      )
    }
    lag_df <- ComputeLigandReceptorLag(
      temporal_scores = lag_scores,
      lineage_weights = lw,
      pseudotime = pt,
      celltypes = celltypes,
      n_bins = n_bins,
      method = lag_method,
      lag_threshold = lag_threshold,
      expr = expr,
      lr_db = lr_db,
      by = by
    )
    out <- MergeLagIntoRoleFeatures(out, lag_df, by = by)
    attr(out, "lag_df") <- lag_df
  }
  out <- .apply_role_timing_source(out, timing_source)
  out
}

# ─── Ligand-receptor timing lag ─────────────────────────────────────────

.lag_side_profile_columns <- function(temporal_scores) {
  if (all(c("ligand_score", "receptor_score") %in% colnames(temporal_scores))) {
    return(c("ligand_score", "receptor_score"))
  }
  NULL
}

.weighted_complex_gene_mean <- function(expr, gene_rows, cells, weights) {
  if (!length(cells)) {
    return(NA_real_)
  }
  w <- as.numeric(weights)
  if (length(w) != length(cells)) {
    stop("`weights` length must equal `cells` length", call. = FALSE)
  }
  ok <- is.finite(w) & w > 0
  if (!any(ok)) {
    return(NA_real_)
  }
  cells <- cells[ok]
  w <- w[ok]
  if (length(gene_rows) == 1L) {
    x <- as.numeric(expr[gene_rows, cells, drop = FALSE])
    return(sum(w * x) / sum(w))
  }
  subunit_means <- vapply(gene_rows, function(gr) {
    x <- as.numeric(expr[gr, cells, drop = FALSE])
    sum(w * x) / sum(w)
  }, numeric(1))
  min(subunit_means, na.rm = TRUE)
}

.compute_lr_side_profiles <- function(traj_df, temporal_scores, expr, lr_db,
                                      celltypes, pseudotime, lineage_weights,
                                      n_bins, ct_levels, src_ix, tgt_ix, lin_id,
                                      wt_sum) {
  n_bins <- as.integer(n_bins)
  score_vec <- rep(NA_real_, n_bins)
  lig_vec <- rep(NA_real_, n_bins)
  rec_vec <- rep(NA_real_, n_bins)
  for (i in seq_len(nrow(traj_df))) {
    b <- as.integer(traj_df$bin[i])
    if (b >= 1L && b <= n_bins) {
      score_vec[b] <- as.numeric(traj_df$score[i])
      if ("ligand_score" %in% colnames(traj_df)) {
        lig_vec[b] <- as.numeric(traj_df$ligand_score[i])
      }
      if ("receptor_score" %in% colnames(traj_df)) {
        rec_vec[b] <- as.numeric(traj_df$receptor_score[i])
      }
    }
  }
  if (all(is.finite(lig_vec)) && all(is.finite(rec_vec))) {
    return(list(
      ligand = lig_vec,
      receptor = rec_vec,
      lag_profile_source = "explicit_columns"
    ))
  }
  if (!is.null(expr) && !is.null(lr_db)) {
    lr_db <- .validate_lr_db(lr_db)
    pair_idx <- as.integer(traj_df$pair_idx[1L])
    lig_rows <- lr_db$ligand_rows[[pair_idx]]
    rec_rows <- lr_db$receptor_rows[[pair_idx]]
    bounds <- .compute_lineage_bin_bounds(pseudotime, lineage_weights, n_bins)
    bnd <- bounds[[lin_id]]
    ct <- .assert_celltypes(celltypes, length(pseudotime))
    n_cells <- length(pseudotime)
    src_label <- as.character(ct_levels[src_ix])
    tgt_label <- as.character(ct_levels[tgt_ix])
    ct_labels <- as.character(attr(ct, "levels")[ct])
    for (b in seq_len(n_bins)) {
      src_cells <- integer(0)
      tgt_cells <- integer(0)
      src_w <- numeric(0)
      tgt_w <- numeric(0)
      for (cell in seq_len(n_cells)) {
        pt <- pseudotime[cell]
        if (!is.finite(pt)) {
          next
        }
        w <- lineage_weights[cell, lin_id]
        if (w < 1e-8) {
          next
        }
        bin_i <- .lineage_cell_bin_index(pt, bnd, n_bins)
        if (is.na(bin_i) || bin_i != b) {
          next
        }
        if (ct_labels[cell] == src_label) {
          src_cells <- c(src_cells, cell)
          src_w <- c(src_w, w)
        }
        if (ct_labels[cell] == tgt_label) {
          tgt_cells <- c(tgt_cells, cell)
          tgt_w <- c(tgt_w, w)
        }
      }
      lig_vec[b] <- .weighted_complex_gene_mean(expr, lig_rows, src_cells, src_w)
      rec_vec[b] <- .weighted_complex_gene_mean(expr, rec_rows, tgt_cells, tgt_w)
    }
    return(list(
      ligand = lig_vec,
      receptor = rec_vec,
      lag_profile_source = "expr_lr_db"
    ))
  }
  src_mass <- wt_sum[src_ix, lin_id, ]
  tgt_mass <- wt_sum[tgt_ix, lin_id, ]
  list(
    ligand = score_vec * src_mass,
    receptor = score_vec * tgt_mass,
    lag_profile_source = "mass_proxy"
  )
}

.enrich_temporal_scores_lr_sides <- function(temporal_scores,
                                              expr,
                                              lr_db,
                                              celltypes,
                                              pseudotime,
                                              lineage_weights,
                                              n_bins,
                                              ct_levels = NULL) {
  EnrichTemporalScoresWithLRSides(
    temporal_scores, expr, lr_db, celltypes, pseudotime, lineage_weights,
    n_bins, ct_levels
  )
}

#' Add gene-level ligand/receptor side scores to temporal LR rows.
#'
#' Populates `ligand_score` and `receptor_score` from `expr` and `lr_db` when
#' those columns are not already complete. Used internally by
#' [ComputeRoleFeatures()] and [ComputeLigandReceptorLag()].
#'
#' @param temporal_scores data.frame with trajectory keys and `bin`.
#' @param expr `dgCMatrix` (genes × cells).
#' @param lr_db ligand/receptor row lists.
#' @param celltypes,pseudotime,lineage_weights see [ComputeTemporalLR()].
#' @param n_bins number of pseudotime bins.
#' @param ct_levels optional cell-type level labels.
#' @return `temporal_scores` with `ligand_score` / `receptor_score` columns.
#' @export
EnrichTemporalScoresWithLRSides <- function(temporal_scores,
                                              expr,
                                              lr_db,
                                              celltypes,
                                              pseudotime,
                                              lineage_weights,
                                              n_bins,
                                              ct_levels = NULL) {
  if (!is.null(.lag_side_profile_columns(temporal_scores))) {
    return(temporal_scores)
  }
  if (is.null(expr) || is.null(lr_db)) {
    return(temporal_scores)
  }
  n_bins <- as.integer(n_bins)
  timing <- .compute_celltype_lineage_bin_timing(
    celltypes, pseudotime, lineage_weights, n_bins,
    ct_levels_full = ct_levels
  )
  ct_levels <- timing$ct_levels
  out <- temporal_scores
  out$ligand_score <- NA_real_
  out$receptor_score <- NA_real_
  by <- c("pair_idx", "source", "target", "lineage")
  key <- .trajectory_key_col(out, by)
  for (k in unique(key)) {
    idx <- which(key == k)
    traj_df <- out[idx, , drop = FALSE]
    src_ix <- .resolve_scores_celltype_index(traj_df$source[1L], ct_levels)
    tgt_ix <- .resolve_scores_celltype_index(traj_df$target[1L], ct_levels)
    lin_id <- .as_integerish(traj_df$lineage[1L], "lineage")
    prof <- .compute_lr_side_profiles(
      traj_df, temporal_scores, expr, lr_db,
      celltypes, pseudotime, lineage_weights, n_bins, ct_levels,
      src_ix, tgt_ix, lin_id, timing$wt_sum
    )
    for (i in idx) {
      b <- as.integer(out$bin[i])
      if (b >= 1L && b <= n_bins) {
        out$ligand_score[i] <- prof$ligand[b]
        out$receptor_score[i] <- prof$receptor[b]
      }
    }
  }
  out
}

.infer_lag_bin_width <- function(timing_src, timing_tgt, n_bins) {
  tv <- c(timing_src, timing_tgt)
  tv <- tv[is.finite(tv)]
  if (length(tv) >= 2L) {
    rg <- diff(range(tv))
    if (is.finite(rg) && rg > 0) {
      return(rg / max(n_bins - 1L, 1L))
    }
  }
  1 / max(n_bins, 1L)
}

.lag_se_heuristic <- function(lag_value, bin_width, n_valid_bins) {
  base <- pmax(abs(lag_value), bin_width, na.rm = TRUE)
  denom <- sqrt(pmax(as.numeric(n_valid_bins), 1))
  pmax(base / denom, bin_width * 0.05, na.rm = TRUE)
}

.lag_from_peak <- function(lig_prof, rec_prof, timing_src, timing_tgt) {
  lig_ok <- is.finite(lig_prof)
  rec_ok <- is.finite(rec_prof)
  if (!any(lig_ok) || !any(rec_ok)) {
    return(list(
      ligand_time = NA_real_, receptor_time = NA_real_, lag = NA_real_,
      ligand_peak_time = NA_real_, receptor_peak_time = NA_real_,
      lag_peak = NA_real_, n_valid_bins = 0L
    ))
  }
  lig_bin <- which.max(replace(lig_prof, !lig_ok, -Inf))
  rec_bin <- which.max(replace(rec_prof, !rec_ok, -Inf))
  lig_t <- timing_src[lig_bin]
  rec_t <- timing_tgt[rec_bin]
  lag <- rec_t - lig_t
  list(
    ligand_time = lig_t,
    receptor_time = rec_t,
    lag = lag,
    ligand_peak_time = lig_t,
    receptor_peak_time = rec_t,
    lag_peak = lag,
    n_valid_bins = sum(lig_ok | rec_ok)
  )
}

.lag_from_centroid <- function(lig_prof, rec_prof, timing_src, timing_tgt) {
  lig_ok <- is.finite(lig_prof)
  rec_ok <- is.finite(rec_prof)
  lig <- replace(as.numeric(lig_prof), !lig_ok, NA_real_)
  rec <- replace(as.numeric(rec_prof), !rec_ok, NA_real_)
  sw_l <- sum(lig, na.rm = TRUE)
  sw_r <- sum(rec, na.rm = TRUE)
  if (!is.finite(sw_l) || !is.finite(sw_r) || sw_l <= 0 || sw_r <= 0) {
    return(list(
      ligand_time = NA_real_, receptor_time = NA_real_, lag = NA_real_,
      ligand_centroid_time = NA_real_, receptor_centroid_time = NA_real_,
      lag_centroid = NA_real_, n_valid_bins = 0L
    ))
  }
  lig_t <- sum(lig * timing_src, na.rm = TRUE) / sw_l
  rec_t <- sum(rec * timing_tgt, na.rm = TRUE) / sw_r
  lag <- rec_t - lig_t
  list(
    ligand_time = lig_t,
    receptor_time = rec_t,
    lag = lag,
    ligand_centroid_time = lig_t,
    receptor_centroid_time = rec_t,
    lag_centroid = lag,
    n_valid_bins = sum(lig_ok | rec_ok)
  )
}

.lag_from_cross_cor <- function(lig_prof, rec_prof, timing_src, timing_tgt,
                                bin_width, max_lag_bins = NULL) {
  lig <- as.numeric(lig_prof)
  rec <- as.numeric(rec_prof)
  n <- length(lig)
  if (n < 2L || sum(is.finite(lig)) < 2L || sum(is.finite(rec)) < 2L) {
    return(list(
      ligand_time = NA_real_, receptor_time = NA_real_, lag = NA_real_,
      lag_cross_cor = NA_real_, cross_cor_max = NA_real_, n_valid_bins = 0L,
      lag_tie = NA, n_equivalent_lags = 0L
    ))
  }
  if (is.null(max_lag_bins)) {
    max_lag_bins <- max(1L, n %/% 2L)
  } else {
    max_lag_bins <- .validate_positive_int(max_lag_bins, "max_lag_bins", min = 1L)
    max_lag_bins <- min(max_lag_bins, n - 1L)
  }
  candidate_shift <- integer(0)
  candidate_cor <- numeric(0)
  candidate_n_pairs <- integer(0)
  for (shift in seq(-max_lag_bins, max_lag_bins)) {
    if (shift >= 0L) {
      lig_idx <- seq_len(n - shift)
      rec_idx <- shift + seq_len(n - shift)
    } else {
      s <- -shift
      lig_idx <- s + seq_len(n - s)
      rec_idx <- seq_len(n - s)
    }
    aligned_ok <- is.finite(lig[lig_idx]) & is.finite(rec[rec_idx])
    if (sum(aligned_ok) < 2L) {
      next
    }
    l1 <- lig[lig_idx][aligned_ok]
    r1 <- rec[rec_idx][aligned_ok]
    if (stats::sd(l1) <= 0 || stats::sd(r1) <= 0) {
      next
    }
    cc <- stats::cor(l1, r1)
    if (is.finite(cc)) {
      candidate_shift <- c(candidate_shift, shift)
      candidate_cor <- c(candidate_cor, cc)
      candidate_n_pairs <- c(candidate_n_pairs, length(l1))
    }
  }
  if (!length(candidate_cor)) {
    return(list(
      ligand_time = NA_real_, receptor_time = NA_real_, lag = NA_real_,
      lag_cross_cor = NA_real_, cross_cor_max = NA_real_, n_valid_bins = 0L,
      lag_tie = NA, n_equivalent_lags = 0L
    ))
  }
  tol <- 1e-10
  best_cor <- max(candidate_cor)
  equivalent <- which(abs(candidate_cor - best_cor) <= tol)
  n_equivalent_lags <- length(equivalent)
  best_n_pairs <- max(candidate_n_pairs[equivalent])
  preferred <- equivalent[candidate_n_pairs[equivalent] == best_n_pairs]
  min_abs_shift <- min(abs(candidate_shift[preferred]))
  preferred <- preferred[abs(candidate_shift[preferred]) == min_abs_shift]
  best_lag_bins <- if (length(preferred) == 1L) candidate_shift[preferred] else 0L
  lag <- as.numeric(best_lag_bins) * bin_width
  timing_ok <- c(
    timing_src[is.finite(lig)],
    timing_tgt[is.finite(rec)]
  )
  timing_ok <- timing_ok[is.finite(timing_ok)]
  mid <- if (length(timing_ok)) stats::median(timing_ok) else 0.5
  if (!is.finite(mid)) {
    mid <- 0.5
  }
  list(
    ligand_time = mid - lag / 2,
    receptor_time = mid + lag / 2,
    lag = lag,
    lag_cross_cor = lag,
    cross_cor_max = if (is.finite(best_cor)) best_cor else NA_real_,
    n_valid_bins = best_n_pairs,
    lag_tie = n_equivalent_lags > 1L,
    n_equivalent_lags = n_equivalent_lags
  )
}

.compute_trajectory_lag_row <- function(traj_df, temporal_scores, timing, wt_sum,
                                        ct_levels, expr, lr_db, celltypes,
                                        pseudotime, lineage_weights, n_bins,
                                        method, lag_threshold, max_lag_bins,
                                        by, compute_lag_null = FALSE,
                                        n_lag_perm = 199L, lag_null_seed = 41L) {
  src_ix <- .resolve_scores_celltype_index(traj_df$source[1L], ct_levels)
  tgt_ix <- .resolve_scores_celltype_index(traj_df$target[1L], ct_levels)
  lin_id <- .as_integerish(traj_df$lineage[1L], "lineage")
  timing_src <- timing[src_ix, lin_id, ]
  timing_tgt <- timing[tgt_ix, lin_id, ]
  profiles <- .compute_lr_side_profiles(
    traj_df, temporal_scores, expr, lr_db,
    celltypes, pseudotime, lineage_weights, n_bins, ct_levels,
    src_ix, tgt_ix, lin_id, wt_sum
  )
  bin_width <- .infer_lag_bin_width(timing_src, timing_tgt, n_bins)
  if (is.null(lag_threshold)) {
    lag_threshold <- bin_width
  }
  method <- match.arg(method, c("peak", "centroid", "cross_cor"))
  est <- switch(
    method,
    peak = .lag_from_peak(profiles$ligand, profiles$receptor, timing_src, timing_tgt),
    centroid = .lag_from_centroid(profiles$ligand, profiles$receptor,
                                  timing_src, timing_tgt),
    cross_cor = .lag_from_cross_cor(
      profiles$ligand, profiles$receptor, timing_src, timing_tgt,
      bin_width, max_lag_bins
    )
  )
  lag_se <- .lag_se_heuristic(est$lag, bin_width, est$n_valid_bins)
  lag_exceeds_threshold <- is.finite(est$lag) && abs(est$lag) > lag_threshold
  lag_null <- if (isTRUE(compute_lag_null)) {
    .lag_permutation_null(
      profiles$ligand, profiles$receptor, timing_src, timing_tgt,
      method, bin_width, max_lag_bins, est$lag, n_lag_perm, lag_null_seed
    )
  } else {
    list(
      p_lag = NA_real_, lag_null_median = NA_real_,
      n_lag_perm_observed = 0L, lag_null_mode = NA_character_
    )
  }
  lag_null_centered <- if (is.finite(est$lag) && is.finite(lag_null$lag_null_median)) {
    est$lag - lag_null$lag_null_median
  } else {
    NA_real_
  }
  key_row <- traj_df[1L, by, drop = FALSE]
  cbind(
    key_row,
    data.frame(
      lag_method = method,
      ligand_time = est$ligand_time,
      receptor_time = est$receptor_time,
      lag = est$lag,
      p_lag = lag_null$p_lag,
      lag_null_median = lag_null$lag_null_median,
      n_lag_perm_observed = lag_null$n_lag_perm_observed,
      lag_null_mode = lag_null$lag_null_mode,
      lag_null_centered = lag_null_centered,
      ligand_peak_time = if ("ligand_peak_time" %in% names(est)) {
        est$ligand_peak_time
      } else {
        NA_real_
      },
      receptor_peak_time = if ("receptor_peak_time" %in% names(est)) {
        est$receptor_peak_time
      } else {
        NA_real_
      },
      lag_peak = if ("lag_peak" %in% names(est)) est$lag_peak else NA_real_,
      ligand_centroid_time = if ("ligand_centroid_time" %in% names(est)) {
        est$ligand_centroid_time
      } else {
        NA_real_
      },
      receptor_centroid_time = if ("receptor_centroid_time" %in% names(est)) {
        est$receptor_centroid_time
      } else {
        NA_real_
      },
      lag_centroid = if ("lag_centroid" %in% names(est)) est$lag_centroid else NA_real_,
      lag_cross_cor = if ("lag_cross_cor" %in% names(est)) est$lag_cross_cor else NA_real_,
      cross_cor_max = if ("cross_cor_max" %in% names(est)) est$cross_cor_max else NA_real_,
      lag_tie = if ("lag_tie" %in% names(est)) est$lag_tie else FALSE,
      n_equivalent_lags = if ("n_equivalent_lags" %in% names(est)) {
        est$n_equivalent_lags
      } else {
        1L
      },
      lag_se = lag_se,
      lag_exceeds_threshold = lag_exceeds_threshold,
      lag_significant = lag_exceeds_threshold,
      lag_threshold = lag_threshold,
      bin_width = bin_width,
      n_valid_bins = est$n_valid_bins,
      lag_profile_source = profiles$lag_profile_source,
      stringsAsFactors = FALSE
    )
  )
}

#' Ligand-receptor timing lag along pseudotime trajectories.
#'
#' Estimates sender-receiver timing separation per trajectory from binned LR
#' scores. Side profiles use `ligand_score` / `receptor_score` columns when
#' present; otherwise gene-level profiles from optional `expr` + `lr_db`, or
#' score weighted by sender/receiver lineage support per bin.
#'
#' @param temporal_scores data.frame with `pair_idx`, `source`, `target`,
#'   `lineage`, `bin`, and `score`.
#' @param lineage_weights numeric matrix `[n_cells × n_lineages]`.
#' @param pseudotime numeric vector of length `n_cells`.
#' @param celltypes integer vector (1-based) of length `n_cells`.
#' @param n_bins number of pseudotime bins (inferred when `NULL`).
#' @param method `"peak"` (separate peak bins per side), `"centroid"`
#'   (score-weighted mean pseudotime), or `"cross_cor"` (discrete shift search).
#' @param lag_threshold minimum `|lag|` for `lag_exceeds_threshold`; defaults to one
#'   bin width. Deprecated alias: `lag_significant`.
#' @param max_lag_bins maximum shift (in bins) for `"cross_cor"`.
#' @param expr optional expression matrix for gene-level side profiles.
#' @param lr_db optional LR database (with `expr`).
#' @param by trajectory key columns.
#' @return data.frame with one row per trajectory: `ligand_peak_time`,
#'   `receptor_peak_time`, `lag_peak` (peak method), centroid / cross-cor
#'   columns, `lag`, `lag_se`, `lag_exceeds_threshold`, and `bin_width`.
#'   Cross-correlation output also reports `lag_tie` and `n_equivalent_lags`;
#'   unresolved equal-strength opposite shifts return a direction-neutral lag 0.
#'   `lag_significant` is a deprecated alias of `lag_exceeds_threshold`.
#' @param compute_lag_null when `TRUE`, run a circular-shift profile null that
#'   preserves each side's temporal shape and autocorrelation.
#' @param n_lag_perm permutation draws for the lag null (default 199).
#' @param lag_null_seed random seed for lag null.
#' @export
ComputeLigandReceptorLag <- function(temporal_scores,
                                   lineage_weights,
                                   pseudotime,
                                   celltypes,
                                   n_bins = NULL,
                                   method = c("peak", "centroid", "cross_cor"),
                                   lag_threshold = NULL,
                                   max_lag_bins = NULL,
                                   expr = NULL,
                                   lr_db = NULL,
                                   compute_lag_null = FALSE,
                                   n_lag_perm = 199L,
                                   lag_null_seed = 41L,
                                   by = c("pair_idx", "source", "target", "lineage")) {
  stopifnot(is.data.frame(temporal_scores))
  req <- c("pair_idx", "source", "target", "lineage", "bin", "score")
  missing <- setdiff(req, colnames(temporal_scores))
  if (length(missing) > 0L) {
    stop("`temporal_scores` must contain columns: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  by <- as.character(by)
  method <- match.arg(method)
  n_cells <- length(pseudotime)
  lw <- .assert_lineage_weights(lineage_weights, n_cells)
  pt <- .assert_pseudotime_1d(pseudotime, n_cells)
  n_bins <- .infer_n_bins(temporal_scores, n_bins)
  if (!is.null(expr)) {
    expr <- .prepare_lr_expr(expr)
  }
  ct_levels_full <- attr(temporal_scores, "celltype_levels")
  timing <- .compute_celltype_lineage_bin_timing(
    celltypes, pt, lw, n_bins, ct_levels_full = ct_levels_full
  )
  ct_levels <- timing$ct_levels
  if (!is.null(ct_levels_full)) {
    ct_levels <- ct_levels_full
  }
  scores_work <- temporal_scores
  scores_work$.traj_key <- .trajectory_key_col(scores_work, by)
  groups <- split(scores_work, scores_work$.traj_key)
  rows <- lapply(seq_along(groups), function(gi) {
    traj_df <- groups[[gi]]
    .compute_trajectory_lag_row(
      traj_df, temporal_scores, timing$timing, timing$wt_sum,
      ct_levels, expr, lr_db, celltypes, pt, lw, n_bins,
      method, lag_threshold, max_lag_bins, by,
      compute_lag_null = compute_lag_null,
      n_lag_perm = n_lag_perm,
      lag_null_seed = lag_null_seed + gi - 1L
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out$.traj_key <- NULL
  out$q_lag <- stats::p.adjust(out$p_lag, method = "BH")
  attr(out, "lag_method") <- method
  attr(out, "n_bins") <- n_bins
  out
}

#' @rdname ComputeLigandReceptorLag
#' @export
ComputeLigandTargetLag <- ComputeLigandReceptorLag

#' Merge trajectory-level lag columns into role features.
#'
#' @param role_features data.frame from [ComputeRoleFeatures()].
#' @param lag_df data.frame from [ComputeLigandReceptorLag()].
#' @param by trajectory key columns.
#' @return `role_features` with lag columns attached (replicated per bin-row).
#' @export
MergeLagIntoRoleFeatures <- function(role_features,
                                     lag_df,
                                     by = c("pair_idx", "source", "target", "lineage")) {
  stopifnot(is.data.frame(role_features), is.data.frame(lag_df))
  by <- as.character(by)
  missing <- setdiff(by, colnames(role_features))
  if (length(missing) > 0L) {
    stop("`role_features` missing trajectory key columns: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  missing_lag <- setdiff(by, colnames(lag_df))
  if (length(missing_lag) > 0L) {
    stop("`lag_df` missing trajectory key columns: ",
         paste(missing_lag, collapse = ", "), call. = FALSE)
  }
  lag_cols <- setdiff(colnames(lag_df), c(by, ".traj_key"))
  key <- .trajectory_key_col(role_features, by)
  lag_key <- .trajectory_key_col(lag_df, by)
  lag_unique <- lag_df[!duplicated(lag_key), , drop = FALSE]
  lag_unique_key <- .trajectory_key_col(lag_unique, by)
  for (col in lag_cols) {
    role_features[[col]] <- lag_unique[[col]][match(key, lag_unique_key)]
  }
  role_features
}

#' Per-donor ligand-receptor lag summary.
#'
#' Computes [ComputeLigandReceptorLag()] independently within each `block_labels`
#' stratum for donor-level meta-analysis.
#'
#' @param temporal_scores pooled temporal LR scores (see [ComputeLigandReceptorLag()]).
#' @param block_labels donor/sample labels of length `n_cells`.
#' @param lineage_weights,pseudotime,celltypes see [ComputeLigandReceptorLag()].
#' @param n_bins,method,lag_threshold,max_lag_bins,expr,lr_db,by passed through.
#' @param compute_lag_null when `TRUE`, run circular-shift profile null per block
#'   (see [ComputeLigandReceptorLag()]).
#' @param n_lag_perm,lag_null_seed passed to [ComputeLigandReceptorLag()].
#' @param min_block_cells minimum cells per block.
#' @return long data.frame with `block_id`, `block_label`, and lag columns per
#'   trajectory × donor.
#' @export
SummarizeLagByDonor <- function(temporal_scores,
                                block_labels,
                                lineage_weights,
                                pseudotime,
                                celltypes,
                                n_bins = NULL,
                                method = c("peak", "centroid", "cross_cor"),
                                lag_threshold = NULL,
                                max_lag_bins = NULL,
                                expr = NULL,
                                lr_db = NULL,
                                by = c("pair_idx", "source", "target", "lineage"),
                                min_block_cells = 10L,
                                compute_lag_null = FALSE,
                                n_lag_perm = 199L,
                                lag_null_seed = 41L) {
  stopifnot(is.data.frame(temporal_scores))
  if (is.null(block_labels)) {
    stop("`block_labels` is required for `SummarizeLagByDonor()`", call. = FALSE)
  }
  n_cells <- length(pseudotime)
  blocks <- .prepare_block_labels(block_labels, n_cells)
  block_levels <- attr(blocks, "levels")
  min_block_cells <- .validate_positive_int(min_block_cells, "min_block_cells", min = 1L)
  method <- match.arg(method)
  if (!is.null(expr)) {
    expr <- .prepare_lr_expr(expr)
  }
  rows <- list()
  row_i <- 0L
  for (b in seq_along(block_levels)) {
    block_cells <- which(blocks == b)
    if (length(block_cells) < min_block_cells) {
      next
    }
    block_lag <- ComputeLigandReceptorLag(
      temporal_scores = temporal_scores,
      lineage_weights = lineage_weights[block_cells, , drop = FALSE],
      pseudotime = pseudotime[block_cells],
      celltypes = celltypes[block_cells],
      n_bins = n_bins,
      method = method,
      lag_threshold = lag_threshold,
      max_lag_bins = max_lag_bins,
      expr = if (is.null(expr)) NULL else expr[, block_cells, drop = FALSE],
      lr_db = lr_db,
      by = by,
      compute_lag_null = compute_lag_null,
      n_lag_perm = n_lag_perm,
      lag_null_seed = lag_null_seed + b - 1L
    )
    if (nrow(block_lag) < 1L) {
      next
    }
    row_i <- row_i + 1L
    block_lag$block_id <- b
    block_lag$block_label <- block_levels[[b]]
    block_lag$n_cells <- length(block_cells)
    rows[[row_i]] <- block_lag
  }
  if (row_i < 1L) {
    stop("No blocks passed `min_block_cells`; check `block_labels`", call. = FALSE)
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

.role_like_column_map <- c(
  instructive = "instructive_like",
  transition_like = "transition_like",
  reinforcing = "reinforcing_like",
  feedback = "feedback_like",
  lateral = "lateral_like",
  bystander = "bystander_like"
)

.wilson_binomial_ci <- function(successes, trials, conf = 0.95) {
  if (!is.finite(successes) || !is.finite(trials) || trials <= 0) {
    return(c(lower = NA_real_, upper = NA_real_))
  }
  z <- stats::qnorm(1 - (1 - conf) / 2)
  phat <- successes / trials
  denom <- 1 + z^2 / trials
  center <- (phat + z^2 / (2 * trials)) / denom
  half <- z * sqrt(phat * (1 - phat) / trials + z^2 / (4 * trials^2)) / denom
  c(
    lower = max(0, center - half),
    upper = min(1, center + half)
  )
}

.resolve_role_timing_delta <- function(role_features, timing_source) {
  timing_source <- match.arg(
    timing_source,
    c("cell_state_tau", "ligand_receptor_lag", "ligand_target_lag")
  )
  if (timing_source == "ligand_target_lag") {
    timing_source <- "ligand_receptor_lag"
  }
  switch(
    timing_source,
    cell_state_tau = {
      if ("delta_tau" %in% colnames(role_features)) {
        role_features$delta_tau
      } else {
        ComputeDeltaTau(role_features$source_timing, role_features$target_timing)
      }
    },
    ligand_receptor_lag = {
      if (all(c("ligand_time", "receptor_time") %in% colnames(role_features))) {
        role_features$receptor_time - role_features$ligand_time
      } else if ("lag_peak" %in% colnames(role_features)) {
        role_features$lag_peak
      } else {
        warning(
          "`timing_source = \"ligand_receptor_lag\"` requires ligand/receptor ",
          "timing columns; falling back to cell_state_tau.",
          call. = FALSE
        )
        role_features$delta_tau
      }
    },
    ligand_target_lag = {
      if ("lag" %in% colnames(role_features)) {
        role_features$lag
      } else {
        warning(
          "`timing_source = \"ligand_target_lag\"` requires `lag`; ",
          "falling back to cell_state_tau.",
          call. = FALSE
        )
        role_features$delta_tau
      }
    }
  )
}

.apply_role_timing_source <- function(role_features, timing_source) {
  timing_source <- match.arg(
    timing_source,
    c("cell_state_tau", "ligand_receptor_lag", "ligand_target_lag")
  )
  if (timing_source == "ligand_target_lag") {
    timing_source <- "ligand_receptor_lag"
  }
  if (timing_source == "cell_state_tau") {
    return(role_features)
  }
  delta <- .resolve_role_timing_delta(role_features, timing_source)
  role_features$delta_tau <- delta
  if ("signed_lag" %in% colnames(role_features)) {
    role_features$signed_lag <- delta
  }
  attr(role_features, "timing_source") <- timing_source
  role_features
}

.role_formal_eligible <- function(role_features) {
  supported <- role_features$fully_supported %in% TRUE
  supported <- supported &
    is.finite(role_features$spec_self) &
    is.finite(role_features$spec_max_other) &
    is.finite(role_features$delta_tau)
  if ("structural_valid" %in% colnames(role_features)) {
    supported <- supported & role_features$structural_valid %in% TRUE
  }
  if ("n_evaluable_lineages" %in% colnames(role_features)) {
    supported <- supported &
      role_features$n_evaluable_lineages >= 2L &
      is.finite(role_features$n_evaluable_lineages)
  }
  supported
}

.dominant_role_like <- function(scores_mat) {
  if (!nrow(scores_mat)) {
    return(character(0))
  }
  main_cols <- c("instructive_like", "transition_like", "reinforcing_like",
                 "feedback_like")
  aux_cols <- c("lateral_like", "bystander_like")
  all_cols <- c(main_cols, aux_cols)
  scores_mat <- scores_mat[, intersect(all_cols, colnames(scores_mat)), drop = FALSE]
  apply(scores_mat, 1L, function(row) {
    row[!is.finite(row)] <- -Inf
    if (all(!is.finite(row) | row == -Inf)) {
      return(NA_character_)
    }
    main <- row[intersect(main_cols, names(row))]
    aux <- row[intersect(aux_cols, names(row))]
    if (length(main) && any(is.finite(main)) && max(main, na.rm = TRUE) > 0) {
      return(names(main)[which.max(main)])
    }
    if (length(aux) && any(is.finite(aux)) && max(aux, na.rm = TRUE) > 0) {
      return(names(aux)[which.max(aux)])
    }
    names(row)[which.max(row)]
  })
}

.role_main_like_cols <- function() {
  c(
    "instructive_like", "transition_like", "reinforcing_like", "feedback_like",
    "lateral_like", "bystander_like"
  )
}

.resolve_role_bif_df <- function(role_features) {
  bif_df <- attr(role_features, "bif_df")
  if (is.null(bif_df)) {
    if (!all(c("b_star", "bif_position") %in% colnames(role_features))) {
      stop("`role_features` must carry attribute `bif_df` or bifurcation columns",
           call. = FALSE)
    }
    bif_df <- unique(role_features[, c("lineage", "b_star"), drop = FALSE])
    if (!"reason" %in% colnames(role_features)) {
      bif_df$reason <- ifelse(is.na(bif_df$b_star),
                              "c2_not_informative", "crossing")
    }
  }
  bif_df
}

.resolve_role_tau_settings <- function(role_features, tau_0 = NULL) {
  tau_by_lin <- attr(role_features, "tau_0_by_lineage")
  use_scalar_tau <- !is.null(tau_0)
  if (use_scalar_tau) {
    list(use_scalar_tau = TRUE, tau_scalar = .validate_positive_num(tau_0, "tau_0"))
  } else if (is.null(tau_by_lin)) {
    list(use_scalar_tau = TRUE, tau_scalar = 1.0)
  } else {
    list(use_scalar_tau = FALSE, tau_by_lin = tau_by_lin)
  }
}

.compute_role_like_scores <- function(role_features,
                                     bif_df,
                                     tau_settings,
                                     h_fork = 0.5,
                                     min_n_eff_cells = 5.0,
                                     transition_window = 0L,
                                     timing_source = NA_character_) {
  scores_df <- role_features[, c("bin", "lineage", "n_eff_cells"), drop = FALSE]
  raw_scores <- data.frame(
    instructive = rep(NA_real_, nrow(role_features)),
    transition_like = rep(NA_real_, nrow(role_features)),
    reinforcing = rep(NA_real_, nrow(role_features)),
    feedback = rep(NA_real_, nrow(role_features)),
    lateral = rep(NA_real_, nrow(role_features)),
    bystander = rep(NA_real_, nrow(role_features)),
    stringsAsFactors = FALSE
  )
  lin_levels <- unique(.as_integerish(role_features$lineage, "role_features$lineage"))
  tau_0_by_lineage <- stats::setNames(rep(NA_real_, length(lin_levels)), lin_levels)
  for (lin in lin_levels) {
    idx <- which(.as_integerish(role_features$lineage, "role_features$lineage") == lin)
    lin_tau <- if (tau_settings$use_scalar_tau) {
      tau_settings$tau_scalar
    } else {
      .validate_positive_num(
        tau_settings$tau_by_lin[lin], sprintf("tau_0_by_lineage[%d]", lin)
      )
    }
    valid_idx <- idx[
      is.finite(role_features$spec_self[idx]) &
        is.finite(role_features$spec_max_other[idx]) &
        is.finite(role_features$delta_tau[idx]) &
        is.finite(scores_df$n_eff_cells[idx])
    ]
    if (length(valid_idx)) {
      part <- ComputeRoleScores(
        scores_df[valid_idx, , drop = FALSE],
        spec_self = role_features$spec_self[valid_idx],
        spec_max_other = role_features$spec_max_other[valid_idx],
        delta_tau = role_features$delta_tau[valid_idx],
        bif_df = bif_df,
        tau_0 = lin_tau,
        h_fork = h_fork,
        min_n_eff_cells = min_n_eff_cells,
        transition_window = transition_window
      )
      raw_scores[valid_idx, ] <- part
    }
    tau_0_by_lineage[as.character(lin)] <- lin_tau
  }
  like_scores <- raw_scores
  names(like_scores) <- unname(.role_like_column_map[names(like_scores)])
  evidence_cols <- intersect(
    c("analysis_data_hash", "scoring_config_hash", "fake_null_config_hash",
      "temporal_null_config_hash", "adaptive_perm_config_hash"),
    colnames(role_features)
  )
  evidence_hashes <- lapply(
    role_features[, evidence_cols, drop = FALSE],
    function(x) unique(as.character(x))
  )
  role_config_hash <- .compute_role_config_hash(
    tau_0_by_lineage, h_fork, min_n_eff_cells, transition_window,
    timing_source = timing_source,
    bif_df = bif_df,
    evidence_hashes = evidence_hashes
  )
  list(like_scores = like_scores, role_config_hash = role_config_hash)
}

.trajectory_max_like_scores <- function(like_scores, traj_key) {
  main_cols <- .role_main_like_cols()
  split_idx <- split(seq_len(nrow(like_scores)), traj_key)
  vapply(split_idx, function(idx) {
    mat <- as.matrix(like_scores[idx, main_cols, drop = FALSE])
    colMaxs <- apply(mat, 2L, function(col) {
      if (!any(is.finite(col))) {
        return(-Inf)
      }
      max(col, na.rm = TRUE)
    })
    colMaxs
  }, numeric(length(main_cols)))
}

.trajectory_dominant_role <- function(max_like_mat) {
  if (!length(max_like_mat)) {
    return(character(0))
  }
  if (is.null(dim(max_like_mat))) {
    max_like_mat <- matrix(max_like_mat, nrow = 1L)
  }
  main_cols <- .role_main_like_cols()
  colnames(max_like_mat) <- main_cols
  .dominant_role_like(max_like_mat)
}

.perturb_role_features_cell_boot <- function(role_features,
                                             cell_idx,
                                             pseudotime,
                                             lineage_weights,
                                             timing_noise_scale = 0.05) {
  out <- role_features
  if (timing_noise_scale <= 0) {
    return(out)
  }
  n_lineages <- ncol(lineage_weights)
  lin_levels <- unique(.as_integerish(out$lineage, "role_features$lineage"))
  pt_boot <- pseudotime[cell_idx]
  lw_boot <- lineage_weights[cell_idx, , drop = FALSE]
  for (lin in lin_levels) {
    w <- lw_boot[, lin]
    ok <- is.finite(pt_boot) & w > 0
    jitter <- if (sum(ok) >= 2L) {
      stats::sd(pt_boot[ok], na.rm = TRUE) * timing_noise_scale
    } else {
      timing_noise_scale * 0.01
    }
    idx <- which(.as_integerish(out$lineage, "role_features$lineage") == lin)
    if (length(idx)) {
      out$delta_tau[idx] <- out$delta_tau[idx] +
        stats::rnorm(length(idx), sd = jitter)
      if ("signed_lag" %in% colnames(out)) {
        out$signed_lag[idx] <- out$delta_tau[idx]
      }
    }
  }
  out
}

.resample_block_cells_with_replacement <- function(block_labels, n_cells) {
  blocks <- .prepare_block_labels(block_labels, n_cells)
  block_levels <- attr(blocks, "levels")
  n_blocks <- length(block_levels)
  if (n_blocks < 1L) {
    return(sample.int(n_cells, n_cells, replace = TRUE))
  }
  draw_blocks <- sample.int(n_blocks, n_blocks, replace = TRUE)
  unlist(lapply(draw_blocks, function(b) which(blocks == b)), use.names = FALSE)
}

.resample_cell_indices <- function(n_cells, block_labels = NULL, n_draw = NULL) {
  n_draw <- if (is.null(n_draw)) n_cells else n_draw
  if (is.null(block_labels)) {
    return(sample.int(n_cells, n_draw, replace = TRUE))
  }
  blocks <- .prepare_block_labels(block_labels, n_cells)
  block_levels <- attr(blocks, "levels")
  if (length(block_levels) < 2L) {
    return(sample.int(n_cells, n_draw, replace = TRUE))
  }
  draw_blocks <- sample.int(length(block_levels), n_draw, replace = TRUE)
  unlist(lapply(draw_blocks, function(b) {
    pool <- which(blocks == b)
    if (!length(pool)) {
      return(sample.int(n_cells, 1L, replace = TRUE))
    }
    sample(pool, 1L, replace = TRUE)
  }), use.names = FALSE)
}

#' End-to-end soft role inference with formal/exploratory gating.
#'
#' Wraps [ComputeRoleScores()] using features from [ComputeRoleFeatures()].
#' Rows on trajectories with `fully_supported == TRUE` receive formal
#' `*_like` role scores and a `dominant_role_like` label; all other rows are
#' marked `role_inference_status = "exploratory"` with `NA` formal role outputs.
#'
#' @param role_features data.frame from `ComputeRoleFeatures()`.
#' @param merged_evidence optional `MergeEvidence()` table; when supplied,
#'   overrides / supplies trajectory-level `fully_supported`.
#' @param tau_0 numeric sigmoid scale; defaults to per-lineage
#'   `0.1 * sd(pseudotime)` stored in `attr(role_features, "tau_0_by_lineage")`.
#' @param h_fork bifurcation threshold recorded with the role-scoring settings;
#'   row-level fork positions are supplied by `role_features`.
#' @param min_n_eff_cells row-level reliability gate (see `ComputeRoleScores()`).
#' @param transition_window non-negative integer; bins with
#'   `abs(bin - b_star) <= transition_window` receive non-zero `transition_like`
#'   in the Rust core and `in_transition_window = TRUE` in the output.
#' @param timing_source timing signal wired into role scoring (see
#'   [ComputeRoleFeatures()]).
#' @param by trajectory key columns.
#' @return data.frame with `*_like` score columns, `role_inference_status`,
#'   `dominant_role_like`, `fully_supported`, `n_total_lineages`,
#'   `n_evaluable_lineages`, `n_supported_lineages`, `in_transition_window`, and
#'   calibration columns when available. Attribute `role_config_hash` fingerprints
#'   role parameters; attribute `timing_source` records the timing signal used.
#' @export
ComputeRoleInference <- function(role_features,
                                 merged_evidence = NULL,
                                 tau_0 = NULL,
                                 h_fork = 0.5,
                                 min_n_eff_cells = 5.0,
                                 transition_window = 0L,
                                 timing_source = c("cell_state_tau",
                                                   "ligand_receptor_lag",
                                                   "ligand_target_lag"),
                                 by = c("pair_idx", "source", "target", "lineage")) {
  stopifnot(is.data.frame(role_features))
  req <- c("spec_self", "spec_max_other", "delta_tau", "bin", "lineage",
           "n_eff_cells")
  missing <- setdiff(req, colnames(role_features))
  if (length(missing) > 0L) {
    stop("`role_features` must contain columns: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  by <- as.character(by)
  timing_source <- match.arg(timing_source)
  transition_window <- .validate_positive_int(transition_window, "transition_window",
                                                min = 0L)
  role_features_scoring <- .apply_role_timing_source(role_features, timing_source)
  timing_source <- if (is.null(attr(role_features_scoring, "timing_source"))) {
    timing_source
  } else {
    attr(role_features_scoring, "timing_source")
  }
  bif_df <- .resolve_role_bif_df(role_features_scoring)
  tau_settings <- .resolve_role_tau_settings(role_features_scoring, tau_0)
  scored <- .compute_role_like_scores(
    role_features_scoring, bif_df, tau_settings, h_fork, min_n_eff_cells,
    transition_window = transition_window,
    timing_source = timing_source
  )
  like_scores <- scored$like_scores
  role_config_hash <- scored$role_config_hash

  out <- cbind(role_features_scoring, like_scores)
  out <- .attach_fully_supported(out, merged_evidence, by)
  if (is.null(merged_evidence) && "fully_supported" %in% colnames(role_features)) {
    out$fully_supported <- role_features$fully_supported
  }

  b_star_row <- out$b_star
  if (is.null(b_star_row) && "b_star" %in% colnames(role_features)) {
    b_star_row <- role_features$b_star
  }
  if (!is.null(b_star_row)) {
    out$in_transition_window <- is.finite(b_star_row) &
      abs(as.numeric(out$bin) - as.numeric(b_star_row)) <= transition_window
  } else {
    out$in_transition_window <- rep(FALSE, nrow(out))
  }

  formal <- .role_formal_eligible(out)
  out$role_inference_status <- ifelse(formal, "formal", "exploratory")
  like_cols <- unname(.role_like_column_map)
  formal_mat <- as.matrix(out[, like_cols, drop = FALSE])
  formal_mat[!formal, ] <- NA_real_
  out[, like_cols] <- formal_mat
  out$dominant_role_like <- NA_character_
  if (any(formal)) {
    out$dominant_role_like[formal] <- .dominant_role_like(formal_mat[formal, , drop = FALSE])
  }

  wanted_cal <- c(
    "q_joint", "p_joint", "q_identity", "p_identity",
    "q_temporal", "p_temporal"
  )
  missing_cal <- setdiff(wanted_cal, colnames(out))
  if (length(missing_cal) && !is.null(merged_evidence)) {
    keep <- c(by, intersect(missing_cal, colnames(merged_evidence)))
    if (length(setdiff(keep, by))) {
      cal <- merged_evidence[, keep, drop = FALSE]
      cal <- cal[!duplicated(.trajectory_key_col(cal, by)), , drop = FALSE]
      key <- .trajectory_key_col(out, by)
      cal_key <- .trajectory_key_col(cal, by)
      for (col in setdiff(keep, by)) {
        out[[col]] <- cal[[col]][match(key, cal_key)]
      }
    }
  }

  attr(out, "role_config_hash") <- role_config_hash
  attr(out, "timing_source") <- timing_source
  out
}

# Rebuild role features after block/cell resampling.
.rebuild_role_features_block_boot <- function(rebuild_context, cell_idx, by,
                                               timing_source) {
  required <- c("expr", "celltypes", "lr_db", "pseudotime", "lineage_weights")
  missing <- setdiff(required, names(rebuild_context))
  if (length(missing) > 0L) {
    stop(
      "`rebuild_context` missing required fields: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  expr_b <- rebuild_context$expr[, cell_idx, drop = FALSE]
  ct_b <- rebuild_context$celltypes[cell_idx]
  pt_b <- rebuild_context$pseudotime[cell_idx]
  lw_b <- rebuild_context$lineage_weights[cell_idx, , drop = FALSE]
  n_bins <- rebuild_context$n_bins
  if (is.null(n_bins)) {
    n_bins <- 12L
  }
  tl_extra <- rebuild_context$temporal_lr_args
  if (is.null(tl_extra)) {
    tl_extra <- list()
  }
  use_intermediate <- isTRUE(rebuild_context$use_intermediate)
  if (use_intermediate) {
    res_b <- do.call(
      ComputeTemporalLRWithIntermediate,
      c(
        list(
          expr = expr_b,
          celltypes = ct_b,
          lr_db = rebuild_context$lr_db,
          pseudotime = pt_b,
          lineage_weights = lw_b,
          n_bins = n_bins
        ),
        tl_extra
      )
    )
    scores_b <- res_b$scores
    intermediate_b <- res_b$intermediate
  } else {
    scores_b <- do.call(
      ComputeTemporalLR,
      c(
        list(
          expr = expr_b,
          celltypes = ct_b,
          lr_db = rebuild_context$lr_db,
          pseudotime = pt_b,
          lineage_weights = lw_b,
          n_bins = n_bins
        ),
        tl_extra
      )
    )
    intermediate_b <- rebuild_context$intermediate
  }
  rf_extra <- rebuild_context$role_feature_args
  if (is.null(rf_extra)) {
    rf_extra <- list()
  }
  rf_base <- list(
    temporal_scores = scores_b,
    lineage_weights = lw_b,
    pseudotime = pt_b,
    celltypes = ct_b,
    by = by,
    timing_source = timing_source
  )
  if (!is.null(intermediate_b)) {
    rf_base$intermediate <- intermediate_b
  }
  if (!is.null(rebuild_context$evidence)) {
    rf_base$evidence <- rebuild_context$evidence
  }
  do.call(ComputeRoleFeatures, c(rf_base, rf_extra))
}

#' Bootstrap sensitivity and stability for trajectory-level role assignment.
#'
#' Resamples uncertainty in role assignment and reports per-trajectory stability
#' metrics (modal role, win rates, Wilson CIs). This quantifies label stability
#' under the documented resampling protocol; it does not replace FakeLR or
#' temporal permutation p-values for LR identity or dynamics. Block rebuilds
#' reuse the original evidence table, so they estimate conditional role
#' stability after discovery rather than full-pipeline uncertainty.
#'
#' @param role_features data.frame from `ComputeRoleFeatures()`.
#' @param merged_evidence optional `MergeEvidence()` table for
#'   `fully_supported` gating.
#' @param n_boot number of bootstrap replicates (default 500).
#' @param seed random seed.
#' @param tau_0,h_fork,min_n_eff_cells,by,timing_source passed to role scoring
#'   (see `ComputeRoleInference()`).
#' @param transition_window passed to `ComputeRoleInference()` for the point
#'   estimate and `in_transition_window` column.
#' @param bootstrap_mode `"row"`, `"cell"`, or `"block"`.
#' @param block_labels optional length-`n_cells` donor/sample labels for
#'   block-stratified resampling in `"cell"` and `"block"` modes.
#' @param pseudotime,lineage_weights required for `"cell"` mode.
#' @param rebuild_context named list for `"block"` mode with `expr`, `celltypes`,
#'   `lr_db`, `pseudotime`, `lineage_weights`, and optional `n_bins`,
#'   `temporal_lr_args`, `role_feature_args` passed to
#'   [ComputeTemporalLR()] / [ComputeRoleFeatures()].
#' @param timing_noise_scale scale for `delta_tau` perturbation in cell mode.
#' @param ci_level two-sided confidence level for Wilson binomial CIs (default
#'   0.95).
#' @return list with `inference` (point estimates) and `stability` (one row per
#'   trajectory with modal-role frequencies, Wilson intervals, and an explicit
#'   sensitivity/conditional-stability status). Attributes record the bootstrap
#'   settings.
#' @export
BootstrapRoleInference <- function(role_features,
                                   merged_evidence = NULL,
                                   n_boot = 500L,
                                   seed = 12345L,
                                   tau_0 = NULL,
                                   h_fork = 0.5,
                                   min_n_eff_cells = 5.0,
                                   transition_window = 0L,
                                   timing_source = c("cell_state_tau",
                                                     "ligand_receptor_lag",
                                                     "ligand_target_lag"),
                                   by = c("pair_idx", "source", "target", "lineage"),
                                   bootstrap_mode = c("row", "cell", "block"),
                                   block_labels = NULL,
                                   pseudotime = NULL,
                                   lineage_weights = NULL,
                                   rebuild_context = NULL,
                                   timing_noise_scale = 0.05,
                                   ci_level = 0.95) {
  stopifnot(is.data.frame(role_features))
  by <- as.character(by)
  timing_source <- match.arg(timing_source)
  bootstrap_mode <- match.arg(bootstrap_mode)
  n_boot <- .validate_positive_int(n_boot, "n_boot", min = 1L)
  seed <- .validate_seed(seed)
  ci_level <- .validate_unit_interval(ci_level, "ci_level")
  if (ci_level <= 0 || ci_level >= 1) {
    stop("`ci_level` must be strictly between 0 and 1", call. = FALSE)
  }
  if (bootstrap_mode == "block") {
    if (is.null(block_labels)) {
      stop("`block_labels` is required when `bootstrap_mode = \"block\"`",
           call. = FALSE)
    }
    if (is.null(rebuild_context)) {
      stop("`rebuild_context` is required when `bootstrap_mode = \"block\"`",
           call. = FALSE)
    }
  }

  inference <- ComputeRoleInference(
    role_features,
    merged_evidence = merged_evidence,
    tau_0 = tau_0,
    h_fork = h_fork,
    min_n_eff_cells = min_n_eff_cells,
    transition_window = transition_window,
    timing_source = timing_source,
    by = by
  )

  role_features_scoring <- .apply_role_timing_source(role_features, timing_source)
  timing_source <- if (is.null(attr(role_features_scoring, "timing_source"))) {
    timing_source
  } else {
    attr(role_features_scoring, "timing_source")
  }
  bif_df <- .resolve_role_bif_df(role_features_scoring)
  tau_settings <- .resolve_role_tau_settings(role_features_scoring, tau_0)

  traj_key <- .trajectory_key_col(inference, by)
  traj_keys <- unique(traj_key)
  formal_mask <- .role_formal_eligible(inference)
  formal_keys <- unique(traj_key[formal_mask])

  n_cells <- NULL
  if (bootstrap_mode == "cell") {
    if (is.null(pseudotime) || is.null(lineage_weights)) {
      stop("`pseudotime` and `lineage_weights` are required when ",
           "`bootstrap_mode = \"cell\"`", call. = FALSE)
    }
    n_cells <- length(pseudotime)
    lw <- .assert_lineage_weights(lineage_weights, n_cells)
    pt <- .assert_pseudotime_1d(pseudotime, n_cells)
    if (!is.null(block_labels)) {
      .prepare_block_labels(block_labels, n_cells)
    }
  } else if (bootstrap_mode == "block") {
    n_cells <- length(rebuild_context$pseudotime)
    .prepare_block_labels(block_labels, n_cells)
  }

  main_cols <- .role_main_like_cols()

  stability_rows <- vector("list", length(traj_keys))
  set.seed(seed)

  for (ti in seq_along(traj_keys)) {
    key <- traj_keys[ti]
    row_idx <- which(traj_key == key)
    g <- role_features_scoring[row_idx, , drop = FALSE]
    is_formal <- key %in% formal_keys

    if (!is_formal) {
      stability_rows[[ti]] <- cbind(
        inference[row_idx[1L], by, drop = FALSE],
        data.frame(
          modal_role = NA_character_,
          dominant_role_like = NA_character_,
          role_win_frequency = I(list(NULL)),
          role_bootstrap_fraction = I(list(NULL)),
          role_stability = NA_real_,
          role_dominant_wilson_ci = I(list(NULL)),
          role_ci = I(list(NULL)),
          role_bootstrap_status = "exploratory",
          n_boot = n_boot,
          n_boot_valid = NA_integer_,
          stringsAsFactors = FALSE
        )
      )
      next
    }

    boot_dom <- character(n_boot)
    boot_frac_mat <- matrix(0, nrow = n_boot, ncol = length(main_cols),
                            dimnames = list(NULL, main_cols))

    for (b in seq_len(n_boot)) {
      g_boot <- if (bootstrap_mode == "block") {
        cell_idx <- .resample_block_cells_with_replacement(block_labels, n_cells)
        feats_boot <- .rebuild_role_features_block_boot(
          rebuild_context, cell_idx, by, timing_source
        )
        key_match <- .trajectory_key_col(feats_boot, by) == key
        feats_boot[key_match, , drop = FALSE]
      } else {
        boot_idx <- sample.int(nrow(g), nrow(g), replace = TRUE)
        g_boot <- g[boot_idx, , drop = FALSE]
        if (bootstrap_mode == "cell") {
          cell_idx <- .resample_cell_indices(n_cells, block_labels)
          g_boot <- .perturb_role_features_cell_boot(
            g_boot, cell_idx, pt, lw, timing_noise_scale
          )
        }
        g_boot
      }

      if (nrow(g_boot) < 1L) {
        boot_dom[b] <- NA_character_
        next
      }

      scored <- .compute_role_like_scores(
        g_boot, bif_df, tau_settings, h_fork, min_n_eff_cells,
        transition_window = transition_window,
        timing_source = timing_source
      )
      max_like <- .trajectory_max_like_scores(scored$like_scores, rep(1L, nrow(g_boot)))
      if (is.null(dim(max_like))) {
        max_like <- matrix(max_like, nrow = 1L, dimnames = list(NULL, main_cols))
      }
      boot_dom[b] <- .trajectory_dominant_role(max_like)
      for (j in seq_along(main_cols)) {
        boot_frac_mat[b, j] <- as.numeric(boot_dom[b] == main_cols[j])
      }
    }

    valid_boot <- !is.na(boot_dom)
    n_valid_boot <- sum(valid_boot)
    if (n_valid_boot < 1L) {
      stability_rows[[ti]] <- cbind(
        inference[row_idx[1L], by, drop = FALSE],
        data.frame(
          modal_role = NA_character_,
          dominant_role_like = NA_character_,
          role_win_frequency = I(list(NULL)),
          role_bootstrap_fraction = I(list(NULL)),
          role_stability = NA_real_,
          role_dominant_wilson_ci = I(list(NULL)),
          role_ci = I(list(NULL)),
          role_bootstrap_status = "insufficient_boot",
          n_boot = n_boot,
          n_boot_valid = 0L,
          stringsAsFactors = FALSE
        )
      )
      next
    }

    frac_named <- stats::setNames(
      colMeans(boot_frac_mat[valid_boot, , drop = FALSE]),
      main_cols
    )
    role_stability <- max(frac_named, na.rm = TRUE)
    modal_role <- main_cols[which.max(frac_named)]
    modal_wins <- sum(boot_dom[valid_boot] == modal_role, na.rm = TRUE)
    dominant_wilson <- .wilson_binomial_ci(modal_wins, n_valid_boot, conf = ci_level)
    ci_named <- stats::setNames(
      lapply(main_cols, function(col) {
        wins <- sum(boot_frac_mat[valid_boot, col], na.rm = TRUE)
        stats::setNames(
          .wilson_binomial_ci(wins, n_valid_boot, conf = ci_level),
          c("lower", "upper")
        )
      }),
      main_cols
    )

    role_stability_mode <- if (bootstrap_mode == "row") {
      "sensitivity"
    } else if (bootstrap_mode == "cell") {
      "timing_sensitivity"
    } else if (isTRUE(rebuild_context$use_intermediate)) {
      "conditional_role_stability"
    } else {
      "rebuild_sensitivity"
    }
    boot_status <- if (n_valid_boot < 2L) {
      "insufficient_boot"
    } else {
      role_stability_mode
    }

    stability_rows[[ti]] <- cbind(
      inference[row_idx[1L], by, drop = FALSE],
      data.frame(
        modal_role = modal_role,
        dominant_role_like = modal_role,
        role_win_frequency = I(list(frac_named)),
        role_bootstrap_fraction = I(list(frac_named)),
        role_stability = role_stability,
        role_dominant_wilson_ci = I(list(dominant_wilson)),
        role_ci = I(list(ci_named)),
        role_bootstrap_status = boot_status,
        n_boot = n_boot,
        n_boot_valid = n_valid_boot,
        stringsAsFactors = FALSE
      )
    )
  }

  stability <- do.call(rbind, stability_rows)
  rownames(stability) <- NULL

  out <- list(
    inference = inference,
    stability = stability
  )
  attr(out, "n_boot") <- n_boot
  attr(out, "bootstrap_mode") <- bootstrap_mode
  attr(out, "seed") <- seed
  attr(out, "role_stability_mode") <- if (bootstrap_mode == "row") {
    "sensitivity"
  } else if (bootstrap_mode == "cell") {
    "timing_sensitivity"
  } else if (isTRUE(rebuild_context$use_intermediate)) {
    "conditional_role_stability"
  } else {
    "rebuild_sensitivity"
  }
  out
}

#' Per-gene Spearman correlation of expression vs. pseudotime, BH-corrected.
#'
#' @param expr_dense numeric matrix \[n_genes × n_cells\] (raw or
#'   normalised counts).
#' @param pseudotime numeric vector of length n_cells; NaN cells excluded.
#' @param min_cells integer, minimum qualifying cells per gene.
#' @return data.frame with `rho`, `pval`, `qval` (length n_genes); genes
#'   with fewer than `min_cells` qualifying cells get rho=0, pval=1, qval=1.
#' @export
TrajectoryDE <- function(expr_dense, pseudotime, min_cells = 20L) {
  if (!is.numeric(expr_dense)) {
    stop("`expr_dense` must be numeric", call. = FALSE)
  }
  m <- as.matrix(expr_dense)
  n_genes <- nrow(m)
  n_cells <- ncol(m)
  if (n_genes < 1L || n_cells < 3L || !all(is.finite(m)) || any(m < 0)) {
    stop("`expr_dense` must be a non-empty finite non-negative matrix with at least three cells",
         call. = FALSE)
  }
  pseudotime <- .assert_pseudotime_1d(pseudotime, n_cells)
  min_cells <- .validate_positive_int(min_cells, "min_cells", min = 3L)
  if (min_cells > n_cells) {
    stop("`min_cells` cannot exceed ncol(expr_dense)", call. = FALSE)
  }
  out <- rs_trajectory_de(
    as.numeric(t(m)),
    pseudotime,
    as.integer(n_genes),
    as.integer(n_cells),
    as.integer(min_cells)
  )
  data.frame(out, stringsAsFactors = FALSE)
}

# ───  Fake LR null model ────────────────────────────────────────────────

#' Construct and score fake LR pairs for null-model calibration.
#'
#' Identifies candidate
#' genes NOT in any known LR database, matches them to real ligand/receptor
#' genes via stratified sampling on marginal expression features, generates
#' fake LR pairs, scores them, and returns the null distribution.
#'
#' @param expr `dgCMatrix` (genes × cells), log-normalised.
#' @param celltypes integer vector (1-based).
#' @param lr_db list with `ligand_rows` and `receptor_rows`
#'   (see `ComputeLRScores`).
#' @param candidate_fake_genes integer vector of 1-based gene row indices
#'   that are **eligible** for fake pairing.  Caller MUST pre-filter to
#'   exclude genes appearing in any LR database (CellChat, CellPhoneDB,
#'   OmniPath, NicheNet, ICELLNET, Ramilowski) AND STRING-confidence
#'   interactors AND BioGRID physical interactors.
#' @param pseudotime numeric vector.
#' @param lineage_weights numeric matrix \[n_cells × n_lineages\].
#' @param n_bins,n_fake_per_real,n_strata,caliper_k integers.
#' @param max_match_distance maximum normalised four-feature Euclidean
#'   distance allowed for either side of a fake pair. Candidates outside this
#'   numeric caliper are unmatched rather than silently accepted. Must be finite.
#' @param fake_null_mode fake-null construction:
#'   `"context_generic"` (matched gene-pair draws),
#'   `"edge_rewire"` (one-endpoint-preserving matched replacement; the legacy
#'   option name is retained for API stability), `"both"` (generic + endpoint
#'   replacement), `"lr_network_rewire"` (optional third null layer:
#'   replace one endpoint from lr_db ligand/receptor gene pools while preserving
#'   context matching), or `"all"` / `"triple"` (all three layers with budget
#'   split across generic, edge, and network). Primary `p_fake` in `"all"` mode
#'   follows `p_network_iu > p_edge_iu > p_fake_generic`. `"both"` does not
#'   include network rewiring; use `"lr_network_rewire"` or `"all"` when
#'   network-constrained nulls are needed.
#' @param complex_null_mode how multi-subunit LR pairs enter FakeLR:
#'   `"skip"` (default; complex pairs receive NA identity p-values) or
#'   `"cardinality_matched"` (generic null draws fake subunit sets with the
#'   same cardinality as the real complex; max 8 subunits per side).
#' @param min_pct numeric.
#' @param seed numeric.
#' @param min_n_eff minimum Kish effective sample size required for a
#'   structurally valid bin.
#' @param min_valid_bins minimum number of structurally valid bins required
#'   for a trajectory.
#' @param min_valid_fraction minimum fraction of structurally valid bins
#'   required for a trajectory.
#' @param return_bin_scores logical; when `FALSE`, skip the potentially large
#'   fake bin-score table while retaining trajectory calibration.
#' @param trajectory_stat trajectory statistic used for the streamed
#'   calibration table.
#' @param top_k number of bins used by `"topk_mean"`.
#' @param min_fake_observed minimum observed fake trajectories per real pair.
#' @param min_fake_coverage minimum `n_fake_observed / n_fake_expected`.
#' @param min_fake_observed_generic,min_fake_observed_edge minimum observed
#'   trajectories required separately for the generic and endpoint-replacement
#'   nulls.
#' @param min_fake_coverage_generic,min_fake_coverage_edge minimum coverage
#'   required separately for the two null families.
#' @param min_unique_fake_draws minimum number of distinct fake LR pairs needed
#'   within the selected null family.
#' @param min_unique_fake_fraction minimum distinct-pair fraction among
#'   observed fake trajectories within the selected null family.
#' @param block_labels optional donor/sample/batch labels. Context matching
#'   for fake draws is performed within each block when multiple blocks exist.
#' @return list with elements:
#'   - `fake_scores`: data.frame (same columns as `ComputeTemporalLR`)
#'   - `fake_pairs`: data.frame of `ligand_gene`, `receptor_gene`,
#'     `shadowing_pair` (1-based), and per-side match distances
#'   - `trajectory_calibration`: trajectory-level FakeLR table
#'   - `matching_audit`: matching mode, coverage, uniqueness, and distance audit
#'   - `n_fake_total`: integer, total fake pairs sampled
#'   - `candidate_pool_size`: integer, number of eligible fake genes
#'   - `n_real_pairs`: integer, |D|, number of real LR pairs
#'   - `analysis_data_hash`: dataset fingerprint (expr, cell types, LR db, pt, lw, blocks)
#'   - `scoring_config_hash`: shared scoring/evaluation configuration (alias
#'     `calibration_config_hash` for backward compatibility)
#'   - `fake_null_config_hash`: FakeLR null-construction configuration
#'   - `hill_scale`: dataset-level Hill saturation scale
#' @export
ComputeFakeLRNull <- function(expr, celltypes, lr_db,
                               candidate_fake_genes,
                               pseudotime, lineage_weights,
                               n_bins = 12L, n_fake_per_real = 10L,
                               n_strata = 5L, caliper_k = 10L,
                               max_match_distance = 0.5,
                               fake_null_mode = c("context_generic", "edge_rewire", "both",
                                                 "lr_network_rewire", "all", "triple"),
                               min_pct = 0.05, seed = 12345L,
                               min_n_eff = 5.0,
                               min_valid_bins = 3L,
                               min_valid_fraction = 0.25,
                               return_bin_scores = TRUE,
                               trajectory_stat = c("amplitude", "max", "mean_score", "auc", "topk_mean"),
                               top_k = 3L,
                               min_fake_observed = 100L,
                               min_fake_coverage = 0.8,
                               block_labels = NULL,
                               min_fake_observed_generic = min_fake_observed,
                               min_fake_coverage_generic = min_fake_coverage,
                               min_fake_observed_edge = min_fake_observed,
                               min_fake_coverage_edge = min_fake_coverage,
                               min_unique_fake_draws = 10L,
                               min_unique_fake_fraction = 0.01,
                               complex_null_mode = c("skip", "cardinality_matched")) {
  stat_info <- .normalize_trajectory_stat(trajectory_stat)
  stat_code <- stat_info$code
  fake_null_mode <- .normalize_fake_null_mode(fake_null_mode)
  fake_null_code <- .fake_null_mode_code(fake_null_mode)
  complex_null_mode <- match.arg(complex_null_mode)
  complex_null_code <- match(complex_null_mode, c("skip", "cardinality_matched")) - 1L
  n_bins <- .validate_positive_int(n_bins, "n_bins", min = 2L)
  n_fake_per_real <- .validate_positive_int(n_fake_per_real, "n_fake_per_real", min = 1L)
  n_strata <- .validate_positive_int(n_strata, "n_strata", min = 1L)
  caliper_k <- .validate_positive_int(caliper_k, "caliper_k", min = 1L)
  max_match_distance <- as.numeric(max_match_distance)
  if (length(max_match_distance) != 1L || is.na(max_match_distance) ||
      !is.finite(max_match_distance) || max_match_distance <= 0) {
    stop("`max_match_distance` must be a finite positive numeric value", call. = FALSE)
  }
  top_k <- .validate_positive_int(top_k, "top_k", min = 1L)
  min_valid_bins <- .validate_positive_int(min_valid_bins, "min_valid_bins", min = 2L)
  min_n_eff <- .validate_positive_num(min_n_eff, "min_n_eff")
  min_pct <- .validate_unit_interval(min_pct, "min_pct")
  min_valid_fraction <- .validate_unit_interval(min_valid_fraction, "min_valid_fraction")
  min_fake_observed <- .validate_non_negative_int(min_fake_observed, "min_fake_observed")
  min_fake_coverage <- .validate_unit_interval(min_fake_coverage, "min_fake_coverage")
  min_fake_observed_generic <- .validate_non_negative_int(
    min_fake_observed_generic, "min_fake_observed_generic"
  )
  min_fake_coverage_generic <- .validate_unit_interval(
    min_fake_coverage_generic, "min_fake_coverage_generic"
  )
  min_fake_observed_edge <- .validate_non_negative_int(
    min_fake_observed_edge, "min_fake_observed_edge"
  )
  min_fake_coverage_edge <- .validate_unit_interval(
    min_fake_coverage_edge, "min_fake_coverage_edge"
  )
  min_unique_fake_draws <- .validate_non_negative_int(
    min_unique_fake_draws, "min_unique_fake_draws"
  )
  min_unique_fake_fraction <- .validate_unit_interval(
    min_unique_fake_fraction, "min_unique_fake_fraction"
  )
  seed <- .validate_seed(seed)
  expr <- .prepare_lr_expr(expr)
  lr_db <- .validate_lr_db(lr_db)
  .assert_lr_db_indices(lr_db, nrow(expr))
  .assert_indices_in_range(candidate_fake_genes, 1L, nrow(expr),
                           "candidate_fake_genes")
  candidate_input <- .as_integerish(candidate_fake_genes, "candidate_fake_genes")
  candidate_unique <- unique(candidate_input)
  real_lr_genes <- unique(unlist(
    c(lr_db$ligand_rows, lr_db$receptor_rows),
    use.names = FALSE
  ))
  candidate_fake_genes <- setdiff(candidate_unique, real_lr_genes)
  if (length(candidate_fake_genes) < 2L) {
    stop("`candidate_fake_genes` must contain at least two non-LR genes after filtering",
         call. = FALSE)
  }
  n_candidate_input <- length(candidate_input)
  n_candidate_unique <- length(candidate_unique)
  n_duplicates_removed <- n_candidate_input - n_candidate_unique
  n_real_lr_removed <- n_candidate_unique - length(candidate_fake_genes)
  lw <- .assert_lineage_weights(lineage_weights, ncol(expr))
  n_lineages <- ncol(lw)
  ct <- .assert_celltypes(celltypes, ncol(expr))
  ct_levels <- attr(ct, "levels")
  n_types <- max(ct)
  pt_vec <- .assert_pseudotime_1d(pseudotime, ncol(expr))
  blocks <- .prepare_block_labels(block_labels, ncol(expr))
  block_levels <- attr(blocks, "levels")
  raw <- rs_compute_fake_lr_null(
    expr,
    ct,
    as.integer(n_types),
    as.integer(candidate_fake_genes),
    lapply(lr_db$ligand_rows,   as.integer),
    lapply(lr_db$receptor_rows, as.integer),
    pt_vec,
    as.numeric(t(lw)),
    as.integer(n_lineages),
    n_bins,
    min_pct,
    n_fake_per_real,
    n_strata,
    caliper_k,
    max_match_distance,
    as.numeric(seed),
    as.logical(return_bin_scores),
    as.integer(stat_code),
    top_k,
    min_n_eff,
    min_valid_bins,
    min_valid_fraction,
    as.integer(fake_null_code),
    as.integer(complex_null_code),
    blocks
  )
  trajectory_calibration <- data.frame(
    raw$trajectory_calibration,
    stringsAsFactors = FALSE
  )
  trajectory_calibration$source <- ct_levels[trajectory_calibration$source]
  trajectory_calibration$target <- ct_levels[trajectory_calibration$target]
  trajectory_calibration$unique_fake_fraction_generic <- ifelse(
    trajectory_calibration$n_fake_observed_generic > 0L,
    trajectory_calibration$n_unique_fake_draws_generic /
      trajectory_calibration$n_fake_observed_generic,
    NA_real_
  )
  trajectory_calibration$unique_fake_fraction_edge <- ifelse(
    trajectory_calibration$n_fake_observed_edge > 0L,
    trajectory_calibration$n_unique_fake_draws_edge /
      trajectory_calibration$n_fake_observed_edge,
    NA_real_
  )
  trajectory_calibration$unique_fake_fraction_replace_ligand <- ifelse(
    trajectory_calibration$n_fake_observed_replace_ligand > 0L,
    trajectory_calibration$n_unique_fake_draws_replace_ligand /
      trajectory_calibration$n_fake_observed_replace_ligand,
    NA_real_
  )
  trajectory_calibration$unique_fake_fraction_replace_receptor <- ifelse(
    trajectory_calibration$n_fake_observed_replace_receptor > 0L,
    trajectory_calibration$n_unique_fake_draws_replace_receptor /
      trajectory_calibration$n_fake_observed_replace_receptor,
    NA_real_
  )
  trajectory_calibration$unique_fake_fraction_network <- ifelse(
    trajectory_calibration$n_fake_observed_network > 0L,
    trajectory_calibration$n_unique_fake_draws_network /
      trajectory_calibration$n_fake_observed_network,
    NA_real_
  )
  trajectory_calibration$unique_fake_fraction_network_rewire_ligand <- ifelse(
    trajectory_calibration$n_fake_observed_network_rewire_ligand > 0L,
    trajectory_calibration$n_unique_fake_draws_network_rewire_ligand /
      trajectory_calibration$n_fake_observed_network_rewire_ligand,
    NA_real_
  )
  trajectory_calibration$unique_fake_fraction_network_rewire_receptor <- ifelse(
    trajectory_calibration$n_fake_observed_network_rewire_receptor > 0L,
    trajectory_calibration$n_unique_fake_draws_network_rewire_receptor /
      trajectory_calibration$n_fake_observed_network_rewire_receptor,
    NA_real_
  )
  generic_ok <- trajectory_calibration$n_fake_expected_generic > 0L &
    trajectory_calibration$n_fake_observed_generic >= min_fake_observed_generic &
    trajectory_calibration$fake_coverage_generic >= min_fake_coverage_generic &
    trajectory_calibration$n_unique_fake_draws_generic >= min_unique_fake_draws &
    trajectory_calibration$unique_fake_fraction_generic >= min_unique_fake_fraction
  replace_ligand_ok <- trajectory_calibration$n_fake_expected_replace_ligand > 0L &
    trajectory_calibration$n_fake_observed_replace_ligand >= min_fake_observed_edge &
    trajectory_calibration$fake_coverage_replace_ligand >= min_fake_coverage_edge &
    trajectory_calibration$n_unique_fake_draws_replace_ligand >= min_unique_fake_draws &
    trajectory_calibration$unique_fake_fraction_replace_ligand >= min_unique_fake_fraction
  replace_receptor_ok <- trajectory_calibration$n_fake_expected_replace_receptor > 0L &
    trajectory_calibration$n_fake_observed_replace_receptor >= min_fake_observed_edge &
    trajectory_calibration$fake_coverage_replace_receptor >= min_fake_coverage_edge &
    trajectory_calibration$n_unique_fake_draws_replace_receptor >= min_unique_fake_draws &
    trajectory_calibration$unique_fake_fraction_replace_receptor >= min_unique_fake_fraction
  edge_ok <- trajectory_calibration$n_fake_expected_edge > 0L &
    trajectory_calibration$n_fake_observed_edge >= min_fake_observed_edge &
    trajectory_calibration$fake_coverage_edge >= min_fake_coverage_edge &
    trajectory_calibration$n_unique_fake_draws_edge >= min_unique_fake_draws &
    trajectory_calibration$unique_fake_fraction_edge >= min_unique_fake_fraction
  network_ligand_ok <- trajectory_calibration$n_fake_expected_network_rewire_ligand > 0L &
    trajectory_calibration$n_fake_observed_network_rewire_ligand >= min_fake_observed_edge &
    trajectory_calibration$fake_coverage_network_rewire_ligand >= min_fake_coverage_edge &
    trajectory_calibration$n_unique_fake_draws_network_rewire_ligand >= min_unique_fake_draws &
    trajectory_calibration$unique_fake_fraction_network_rewire_ligand >= min_unique_fake_fraction
  network_receptor_ok <- trajectory_calibration$n_fake_expected_network_rewire_receptor > 0L &
    trajectory_calibration$n_fake_observed_network_rewire_receptor >= min_fake_observed_edge &
    trajectory_calibration$fake_coverage_network_rewire_receptor >= min_fake_coverage_edge &
    trajectory_calibration$n_unique_fake_draws_network_rewire_receptor >= min_unique_fake_draws &
    trajectory_calibration$unique_fake_fraction_network_rewire_receptor >= min_unique_fake_fraction
  network_ok <- trajectory_calibration$n_fake_expected_network > 0L &
    trajectory_calibration$n_fake_observed_network >= min_fake_observed_edge &
    trajectory_calibration$fake_coverage_network >= min_fake_coverage_edge &
    trajectory_calibration$n_unique_fake_draws_network >= min_unique_fake_draws &
    trajectory_calibration$unique_fake_fraction_network >= min_unique_fake_fraction
  trajectory_calibration$generic_null_evaluable <- generic_ok
  trajectory_calibration$replace_ligand_null_evaluable <- replace_ligand_ok
  trajectory_calibration$replace_receptor_null_evaluable <- replace_receptor_ok
  trajectory_calibration$edge_null_evaluable <- edge_ok
  trajectory_calibration$network_ligand_null_evaluable <- network_ligand_ok
  trajectory_calibration$network_receptor_null_evaluable <- network_receptor_ok
  trajectory_calibration$network_null_evaluable <- network_ok
  trajectory_calibration$p_fake_generic[!generic_ok] <- NA_real_
  trajectory_calibration$p_replace_ligand[!replace_ligand_ok] <- NA_real_
  trajectory_calibration$p_replace_receptor[!replace_receptor_ok] <- NA_real_
  trajectory_calibration$p_fake_edge[!edge_ok] <- NA_real_
  trajectory_calibration$p_fake_network[!network_ok] <- NA_real_
  trajectory_calibration$p_network_rewire_ligand[!network_ligand_ok] <- NA_real_
  trajectory_calibration$p_network_rewire_receptor[!network_receptor_ok] <- NA_real_
  edge_iu_evaluable <- replace_ligand_ok & replace_receptor_ok
  network_iu_evaluable <- network_ligand_ok & network_receptor_ok
  trajectory_calibration$edge_iu_null_evaluable <- edge_iu_evaluable
  trajectory_calibration$network_iu_null_evaluable <- network_iu_evaluable
  trajectory_calibration$p_edge_iu <- .intersection_union_p(
    trajectory_calibration$p_replace_ligand,
    trajectory_calibration$p_replace_receptor,
    edge_iu_evaluable
  )
  trajectory_calibration$p_network_iu <- .intersection_union_p(
    trajectory_calibration$p_network_rewire_ligand,
    trajectory_calibration$p_network_rewire_receptor,
    network_iu_evaluable
  )
  trajectory_calibration$p_fake <- .primary_p_fake(
    fake_null_mode,
    trajectory_calibration$p_network_iu,
    trajectory_calibration$p_edge_iu,
    trajectory_calibration$p_fake_edge,
    trajectory_calibration$p_fake_generic,
    network_ok,
    edge_ok,
    generic_ok
  )
  trajectory_calibration$identity_inference_mode <- ifelse(
    fake_null_mode %in% c("all", "lr_network_rewire") &
      (is.finite(trajectory_calibration$p_network_iu) | network_ok),
    "lr_network_rewire",
    ifelse(
      is.finite(trajectory_calibration$p_edge_iu) | edge_ok,
      "one_endpoint_replacement",
      ifelse(generic_ok, "context_generic", NA_character_)
    )
  )
  trajectory_calibration$trajectory_stat_name <- stat_info$stat
  trajectory_calibration$fake_null_mode <- fake_null_mode
  trajectory_calibration$calibration_hypothesis <- "lr_identity"
  provenance <- list(
    analysis_data_hash = .compute_analysis_data_hash(
      expr, ct, lr_db, pt_vec, lw, block_labels
    ),
    scoring_config_hash = .compute_scoring_config_hash(
      n_bins, min_pct, min_n_eff, min_valid_bins, min_valid_fraction,
      stat_code, top_k
    ),
    fake_null_config_hash = .compute_fake_null_config_hash(
      candidate_fake_genes, fake_null_mode, n_fake_per_real, seed,
      max_match_distance, block_labels, ncol(expr), complex_null_mode,
      n_strata, caliper_k,
      min_fake_observed_generic, min_fake_coverage_generic,
      min_fake_observed_edge, min_fake_coverage_edge,
      min_unique_fake_draws, min_unique_fake_fraction
    )
  )
  trajectory_calibration <- .attach_provenance_columns(
    trajectory_calibration, provenance
  )
  trajectory_calibration$calibration_config_hash <-
    trajectory_calibration$scoring_config_hash
  trajectory_calibration$hill_scale <- as.numeric(raw$hill_scale)
  complex_pair <- lengths(lr_db$ligand_rows) != 1L |
    lengths(lr_db$receptor_rows) != 1L
  trajectory_calibration$complex_pair <- complex_pair[
    trajectory_calibration$pair_idx
  ]
  trajectory_calibration$complex_null_mode <- complex_null_mode
  if (complex_null_mode == "skip") {
    trajectory_calibration$p_fake[
      trajectory_calibration$complex_pair
    ] <- NA_real_
    trajectory_calibration$p_fake_generic[
      trajectory_calibration$complex_pair
    ] <- NA_real_
    trajectory_calibration$p_fake_complex[
      trajectory_calibration$complex_pair
    ] <- NA_real_
    trajectory_calibration$p_fake_edge[
      trajectory_calibration$complex_pair
    ] <- NA_real_
    trajectory_calibration$p_replace_ligand[
      trajectory_calibration$complex_pair
    ] <- NA_real_
    trajectory_calibration$p_replace_receptor[
      trajectory_calibration$complex_pair
    ] <- NA_real_
    trajectory_calibration$p_edge_iu[
      trajectory_calibration$complex_pair
    ] <- NA_real_
    trajectory_calibration$p_fake_network[
      trajectory_calibration$complex_pair
    ] <- NA_real_
    trajectory_calibration$p_network_rewire_ligand[
      trajectory_calibration$complex_pair
    ] <- NA_real_
    trajectory_calibration$p_network_rewire_receptor[
      trajectory_calibration$complex_pair
    ] <- NA_real_
    trajectory_calibration$p_network_iu[
      trajectory_calibration$complex_pair
    ] <- NA_real_
    trajectory_calibration$identity_calibration_evaluable <-
      !trajectory_calibration$complex_pair &
      (is.finite(trajectory_calibration$p_network_iu) |
       is.finite(trajectory_calibration$p_fake_network) |
       is.finite(trajectory_calibration$p_network_rewire_ligand) |
       is.finite(trajectory_calibration$p_network_rewire_receptor) |
       is.finite(trajectory_calibration$p_edge_iu) |
       is.finite(trajectory_calibration$p_fake_edge) |
       is.finite(trajectory_calibration$p_replace_ligand) |
       is.finite(trajectory_calibration$p_replace_receptor) |
       is.finite(trajectory_calibration$p_fake_generic) |
       is.finite(trajectory_calibration$p_fake))
  } else {
    if (any(complex_pair)) {
      unsupported <- vapply(lr_db$ligand_rows, length, integer(1L)) > 8L |
        vapply(lr_db$receptor_rows, length, integer(1L)) > 8L
      if (any(unsupported)) {
        warning(
          "Cardinality-matched FakeLR supports at most 8 subunits per side; ",
          "pairs exceeding this limit receive NA p_fake_complex.",
          call. = FALSE
        )
      }
    }
    complex_idx <- trajectory_calibration$complex_pair
    trajectory_calibration$p_fake_complex[!generic_ok] <- NA_real_
    trajectory_calibration$p_fake_edge[!edge_ok] <- NA_real_
    trajectory_calibration$p_replace_ligand[!replace_ligand_ok] <- NA_real_
    trajectory_calibration$p_replace_receptor[!replace_receptor_ok] <- NA_real_
    trajectory_calibration$p_fake_network[!network_ok] <- NA_real_
    trajectory_calibration$p_network_rewire_ligand[!network_ligand_ok] <- NA_real_
    trajectory_calibration$p_network_rewire_receptor[!network_receptor_ok] <- NA_real_
    trajectory_calibration$p_edge_iu <- .intersection_union_p(
      trajectory_calibration$p_replace_ligand,
      trajectory_calibration$p_replace_receptor,
      replace_ligand_ok & replace_receptor_ok
    )
    trajectory_calibration$p_network_iu <- .intersection_union_p(
      trajectory_calibration$p_network_rewire_ligand,
      trajectory_calibration$p_network_rewire_receptor,
      network_ligand_ok & network_receptor_ok
    )
    complex_generic_p <- ifelse(
      is.finite(trajectory_calibration$p_fake_complex),
      trajectory_calibration$p_fake_complex,
      trajectory_calibration$p_fake_generic
    )
    trajectory_calibration$p_fake[complex_idx] <- .primary_p_fake(
      fake_null_mode,
      trajectory_calibration$p_network_iu[complex_idx],
      trajectory_calibration$p_edge_iu[complex_idx],
      trajectory_calibration$p_fake_edge[complex_idx],
      complex_generic_p[complex_idx],
      network_ok[complex_idx],
      edge_ok[complex_idx],
      generic_ok[complex_idx]
    )
    trajectory_calibration$identity_inference_mode[complex_idx] <- ifelse(
      fake_null_mode %in% c("all", "lr_network_rewire") &
        (is.finite(trajectory_calibration$p_network_iu[complex_idx]) |
         network_ok[complex_idx]),
      "lr_network_rewire_complex",
      ifelse(
        is.finite(trajectory_calibration$p_edge_iu[complex_idx]) |
          edge_ok[complex_idx],
        "one_endpoint_replacement_complex",
        ifelse(
          generic_ok[complex_idx], "context_generic_complex", NA_character_
        )
      )
    )
    simple_evaluable <- !complex_idx &
      (is.finite(trajectory_calibration$p_network_iu) |
       is.finite(trajectory_calibration$p_fake_network) |
       is.finite(trajectory_calibration$p_network_rewire_ligand) |
       is.finite(trajectory_calibration$p_network_rewire_receptor) |
       is.finite(trajectory_calibration$p_edge_iu) |
       is.finite(trajectory_calibration$p_fake_edge) |
       is.finite(trajectory_calibration$p_replace_ligand) |
       is.finite(trajectory_calibration$p_replace_receptor) |
       is.finite(trajectory_calibration$p_fake_generic) |
       is.finite(trajectory_calibration$p_fake))
    complex_evaluable <- complex_idx &
      (
        (generic_ok & is.finite(trajectory_calibration$p_fake_complex)) |
          (edge_ok & is.finite(trajectory_calibration$p_edge_iu)) |
          (network_ok & is.finite(trajectory_calibration$p_network_iu)) |
          is.finite(trajectory_calibration$p_fake)
      )
    trajectory_calibration$identity_calibration_evaluable <-
      simple_evaluable | complex_evaluable
  }
  exploratory <- isTRUE(n_fake_per_real < 100L)
  if (exploratory) {
    stop(
      "ComputeFakeLRNull requires `n_fake_per_real >= 100` for formal BH q-values; ",
      "got ", n_fake_per_real, ". Use ComputeTrajectoryFakeLRNull(n_fake_per_real >= 1000) ",
      "for trajectory-level inference.",
      call. = FALSE
    )
  } else {
    trajectory_calibration$q_fake <- stats::p.adjust(
      trajectory_calibration$p_fake,
      method = "BH"
    )
    trajectory_calibration$q_fake_edge <- stats::p.adjust(
      trajectory_calibration$p_fake_edge,
      method = "BH"
    )
    trajectory_calibration$q_fake_generic <- stats::p.adjust(
      trajectory_calibration$p_fake_generic,
      method = "BH"
    )
    trajectory_calibration$q_replace_ligand <- stats::p.adjust(
      trajectory_calibration$p_replace_ligand,
      method = "BH"
    )
    trajectory_calibration$q_replace_receptor <- stats::p.adjust(
      trajectory_calibration$p_replace_receptor,
      method = "BH"
    )
    trajectory_calibration$q_edge_iu <- stats::p.adjust(
      trajectory_calibration$p_edge_iu,
      method = "BH"
    )
    trajectory_calibration$q_fake_network <- stats::p.adjust(
      trajectory_calibration$p_fake_network,
      method = "BH"
    )
    trajectory_calibration$q_network_rewire_ligand <- stats::p.adjust(
      trajectory_calibration$p_network_rewire_ligand,
      method = "BH"
    )
    trajectory_calibration$q_network_rewire_receptor <- stats::p.adjust(
      trajectory_calibration$p_network_rewire_receptor,
      method = "BH"
    )
    trajectory_calibration$q_network_iu <- stats::p.adjust(
      trajectory_calibration$p_network_iu,
      method = "BH"
    )
    trajectory_calibration$q_fake_complex <- stats::p.adjust(
      trajectory_calibration$p_fake_complex,
      method = "BH"
    )
  }
  attr(trajectory_calibration, "exploratory") <- FALSE
  fake_scores <- data.frame(raw$fake_scores, stringsAsFactors = FALSE)
  fake_scores$source <- ct_levels[fake_scores$source]
  fake_scores$target <- ct_levels[fake_scores$target]
  attr(fake_scores, "celltype_levels") <- ct_levels
  attr(trajectory_calibration, "celltype_levels") <- ct_levels
  fake_pairs <- data.frame(raw$fake_pairs, stringsAsFactors = FALSE)
  fake_pairs$source <- ct_levels[fake_pairs$source]
  fake_pairs$target <- ct_levels[fake_pairs$target]
  matching_balance <- .compute_fake_matching_balance(fake_pairs)
  if (nrow(matching_balance) > 0L && nrow(trajectory_calibration) > 0L) {
    trajectory_calibration <- merge(
      trajectory_calibration,
      matching_balance,
      by = c("pair_idx", "source", "target", "lineage"),
      all.x = TRUE,
      sort = FALSE
    )
  }
  matching_mode <- "context_aware_caliper"
  n_traj_keys <- as.integer(raw$n_traj_keys)
  n_expected_total <- n_traj_keys * n_fake_per_real
  matched_fraction <- as.integer(raw$n_fake_total) / n_expected_total
  unique_pair_fraction <- ifelse(
    nrow(fake_pairs) > 0L,
    nrow(unique(fake_pairs[c("ligand_gene", "receptor_gene", "shadowing_pair")])) /
      nrow(fake_pairs),
    NA_real_
  )
  matching_audit <- data.frame(
    matching_mode = matching_mode,
    fake_null_mode = fake_null_mode,
    n_traj_keys = n_traj_keys,
    max_match_distance = max_match_distance,
    n_matching_blocks = as.integer(raw$n_matching_blocks),
    matching_block_mode = if (as.integer(raw$n_matching_blocks) > 1L) {
      "block_x_context"
    } else if (is.null(block_labels)) {
      "pooled"
    } else {
      "single_block"
    },
    n_expected_total = n_expected_total,
    n_matched_total = as.integer(raw$n_fake_total),
    matched_fraction = matched_fraction,
    unmatched_fraction = 1 - matched_fraction,
    unique_pair_fraction = unique_pair_fraction,
    mean_ligand_match_distance = ifelse(
      nrow(fake_pairs) > 0L,
      mean(fake_pairs$ligand_match_distance),
      NA_real_
    ),
    mean_receptor_match_distance = ifelse(
      nrow(fake_pairs) > 0L,
      mean(fake_pairs$receptor_match_distance),
      NA_real_
    ),
    stringsAsFactors = FALSE
  )
  list(
    fake_scores         = fake_scores,
    fake_pairs          = fake_pairs,
    trajectory_calibration = trajectory_calibration,
    matching_audit      = matching_audit,
    matching_balance    = matching_balance,
    n_fake_total        = as.integer(raw$n_fake_total),
    candidate_pool_size = as.integer(raw$candidate_pool_size),
    n_candidate_input   = as.integer(n_candidate_input),
    n_candidate_unique  = as.integer(n_candidate_unique),
    n_duplicates_removed = as.integer(n_duplicates_removed),
    n_real_lr_removed   = as.integer(n_real_lr_removed),
    celltype_levels     = ct_levels,
    n_real_pairs        = as.integer(raw$n_real_pairs),
    analysis_data_hash  = provenance$analysis_data_hash,
    scoring_config_hash = provenance$scoring_config_hash,
    calibration_config_hash = provenance$scoring_config_hash,
    fake_null_config_hash = provenance$fake_null_config_hash,
    hill_scale          = as.numeric(raw$hill_scale)
  )
}

#' Low-memory trajectory-level FakeLR null.
#'
#' Recommended publication path.  Generates many matched fake pairs in Rust,
#' samples with replacement inside each feature-matched caliper, streams their
#' trajectory statistics on the same Hill scale as the real scores, and returns
#' calibrated real trajectories without materialising the full fake bin table.
#'
#' @inheritParams ComputeFakeLRNull
#' @param n_bins number of pseudotime bins.
#' @param n_fake_per_real number of matched fake trajectories per real LR pair.
#' @param n_strata number of matching strata per feature axis.
#' @param caliper_k nearest matched candidate genes retained per side.
#' @param min_fake_observed minimum observed fake trajectories required.
#' @param min_fake_coverage minimum observed/expected fake coverage.
#' @return list with `trajectories`, `fake_pairs`, and audit counts.
#' @export
ComputeTrajectoryFakeLRNull <- function(expr, celltypes, lr_db,
                                        candidate_fake_genes,
                                        pseudotime, lineage_weights,
                                        n_bins = 12L,
                                        n_fake_per_real = 1000L,
                                        n_strata = 5L, caliper_k = 50L,
                                        max_match_distance = 0.5,
                                        fake_null_mode = c("context_generic", "edge_rewire", "both",
                                                 "lr_network_rewire", "all", "triple"),
                                        min_pct = 0.05, seed = 12345L,
                                        min_n_eff = 5.0,
                                        min_valid_bins = 3L,
                                        min_valid_fraction = 0.25,
                                        trajectory_stat = c("amplitude", "max", "mean_score", "auc", "topk_mean"),
                                        top_k = 3L,
                                        min_fake_observed = 100L,
                                        min_fake_coverage = 0.8,
                                        block_labels = NULL,
                                        min_fake_observed_generic = min_fake_observed,
                                        min_fake_coverage_generic = min_fake_coverage,
                                        min_fake_observed_edge = min_fake_observed,
                                        min_fake_coverage_edge = min_fake_coverage,
                                        min_unique_fake_draws = 10L,
                                        min_unique_fake_fraction = 0.01,
                                        complex_null_mode = c("skip", "cardinality_matched")) {
  trajectory_stat <- match.arg(trajectory_stat)
  raw <- ComputeFakeLRNull(
    expr = expr,
    celltypes = celltypes,
    lr_db = lr_db,
    candidate_fake_genes = candidate_fake_genes,
    pseudotime = pseudotime,
    lineage_weights = lineage_weights,
    n_bins = n_bins,
    n_fake_per_real = n_fake_per_real,
    n_strata = n_strata,
    caliper_k = caliper_k,
    max_match_distance = max_match_distance,
    fake_null_mode = fake_null_mode,
    min_pct = min_pct,
    seed = seed,
    return_bin_scores = FALSE,
    trajectory_stat = trajectory_stat,
    top_k = top_k,
    min_n_eff = min_n_eff,
    min_valid_bins = min_valid_bins,
    min_valid_fraction = min_valid_fraction,
    min_fake_observed = min_fake_observed,
    min_fake_coverage = min_fake_coverage,
    block_labels = block_labels,
    min_fake_observed_generic = min_fake_observed_generic,
    min_fake_coverage_generic = min_fake_coverage_generic,
    min_fake_observed_edge = min_fake_observed_edge,
    min_fake_coverage_edge = min_fake_coverage_edge,
    min_unique_fake_draws = min_unique_fake_draws,
    min_unique_fake_fraction = min_unique_fake_fraction,
    complex_null_mode = complex_null_mode
  )
  trajectories <- raw$trajectory_calibration
  if (isTRUE(n_fake_per_real < 100L)) {
    stop(
      "ComputeTrajectoryFakeLRNull requires `n_fake_per_real >= 100` for formal ",
      "BH q-values; got ", n_fake_per_real, ". Recommended: >= 1000.",
      call. = FALSE
    )
  }
  trajectories$q_fake <- stats::p.adjust(trajectories$p_fake, method = "BH")
  if ("p_fake_edge" %in% colnames(trajectories)) {
    trajectories$q_fake_edge <- stats::p.adjust(trajectories$p_fake_edge, method = "BH")
  }
  if ("p_fake_generic" %in% colnames(trajectories)) {
    trajectories$q_fake_generic <- stats::p.adjust(
      trajectories$p_fake_generic,
      method = "BH"
    )
  }
  if ("p_replace_ligand" %in% colnames(trajectories)) {
    trajectories$q_replace_ligand <- stats::p.adjust(
      trajectories$p_replace_ligand,
      method = "BH"
    )
  }
  if ("p_replace_receptor" %in% colnames(trajectories)) {
    trajectories$q_replace_receptor <- stats::p.adjust(
      trajectories$p_replace_receptor,
      method = "BH"
    )
  }
  if ("p_edge_iu" %in% colnames(trajectories)) {
    trajectories$q_edge_iu <- stats::p.adjust(trajectories$p_edge_iu, method = "BH")
  }
  if ("p_network_iu" %in% colnames(trajectories)) {
    trajectories$q_network_iu <- stats::p.adjust(trajectories$p_network_iu, method = "BH")
  }
  if ("p_fake_network" %in% colnames(trajectories)) {
    trajectories$q_fake_network <- stats::p.adjust(trajectories$p_fake_network, method = "BH")
  }
  if ("p_network_rewire_ligand" %in% colnames(trajectories)) {
    trajectories$q_network_rewire_ligand <- stats::p.adjust(
      trajectories$p_network_rewire_ligand,
      method = "BH"
    )
  }
  if ("p_network_rewire_receptor" %in% colnames(trajectories)) {
    trajectories$q_network_rewire_receptor <- stats::p.adjust(
      trajectories$p_network_rewire_receptor,
      method = "BH"
    )
  }
  list(
    trajectories = trajectories,
    fake_pairs = raw$fake_pairs,
    matching_audit = raw$matching_audit,
    matching_balance = raw$matching_balance,
    n_fake_total = raw$n_fake_total,
    candidate_pool_size = raw$candidate_pool_size,
    n_candidate_input = raw$n_candidate_input,
    n_candidate_unique = raw$n_candidate_unique,
    n_duplicates_removed = raw$n_duplicates_removed,
    n_real_lr_removed = raw$n_real_lr_removed,
    celltype_levels = raw$celltype_levels,
    n_real_pairs = raw$n_real_pairs,
    analysis_data_hash = raw$analysis_data_hash,
    scoring_config_hash = raw$scoring_config_hash,
    calibration_config_hash = raw$scoring_config_hash,
    fake_null_config_hash = raw$fake_null_config_hash,
    hill_scale = raw$hill_scale
  )
}

.merge_adaptive_fake_null <- function(screen_traj, confirm_traj,
                                      candidate_keys,
                                      n_fake_screen, n_fake_confirm,
                                      p_col = "p_fake",
                                      by = c("pair_idx", "source", "target", "lineage")) {
  traj <- screen_traj
  traj[[paste0(p_col, "_screen")]] <- traj[[p_col]]
  q_col <- sub("^p_", "q_", p_col)
  if (q_col %in% colnames(traj)) {
    traj[[paste0(q_col, "_screen")]] <- traj[[q_col]]
  }
  screen_observed <- if ("n_fake_observed" %in% colnames(traj)) {
    traj$n_fake_observed
  } else {
    n_fake_screen
  }
  traj <- .adaptive_perm_attach_diagnostics(
    traj, p_col, screen_observed, "screen"
  )
  keys <- .trajectory_key_col(traj, by)
  is_candidate <- keys %in% candidate_keys
  traj$adaptive_perm_status <- ifelse(is_candidate, "pending_confirm", "screened_out")

  if (length(candidate_keys) > 0L && !is.null(confirm_traj)) {
    conf_keys <- .trajectory_key_col(confirm_traj, by)
    hit <- match(keys[is_candidate], conf_keys)
    idx <- which(is_candidate)
    ok <- !is.na(hit)
    if (any(ok)) {
      replace_cols <- setdiff(
        intersect(colnames(confirm_traj), colnames(traj)),
        by
      )
      for (col in replace_cols) {
        traj[[col]][idx[ok]] <- confirm_traj[[col]][hit[ok]]
      }
      traj$perm_stage[idx[ok]] <- "confirm"
      confirm_observed <- if ("n_fake_observed" %in% colnames(confirm_traj)) {
        confirm_traj$n_fake_observed[hit[ok]]
      } else {
        rep(n_fake_confirm, sum(ok))
      }
      traj$n_perm_used[idx[ok]] <- as.integer(confirm_observed)
      traj$adaptive_perm_status[idx[ok]] <- "confirmed"
    }
  }

  traj$adaptive_perm_status[traj$adaptive_perm_status == "pending_confirm"] <- "screened_out"
  if ("n_fake_observed" %in% colnames(traj)) {
    use_observed <- is.finite(traj$n_fake_observed) & traj$n_fake_observed >= 0L
    traj$n_perm_used[use_observed] <- as.integer(traj$n_fake_observed[use_observed])
  }
  diag_all <- .adaptive_perm_diagnostics(traj[[p_col]], traj$n_perm_used)
  traj$min_achievable_p <- diag_all$min_achievable_p
  traj$monte_carlo_se <- diag_all$monte_carlo_se
  .recompute_all_bh_q(traj)
}

#' Adaptive two-stage trajectory-level FakeLR null.
#'
#' Two-stage FakeLR with a locked null universe: both stages use the complete LR
#' database and fake-candidate universe; confirm-stage results replace every
#' null sub-layer, quality audit, median, and stage hash for candidate rows.
#' Final `q_fake` uses BH on the full trajectory family. Returned matching
#' tables retain both stages through `adaptive_stage`. Confirmation currently
#' recomputes the full universe, so this allocates precision adaptively but is
#' not a computational acceleration.
#'
#' @inheritParams ComputeTrajectoryFakeLRNull
#' @param n_fake_screen matched fakes per real pair for screening (default 199).
#' @param n_fake_confirm matched fakes per real pair for confirmation (default 4999).
#' @param q_confirm_threshold screen-stage `p_fake` cutoff (default 0.1).
#' @param confirm_top_k additionally confirm the smallest `confirm_top_k`
#'   finite screen p-values (default 20; set 0 to disable).
#' @return list as `ComputeTrajectoryFakeLRNull()` plus `adaptive_perm_audit`.
#' @export
ComputeAdaptiveTrajectoryFakeLRNull <- function(expr, celltypes, lr_db,
                                                candidate_fake_genes,
                                                pseudotime, lineage_weights,
                                                n_bins = 12L,
                                                n_fake_screen = 199L,
                                                n_fake_confirm = 4999L,
                                                q_confirm_threshold = 0.1,
                                                confirm_top_k = 20L,
                                                n_strata = 5L, caliper_k = 50L,
                                                max_match_distance = 0.5,
                                                fake_null_mode = c("context_generic", "edge_rewire", "both",
                                                 "lr_network_rewire", "all", "triple"),
                                                min_pct = 0.05, seed = 12345L,
                                                min_n_eff = 5.0,
                                                min_valid_bins = 3L,
                                                min_valid_fraction = 0.25,
                                                trajectory_stat = c("amplitude", "max", "mean_score", "auc", "topk_mean"),
                                                top_k = 3L,
                                                min_fake_observed = 100L,
                                                min_fake_coverage = 0.8,
                                                block_labels = NULL,
                                                min_fake_observed_generic = min_fake_observed,
                                                min_fake_coverage_generic = min_fake_coverage,
                                                min_fake_observed_edge = min_fake_observed,
                                                min_fake_coverage_edge = min_fake_coverage,
                                                min_unique_fake_draws = 10L,
                                                min_unique_fake_fraction = 0.01) {
  n_fake_screen <- .validate_positive_int(n_fake_screen, "n_fake_screen", min = 1L)
  n_fake_confirm <- .validate_positive_int(n_fake_confirm, "n_fake_confirm", min = 1L)
  q_confirm_threshold <- .validate_unit_interval(q_confirm_threshold, "q_confirm_threshold")
  confirm_top_k <- .validate_non_negative_int(confirm_top_k, "confirm_top_k")
  seed <- .validate_seed(seed)
  fake_args <- list(
    expr = expr,
    celltypes = celltypes,
    lr_db = lr_db,
    candidate_fake_genes = candidate_fake_genes,
    pseudotime = pseudotime,
    lineage_weights = lineage_weights,
    n_bins = n_bins,
    n_strata = n_strata,
    caliper_k = caliper_k,
    max_match_distance = max_match_distance,
    fake_null_mode = fake_null_mode,
    min_pct = min_pct,
    min_n_eff = min_n_eff,
    min_valid_bins = min_valid_bins,
    min_valid_fraction = min_valid_fraction,
    trajectory_stat = trajectory_stat,
    top_k = top_k,
    min_fake_observed = min_fake_observed,
    min_fake_coverage = min_fake_coverage,
    block_labels = block_labels,
    min_fake_observed_generic = min_fake_observed_generic,
    min_fake_coverage_generic = min_fake_coverage_generic,
    min_fake_observed_edge = min_fake_observed_edge,
    min_fake_coverage_edge = min_fake_coverage_edge,
    min_unique_fake_draws = min_unique_fake_draws,
    min_unique_fake_fraction = min_unique_fake_fraction
  )
  screen_out <- do.call(
    ComputeTrajectoryFakeLRNull,
    c(fake_args, list(n_fake_per_real = n_fake_screen, seed = seed))
  )
  by <- c("pair_idx", "source", "target", "lineage")
  candidate_keys <- .adaptive_perm_select_candidates(
    screen_out$trajectories,
    p_col = "p_fake",
    threshold = q_confirm_threshold,
    top_k = confirm_top_k,
    by = by
  )
  confirm_traj <- NULL
  confirm_out <- NULL
  if (length(candidate_keys) > 0L) {
    confirm_out <- do.call(
      ComputeTrajectoryFakeLRNull,
      utils::modifyList(fake_args, list(
        n_fake_per_real = n_fake_confirm,
        seed = seed + 1
      ))
    )
    confirm_traj <- confirm_out$trajectories
  }
  trajectories <- .merge_adaptive_fake_null(
    screen_out$trajectories,
    confirm_traj,
    candidate_keys,
    n_fake_screen,
    n_fake_confirm,
    p_col = "p_fake",
    by = by
  )
  adaptive_hash <- .compute_adaptive_perm_config_hash(
    n_fake_screen, n_fake_confirm, q_confirm_threshold, confirm_top_k, seed,
    kind = "fake_lr"
  )
  trajectories$adaptive_perm_config_hash <- adaptive_hash
  n_confirmed <- sum(trajectories$adaptive_perm_status == "confirmed")
  fake_pairs <- screen_out$fake_pairs
  fake_pairs$adaptive_stage <- rep("screen", nrow(fake_pairs))
  matching_audit <- screen_out$matching_audit
  matching_audit$adaptive_stage <- rep("screen", nrow(matching_audit))
  matching_balance <- screen_out$matching_balance
  matching_balance$adaptive_stage <- rep("screen", nrow(matching_balance))
  if (!is.null(confirm_out)) {
    confirm_pairs <- confirm_out$fake_pairs
    confirm_pairs$adaptive_stage <- rep("confirm", nrow(confirm_pairs))
    fake_pairs <- rbind(fake_pairs, confirm_pairs)
    confirm_audit <- confirm_out$matching_audit
    confirm_audit$adaptive_stage <- rep("confirm", nrow(confirm_audit))
    matching_audit <- rbind(matching_audit, confirm_audit)
    confirm_balance <- confirm_out$matching_balance
    confirm_balance$adaptive_stage <- rep("confirm", nrow(confirm_balance))
    matching_balance <- rbind(matching_balance, confirm_balance)
  }
  list(
    trajectories = trajectories,
    fake_pairs = fake_pairs,
    matching_audit = matching_audit,
    matching_balance = matching_balance,
    n_fake_total = screen_out$n_fake_total + if (is.null(confirm_out)) {
      0L
    } else {
      confirm_out$n_fake_total
    },
    candidate_pool_size = screen_out$candidate_pool_size,
    n_candidate_input = screen_out$n_candidate_input,
    n_candidate_unique = screen_out$n_candidate_unique,
    n_duplicates_removed = screen_out$n_duplicates_removed,
    n_real_lr_removed = screen_out$n_real_lr_removed,
    celltype_levels = screen_out$celltype_levels,
    n_real_pairs = screen_out$n_real_pairs,
    analysis_data_hash = screen_out$analysis_data_hash,
    scoring_config_hash = screen_out$scoring_config_hash,
    calibration_config_hash = screen_out$scoring_config_hash,
    fake_null_config_hash = adaptive_hash,
    stage_fake_null_config_hashes = unique(trajectories$fake_null_config_hash),
    hill_scale = screen_out$hill_scale,
    adaptive_perm_config_hash = adaptive_hash,
    adaptive_perm_audit = list(
      n_trajectories = nrow(trajectories),
      n_candidates = length(candidate_keys),
      n_confirmed = n_confirmed,
      n_screened_out = sum(trajectories$adaptive_perm_status == "screened_out"),
      n_fake_screen = n_fake_screen,
      n_fake_confirm = n_fake_confirm,
      q_confirm_threshold = q_confirm_threshold,
      confirm_top_k = confirm_top_k
    )
  )
}

#' Compute per-row p_fake values from real and fake LR scores.
#'
#' For each real LR row, counts how many fake LR rows at the same
#' `(pair_idx, source, target, lineage, bin)` position have score >= the
#' real score, then applies Laplace smoothing:
#' `p_fake = (count_geq + 1) / (n_fake_observed + 1)`.
#'
#' @param real_scores data.frame from `ComputeTemporalLR` /
#'   `ComputeTemporalLRWithIntermediate` (its `scores` element).
#' @param fake_scores data.frame, the `fake_scores` element from
#'   `ComputeFakeLRNull`.
#' @param n_fake_total deprecated; retained for source compatibility.  The
#'   denominator is the number of actually matched fake rows at the row's
#'   comparison key, not the global fake-pair count.
#' @return numeric vector of length `nrow(real_scores)`, each in (0, 1].
#' @export
ComputeFakePValues <- function(real_scores, fake_scores, n_fake_total = NULL) {
  detail <- ComputeFakePValuesDetailed(real_scores, fake_scores)
  p <- detail$p_fake
  attr(p, "n_fake_observed") <- detail$n_fake_observed
  attr(p, "n_fake_expected") <- detail$n_fake_expected
  attr(p, "fake_coverage") <- detail$fake_coverage
  p
}

#' Detailed per-row p_fake values and matched-fake coverage.
#'
#' @inheritParams ComputeFakePValues
#' @param n_fake_expected optional scalar or vector giving the expected number
#'   of matched fake rows per real row.  Used only for `fake_coverage`.
#' @param min_fake_observed minimum matched fake rows required to report p_fake.
#' @param min_fake_coverage minimum observed/expected coverage required when
#'   `n_fake_expected` is supplied.
#' @return data.frame with `p_fake`, `n_fake_expected`, `n_fake_observed`,
#'   and `fake_coverage`.
#' @export
ComputeFakePValuesDetailed <- function(real_scores, fake_scores,
                                       n_fake_expected = NA_integer_,
                                       min_fake_observed = 0L,
                                       min_fake_coverage = 0) {
  stopifnot(is.data.frame(real_scores), is.data.frame(fake_scores),
            all(c("pair_idx", "source", "target", "lineage", "bin", "score")
                %in% colnames(real_scores)),
            all(c("pair_idx", "source", "target", "lineage", "bin", "score")
                %in% colnames(fake_scores)))

  # Build only observed keys; avoid a Cartesian product of factor levels.
  fake_key <- do.call(
    paste,
    c(fake_scores[c("pair_idx", "source", "target", "lineage", "bin")],
      sep = ".")
  )
  fake_list <- split(fake_scores$score, fake_key)

  expected <- rep_len(as.numeric(n_fake_expected), nrow(real_scores))
  observed <- vapply(seq_len(nrow(real_scores)), function(i) {
    key <- paste(real_scores$pair_idx[i], real_scores$source[i],
                 real_scores$target[i], real_scores$lineage[i],
                 real_scores$bin[i], sep = ".")
    f_scores <- fake_list[[key]]
    if (is.null(f_scores)) 0L else length(f_scores)
  }, integer(1))
  p_fake <- vapply(seq_len(nrow(real_scores)), function(i) {
    key <- paste(real_scores$pair_idx[i], real_scores$source[i],
                 real_scores$target[i], real_scores$lineage[i],
                 real_scores$bin[i], sep = ".")
    f_scores <- fake_list[[key]]
    if (is.null(f_scores) || !length(f_scores)) {
      return(NA_real_)
    }
    count_geq <- sum(f_scores >= real_scores$score[i])
    (count_geq + 1.0) / (length(f_scores) + 1.0)
  }, numeric(1))
  coverage <- ifelse(is.finite(expected) & expected > 0,
                     observed / expected,
                     NA_real_)
  p_fake[observed < as.integer(min_fake_observed)] <- NA_real_
  p_fake[is.finite(coverage) & coverage < min_fake_coverage] <- NA_real_
  data.frame(
    p_fake = p_fake,
    n_fake_expected = expected,
    n_fake_observed = observed,
    fake_coverage = coverage,
    stringsAsFactors = FALSE
  )
}

#' Trajectory-level FakeLR support (deprecated).
#'
#' **Hard error — do not use for formal inference.** Missing bins are
#' zero-filled, which inflates amplitude statistics. Use
#' [ComputeTrajectoryFakeLRNull()] for publication analyses.
#'
#' @param real_scores data.frame from `ComputeTemporalLR` /
#'   `ComputeTemporalLRWithIntermediate`.
#' @param fake_scores data.frame from `ComputeFakeLRNull$fake_scores`.
#' @param n_bins optional integer number of bins.
#' @param stat trajectory statistic.
#' @param top_k number of bins used by `"topk_mean"`.
#' @param n_fake_expected optional scalar or vector for coverage reporting.
#' @param min_fake_observed,min_fake_coverage ignored; retained for API stability.
#' @return Never returns; always errors.
#' @export
ComputeTrajectoryFakePValues <- function(real_scores, fake_scores,
                                         n_bins = NULL,
                                         stat = c("amplitude", "max", "auc", "topk_mean"),
                                         top_k = 3L,
                                         n_fake_expected = NA_integer_,
                                         min_fake_observed = 100L,
                                         min_fake_coverage = 0) {
  .Deprecated(
    msg = paste0(
      "ComputeTrajectoryFakePValues() is deprecated: it zero-fills missing bins ",
      "and produces biased trajectory statistics. ",
      "Use ComputeTrajectoryFakeLRNull() for publication inference."
    )
  )
  stop(
    "ComputeTrajectoryFakePValues() must not be used for formal inference ",
    "(missing bins are zero-filled, inflating amplitude). ",
    "Use ComputeTrajectoryFakeLRNull() instead.",
    call. = FALSE
  )
}

#' @keywords internal
.ComputeTrajectoryFakePValues_legacy <- function(real_scores, fake_scores,
                                         n_bins = NULL,
                                         stat = c("amplitude", "max", "auc", "topk_mean"),
                                         top_k = 3L,
                                         n_fake_expected = NA_integer_,
                                         min_fake_observed = 100L,
                                         min_fake_coverage = 0) {
  stat <- match.arg(stat)
  stopifnot(is.data.frame(real_scores), is.data.frame(fake_scores),
            all(c("pair_idx", "source", "target", "lineage", "bin", "score")
                %in% colnames(real_scores)),
            all(c("fake_id", "pair_idx", "source", "target", "lineage", "bin", "score")
                %in% colnames(fake_scores)))
  n_bins <- if (is.null(n_bins)) {
    max(c(real_scores$bin, fake_scores$bin), na.rm = TRUE)
  } else {
    as.integer(n_bins)
  }

  key_cols <- c("pair_idx", "source", "target", "lineage")
  real_key <- do.call(paste, c(real_scores[key_cols], sep = "."))
  fake_key <- do.call(paste, c(fake_scores[c("fake_id", key_cols)], sep = "."))
  real_split <- split(real_scores, real_key)
  fake_split <- split(fake_scores, fake_key)

  stat_one <- function(df) {
    y <- numeric(n_bins)
    b <- as.integer(df$bin)
    y[b] <- pmax(y[b], as.numeric(df$score))
    if (stat == "amplitude") {
      return(max(y) - min(y))
    }
    if (stat == "max") {
      return(max(y))
    }
    if (stat == "auc") {
      return(mean(y))
    }
    kk <- min(as.integer(top_k), length(y))
    mean(sort(y, decreasing = TRUE)[seq_len(kk)])
  }

  real_stat <- vapply(real_split, stat_one, numeric(1))
  real_keys <- names(real_split)
  key_df <- do.call(rbind, strsplit(real_keys, ".", fixed = TRUE))
  key_df <- data.frame(
    pair_idx = as.integer(key_df[, 1]),
    source = as.integer(key_df[, 2]),
    target = as.integer(key_df[, 3]),
    lineage = as.integer(key_df[, 4]),
    stringsAsFactors = FALSE
  )
  fake_stat <- vapply(fake_split, stat_one, numeric(1))
  fake_key_df <- do.call(rbind, strsplit(names(fake_split), ".", fixed = TRUE))
  fake_group_key <- paste(fake_key_df[, 2], fake_key_df[, 3],
                          fake_key_df[, 4], fake_key_df[, 5], sep = ".")
  fake_stat_list <- split(fake_stat, fake_group_key)

  expected <- rep_len(as.numeric(n_fake_expected), length(real_keys))
  observed <- vapply(real_keys, function(key) {
    x <- fake_stat_list[[key]]
    if (is.null(x)) 0L else length(x)
  }, integer(1))
  p_fake <- vapply(seq_along(real_keys), function(i) {
    fs <- fake_stat_list[[real_keys[[i]]]]
    if (is.null(fs) || !length(fs)) {
      return(NA_real_)
    }
    (sum(fs >= real_stat[[i]]) + 1.0) / (length(fs) + 1.0)
  }, numeric(1))
  coverage <- ifelse(is.finite(expected) & expected > 0,
                     observed / expected,
                     NA_real_)
  p_fake[observed < as.integer(min_fake_observed)] <- NA_real_
  p_fake[is.finite(coverage) & coverage < min_fake_coverage] <- NA_real_
  key_df$trajectory_stat <- as.numeric(real_stat)
  key_df$p_fake <- p_fake
  key_df$q_fake <- stats::p.adjust(p_fake, method = "BH")
  key_df$n_fake_expected <- expected
  key_df$n_fake_observed <- observed
  key_df$fake_coverage <- coverage
  key_df
}

# ───  Evidence decomposition ────────────────────────────────────────────

#' Per-stratum aligned profile-correlation evidence decomposition.
#'
#' For each row, correlates aligned ligand/receptor micro-bin profiles within
#' the `(source, target, lineage, bin)` stratum.  This is a descriptive
#' profile-concordance layer, not a formal significance test.
#'
#' @param scores_df data.frame from `ComputeTemporalLR` or
#'   `ComputeTemporalLRWithIntermediate$scores`.
#' @param expr `dgCMatrix` (genes × cells), log-normalised.
#' @param celltypes integer vector (1-based).
#' @param lr_db list with `ligand_rows` and `receptor_rows`
#'   (see `ComputeLRScores`).
#' @param pseudotime numeric vector.
#' @param lineage_weights numeric matrix \[n_cells × n_lineages\].
#' @param n_bins integer.
#' @param min_pct numeric.
#' @return data.frame with columns `aligned_r`, `profile_concordance`,
#'   deprecated `z_score` (always `NA`), `n_eff`, `n_microbins`, and
#'   `concordance_available`, parallel to `scores_df` rows.
#' @export
EvidenceDecomp <- function(scores_df, expr, celltypes, lr_db,
                           pseudotime, lineage_weights,
                           n_bins = 12L, min_pct = 0.05) {
  n_bins <- .validate_positive_int(n_bins, "n_bins", min = 2L)
  min_pct <- .validate_unit_interval(min_pct, "min_pct")
  expr <- .prepare_lr_expr(expr)
  lr_db <- .validate_lr_db(lr_db)
  .assert_lr_db_indices(lr_db, nrow(expr))
  stopifnot(is.data.frame(scores_df),
            all(c("pair_idx", "source", "target", "lineage", "bin", "n_eff_cells")
                %in% colnames(scores_df)))

  lw         <- .assert_lineage_weights(lineage_weights, ncol(expr))
  n_lineages <- ncol(lw)
  ct         <- .assert_celltypes(celltypes, ncol(expr))
  ct_levels  <- attr(ct, "levels")
  n_types    <- max(ct)
  pt_vec     <- .assert_pseudotime_1d(pseudotime, ncol(expr))
  score_source <- match(as.integer(scores_df$source), ct_levels)
  score_target <- match(as.integer(scores_df$target), ct_levels)
  if (anyNA(score_source) || anyNA(score_target)) {
    stop("`scores_df` source/target labels must occur in `celltypes`", call. = FALSE)
  }

  raw <- rs_evidence_decomp(
    as.integer(scores_df$pair_idx),
    as.integer(score_source),
    as.integer(score_target),
    as.integer(scores_df$lineage),
    as.integer(scores_df$bin),
    as.numeric(scores_df$n_eff_cells),
    expr,
    ct,
    as.integer(n_types),
    lapply(lr_db$ligand_rows,   as.integer),
    lapply(lr_db$receptor_rows, as.integer),
    pt_vec,
    as.numeric(t(lw)),
    as.integer(n_lineages),
    n_bins,
    min_pct
  )
  data.frame(
    aligned_r            = raw$aligned_r,
    profile_concordance  = raw$aligned_r,
    z_score              = raw$z_score,
    n_eff                = raw$n_eff,
    n_microbins          = raw$n_microbins,
    concordance_available = raw$evidence_sufficient == "TRUE",
    evidence_sufficient  = raw$evidence_sufficient == "TRUE",
    stringsAsFactors = FALSE
  )
}

# ─── Spatial adjacency for ligand-receptor pairs ─────────────────────────

.prepare_spatial_coords <- function(spatial_coords, n_cells) {
  if (!is.numeric(spatial_coords)) {
    stop("`spatial_coords` must be numeric", call. = FALSE)
  }
  coords <- as.matrix(spatial_coords)
  if (nrow(coords) != n_cells) {
    stop("`spatial_coords` must have nrow equal to length(celltypes)",
         call. = FALSE)
  }
  if (ncol(coords) < 2L || ncol(coords) > 3L) {
    stop("`spatial_coords` must have 2 or 3 columns", call. = FALSE)
  }
  if (!all(is.finite(coords))) {
    stop("`spatial_coords` must contain only finite values", call. = FALSE)
  }
  coords
}

.spatial_knn_neighbor_indices <- function(coords, k) {
  n <- nrow(coords)
  k_use <- min(as.integer(k), n - 1L)
  if (k_use < 1L) {
    stop("kNN spatial adjacency requires at least two cells", call. = FALSE)
  }
  idx <- matrix(0L, n, k_use)
  for (i in seq_len(n)) {
    d2 <- colSums((t(coords) - coords[i, ])^2)
    d2[i] <- Inf
    ord <- order(d2)
    idx[i, ] <- ord[seq_len(k_use)]
  }
  idx
}

.spatial_adjacency_dimnames <- function(adjacency) {
  row_ids <- rownames(adjacency)
  col_ids <- colnames(adjacency)
  if (is.null(row_ids)) row_ids <- as.character(seq_len(nrow(adjacency)))
  if (is.null(col_ids)) col_ids <- as.character(seq_len(ncol(adjacency)))
  list(row_ids = row_ids, col_ids = col_ids)
}

.spatial_support_for_pairs <- function(scores_df, adjacency) {
  if (!is.matrix(adjacency) || nrow(adjacency) < 1L || ncol(adjacency) < 1L) {
    stop("`adjacency` must be a non-empty matrix", call. = FALSE)
  }
  if (!"source" %in% colnames(scores_df) || !"target" %in% colnames(scores_df)) {
    stop("`scores_df` must contain `source` and `target` columns",
         call. = FALSE)
  }
  dn <- .spatial_adjacency_dimnames(adjacency)
  src_match <- match(as.character(scores_df$source), dn$row_ids)
  tgt_match <- match(as.character(scores_df$target), dn$col_ids)
  support <- rep(NA_real_, nrow(scores_df))
  ok <- !is.na(src_match) & !is.na(tgt_match)
  if (any(ok)) {
    support[ok] <- mapply(
      function(r, c) adjacency[r, c],
      src_match[ok],
      tgt_match[ok],
      USE.NAMES = FALSE
    )
  }
  support
}

.attach_spatial_columns <- function(scores_df, spatial_adjacency, min_spatial_support) {
  if (is.null(spatial_adjacency)) {
    return(scores_df)
  }
  min_spatial_support <- .validate_unit_interval(
    min_spatial_support, "min_spatial_support"
  )
  scores_df$spatial_support <- .spatial_support_for_pairs(scores_df, spatial_adjacency)
  scores_df$spatially_eligible <- is.finite(scores_df$spatial_support) &
    scores_df$spatial_support >= min_spatial_support
  scores_df
}

#' Compute cell-type spatial adjacency from coordinates.
#'
#' Returns an `n_types` by `n_types` matrix with sender types as rows and
#' receiver types as columns. Each entry is the fraction of cells in the
#' source type that lie within the `k`-nearest neighbourhood of at least one
#' target-type cell (`method = "knn"`) or within `radius` of a target-type cell
#' (`method = "radius"`). Self-pairs (`source == target`) exclude the same cell
#' index so autocrine adjacency does not inflate to 1 at any positive radius.
#' With `binary = TRUE`, entries are `0` or `1` indicating whether any source
#' cell meets the neighbourhood criterion.
#'
#' @param celltypes length-`n_cells` vector of positive 1-based cell-type ids.
#' @param spatial_coords `n_cells` by 2 or 3 numeric matrix of spatial
#'   coordinates.
#' @param method `"knn"` or `"radius"`.
#' @param k number of neighbours per cell for the kNN method.
#' @param radius distance threshold for the radius method (required when
#'   `method = "radius"`).
#' @param binary logical; return binary adjacency instead of fractions.
#' @return numeric matrix with row/column names equal to the observed cell-type
#'   ids.
#' @export
ComputeSpatialAdjacency <- function(celltypes,
                                    spatial_coords,
                                    method = c("knn", "radius"),
                                    k = 6L,
                                    radius = NULL,
                                    binary = FALSE) {
  method <- match.arg(method)
  n_cells <- length(celltypes)
  ct <- .assert_celltypes(celltypes, n_cells)
  type_levels <- attr(ct, "levels")
  n_types <- length(type_levels)
  coords <- .prepare_spatial_coords(spatial_coords, n_cells)
  adj <- matrix(0, n_types, n_types)
  type_cells <- lapply(seq_len(n_types), function(t) which(ct == t))

  if (method == "knn") {
    k <- as.integer(k)
    if (k < 1L) {
      stop("`k` must be at least 1", call. = FALSE)
    }
    nn_idx <- .spatial_knn_neighbor_indices(coords, k)
    for (src in seq_len(n_types)) {
      src_cells <- type_cells[[src]]
      if (!length(src_cells)) next
      for (tgt in seq_len(n_types)) {
        tgt_cells <- type_cells[[tgt]]
        if (!length(tgt_cells)) {
          adj[src, tgt] <- 0
          next
        }
        near <- vapply(
          src_cells,
          function(i) any(nn_idx[i, ] %in% tgt_cells),
          logical(1L)
        )
        adj[src, tgt] <- if (isTRUE(binary)) {
          as.numeric(any(near))
        } else {
          mean(near)
        }
      }
    }
  } else {
    if (is.null(radius) || !is.finite(radius) || radius <= 0) {
      stop("`radius` must be a positive finite value for method = \"radius\"",
           call. = FALSE)
    }
    r2 <- radius * radius
    for (src in seq_len(n_types)) {
      src_cells <- type_cells[[src]]
      if (!length(src_cells)) next
      for (tgt in seq_len(n_types)) {
        tgt_cells <- type_cells[[tgt]]
        if (!length(tgt_cells)) {
          adj[src, tgt] <- 0
          next
        }
        tgt_coords <- coords[tgt_cells, , drop = FALSE]
        near <- vapply(
          src_cells,
          function(i) {
            if (src == tgt) {
              other_tgt <- tgt_cells[tgt_cells != i]
              if (!length(other_tgt)) {
                return(FALSE)
              }
              tgt_coords <- coords[other_tgt, , drop = FALSE]
            }
            d2 <- colSums((t(tgt_coords) - coords[i, ])^2)
            any(d2 <= r2)
          },
          logical(1L)
        )
        adj[src, tgt] <- if (isTRUE(binary)) {
          as.numeric(any(near))
        } else {
          mean(near)
        }
      }
    }
  }

  dimnames(adj) <- list(as.character(type_levels), as.character(type_levels))
  attr(adj, "method") <- method
  if (method == "knn") attr(adj, "k") <- as.integer(k)
  if (method == "radius") attr(adj, "radius") <- radius
  attr(adj, "binary") <- isTRUE(binary)
  adj
}

.permute_celltypes_within_blocks <- function(celltypes, block_labels) {
  n_cells <- length(celltypes)
  ct <- .assert_celltypes(celltypes, n_cells)
  blocks <- .prepare_block_labels(block_labels, n_cells)
  out <- ct
  for (b in unique(blocks)) {
    idx <- which(blocks == b)
    if (length(idx) >= 2L) {
      out[idx] <- ct[sample(idx)]
    }
  }
  attr(out, "levels") <- attr(ct, "levels")
  out
}

.permute_celltypes_preserving_counts <- function(celltypes) {
  ct <- .assert_celltypes(celltypes, length(celltypes))
  perm <- sample.int(length(ct))
  shuffled <- ct[perm]
  attr(shuffled, "levels") <- attr(ct, "levels")
  shuffled
}

.permute_celltypes_size_matched <- function(celltypes,
                                            spatial_coords,
                                            n_strata = 4L) {
  n_cells <- length(celltypes)
  ct <- .assert_celltypes(celltypes, n_cells)
  coords <- .prepare_spatial_coords(spatial_coords, n_cells)
  n_strata <- as.integer(n_strata)
  if (n_strata < 2L) {
    stop("`size_match_strata` must be at least 2", call. = FALSE)
  }
  x_breaks <- unique(stats::quantile(
    coords[, 1L], probs = seq(0, 1, length.out = n_strata + 1L), na.rm = TRUE
  ))
  y_breaks <- unique(stats::quantile(
    coords[, 2L], probs = seq(0, 1, length.out = n_strata + 1L), na.rm = TRUE
  ))
  if (length(x_breaks) < 2L || length(y_breaks) < 2L) {
    return(.permute_celltypes_preserving_counts(celltypes))
  }
  x_bin <- cut(coords[, 1L], breaks = x_breaks, include.lowest = TRUE, labels = FALSE)
  y_bin <- cut(coords[, 2L], breaks = y_breaks, include.lowest = TRUE, labels = FALSE)
  stratum <- paste(x_bin, y_bin, sep = ":")
  out <- ct
  for (s in unique(stratum)) {
    idx <- which(stratum == s)
    if (length(idx) >= 2L) {
      out[idx] <- ct[sample(idx)]
    }
  }
  attr(out, "levels") <- attr(ct, "levels")
  out
}

.spatial_adjacency_params_from_matrix <- function(adjacency) {
  method <- attr(adjacency, "method")
  if (is.null(method)) {
    stop("`adjacency` must carry a `method` attribute from ",
         "`ComputeSpatialAdjacency()`", call. = FALSE)
  }
  list(
    method = method,
    k = attr(adjacency, "k"),
    radius = attr(adjacency, "radius"),
    binary = isTRUE(attr(adjacency, "binary"))
  )
}

#' Permutation null for cell-type spatial adjacency.
#'
#' Shuffles cell-type labels while preserving coordinates and recomputes
#' adjacency under `null_mode`. Returns observed adjacency plus null summary
#' statistics for downstream normalization.
#'
#' @param celltypes length-`n_cells` vector of positive 1-based cell-type ids.
#' @param spatial_coords `n_cells` by 2 or 3 numeric matrix of spatial
#'   coordinates.
#' @param method,k,radius,binary passed to [ComputeSpatialAdjacency()].
#' @param null_mode `"label_permutation"` shuffles labels globally preserving
#'   type counts; `"size_matched"` shuffles within spatial quantile bins;
#'   `"block_label_permutation"` shuffles labels within each `block_labels`
#'   stratum (requires `block_labels`).
#' @param n_perm number of null draws (default 199).
#' @param seed optional random seed.
#' @param size_match_strata number of quantile bins per axis for
#'   `null_mode = "size_matched"`.
#' @param block_labels length-`n_cells` donor/sample labels; required when
#'   `null_mode = "block_label_permutation"`.
#' @return list with `observed`, `null_mean`, `null_sd`, `null_ge_count`,
#'   upper-tail `p_spatial`, family-wise BH `q_spatial`, legacy alias
#'   `null_percentile`, `null_mode`, `n_perm`, and adjacency metadata.
#' @export
ComputeSpatialAdjacencyNull <- function(celltypes,
                                        spatial_coords,
                                        method = c("knn", "radius"),
                                        k = 6L,
                                        radius = NULL,
                                        binary = FALSE,
                                        null_mode = c(
                                          "label_permutation", "size_matched",
                                          "block_label_permutation"
                                        ),
                                        n_perm = 199L,
                                        seed = NULL,
                                        size_match_strata = 4L,
                                        block_labels = NULL) {
  method <- match.arg(method)
  null_mode <- match.arg(null_mode)
  n_perm <- as.integer(n_perm)
  if (n_perm < 1L) {
    stop("`n_perm` must be at least 1", call. = FALSE)
  }
  if (!is.null(seed)) {
    seed <- .validate_seed(seed)
    set.seed(seed)
  }

  observed <- ComputeSpatialAdjacency(
    celltypes, spatial_coords, method = method, k = k,
    radius = radius, binary = binary
  )
  n_types <- nrow(observed)
  null_sum <- matrix(0, n_types, n_types)
  null_sq_sum <- matrix(0, n_types, n_types)
  null_ge <- matrix(0L, n_types, n_types)
  ct_base <- .assert_celltypes(celltypes, length(celltypes))
  if (null_mode == "block_label_permutation") {
    if (is.null(block_labels)) {
      stop("`block_labels` is required when `null_mode = \"block_label_permutation\"`",
           call. = FALSE)
    }
    .prepare_block_labels(block_labels, length(celltypes))
  }
  perm_fn <- if (null_mode == "label_permutation") {
    function(ct) .permute_celltypes_preserving_counts(ct)
  } else if (null_mode == "block_label_permutation") {
    function(ct) .permute_celltypes_within_blocks(ct, block_labels)
  } else {
    function(ct) {
      .permute_celltypes_size_matched(
        ct, spatial_coords, n_strata = size_match_strata
      )
    }
  }

  for (p in seq_len(n_perm)) {
    perm_ct <- perm_fn(ct_base)
    adj_p <- ComputeSpatialAdjacency(
      perm_ct, spatial_coords, method = method, k = k,
      radius = radius, binary = binary
    )
    null_sum <- null_sum + adj_p
    null_sq_sum <- null_sq_sum + adj_p^2
    null_ge <- null_ge + (adj_p >= observed)
  }

  null_mean <- null_sum / n_perm
  null_var <- pmax(null_sq_sum / n_perm - null_mean^2, 0)
  null_sd <- sqrt(null_var)
  p_spatial <- (null_ge + 1L) / (n_perm + 1L)
  q_spatial <- matrix(
    stats::p.adjust(as.numeric(p_spatial), method = "BH"),
    nrow = nrow(p_spatial), ncol = ncol(p_spatial)
  )
  dimnames(null_mean) <- dimnames(observed)
  dimnames(null_sd) <- dimnames(null_mean)
  dimnames(p_spatial) <- dimnames(null_mean)
  dimnames(q_spatial) <- dimnames(null_mean)

  list(
    observed = observed,
    null_mean = null_mean,
    null_sd = null_sd,
    null_ge_count = null_ge,
    p_spatial = p_spatial,
    q_spatial = q_spatial,
    null_percentile = p_spatial,
    null_mode = null_mode,
    n_perm = n_perm,
    method = method,
    k = if (method == "knn") as.integer(k) else NULL,
    radius = if (method == "radius") radius else NULL,
    binary = isTRUE(binary),
    size_match_strata = if (null_mode == "size_matched") {
      as.integer(size_match_strata)
    } else {
      NULL
    }
  )
}

#' Normalize spatial adjacency against a permutation null.
#'
#' Converts observed adjacency to z-scores and upper-tail spatial p/q values using
#' output from [ComputeSpatialAdjacencyNull()] or by computing a fresh null when
#' `celltypes` and `spatial_coords` are supplied.
#'
#' @param adjacency observed adjacency matrix (defaults to `null_result$observed`
#'   when present).
#' @param null_result optional list from [ComputeSpatialAdjacencyNull()].
#' @param celltypes,spatial_coords required when `null_result` is NULL.
#' @param method,k,radius,binary,null_mode,n_perm,seed,size_match_strata passed
#'   to [ComputeSpatialAdjacencyNull()] when a null is computed on the fly.
#' @param block_labels optional block labels used by the block-restricted null.
#' @return list with `adjacency`, `z_score`, `p_spatial`, `q_spatial`, legacy
#'   alias `null_percentile`, `null_mean`, `null_sd`, and `null_mode`.
#' @export
NormalizeSpatialAdjacency <- function(adjacency = NULL,
                                      null_result = NULL,
                                      celltypes = NULL,
                                      spatial_coords = NULL,
                                      method = c("knn", "radius"),
                                      k = 6L,
                                      radius = NULL,
                                      binary = FALSE,
                                      null_mode = c(
                                        "label_permutation", "size_matched",
                                        "block_label_permutation"
                                      ),
                                      n_perm = 199L,
                                      seed = NULL,
                                      size_match_strata = 4L,
                                      block_labels = NULL) {
  if (is.null(null_result)) {
    if (is.null(celltypes) || is.null(spatial_coords)) {
      stop(
        "Provide `null_result` or both `celltypes` and `spatial_coords`",
        call. = FALSE
      )
    }
    if (!is.null(adjacency)) {
      params <- .spatial_adjacency_params_from_matrix(adjacency)
      method <- params$method
      k <- params$k
      radius <- params$radius
      binary <- params$binary
    }
    null_result <- ComputeSpatialAdjacencyNull(
      celltypes = celltypes,
      spatial_coords = spatial_coords,
      method = method,
      k = k,
      radius = radius,
      binary = binary,
      null_mode = null_mode,
      n_perm = n_perm,
      seed = seed,
      size_match_strata = size_match_strata,
      block_labels = block_labels
    )
  }
  obs <- if (!is.null(null_result$observed)) {
    null_result$observed
  } else {
    if (is.null(adjacency)) {
      stop("`adjacency` is required when `null_result$observed` is missing",
           call. = FALSE)
    }
    adjacency
  }
  z <- (obs - null_result$null_mean) /
    pmax(null_result$null_sd, .Machine$double.eps)
  list(
    adjacency = obs,
    z_score = z,
    p_spatial = if (!is.null(null_result$p_spatial)) {
      null_result$p_spatial
    } else {
      null_result$null_percentile
    },
    q_spatial = if (!is.null(null_result$q_spatial)) {
      null_result$q_spatial
    } else {
      matrix(
        stats::p.adjust(as.numeric(null_result$null_percentile), method = "BH"),
        nrow = nrow(null_result$null_percentile),
        ncol = ncol(null_result$null_percentile),
        dimnames = dimnames(null_result$null_percentile)
      )
    },
    null_percentile = if (!is.null(null_result$p_spatial)) {
      null_result$p_spatial
    } else {
      null_result$null_percentile
    },
    null_mean = null_result$null_mean,
    null_sd = null_result$null_sd,
    null_mode = null_result$null_mode,
    n_perm = null_result$n_perm
  )
}

.spatial_null_metrics_for_pairs <- function(scores_df, norm_result) {
  z_mat <- norm_result$z_score
  p_mat <- norm_result$p_spatial
  q_mat <- norm_result$q_spatial
  dn <- .spatial_adjacency_dimnames(z_mat)
  src_match <- match(as.character(scores_df$source), dn$row_ids)
  tgt_match <- match(as.character(scores_df$target), dn$col_ids)
  z_out <- rep(NA_real_, nrow(scores_df))
  p_out <- rep(NA_real_, nrow(scores_df))
  q_out <- rep(NA_real_, nrow(scores_df))
  ok <- !is.na(src_match) & !is.na(tgt_match)
  if (any(ok)) {
    z_out[ok] <- mapply(
      function(r, c) z_mat[r, c],
      src_match[ok],
      tgt_match[ok],
      USE.NAMES = FALSE
    )
    p_out[ok] <- mapply(
      function(r, c) p_mat[r, c],
      src_match[ok],
      tgt_match[ok],
      USE.NAMES = FALSE
    )
    q_out[ok] <- mapply(
      function(r, c) q_mat[r, c],
      src_match[ok],
      tgt_match[ok],
      USE.NAMES = FALSE
    )
  }
  list(z_score = z_out, p_spatial = p_out, q_spatial = q_out)
}

#' Flag or filter LR pairs by spatial adjacency support.
#'
#' Looks up `adjacency[source, target]` for each row (sender row, receiver
#' column) and adds `spatial_support` and `spatially_supported`. When
#' `null_mode != "raw"`, also computes permutation-based `spatial_support_z`
#' and `spatial_support_null_pctl` via [ComputeSpatialAdjacencyNull()].
#'
#' @param scores_df data.frame with `source` and `target` cell-type columns.
#' @param adjacency matrix from `ComputeSpatialAdjacency()`.
#' @param min_spatial_support minimum adjacency fraction for support.
#' @param filter when `TRUE`, drop rows that fail the spatial support gate.
#' @param null_mode `"raw"` (default), `"label_permutation"`, or
#'   `"size_matched"`.
#' @param celltypes,spatial_coords required when `null_mode != "raw"`.
#' @param method,k,radius,binary adjacency settings for null computation; when
#'   omitted, inferred from `adjacency` attributes.
#' @param n_perm,seed,size_match_strata passed to [ComputeSpatialAdjacencyNull()].
#' @param block_labels optional block labels used by the block-restricted null.
#' @param min_spatial_z optional z-score floor; when set, `spatially_supported`
#'   also requires `spatial_support_z >= min_spatial_z`.
#' @param max_spatial_q BH threshold for calibrated spatial support in non-raw
#'   null modes.
#' @return data.frame with spatial support columns; filtered when requested.
#' @export
FilterLRBySpatialSupport <- function(scores_df,
                                     adjacency,
                                     min_spatial_support = 0.05,
                                     filter = FALSE,
                                     null_mode = c("raw", "label_permutation",
                                                   "size_matched",
                                                   "block_label_permutation"),
                                     celltypes = NULL,
                                     spatial_coords = NULL,
                                     method = c("knn", "radius"),
                                     k = 6L,
                                     radius = NULL,
                                     binary = FALSE,
                                     n_perm = 199L,
                                     seed = NULL,
                                     size_match_strata = 4L,
                                     block_labels = NULL,
                                     min_spatial_z = NULL,
                                     max_spatial_q = 0.05) {
  stopifnot(is.data.frame(scores_df))
  null_mode <- match.arg(null_mode)
  min_spatial_support <- .validate_unit_interval(
    min_spatial_support, "min_spatial_support"
  )
  if (!is.null(min_spatial_z) && !is.finite(min_spatial_z)) {
    stop("`min_spatial_z` must be finite when supplied", call. = FALSE)
  }
  out <- scores_df
  out$spatial_support <- .spatial_support_for_pairs(out, adjacency)
  out$spatially_supported <- is.finite(out$spatial_support) &
    out$spatial_support >= min_spatial_support

  if (null_mode != "raw") {
    if (is.null(celltypes) || is.null(spatial_coords)) {
      stop(
        "`celltypes` and `spatial_coords` are required when ",
        "`null_mode != \"raw\"`", call. = FALSE
      )
    }
    params <- .spatial_adjacency_params_from_matrix(adjacency)
    null_result <- ComputeSpatialAdjacencyNull(
      celltypes = celltypes,
      spatial_coords = spatial_coords,
      method = params$method,
      k = params$k,
      radius = params$radius,
      binary = params$binary,
      null_mode = null_mode,
      n_perm = n_perm,
      seed = seed,
      size_match_strata = size_match_strata,
      block_labels = block_labels
    )
    norm <- NormalizeSpatialAdjacency(null_result = null_result)
    pair_metrics <- .spatial_null_metrics_for_pairs(out, norm)
    out$spatial_support_z <- pair_metrics$z_score
    out$p_spatial <- pair_metrics$p_spatial
    out$q_spatial <- pair_metrics$q_spatial
    out$spatial_support_null_pctl <- out$p_spatial
    out$spatial_null_mode <- null_mode
    max_spatial_q <- .validate_unit_interval(max_spatial_q, "max_spatial_q")
    out$spatially_supported <- out$spatially_supported &
      is.finite(out$q_spatial) & out$q_spatial <= max_spatial_q
    if (!is.null(min_spatial_z)) {
      out$spatially_supported <- out$spatially_supported &
        is.finite(out$spatial_support_z) &
        out$spatial_support_z >= min_spatial_z
    }
  }

  if (isTRUE(filter)) {
    out <- out[out$spatially_supported, , drop = FALSE]
  }
  out
}

# ─── MergeEvidence helpers ─────────────────────────────────────────────

.validate_trajectory_keys <- function(df, name, by) {
  stopifnot(is.data.frame(df))
  missing <- setdiff(by, colnames(df))
  if (length(missing) > 0L) {
    stop(name, " must contain columns: ", paste(missing, collapse = ", "))
  }
  key_df <- df[, by, drop = FALSE]
  if (anyDuplicated(key_df)) {
    stop(name, " contains duplicate trajectory keys; keys must be unique", call. = FALSE)
  }
  invisible(df)
}

.merge_calibration_table <- function(scores_df, calibration, by, layer) {
  p_col <- if (layer == "fake") "p_fake" else "p_temporal"
  q_col <- if (layer == "fake") "q_fake" else "q_temporal"
  audit_cols <- if (layer == "fake") {
    c(
      "n_fake_expected", "n_fake_observed", "fake_coverage",
      "n_fake_expected_generic", "n_fake_expected_edge",
      "fake_coverage_generic", "fake_coverage_edge",
      "identity_calibration_evaluable", "complex_pair",
      .fake_identity_detail_cols()
    )
  } else {
    c("n_perm_expected", "n_perm_observed", "perm_coverage")
  }
  quality_cols <- c(
    "trajectory_evaluable", "n_valid_bins", "valid_bin_fraction",
    "n_expression_supported_bins", "n_eff_at_max", "n_eff_at_min",
    "valid_at_max", "valid_at_min"
  )
  config_cols <- c(
    "analysis_data_hash", "scoring_config_hash", "calibration_config_hash",
    "fake_null_config_hash", "temporal_null_config_hash",
    "hill_scale", "trajectory_stat_name"
  )
  layer_meta_cols <- "calibration_hypothesis"
  extra <- intersect(
    colnames(calibration),
    c(p_col, q_col, quality_cols, audit_cols, config_cols, layer_meta_cols)
  )
  cal <- calibration[, c(by, extra), drop = FALSE]
  cal <- .validate_trajectory_keys(cal, paste0(layer, "_calibration"), by)
  if (!q_col %in% colnames(cal) && p_col %in% colnames(cal)) {
    cal[[q_col]] <- stats::p.adjust(cal[[p_col]], method = "BH")
  }
  rename_cols <- setdiff(colnames(cal), by)
  colnames(cal)[match(rename_cols, colnames(cal))] <- paste0(
    rename_cols, ".", layer
  )
  p_col2 <- paste0(p_col, ".", layer)
  q_col2 <- paste0(q_col, ".", layer)
  n_before <- nrow(scores_df)
  merged <- merge(scores_df, cal, by = by, all.x = TRUE, sort = FALSE)
  if (nrow(merged) != n_before) {
    stop(
      "MergeEvidence: merge with ", layer, " calibration changed row count ",
      "(", n_before, " -> ", nrow(merged), "); check for duplicate keys",
      call. = FALSE
    )
  }
  promote_cols <- c(p_col, q_col, audit_cols, config_cols)
  for (col in promote_cols) {
    col2 <- paste0(col, ".", layer)
    if (col %in% colnames(merged) && col2 %in% colnames(merged)) {
      merged[[col]] <- merged[[col2]]
      merged[[col2]] <- NULL
    } else if (col2 %in% colnames(merged)) {
      names(merged)[names(merged) == col2] <- col
    }
  }
  merged
}

# ───  Multi-evidence merge (intersection gating) ────────────────────────

#' Merge structural, FakeLR, and temporal-null evidence layers.
#'
#' Structural gates (`n_eff_cells`, `assignment_confidence`) are quality
#' criteria for lineage-assignment support. `fate_entropy` / `fate_ambiguity`
#' (alias of upstream `lineage_entropy`) retain the biological commitment signal
#' used for role scoring and bifurcation context. When no explicit
#' `assignment_confidence` or `lineage_weight_max` is supplied,
#' `assignment_confidence` is `NA` with `assignment_confidence_source = "missing"`
#' and `pass_assignment` is `NA` (structural gate not applied). `pass_entropy`
#' remains a separate biological fate-commitment check (`fate_entropy <=
#' max_entropy`) and does not substitute for assignment confidence.
#'
#' FakeLR tests LR identity support. Temporal permutation tests whether the
#' trajectory pattern exceeds a within-cell-type/lineage time null.
#'
#' Calibration tables are merged by trajectory key (`pair_idx`, `source`,
#' `target`, `lineage` by default). Trajectory-level q-values are computed
#' once on the calibration table and copied to bin rows; they are never
#' recomputed after bin expansion.
#'
#' `prioritized_candidate` is TRUE when structurally (and spatially, if requested)
#' eligible and at least one calibration layer is supported (OR over layers).
#' `identity_prioritized` and `temporal_prioritized` carry layer-specific
#' prioritization. `high_confidence` is a deprecated alias of `fully_supported`.
#' `fully_supported` uses the conjunction p-value
#' `p_joint = max(p_identity, p_temporal)` where `p_identity` prefers
#' edge-rewiring FakeLR (`p_fake_edge`) over generic context matching
#' (`p_fake_generic`) when both are available. `q_identity` applies one
#' trajectory-level BH correction to that mixed primary family; sub-family
#' `q_fake_edge` / `q_fake_generic` remain diagnostic. One trajectory-level BH
#' correction yields `q_joint`.
#'
#' @param scores_df data.frame from `ComputeTemporalLR` or a trajectory-level
#'   table. Must contain `lineage_entropy` or `fate_entropy`, and optionally
#'   an explicit `assignment_confidence` column.
#' @param fake_calibration optional trajectory-level FakeLR calibration table.
#' @param temporal_calibration optional trajectory-level temporal-null table.
#' @param p_fake deprecated numeric vector; use `fake_calibration` instead.
#' @param p_temporal deprecated numeric vector; use `temporal_calibration`.
#' @param by character vector of trajectory key columns.
#' @param min_n_eff structural gate on Kish effective sample size.
#' @param max_entropy deprecated; used only when `min_assignment_confidence` is
#'   NULL to set the structural confidence floor (`1 - max_entropy`).
#' @param min_assignment_confidence minimum structural assignment confidence;
#'   defaults to `1 - max_entropy` when NULL. Rows with missing
#'   `assignment_confidence` yield `pass_assignment = NA`.
#' @param spatial_adjacency optional `n_types` by `n_types` adjacency matrix
#'   from `ComputeSpatialAdjacency()` (sender row, receiver column). When
#'   supplied, adds `spatial_support` and `spatially_eligible` as orthogonal
#'   annotations; raw adjacency does not alter calibrated discovery status.
#' @param min_spatial_support minimum `adjacency[source, target]` for spatial
#'   eligibility when `spatial_adjacency` is provided.
#' @param alpha q-value threshold.
#' @return data.frame with explicit evidence columns including
#'   `assignment_confidence`, `assignment_confidence_source`, `fate_entropy`,
#'   `fate_ambiguity`, `pass_assignment`, provenance hashes when calibration
#'   `spatial_support`, `spatially_eligible` (when `spatial_adjacency` is given),
#'   `identity_prioritized`, `temporal_variation_prioritized`,
#'   `trajectory_strength_prioritized`, `prioritized_candidate`,
#'   `high_confidence`, `p_joint`, `q_joint`, `fully_supported`,
#'   `p_identity`, `q_identity`, `identity_inference_mode`,
#'   `lr_identity_edge_supported`, `lr_identity_generic_supported`,
#'   `identity_and_strength_supported`, and `calibration_level`.
#' @export
MergeEvidence <- function(scores_df,
                          fake_calibration = NULL,
                          temporal_calibration = NULL,
                          p_fake = NULL,
                          p_temporal = NULL,
                          by = c("pair_idx", "source", "target", "lineage"),
                          min_n_eff = 5.0,
                          max_entropy = 0.8,
                          min_assignment_confidence = NULL,
                          spatial_adjacency = NULL,
                          min_spatial_support = 0.05,
                          alpha = 0.05) {
  stopifnot(is.data.frame(scores_df))
  by <- as.character(by)
  missing <- setdiff(by, colnames(scores_df))
  if (length(missing) > 0L) {
    stop("`scores_df` must contain trajectory key columns: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  required <- c("n_eff_cells")
  if (!all(required %in% colnames(scores_df))) {
    stop("`scores_df` must contain columns: ",
         paste(required, collapse = ", "), call. = FALSE)
  }
  if (!"lineage_entropy" %in% colnames(scores_df) &&
      !"fate_entropy" %in% colnames(scores_df)) {
    stop("`scores_df` must contain `lineage_entropy` or `fate_entropy`",
         call. = FALSE)
  }

  out <- scores_df
  entropy_sem <- .derive_entropy_semantics(out)
  out$assignment_confidence <- entropy_sem$assignment_confidence
  out$assignment_confidence_source <- entropy_sem$assignment_confidence_source
  out$fate_entropy <- entropy_sem$fate_entropy
  out$fate_ambiguity <- entropy_sem$fate_ambiguity
  if (!"lineage_entropy" %in% colnames(out)) {
    out$lineage_entropy <- out$fate_entropy
  }

  if (is.null(min_assignment_confidence)) {
    min_assignment_confidence <- 1 - max_entropy
  } else {
    min_assignment_confidence <- .validate_unit_interval(
      min_assignment_confidence, "min_assignment_confidence"
    )
  }

  out$.original_row_id <- seq_len(nrow(out))

  fake_layer_requested <- !is.null(fake_calibration) || !is.null(p_fake) ||
    "p_fake" %in% colnames(scores_df) || "q_fake" %in% colnames(scores_df)
  temporal_layer_requested <- !is.null(temporal_calibration) || !is.null(p_temporal) ||
    "p_temporal" %in% colnames(scores_df) || "q_temporal" %in% colnames(scores_df)

  fake_hash_source <- if (!is.null(fake_calibration)) fake_calibration else scores_df
  temporal_hash_source <- if (!is.null(temporal_calibration)) temporal_calibration else scores_df
  fake_scoring_hash <- .scoring_hash_from_df(fake_hash_source)
  temporal_scoring_hash <- .scoring_hash_from_df(temporal_hash_source)
  fake_analysis_hash <- .analysis_hash_from_df(fake_hash_source)
  temporal_analysis_hash <- .analysis_hash_from_df(temporal_hash_source)
  if (length(fake_scoring_hash) > 1L || anyNA(fake_scoring_hash)) {
    stop("FakeLR calibration must contain exactly one non-missing scoring configuration hash",
         call. = FALSE)
  }
  if (length(temporal_scoring_hash) > 1L || anyNA(temporal_scoring_hash)) {
    stop("Temporal calibration must contain exactly one non-missing scoring configuration hash",
         call. = FALSE)
  }
  if (fake_layer_requested && temporal_layer_requested &&
      xor(length(fake_scoring_hash) == 1L, length(temporal_scoring_hash) == 1L)) {
    stop(
      "FakeLR and temporal calibration must either both provide scoring configuration hashes or both be legacy tables",
      call. = FALSE
    )
  }
  if (length(fake_scoring_hash) == 1L && length(temporal_scoring_hash) == 1L &&
      !identical(fake_scoring_hash, temporal_scoring_hash)) {
    stop("FakeLR and temporal calibration scoring configuration hashes do not match",
         call. = FALSE)
  }
  if (length(fake_analysis_hash) > 1L || anyNA(fake_analysis_hash)) {
    stop("FakeLR calibration must contain exactly one non-missing analysis data hash",
         call. = FALSE)
  }
  if (length(temporal_analysis_hash) > 1L || anyNA(temporal_analysis_hash)) {
    stop("Temporal calibration must contain exactly one non-missing analysis data hash",
         call. = FALSE)
  }
  if (fake_layer_requested && temporal_layer_requested &&
      length(fake_analysis_hash) == 1L && length(temporal_analysis_hash) == 1L &&
      !identical(fake_analysis_hash, temporal_analysis_hash)) {
    stop("FakeLR and temporal calibration analysis data hashes do not match",
         call. = FALSE)
  }
  if (length(fake_scoring_hash) == 1L) {
    out$scoring_config_hash <- fake_scoring_hash
    out$calibration_config_hash <- fake_scoring_hash
  }
  if (length(fake_analysis_hash) == 1L) {
    out$analysis_data_hash <- fake_analysis_hash
  }
  if (length(fake_scoring_hash) == 0L && length(temporal_scoring_hash) == 1L) {
    out$scoring_config_hash <- temporal_scoring_hash
    out$calibration_config_hash <- temporal_scoring_hash
  }
  if (length(fake_analysis_hash) == 0L && length(temporal_analysis_hash) == 1L) {
    out$analysis_data_hash <- temporal_analysis_hash
  }

  out$pass_n_eff <- out$n_eff_cells >= min_n_eff
  out$pass_assignment <- ifelse(
    out$assignment_confidence_source == "missing",
    NA_real_,
    is.finite(out$assignment_confidence) &
      out$assignment_confidence >= min_assignment_confidence
  )
  out$pass_entropy <- out$fate_entropy <= max_entropy
  assignment_structural_ok <- ifelse(
    out$assignment_confidence_source == "missing",
    TRUE,
    is.finite(out$pass_assignment) & out$pass_assignment
  )
  out$structurally_eligible <- out$pass_n_eff & assignment_structural_ok

  spatial_layer_requested <- !is.null(spatial_adjacency)
  if (spatial_layer_requested) {
    min_spatial_support <- .validate_unit_interval(
      min_spatial_support, "min_spatial_support"
    )
    out$spatial_support <- .spatial_support_for_pairs(out, spatial_adjacency)
    out$spatially_eligible <- is.finite(out$spatial_support) &
      out$spatial_support >= min_spatial_support
  }

  if (!is.null(p_fake) && !is.null(fake_calibration)) {
    stop("Provide either `fake_calibration` or deprecated `p_fake`, not both",
         call. = FALSE)
  }
  if (!is.null(p_temporal) && !is.null(temporal_calibration)) {
    stop("Provide either `temporal_calibration` or deprecated `p_temporal`, not both",
         call. = FALSE)
  }

  if (!is.null(p_fake)) {
    .Deprecated(msg = "MergeEvidence(p_fake=) is deprecated; use fake_calibration=")
    stopifnot(length(p_fake) == nrow(out))
    out$p_fake <- p_fake
    traj_keys <- unique(out[, c(by, "p_fake"), drop = FALSE])
    traj_keys <- .validate_trajectory_keys(traj_keys, "deprecated p_fake keys", by)
    traj_keys$q_fake <- stats::p.adjust(traj_keys$p_fake, method = "BH")
    out$q_fake <- traj_keys$q_fake[match(
      do.call(paste, c(out[, by, drop = FALSE], sep = "\001")),
      do.call(paste, c(traj_keys[, by, drop = FALSE], sep = "\001"))
    )]
  } else if (!is.null(fake_calibration)) {
    out <- .merge_calibration_table(out, fake_calibration, by, "fake")
    out <- .resolve_fake_identity_layer(out, by)
  } else if ("p_fake" %in% colnames(out) || "q_fake" %in% colnames(out)) {
    if ("p_fake" %in% colnames(out) && !"q_fake" %in% colnames(out)) {
      traj_keys <- unique(out[, c(by, "p_fake"), drop = FALSE])
      traj_keys <- .validate_trajectory_keys(traj_keys, "scores_df trajectory p_fake", by)
      traj_keys$q_fake <- stats::p.adjust(traj_keys$p_fake, method = "BH")
      out <- .merge_calibration_table(out, traj_keys, by, "fake")
    } else if ("q_fake" %in% colnames(out)) {
      traj_keys <- unique(out[, intersect(c(by, "p_fake", "q_fake"), colnames(out)), drop = FALSE])
      traj_keys <- .validate_trajectory_keys(traj_keys, "scores_df trajectory q_fake", by)
      out <- .merge_calibration_table(out, traj_keys, by, "fake")
    }
    if (fake_layer_requested) {
      out <- .resolve_fake_identity_layer(out, by)
    }
  }

  if (!is.null(p_temporal)) {
    .Deprecated(msg = "MergeEvidence(p_temporal=) is deprecated; use temporal_calibration=")
    stopifnot(length(p_temporal) == nrow(out))
    out$p_temporal <- p_temporal
    traj_keys <- unique(out[, c(by, "p_temporal"), drop = FALSE])
    traj_keys <- .validate_trajectory_keys(traj_keys, "deprecated p_temporal keys", by)
    traj_keys$q_temporal <- stats::p.adjust(traj_keys$p_temporal, method = "BH")
    out$q_temporal <- traj_keys$q_temporal[match(
      do.call(paste, c(out[, by, drop = FALSE], sep = "\001")),
      do.call(paste, c(traj_keys[, by, drop = FALSE], sep = "\001"))
    )]
  } else if (!is.null(temporal_calibration)) {
    out <- .merge_calibration_table(out, temporal_calibration, by, "temporal")
  } else if ("p_temporal" %in% colnames(out) || "q_temporal" %in% colnames(out)) {
    if ("p_temporal" %in% colnames(out) && !"q_temporal" %in% colnames(out)) {
      traj_keys <- unique(out[, c(by, "p_temporal"), drop = FALSE])
      traj_keys <- .validate_trajectory_keys(traj_keys, "scores_df trajectory p_temporal", by)
      traj_keys$q_temporal <- stats::p.adjust(traj_keys$p_temporal, method = "BH")
      out <- .merge_calibration_table(out, traj_keys, by, "temporal")
    } else if ("q_temporal" %in% colnames(out)) {
      traj_keys <- unique(out[, intersect(c(by, "p_temporal", "q_temporal"), colnames(out)), drop = FALSE])
      traj_keys <- .validate_trajectory_keys(traj_keys, "scores_df trajectory q_temporal", by)
      out <- .merge_calibration_table(out, traj_keys, by, "temporal")
    }
  }

  if (!"p_fake" %in% colnames(out) && !"q_fake" %in% colnames(out)) {
    out$p_fake <- NA_real_
    out$q_fake <- NA_real_
  } else {
    if (!"p_fake" %in% colnames(out)) out$p_fake <- NA_real_
    if (!"q_fake" %in% colnames(out)) out$q_fake <- NA_real_
  }
  if (fake_layer_requested && !"p_identity" %in% colnames(out)) {
    out <- .resolve_fake_identity_layer(out, by)
  }
  if (fake_layer_requested && "identity_calibration_evaluable" %in% colnames(out)) {
    bad <- is.na(out$identity_calibration_evaluable) |
      !out$identity_calibration_evaluable
    out$p_identity[bad] <- NA_real_
    out$q_identity[bad] <- NA_real_
    if ("p_fake_edge" %in% colnames(out)) out$p_fake_edge[bad] <- NA_real_
    if ("q_fake_edge" %in% colnames(out)) out$q_fake_edge[bad] <- NA_real_
    if ("p_replace_ligand" %in% colnames(out)) out$p_replace_ligand[bad] <- NA_real_
    if ("q_replace_ligand" %in% colnames(out)) out$q_replace_ligand[bad] <- NA_real_
    if ("p_replace_receptor" %in% colnames(out)) out$p_replace_receptor[bad] <- NA_real_
    if ("q_replace_receptor" %in% colnames(out)) out$q_replace_receptor[bad] <- NA_real_
    if ("p_edge_iu" %in% colnames(out)) out$p_edge_iu[bad] <- NA_real_
    if ("q_edge_iu" %in% colnames(out)) out$q_edge_iu[bad] <- NA_real_
    if ("p_fake_generic" %in% colnames(out)) out$p_fake_generic[bad] <- NA_real_
    if ("q_fake_generic" %in% colnames(out)) out$q_fake_generic[bad] <- NA_real_
    if ("p_fake_complex" %in% colnames(out)) out$p_fake_complex[bad] <- NA_real_
    if ("q_fake_complex" %in% colnames(out)) out$q_fake_complex[bad] <- NA_real_
    out$p_fake[bad] <- NA_real_
    out$q_fake[bad] <- NA_real_
  }
  if (!"p_temporal" %in% colnames(out) && !"q_temporal" %in% colnames(out)) {
    out$p_temporal <- NA_real_
    out$q_temporal <- NA_real_
  } else {
    if (!"p_temporal" %in% colnames(out)) out$p_temporal <- NA_real_
    if (!"q_temporal" %in% colnames(out)) out$q_temporal <- NA_real_
  }

  p_fake_joint <- .fake_p_for_joint(out)
  out$p_joint <- ifelse(
    is.finite(p_fake_joint) & is.finite(out$p_temporal),
    pmax(p_fake_joint, out$p_temporal),
    NA_real_
  )
  joint_keys <- unique(out[, c(by, "p_joint"), drop = FALSE])
  joint_keys <- .validate_trajectory_keys(joint_keys, "joint calibration", by)
  joint_keys$q_joint <- stats::p.adjust(joint_keys$p_joint, method = "BH")
  out$q_joint <- joint_keys$q_joint[match(
    do.call(paste, c(out[, by, drop = FALSE], sep = "\001")),
    do.call(paste, c(joint_keys[, by, drop = FALSE], sep = "\001"))
  )]
  out$joint_evaluable <- is.finite(out$p_joint) & is.finite(out$q_joint)
  out$joint_supported <- ifelse(
    out$joint_evaluable,
    out$structurally_eligible & out$q_joint <= alpha,
    NA
  )

  out$identity_layer_requested <- fake_layer_requested
  out$temporal_layer_requested <- temporal_layer_requested
  out$identity_evaluable <- if (fake_layer_requested) {
    edge_ev <- .calibration_evaluable(
      if ("p_fake_edge" %in% colnames(out)) out$p_fake_edge else NA_real_,
      if ("q_fake_edge" %in% colnames(out)) out$q_fake_edge else NA_real_
    )
    generic_ev <- .calibration_evaluable(
      if ("p_fake_generic" %in% colnames(out)) out$p_fake_generic else NA_real_,
      if ("q_fake_generic" %in% colnames(out)) out$q_fake_generic else NA_real_
    )
    legacy_ev <- .calibration_evaluable(out$p_fake, out$q_fake)
    edge_ev | generic_ev | legacy_ev
  } else {
    NA
  }
  out$temporal_evaluable <- if (temporal_layer_requested) {
    .calibration_evaluable(out$p_temporal, out$q_temporal)
  } else {
    NA
  }

  out$lr_identity_edge_supported <- if (fake_layer_requested) {
    .calibration_supported(
      .calibration_evaluable(out$p_fake_edge, out$q_fake_edge),
      out$p_fake_edge,
      out$q_fake_edge,
      alpha
    )
  } else {
    NA
  }
  out$lr_identity_generic_supported <- if (fake_layer_requested) {
    .calibration_supported(
      .calibration_evaluable(out$p_fake_generic, out$q_fake_generic),
      out$p_fake_generic,
      out$q_fake_generic,
      alpha
    )
  } else {
    NA
  }
  out$lr_identity_supported <- if (fake_layer_requested) {
    .calibration_supported(
      out$identity_evaluable,
      out$p_identity,
      out$q_identity,
      alpha
    )
  } else {
    NA
  }
  out$pass_fake <- out$lr_identity_supported
  out$temporal_evidence_supported <- .calibration_supported(
    out$temporal_evaluable, out$p_temporal, out$q_temporal, alpha
  )
  temporal_hypothesis <- rep("legacy_unspecified", nrow(out))
  if ("calibration_hypothesis.temporal" %in% colnames(out)) {
    temporal_hypothesis <- as.character(out$calibration_hypothesis.temporal)
  } else if ("calibration_hypothesis" %in% colnames(out)) {
    temporal_hypothesis <- as.character(out$calibration_hypothesis)
  } else if ("trajectory_stat_name" %in% colnames(out)) {
    temporal_hypothesis <- ifelse(
      out$trajectory_stat_name == "amplitude",
      "temporal_variation",
      "trajectory_strength"
    )
  }
  temporal_hypothesis[is.na(temporal_hypothesis)] <- "legacy_unspecified"
  if (any(temporal_hypothesis == "legacy_unspecified") &&
      is.null(temporal_calibration) && !is.null(p_temporal)) {
    temporal_hypothesis[temporal_hypothesis == "legacy_unspecified"] <- "temporal_variation"
  }
  out$temporal_hypothesis <- temporal_hypothesis
  out$temporal_variation_evaluable <- out$temporal_evaluable &
    out$temporal_hypothesis == "temporal_variation"
  out$trajectory_strength_evaluable <- out$temporal_evaluable &
    out$temporal_hypothesis == "trajectory_strength"
  out$temporal_variation_supported <- .calibration_supported(
    out$temporal_variation_evaluable, out$p_temporal, out$q_temporal, alpha
  )
  out$trajectory_strength_supported <- .calibration_supported(
    out$trajectory_strength_evaluable, out$p_temporal, out$q_temporal, alpha
  )
  out$temporal_specificity_supported <- out$temporal_variation_supported
  out$pass_temporal <- out$temporal_evidence_supported
  out$fully_supported <- ifelse(
    out$joint_evaluable & out$temporal_variation_evaluable,
    out$joint_supported,
    NA
  )
  out$identity_and_strength_supported <- ifelse(
    out$joint_evaluable & out$trajectory_strength_evaluable,
    out$joint_supported,
    NA
  )

  structural_gate <- out$structurally_eligible

  out$lr_identity_calibrated <- out$identity_evaluable
  out$temporal_evidence_calibrated <- out$temporal_evaluable
  out$temporal_specificity_calibrated <- out$temporal_variation_evaluable
  out$trajectory_strength_calibrated <- out$trajectory_strength_evaluable

  out$identity_prioritized <- .layer_prioritized(
    structural_gate, out$identity_evaluable, out$lr_identity_supported
  )
  out$temporal_prioritized <- .layer_prioritized(
    structural_gate, out$temporal_evaluable, out$temporal_evidence_supported
  )
  out$temporal_variation_prioritized <- .layer_prioritized(
    structural_gate,
    out$temporal_variation_evaluable,
    out$temporal_variation_supported
  )
  out$trajectory_strength_prioritized <- .layer_prioritized(
    structural_gate,
    out$trajectory_strength_evaluable,
    out$trajectory_strength_supported
  )

  both_evaluable <- fake_layer_requested & temporal_layer_requested &
    out$identity_evaluable & out$temporal_evaluable
  out$dual_layer_evaluable <- both_evaluable
  out$high_confidence <- out$fully_supported
  attr(out$high_confidence, "deprecated_alias_of") <- "fully_supported"

  identity_support <- fake_layer_requested & out$identity_evaluable &
    out$lr_identity_supported %in% TRUE
  variation_support <- temporal_layer_requested & out$temporal_variation_evaluable &
    out$temporal_variation_supported %in% TRUE
  strength_support <- temporal_layer_requested & out$trajectory_strength_evaluable &
    out$trajectory_strength_supported %in% TRUE
  no_calibration_requested <- !fake_layer_requested & !temporal_layer_requested
  out$prioritized_candidate <- structural_gate &
    (identity_support | variation_support | strength_support)

  n_out <- nrow(out)
  fake_only <- rep(fake_layer_requested & !temporal_layer_requested, n_out)
  temporal_only <- rep(!fake_layer_requested & temporal_layer_requested, n_out)
  both_requested <- rep(fake_layer_requested & temporal_layer_requested, n_out)
  out$calibration_level <- rep("attempted_but_unevaluable", n_out)
  out$calibration_level[no_calibration_requested] <- "structural_only"
  out$calibration_level[fake_only & out$identity_evaluable] <-
    "identity_calibrated"
  out$calibration_level[temporal_only & out$temporal_evaluable] <-
    "temporal_calibrated"
  out$calibration_level[both_requested & out$identity_evaluable &
                          out$temporal_evaluable] <- "fully_calibrated"
  out$calibration_level[both_requested &
                          xor(out$identity_evaluable, out$temporal_evaluable)] <-
    "partially_calibrated"

  if (fake_layer_requested) {
    te_col <- "trajectory_evaluable.fake"
    miss_fake <- !out$identity_evaluable
    if (te_col %in% colnames(out)) {
      uneval_fake <- !is.na(out[[te_col]]) & !out[[te_col]] & !out$identity_evaluable
      miss_fake <- is.na(out[[te_col]]) & !out$identity_evaluable
      out$calibration_level[uneval_fake &
                              (!temporal_layer_requested | !out$temporal_evaluable)] <-
        "attempted_but_unevaluable"
      out$calibration_level[miss_fake & !temporal_layer_requested] <-
        "missing_calibration_key"
      out$calibration_level[miss_fake & temporal_layer_requested &
                              out$temporal_evaluable] <- "temporal_calibrated"
    } else {
      out$calibration_level[miss_fake & !temporal_layer_requested] <-
        "missing_calibration_key"
      out$calibration_level[miss_fake & temporal_layer_requested &
                              out$temporal_evaluable] <- "temporal_calibrated"
    }
  }
  if (temporal_layer_requested) {
    te_col <- "trajectory_evaluable.temporal"
    miss_temp <- !out$temporal_evaluable
    if (te_col %in% colnames(out)) {
      uneval_temp <- !is.na(out[[te_col]]) & !out[[te_col]] & !out$temporal_evaluable
      miss_temp <- is.na(out[[te_col]]) & !out$temporal_evaluable
      out$calibration_level[uneval_temp &
                              (!fake_layer_requested | !out$identity_evaluable)] <-
        "attempted_but_unevaluable"
      out$calibration_level[miss_temp & !fake_layer_requested] <-
        "missing_calibration_key"
      out$calibration_level[miss_temp & fake_layer_requested &
                              out$identity_evaluable] <- "identity_calibrated"
    } else {
      out$calibration_level[miss_temp & !fake_layer_requested] <-
        "missing_calibration_key"
      out$calibration_level[miss_temp & fake_layer_requested &
                              out$identity_evaluable] <- "identity_calibrated"
    }
  }

  out <- out[order(out$.original_row_id), , drop = FALSE]
  out$.original_row_id <- NULL
  rownames(out) <- NULL

  out
}

# ─── Sample-level temporal effects ──────────────────────────────────────

.trajectory_key_col <- function(df, by) {
  do.call(paste, c(df[, by, drop = FALSE], sep = "\001"))
}

.pooled_lr_context_vectors <- function(scores_df, ct_levels = NULL) {
  ctx_cols <- c("pair_idx", "source", "target")
  missing <- setdiff(ctx_cols, colnames(scores_df))
  if (length(missing) > 0L) {
    stop("`scores_df` must contain ", paste(ctx_cols, collapse = ", "),
         " to define pooled LR contexts", call. = FALSE)
  }
  ctx <- unique(as.data.frame(scores_df)[, ctx_cols, drop = FALSE])
  if (is.character(ctx$source) || is.factor(ctx$source)) {
    if (is.null(ct_levels)) {
      ct_levels <- unique(c(as.character(ctx$source), as.character(ctx$target)))
    }
    ctx$source <- match(as.character(ctx$source), ct_levels)
    ctx$target <- match(as.character(ctx$target), ct_levels)
  } else {
    ctx$source <- as.integer(ctx$source)
    ctx$target <- as.integer(ctx$target)
  }
  if (any(!is.finite(ctx$source)) || any(!is.finite(ctx$target))) {
    stop("requested contexts contain unknown source/target labels", call. = FALSE)
  }
  list(
    pair_idx = as.integer(ctx$pair_idx),
    source = as.integer(ctx$source),
    target = as.integer(ctx$target)
  )
}

.pooled_trajectory_key_vectors <- function(scores_df, by, ct_levels = NULL) {
  missing <- setdiff(by, colnames(scores_df))
  if (length(missing) > 0L) {
    stop("`scores_df` must contain trajectory key columns: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  traj <- unique(as.data.frame(scores_df)[, by, drop = FALSE])
  if (any(c("source", "target") %in% by) &&
      (is.character(traj$source) || is.factor(traj$source))) {
    if (is.null(ct_levels)) {
      ct_levels <- unique(c(as.character(traj$source), as.character(traj$target)))
    }
    if ("source" %in% by) {
      traj$source <- match(as.character(traj$source), ct_levels)
    }
    if ("target" %in% by) {
      traj$target <- match(as.character(traj$target), ct_levels)
    }
  }
  out <- lapply(traj, function(col) as.integer(col))
  names(out) <- by
  out
}

.donor_bins_from_block_scores <- function(bdf, n_bins) {
  y <- rep(NA_real_, n_bins)
  if (is.null(bdf) || !nrow(bdf)) {
    return(list(y = y, n_valid = 0L, n_eff = NA_real_))
  }
  for (i in seq_len(nrow(bdf))) {
    bin_i <- as.integer(bdf$bin[i])
    if (bin_i < 1L || bin_i > n_bins) {
      next
    }
    if ("structural_valid" %in% colnames(bdf)) {
      if (isTRUE(bdf$structural_valid[i])) {
        y[bin_i] <- as.numeric(bdf$score[i])
      }
    } else if (is.finite(bdf$score[i])) {
      y[bin_i] <- as.numeric(bdf$score[i])
    }
  }
  n_valid <- sum(is.finite(y))
  n_eff <- if ("n_eff_cells" %in% colnames(bdf)) {
    valid_eff <- as.numeric(bdf$n_eff_cells[is.finite(y)])
    if (length(valid_eff)) min(valid_eff, na.rm = TRUE) else NA_real_
  } else {
    NA_real_
  }
  list(y = y, n_valid = n_valid, n_eff = n_eff)
}

.trajectory_stat_from_bins <- function(bin_scores, stat, top_k = 3L) {
  y <- as.numeric(bin_scores)
  if (!length(y) || !any(is.finite(y))) {
    return(NA_real_)
  }
  y[!is.finite(y)] <- NA_real_
  if (stat == "amplitude") {
    return(max(y, na.rm = TRUE) - min(y, na.rm = TRUE))
  }
  if (stat == "max") {
    return(max(y, na.rm = TRUE))
  }
  if (stat == "mean_score") {
    return(mean(y, na.rm = TRUE))
  }
  kk <- min(as.integer(top_k), sum(is.finite(y)))
  if (kk < 1L) {
    return(NA_real_)
  }
  mean(sort(y, decreasing = TRUE, na.last = TRUE)[seq_len(kk)], na.rm = TRUE)
}

.infer_n_bins <- function(scores_df, n_bins = NULL) {
  if (!is.null(n_bins)) {
    return(.validate_positive_int(n_bins, "n_bins", min = 2L))
  }
  if (!"bin" %in% colnames(scores_df)) {
    stop("`scores_df` must contain `bin` or pass explicit `n_bins`", call. = FALSE)
  }
  max(.as_integerish(scores_df$bin, "scores_df$bin"), na.rm = TRUE)
}

.block_lineage_bin_mass <- function(block_cells, pseudotime, lineage_weights,
                                     n_bins, min_pct = 0) {
  lw <- lineage_weights[block_cells, , drop = FALSE]
  pt <- pseudotime[block_cells]
  n_lineages <- ncol(lw)
  bounds <- .compute_lineage_bin_bounds(pt, lw, n_bins)
  mass <- matrix(0, nrow = n_lineages, ncol = n_bins)
  for (l in seq_len(n_lineages)) {
    w <- lw[, l]
    ok <- is.finite(pt) & w >= 1e-8
    if (!any(ok)) {
      next
    }
    pt_l <- pt[ok]
    w_l <- w[ok]
    bin_idx <- vapply(
      pt_l,
      .lineage_cell_bin_index,
      integer(1),
      bounds = bounds[[l]],
      n_bins = n_bins
    )
    for (b in seq_len(n_bins)) {
      mass[l, b] <- sum(w_l[bin_idx == b])
    }
  }
  mass
}

.block_peak_time <- function(block_cells, pseudotime, lineage_weights,
                            lineage_id, peak_bin, n_bins) {
  lw_block <- lineage_weights[block_cells, , drop = FALSE]
  lw <- lw_block[, lineage_id, drop = TRUE]
  pt <- pseudotime[block_cells]
  ok <- is.finite(pt) & lw >= 1e-8
  if (!any(ok)) {
    return(NA_real_)
  }
  pt_l <- pt[ok]
  w_l <- lw[ok]
  bounds <- .compute_lineage_bin_bounds(pt, lw_block, n_bins)[[lineage_id]]
  bin_idx <- vapply(
    pt_l,
    .lineage_cell_bin_index,
    integer(1),
    bounds = bounds,
    n_bins = n_bins
  )
  in_peak <- bin_idx == peak_bin
  if (!any(in_peak)) {
    return(NA_real_)
  }
  stats::weighted.mean(pt_l[in_peak], w_l[in_peak])
}

.merge_trajectory_calibration <- function(by, merged_df = NULL,
                                          fake_calibration = NULL,
                                          temporal_calibration = NULL) {
  keys <- NULL
  if (!is.null(merged_df)) {
    stopifnot(is.data.frame(merged_df))
    missing <- setdiff(by, colnames(merged_df))
    if (length(missing) > 0L) {
      stop("`merged_df` must contain trajectory key columns: ",
           paste(missing, collapse = ", "), call. = FALSE)
    }
    keep <- unique(c(
      by,
      intersect(
        c("p_identity", "q_identity", "p_fake", "q_fake",
          "p_temporal", "q_temporal", "trajectory_evaluable",
          "identity_evaluable", "temporal_evaluable",
          "n_eff_cells", "lineage_entropy"),
        colnames(merged_df)
      )
    ))
    keys <- unique(merged_df[, keep, drop = FALSE])
    keys <- .validate_trajectory_keys(keys, "merged_df", by)
  }
  if (!is.null(fake_calibration)) {
    cal <- fake_calibration[, unique(c(by, intersect(
      c("p_identity", "q_identity", "p_fake", "q_fake",
        "identity_calibration_evaluable"),
      colnames(fake_calibration)
    ))), drop = FALSE]
    cal <- .validate_trajectory_keys(cal, "fake_calibration", by)
    keys <- if (is.null(keys)) cal else merge(keys, cal, by = by, all = TRUE, sort = FALSE)
  }
  if (!is.null(temporal_calibration)) {
    cal <- temporal_calibration[, unique(c(by, intersect(
      c("p_temporal", "q_temporal", "trajectory_evaluable",
        "trajectory_stat", "trajectory_stat_name"),
      colnames(temporal_calibration)
    ))), drop = FALSE]
    cal <- .validate_trajectory_keys(cal, "temporal_calibration", by)
    keys <- if (is.null(keys)) cal else merge(keys, cal, by = by, all = TRUE, sort = FALSE)
  }
  if (!is.null(keys) && "p_identity" %in% colnames(keys)) {
    keys
  } else if (!is.null(keys)) {
    keys$p_identity <- if ("p_fake" %in% colnames(keys)) keys$p_fake else NA_real_
    keys
  } else {
    NULL
  }
}

.donor_rescore_formal_mode <- function(expr, celltypes, pseudotime,
                                        lineage_weights, lr_db, block_labels) {
  !is.null(expr) && !is.null(celltypes) && !is.null(pseudotime) &&
    !is.null(lineage_weights) && !is.null(lr_db) && !is.null(block_labels)
}

.stouffer_combine_p <- function(p, sided = c("one", "two"),
                                signs = NULL, weights = NULL) {
  p <- as.numeric(p)
  sided <- match.arg(sided)
  ok <- is.finite(p) & p > 0 & p <= 1
  k <- sum(ok)
  if (k < 2L) {
    return(NA_real_)
  }
  w <- if (is.null(weights)) {
    rep(1, k)
  } else {
    as.numeric(weights)[ok]
  }
  if (length(w) != k || any(!is.finite(w)) || any(w <= 0)) {
    stop("`weights` must be finite positive values aligned with `p`", call. = FALSE)
  }
  if (sided == "one") {
    q <- stats::qnorm(1 - p[ok])
    z <- sum(w * q) / sqrt(sum(w^2))
    stats::pnorm(z, lower.tail = FALSE)
  } else {
    q <- stats::qnorm(1 - p[ok] / 2)
    if (is.null(signs)) {
      stop("two-sided Stouffer combination requires `signs`", call. = FALSE)
    }
    sgn <- sign(as.numeric(signs)[ok])
    if (any(!is.finite(sgn) | sgn == 0)) {
      stop("two-sided Stouffer `signs` must be finite and non-zero", call. = FALSE)
    }
    z <- sum(w * sgn * q) / sqrt(sum(w^2))
    2 * stats::pnorm(-abs(z))
  }
}

.fisher_combine_p <- function(p) {
  p <- as.numeric(p)
  ok <- is.finite(p) & p > 0 & p <= 1
  k <- sum(ok)
  if (k < 2L) {
    return(NA_real_)
  }
  stats::pchisq(-2 * sum(log(p[ok])), df = 2 * k, lower.tail = FALSE)
}

.block_trajectory_null_median <- function(bin_scores, stat, top_k,
                                         n_draws = 199L, seed = 1L) {
  y <- as.numeric(bin_scores)
  if (!length(y) || !any(is.finite(y))) {
    return(NA_real_)
  }
  n_draws <- .validate_positive_int(n_draws, "n_null_median_draws", min = 1L)
  seed <- .validate_seed(seed)
  fin_idx <- which(is.finite(y))
  n_fin <- length(fin_idx)
  if (n_fin < 1L) {
    return(NA_real_)
  }
  set.seed(seed)
  stats <- vapply(seq_len(n_draws), function(i) {
    perm <- y
    perm[fin_idx] <- y[fin_idx][sample.int(n_fin, n_fin, replace = FALSE)]
    .trajectory_stat_from_bins(perm, stat, top_k)
  }, numeric(1))
  stats <- stats[is.finite(stats)]
  if (length(stats) < 1L) {
    return(NA_real_)
  }
  stats::median(stats)
}

.compute_block_temporal_permutation_table <- function(expr, celltypes, lr_db,
                                                        pseudotime, lineage_weights,
                                                        n_bins, stat, top_k,
                                                        n_perms, seed,
                                                        min_pct,
                                                        min_n_eff,
                                                        min_valid_bins,
                                                        min_valid_fraction,
                                                        permutation_stratum = c(
                                                          "block_celltype_lineage",
                                                          "block_celltype"
                                                        ),
                                                        min_perm_observed = 0L,
                                                        min_perm_coverage = 0,
                                                        requested_trajectories = NULL) {
  permutation_stratum <- match.arg(permutation_stratum)
  out <- ComputeTemporalPermutation(
    expr = expr,
    celltypes = celltypes,
    lr_db = lr_db,
    pseudotime = pseudotime,
    lineage_weights = lineage_weights,
    n_bins = n_bins,
    n_perms = n_perms,
    seed = seed,
    stat = stat,
    top_k = top_k,
    min_pct = min_pct,
    min_n_eff = min_n_eff,
    min_valid_bins = min_valid_bins,
    min_valid_fraction = min_valid_fraction,
    min_perm_observed = min_perm_observed,
    min_perm_coverage = min_perm_coverage,
    block_labels = NULL,
    requested_trajectories = requested_trajectories,
    permutation_stratum = permutation_stratum
  )
  traj <- out$trajectories
  by <- c("pair_idx", "source", "target", "lineage")
  traj$.traj_key <- .trajectory_key_col(traj, by)
  if (!"temporal_null_median" %in% colnames(traj)) {
    traj$temporal_null_median <- NA_real_
  }
  traj
}

.lag_permutation_null <- function(lig, rec, timing_src, timing_tgt, method,
                                  bin_width, max_lag_bins, observed_lag,
                                  n_perm, seed) {
  if (!is.finite(observed_lag)) {
    return(list(
      p_lag = NA_real_, lag_null_median = NA_real_,
      n_lag_perm_observed = 0L, lag_null_mode = "circular_shift"
    ))
  }
  n <- length(lig)
  if (n < 3L || sum(is.finite(lig)) < 2L || sum(is.finite(rec)) < 2L) {
    return(list(
      p_lag = NA_real_, lag_null_median = NA_real_,
      n_lag_perm_observed = 0L, lag_null_mode = "circular_shift"
    ))
  }
  n_perm <- .validate_positive_int(n_perm, "n_lag_perm", min = 1L)
  seed <- .validate_seed(seed)
  set.seed(seed)
  available_shifts <- seq_len(n - 1L)
  sampled_shifts <- if (n_perm >= length(available_shifts)) {
    available_shifts
  } else {
    sample(available_shifts, n_perm, replace = FALSE)
  }
  null_lags <- vapply(sampled_shifts, function(shift) {
    rec_idx <- ((seq_len(n) - 1L + shift) %% n) + 1L
    rec_p <- rec[rec_idx]
    est <- switch(
      method,
      peak = .lag_from_peak(lig, rec_p, timing_src, timing_tgt),
      centroid = .lag_from_centroid(lig, rec_p, timing_src, timing_tgt),
      cross_cor = .lag_from_cross_cor(
        lig, rec_p, timing_src, timing_tgt, bin_width, max_lag_bins
      )
    )
    est$lag
  }, numeric(1))
  null_lags <- null_lags[is.finite(null_lags)]
  if (!length(null_lags)) {
    return(list(
      p_lag = NA_real_, lag_null_median = NA_real_,
      n_lag_perm_observed = 0L, lag_null_mode = "circular_shift"
    ))
  }
  n_ge <- sum(abs(null_lags) >= abs(observed_lag))
  n_null_obs <- length(null_lags)
  list(
    p_lag = (n_ge + 1L) / (n_null_obs + 1L),
    lag_null_median = stats::median(null_lags),
    n_lag_perm_observed = n_null_obs,
    lag_null_mode = "circular_shift"
  )
}

.resolve_identity_meta_defaults <- function(sample_effects, effect_col, method, p_col) {
  if (is.null(effect_col)) {
    if ("identity_effect_null_centered" %in% colnames(sample_effects) &&
        any(is.finite(sample_effects$identity_effect_null_centered))) {
      effect_col <- "identity_effect_null_centered"
    } else {
      effect_col <- "identity_effect"
    }
  }
  list(effect_col = effect_col, method = method, p_col = p_col)
}

.resolve_lag_meta_defaults <- function(sample_effects, effect_col, method, p_col) {
  if (is.null(effect_col)) {
    if ("lag_null_centered" %in% colnames(sample_effects) &&
        any(is.finite(sample_effects$lag_null_centered))) {
      effect_col <- "lag_null_centered"
    } else if ("lag" %in% colnames(sample_effects)) {
      effect_col <- "lag"
    } else if ("signed_lag" %in% colnames(sample_effects)) {
      effect_col <- "signed_lag"
    } else {
      effect_col <- "lag"
    }
  }
  list(effect_col = effect_col, method = method, p_col = p_col)
}

.resolve_meta_effect_col <- function(sample_effects, effect_col) {
  if (effect_col %in% colnames(sample_effects)) {
    vals <- sample_effects[[effect_col]]
    if (any(is.finite(as.numeric(vals)))) {
      return(effect_col)
    }
  }
  if ("temporal_effect_null_centered" %in% colnames(sample_effects) &&
      any(is.finite(sample_effects$temporal_effect_null_centered))) {
    return("temporal_effect_null_centered")
  }
  if ("temporal_effect" %in% colnames(sample_effects) &&
      any(is.finite(sample_effects$temporal_effect))) {
    return("temporal_effect")
  }
  effect_col
}

.warn_random_effects_donor_count <- function(groups, effect_col, se_col, method) {
  if (!identical(method, "random")) {
    return(invisible(NULL))
  }
  n_by_traj <- vapply(groups, function(g) {
    y <- g[[effect_col]]
    se <- g[[se_col]]
    sum(is.finite(as.numeric(y)) & is.finite(as.numeric(se)) & se > 0)
  }, integer(1))
  if (any(n_by_traj > 0L & n_by_traj < 3L, na.rm = TRUE)) {
    warning(
      "Random-effects donor meta is best with n_samples >= 3; ",
      "some trajectories have fewer evaluable donors.",
      call. = FALSE
    )
  }
  invisible(NULL)
}

.fixed_effect_meta <- function(y, se) {
  y <- as.numeric(y)
  se <- as.numeric(se)
  ok <- is.finite(y) & is.finite(se) & se > 0
  k <- sum(ok)
  if (k < 1L) {
    return(list(meta_effect = NA_real_, meta_se = NA_real_,
                meta_z = NA_real_, meta_p = NA_real_,
                n_samples = 0L))
  }
  if (k == 1L) {
    return(list(
      meta_effect = y[ok][1L],
      meta_se = se[ok][1L],
      meta_z = y[ok][1L] / se[ok][1L],
      meta_p = 2 * stats::pnorm(-abs(y[ok][1L] / se[ok][1L])),
      n_samples = 1L
    ))
  }
  w <- 1 / se[ok]^2
  est <- sum(w * y[ok]) / sum(w)
  se_est <- sqrt(1 / sum(w))
  z <- est / se_est
  list(
    meta_effect = est,
    meta_se = se_est,
    meta_z = z,
    meta_p = 2 * stats::pnorm(-abs(z)),
    n_samples = k
  )
}

.heterogeneity_stats <- function(y, se, meta_effect) {
  y <- as.numeric(y)
  se <- as.numeric(se)
  ok <- is.finite(y) & is.finite(se) & se > 0
  k <- sum(ok)
  if (k < 2L || !is.finite(meta_effect)) {
    return(list(Q = NA_real_, I2 = NA_real_, df = NA_integer_))
  }
  Q <- sum(((y[ok] - meta_effect) / se[ok])^2)
  df <- k - 1L
  I2 <- max(0, (Q - df) / Q)
  list(Q = Q, I2 = I2, df = df)
}

.dersimonian_laird_tau2 <- function(y, se, fixed_effect = NULL) {
  y <- as.numeric(y)
  se <- as.numeric(se)
  ok <- is.finite(y) & is.finite(se) & se > 0
  k <- sum(ok)
  if (k < 2L) {
    return(NA_real_)
  }
  vi <- se[ok]^2
  w <- 1 / vi
  if (is.null(fixed_effect) || !is.finite(fixed_effect)) {
    fixed_effect <- sum(w * y[ok]) / sum(w)
  }
  Q <- sum(w * (y[ok] - fixed_effect)^2)
  c_val <- sum(w) - sum(w^2) / sum(w)
  if (!is.finite(c_val) || c_val <= 0) {
    return(0)
  }
  max(0, (Q - (k - 1L)) / c_val)
}

.random_effect_meta <- function(y, se, tau2 = NULL) {
  y <- as.numeric(y)
  se <- as.numeric(se)
  ok <- is.finite(y) & is.finite(se) & se > 0
  k <- sum(ok)
  if (k < 1L) {
    return(list(
      meta_effect = NA_real_, meta_se = NA_real_,
      meta_z = NA_real_, meta_p = NA_real_,
      tau2 = NA_real_, n_samples = 0L
    ))
  }
  if (is.null(tau2) || !is.finite(tau2)) {
    fe <- .fixed_effect_meta(y, se)
    tau2 <- .dersimonian_laird_tau2(y, se, fe$meta_effect)
  }
  tau2 <- max(0, as.numeric(tau2))
  if (k == 1L) {
    return(list(
      meta_effect = y[ok][1L],
      meta_se = sqrt(se[ok][1L]^2 + tau2),
      meta_z = y[ok][1L] / sqrt(se[ok][1L]^2 + tau2),
      meta_p = 2 * stats::pnorm(-abs(y[ok][1L] / sqrt(se[ok][1L]^2 + tau2))),
      tau2 = tau2,
      n_samples = 1L
    ))
  }
  vi <- se[ok]^2
  w <- 1 / (vi + tau2)
  est <- sum(w * y[ok]) / sum(w)
  se_est <- sqrt(1 / sum(w))
  z <- est / se_est
  list(
    meta_effect = est,
    meta_se = se_est,
    meta_z = z,
    meta_p = 2 * stats::pnorm(-abs(z)),
    tau2 = tau2,
    n_samples = k
  )
}

.sign_consistency <- function(y) {
  y <- as.numeric(y)
  y <- y[is.finite(y) & y != 0]
  if (length(y) < 1L) {
    return(NA_real_)
  }
  pos <- sum(y > 0)
  neg <- sum(y < 0)
  max(pos, neg) / length(y)
}

# Sample-level effect standard error from bin scores or heuristic fallback.
.sample_effect_se <- function(temporal_effect, n_eff_cells, n_valid_bins,
                              bin_scores = NULL) {
  if (!is.null(bin_scores)) {
    valid <- as.numeric(bin_scores)
    valid <- valid[is.finite(valid)]
    if (length(valid) >= 2L) {
      bin_var <- stats::var(valid)
      n_v <- length(valid)
      if (is.finite(bin_var) && n_v >= 2L) {
        return(pmax(sqrt(bin_var / n_v), 1e-3))
      }
    }
  }
  base <- pmax(abs(temporal_effect), 1e-6)
  denom <- pmax(1, sqrt(pmax(n_eff_cells, 1) * pmax(n_valid_bins, 1)))
  pmax(base / denom, 1e-3)
}

.extract_calibration_p_identity <- function(cal_df, traj_key, by) {
  if (is.null(cal_df) || nrow(cal_df) == 0L) {
    return(NA_real_)
  }
  cal_match <- .trajectory_key_col(cal_df, by)
  idx <- match(traj_key, cal_match)
  if (is.na(idx)) {
    return(NA_real_)
  }
  row <- cal_df[idx, , drop = FALSE]
  resolved <- .resolve_fake_identity_layer(row, by)
  if ("p_identity" %in% colnames(resolved) && is.finite(resolved$p_identity[1L])) {
    return(resolved$p_identity[1L])
  }
  if ("p_fake" %in% colnames(row) && is.finite(row$p_fake[1L])) {
    return(row$p_fake[1L])
  }
  if ("p_fake_complex" %in% colnames(row) && is.finite(row$p_fake_complex[1L])) {
    return(row$p_fake_complex[1L])
  }
  NA_real_
}

.extract_calibration_fake_null_median <- function(cal_df, traj_key, by) {
  if (is.null(cal_df) || nrow(cal_df) == 0L) {
    return(NA_real_)
  }
  cal_match <- .trajectory_key_col(cal_df, by)
  idx <- match(traj_key, cal_match)
  if (is.na(idx)) {
    return(NA_real_)
  }
  row <- cal_df[idx, , drop = FALSE]
  if ("fake_null_median" %in% colnames(row) &&
      is.finite(row$fake_null_median[1L])) {
    return(as.numeric(row$fake_null_median[1L]))
  }
  NA_real_
}

.extract_calibration_trajectory_stat <- function(cal_df, traj_key, by) {
  if (is.null(cal_df) || nrow(cal_df) == 0L) {
    return(NA_real_)
  }
  cal_match <- .trajectory_key_col(cal_df, by)
  idx <- match(traj_key, cal_match)
  if (is.na(idx)) {
    return(NA_real_)
  }
  row <- cal_df[idx, , drop = FALSE]
  if ("trajectory_stat" %in% colnames(row) &&
      is.finite(row$trajectory_stat[1L])) {
    return(as.numeric(row$trajectory_stat[1L]))
  }
  NA_real_
}

.compute_block_fake_calibration <- function(expr, celltypes, lr_db,
                                            candidate_fake_genes,
                                            pseudotime, lineage_weights,
                                            n_bins, n_fake_per_real,
                                            complex_null_mode,
                                            fake_null_mode, seed,
                                            n_strata, caliper_k,
                                            max_match_distance, min_pct,
                                            min_n_eff, min_valid_bins,
                                            min_valid_fraction,
                                            min_fake_observed,
                                            min_fake_coverage,
                                            min_unique_fake_draws,
                                            min_unique_fake_fraction,
                                            ...) {
  args <- c(
    list(
      expr = expr,
      celltypes = celltypes,
      lr_db = lr_db,
      candidate_fake_genes = candidate_fake_genes,
      pseudotime = pseudotime,
      lineage_weights = lineage_weights,
      n_bins = n_bins,
      n_fake_per_real = n_fake_per_real,
      n_strata = n_strata,
      caliper_k = caliper_k,
      max_match_distance = max_match_distance,
      min_pct = min_pct,
      min_n_eff = min_n_eff,
      min_valid_bins = min_valid_bins,
      min_valid_fraction = min_valid_fraction,
      min_fake_observed = min_fake_observed,
      min_fake_coverage = min_fake_coverage,
      min_fake_observed_generic = min_fake_observed,
      min_fake_coverage_generic = min_fake_coverage,
      min_fake_observed_edge = min_fake_observed,
      min_fake_coverage_edge = min_fake_coverage,
      min_unique_fake_draws = min_unique_fake_draws,
      min_unique_fake_fraction = min_unique_fake_fraction,
      seed = seed,
      complex_null_mode = complex_null_mode,
      fake_null_mode = fake_null_mode
    ),
    list(...)
  )
  do.call(ComputeFakeLRNull, args)$trajectory_calibration
}

#' Per-donor / sample temporal and identity effects (formal block rescoring).
#'
#' Requires `expr`, `celltypes`, `lr_db`, `block_labels`, `pseudotime`, and
#' `lineage_weights`. Each block is rescored with [ComputeTemporalLR()]; optional
#' block-internal FakeLR and temporal permutation supply `p_identity` and
#' `p_temporal_block`. All rows use `effect_mode = "donor_rescore"`.
#'
#' Splits cells by `block_labels` and computes trajectory-level temporal
#' effect size, identity effect, peak pseudotime, and evaluability per block.
#' Bin scores from the pooled `scores_df` are localised using block-specific
#' lineage-weighted bin occupancy; optional `merged_df` or calibration tables
#' supply trajectory-level identity p-values.
#'
#' @param scores_df data.frame from `ComputeTemporalLR` or `MergeEvidence`.
#' @param block_labels length-`n_cells` donor/sample/batch labels (see
#'   `ComputeTemporalPermutation()`).
#' @param pseudotime numeric vector of length `n_cells`.
#' @param lineage_weights numeric matrix \[n_cells × n_lineages\].
#' @param merged_df optional merged evidence table from `MergeEvidence()`.
#' @param fake_calibration optional trajectory-level FakeLR table.
#' @param temporal_calibration optional trajectory-level temporal-null table.
#' @param celltypes optional integer vector; used only when `expr` is supplied
#'   for per-block rescoring.
#' @param expr optional `dgCMatrix` for per-block `ComputeTemporalLR()` rescoring.
#' @param lr_db optional LR database when `expr` is supplied.
#' @param candidate_fake_genes integer vector of fake-gene candidates; required
#'   when `fake_null_per_block = TRUE`.
#' @param fake_null_per_block when `TRUE` and `expr`/`lr_db`/`candidate_fake_genes`
#'   are supplied, run block-stratified `ComputeFakeLRNull()` for per-donor
#'   `p_identity`. Pooled calibration is returned only as `p_identity_pooled`
#'   annotation and never enters donor identity meta-analysis.
#' @param temporal_null_per_block when `TRUE` and formal donor-rescore inputs
#'   are supplied (`expr`, `celltypes`, `pseudotime`, `lineage_weights`, `lr_db`,
#'   `block_labels`), run block-internal temporal permutation for `p_temporal_block`
#'   and bin-shuffle null medians for `temporal_effect_null_centered`.
#' @param n_perms_block temporal permutations per block (default 199).
#' @param n_null_median_draws legacy bin-shuffle draws when permutation median
#'   unavailable (default 199).
#' @param temporal_null_seed base seed for per-block temporal null (incremented per block).
#' @param permutation_stratum block temporal-permutation stratum.
#' @param min_perm_observed_block,min_perm_coverage_block block temporal-null
#'   support and coverage gates.
#' @param n_fake_per_real_block FakeLR draws per trajectory when
#'   `fake_null_per_block = TRUE` (default 199; minimum 100).
#' @param fake_null_mode passed to per-block `ComputeFakeLRNull()`.
#' @param complex_null_mode passed to per-block `ComputeFakeLRNull()`.
#' @param fake_null_seed base seed for per-block FakeLR (incremented per block).
#' @param n_strata_block,caliper_k_block,max_match_distance_block block FakeLR
#'   matching parameters.
#' @param min_fake_observed_block,min_fake_coverage_block block FakeLR quality
#'   gates applied before a donor identity p-value is evaluable.
#' @param min_unique_fake_draws_block,min_unique_fake_fraction_block block
#'   FakeLR distinct-null support gates.
#' @param min_pct_block,min_valid_fraction_block scoring gates used consistently
#'   for donor rescoring, block temporal permutation, and block FakeLR.
#' @param ... additional arguments passed to `ComputeFakeLRNull()` when
#'   `fake_null_per_block = TRUE` (e.g. `n_strata`, `min_fake_observed`).
#' @param by trajectory key columns (default `pair_idx/source/target/lineage`).
#' @param stat trajectory statistic for temporal effect (`amplitude` default).
#' @param top_k bins for `"topk_mean"`.
#' @param min_n_eff minimum Kish ESS per bin for evaluability.
#' @param min_valid_bins minimum structurally valid bins per block trajectory.
#' @param min_block_cells minimum cells required in a block.
#' @param min_bin_mass minimum lineage-weight mass in a bin to borrow pooled score.
#' @param n_bins number of pseudotime bins (inferred from `scores_df` when NULL).
#' @return data.frame with one row per trajectory × block, including
#'   `temporal_effect`, `identity_effect`, `peak_time`, `evaluable`,
#'   `effect_se` (bin-variance SE when \eqn{\ge 2} valid bins, else heuristic),
#'   `effect_mode` (always `donor_rescore`), optional `temporal_null_median`,
#'   `temporal_effect_null_centered`, `p_temporal_block`, and audit columns.
#'
#' Raw `temporal_effect` without null centering is descriptive only; use
#' `temporal_effect_null_centered` and [MetaAnalyzeTemporalEffects()] for inference.
#' @export
ComputeTemporalEffectsBySample <- function(scores_df,
                                           block_labels,
                                           pseudotime,
                                           lineage_weights,
                                           merged_df = NULL,
                                           fake_calibration = NULL,
                                           temporal_calibration = NULL,
                                           celltypes = NULL,
                                           expr = NULL,
                                           lr_db = NULL,
                                           candidate_fake_genes = NULL,
                                           fake_null_per_block = FALSE,
                                           temporal_null_per_block = NULL,
                                           n_perms_block = 199L,
                                           n_null_median_draws = 199L,
                                           temporal_null_seed = 23L,
                                           permutation_stratum = c(
                                             "block_celltype_lineage", "block_celltype"
                                           ),
                                           min_perm_observed_block = 50L,
                                           min_perm_coverage_block = 0.5,
                                           n_fake_per_real_block = 199L,
                                           fake_null_mode = "context_generic",
                                           complex_null_mode = c("skip", "cardinality_matched"),
                                           fake_null_seed = 17L,
                                           n_strata_block = 3L,
                                           caliper_k_block = 8L,
                                           max_match_distance_block = 2,
                                           min_fake_observed_block = 100L,
                                           min_fake_coverage_block = 0.8,
                                           min_unique_fake_draws_block = 10L,
                                           min_unique_fake_fraction_block = 0.05,
                                           min_pct_block = 0.05,
                                           min_valid_fraction_block = 0.1,
                                           ...,
                                           by = c("pair_idx", "source", "target", "lineage"),
                                           stat = "amplitude",
                                           top_k = 3L,
                                           min_n_eff = 5.0,
                                           min_valid_bins = 2L,
                                           min_block_cells = 10L,
                                           min_bin_mass = 1.0,
                                           n_bins = NULL) {
  stopifnot(is.data.frame(scores_df))
  by <- as.character(by)
  missing <- setdiff(by, colnames(scores_df))
  if (length(missing) > 0L) {
    stop("`scores_df` must contain trajectory key columns: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  if (!"bin" %in% colnames(scores_df) || !"score" %in% colnames(scores_df)) {
    stop("`scores_df` must contain `bin` and `score` columns", call. = FALSE)
  }
  stat_info <- .normalize_trajectory_stat(stat)
  stat <- stat_info$stat
  top_k <- .validate_positive_int(top_k, "top_k", min = 1L)
  min_valid_bins <- .validate_positive_int(min_valid_bins, "min_valid_bins", min = 1L)
  min_block_cells <- .validate_positive_int(min_block_cells, "min_block_cells", min = 1L)
  min_n_eff <- as.numeric(min_n_eff)
  min_bin_mass <- as.numeric(min_bin_mass)
  min_pct_block <- .validate_unit_interval(min_pct_block, "min_pct_block")
  min_valid_fraction_block <- .validate_unit_interval(
    min_valid_fraction_block, "min_valid_fraction_block"
  )
  if (!is.finite(min_n_eff) || min_n_eff < 0 ||
      !is.finite(min_bin_mass) || min_bin_mass < 0) {
    stop("`min_n_eff` and `min_bin_mass` must be finite non-negative numerics",
         call. = FALSE)
  }

  n_cells <- length(pseudotime)
  blocks <- .prepare_block_labels(block_labels, n_cells)
  block_levels <- attr(blocks, "levels")
  lw <- .assert_lineage_weights(lineage_weights, n_cells)
  pt <- .assert_pseudotime_1d(pseudotime, n_cells)
  n_bins <- .infer_n_bins(scores_df, n_bins)

  if (!.donor_rescore_formal_mode(
    expr, celltypes, pseudotime, lineage_weights, lr_db, block_labels
  )) {
    stop(
      "`ComputeTemporalEffectsBySample()` requires `expr`, `celltypes`, `lr_db`, ",
      "`block_labels`, `pseudotime`, and `lineage_weights` for formal per-donor ",
      "inference.",
      call. = FALSE
    )
  }
  formal_donor_mode <- TRUE
  rescore_blocks <- TRUE
  effect_mode <- "donor_rescore"
  ct <- .assert_celltypes(celltypes, n_cells)
  ct_levels_scoring <- attr(ct, "levels")
  pooled_contexts_df <- unique(scores_df[, c("pair_idx", "source", "target"), drop = FALSE])
  pooled_trajectories_df <- unique(scores_df[, by, drop = FALSE])
  if (is.null(temporal_null_per_block)) {
    temporal_null_per_block <- TRUE
  } else {
    temporal_null_per_block <- isTRUE(temporal_null_per_block)
  }
  permutation_stratum <- match.arg(permutation_stratum)
  min_perm_observed_block <- .validate_non_negative_int(
    min_perm_observed_block, "min_perm_observed_block"
  )
  min_perm_coverage_block <- .validate_unit_interval(
    min_perm_coverage_block, "min_perm_coverage_block"
  )
  temporal_null_seed <- .validate_seed(temporal_null_seed)
  if (temporal_null_per_block) {
    if (!formal_donor_mode) {
      stop(
        "`temporal_null_per_block = TRUE` requires formal donor-rescore inputs ",
        "(`expr`, `celltypes`, `pseudotime`, `lineage_weights`, `lr_db`, ",
        "`block_labels`)",
        call. = FALSE
      )
    }
    n_perms_block <- .validate_positive_int(n_perms_block, "n_perms_block", min = 1L)
    n_null_median_draws <- .validate_positive_int(
      n_null_median_draws, "n_null_median_draws", min = 1L
    )
  }
  fake_null_per_block <- isTRUE(fake_null_per_block)
  if (fake_null_per_block) {
    if (n_fake_per_real_block < 100L) {
      stop("`n_fake_per_real_block` must be >= 100 for formal FakeLR per block",
           call. = FALSE)
    }
    if (is.null(expr) || is.null(lr_db) || is.null(candidate_fake_genes)) {
      stop("`fake_null_per_block = TRUE` requires `expr`, `lr_db`, and ",
           "`candidate_fake_genes`", call. = FALSE)
    }
    n_fake_per_real_block <- .validate_positive_int(
      n_fake_per_real_block, "n_fake_per_real_block", min = 100L
    )
    fake_null_mode <- .normalize_fake_null_mode(fake_null_mode)
    complex_null_mode <- match.arg(complex_null_mode)
    fake_null_seed <- .validate_seed(fake_null_seed)
    n_strata_block <- .validate_positive_int(n_strata_block, "n_strata_block", min = 1L)
    caliper_k_block <- .validate_positive_int(caliper_k_block, "caliper_k_block", min = 1L)
    max_match_distance_block <- .validate_positive_num(
      max_match_distance_block, "max_match_distance_block"
    )
    min_fake_observed_block <- .validate_non_negative_int(
      min_fake_observed_block, "min_fake_observed_block"
    )
    min_fake_coverage_block <- .validate_unit_interval(
      min_fake_coverage_block, "min_fake_coverage_block"
    )
    min_unique_fake_draws_block <- .validate_non_negative_int(
      min_unique_fake_draws_block, "min_unique_fake_draws_block"
    )
    min_unique_fake_fraction_block <- .validate_unit_interval(
      min_unique_fake_fraction_block, "min_unique_fake_fraction_block"
    )
  }

  cal <- if (fake_null_per_block) {
    NULL
  } else {
    .merge_trajectory_calibration(
      by, merged_df, fake_calibration, temporal_calibration
    )
  }
  scores_df <- scores_df[, c(by, "bin", "lineage", "score",
                             intersect(c("n_eff_cells", "lineage_entropy"),
                                       colnames(scores_df))), drop = FALSE]
  scores_df$bin <- .as_integerish(scores_df$bin, "scores_df$bin")
  scores_df$lineage <- .as_integerish(scores_df$lineage, "scores_df$lineage")
  scores_df$.traj_key <- .trajectory_key_col(scores_df, by)
  score_split <- split(scores_df, scores_df$.traj_key)

  rows <- list()
  row_i <- 0L
  for (b in seq_along(block_levels)) {
    block_cells <- which(blocks == b)
    n_block <- length(block_cells)
    if (n_block < min_block_cells) {
      next
    }
    block_scores <- if (rescore_blocks) {
      block_expr <- expr[, block_cells, drop = FALSE]
      block_ct <- celltypes[block_cells]
      block_pt <- pt[block_cells]
      block_lw <- lw[block_cells, , drop = FALSE]
      ComputeTemporalLR(
        block_expr, block_ct, lr_db, block_pt, block_lw,
        n_bins = n_bins, min_pct = min_pct_block,
        min_n_eff = max(min_n_eff, 0.1), return_complete_grid = TRUE,
        requested_contexts = pooled_contexts_df
      )
    } else {
      NULL
    }
    if (rescore_blocks && nrow(block_scores) > 0L) {
      block_scores$.traj_key <- .trajectory_key_col(block_scores, by)
      block_split <- split(block_scores, block_scores$.traj_key)
    } else {
      block_split <- NULL
    }
    mass_mat <- .block_lineage_bin_mass(block_cells, pt, lw, n_bins)
    block_fake_cal <- NULL
    if (fake_null_per_block) {
      block_fake_cal <- .compute_block_fake_calibration(
          expr[, block_cells, drop = FALSE],
          celltypes[block_cells],
          lr_db,
          candidate_fake_genes = as.integer(candidate_fake_genes),
          pseudotime = pt[block_cells],
          lineage_weights = lw[block_cells, , drop = FALSE],
          n_bins = n_bins,
          n_fake_per_real = n_fake_per_real_block,
          complex_null_mode = complex_null_mode,
          fake_null_mode = fake_null_mode,
          seed = fake_null_seed + b - 1L,
          n_strata = n_strata_block,
          caliper_k = caliper_k_block,
          max_match_distance = max_match_distance_block,
          min_pct = min_pct_block,
          min_n_eff = max(min_n_eff, 0.1),
          min_valid_bins = max(min_valid_bins, 2L),
          min_valid_fraction = min_valid_fraction_block,
          min_fake_observed = min_fake_observed_block,
          min_fake_coverage = min_fake_coverage_block,
          min_unique_fake_draws = min_unique_fake_draws_block,
          min_unique_fake_fraction = min_unique_fake_fraction_block,
          ...
      )
    }
    block_temporal_cal <- NULL
    if (temporal_null_per_block) {
      block_temporal_cal <- .compute_block_temporal_permutation_table(
          expr = expr[, block_cells, drop = FALSE],
          celltypes = celltypes[block_cells],
          lr_db = lr_db,
          pseudotime = pt[block_cells],
          lineage_weights = lw[block_cells, , drop = FALSE],
          n_bins = n_bins,
          stat = stat,
          top_k = top_k,
          n_perms = n_perms_block,
          seed = temporal_null_seed + b - 1L,
          min_pct = min_pct_block,
          min_n_eff = max(min_n_eff, 0.1),
          min_valid_bins = max(min_valid_bins, 2L),
          min_valid_fraction = min_valid_fraction_block,
          permutation_stratum = permutation_stratum,
          min_perm_observed = min_perm_observed_block,
          min_perm_coverage = min_perm_coverage_block,
          requested_trajectories = pooled_trajectories_df
      )
    }

    for (traj_key in names(score_split)) {
      traj_df <- score_split[[traj_key]]
      lineage_id <- traj_df$lineage[1L]
      if (rescore_blocks && !is.null(block_split) && traj_key %in% names(block_split)) {
        bdf <- block_split[[traj_key]]
        bin_info <- .donor_bins_from_block_scores(bdf, n_bins)
        y <- bin_info$y
        n_valid <- bin_info$n_valid
        n_eff <- bin_info$n_eff
      } else if (rescore_blocks) {
        y <- rep(NA_real_, n_bins)
        n_valid <- 0L
        n_eff <- NA_real_
      } else {
        y <- rep(NA_real_, n_bins)
        for (i in seq_len(nrow(traj_df))) {
          bin_i <- traj_df$bin[i]
          if (bin_i >= 1L && bin_i <= n_bins &&
              mass_mat[lineage_id, bin_i] >= min_bin_mass) {
            y[bin_i] <- as.numeric(traj_df$score[i])
          }
        }
        n_valid <- sum(is.finite(y))
        n_eff <- if ("n_eff_cells" %in% colnames(traj_df)) {
          stats::median(as.numeric(traj_df$n_eff_cells), na.rm = TRUE) *
            (n_block / n_cells)
        } else {
          NA_real_
        }
      }
      temporal_effect <- .trajectory_stat_from_bins(y, stat, top_k)
      peak_bin <- if (any(is.finite(y))) which.max(replace(y, !is.finite(y), -Inf)) else NA_integer_
      peak_time <- if (is.na(peak_bin)) {
        NA_real_
      } else {
        .block_peak_time(block_cells, pt, lw, lineage_id, peak_bin, n_bins)
      }
      evaluable <- is.finite(temporal_effect) &&
        n_valid >= min_valid_bins &&
        (is.na(n_eff) || n_eff >= min_n_eff)

      p_identity <- NA_real_
      p_identity_pooled <- NA_real_
      if (fake_null_per_block) {
        p_identity <- .extract_calibration_p_identity(
          block_fake_cal, traj_key, by
        )
      } else if (!is.null(cal)) {
        p_identity_pooled <- .extract_calibration_p_identity(
          cal, traj_key, by
        )
      }
      identity_effect <- NA_real_
      identity_effect_null_centered <- NA_real_
      fake_null_median <- NA_real_
      trajectory_stat_cal <- NA_real_
      cal_for_identity <- if (fake_null_per_block) block_fake_cal else NULL
      if (!is.null(cal_for_identity)) {
        fake_null_median <- .extract_calibration_fake_null_median(
          cal_for_identity, traj_key, by
        )
        trajectory_stat_cal <- .extract_calibration_trajectory_stat(
          cal_for_identity, traj_key, by
        )
      }
      if (is.finite(trajectory_stat_cal) && is.finite(fake_null_median)) {
        identity_effect_null_centered <- trajectory_stat_cal - fake_null_median
        identity_effect <- identity_effect_null_centered
      }
      evaluable_identity <- fake_null_per_block && is.finite(p_identity) &&
        is.finite(identity_effect_null_centered)

      effect_se <- if (is.finite(temporal_effect)) {
        .sample_effect_se(temporal_effect, n_eff, n_valid, bin_scores = y)
      } else {
        NA_real_
      }

      temporal_null_median <- NA_real_
      temporal_null_source <- NA_character_
      p_temporal_block <- NA_real_
      if (temporal_null_per_block && !is.null(block_temporal_cal)) {
        cal_idx <- match(traj_key, block_temporal_cal$.traj_key)
        if (!is.na(cal_idx)) {
          if ("temporal_null_median" %in% colnames(block_temporal_cal) &&
              is.finite(block_temporal_cal$temporal_null_median[cal_idx])) {
            temporal_null_median <- block_temporal_cal$temporal_null_median[cal_idx]
            temporal_null_source <- "permutation"
          } else if ("temporal_stat" %in% colnames(block_temporal_cal) &&
                     is.finite(block_temporal_cal$temporal_stat[cal_idx])) {
            temporal_null_median <- block_temporal_cal$temporal_stat[cal_idx]
            temporal_null_source <- "permutation_observed_only"
          }
          if ("p_temporal" %in% colnames(block_temporal_cal)) {
            p_temporal_block <- block_temporal_cal$p_temporal[cal_idx]
          }
        }
      } else if (formal_donor_mode && any(is.finite(y))) {
        traj_idx <- match(traj_key, names(score_split))
        temporal_null_median <- .block_trajectory_null_median(
          y, stat, top_k,
          n_draws = n_null_median_draws,
          seed = temporal_null_seed + b * 1009L + traj_idx
        )
        temporal_null_source <- "bin_shuffle_exploratory"
      }
      temporal_effect_null_centered <- if (is.finite(temporal_effect) &&
                                          is.finite(temporal_null_median)) {
        temporal_effect - temporal_null_median
      } else {
        NA_real_
      }

      row_i <- row_i + 1L
      key_row <- traj_df[1L, by, drop = FALSE]
      rows[[row_i]] <- data.frame(
        key_row,
        block_id = b,
        block_label = block_levels[[b]],
        effect_mode = effect_mode,
        temporal_effect = temporal_effect,
        temporal_null_median = temporal_null_median,
        temporal_null_source = temporal_null_source,
        temporal_effect_null_centered = temporal_effect_null_centered,
        p_temporal_block = p_temporal_block,
        identity_effect = identity_effect,
        identity_effect_null_centered = identity_effect_null_centered,
        fake_null_median = fake_null_median,
        peak_time = peak_time,
        evaluable = evaluable,
        n_cells = n_block,
        n_eff_cells = n_eff,
        n_valid_bins = n_valid,
        effect_se = effect_se,
        p_identity = p_identity,
        p_identity_pooled = p_identity_pooled,
        evaluable_identity = evaluable_identity,
        identity_source = if (fake_null_per_block) {
          "block_fake_null"
        } else if (!is.null(cal)) {
          "pooled_annotation"
        } else {
          NA_character_
        },
        trajectory_stat = stat,
        stringsAsFactors = FALSE
      )
    }
  }

  if (row_i < 1L) {
    stop("No blocks passed `min_block_cells`; check `block_labels` and thresholds",
         call. = FALSE)
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' Donor-level meta-analysis of sample temporal effects.
#'
#' Combines per-block effects from [ComputeTemporalEffectsBySample()]. Prefer
#' `method = "stouffer"` / `"fisher"` on `p_temporal_block` with
#' `temporal_effect_null_centered`.
#'
#' @param sample_effects data.frame from `ComputeTemporalEffectsBySample()`.
#' @param effect_col column holding sample-level effect sizes (default
#'   `"temporal_effect_null_centered"` when present, else `"temporal_effect"`).
#' @param se_col column holding sample-level standard errors (used for
#'   `"fixed"` / `"random"` methods).
#' @param p_col column holding block-level p-values for `"stouffer"` / `"fisher"`
#'   (default `"p_temporal_block"`).
#' @param by trajectory key columns.
#' @param method `"fisher"` (default) or one-sided `"stouffer"` p-value
#'   combination across donors. `"fixed"` and `"random"` are exploratory
#'   effect summaries unless validated donor-level sampling SEs are supplied.
#' @param include_lodo include leave-one-donor-out stability columns (fixed meta).
#' @param alpha significance threshold for LOO stability audit.
#' @return data.frame with `meta_effect`, `meta_p`, `q_meta` (BH on the full
#'   trajectory family), `heterogeneity_I2`, `sign_consistency`, `meta_method`,
#'   `tau2`, `meta_p_random`, `meta_se_random`, and optional `lodo_stability`.
#'   Requires `n_samples >= 2` evaluable donors for finite `meta_p`; random
#'   effects warn when `n_samples < 3`.
#' @export
MetaAnalyzeTemporalEffects <- function(sample_effects,
                                       effect_col = NULL,
                                       se_col = "effect_se",
                                       p_col = "p_temporal_block",
                                       by = c("pair_idx", "source", "target", "lineage"),
                                       method = c("fisher", "stouffer", "fixed", "random"),
                                       include_lodo = TRUE,
                                       alpha = 0.05) {
  if (is.null(effect_col)) {
    effect_col <- if ("temporal_effect_null_centered" %in% colnames(sample_effects)) {
      "temporal_effect_null_centered"
    } else {
      "temporal_effect"
    }
  }
  .meta_analyze_sample_effects(
    sample_effects = sample_effects,
    effect_col = effect_col,
    se_col = se_col,
    p_col = p_col,
    by = by,
    method = match.arg(method),
    p_sided = "one",
    include_lodo = include_lodo,
    alpha = alpha
  )
}

#' Donor-level meta-analysis of sample identity effects.
#'
#' Parallel to [MetaAnalyzeTemporalEffects()] for identity columns. Formal input
#' is restricted to evaluable rows with `identity_source = "block_fake_null"`;
#' pooled identity calibration is never treated as replicated donor evidence.
#'
#' @param sample_effects data.frame from `ComputeTemporalEffectsBySample()`.
#' @param effect_col column holding identity effect sizes (default
#'   `"identity_effect_null_centered"` when present, else `"identity_effect"`).
#' @param se_col column holding sample-level standard errors.
#' @param p_col block-level identity p-values (default `"p_identity"`).
#' @param by trajectory key columns.
#' @param method one-sided `"stouffer"` / `"fisher"` recommended for formal inference;
#'   see [MetaAnalyzeTemporalEffects()].
#' @param include_lodo include leave-one-donor-out stability columns.
#' @param alpha significance threshold for LOO stability audit.
#' @return data.frame with the same meta columns as
#'   [MetaAnalyzeTemporalEffects()], including `q_meta`.
#' @export
MetaAnalyzeIdentityEffects <- function(sample_effects,
                                       effect_col = NULL,
                                       se_col = "effect_se",
                                       p_col = "p_identity",
                                       by = c("pair_idx", "source", "target", "lineage"),
                                       method = c("stouffer", "fisher", "fixed", "random"),
                                       include_lodo = TRUE,
                                       alpha = 0.05) {
  if (!"identity_source" %in% colnames(sample_effects)) {
    stop(
      "Formal identity meta-analysis requires `identity_source = \"block_fake_null\"`",
      call. = FALSE
    )
  }
  if (!"evaluable_identity" %in% colnames(sample_effects)) {
    sample_effects$evaluable_identity <- is.finite(sample_effects[[p_col]])
  }
  sample_effects <- sample_effects[
    sample_effects$identity_source == "block_fake_null" &
      sample_effects$evaluable_identity %in% TRUE,
    ,
    drop = FALSE
  ]
  if (nrow(sample_effects) < 1L) {
    stop("No evaluable block-specific FakeLR rows for identity meta-analysis",
         call. = FALSE)
  }
  sample_effects$evaluable <- sample_effects$evaluable_identity
  defaults <- .resolve_identity_meta_defaults(
    sample_effects, effect_col, method[1L], p_col
  )
  .meta_analyze_sample_effects(
    sample_effects = sample_effects,
    effect_col = defaults$effect_col,
    se_col = se_col,
    p_col = defaults$p_col,
    by = by,
    method = match.arg(method, c("stouffer", "fisher", "fixed", "random")),
    p_sided = "one",
    include_lodo = include_lodo,
    alpha = alpha
  )
}

#' Meta-analyse per-donor lag effects.
#'
#' Parallel to [MetaAnalyzeTemporalEffects()] for lag columns from
#' [SummarizeLagByDonor()]. Prefer `method = "stouffer"` / `"fisher"` on
#' `p_lag` when `compute_lag_null = TRUE` was used upstream.
#'
#' @param sample_effects data.frame from [SummarizeLagByDonor()].
#' @param effect_col lag effect column (default `"lag_null_centered"` when
#'   present, else `"lag"`).
#' @param se_col column holding sample-level standard errors (optional; many lag
#'   rows omit SE — fixed/random meta may be NA).
#' @param p_col block-level lag p-values (default `"p_lag"`).
#' @param by trajectory key columns.
#' @param method signed two-sided `"stouffer"` / `"fisher"` recommended for
#'   lag inference.
#' @param include_lodo include leave-one-donor-out stability columns.
#' @param alpha significance threshold for LOO stability audit.
#' @return data.frame with meta columns including `q_meta`.
#' @export
MetaAnalyzeLagEffects <- function(sample_effects,
                                  effect_col = NULL,
                                  se_col = "effect_se",
                                  p_col = "p_lag",
                                  by = c("pair_idx", "source", "target", "lineage"),
                                  method = c("stouffer", "fisher", "fixed", "random"),
                                  include_lodo = TRUE,
                                  alpha = 0.05) {
  defaults <- .resolve_lag_meta_defaults(
    sample_effects, effect_col, method[1L], p_col
  )
  if (!defaults$effect_col %in% colnames(sample_effects)) {
    stop(
      "`sample_effects` missing lag effect column `",
      defaults$effect_col, "`",
      call. = FALSE
    )
  }
  se_use <- se_col
  if (!se_use %in% colnames(sample_effects)) {
    sample_effects[[se_use]] <- NA_real_
  }
  if ("evaluable" %in% colnames(sample_effects)) {
    sample_effects <- .filter_formal_evaluable_rows(sample_effects)
  } else if ("p_lag" %in% colnames(sample_effects)) {
    sample_effects <- sample_effects[is.finite(sample_effects$p_lag), , drop = FALSE]
  }
  .meta_analyze_sample_effects(
    sample_effects = sample_effects,
    effect_col = defaults$effect_col,
    se_col = se_use,
    p_col = defaults$p_col,
    by = by,
    method = match.arg(method, c("stouffer", "fisher", "fixed", "random")),
    p_sided = "two",
    include_lodo = include_lodo,
    alpha = alpha
  )
}

.meta_analyze_sample_effects <- function(sample_effects,
                                          effect_col,
                                          se_col,
                                          p_col,
                                          by,
                                          method,
                                          p_sided,
                                          include_lodo,
                                          alpha) {
  stopifnot(is.data.frame(sample_effects))
  by <- as.character(by)
  method <- match.arg(method, c("fixed", "random", "stouffer", "fisher"))
  effect_col <- .resolve_meta_effect_col(sample_effects, effect_col)
  if (!all(c(by, effect_col, se_col) %in% colnames(sample_effects))) {
    stop("`sample_effects` missing required columns", call. = FALSE)
  }
  if (method %in% c("stouffer", "fisher") && !p_col %in% colnames(sample_effects)) {
    stop(
      sprintf("`method = \"%s\"` requires `%s` in `sample_effects`", method, p_col),
      call. = FALSE
    )
  }
  alpha <- as.numeric(alpha)
  if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be a finite numeric in (0, 1)", call. = FALSE)
  }

  sample_effects <- .filter_formal_evaluable_rows(sample_effects)
  sample_effects$.traj_key <- .trajectory_key_col(sample_effects, by)
  groups <- split(sample_effects, sample_effects$.traj_key)
  .warn_random_effects_donor_count(groups, effect_col, se_col, method)

  meta_rows <- lapply(names(groups), function(key) {
    g <- groups[[key]]
    y <- g[[effect_col]]
    se <- g[[se_col]]
    fe <- .fixed_effect_meta(y, se)
    tau2 <- if (fe$n_samples >= 2L) {
      .dersimonian_laird_tau2(y, se, fe$meta_effect)
    } else {
      NA_real_
    }
    re <- if (method == "random") {
      .random_effect_meta(y, se, tau2 = tau2)
    } else {
      NULL
    }
    het <- .heterogeneity_stats(y, se, fe$meta_effect)
    sc <- .sign_consistency(y)
    lodo <- if (include_lodo && fe$n_samples >= 2L &&
                method %in% c("fixed", "random")) {
      LeaveOneSampleOut(
        g, effect_col = effect_col, se_col = se_col, by = by, alpha = alpha
      )
    } else {
      NULL
    }
    lodo_stability <- NA_real_
    lodo_min_sign <- NA_real_
    if (!is.null(lodo) && nrow(lodo) > 0L) {
      lodo_stability <- mean(lodo$sign_stable & lodo$meta_p <= alpha, na.rm = TRUE)
      lodo_min_sign <- if (any(is.finite(lodo$sign_consistency))) {
        min(lodo$sign_consistency, na.rm = TRUE)
      } else {
        NA_real_
      }
    }
    meta_p_random <- if (method == "random" && !is.null(re)) re$meta_p else NA_real_
    meta_se_random <- if (method == "random" && !is.null(re)) re$meta_se else NA_real_
    meta_effect_out <- if (method == "random" && !is.null(re)) re$meta_effect else fe$meta_effect
    meta_se_out <- if (method == "random" && !is.null(re)) re$meta_se else fe$meta_se
    meta_p_out <- if (method == "random" && !is.null(re)) re$meta_p else fe$meta_p
    meta_z_out <- if (method == "random" && !is.null(re)) re$meta_z else fe$meta_z
    if (method == "stouffer") {
      meta_p_out <- .stouffer_combine_p(
        g[[p_col]],
        sided = p_sided,
        signs = if (p_sided == "two") y else NULL
      )
      meta_z_out <- NA_real_
      meta_se_out <- NA_real_
    } else if (method == "fisher") {
      meta_p_out <- .fisher_combine_p(g[[p_col]])
      meta_z_out <- NA_real_
      meta_se_out <- NA_real_
    }
    n_p_ok <- sum(is.finite(g[[p_col]]) & g[[p_col]] > 0 & g[[p_col]] <= 1)
    n_meta_ok <- if (method %in% c("stouffer", "fisher")) n_p_ok else fe$n_samples
    if (n_meta_ok < 2L) {
      meta_p_out <- NA_real_
      meta_z_out <- NA_real_
      meta_p_random <- NA_real_
    }
    cbind(
      g[1L, by, drop = FALSE],
      data.frame(
        meta_effect = meta_effect_out,
        meta_se = meta_se_out,
        meta_z = meta_z_out,
        meta_p = meta_p_out,
        meta_method = method,
        meta_inference_status = if (method %in% c("fisher", "stouffer")) {
          "formal_p_combination"
        } else {
          "exploratory_effect_summary"
        },
        meta_effect_fixed = fe$meta_effect,
        meta_se_fixed = fe$meta_se,
        meta_p_fixed = if (fe$n_samples >= 2L) fe$meta_p else NA_real_,
        tau2 = if (method == "random" && !is.null(re)) re$tau2 else tau2,
        meta_p_random = if (fe$n_samples >= 2L) meta_p_random else NA_real_,
        meta_se_random = meta_se_random,
        heterogeneity_Q = het$Q,
        heterogeneity_I2 = het$I2,
        heterogeneity_df = het$df,
        sign_consistency = sc,
        n_samples = if (method %in% c("stouffer", "fisher")) n_p_ok else fe$n_samples,
        n_blocks = length(unique(g$block_id)),
        lodo_stability = lodo_stability,
        lodo_min_sign_consistency = lodo_min_sign,
        stringsAsFactors = FALSE
      )
    )
  })
  out <- do.call(rbind, meta_rows)
  rownames(out) <- NULL
  out$q_meta <- stats::p.adjust(out$meta_p, method = "BH")
  attr(out, "hypothesis_family") <- "donor_meta_trajectories"
  out
}

#' Leave-one-donor-out meta-analysis stability.
#'
#' For each trajectory, recomputes the fixed-effect meta statistic after
#' dropping one block at a time.
#'
#' @param sample_effects data.frame from `ComputeTemporalEffectsBySample()`.
#' @param effect_col,se_col,by see [MetaAnalyzeTemporalEffects()].
#' @param alpha significance threshold used for stability flags.
#' @return long data.frame with one row per trajectory × left-out block.
#' @export
LeaveOneSampleOut <- function(sample_effects,
                              effect_col = "temporal_effect",
                              se_col = "effect_se",
                              by = c("pair_idx", "source", "target", "lineage"),
                              alpha = 0.05) {
  stopifnot(is.data.frame(sample_effects))
  by <- as.character(by)
  if (!all(c(by, effect_col, se_col, "block_id") %in% colnames(sample_effects))) {
    stop("`sample_effects` must contain effect, se, and `block_id` columns",
         call. = FALSE)
  }
  alpha <- as.numeric(alpha)
  sample_effects <- .filter_formal_evaluable_rows(sample_effects)
  sample_effects$.traj_key <- .trajectory_key_col(sample_effects, by)
  groups <- split(sample_effects, sample_effects$.traj_key)

  rows <- lapply(names(groups), function(key) {
    g <- groups[[key]]
    full <- .fixed_effect_meta(g[[effect_col]], g[[se_col]])
    blocks <- unique(g$block_id)
    lapply(blocks, function(left_out) {
      h <- g[g$block_id != left_out, , drop = FALSE]
      fe <- .fixed_effect_meta(h[[effect_col]], h[[se_col]])
      het <- .heterogeneity_stats(h[[effect_col]], h[[se_col]], fe$meta_effect)
      sc <- .sign_consistency(h[[effect_col]])
      sign_stable <- is.finite(full$meta_effect) && is.finite(fe$meta_effect) &&
        sign(full$meta_effect) == sign(fe$meta_effect)
      cbind(
        g[1L, by, drop = FALSE],
        data.frame(
          left_out_block_id = left_out,
          left_out_block_label = g$block_label[match(left_out, g$block_id)],
          meta_effect = fe$meta_effect,
          meta_p = fe$meta_p,
          sign_consistency = sc,
          heterogeneity_I2 = het$I2,
          sign_stable = sign_stable,
          meta_p_significant = is.finite(fe$meta_p) && fe$meta_p <= alpha,
          stringsAsFactors = FALSE
        )
      )
    })
  })
  out <- do.call(rbind, unlist(rows, recursive = FALSE))
  rownames(out) <- NULL
  out
}

#' Bootstrap donor-level meta-analysis over blocks.
#'
#' Resamples blocks with replacement and recomputes fixed-effect meta statistics.
#'
#' @param sample_effects data.frame from `ComputeTemporalEffectsBySample()`.
#' @param n_boot number of bootstrap replicates.
#' @param seed random seed.
#' @param effect_col,se_col,by see [MetaAnalyzeTemporalEffects()].
#' @return list with `summary` (per-trajectory bootstrap quantiles) and
#'   `replicates` (long replicate table when `n_boot <= 200`).
#' @export
BootstrapSamples <- function(sample_effects,
                             n_boot = 999L,
                             seed = 12345L,
                             effect_col = "temporal_effect",
                             se_col = "effect_se",
                             by = c("pair_idx", "source", "target", "lineage")) {
  stopifnot(is.data.frame(sample_effects))
  by <- as.character(by)
  n_boot <- .validate_positive_int(n_boot, "n_boot", min = 1L)
  seed <- .validate_seed(seed)
  if (!all(c(by, effect_col, se_col, "block_id") %in% colnames(sample_effects))) {
    stop("`sample_effects` missing required columns", call. = FALSE)
  }

  sample_effects <- .filter_formal_evaluable_rows(sample_effects)
  sample_effects$.traj_key <- .trajectory_key_col(sample_effects, by)
  groups <- split(sample_effects, sample_effects$.traj_key)

  set.seed(seed)
  summary_rows <- list()
  replicate_rows <- if (n_boot <= 200L) list() else NULL

  for (key in names(groups)) {
    g <- groups[[key]]
    blocks <- unique(g$block_id)
    if (length(blocks) < 2L) {
      next
    }
    boot_effects <- numeric(n_boot)
    boot_p <- numeric(n_boot)
    for (r in seq_len(n_boot)) {
      draw <- sample(blocks, length(blocks), replace = TRUE)
      h <- do.call(
        rbind,
        lapply(draw, function(b) g[g$block_id == b, , drop = FALSE])
      )
      fe <- .fixed_effect_meta(h[[effect_col]], h[[se_col]])
      boot_effects[r] <- fe$meta_effect
      boot_p[r] <- fe$meta_p
      if (!is.null(replicate_rows)) {
        replicate_rows[[length(replicate_rows) + 1L]] <- data.frame(
          g[1L, by, drop = FALSE],
          replicate = r,
          meta_effect = fe$meta_effect,
          meta_p = fe$meta_p,
          stringsAsFactors = FALSE
        )
      }
    }
    summary_rows[[key]] <- cbind(
      g[1L, by, drop = FALSE],
      data.frame(
        meta_effect_mean = mean(boot_effects, na.rm = TRUE),
        meta_effect_q025 = stats::quantile(boot_effects, 0.025, na.rm = TRUE),
        meta_effect_q975 = stats::quantile(boot_effects, 0.975, na.rm = TRUE),
        meta_p_median = stats::median(boot_p, na.rm = TRUE),
        n_boot = n_boot,
        stringsAsFactors = FALSE
      )
    )
  }

  list(
    summary = if (length(summary_rows)) do.call(rbind, summary_rows) else NULL,
    replicates = if (length(replicate_rows)) do.call(rbind, replicate_rows) else NULL
  )
}
