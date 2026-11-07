# Fork-local sibling-lineage evidence.

#' Map each lineage to its nearest sibling lineage or lineages.
#'
#' Nearest siblings are the other root-to-terminal paths sharing the longest
#' prefix with a focal path in the same rooted component. Ties are retained as
#' separate rows; unrelated roots do not define a sibling contrast.
#'
#' @param lineage_paths list of positive integer root-to-terminal cluster paths.
#' @return Data frame with focal and sibling lineage indices, fork cluster, and
#'   shared-prefix length.
#' @export
BuildSiblingContrastMap <- function(lineage_paths) {
  if (!is.list(lineage_paths) || length(lineage_paths) < 2L) {
    stop("`lineage_paths` must contain at least two paths", call. = FALSE)
  }
  paths <- lapply(seq_along(lineage_paths), function(i) {
    path <- .as_integerish(lineage_paths[[i]], sprintf("lineage_paths[[%d]]", i))
    if (!length(path) || any(path < 1L)) {
      stop("lineage paths must be non-empty positive integer vectors",
           call. = FALSE)
    }
    path
  })
  roots <- vapply(paths, `[[`, integer(1L), 1L)
  keys <- vapply(paths, paste, collapse = ">", character(1L))
  if (anyDuplicated(keys)) {
    stop("duplicate lineage path", call. = FALSE)
  }

  rows <- lapply(seq_along(paths), function(focal) {
    siblings <- which(roots == roots[focal] & seq_along(paths) != focal)
    prefix <- vapply(siblings, function(sibling) {
      .common_prefix_length(paths[[focal]], paths[[sibling]])
    }, integer(1L))
    keep <- siblings[prefix == max(c(0L, prefix))]
    kept_prefix <- prefix[prefix == max(c(0L, prefix))]
    data.frame(
      focal_lineage = rep.int(as.integer(focal), length(keep)),
      sibling_lineage = as.integer(keep),
      fork_cluster = as.integer(paths[[focal]][kept_prefix]),
      shared_prefix_length = as.integer(kept_prefix),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

.common_prefix_length <- function(a, b) {
  n <- min(length(a), length(b))
  mismatch <- which(a[seq_len(n)] != b[seq_len(n)])
  if (length(mismatch) == 0L) n else mismatch[[1L]] - 1L
}

#' Compute fork-local contrasts against nearest sibling lineages.
#'
#' Two estimators are available.
#'
#' Shared-window estimator (used when `fork_bins` is supplied or the input
#' declares a common temporal coordinate): every compared branch is scored on
#' one shared grid, one common fork bin `b*` anchors the pre/post windows of
#' the focal lineage and all of its sisters, and the effect is
#' `theta = mean_post(D) - mean_pre(D)` with `D(t) = focal(t) - mean_sister(t)`
#' evaluated only at bins where the focal and at least one sister both have
#' support.  All branches enter `D(t)` with identical (uniform) time weights
#' and sisters enter with equal weights, so a common smooth temporal trend
#' cancels exactly.  `fork_bins` gives the topology-derived anchor; without it
#' the focal lineage's own entropy `b_star` is used and recorded as
#' `fork_time_source = "entropy_shared"`.
#'
#' Legacy estimator (default when no shared coordinate is declared): each
#' lineage is windowed around its own `b_star` and per-lineage changes are
#' subtracted afterwards, which does not guarantee a common time support.
#'
#' @param temporal_scores Temporal LR score table.
#' @param bif_df Lineage-specific fork bins from `FindBifurcationBins()`.
#' @param sibling_map Output from `BuildSiblingContrastMap()`.
#' @param sample_col Optional sample identifier column.
#' @param window_bins Number of bins on each side of the fork.
#' @param activity_floor Minimum activity required for evaluation.
#' @param min_n_eff Minimum row-level effective cell count.
#' @param fork_bins Optional data.frame with columns `lineage` and `fork_bin`
#'   giving the common fork bin of each focal lineage (see
#'   `LocateForkBinsOnGrid()`).  Presence of this argument selects the
#'   shared-window estimator.
#' @param common_grid `NULL` (default) infers the coordinate from the
#'   `temporal_coordinate` attribute written by `ComputeTemporalLR*`;
#'   `TRUE`/`FALSE` override the inference.
#' @return One row per LR trajectory and focal lineage.
#' @export
ComputeForkContrasts <- function(temporal_scores, bif_df, sibling_map,
                                 sample_col = NULL, window_bins = 2L,
                                 activity_floor = 1e-4, min_n_eff = 5,
                                 fork_bins = NULL, common_grid = NULL) {
  score_required <- c(
    "pair_idx", "source", "target", "lineage", "bin", "score"
  )
  bif_required <- c("lineage", "b_star", "reason")
  sibling_required <- c(
    "focal_lineage", "sibling_lineage", "fork_cluster",
    "shared_prefix_length"
  )
  if (!is.data.frame(temporal_scores) ||
      !all(score_required %in% colnames(temporal_scores))) {
    stop("`temporal_scores` is missing required columns", call. = FALSE)
  }
  if (!is.data.frame(bif_df) ||
      !all(bif_required %in% colnames(bif_df))) {
    stop("`bif_df` is missing required columns", call. = FALSE)
  }
  if (!is.data.frame(sibling_map) ||
      !all(sibling_required %in% colnames(sibling_map))) {
    stop("`sibling_map` is missing required columns", call. = FALSE)
  }
  if (!is.null(sample_col) &&
      (length(sample_col) != 1L ||
       !sample_col %in% colnames(temporal_scores))) {
    stop("`sample_col` must name one column in `temporal_scores`",
         call. = FALSE)
  }
  window_bins <- .as_integerish(window_bins, "window_bins")
  if (length(window_bins) != 1L || window_bins < 1L) {
    stop("`window_bins` must be one positive integer", call. = FALSE)
  }
  activity_floor <- as.numeric(activity_floor)
  min_n_eff <- as.numeric(min_n_eff)
  if (length(activity_floor) != 1L || !is.finite(activity_floor) ||
      activity_floor < 0) {
    stop("`activity_floor` must be one finite value >= 0", call. = FALSE)
  }
  if (length(min_n_eff) != 1L || !is.finite(min_n_eff) || min_n_eff < 0) {
    stop("`min_n_eff` must be one finite value >= 0", call. = FALSE)
  }
  if (anyDuplicated(bif_df$lineage)) {
    stop("`bif_df$lineage` must be unique", call. = FALSE)
  }
  if (!is.null(fork_bins)) {
    if (!is.data.frame(fork_bins) ||
        !all(c("lineage", "fork_bin") %in% colnames(fork_bins)) ||
        !is.numeric(fork_bins$fork_bin)) {
      stop("`fork_bins` must be a data.frame with numeric `lineage` and `fork_bin`",
           call. = FALSE)
    }
    if (anyDuplicated(fork_bins$lineage)) {
      stop("`fork_bins$lineage` must be unique", call. = FALSE)
    }
  }
  if (is.null(common_grid)) {
    common_grid <- identical(
      attr(temporal_scores, "temporal_coordinate"),
      "common_pseudotime_grid"
    )
  }
  if (!is.logical(common_grid) || length(common_grid) != 1L ||
      is.na(common_grid)) {
    stop("`common_grid` must be TRUE, FALSE or NULL", call. = FALSE)
  }
  coordinate <- attr(temporal_scores, "temporal_coordinate")
  if (is.null(coordinate) || !is.character(coordinate) ||
      length(coordinate) != 1L) {
    coordinate <- "unspecified"
  }
  shared_windows <- !is.null(fork_bins) || isTRUE(common_grid)

  scores <- data.table::as.data.table(data.table::copy(temporal_scores))
  n_eff_source <- if ("n_eff_cells" %in% colnames(scores)) {
    "reported"
  } else {
    scores$n_eff_cells <- 1
    "unavailable"
  }
  if ("structural_valid" %in% colnames(scores)) {
    scores$structurally_valid <- as.logical(scores$structural_valid)
  } else if (!"structurally_valid" %in% colnames(scores)) {
    scores$structurally_valid <- TRUE
  }
  scores[, b_star := bif_df$b_star[match(lineage, bif_df$lineage)]]
  scores[, bifurcation_reason := bif_df$reason[
    match(lineage, bif_df$lineage)
  ]]
  group_cols <- c(
    "pair_idx", "source", "target", "lineage",
    if (!is.null(sample_col)) sample_col
  )
  scores[, valid_bif := !is.na(b_star) & bifurcation_reason == "crossing"]
  scores[, valid := is.finite(score) & is.finite(n_eff_cells) &
           n_eff_cells >= min_n_eff & !is.na(structurally_valid) &
           structurally_valid]
  scores[, cap := stats::quantile(
    n_eff_cells[valid], 0.95, na.rm = TRUE, names = FALSE
  ), by = group_cols]
  scores[, w_eff := pmin(n_eff_cells, cap)]
  scores[, pre := valid & valid_bif &
           bin >= b_star - window_bins & bin <= b_star - 1L]
  scores[, post := valid & valid_bif &
           bin >= b_star + 1L & bin <= b_star + window_bins]
  scores[, fork := valid & valid_bif & bin == b_star]

  if (shared_windows) {
    anchor <- if (!is.null(fork_bins)) {
      data.table::data.table(
        lineage = as.integer(fork_bins$lineage),
        anchor_bin = as.integer(round(as.numeric(fork_bins$fork_bin))),
        anchor_valid = is.finite(as.numeric(fork_bins$fork_bin))
      )
    } else {
      data.table::data.table(
        lineage = as.integer(bif_df$lineage),
        anchor_bin = as.integer(bif_df$b_star),
        anchor_valid = !is.na(bif_df$b_star) & bif_df$reason == "crossing"
      )
    }
    out <- .fork_contrasts_shared_windows(
      scores, anchor, sibling_map, group_cols, window_bins,
      activity_floor, min_n_eff, n_eff_source, sample_col,
      coordinate,
      fork_time_source = if (!is.null(fork_bins)) {
        "topology"
      } else {
        "entropy_shared"
      }
    )
    return(data.frame(out, stringsAsFactors = FALSE))
  }

  focal <- scores[, .(
    b_star = data.table::first(data.table::fifelse(
      valid_bif, as.integer(b_star), NA_integer_
    )),
    valid_bif = data.table::first(valid_bif),
    pre_count = sum(pre),
    pre_numerator = sum(data.table::fifelse(pre, score * w_eff, 0)),
    pre_denominator = sum(data.table::fifelse(pre, w_eff, 0)),
    fork_count = sum(fork),
    fork_numerator = sum(data.table::fifelse(fork, score * w_eff, 0)),
    fork_denominator = sum(data.table::fifelse(fork, w_eff, 0)),
    post_count = sum(post),
    post_numerator = sum(data.table::fifelse(post, score * w_eff, 0)),
    post_denominator = sum(data.table::fifelse(post, w_eff, 0)),
    post_score_sum = sum(data.table::fifelse(post, score, 0)),
    post_score_sum_squares = sum(data.table::fifelse(post, score * score, 0)),
    n_pre_bins = data.table::uniqueN(bin[pre]),
    n_post_bins = data.table::uniqueN(bin[post]),
    effective_weight = sum(data.table::fifelse(valid, w_eff, 0))
  ), by = group_cols]
  focal[, pre_activity := data.table::fifelse(
    pre_count > 0L, pre_numerator / pre_denominator, NA_real_
  )]
  focal[, fork_activity := data.table::fifelse(
    fork_count > 0L, fork_numerator / fork_denominator, NA_real_
  )]
  focal[, post_activity := data.table::fifelse(
    post_count > 0L, post_numerator / post_denominator, NA_real_
  )]
  focal[, post_mean := data.table::fifelse(
    post_count > 0L, post_score_sum / post_count, NA_real_
  )]
  focal[, post_sd := data.table::fifelse(
    post_count > 1L,
    sqrt(pmax(
      0,
      (post_score_sum_squares - post_score_sum * post_score_sum / post_count) /
        (post_count - 1L)
    )),
    NA_real_
  )]
  focal[, activity_strength := pmax(
    pre_activity, fork_activity, post_activity, na.rm = TRUE
  )]
  focal[!is.finite(activity_strength), activity_strength := NA_real_]
  focal[, reason := data.table::fifelse(
    !valid_bif,
    "invalid_bifurcation",
    data.table::fifelse(
      n_pre_bins != window_bins,
      "insufficient_pre_bins",
      data.table::fifelse(
        n_post_bins != window_bins,
        "insufficient_post_bins",
        data.table::fifelse(
          !is.finite(activity_strength) | activity_strength < activity_floor,
          "below_activity_floor",
          "evaluable"
        )
      )
    )
  )]
  focal[, evaluable := reason == "evaluable"]
  focal[, post_fork_change := data.table::fifelse(
    evaluable, post_activity - pre_activity, NA_real_
  )]
  focal[, pre_fork_enrichment := pre_activity /
          (pre_activity + post_activity + 1e-12)]
  focal[, fork_locality := fork_activity /
          (pre_activity + fork_activity + post_activity + 1e-12)]
  focal[, post_fork_persistence := data.table::fifelse(
    is.finite(post_sd),
    1 - post_sd / (abs(post_mean) + post_sd + 1e-12),
    NA_real_
  )]
  focal[, `:=`(
    n_eff_source = n_eff_source,
    inference_unit = data.table::fifelse(
      is.null(sample_col),
      "pooled_cells_descriptive",
      "sample_rescore_descriptive"
    )
  )]
  focal[, c(
    "valid_bif", "pre_count", "pre_numerator", "pre_denominator",
    "fork_count", "fork_numerator", "fork_denominator", "post_count",
    "post_numerator", "post_denominator", "post_score_sum",
    "post_score_sum_squares", "post_mean", "post_sd"
  ) := NULL]
  data.table::setorderv(focal, group_cols)

  mapped <- merge(
    focal, data.table::as.data.table(data.table::copy(sibling_map)),
    by.x = "lineage", by.y = "focal_lineage", all.x = TRUE,
    sort = FALSE, allow.cartesian = TRUE
  )
  sibling_keys <- c("pair_idx", "source", "target", "lineage")
  if (!is.null(sample_col)) sibling_keys <- c(sibling_keys, sample_col)
  sibling <- focal[, c(
    sibling_keys, "post_fork_change", "effective_weight", "evaluable"
  ), with = FALSE]
  data.table::setnames(
    sibling,
    c("lineage", "post_fork_change", "effective_weight", "evaluable"),
    c(
      "sibling_lineage", "sibling_post_fork_change",
      "sibling_effective_weight", "sibling_evaluable"
    )
  )
  join_keys <- c("pair_idx", "source", "target", "sibling_lineage")
  if (!is.null(sample_col)) join_keys <- c(join_keys, sample_col)
  mapped <- merge(mapped, sibling, by = join_keys, all.x = TRUE, sort = FALSE)
  focal_id_cols <- c(group_cols, "reason")
  mapped[, sibling_valid := !is.na(sibling_evaluable) & sibling_evaluable &
           is.finite(sibling_post_fork_change) &
           is.finite(sibling_effective_weight) & sibling_effective_weight > 0]
  mapped[, sibling_cap := stats::quantile(
    sibling_effective_weight[sibling_valid], 0.95,
    names = FALSE, na.rm = TRUE
  ), by = focal_id_cols]
  mapped[, sibling_weight := pmin(sibling_effective_weight, sibling_cap)]
  sibling_summary <- mapped[, .(
    sibling_valid_count = sum(sibling_valid),
    sibling_numerator = sum(data.table::fifelse(
      sibling_valid, sibling_post_fork_change * sibling_weight, 0
    )),
    sibling_denominator = sum(data.table::fifelse(
      sibling_valid, sibling_weight, 0
    ))
  ), by = focal_id_cols]
  sibling_summary[, sibling_change := data.table::fifelse(
    sibling_valid_count > 0L,
    sibling_numerator / sibling_denominator,
    NA_real_
  )]
  sibling_summary[, c(
    "sibling_valid_count", "sibling_numerator", "sibling_denominator"
  ) := NULL]
  out <- merge(
    focal, sibling_summary,
    by = focal_id_cols, all.x = TRUE, sort = FALSE
  )
  out[, sibling_contrast := data.table::fifelse(
    evaluable & is.finite(sibling_change),
    post_fork_change - sibling_change,
    NA_real_
  )]
  out[evaluable & !is.finite(sibling_change), `:=`(
    evaluable = FALSE,
    reason = "no_evaluable_sibling"
  )]
  fork_cluster_tbl <- data.table::as.data.table(sibling_map)[
    , .(fork_cluster = data.table::first(fork_cluster)),
    by = focal_lineage
  ]
  out <- merge(
    out, fork_cluster_tbl,
    by.x = "lineage", by.y = "focal_lineage", all.x = TRUE, sort = FALSE
  )
  out[, `:=`(
    fork_time_source = "per_lineage",
    temporal_coordinate = coordinate,
    n_pre_bins_shared = NA_integer_,
    n_post_bins_shared = NA_integer_,
    n_sisters_mapped = NA_integer_
  )]
  data.table::setorderv(out, group_cols)
  out <- out[, c(
    group_cols, "fork_cluster", "b_star", "pre_activity", "fork_activity",
    "post_activity", "post_fork_change", "sibling_change",
    "sibling_contrast", "activity_strength", "pre_fork_enrichment",
    "fork_locality", "post_fork_persistence", "n_pre_bins",
    "n_post_bins", "n_pre_bins_shared", "n_post_bins_shared",
    "n_sisters_mapped", "effective_weight", "evaluable", "reason",
    "fork_time_source", "temporal_coordinate",
    "n_eff_source", "inference_unit"
  ), with = FALSE]
  data.frame(out, stringsAsFactors = FALSE)
}

#' Shared-window fork contrast estimator (internal).
#'
#' One common fork bin `anchor_bin` per focal lineage anchors pre/post windows
#' for the focal and every mapped sister.  The effect is
#' `theta = mean_post(D) - mean_pre(D)` where `D(t) = focal(t) - mean_sister(t)`
#' is evaluated on bins where the focal and at least one sister are both
#' supported; all branches use identical uniform time weights and sisters use
#' equal weights.
.fork_contrasts_shared_windows <- function(scores, anchor, sibling_map,
                                           group_cols, window_bins,
                                           activity_floor, min_n_eff,
                                           n_eff_source, sample_col,
                                           coordinate, fork_time_source) {
  unit_cols <- c(
    "pair_idx", "source", "target",
    if (!is.null(sample_col)) sample_col
  )

  long <- scores[, c(
    unit_cols, "lineage", "bin", "score", "valid", "n_eff_cells"
  ), with = FALSE]

  focal <- unique(long[, c(unit_cols, "lineage"), with = FALSE])
  focal <- merge(
    focal, anchor, by = "lineage", all.x = TRUE, sort = FALSE
  )
  focal[is.na(anchor_valid), anchor_valid := FALSE]

  smap <- data.table::as.data.table(sibling_map)[
    , .(focal_lineage = as.integer(focal_lineage),
        sibling_lineage = as.integer(sibling_lineage),
        fork_cluster = as.integer(fork_cluster))
  ]
  sister_counts <- smap[
    , .(n_sisters_mapped = .N, fork_cluster = data.table::first(fork_cluster)),
    by = focal_lineage
  ]
  focal <- merge(
    focal, sister_counts,
    by.x = "lineage", by.y = "focal_lineage",
    all.x = TRUE, sort = FALSE
  )
  focal[is.na(n_sisters_mapped), n_sisters_mapped := 0L]

  # Expand each valid anchor into its window bins without evaluating `seq()`
  # inside a grouped `j`: with zero valid anchors (e.g. every lineage ends in a
  # non-crossing bifurcation) data.table still probes `j` once on the empty
  # subset, which made `seq(NA, NA)` abort instead of returning no-call rows.
  win_keys <- focal[
    anchor_valid %in% TRUE,
    c(unit_cols, "lineage", "anchor_bin"),
    with = FALSE
  ]
  win <- win_keys[rep(seq_len(nrow(win_keys)), each = length(seq.int(
    -window_bins, window_bins
  )))]
  win[, bin := anchor_bin + rep(
    seq.int(-window_bins, window_bins),
    times = nrow(win_keys)
  )]

  win <- merge(
    win, long,
    by = c(unit_cols, "lineage", "bin"),
    all.x = TRUE, sort = FALSE
  )
  win[, focal_valid := !is.na(valid) & valid]
  win[, window_side := data.table::fcase(
    bin < anchor_bin, "pre",
    bin > anchor_bin, "post",
    default = "fork"
  )]
  win[, cap := stats::quantile(
    n_eff_cells[focal_valid], 0.95, na.rm = TRUE, names = FALSE
  ), by = c(unit_cols, "lineage")]

  sis <- merge(
    win[, c(unit_cols, "lineage", "bin"), with = FALSE],
    smap,
    by.x = "lineage", by.y = "focal_lineage",
    allow.cartesian = TRUE, sort = FALSE
  )
  sis <- merge(
    sis, long,
    by.x = c(unit_cols, "sibling_lineage", "bin"),
    by.y = c(unit_cols, "lineage", "bin"),
    all.x = TRUE, sort = FALSE
  )
  sis[, sister_valid := !is.na(valid) & valid]
  sis_bin <- sis[, .(
    n_sisters_valid = sum(sister_valid),
    sister_mean = if (sum(sister_valid) > 0L) {
      mean(score[sister_valid])
    } else {
      NA_real_
    }
  ), by = c(unit_cols, "lineage", "bin")]

  win <- merge(
    win, sis_bin,
    by = c(unit_cols, "lineage", "bin"),
    all.x = TRUE, sort = FALSE
  )
  win[is.na(n_sisters_valid), n_sisters_valid := 0L]
  win[, shared := focal_valid & n_sisters_valid > 0L]

  agg <- win[, .(
    n_pre_bins = sum(focal_valid & window_side == "pre"),
    n_post_bins = sum(focal_valid & window_side == "post"),
    n_pre_bins_shared = sum(shared & window_side == "pre"),
    n_post_bins_shared = sum(shared & window_side == "post"),
    pre_activity = if (sum(focal_valid & window_side == "pre") > 0L) {
      mean(score[focal_valid & window_side == "pre"])
    } else {
      NA_real_
    },
    fork_activity = if (sum(focal_valid & window_side == "fork") > 0L) {
      mean(score[focal_valid & window_side == "fork"])
    } else {
      NA_real_
    },
    post_activity = if (sum(focal_valid & window_side == "post") > 0L) {
      mean(score[focal_valid & window_side == "post"])
    } else {
      NA_real_
    },
    post_sd = if (sum(focal_valid & window_side == "post") > 1L) {
      stats::sd(score[focal_valid & window_side == "post"])
    } else {
      NA_real_
    },
    pre_f_shared = if (sum(shared & window_side == "pre") > 0L) {
      mean(score[shared & window_side == "pre"])
    } else {
      NA_real_
    },
    post_f_shared = if (sum(shared & window_side == "post") > 0L) {
      mean(score[shared & window_side == "post"])
    } else {
      NA_real_
    },
    pre_s_shared = if (sum(shared & window_side == "pre") > 0L) {
      mean(sister_mean[shared & window_side == "pre"])
    } else {
      NA_real_
    },
    post_s_shared = if (sum(shared & window_side == "post") > 0L) {
      mean(sister_mean[shared & window_side == "post"])
    } else {
      NA_real_
    },
    effective_weight = sum(data.table::fifelse(
      focal_valid, pmin(n_eff_cells, cap), 0
    ), na.rm = TRUE)
  ), by = c(unit_cols, "lineage")]

  out <- merge(
    focal, agg,
    by = c(unit_cols, "lineage"),
    all.x = TRUE, sort = FALSE
  )
  count_cols <- c(
    "n_pre_bins", "n_post_bins",
    "n_pre_bins_shared", "n_post_bins_shared"
  )
  out[, (count_cols) := lapply(.SD, function(x) {
    data.table::fifelse(is.na(x), 0L, x)
  }), .SDcols = count_cols]

  out[, activity_strength := pmax(
    pre_activity, fork_activity, post_activity, na.rm = TRUE
  )]
  out[!is.finite(activity_strength), activity_strength := NA_real_]
  out[, post_fork_change := post_f_shared - pre_f_shared]
  out[, sibling_change := post_s_shared - pre_s_shared]
  out[, sibling_contrast := post_fork_change - sibling_change]
  out[, pre_fork_enrichment := pre_activity /
        (pre_activity + post_activity + 1e-12)]
  out[, fork_locality := fork_activity /
        (pre_activity + fork_activity + post_activity + 1e-12)]
  out[, post_fork_persistence := data.table::fifelse(
    is.finite(post_sd),
    1 - post_sd / (abs(post_activity) + post_sd + 1e-12),
    NA_real_
  )]

  out[, reason := data.table::fcase(
    !anchor_valid, if (fork_time_source == "topology") {
      "fork_bin_unknown"
    } else {
      "invalid_bifurcation"
    },
    n_pre_bins < window_bins, "insufficient_pre_bins",
    n_post_bins < window_bins, "insufficient_post_bins",
    n_pre_bins_shared < window_bins, "insufficient_shared_pre_bins",
    n_post_bins_shared < window_bins, "insufficient_shared_post_bins",
    !is.finite(activity_strength) | activity_strength < activity_floor,
    "below_activity_floor",
    default = "evaluable"
  )]
  out[, evaluable := reason == "evaluable"]
  out[evaluable == FALSE, `:=`(
    post_fork_change = NA_real_,
    sibling_change = NA_real_,
    sibling_contrast = NA_real_
  )]
  out[, `:=`(
    n_eff_source = n_eff_source,
    inference_unit = data.table::fifelse(
      is.null(sample_col),
      "pooled_cells_descriptive",
      "sample_rescore_descriptive"
    ),
    fork_time_source = fork_time_source,
    temporal_coordinate = coordinate
  )]

  data.table::setnames(out, "anchor_bin", "b_star")
  out[, c(
    "anchor_valid", "post_sd", "pre_f_shared", "post_f_shared",
    "pre_s_shared", "post_s_shared"
  ) := NULL]
  data.table::setorderv(out, group_cols)
  out[, c(
    group_cols, "fork_cluster", "b_star", "pre_activity", "fork_activity",
    "post_activity", "post_fork_change", "sibling_change",
    "sibling_contrast", "activity_strength", "pre_fork_enrichment",
    "fork_locality", "post_fork_persistence", "n_pre_bins",
    "n_post_bins", "n_pre_bins_shared", "n_post_bins_shared",
    "n_sisters_mapped", "effective_weight", "evaluable", "reason",
    "fork_time_source", "temporal_coordinate",
    "n_eff_source", "inference_unit"
  ), with = FALSE]
}

#' Summarize fork contrasts across biological samples.
#'
#' Samples, not cells, are the resampling unit. The function reports sign
#' fractions and a bootstrap interval for the sample median, but no cell-level
#' inferential p-value.
#'
#' @param sample_contrasts One evaluable contrast row per trajectory and sample.
#' @param n_boot Number of sample bootstrap replicates.
#' @param seed Random seed.
#' @return Sample-level descriptive summary.
#' @export
SummarizeForkContrastsBySample <- function(sample_contrasts,
                                           n_boot = 999L,
                                           seed = 12345L,
                                           by = NULL) {
  if (is.null(by)) {
    by <- c("pair_idx", "source", "target", "lineage")
    event_key <- c("event_node", "child")
    if (all(event_key %in% colnames(sample_contrasts))) {
      by <- c(by, event_key)
    }
  }
  required <- c(
    by, "block_label", "sibling_contrast", "evaluable"
  )
  if (!is.data.frame(sample_contrasts) ||
      !all(required %in% colnames(sample_contrasts))) {
    stop("`sample_contrasts` is missing required columns", call. = FALSE)
  }
  n_boot <- .validate_positive_int(n_boot, "n_boot", min = 1L)
  seed <- .validate_seed(seed)
  x <- sample_contrasts[
    sample_contrasts$evaluable %in% TRUE &
      is.finite(sample_contrasts$sibling_contrast),
    , drop = FALSE
  ]
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  if (!nrow(x)) {
    stop("no evaluable sample contrasts", call. = FALSE)
  }
  key_parts <- lapply(x[by], function(v) {
    v <- as.character(v)
    v[is.na(v)] <- "\x01NA"
    v
  })
  key <- do.call(paste, c(key_parts, sep = "\x01"))
  groups <- split(seq_len(nrow(x)), key)
  duplicated_samples <- vapply(groups, function(idx) {
    anyDuplicated(x$block_label[idx]) > 0L
  }, logical(1L))
  if (any(duplicated_samples)) {
    stop("each trajectory must have at most one row per `block_label`",
         call. = FALSE)
  }
  set.seed(seed)
  rows <- lapply(groups, function(idx) {
    g <- x[idx, , drop = FALSE]
    values <- as.numeric(g$sibling_contrast)
    n_samples <- length(values)
    median_contrast <- stats::median(values)
    positive_fraction <- mean(values > 0)
    negative_fraction <- mean(values < 0)
    sample_support <- if (median_contrast > 0) {
      positive_fraction
    } else if (median_contrast < 0) {
      negative_fraction
    } else {
      max(positive_fraction, negative_fraction)
    }
    boot <- if (n_samples >= 2L) {
      # Consecutive blocks retain the original draws and bootstrap statistic.
      draws <- sample(values, n_samples * n_boot, replace = TRUE)
      .Call(wrap__rs_bootstrap_medians, draws, as.integer(n_samples))
    } else {
      numeric()
    }
    ci <- if (length(boot)) {
      stats::quantile(
        boot, c(0.025, 0.975), names = FALSE, na.rm = TRUE
      )
    } else {
      c(NA_real_, NA_real_)
    }
    cbind(
      g[1L, by, drop = FALSE],
      data.frame(
        median_contrast = median_contrast,
        sample_support = sample_support,
        positive_sample_fraction = positive_fraction,
        negative_sample_fraction = negative_fraction,
        n_samples = n_samples,
        ci_low = ci[[1L]],
        ci_high = ci[[2L]],
        inference_unit = if (n_samples >= 2L) {
          "biological_sample_bootstrap"
        } else {
          "single_sample_descriptive"
        },
        stringsAsFactors = FALSE
      )
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' Compose fork-local evidence without assigning mechanistic roles.
#'
#' Missing response or sample support remains missing rather than being imputed
#' as zero. All merges use the complete LR trajectory key.
#'
#' @param fork_contrasts Output from `ComputeForkContrasts()`.
#' @param response_support Optional receiver-response evidence table.
#' @param sample_support Optional output from
#'   `SummarizeForkContrastsBySample()`.
#' @param nomination_groups Optional query table with context, target and nonnegative
#'   weight. Weights are normalized within each prespecified receiver population.
#' @param pair_gene_map Molecular map with pair_idx, ligand_genes and receptor_genes.
#' @param nomination_candidates Optional fixed context, ligand, receptor panel.
#' @return Fork evidence table with an explicit completeness label. Optional
#'   ligand_nominations and nomination_receiver_evidence attributes retain every
#'   planned receiver. Nomination orders Pareto fronts of peak and population
#'   effects, then the original peak effect within each front. Missing receiver
#'   evidence contributes no nomination mass and retains an explicit missing value;
#'   missing response and sample support are unchanged.
#' @export
ComposeForkEvidence <- function(fork_contrasts, response_support = NULL,
                                sample_support = NULL, nomination_groups = NULL,
                                pair_gene_map = NULL, nomination_candidates = NULL) {
  key <- c("pair_idx", "source", "target", "lineage")
  event_key <- c("event_node", "child")
  has_event_identity <- all(event_key %in% colnames(fork_contrasts))
  if (has_event_identity) {
    key <- c(key, event_key)
  }
  fork_required <- c(
    key, "activity_strength", "pre_fork_enrichment", "fork_locality",
    "post_fork_change", "sibling_contrast", "post_fork_persistence",
    "evaluable"
  )
  if (!is.data.frame(fork_contrasts) ||
      !all(fork_required %in% colnames(fork_contrasts))) {
    stop("`fork_contrasts` is missing required columns", call. = FALSE)
  }
  out <- fork_contrasts
  out$.input_order <- seq_len(nrow(out))
  if (!is.null(response_support)) {
    response_required <- c(key, "receiver_response_support")
    if (!is.data.frame(response_support) ||
        !all(c("pair_idx", "source", "target", "lineage",
               "receiver_response_support") %in%
             colnames(response_support))) {
      stop("`response_support` is missing required columns", call. = FALSE)
    }
    if (has_event_identity &&
        !all(event_key %in% colnames(response_support))) {
      warning(
        "`response_support` lacks the fork event identity; response ",
        "evidence stays NA instead of merging on an incomplete key",
        call. = FALSE
      )
      out$receiver_response_support <- NA_real_
    } else {
      if (anyDuplicated(response_support[key])) {
        stop("`response_support` trajectory keys must be unique",
             call. = FALSE)
      }
      out <- merge(
        out,
        response_support[, response_required, drop = FALSE],
        by = key, all.x = TRUE, sort = FALSE
      )
    }
  } else {
    out$receiver_response_support <- NA_real_
  }
  if (!is.null(sample_support)) {
    sample_required <- c(key, "sample_support")
    if (!is.data.frame(sample_support) ||
        !all(c("pair_idx", "source", "target", "lineage",
               "sample_support") %in% colnames(sample_support))) {
      stop("`sample_support` is missing required columns", call. = FALSE)
    }
    if (has_event_identity &&
        !all(event_key %in% colnames(sample_support))) {
      warning(
        "`sample_support` lacks the fork event identity; sample evidence ",
        "stays NA instead of merging on an incomplete key",
        call. = FALSE
      )
      if (!"sample_support" %in% colnames(out)) {
        out$sample_support <- NA_real_
      }
    } else {
      if (anyDuplicated(sample_support[key])) {
        stop("`sample_support` trajectory keys must be unique", call. = FALSE)
      }
      sample_data <- sample_support
      sample_inference_index <- match(
        "inference_unit", colnames(sample_data), nomatch = 0L
      )
      if (sample_inference_index > 0L) {
        colnames(sample_data)[sample_inference_index] <-
          "sample_inference_unit"
      }
      sample_cols <- unique(c(
        key,
        intersect(
          c(
            "median_contrast", "sample_support",
            "positive_sample_fraction", "negative_sample_fraction",
            "n_samples", "ci_low", "ci_high", "sample_inference_unit"
          ),
          colnames(sample_data)
        )
      ))
      collisions <- intersect(
        setdiff(sample_cols, key),
        colnames(out)
      )
      if (length(collisions)) {
        stop("sample evidence columns already exist in `fork_contrasts`",
             call. = FALSE)
      }
      out <- merge(
        out,
        sample_data[, sample_cols, drop = FALSE],
        by = key, all.x = TRUE, sort = FALSE
      )
    }
  } else if (!"sample_support" %in% colnames(out)) {
    out$sample_support <- NA_real_
  }
  response_present <- is.finite(out$receiver_response_support)
  sample_present <- is.finite(out$sample_support)
  out$evidence_completeness <- ifelse(
    response_present & sample_present,
    "fork+response+sample",
    ifelse(
      response_present,
      "fork+response",
      ifelse(sample_present, "fork+sample", "fork_only")
    )
  )
  out$sample_supported_score <- switch(
    as.character(all(c("median_contrast", "n_samples") %in% colnames(out))),
    "TRUE" = out$median_contrast * out$n_samples,
    "FALSE" = rep(NA_real_, nrow(out))
  )
  out <- out[order(out$.input_order), , drop = FALSE]
  out$.input_order <- NULL
  rownames(out) <- NULL
  if (!is.null(nomination_groups) && nrow(out) > 0L) {
    groups <- data.table::copy(data.table::as.data.table(nomination_groups))
    groups[, weight := weight / sum(weight), by = context]
    mapping <- data.table::as.data.table(pair_gene_map)[,
      .(pair_idx, ligand = ligand_genes, receptor = receptor_genes)]
    candidates <- switch(as.character(is.null(nomination_candidates)),
      "TRUE" = merge(unique(groups[, .(context, join_key = 1L)]),
        unique(mapping[, .(ligand, receptor, join_key = 1L)]),
        by = "join_key", allow.cartesian = TRUE)[, join_key := NULL],
      "FALSE" = data.table::copy(data.table::as.data.table(nomination_candidates)))
    fields <- intersect(c("event_ranking_score", "no_fork_ranking_score",
      "sample_unadjusted_ranking_score"), names(out))
    events <- merge(data.table::as.data.table(out)[,
      c("source", "target", "pair_idx", "evaluable", fields), with = FALSE],
      mapping, by = "pair_idx", sort = FALSE)
    events <- merge(events, groups[, .(context, target, weight)], by = "target",
      allow.cartesian = TRUE)
    events <- merge(events, candidates, by = c("context", "ligand", "receptor"))
    values <- data.table::melt(events, id.vars = c("source", "target", "context",
      "ligand", "receptor", "evaluable", "weight"), measure.vars = fields,
      variable.name = "method", value.name = "effect", variable.factor = FALSE)
    receiver <- values[evaluable %in% TRUE & is.finite(effect),
      .(effect = max(effect)), by = .(source, method, context, ligand, target)]
    grid <- merge(data.table::CJ(source = unique(out$source), method = fields)[, join_key := 1L],
      unique(candidates[, .(context, ligand, join_key = 1L)]),
      by = "join_key", allow.cartesian = TRUE)[, join_key := NULL]
    grid <- merge(grid, groups, by = "context", allow.cartesian = TRUE)
    planned <- merge(grid, receiver, by = c("source", "method", "context", "ligand", "target"),
      all.x = TRUE, sort = FALSE)
    nomination <- planned[, .(
      peak_effect = data.table::fifelse(any(is.finite(effect)), max(data.table::fcoalesce(effect, -Inf)), NA_real_),
      population_effect = data.table::fifelse(any(is.finite(effect)),
        sum(data.table::fcoalesce(effect, 0) * weight), NA_real_),
      n_planned_receivers = .N, n_observed_receivers = sum(is.finite(effect)),
      observed_weight = sum(weight[is.finite(effect)])), by = .(source, method, context, ligand)]
    nomination[, c("nomination_front", "nomination_order_score") := {
      peak <- data.table::fcoalesce(peak_effect, -Inf)
      population <- data.table::fcoalesce(population_effect, -Inf)
      dominates <- outer(peak, peak, ">=") & outer(population, population, ">=") &
        (outer(peak, peak, ">") | outer(population, population, ">"))
      graph <- igraph::graph_from_adjacency_matrix(-dominates,
        mode = "directed", weighted = TRUE, diag = FALSE)
      graph <- igraph::add_vertices(graph, 1L)
      graph <- igraph::add_edges(graph, as.vector(rbind(rep(.N + 1L, .N), seq_len(.N))),
        attr = list(weight = rep(0, .N)))
      depth <- as.numeric(igraph::distances(graph, v = .N + 1L, to = seq_len(.N),
        weights = igraph::edge_attr(graph, "weight"), mode = "out", algorithm = "bellman-ford"))
      score <- depth + (data.table::frank(peak_effect, ties.method = "average") - 1) / .N
      list(data.table::fifelse(is.finite(peak_effect), 1 - depth, NA_real_),
        data.table::fifelse(is.finite(peak_effect), score, NA_real_))
    }, by = .(source, method, context)]
    nomination[, nomination_rank := data.table::frank(
      -data.table::fcoalesce(nomination_order_score, -Inf), ties.method = "min"),
      by = .(source, method, context)]
    nomination[, prediction_status := data.table::fifelse(n_observed_receivers == 0L,
      "no_observed_receiver", data.table::fifelse(n_observed_receivers < n_planned_receivers,
        "partial_receiver_coverage", "complete_receiver_coverage"))]
    attr(out, "ligand_nominations") <- as.data.frame(nomination)
    attr(out, "nomination_receiver_evidence") <- as.data.frame(planned)
  }
  out
}

#' Locate topology fork clusters on a shared pseudotime grid.
#'
#' Each focal lineage's fork event is the last shared cluster of its
#' longest-prefix sibling group (`sibling_map$fork_cluster`).  The fork bin is
#' the bin containing the median pseudotime of the cells assigned to that
#' cluster, evaluated on the common boundary vector — so the returned anchor
#' denotes the same topology event for the focal and its sisters.
#'
#' @param sibling_map Output from `BuildSiblingContrastMap()`.
#' @param cluster_labels Cluster label per cell (values matching the path
#'   entries used to build `sibling_map`).
#' @param pseudotime numeric vector of length `length(cluster_labels)` on the
#'   same `[0,1]` graph-time coordinate the grid was built from.
#' @param bin_bounds Common boundary vector of length `n_bins + 1`
#'   (`ComputeCommonBinBounds()` output).  A matrix or list of identical
#'   per-lineage rows is accepted and checked for equality. Repeated boundaries
#'   represent empty bins and retain the scorer's rightmost-bin convention.
#' @return data.frame with columns `lineage`, `fork_cluster`,
#'   `fork_pseudotime`, `fork_bin` (1-based; `NA` when the fork cluster has no
#'   cells with finite pseudotime).  Pass `lineage`/`fork_bin` as the
#'   `fork_bins` argument of `ComputeForkContrasts()`.
#' @export
LocateForkBinsOnGrid <- function(sibling_map, cluster_labels, pseudotime,
                                 bin_bounds) {
  if (!is.data.frame(sibling_map) ||
      !all(c("focal_lineage", "fork_cluster") %in% colnames(sibling_map))) {
    stop("`sibling_map` is missing required columns", call. = FALSE)
  }
  if (is.matrix(bin_bounds) || is.data.frame(bin_bounds)) {
    bound_rows <- lapply(seq_len(nrow(bin_bounds)), function(i) {
      as.numeric(bin_bounds[i, ])
    })
  } else if (is.list(bin_bounds)) {
    bound_rows <- lapply(bin_bounds, as.numeric)
  } else {
    bound_rows <- list(as.numeric(bin_bounds))
  }
  if (length(bound_rows) > 1L &&
      !all(vapply(bound_rows[-1L], function(row) {
        isTRUE(all.equal(row, bound_rows[[1L]], tolerance = 1e-12))
      }, logical(1L)))) {
    stop("`bin_bounds` must be a common grid (identical rows)", call. = FALSE)
  }
  bounds <- bound_rows[[1L]]
  n_bins <- length(bounds) - 1L
  if (length(bounds) < 3L || any(!is.finite(bounds)) || any(diff(bounds) < 0)) {
    stop("`bin_bounds` must be a finite non-decreasing boundary vector",
         call. = FALSE)
  }
  if (length(cluster_labels) != length(pseudotime)) {
    stop("`cluster_labels` and `pseudotime` lengths differ", call. = FALSE)
  }
  pt <- as.numeric(pseudotime)

  forks <- unique(sibling_map[, c("focal_lineage", "fork_cluster")])
  rows <- lapply(seq_len(nrow(forks)), function(i) {
    cells <- which(cluster_labels == forks$fork_cluster[i] & is.finite(pt))
    fork_pt <- if (length(cells)) stats::median(pt[cells]) else NA_real_
    fork_bin <- if (is.finite(fork_pt)) {
      b <- findInterval(fork_pt, bounds, rightmost.closed = TRUE)
      as.integer(max(1L, min(b, n_bins)))
    } else {
      NA_integer_
    }
    data.frame(
      lineage = as.integer(forks$focal_lineage[i]),
      fork_cluster = as.integer(forks$fork_cluster[i]),
      fork_pseudotime = fork_pt,
      fork_bin = fork_bin,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' Re-estimate no-call shared-window fork contrasts on merged adjacent bins.
#'
#' Fixed per-side windows reject a side when it holds fewer than `window_bins`
#' supported bins even though a shorter continuous interval has enough cells.
#' This function keeps the same fork, global time grid, Hill scale and fixed
#' sister set, and for a thin side merges the supported bins of that side
#' (never crossing the fork and never extending beyond the original window)
#' into one explicit interval. Ligand/receptor means are recomputed from the
#' raw weighted expression on that interval via `ComputeTemporalLRFixedBounds`
#' (the interval is expressed as the middle bin of the 3-bin grid
#' `[0, lo, hi, 1]`, so the pooled Hill scale and support gates are reused
#' unchanged). Boundaries depend only on cell time and support - never on LR
#' effect, truth or significance. If any fixed sister lacks support on a side,
#' theta stays NA; a resolvable post side is reported separately as
#' `post_branch_difference`.
#'
#' @param temporal_scores Common-grid temporal LR score table.
#' @param bif_df Bifurcation table (passed through to the fixed-window path).
#' @param sibling_map Focal-to-sibling lineage map with `fork_cluster`.
#' @param expr Log-normalized expression matrix used for the scores.
#' @param celltypes Integer cell-type vector.
#' @param lr_db Ligand-receptor database rows.
#' @param pseudotime Pseudotime vector used for the scores.
#' @param lineage_weights Cell-by-lineage weight matrix.
#' @param fork_bins One topology fork bin per focal lineage.
#' @param window_bins Fixed window half-width in bins.
#' @param min_n_eff Minimum effective cells per side.
#' @param min_pct Minimum weighted expressing fraction.
#' @param activity_floor Activity floor used by the fixed-window path.
#' @param hill_scale Pooled Hill scale; defaults to the scores attribute.
#' @return Fixed-window contrast table with `estimation_mode`, merged-interval
#'   metadata and `post_branch_difference` columns.
#' @export
ComputeMergedIntervalForkContrasts <- function(
    temporal_scores, bif_df, sibling_map, expr, celltypes, lr_db,
    pseudotime, lineage_weights, fork_bins,
    window_bins = 2L, min_n_eff = 5.0, min_pct = 0.05,
    activity_floor = 1e-4, hill_scale = NULL, sample_col = NULL) {
  base <- data.table::as.data.table(ComputeForkContrasts(
    temporal_scores, bif_df, sibling_map,
    window_bins = window_bins, activity_floor = activity_floor,
    min_n_eff = min_n_eff, fork_bins = fork_bins, common_grid = TRUE,
    sample_col = sample_col
  ))
  base[, `:=`(
    estimation_mode = "not_evaluable",
    merged_pre_edges = NA_character_,
    merged_post_edges = NA_character_,
    merged_pre_bins = 0L,
    merged_post_bins = 0L,
    post_branch_difference = NA_real_
  )]
  base[evaluable %in% TRUE, estimation_mode := "fixed_window"]

  bounds <- attr(temporal_scores, "bin_bounds")
  if (is.null(hill_scale)) hill_scale <- attr(temporal_scores, "hill_scale")
  ct_levels <- attr(temporal_scores, "celltype_levels")
  n_lin <- ncol(lineage_weights)
  fb <- data.table::as.data.table(fork_bins)
  smap <- data.table::as.data.table(sibling_map)
  scores <- data.table::as.data.table(data.table::copy(temporal_scores))
  if (!"structurally_valid" %in% colnames(scores)) {
    scores[, structurally_valid := TRUE]
  }
  scores[, valid := is.finite(score) & is.finite(n_eff_cells) &
           n_eff_cells >= min_n_eff & !is.na(structurally_valid) &
           structurally_valid]
  scores[, anchor_bin := as.integer(round(
    fb$fork_bin[match(lineage, fb$lineage)]
  ))]

  group_cols <- c("pair_idx", "source", "target", "lineage",
                  if (!is.null(sample_col)) sample_col)
  side_bins <- scores[!is.na(anchor_bin) & valid %in% TRUE, {
    a <- anchor_bin[1L]
    .(pre_bins = list(sort(bin[bin >= a - window_bins & bin <= a - 1L])),
      post_bins = list(sort(bin[bin >= a + 1L & bin <= a + window_bins])))
  }, by = group_cols]
  base <- merge(base, side_bins, by = group_cols, all.x = TRUE, sort = FALSE)

  needs_merge <- !(base$evaluable %in% TRUE) &
    base$reason %in% c("insufficient_pre_bins", "insufficient_post_bins",
                       "insufficient_shared_pre_bins",
                       "insufficient_shared_post_bins")
  merge_rows <- which(needs_merge)

  to_interval <- function(lin, bins) {
    edges <- bounds[[lin]]
    c(edges[min(bins)], edges[max(bins) + 1L])
  }
  intervals <- rbindlist(lapply(merge_rows, function(i) {
    lin <- base$lineage[i]
    out <- list()
    if (length(base$pre_bins[[i]]) >= 1L &&
        !(is.finite(base$n_pre_bins_shared[i]) &&
          base$n_pre_bins_shared[i] >= window_bins)) {
      iv <- to_interval(lin, base$pre_bins[[i]])
      out[[length(out) + 1L]] <- data.table(
        lineage = lin, lo = iv[1L], hi = iv[2L], side = "pre", row = i
      )
    }
    if (length(base$post_bins[[i]]) >= 1L) {
      iv <- to_interval(lin, base$post_bins[[i]])
      out[[length(out) + 1L]] <- data.table(
        lineage = lin, lo = iv[1L], hi = iv[2L], side = "post", row = i
      )
    }
    if (length(out)) rbindlist(out) else NULL
  }), fill = TRUE)

  if (is.null(intervals) || !nrow(intervals)) {
    base[, c("pre_bins", "post_bins") := NULL]
    return(as.data.frame(base, stringsAsFactors = FALSE))
  }

  # An interval is identified by its (lo, hi) edges alone: the common grid is
  # shared by all lineages and one rescore returns rows for every lineage.
  interval_key <- unique(intervals[, .(lo, hi)])
  interval_key[, interval_id := seq_len(.N)]
  lookup <- rbindlist(lapply(seq_len(nrow(interval_key)), function(k) {
    lo <- interval_key$lo[k]
    hi <- interval_key$hi[k]
    if (!(hi > lo)) return(NULL)
    bm <- matrix(rep(c(0, lo, hi, 1), each = n_lin), nrow = n_lin)
    tb <- ComputeTemporalLRFixedBounds(
      expr, celltypes, lr_db, pseudotime, lineage_weights,
      bin_bounds = bm, hill_scale = hill_scale, n_bins = 3L,
      min_pct = min_pct, min_n_eff = min_n_eff,
      return_complete_grid = TRUE
    )
    st <- data.table::as.data.table(tb)[bin == 2L]
    st[, source_int := match(source, ct_levels)]
    st[, target_int := match(target, ct_levels)]
    data.table::data.table(
      interval_id = as.integer(interval_key$interval_id[k]),
      lineage = as.integer(st$lineage),
      pair_idx = as.integer(st$pair_idx),
      source = as.integer(st$source_int),
      target = as.integer(st$target_int),
      score = as.numeric(st$score)
    )
  }), fill = TRUE)
  intervals <- merge(intervals, interval_key, by = c("lo", "hi"))

  lookup_key <- function(iv_id, pair, src, tgt, lin) {
    hit <- lookup[interval_id == iv_id & pair_idx == pair & source == src &
                    target == tgt & lineage == lin, score]
    if (length(hit) == 1L) hit else NA_real_
  }
  pre_iv <- intervals[side == "pre", .(row, interval_id)]
  post_iv <- intervals[side == "post", .(row, interval_id)]

  refined <- rbindlist(lapply(merge_rows, function(i) {
    pair <- base$pair_idx[i]; src <- base$source[i]; tgt <- base$target[i]
    lin <- base$lineage[i]
    sisters <- smap$sibling_lineage[smap$focal_lineage == lin]
    pre_id <- pre_iv$interval_id[pre_iv$row == i]
    post_id <- post_iv$interval_id[post_iv$row == i]
    pre_id <- if (length(pre_id)) pre_id[1L] else NA_integer_
    post_id <- if (length(post_id)) post_id[1L] else NA_integer_

    focal_pre <- if (!is.na(pre_id)) lookup_key(pre_id, pair, src, tgt, lin) else NA_real_
    focal_post <- if (!is.na(post_id)) lookup_key(post_id, pair, src, tgt, lin) else NA_real_
    sister_pre <- if (!is.na(pre_id) && length(sisters)) {
      vapply(sisters, function(s) lookup_key(pre_id, pair, src, tgt, s), numeric(1L))
    } else numeric()
    sister_post <- if (!is.na(post_id) && length(sisters)) {
      vapply(sisters, function(s) lookup_key(post_id, pair, src, tgt, s), numeric(1L))
    } else numeric()

    pre_ok <- is.finite(focal_pre) && length(sister_pre) > 0L &&
      all(is.finite(sister_pre))
    post_ok <- is.finite(focal_post) && length(sister_post) > 0L &&
      all(is.finite(sister_post))
    post_diff <- if (post_ok) focal_post - mean(sister_post) else NA_real_
    theta <- if (pre_ok && post_ok) {
      (focal_post - focal_pre) - (mean(sister_post) - mean(sister_pre))
    } else NA_real_
    data.table(
      row = i,
      merged_pre = if (!is.na(pre_id)) {
        sprintf("%.4f|%.4f", intervals$lo[intervals$row == i & intervals$side == "pre"][1L],
                intervals$hi[intervals$row == i & intervals$side == "pre"][1L])
      } else NA_character_,
      merged_post = if (!is.na(post_id)) {
        sprintf("%.4f|%.4f", intervals$lo[intervals$row == i & intervals$side == "post"][1L],
                intervals$hi[intervals$row == i & intervals$side == "post"][1L])
      } else NA_character_,
      merged_pre_bins = if (!is.na(pre_id)) length(base$pre_bins[[i]]) else 0L,
      merged_post_bins = if (!is.na(post_id)) length(base$post_bins[[i]]) else 0L,
      post_fork_change = if (pre_ok && post_ok) focal_post - focal_pre else NA_real_,
      sibling_change = if (pre_ok && post_ok) mean(sister_post) - mean(sister_pre) else NA_real_,
      sibling_contrast = theta,
      post_branch_difference = post_diff,
      evaluable_merged = is.finite(theta)
    )
  }), fill = TRUE)

  base[refined$row, `:=`(
    merged_pre_edges = refined$merged_pre,
    merged_post_edges = refined$merged_post,
    merged_pre_bins = refined$merged_pre_bins,
    merged_post_bins = refined$merged_post_bins,
    post_branch_difference = refined$post_branch_difference
  )]
  base[refined$row, estimation_mode := data.table::fifelse(
    refined$evaluable_merged, "merged_interval", "not_evaluable"
  )]
  merged_rows <- refined$row[refined$evaluable_merged]
  if (length(merged_rows)) {
    base[merged_rows, `:=`(
      post_fork_change = refined$post_fork_change[refined$evaluable_merged],
      sibling_change = refined$sibling_change[refined$evaluable_merged],
      sibling_contrast = refined$sibling_contrast[refined$evaluable_merged],
      evaluable = TRUE,
      reason = "evaluable_merged_interval"
    )]
  }
  base[, c("pre_bins", "post_bins") := NULL]
  as.data.frame(base, stringsAsFactors = FALSE)
}

#' Combine geometric support with local branch assignment.
#'
#' Normalize the componentwise maximum of the two row-normalized inputs.
#' The result preserves their union of positive support, and retains at least
#' half of each input's mass. It is a branch propensity, not a calibrated
#' individual-cell fate posterior. Geometric memberships remain the reference
#' population weights for communication contrasts.
#'
#' @param reference_weights Cell-by-lineage geometric memberships.
#' @param assignment_weights Matching local branch memberships.
#' @return Cell-by-lineage normalized propensities in the original path order.
#' @export
CombineLineageMemberships <- function(reference_weights, assignment_weights) {
  reference <- .assert_lineage_weights(reference_weights, nrow(reference_weights))
  assignment <- .assert_lineage_weights(assignment_weights, nrow(reference))
  stopifnot(identical(dim(reference), dim(assignment)))
  weights <- pmax(reference, assignment)
  weights / rowSums(weights)
}
