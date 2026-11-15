# Fork-event-based contrast computation.
#
# Enumerates all fork events (event_scope = "all") or nearest-sibling units
# ("nearest_sibling") and scores each event child on shared pre/post windows
# with the fixed-boundary scorer `ComputeTemporalLRFixedBounds`,
# so complex activity uses the minimum subunit mean, percent-expressed is
# the subunit minimum, and the dataset-level Hill scale is reused unchanged.
#
# Event-child weights: each child column is the SUM of the original
# `lineage_weights` over that child's descendant terminal lineages. One extra
# remainder column carries the residual assignment mass so each cell row
# keeps the same total weight as the pooled call - the scorer row-normalizes,
# and a lone child column would collapse to all-1. The remainder preserves
# mass but is never a contrast candidate.
#
# Output rows carry the full event identity (`event_node`, `child`,
# `child_descendants`). `lineage` is only set for nearest-sibling mode, where
# the child IS a terminal lineage; for "all" it stays NA because a child is a
# subtree, not a terminal. Fork-bin activity on exact boundaries and
# post-fork persistence are undefined for this estimator and remain NA.

#' Compute fork-local contrasts using fork-event enumeration or nearest
#' siblings, with fixed-bin or exact-boundary windows.
#'
#' One scorer call per comparison unit evaluates all children on the same partition;
#' units only differ in which children are emitted and contrasted.
#'
#' @param expr Sparse expression matrix (genes x cells), log-normalized.
#' @param celltypes Integer cell type vector (one per cell).
#' @param lr_db LR database (from `.validate_lr_db`).
#' @param pseudotime Numeric pseudotime vector (one per cell).
#' @param lineage_weights Matrix (cells x lineages) of soft assignments.
#' @param lineage_paths List of root-to-terminal cluster paths.
#' @param sibling_map Output from `BuildSiblingContrastMap()` (for
#'   `event_scope = "nearest_sibling"`).
#' @param lineage_clusters Integer cluster labels per cell.
#' @param fork_bins Data frame with `lineage`, `fork_bin`,
#'   `fork_pseudotime`.
#' @param bin_bounds List of per-lineage bin edge vectors. Modes that need a
#'   single grid (`boundary = "exact"` or `event_scope = "all"`) require all
#'   rows to be identical (common grid).
#' @param event_scope `"all"` or `"nearest_sibling"`.
#' @param boundary `"exact"`, `"fixed_bin"`, or `"tree_segment_interval"`.
#' @param window_bins Bins on each side of the fork.
#' @param hill_scale Hill saturation scale.
#' @param min_n_eff Minimum effective cell count.
#' @param min_pct Minimum percent expressed.
#' @param activity_floor Minimum activity required for evaluation.
#' @param activity_transform `"hill"` retains the native activity;
#'   `"log_exposure"` uses log(1 + geometric mean / Hill scale) for contrasts.
#'   The activity floor remains on the original Hill scale.
#' @param prefork_weighting In native tree-segment intervals, `shared` gives
#'   every child the same observed parent population. `descendant` retains
#'   retrospective descendant-specific weights on the shared parent cells.
#' @param measurement_error Compute conditional delete-one-cell error of the
#'   native event contrast on its selected activity scale. Requires tree-segment intervals.
#' @param measurement_median_bounds Optional per-event, per-pair bounds for
#'   propagating cell measurement error through a fixed sample median.
#' @param bif_df Optional data frame with `lineage` and `b_star` columns
#'   (per-lineage bifurcation bins from the upstream detector). Used to
#'   anchor each event at the earliest of its children's own crossings.
#' @param anchor_rule How child crossings combine into the event anchor:
#'   `"min"` (default; earliest child crossing, the operational anchor) or
#'   `"mean"` (sensitivity comparison).
#' @param sibling_selection `"all"` averages every supported sister;
#'   `"time_matched"` uses the supported sister closest in post-fork
#'   pseudotime to the focal child.
#' @param units_override Optional `attr(result, "event_units")` of a pooled
#'   call. Reuses the same events, children, windows and boundary mode so a
#'   cell subset (e.g. one biological sample) is scored on identical objects.
#' @param descendant_window_children Child nodes whose post window includes
#'   their descendants beyond the next fork after a local transition split.
#' @return Data frame with one row per unit x emitted child x LR trajectory.
#'   `attr(out, "event_units")` stores the unit definitions for
#'   `units_override` reuse.
#' @export
ComputeForkContrastsEvent <- function(
    expr, celltypes, lr_db, pseudotime, lineage_weights,
    lineage_paths = NULL, sibling_map = NULL, lineage_clusters = NULL,
    fork_bins = NULL, bin_bounds = NULL,
    event_scope = c("all", "nearest_sibling"),
    boundary = c("exact", "fixed_bin", "tree_segment_interval"),
    window_bins = 2L, hill_scale = NULL,
    min_n_eff = 5.0, min_pct = 0.05, activity_floor = 1e-4,
    bif_df = NULL, units_override = NULL,
    descendant_window_children = integer(),
    anchor_rule = c("min", "mean", "child_mass_entropy"),
    scoring_mode = c("bin", "interval"),
    expression_mean_mode = c("all_cells", "positive_cells"),
    activity_mode = c("thresholded", "continuous"),
    sibling_selection = c("all", "time_matched", "program_matched"),
    sibling_program_profiles = NULL,
    activity_transform = c("hill", "log_exposure"),
    measurement_error = FALSE, measurement_median_bounds = NULL,
    prefork_weighting = c("shared", "descendant"), fork_anchors = NULL) {
  pt <- as.numeric(pseudotime)
  if (length(hill_scale) != 1L || !is.finite(hill_scale) || hill_scale <= 0) {
    stop("`hill_scale` is required (from compute_dataset_hill_s)", call. = FALSE)
  }
  lw <- as.matrix(lineage_weights)
  activity_transform <- match.arg(activity_transform)
  prefork_weighting <- match.arg(prefork_weighting)

  if (!is.null(units_override)) {
    event_scope <- units_override$event_scope
    boundary <- units_override$boundary
    scoring_mode <- if (!is.null(units_override$scoring_mode)) {
      units_override$scoring_mode
    } else {
      "bin"
    }
    expression_mean_mode <- if (!is.null(units_override$expression_mean_mode)) {
      units_override$expression_mean_mode
    } else {
      "all_cells"
    }
    activity_mode <- if (!is.null(units_override$activity_mode)) {
      units_override$activity_mode
    } else {
      "thresholded"
    }
    if (!is.null(units_override$activity_transform)) {
      activity_transform <- units_override$activity_transform
    }
    sibling_selection <- if (!is.null(units_override$sibling_selection)) {
      units_override$sibling_selection
    } else {
      "all"
    }
    if (!is.null(units_override$prefork_weighting)) {
      prefork_weighting <- units_override$prefork_weighting
    }
    units <- units_override$units
    if (!is.list(units)) {
      stop("`units_override` must be `attr(result, 'event_units')`", call. = FALSE)
    }
  } else {
    event_scope <- match.arg(event_scope)
    boundary <- match.arg(boundary)
    anchor_rule <- match.arg(anchor_rule)
    scoring_mode <- match.arg(scoring_mode)
    expression_mean_mode <- match.arg(expression_mean_mode)
    activity_mode <- match.arg(activity_mode)
    sibling_selection <- match.arg(sibling_selection)
    units <- .fork_event_units(
      event_scope,
      if (boundary == "tree_segment_interval") "fixed_bin" else boundary,
      lineage_paths, sibling_map, lineage_clusters,
      fork_bins, bin_bounds, pt, window_bins, lw, bif_df, anchor_rule
    )
  }
  if (is.null(units_override) && !is.null(fork_anchors)) {
    units <- lapply(units, function(u) {
      index <- match(u$event_node, fork_anchors$event_node)
      if (is.na(index) || !is.finite(fork_anchors$fork_time[index])) return(u)
      u$fork_pt <- fork_anchors$fork_time[index]
      u$b_star <- findInterval(u$fork_pt, u$grid, all.inside = TRUE)
      u$left_bound <- u$grid[u$b_star]
      u$right_bound <- u$grid[u$b_star + 1L]
      u$pre_lo <- u$grid[max(1L, u$b_star - window_bins)]
      u$pre_hi <- u$left_bound
      u$post_lo <- u$right_bound
      u$post_hi <- u$grid[min(length(u$grid), u$b_star + window_bins + 1L)]
      u$fork_time_source <- "observed_graph_branch_boundary"
      u
    })
  }
  scoring_floor <- switch(activity_transform,
    hill = activity_floor, log_exposure = -log1p(-activity_floor))
  if (!length(units)) return(data.frame())
  segment_mode <- boundary == "tree_segment_interval"
  if (measurement_error && !segment_mode) {
    stop("measurement error requires native tree segments", call. = FALSE)
  }
  if (segment_mode) {
    if (event_scope != "all" || scoring_mode != "interval" ||
        is.null(lineage_clusters)) {
      stop("tree segments require all events, interval scoring, and cell nodes",
           call. = FALSE)
    }
    if (is.null(units_override)) {
      paths <- data.table::rbindlist(lapply(seq_along(lineage_paths), function(i) {
        p <- as.integer(lineage_paths[[i]])
        data.table::data.table(lineage = i, position = seq_along(p), node = p)
      }))
      paths[, parent := data.table::shift(node), by = lineage]
      edges <- unique(paths[!is.na(parent), .(parent, child = node)])
      forks <- edges[, .N, by = parent][N >= 2L, parent]
      units <- lapply(units, function(u) {
        event_path <- paths[node == u$event_node][1L]
        prefix <- paths[lineage == event_path$lineage &
                          position <= event_path$position]
        previous_fork <- max(c(0L, prefix[
          node %in% forks & position < event_path$position, position]))
        u$segment_pre_nodes <- prefix[position > previous_fork, node]
        u$segment_post_nodes <- lapply(u$children, function(child) {
          child_positions <- paths[node == child,
                                   .(lineage, child_pos = position)]
          portion <- merge(paths, child_positions, by = "lineage",
                           allow.cartesian = TRUE)[position >= child_pos]
          portion[, end := min(c(position[
            node %in% forks & !(child %in% descendant_window_children)],
            max(position))),
                  by = lineage]
          unique(portion[position <= end, node])
        })
        if (sibling_selection == "time_matched") {
          post_nodes <- data.table::as.data.table(stack(setNames(
            u$segment_post_nodes, seq_along(u$children))))
          data.table::setnames(post_nodes, c("values", "ind"),
                               c("node", "child_index"))
          post_nodes[, child_index := as.integer(as.character(child_index))]
          cell_times <- data.table::data.table(
            node = as.integer(lineage_clusters), time = pt)
          post_times <- merge(cell_times, post_nodes, by = "node",
                              allow.cartesian = TRUE)[
                                , .(time = median(time[is.finite(time)])),
                                by = child_index]
          u$child_post_time <- post_times$time[
            match(seq_along(u$children), post_times$child_index)]
        }
        u
      })
    }
  }

  switch(sibling_selection, program_matched = {
    if (is.null(units_override)) {
      if (is.null(sibling_program_profiles)) {
        group <- data.table::data.table(node=as.integer(lineage_clusters),target=celltypes)
        group[,profile:=.GRP,by=.(node,target)]
        profiles <- unique(group)
        data.table::setorder(profiles,profile)
        count <- tabulate(group$profile)
        membership <- Matrix::sparseMatrix(i=group$profile,j=seq_len(nrow(group)),x=1,
          dims=c(nrow(profiles),nrow(group)))
        program_db <- .validate_lr_db(lr_db)
        non_lr <- setdiff(seq_len(nrow(expr)),unique(unlist(c(program_db$ligand_rows,program_db$receptor_rows))))
        values <- as.matrix(membership%*%Matrix::t(expr[non_lr,,drop=FALSE]))/count
        sibling_program_profiles <- lapply(split(seq_len(nrow(profiles)),profiles$target),function(index) {
          value <- values[index,,drop=FALSE]
          rownames(value) <- as.character(profiles$node[index])
          attr(value,"cell_count") <- count[index]
          value
        })
      }
      profile_list <- if (is.matrix(sibling_program_profiles)) list(all=sibling_program_profiles) else sibling_program_profiles
      program_paths <- data.table::as.data.table(stack(setNames(lineage_paths,seq_along(lineage_paths))))
      data.table::setnames(program_paths,c("values","ind"),c("node","lineage"))
      program_paths[,lineage:=as.integer(as.character(lineage))]
      program_paths[,position:=seq_len(.N),by=lineage]
      units <- lapply(units,function(u) {
        post_nodes <- switch(as.character(segment_mode),"TRUE"={
          nodes <- data.table::as.data.table(stack(setNames(u$segment_post_nodes,seq_along(u$children))))
          data.table::setnames(nodes,c("values","ind"),c("node","child_index"))
          nodes[,child_index:=as.integer(as.character(child_index))]
          nodes
        },"FALSE"={
          starts <- program_paths[node%in%u$children,.(lineage,child=node,child_pos=position)]
          unique(merge(program_paths,starts,by="lineage",allow.cartesian=TRUE)[position>=child_pos,
            .(node,child_index=match(child,u$children))])
        })
        pre_nodes <- if (segment_mode) u$segment_pre_nodes else u$event_node
        u$child_program_overlap <- lapply(profile_list,function(profile) {
          sizes <- attr(profile,"cell_count")
          if (is.null(sizes)) sizes <- rep.int(1,nrow(profile))
          parent <- match(as.character(pre_nodes),rownames(profile))
          parent <- parent[!is.na(parent)]
          baseline <- colSums(profile[parent,,drop=FALSE]*sizes[parent])/sum(sizes[parent])
          post <- data.table::copy(post_nodes)
          post[,profile_index:=match(as.character(node),rownames(profile))]
          post <- post[!is.na(profile_index)]
          weights <- Matrix::sparseMatrix(i=post$child_index,j=post$profile_index,
            x=sizes[post$profile_index],dims=c(length(u$children),nrow(profile)))
          child_profiles <- as.matrix(weights%*%profile)/as.numeric(Matrix::rowSums(weights))
          delta <- sweep(child_profiles,2L,baseline,"-")
          change <- cbind(pmax(delta,0),pmax(-delta,0))
          pairs <- data.table::CJ(focal=seq_along(u$children),sibling=seq_along(u$children))
          overlap <- rowSums(pmin(change[pairs$focal,,drop=FALSE],change[pairs$sibling,,drop=FALSE])^2)
          post_count <- as.numeric(Matrix::rowSums(weights))
          precision <- 1/(1/post_count[pairs$focal]+1/post_count[pairs$sibling]+2/sum(sizes[parent]))
          overlap <- overlap*precision
          matrix(overlap,nrow=length(u$children),byrow=TRUE)
        })
        u
      })
    }
  })

  # --- Score units ---------------------------------------------------------
  # Per-bin scoring on the unit's own grid (plus the exact boundary for D1),
  # then the shared-window aggregation: each side averages the
  # `window_bins` bins next to the fork over bins where the focal child and at
  # least one sibling child are both supported. Units sharing one scoring
  # matrix and one partition collapse into a single fixed-bounds call:
  # nearest-sibling units always share the full weight matrix and group by
  # window+grid signature; event units carry their own child-descendant
  # matrix and score individually.
  job_key <- vapply(seq_along(units), function(ui) {
    u <- units[[ui]]
    wsig <- paste(formatC(c(u$pre_lo, u$pre_hi, u$post_lo, u$post_hi),
                          digits = 12L, format = "g"), collapse = ",")
    gsig <- paste(formatC(u$grid, digits = 12L, format = "g"), collapse = ",")
    if (event_scope == "all") {
      pref <- if (scoring_mode == "interval") "u_int" else "u"
      paste0(pref, ui, "_", activity_mode)
    } else {
      paste(wsig, gsig, scoring_mode, activity_mode, sep = "|")
    }
  }, character(1L))

  contrasts <- data.table::rbindlist(lapply(unique(job_key), function(jk) {
    ui_set <- which(job_key == jk)
    u0 <- units[[ui_set[1L]]]
    w_score <- if (event_scope == "all") {
      .event_child_weights(u0, lw)
    } else {
      lw
    }
    switch(as.character(!is.null(units_override) && !segment_mode &&
                        event_scope == "all" && !is.null(lineage_paths) &&
                        !is.null(lineage_clusters)), "TRUE" = {
      post_cells <- is.finite(pt) & pt >= u0$post_lo & pt <= u0$post_hi
      child_nodes <- lapply(seq_along(u0$children), function(k) {
        unique(unlist(lapply(lineage_paths[u0$desc[[k]]], function(path) {
          path[seq.int(match(u0$children[k], path), length(path))]
        }), use.names = FALSE))
      })
      allowed <- vapply(child_nodes, function(nodes) {
        !post_cells | lineage_clusters %in% nodes
      }, logical(length(pt)))
      child_weights <- w_score[, seq_along(u0$children), drop = FALSE] * allowed
      w_score <- cbind(child_weights,
                       pmax(0, rowSums(lw) - rowSums(child_weights)))
    })
    if (is.null(u0$grid)) {
      stop("event unit lacks `grid`; rebuild units with this version",
           call. = FALSE)
    }
    bounds <- if (scoring_mode == "interval") {
      sort(unique(c(0, u0$pre_lo, u0$pre_hi, u0$post_lo, u0$post_hi, 1)))
    } else {
      sort(unique(c(
        u0$grid, u0$pre_lo, u0$pre_hi, u0$post_lo, u0$post_hi
      )))
    }
    bins <- .locate_window_bins(bounds, u0)
    score_pt <- pt
    if (segment_mode) {
      w_score <- .event_segment_weights(u0, lw, lineage_clusters, prefork_weighting)
      score_pt <- rep(0.5, length(pt))
      score_pt[lineage_clusters %in% u0$segment_pre_nodes] <- 1 / 6
      score_pt[lineage_clusters %in% unlist(u0$segment_post_nodes)] <- 5 / 6
      bounds <- seq(0, 1, length.out = 4L)
      bins <- list(pre = 1L, fork = 2L, post = 3L)
    }
    tb <- ComputeTemporalLRFixedBounds(
      expr, celltypes, lr_db, score_pt, w_score,
      bin_bounds = replicate(ncol(w_score), bounds, simplify = FALSE),
      hill_scale = hill_scale, min_pct = min_pct, min_n_eff = min_n_eff,
      return_complete_grid = TRUE
    )
    sdt <- data.table::as.data.table(tb)
    if (expression_mean_mode == "positive_cells") {
      sdt[, score := score_positive]
    }
    if (activity_mode == "continuous" &&
        expression_mean_mode == "all_cells") {
      tb_cont <- data.table::as.data.table(ComputeTemporalLRFixedBounds(
        expr, celltypes, lr_db, score_pt, w_score,
        bin_bounds = replicate(ncol(w_score), bounds, simplify = FALSE),
        hill_scale = hill_scale, min_pct = 0.0, min_n_eff = min_n_eff,
        return_complete_grid = TRUE
      ))
      sdt <- merge(
        sdt[, !"score"],
        tb_cont[, .(pair_idx, source, target, lineage, bin, score)],
        by = c("pair_idx", "source", "target", "lineage", "bin"),
        all.x = TRUE, sort = FALSE
      )
    }
    if (activity_transform == "log_exposure") {
      sdt[, score := -log1p(-score)]
    }
    sdt[, row_valid := is.finite(score) & is.finite(n_eff_cells) &
          n_eff_cells >= min_n_eff & structural_valid %in% TRUE]
    sdt <- sdt[bin %in% c(bins$pre, bins$post, bins$fork)]
    sdt[, comp := data.table::fcase(
      bin %in% bins$pre, "pre",
      bin %in% bins$post, "post",
      bin %in% bins$fork, "fork"
    )]
    sdt[, cap := stats::quantile(
      n_eff_cells[row_valid], 0.95, na.rm = TRUE, names = FALSE
    ), by = c("pair_idx", "source", "target", "lineage")]
    sdt[, w_eff := pmin(n_eff_cells, cap)]
    req_bins <- if (scoring_mode == "interval") 1L else window_bins
    data.table::rbindlist(lapply(ui_set, function(ui) {
      .event_unit_contrasts(
        sdt, units[[ui]], ui, event_scope, req_bins, scoring_floor,
        bins = bins, bounds = bounds, activity_mode = activity_mode,
        sibling_selection = sibling_selection,
        measurement = if (measurement_error) list(expr = expr, lr_db = lr_db,
          weights = w_score / rowSums(w_score), celltypes = celltypes,
          cell_bin = findInterval(score_pt, bounds, rightmost.closed = TRUE, all.inside = TRUE),
          hill_scale = hill_scale, expression_mean_mode = expression_mean_mode,
          activity_transform = activity_transform,
          median_bounds = if (is.null(measurement_median_bounds)) NULL else
            measurement_median_bounds[measurement_median_bounds$unit_id == ui, , drop = FALSE]) else NULL
      )
    }), fill = TRUE)
  }), fill = TRUE)

  if (!nrow(contrasts)) return(data.frame())
  if (segment_mode) {
    contrasts[, c("pre_lo", "pre_hi", "post_lo", "post_hi") :=
                rep(list(NA_real_), 4L)]
    contrasts[, temporal_coordinate := "tree_state_ordinal_stages"]
  }

  data.table::setcolorder(contrasts, c(
    "unit_id", "pair_idx", "source", "target", "lineage", "fork_cluster",
    "event_node", "child", "child_descendants", "fork_time", "n_children",
    "pre_lo", "pre_hi", "post_lo", "post_hi",
    "focal_pre", "focal_post", "other_pre", "other_post",
    "pre_activity", "fork_activity", "post_activity",
    "post_fork_change", "sibling_change", "theta", "sibling_contrast",
    "activity_strength", "pre_fork_enrichment", "fork_locality",
    "post_fork_persistence", "n_pre_bins", "n_post_bins",
    "n_pre_bins_shared", "n_post_bins_shared", "n_sisters_mapped",
    "effective_weight", "evaluable", "reason",
    "b_star", "fork_time_source", "temporal_coordinate", "n_eff_source",
    "inference_unit"
  ))
  if (sibling_selection=="program_matched" && is.null(units_override)) {
    selected <- unique(contrasts[!is.na(sibling_children)&nzchar(sibling_children),
      .(unit_id,source,target,child,peer=as.integer(sibling_children))])
    units <- lapply(seq_along(units),function(index) {
      u <- units[[index]]
      u$program_sibling_choices <- selected[unit_id==index,.(source,target,child,peer)]
      u
    })
  }
  out <- as.data.frame(contrasts, stringsAsFactors = FALSE)
  attr(out, "event_units") <- list(
    event_scope = event_scope, boundary = boundary,
    scoring_mode = scoring_mode, expression_mean_mode = expression_mean_mode,
    activity_mode = activity_mode, activity_transform = activity_transform,
    sibling_selection = sibling_selection, prefork_weighting = prefork_weighting,
    units = units
  )
  out
}

#' Build the scoring weight matrix for one "all"-scope event unit: one column
#' per child holding the summed descendant-terminal weights, plus a remainder
#' column preserving residual mass (never a candidate).
#' @keywords internal
.event_child_weights <- function(u, lw) {
  child_mat <- vapply(seq_along(u$children), function(k) {
    d <- u$desc[[k]]
    if (length(d)) rowSums(lw[, d, drop = FALSE]) else numeric(nrow(lw))
  }, numeric(nrow(lw)))
  remainder <- pmax(0, rowSums(lw) - rowSums(child_mat))
  cbind(child_mat, remainder)
}

#' Locate the pre/post/fork bin sets inside a unit's partition. `bounds`
#' spans [0, 1] and contains all four window edges, so pre bins are the
#' partition bins inside [pre_lo, pre_hi], post bins inside [post_lo,
#' post_hi], and fork bins inside [pre_hi, post_lo] (fixed-bin mode only).
#' @keywords internal
.locate_window_bins <- function(bounds, u) {
  lo <- bounds[-length(bounds)]
  hi <- bounds[-1L]
  list(
    pre = which(lo >= u$pre_lo & hi <= u$pre_hi),
    post = which(lo >= u$post_lo & hi <= u$post_hi),
    fork = which(lo >= u$pre_hi & hi <= u$post_lo)
  )
}

#' Mean of a vector, or NA when empty.
#' @keywords internal
.mean_or_na <- function(x) if (length(x)) mean(x) else NA_real_

#' Build contrast rows for one unit from its scored partition table.
#' `sdt` contains one row per (pair, source, target, scoring column, comp,
#' bin). A sister enters the fixed panel only when it has pre and post
#' support. For each panel sister, focal and sister scores are then joined on
#' the bins valid for both, their paired pre/post changes are contrasted, and
#' those sister-level contrasts are averaged with equal weight. Evaluability
#' keeps the existing gates: `window_bins` supported bins on each side,
#' `window_bins` shared bins on each side, activity above floor.
#' @keywords internal
.event_unit_contrasts <- function(sdt, u, unit_id, event_scope, window_bins,
                                  activity_floor, bins = NULL,
                                  bounds = NULL, activity_mode = "thresholded",
                                  sibling_selection = "all", measurement = NULL) {
  kids <- u$children
  col_of <- if (event_scope == "all") seq_along(kids) else kids
  emit <- if (event_scope == "all") seq_along(kids) else 1L
  key <- c("pair_idx", "source", "target")
  out <- data.table::rbindlist(lapply(emit, function(k) {
    post_bins_k <- if (is.null(bins)) NULL else bins$post
    f <- sdt[lineage == col_of[k] &
               (comp %in% c("pre", "fork") |
                  (comp == "post" &
                     (is.null(post_bins_k) | bin %in% post_bins_k)))]
    if (!nrow(f)) return(NULL)
    o <- sdt[lineage %in% col_of[setdiff(seq_along(kids), k)] &
               (comp == "pre" |
                  (comp == "post" &
                     (is.null(post_bins_k) | bin %in% post_bins_k))) &
               row_valid %in% TRUE]
    # The control panel must be identical on both sides of the fork: a
    # sibling enters only when it is supported in BOTH the pre and the
    # post window for that contrast object. Without this, swapping which
    # sister drives each side fabricates a contrast from unchanged focal
    # activity.
    sib_ok <- o[, .(has_pre = any(comp == "pre"),
                    has_post = any(comp == "post")),
                by = c(key, "lineage")]
    panel <- sib_ok[has_pre & has_post, c(key, "lineage"), with = FALSE]
    if (sibling_selection == "time_matched" && nrow(panel)) {
      panel[, age_distance := abs(
        u$child_post_time[match(lineage, col_of)] - u$child_post_time[k])]
      data.table::setorderv(panel, c(key, "age_distance", "lineage"))
      panel <- panel[, .SD[1L], by = key][, c(key, "lineage"), with = FALSE]
    }
    if (sibling_selection == "program_matched" && nrow(panel)) {
      if (is.null(u$program_sibling_choices)) {
        choices <- unique(panel[,.(source,target,lineage)])
        choices[,program_overlap:={
          overlap <- u$child_program_overlap
          if (is.list(overlap)) overlap <- overlap[[if ("all"%in%names(overlap)) "all" else as.character(target[1L])]]
          overlap[k,match(lineage,col_of)]
        },by=target]
        choices <- choices[is.finite(program_overlap)]
        data.table::setorder(choices,source,target,-program_overlap,lineage)
        choices <- choices[,.SD[1L],by=.(source,target)][,.(source,target,lineage)]
      } else {
        choices <- u$program_sibling_choices[child==kids[k],.(source,target,lineage=col_of[match(peer,kids)])]
      }
      panel <- merge(panel,choices,by=c("source","target","lineage"),sort=FALSE)
    }
    o <- merge(o, panel, by = c(key, "lineage"), sort = FALSE)
    # Align the focal and each sister on their own jointly valid bins before
    # taking either branch's pre/post change. This prevents a common temporal
    # shift from becoming branch-specific when one branch loses a bin.
    f_pair <- f[comp %in% c("pre", "post") & row_valid %in% TRUE,
                c(key, "comp", "bin", "score"), with = FALSE]
    data.table::setnames(f_pair, "score", "focal_score")
    s_pair <- o[, c(key, "lineage", "comp", "bin", "score"), with = FALSE]
    paired <- merge(
      f_pair, s_pair, by = c(key, "comp", "bin"),
      all = FALSE, sort = FALSE, allow.cartesian = TRUE
    )
    pair_stats <- paired[, .(
      focal_pre_pair = .mean_or_na(focal_score[comp == "pre"]),
      focal_post_pair = .mean_or_na(focal_score[comp == "post"]),
      sister_pre_pair = .mean_or_na(score[comp == "pre"]),
      sister_post_pair = .mean_or_na(score[comp == "post"])
    ), by = c(key, "lineage")]
    pair_stats[, `:=`(
      focal_delta_pair = focal_post_pair - focal_pre_pair,
      sister_delta_pair = sister_post_pair - sister_pre_pair
    )]
    pair_stats[, contrast_pair := focal_delta_pair - sister_delta_pair]
    pair_stats <- merge(
      panel, pair_stats, by = c(key, "lineage"), all.x = TRUE,
      sort = FALSE
    )
    sdiff <- pair_stats[, .(
      paired_pre_focal = mean(focal_pre_pair),
      paired_post_focal = mean(focal_post_pair),
      paired_pre_sister = mean(sister_pre_pair),
      paired_post_sister = mean(sister_post_pair),
      focal_delta = mean(focal_delta_pair),
      sib_delta = mean(sister_delta_pair),
      paired_contrast = mean(contrast_pair),
      n_sisters_mapped = .N,
      sibling_children = paste(sort(kids[match(lineage, col_of)]), collapse = ",")
    ), by = key]
    ob <- o[, .(n_sis = .N), by = c(key, "comp", "bin")]
    fm <- merge(f, ob, by = c(key, "comp", "bin"), all.x = TRUE, sort = FALSE)
    fm[is.na(n_sis), n_sis := 0L]
    fm[, shared := row_valid %in% TRUE & n_sis > 0L]
    a <- fm[, .(
      n_pre_bins = sum(row_valid & comp == "pre"),
      n_post_bins = sum(row_valid & comp == "post"),
      n_pre_bins_shared = sum(shared & comp == "pre"),
      n_post_bins_shared = sum(shared & comp == "post"),
      pre_activity = .mean_or_na(score[row_valid & comp == "pre"]),
      fork_activity = .mean_or_na(score[row_valid & comp == "fork"]),
      post_activity = .mean_or_na(score[row_valid & comp == "post"]),
      has_expr_support = if ("expression_supported" %in% names(fm)) {
        any(expression_supported %in% TRUE & row_valid %in% TRUE & comp %in% c("pre", "post"))
      } else {
        TRUE
      },
      post_sd = if (sum(row_valid & comp == "post") > 1L) {
        stats::sd(score[row_valid & comp == "post"])
      } else {
        NA_real_
      },
      pre_f = .mean_or_na(score[shared & comp == "pre"]),
      post_f = .mean_or_na(score[shared & comp == "post"]),
      effective_weight = sum(
        data.table::fifelse(row_valid, w_eff, 0), na.rm = TRUE
      )
    ), by = key]
    a <- merge(a, sdiff, by = key, all.x = TRUE, sort = FALSE)
    a[, `:=`(
      pre_f = paired_pre_focal,
      post_f = paired_post_focal,
      pre_s = paired_pre_sister,
      post_s = paired_post_sister,
      post_fork_change = focal_delta,
      sibling_change = sib_delta,
      sibling_contrast = paired_contrast
    )]
    a[, activity_strength := pmax(
      pre_activity, fork_activity, post_activity, na.rm = TRUE
    )]
    a[!is.finite(activity_strength), activity_strength := NA_real_]
    a[, reason := data.table::fcase(
      n_pre_bins < window_bins, "insufficient_pre_bins",
      n_post_bins < window_bins, "insufficient_post_bins",
      n_pre_bins_shared < window_bins, "insufficient_shared_pre_bins",
      n_post_bins_shared < window_bins, "insufficient_shared_post_bins",
      activity_mode == "continuous" & has_expr_support %in% FALSE,
      "below_expression_threshold",
      !is.finite(activity_strength) | activity_strength < activity_floor,
      "below_activity_floor",
      default = "evaluable"
    )]

    # A panel sister that shares zero valid bins with the focal contributes
    # all-NA paired stats; mark the row non-evaluable instead of emitting
    # evaluable = TRUE with a missing contrast.
    a[is.finite(sibling_contrast) == FALSE & reason == "evaluable",
      reason := "no_shared_bins_with_panel"]
    a[, evaluable := reason == "evaluable"]
    a[evaluable == FALSE, `:=`(
      post_fork_change = NA_real_, sibling_change = NA_real_,
      sibling_contrast = NA_real_
    )]
    desc_label <- if (event_scope == "all") {
      paste(u$desc[[k]], collapse = ",")
    } else {
      as.character(kids[k])
    }
    data.table::data.table(
      unit_id = unit_id,
      pair_idx = a$pair_idx, source = a$source, target = a$target,
      lineage = if (event_scope == "all") NA_integer_ else as.integer(kids[k]),
      fork_cluster = u$event_node,
      event_node = u$event_node,
      child = kids[k],
      child_descendants = desc_label,
      child_deep_fork = if (event_scope == "all" && !is.null(u$deep_fork)) {
        u$deep_fork[k]
      } else {
        NA_real_
      },
      fork_time = u$fork_pt,
      n_children = length(kids),
      pre_lo = u$pre_lo, pre_hi = u$pre_hi,
      post_lo = u$post_lo, post_hi = u$post_hi,
      focal_pre = a$pre_f, focal_post = a$post_f,
      other_pre = a$pre_s, other_post = a$post_s,
      pre_activity = a$pre_activity, fork_activity = a$fork_activity,
      post_activity = a$post_activity,
      post_fork_change = a$post_fork_change,
      sibling_change = a$sibling_change,
      theta = a$sibling_contrast,
      sibling_contrast = a$sibling_contrast,
      activity_strength = a$activity_strength,
      pre_fork_enrichment = a$pre_activity /
        (a$pre_activity + a$post_activity + 1e-12),
      fork_locality = a$fork_activity /
        (a$pre_activity + a$fork_activity + a$post_activity + 1e-12),
      post_fork_persistence = data.table::fifelse(
        is.finite(a$post_sd),
        1 - a$post_sd / (abs(a$post_activity) + a$post_sd + 1e-12),
        NA_real_
      ),
      n_pre_bins = a$n_pre_bins, n_post_bins = a$n_post_bins,
      n_pre_bins_shared = a$n_pre_bins_shared,
      n_post_bins_shared = a$n_post_bins_shared,
      n_sisters_mapped = a$n_sisters_mapped,
      sibling_children = a$sibling_children,
      effective_weight = a$effective_weight,
      evaluable = a$evaluable,
      reason = a$reason,
      b_star = if (is.null(u$b_star)) NA_integer_ else as.integer(u$b_star),
      fork_time_source = if (is.null(u$fork_time_source)) {
        "topology"
      } else {
        u$fork_time_source
      },
      temporal_coordinate = "common_pseudotime_grid",
      n_eff_source = "reported",
      inference_unit = "pooled_cells_descriptive"
    )
  }), fill = TRUE)
  if (is.null(measurement)) return(out)
  if (!any(out$evaluable %in% TRUE)) {
    out[, `:=`(theta_measurement_se = NA_real_, median_measurement_variance = NA_real_,
      theta_recomputed = NA_real_, measurement_score = NA_real_)]
    return(out)
  }
  # All deletion calculations condition on the saved geometry, support panel,
  # and expression gate. Shared cells and complex-subunit switches are retained.
  error <- out[evaluable %in% TRUE, {
    current <- data.table::copy(.SD)
    source_id <- .BY$source
    target_id <- .BY$target
    pair_ids <- sort(unique(current$pair_idx))
    scored <- sdt[source == source_id & target == target_id & pair_idx %in% pair_ids]
    active_children <- unique(c(current$child,
      as.integer(unlist(strsplit(current$sibling_children, ",", fixed = TRUE)))))
    contexts <- unique(scored[comp %in% c("pre", "post"), .(lineage, bin, comp)])
    contexts[, child := kids[match(lineage, col_of)]]
    contexts <- contexts[child %in% active_children]
    if (is.null(measurement$source_weights)) {
      weights <- measurement$weights[, contexts$lineage, drop = FALSE] *
        outer(measurement$cell_bin, contexts$bin, "==")
      weights[weights < 1e-8] <- 0
      source_weights <- weights * (measurement$celltypes == source_id)
      receiver_weights <- weights * (measurement$celltypes == target_id)
    } else {
      context_index <- match(paste(contexts$lineage, contexts$bin),
        paste(measurement$contexts$lineage, measurement$contexts$bin))
      source_weights <- measurement$source_weights[, context_index, drop = FALSE]
      receiver_weights <- measurement$receiver_weights[, context_index, drop = FALSE]
    }
    source_member <- Matrix::rowSums(source_weights) > 0
    receiver_member <- Matrix::rowSums(receiver_weights) > 0
    cells <- which(source_member | receiver_member)
    populations <- data.table::data.table(cell_column = seq_along(cells),
      population = as.integer(source_member[cells]) + 2L * as.integer(receiver_member[cells]))
    source_weights <- as.matrix(source_weights[cells, , drop = FALSE])
    receiver_weights <- as.matrix(receiver_weights[cells, , drop = FALSE])
    ligand_rows <- measurement$lr_db$ligand_rows[pair_ids]
    receptor_rows <- measurement$lr_db$receptor_rows[pair_ids]
    ligand_genes <- sort(unique(unlist(ligand_rows)))
    receptor_genes <- sort(unique(unlist(receptor_rows)))
    ligand_rows <- lapply(ligand_rows, match, table = ligand_genes)
    receptor_rows <- lapply(receptor_rows, match, table = receptor_genes)
    ligand_values <- as.matrix(measurement$expr[ligand_genes, cells, drop = FALSE])
    receptor_values <- as.matrix(measurement$expr[receptor_genes, cells, drop = FALSE])
    ligand_first <- vapply(ligand_rows, `[`, integer(1L), 1L)
    receptor_first <- vapply(receptor_rows, `[`, integer(1L), 1L)
    ligand_complex <- which(lengths(ligand_rows) > 1L)
    receptor_complex <- which(lengths(receptor_rows) > 1L)
    contexts[, context_column := .I]
    deleted <- contexts[, {
      sw <- source_weights[, context_column]
      rw <- receiver_weights[, context_column]
      l_sum <- as.numeric(ligand_values %*% sw)
      r_sum <- as.numeric(receptor_values %*% rw)
      l_mass <- rep(sum(sw), nrow(ligand_values))
      r_mass <- rep(sum(rw), nrow(receptor_values))
      l_den <- l_mass - rep(sw, each = nrow(ligand_values))
      r_den <- r_mass - rep(rw, each = nrow(receptor_values))
      if (measurement$expression_mean_mode == "positive_cells") {
        l_mass <- as.numeric((ligand_values > 0) %*% sw)
        r_mass <- as.numeric((receptor_values > 0) %*% rw)
        l_den <- l_mass - sweep(ligand_values > 0, 2L, sw, "*")
        r_den <- r_mass - sweep(receptor_values > 0, 2L, rw, "*")
      }
      l_mean <- l_sum / pmax(l_mass, .Machine$double.eps)
      r_mean <- r_sum / pmax(r_mass, .Machine$double.eps)
      l_deleted <- (l_sum - sweep(ligand_values, 2L, sw, "*")) / pmax(l_den, .Machine$double.eps)
      r_deleted <- (r_sum - sweep(receptor_values, 2L, rw, "*")) / pmax(r_den, .Machine$double.eps)
      l_deleted[l_den <= 0] <- 0
      r_deleted[r_den <= 0] <- 0
      ligand <- l_deleted[ligand_first, , drop = FALSE]
      receptor <- r_deleted[receptor_first, , drop = FALSE]
      ligand[ligand_complex, ] <- t(vapply(ligand_rows[ligand_complex], function(i)
        matrixStats::colMins(l_deleted[i, , drop = FALSE]), numeric(length(cells))))
      receptor[receptor_complex, ] <- t(vapply(receptor_rows[receptor_complex], function(i)
        matrixStats::colMins(r_deleted[i, , drop = FALSE]), numeric(length(cells))))
      geo <- sqrt(pmax(ligand, 0) * pmax(receptor, 0))
      activity <- geo / (measurement$hill_scale + geo)
      full_geo <- sqrt(vapply(ligand_rows, function(i) min(l_mean[i]), numeric(1L)) *
        vapply(receptor_rows, function(i) min(r_mean[i]), numeric(1L)))
      full_activity <- full_geo / (measurement$hill_scale + full_geo)
      switch(measurement$activity_transform, log_exposure = {
        activity <- log1p(geo / measurement$hill_scale)
        full_activity <- log1p(full_geo / measurement$hill_scale)
      })
      if (measurement$expression_mean_mode == "positive_cells" || activity_mode != "continuous") {
        fixed_gate <- scored[lineage == .BY$lineage & bin == .BY$bin][match(pair_ids, pair_idx),
          expression_supported %in% TRUE]
        activity[!fixed_gate, ] <- 0
        full_activity[!fixed_gate] <- 0
      }
      .(deleted_activity = list(activity), full_activity = list(full_activity))
    }, by = .(child, lineage, bin, comp, context_column)]
    changes <- deleted[, .(
      deleted_change = list(deleted_activity[comp == "post"][[1L]] - deleted_activity[comp == "pre"][[1L]]),
      full_change = list(full_activity[comp == "post"][[1L]] - full_activity[comp == "pre"][[1L]])), by = child]
    changes <- as.list(changes[, .(child, deleted_change, full_change)])
    current[, {
      sisters <- as.integer(strsplit(.BY$sibling_children, ",", fixed = TRUE)[[1L]])
      focal <- match(.BY$child, changes$child)
      other <- match(sisters, changes$child)
      deleted_theta <- changes$deleted_change[[focal]] - Reduce("+", changes$deleted_change[other]) / length(other)
      full_theta <- changes$full_change[[focal]] - Reduce("+", changes$full_change[other]) / length(other)
      row_index <- match(pair_idx, pair_ids)
      stopifnot(max(abs(full_theta[row_index] - theta)) < 1e-5)
      variance <- populations[, {
        values <- deleted_theta[row_index, cell_column, drop = FALSE]
        .(variance = list((.N - 1) / .N * rowSums((values - rowMeans(values))^2)))
      }, by = population]
      theta_variance <- Reduce("+", variance$variance)
      median_variance <- rep(NA_real_, length(row_index))
      if (!is.null(measurement$median_bounds)) {
        median_bounds <- data.table::as.data.table(measurement$median_bounds)[
          source == source_id & target == target_id & child == .BY$child]
        median_index <- match(pair_idx, median_bounds$pair_idx)
        values <- deleted_theta[row_index, , drop = FALSE]
        values <- pmin(pmax(values, median_bounds$lower[median_index]), median_bounds$upper[median_index]) *
          median_bounds$scale[median_index]
        median_parts <- populations[, {
          component <- values[, cell_column, drop = FALSE]
          .(variance = list((.N - 1) / .N * rowSums((component - rowMeans(component))^2)))
        }, by = population]
        median_variance <- Reduce("+", median_parts$variance)
      }
      .(pair_idx, theta_measurement_se = sqrt(theta_variance),
        median_measurement_variance = median_variance,
        theta_recomputed = full_theta[row_index])
    }, by = .(child, sibling_children)]
  }, by = .(source, target)]
  output_columns <- names(out)
  out <- merge(out, error[, !"sibling_children"], by = c("pair_idx", "source", "target", "child"), all.x = TRUE)
  data.table::setcolorder(out, c(output_columns, setdiff(names(out), output_columns)))
  out[, measurement_score := data.table::fifelse(evaluable %in% TRUE & theta_measurement_se > 0,
    theta / theta_measurement_se, NA_real_)]
  out
}

#' Enumerate comparison units and attach fork times and windows.
#' Each unit carries: event_node, children (emitted child ids), desc
#' (descendant terminal indices per child, "all" scope only), lineage
#' (focal lineage, nearest scope only), fork_pt, and the four window edges.
#' The event fork time is the earliest of the children's own bifurcation
#' points: each child contributes the `b_star` edge of its descendant
#' terminal lineages from `bif_df` (the per-lineage crossing detected
#' upstream). Anchoring on the ancestor cluster's median pseudotime can
#' land the shared windows past a short branch's support end, and a
#' weight-quantile onset is swamped by diffuse early mass; the detector's
#' crossings locate the actual separation of each child, and the earliest
#' one marks where the branches still share a trajectory: the mean places
#' the pre-window after an already-separated child. Falls back to the
#' ancestor median / `fork_bins` when `bif_df` is absent.
#' @keywords internal
#' Locate an event bifurcation bin from direct child-mass normalized entropy.
#'
#' For each cell, computes its relative share across the event's direct children
#' and the normalized Shannon entropy (divided by log2(n_children)). A fork is
#' detected as the first bin where normalized entropy falls below `h_sep` (for
#' 2 consecutive bins), provided a mixed state (entropy >= `h_mixed`) was
#' observed prior to that bin.
#' @keywords internal
.event_child_mass_fork_bin <- function(u, lw, pt, grid, h_mixed = 0.5, h_sep = 0.5) {
  if (length(u$children) < 2L) return(NULL)
  desc <- u$desc
  cm <- vapply(desc, function(d) rowSums(lw[, d, drop = FALSE]), numeric(nrow(lw)))
  keep <- rowSums(cm) > 1e-6
  if (sum(keep) == 0L) return(NULL)

  cell_bin <- findInterval(pt[keep], grid, rightmost.closed = TRUE, all.inside = TRUE)
  share <- cm[keep, , drop = FALSE] / rowSums(cm[keep, , drop = FALSE])
  h_cell <- -rowSums(ifelse(share > 0, share * log2(share), 0)) / log2(length(u$children))
  w <- rowSums(cm[keep, , drop = FALSE])

  h_by_bin <- tapply(h_cell * w, cell_bin, sum) / tapply(w, cell_bin, sum)
  h_vec <- rep(NA_real_, length(grid) - 1L)
  h_vec[as.integer(names(h_by_bin))] <- as.numeric(h_by_bin)

  b_cand <- which(h_vec[-length(h_vec)] < h_sep & h_vec[-1L] < h_sep)
  has_prior_mixed <- vapply(b_cand, function(b) {
    if (b <= 1L) return(FALSE)
    any(h_vec[seq_len(b - 1L)] >= h_mixed, na.rm = TRUE)
  }, logical(1L))
  b_valid <- b_cand[has_prior_mixed]

  if (!length(b_valid)) {
    reason <- if (all(h_vec < h_sep, na.rm = TRUE)) "c4_boundary_onset" else "c3_no_low_entropy"
    return(list(b_star = NA_integer_, reason = reason, h_vec = h_vec))
  }

  b_star <- b_valid[1L]
  list(
    b_star = b_star,
    left_bound = grid[b_star],
    right_bound = grid[b_star + 1L],
    fork_pt = grid[b_star],
    reason = "crossing",
    h_vec = h_vec
  )
}

.fork_event_units <- function(event_scope, boundary, lineage_paths,
                              sibling_map, lineage_clusters, fork_bins,
                              bin_bounds, pt, window_bins, lw,
                              bif_df = NULL, anchor_rule = "min") {
  base <- if (event_scope == "all") {
    .enumerate_fork_events(lineage_paths)
  } else {
    .enumerate_sibling_units(sibling_map)
  }
  if (!length(base)) return(list())

  bb <- NULL
  if (!is.null(bin_bounds)) {
    bb <- if (is.list(bin_bounds) && is.null(dim(bin_bounds))) {
      lapply(bin_bounds, as.numeric)
    } else {
      lapply(seq_len(nrow(bin_bounds)),
             function(i) as.numeric(bin_bounds[i, ]))
    }
    uses_first <- boundary == "exact" || event_scope == "all"
    if (uses_first && length(bb) > 1L &&
        !all(vapply(bb[-1L], function(x) {
          isTRUE(all.equal(x, bb[[1L]], tolerance = 1e-12))
        }, logical(1L)))) {
      stop("`bin_bounds` must be a common grid (identical rows) for this mode",
           call. = FALSE)
    }
  }

  fb <- if (!is.null(fork_bins)) data.table::as.data.table(fork_bins) else NULL
  bst <- if (!is.null(bif_df) &&
             all(c("lineage", "b_star") %in% names(bif_df))) {
    data.table::as.data.table(bif_df)
  } else {
    NULL
  }

  units <- lapply(base, function(u) {
    b0 <- if (!is.null(bb)) {
      if (event_scope == "all") bb[[1L]] else bb[[u$lineage]]
    } else {
      NULL
    }

    anchor_bstar <- NA_real_
    fork_source <- "topology"

    if (anchor_rule == "child_mass_entropy" && event_scope == "all" && !is.null(b0)) {
      det <- .event_child_mass_fork_bin(u, lw, pt, b0)
      if (!is.null(det) && is.finite(det$b_star)) {
        anchor_bstar <- as.numeric(det$b_star)
        fork_source <- "child_mass_entropy"
      } else if (!is.null(det) && det$reason == "c4_boundary_onset") {
        anchor_bstar <- 1
        fork_source <- "boundary_onset"
      }
    } else if (!is.null(bst)) {
      child_bstars <- vapply(seq_along(u$children), function(k) {
        lins <- if (event_scope == "all") u$desc[[k]] else u$children[k]
        sel <- bst$lineage %in% lins
        if ("reason" %in% names(bst)) {
          sel <- sel & bst$reason %in% "crossing"
        }
        v <- bst$b_star[sel]
        if (length(v)) mean(v) else NA_real_
      }, numeric(1L))
      if (any(is.finite(child_bstars))) {
        cb <- child_bstars[is.finite(child_bstars)]
        anchor_bstar <- if (anchor_rule == "mean") mean(cb) else min(cb)
        fork_source <- "bifurcation_crossing"
      }
    }

    fork_pt <- NA_real_
    if (is.finite(anchor_bstar) && !is.null(b0)) {
      fork_pt <- if (boundary == "exact") {
        b0[max(1L, min(length(b0), round(anchor_bstar) + 1L))]
      } else {
        b0[max(1L, min(length(b0) - 1L, round(anchor_bstar)))]
      }
    }
    if (!is.finite(fork_pt) && event_scope == "all" &&
        !is.null(lineage_clusters)) {
      node_cells <- lineage_clusters == u$event_node
      if (sum(node_cells) > 0) fork_pt <- stats::median(pt[node_cells])
    }
    if (!is.finite(fork_pt) && !is.null(fb) && !is.null(u$lineage) &&
        "fork_pseudotime" %in% names(fb)) {
      row <- fb[fb$lineage == u$lineage]
      if (nrow(row) > 0 && is.finite(row$fork_pseudotime[1])) {
        fork_pt <- row$fork_pseudotime[1]
      }
    }
    if (!is.finite(fork_pt)) return(NULL)

    b_star <- NA_integer_
    if (boundary == "exact") {
      if (is.null(bb)) {
        stop("`boundary='exact'` requires `bin_bounds`", call. = FALSE)
      }
      b <- bb[[1L]]
      n_b <- length(b) - 1L
      fb_idx <- max(1L, min(n_b, sum(b < fork_pt)))
      on_edge <- isTRUE(all.equal(b[fb_idx + 1L], fork_pt))
      left_bound <- b[fb_idx]
      right_bound <- b[min(n_b + 1L, fb_idx + 1L)]
      pre_lo <- b[max(1L, fb_idx - window_bins + 1L)]
      pre_hi <- fork_pt
      post_lo <- fork_pt
      post_hi <- b[min(n_b + 1L, fb_idx + window_bins + as.integer(on_edge))]
    } else {
      if (is.null(bb)) {
        stop("`boundary='fixed_bin'` requires `bin_bounds`", call. = FALSE)
      }
      b <- if (event_scope == "all") bb[[1L]] else bb[[u$lineage]]
      fb_bin <- if (is.finite(anchor_bstar)) {
        round(anchor_bstar)
      } else {
        sum(b <= fork_pt)
      }
      fb_bin <- max(1L, min(length(b) - 1L, fb_bin))
      left_bound <- b[fb_bin]
      right_bound <- b[min(length(b), fb_bin + 1L)]
      pre_lo <- b[max(1L, fb_bin - window_bins)]
      pre_hi <- left_bound
      post_lo <- right_bound
      post_hi <- b[min(length(b), fb_bin + window_bins + 1L)]
      fork_pt <- left_bound
      b_star <- fb_bin
    }
    u$fork_pt <- fork_pt
    u$left_bound <- left_bound
    u$right_bound <- right_bound
    u$pre_lo <- pre_lo
    u$pre_hi <- pre_hi
    u$post_lo <- post_lo
    u$post_hi <- post_hi
    u$grid <- b
    u$b_star <- b_star
    u$fork_time_source <- fork_source
    u
  })
  Filter(Negate(is.null), units)
}

#' Enumerate all fork events from lineage paths. Each child carries its
#' descendant terminal lineage indices (`desc`).
#' @keywords internal
.enumerate_fork_events <- function(lineage_paths) {
  paths <- lineage_paths
  n_paths <- length(paths)

  edge_set <- unique(data.table::rbindlist(lapply(seq_len(n_paths), function(i) {
    p <- paths[[i]]
    if (length(p) < 2L) return(NULL)
    data.table::data.table(parent = p[-length(p)], child = p[-1L])
  })))
  children_of <- tapply(edge_set$child, edge_set$parent,
                        function(z) sort(unique(z)))
  fork_nodes <- names(children_of)[
    vapply(children_of, length, integer(1L)) >= 2L
  ]

  lapply(fork_nodes, function(fn) {
    kids <- children_of[[fn]]
    term_map <- data.table::rbindlist(lapply(seq_len(n_paths), function(pi) {
      p <- paths[[pi]]
      pos <- match(as.integer(fn), p)
      if (is.na(pos) || pos == length(p)) return(NULL)
      data.table::data.table(child = p[pos + 1L], terminal = pi)
    }))
    desc <- lapply(kids, function(k) sort(term_map[child == k, terminal]))
    list(event_node = as.integer(fn), children = kids,
         desc = desc, lineage = NULL)
  })
}

#' Enumerate nearest-sibling comparison units: one unit per focal lineage
#' whose children are the focal plus all its nearest siblings. Only the focal
#' child is emitted; the siblings form the contrast set.
#' @keywords internal
.enumerate_sibling_units <- function(sibling_map) {
  smap <- data.table::as.data.table(sibling_map)
  foc_groups <- split(seq_len(nrow(smap)), smap$focal_lineage)
  lapply(names(foc_groups), function(foc_name) {
    rows <- smap[foc_groups[[foc_name]]]
    foc <- as.integer(foc_name)
    sibs <- sort(unique(rows$sibling_lineage))
    list(event_node = rows$fork_cluster[1L],
         children = c(foc, sibs), desc = NULL, lineage = foc)
  })
}

#' Attach signed late support within the fixed local sibling panel.
#'
#' A short post-fork increase is capped by the focal-minus-other contrast in
#' the later half of their shared, supported trajectory. The original
#' `sibling_contrast` remains the directional effect estimate.
#' Phase summaries additionally compare every evaluable event with its exact
#' native sister panel. They use equally weighted sisters, retain observed
#' zeros, and subtract the corresponding pre-event activities. Their signs
#' describe early and late branch association, not causal initiation or feedback.
#'
#' @param fork_contrasts Event-based fork contrast table.
#' @param temporal_lr Output of `ComputeTemporalLRWithIntermediate()` on the
#'   same lineage weights and pseudotime as `fork_contrasts`.
#' @param lineage_paths Root-to-terminal cluster paths for that run.
#' @param sibling_scope Pool observed cells in the fixed local sibling panel
#'   (default), require a strict common horizon, use historical unrestricted
#'   descendants, or compare pairwise observed horizons as a mechanism control.
#' @param fallback_temporal_lr Optional adjacent-bin grid on the same cells,
#'   scale and support threshold, used only when fine late support is absent.
#' @return `fork_contrasts` with `late_contrast`, `n_late_bins`, and
#'   `persistent_positive_score` columns, together with native-panel late
#'   activities, baseline-relative changes, persistence_support and phase_pattern.
#' @export
AddLateForkPersistence <- function(fork_contrasts, temporal_lr,
                                   lineage_paths,
                                   expression_mean_mode = c("all_cells", "positive_cells"),
                                   sibling_scope = c("native", "pooled_native", "all_descendants", "pairwise"),
                                   persistence_mode = c("signed", "positive"),
                                   fallback_temporal_lr = temporal_lr$late_fallback) {
  expression_mean_mode <- match.arg(expression_mean_mode)
  sibling_scope <- match.arg(sibling_scope)
  persistence_mode <- match.arg(persistence_mode)
  score_field <- if (expression_mean_mode == "positive_cells") {
    "score_positive"
  } else {
    "score"
  }
  out <- fork_contrasts
  out$late_contrast <- rep(NA_real_, nrow(out))
  out$n_late_bins <- rep(NA_integer_, nrow(out))
  out$persistence_phase_start <- out$persistence_phase_end <- rep(NA_real_, nrow(out))
  out$persistence_focal_activity <- out$persistence_other_activity <- rep(NA_real_,nrow(out))
  # A shared parent can make a local contrast observable before either
  # branch has an observable late horizon. That absence is not persistence:
  # retain only local change supported in the focal direction.
  unobserved_score <- switch(persistence_mode,
    signed = sign(out$sibling_contrast) * pmin(abs(out$sibling_contrast),
      pmax(0, sign(out$sibling_contrast) * out$post_fork_change)),
    positive = out$sibling_contrast)
  out$persistent_positive_score <- unobserved_score
  out$persistence_support <- ifelse(out$evaluable & is.finite(out$sibling_contrast),
    ifelse(persistence_mode == "signed", "local_focal_only", "local_contrast_only"),
    "not_evaluable")
  out$persistence_change_contrast <- rep(NA_real_, nrow(out))
  phase_fields <- c("late_focal_activity_native", "late_other_activity_native",
    "late_change_focal", "late_change_other", "late_change_contrast",
    "late_phase_start", "late_phase_end", "n_late_phase_bins")
  out[, phase_fields] <- rep(list(rep(NA_real_, nrow(out))), length(phase_fields))
  out$phase_pattern <- ifelse(out$evaluable & is.finite(out$sibling_contrast),
    "late_unobserved", "not_evaluable")
  events <- data.table::as.data.table(fork_contrasts)
  events[, event_row := .I]
  events <- events[
    evaluable %in% TRUE & is.finite(sibling_contrast)
  ]
  if (!nrow(events)) return(out)
  # The optional adjacent-bin grid is scored from the same original cells,
  # at the same scale and effective-count threshold. Fine observations win.
  fallback <- switch(as.character(is.null(fallback_temporal_lr)),
    "TRUE" = NULL, "FALSE" = {
      coarse <- AddLateForkPersistence(fork_contrasts, fallback_temporal_lr,
        lineage_paths, expression_mean_mode, sibling_scope, persistence_mode,
        fallback_temporal_lr = NULL)
      observed <- is.finite(coarse$persistence_change_contrast) & coarse$n_late_bins >= 1L
      coarse$persistence_support[observed] <- "late_observed_adjacent_bins"
      coarse
    })

  paths <- data.table::as.data.table(stack(setNames(
    lineage_paths, seq_along(lineage_paths)
  )))
  data.table::setnames(paths, c("values", "ind"), c("node", "lineage"))
  paths[, lineage := as.integer(as.character(lineage))]
  paths[, next_node := data.table::shift(node, type = "lead"), by = lineage]
  members <- merge(
    events[, .(event_row, pair_idx, source, target, event_node,
               child, fork_time, sibling_children)],
    paths, by.x = "event_node", by.y = "node", allow.cartesian = TRUE
  )
  if (sibling_scope %in% c("native", "pooled_native", "pairwise")) {
    selected <- events[, .(next_node = as.integer(unlist(strsplit(
      sibling_children, ",", fixed = TRUE)))), by = event_row]
    selected <- unique(data.table::rbindlist(list(
      selected, events[, .(event_row, next_node = child)]
    )))
    members <- merge(members, selected, by = c("event_row", "next_node"))
    members[, branch_group := data.table::fifelse(
      next_node == child, "focal", "other")]
  } else if (attr(fork_contrasts, "event_units")$sibling_selection %in%
      c("time_matched", "program_matched")) {
    members[, selected_child := as.integer(sibling_children)]
    members <- members[next_node == child | next_node == selected_child]
    members[, branch_group := data.table::fifelse(
      next_node == child, "focal", "other")]
  } else {
    focal_paths <- unique(paths[, .(child = node, lineage)])
    focal_paths[, focal := TRUE]
    members <- merge(members, focal_paths,
                     by = c("child", "lineage"), all.x = TRUE)
    members[, branch_group := data.table::fifelse(
      focal %in% TRUE, "focal", "other"
    )]
  }

  bounds <- data.table::as.data.table(stack(setNames(
    temporal_lr$intermediate$bin_bounds,
    seq_along(temporal_lr$intermediate$bin_bounds)
  )))
  data.table::setnames(bounds, c("values", "ind"),
                       c("edge", "lineage"))
  bounds[, lineage := as.integer(as.character(lineage))]
  bounds[, bin := seq_len(.N) - 1L, by = lineage]
  bounds[, mid := (edge + data.table::shift(edge)) / 2, by = lineage]
  bounds <- bounds[bin > 0L, .(lineage, bin, mid)]

  temporal <- data.table::as.data.table(temporal_lr$scores)[
    pair_idx %in% events$pair_idx,
    c("pair_idx", "source", "target", "lineage", "bin", score_field,
      "weight_mass_cells", "structural_valid", "expression_supported",
      intersect("sample_component", names(temporal_lr$scores))), with = FALSE
  ]
  data.table::setnames(temporal, score_field, "score")
  if (identical(attr(fork_contrasts, "event_units")$activity_transform,
                "log_exposure")) {
    temporal[, score := -log1p(-score)]
  }
  joined <- merge(members, temporal,
                  by = c("pair_idx", "source", "target", "lineage"),
                  allow.cartesian = TRUE)
  joined <- merge(joined, bounds, by = c("lineage", "bin"))
  # Expression-gated zeros are observed activity. Only missing structural
  # support removes a bin from the shared late window.
  joined <- joined[
    structural_valid %in% TRUE &
      is.finite(score) & is.finite(weight_mass_cells) &
      weight_mass_cells > 0 & is.finite(fork_time) & mid > fork_time
  ]
  # A phase contrast needs the same sister identities as the local contrast.
  # Missing late support for one of those sisters is not a zero measurement.
  native_panel <- events[, .(next_node = as.integer(unlist(strsplit(
    sibling_children, ",", fixed = TRUE)))), by = event_row]
  native_panel <- unique(data.table::rbindlist(list(native_panel,
    events[, .(event_row, next_node = child)])))
  native_panel[, n_panel := .N, by = event_row]
  phase <- merge(joined, native_panel, by = c("event_row", "next_node"))
  phase_ends <- phase[, .(group_end = max(mid), n_panel = n_panel[1L]),
    by = .(event_row, next_node)][,
      .(common_end = min(group_end), complete = .N == n_panel[1L]),
      by = event_row][complete %in% TRUE]
  phase <- merge(phase, phase_ends[, .(event_row, common_end)], by = "event_row")
  phase <- phase[mid <= common_end & mid >= fork_time + (common_end - fork_time) / 2]
  phase_cells <- phase[, .(activity = sum(score * weight_mass_cells) / sum(weight_mass_cells),
    focal = next_node[1L] == child[1L], n_panel = n_panel[1L]),
    by = .(event_row, mid, next_node)]
  phase_bins <- phase_cells[, .(focal = sum(activity[focal]),
    other = mean(activity[!focal]), complete = .N == n_panel[1L] & sum(focal) == 1L),
    by = .(event_row, mid)][complete %in% TRUE]
  switch(as.character("sample_component" %in% names(phase)), "TRUE" = {
    phase_components <- phase[, .(mass = sum(weight_mass_cells), child = child[1L],
      n_panel = n_panel[1L]), by = .(event_row, mid, next_node, sample_component)]
    phase_components[, coefficient := data.table::fifelse(next_node == child,
      1, -1 / (n_panel - 1)) * mass / sum(mass), by = .(event_row, mid, next_node)]
    phase_identifiable <- phase_components[, .(balance = sum(coefficient)),
      by = .(event_row, mid, sample_component)][,
        .(identifiable = all(is.finite(sample_component)) &
          max(abs(balance)) < sqrt(.Machine$double.eps)), by = .(event_row, mid)]
    phase_bins <- merge(phase_bins,
      phase_identifiable[identifiable %in% TRUE, .(event_row, mid)], by = c("event_row", "mid"))
  })
  phase_summary <- phase_bins[, .(late_focal_activity_native = mean(focal),
    late_other_activity_native = mean(other), late_phase_start = sort(mid)[1L],
    late_phase_end = sort(mid, decreasing = TRUE)[1L], n_late_phase_bins = .N), by = event_row]
  # Each fixed sibling is compared on its own horizon shared with the focal
  # child. A short sibling must not censor late observations for every pair.
  switch(sibling_scope, pairwise = {
    pair_panel <- events[, .(sibling_child = as.integer(unlist(strsplit(
      sibling_children, ",", fixed = TRUE)))), by = event_row]
    pair_rows <- merge(joined[, c("event_row", "next_node", "child", "fork_time",
      "mid", "score", "weight_mass_cells", intersect("sample_component", names(joined))),
      with = FALSE], pair_panel, by = "event_row", allow.cartesian = TRUE)[
        next_node == child | next_node == sibling_child]
    pair_ends <- pair_rows[, .(group_end = max(mid)),
      by = .(event_row, sibling_child, next_node)][,
        .(pair_end = min(group_end), n_children = .N),
        by = .(event_row, sibling_child)][n_children == 2L]
    pair_rows <- merge(pair_rows, pair_ends[, .(event_row, sibling_child, pair_end)],
      by = c("event_row", "sibling_child"))
    pair_rows <- pair_rows[mid <= pair_end &
      mid >= fork_time + (pair_end - fork_time) / 2]
    pair_cells <- pair_rows[, .(activity = sum(score * weight_mass_cells) / sum(weight_mass_cells),
      focal = next_node[1L] == child[1L]), by = .(event_row, sibling_child, mid, next_node)]
    pair_bins <- pair_cells[, .(focal = sum(activity[focal]),
      other = sum(activity[!focal]), complete = .N == 2L & sum(focal) == 1L),
      by = .(event_row, sibling_child, mid)][complete %in% TRUE]
    switch(as.character("sample_component" %in% names(pair_rows)), "TRUE" = {
      pair_components <- pair_rows[, .(mass = sum(weight_mass_cells), child = child[1L]),
        by = .(event_row, sibling_child, mid, next_node, sample_component)]
      pair_components[, coefficient := data.table::fifelse(next_node == child, 1, -1) *
        mass / sum(mass), by = .(event_row, sibling_child, mid, next_node)]
      pair_identifiable <- pair_components[, .(balance = sum(coefficient)),
        by = .(event_row, sibling_child, mid, sample_component)][,
          .(identifiable = all(is.finite(sample_component)) &
            max(abs(balance)) < sqrt(.Machine$double.eps)),
          by = .(event_row, sibling_child, mid)]
      pair_bins <- merge(pair_bins, pair_identifiable[identifiable %in% TRUE,
        .(event_row, sibling_child, mid)], by = c("event_row", "sibling_child", "mid"))
    })
    pair_windows <- pair_bins[, .(focal = mean(focal), other = mean(other),
      phase_start = min(mid), phase_end = max(mid), n_bins = .N),
      by = .(event_row, sibling_child)]
    pair_required <- pair_panel[, .(n_required = .N), by = event_row]
    pair_summary <- pair_windows[, .(late_focal_activity_native = mean(focal),
      late_other_activity_native = mean(other), late_phase_start = min(phase_start),
      late_phase_end = max(phase_end), n_late_phase_bins = min(n_bins), n_available = .N),
      by = event_row]
    phase_summary <- merge(pair_summary, pair_required, by = "event_row")[
      n_available == n_required, .(event_row, late_focal_activity_native,
        late_other_activity_native, late_phase_start, late_phase_end, n_late_phase_bins)]
    attr(out, "late_pair_windows") <- as.data.frame(merge(pair_windows,
      events[, .(event_row, unit_id, pair_idx, source, target, child)], by = "event_row"))
  })
  phase_summary <- merge(phase_summary,
    events[, .(event_row, focal_pre, other_pre, sibling_contrast)], by = "event_row")
  phase_summary[, `:=`(late_change_focal = late_focal_activity_native - focal_pre,
    late_change_other = late_other_activity_native - other_pre)]
  phase_summary[, late_change_contrast := late_change_focal - late_change_other]
  out[phase_summary$event_row, phase_fields] <- as.data.frame(phase_summary[, phase_fields, with = FALSE])
  # Native persistence uses the complete fixed-child phase summary below.
  # The pooled-descendant calculation is retained for the explicit ablation.
  if (sibling_scope %in% c("native","pairwise")) {
    # The persistence reference is the same fixed sibling panel and equal
    # child weights as the local contrast, on bins observed for every child.
    late <- phase_summary[, .(event_row,
      late_contrast = late_focal_activity_native - late_other_activity_native,
      persistence_focal_activity = late_focal_activity_native,
      persistence_other_activity = late_other_activity_native,
      n_late_bins = n_late_phase_bins,
      persistence_phase_start = late_phase_start, persistence_phase_end = late_phase_end)]
  } else {
    if (!nrow(joined)) return(switch(as.character(is.null(fallback)), "TRUE" = out, "FALSE" = fallback))
    ends <- joined[, .(group_end = max(mid)),
                   by = .(event_row, branch_group)]
    ends <- ends[, .(common_end = min(group_end), n_groups = .N),
                 by = event_row][n_groups == 2L,
                                 .(event_row, common_end)]
    joined <- merge(joined, ends, by = "event_row")
    joined <- joined[
      mid <= common_end & mid >=
        fork_time + (common_end - fork_time) / 2
    ]
    bin_score <- joined[, .(
      score = sum(score * weight_mass_cells) / sum(weight_mass_cells)
    ), by = .(event_row, bin, branch_group)]
    if (!all(c("focal", "other") %in% bin_score$branch_group)) return(switch(as.character(is.null(fallback)), "TRUE" = out, "FALSE" = fallback))
    matched <- data.table::dcast(bin_score, event_row + bin ~ branch_group,
                                 value.var = "score")
    matched <- matched[is.finite(focal) & is.finite(other)]
    switch(as.character("sample_component" %in% names(joined)), "TRUE" = {
      component_weights <- joined[, .(mass = sum(weight_mass_cells)),
        by = .(event_row, bin, branch_group, sample_component)]
      component_weights[, coefficient := data.table::fifelse(branch_group == "focal", 1, -1) *
        mass / sum(mass), by = .(event_row, bin, branch_group)]
      identifiable <- component_weights[, .(balance = sum(coefficient)),
        by = .(event_row, bin, sample_component)][,
          .(identifiable = all(is.finite(sample_component)) &
            max(abs(balance)) < sqrt(.Machine$double.eps)), by = .(event_row, bin)]
      matched <- merge(matched, identifiable[identifiable %in% TRUE, .(event_row, bin)],
        by = c("event_row", "bin"))
    })
    late <- matched[, .(late_contrast = mean(focal - other),
      persistence_focal_activity=mean(focal),persistence_other_activity=mean(other),
      n_late_bins = .N), by = event_row]
    used_phase <- merge(unique(joined[, .(event_row, bin, mid)]),
      matched[, .(event_row, bin)], by = c("event_row", "bin"))[,
        .(persistence_phase_start = min(mid), persistence_phase_end = max(mid)), by = event_row]
    late <- merge(late, used_phase, by = "event_row")
  }
  out$late_contrast[late$event_row] <- late$late_contrast
  out$n_late_bins[late$event_row] <- late$n_late_bins
  out$persistence_phase_start[late$event_row] <- late$persistence_phase_start
  out$persistence_phase_end[late$event_row] <- late$persistence_phase_end
  out$persistence_focal_activity[late$event_row] <- late$persistence_focal_activity
  out$persistence_other_activity[late$event_row] <- late$persistence_other_activity
  out$persistence_change_contrast <-
    (out$persistence_focal_activity - out$focal_pre) -
    (out$persistence_other_activity - out$other_pre)
  known_late <- is.finite(out$persistence_change_contrast) & out$n_late_bins >= 1L
  out$persistence_support[known_late] <- "late_observed"
  out$phase_pattern[known_late] <- data.table::fcase(
    out$sibling_contrast[known_late] > 0 & out$persistence_change_contrast[known_late] > 0,
      "early_and_late_positive",
    out$sibling_contrast[known_late] > 0 & out$persistence_change_contrast[known_late] <= 0,
      "early_positive_only",
    out$sibling_contrast[known_late] <= 0 & out$persistence_change_contrast[known_late] > 0,
      "late_positive_only", default = "no_positive_phase")
  out$persistent_positive_score <- data.table::fifelse(
    out$sibling_contrast > 0 & is.finite(out$late_contrast) &
      out$n_late_bins >= 1L,
    pmin(out$sibling_contrast, out$late_contrast),
    out$sibling_contrast
  )
  switch(persistence_mode, signed = {
    # Failure to persist removes evidence in that direction; it does not
    # reverse the local effect or exempt decreases from the same operator.
    out$persistent_positive_score <- data.table::fifelse(
      is.finite(out$persistence_change_contrast) & out$n_late_bins >= 1L,
      sign(out$sibling_contrast) * pmin(abs(out$sibling_contrast),
        pmax(0, sign(out$sibling_contrast) * out$persistence_change_contrast)),
      unobserved_score
    )
  })
  switch(as.character(!is.null(fallback)), "TRUE" = {
    take <- !known_late & is.finite(fallback$persistence_change_contrast) &
      fallback$n_late_bins >= 1L
    fields <- c("late_contrast", "n_late_bins", "persistence_phase_start",
      "persistence_phase_end", "persistence_focal_activity",
      "persistence_other_activity", "persistent_positive_score",
      "persistence_change_contrast", phase_fields, "phase_pattern", "persistence_support")
    out[take, fields] <- fallback[take, fields]
  })
  out
}

#' Describe branch-associated temporal shapes on the observed native panel
#'
#' Child intercepts and shared time effects separate branch-associated changes
#' from common temporal trends. Each sister contributes for its observed horizon;
#' bins require the focal child and at least one observed native sister. Constant,
#' linear, plateau, delayed and two-slope models compete by weighted BIC. These
#' are descriptive RNA activity shapes, not calibrated tests of communication
#' causality. Existing event and persistence rankings remain unchanged.
#' Focal activity timing is estimated separately because a branch contrast can
#' peak after focal activity when sister activity declines. Focal timing alone
#' does not establish branch specificity; its observed horizon is reported.
#' @param fork_contrasts Native event contrast table.
#' @param temporal_lr Temporal LR scores and their original bin bounds.
#' @param lineage_paths Root-to-terminal paths for those scores.
#' @param expression_mean_mode Activity field used by the native event analysis.
#' @param background Use child and shared-time effects, or a focal-only ablation.
#' @param sample_labels Biological sample for each row of temporal_lr$scores.
#' @param condition_on_sample Retain within-sample backgrounds in state mode.
#' @param time_coordinate State pseudotime (default), or observed calendar time.
#'   Calendar scores require an observed_time column and one time per biological
#'   sample. They compare independent sample observations with sample plus child
#'   intercepts; state scores retain the original within-sample conditioning.
#'   Calendar profile support is assessed independently of state-coordinate
#'   pre/post event eligibility; neither mode changes native event rankings.
#' @param condition_labels Optional observed-condition labels, one per temporal
#'   score row. Conditions are fitted separately, retaining event identities.
#' @param reference_conditions Shared reference labels included in each condition
#'   fit, without treating reuse of a reference as additional replication.
#'   Calendar fits also compare one-interval steps and finite-duration pulses.
#'   Their transition brackets use adjacent sampled times and are distinct from
#'   slope-change and peak estimates; they do not change pulse ranking evidence.
#' @return The input table with relative dynamic shape, support, slopes and BIC
#' gains. The branch-conditioned mode also returns focal_dynamic_* summaries,
#' focal_pulse_bic_gain and focal_peak_time, defined only for a focal pulse.
#' @export
SummarizeForkDynamics <- function(fork_contrasts, temporal_lr, lineage_paths,
                                  expression_mean_mode = c("all_cells", "positive_cells"),
                                  background = c("branch_conditioned", "focal_only"),
                                  sample_labels = rep.int("pooled", nrow(temporal_lr$scores)),
                                  condition_on_sample = TRUE,
                                  time_coordinate = c("state_pseudotime", "observed_calendar"),
                                  condition_labels = NULL,
                                  reference_conditions = character(),
                                  sample_contrasts = NULL,transfer_pulse=TRUE,include_level=TRUE) {
  expression_mean_mode <- match.arg(expression_mean_mode)
  background <- match.arg(background)
  time_coordinate <- match.arg(time_coordinate)
  switch(time_coordinate,observed_calendar={
    switch(as.character(!is.null(temporal_lr$lineage_paths)),"TRUE"={
      lineage_paths <- temporal_lr$lineage_paths
    })
  })
  switch(as.character(!is.null(condition_labels)), "TRUE" = {
    conditions <- as.character(condition_labels)
    stopifnot(length(conditions) == nrow(temporal_lr$scores), !anyNA(conditions))
    strata <- setdiff(unique(conditions), reference_conditions)
    conditioned <- lapply(strata, function(stratum) {
      keep <- conditions %in% c(stratum, reference_conditions)
      observed <- temporal_lr
      observed$scores <- data.table::as.data.table(temporal_lr$scores)[keep]
      value <- SummarizeForkDynamics(fork_contrasts, observed, lineage_paths,
        expression_mean_mode = expression_mean_mode, background = background,
        sample_labels = sample_labels[keep], condition_on_sample = condition_on_sample,
        time_coordinate = time_coordinate,sample_contrasts=sample_contrasts,transfer_pulse=transfer_pulse,include_level=include_level)
      value$dynamic_condition <- stratum
      value
    })
    result <- as.data.frame(data.table::rbindlist(conditioned, use.names = TRUE, fill = TRUE))
    attr(result, "event_units") <- attr(fork_contrasts, "event_units")
    return(result)
  })
  out <- fork_contrasts
  out$cohort_effect <- rep(NA_real_,nrow(out))
  out$cohort_pulse_weight <- rep(0,nrow(out))
  if (!is.null(sample_contrasts)) {
    local_field <- switch(background,branch_conditioned="sibling_contrast",focal_only="post_fork_change")
    ids <- intersect(c("unit_id","pair_idx","source","target","event_node","child"),names(out))
    native_start <- data.table::as.data.table(sample_contrasts)[evaluable%in%TRUE&is.finite(get(local_field)),
      .(effect=stats::median(get(local_field))),by=ids]
    locations <- merge(data.table::as.data.table(out)[,c(ids),with=FALSE][,event_row:=.I],
      native_start,by=ids)
    out$cohort_effect <- out[[local_field]]
    out$cohort_effect[locations$event_row] <- locations$effect
  }
  out$dynamic_time_coordinate <- time_coordinate
  switch(time_coordinate, observed_calendar = {
    out$dynamic_observation_scope <- temporal_lr$observation_scope
  })
  fields <- c("dynamic_time_start", "dynamic_time_end", "dynamic_change_time",
    "dynamic_slope_left", "dynamic_slope_right", "dynamic_bic_gain",
    "pulse_bic_gain", "dynamic_n_bins", "dynamic_n_observations", "dynamic_n_samples",
    "dynamic_peak_lower", "dynamic_peak_upper", "dynamic_effect", "dynamic_effect_time", "dynamic_level_effect",
    "dynamic_observed_peak_time","dynamic_observed_peak_lower","dynamic_observed_peak_upper")
  out[, fields] <- rep(list(rep(NA_real_, nrow(out))), length(fields))
  out$dynamic_sample_ids <- rep(NA_character_,nrow(out))
  transition_fields <- c("dynamic_transition_lower", "dynamic_transition_upper",
    "dynamic_transition_effect", "dynamic_transition_bic_gain",
    "dynamic_offset_lower", "dynamic_offset_upper")
  out[, transition_fields] <- rep(list(rep(NA_real_, nrow(out))), length(transition_fields))
  out$dynamic_transition_model <- rep(NA_character_, nrow(out))
  out$dynamic_transition_status <- rep(switch(time_coordinate,
    observed_calendar = "insufficient_observed_overlap",
    state_pseudotime = "not_estimated_on_state_coordinate"), nrow(out))
  out$dynamic_model <- out$dynamic_pattern <- rep(NA_character_, nrow(out))
  out$dynamic_status <- ifelse(time_coordinate == "observed_calendar" | (out$evaluable & is.finite(out$sibling_contrast)),
    "insufficient_observed_overlap", "not_evaluable")
  switch(background, branch_conditioned = {
    focal <- SummarizeForkDynamics(fork_contrasts, temporal_lr, lineage_paths,
      expression_mean_mode = expression_mean_mode, background = "focal_only",
      sample_labels = sample_labels, condition_on_sample = condition_on_sample,
      time_coordinate = time_coordinate,sample_contrasts=sample_contrasts,transfer_pulse=transfer_pulse,include_level=include_level)
    focal_fields <- c("dynamic_bic_gain", "dynamic_time_start", "dynamic_time_end", "dynamic_change_time",
      "dynamic_n_bins", "dynamic_model", "dynamic_pattern", "dynamic_status", "pulse_bic_gain",
      "dynamic_n_samples", "dynamic_peak_lower", "dynamic_peak_upper", "dynamic_effect", "dynamic_effect_time", "dynamic_level_effect",
      "dynamic_observed_peak_time","dynamic_observed_peak_lower","dynamic_observed_peak_upper",
      transition_fields, "dynamic_transition_model",
      "dynamic_transition_status","cohort_effect","cohort_pulse_weight")
    out[, paste0("focal_", focal_fields)] <- as.data.frame(focal)[, focal_fields, drop = FALSE]
    out$focal_peak_time <- ifelse(focal$dynamic_pattern %in% "pulse",
      (focal$dynamic_peak_lower + focal$dynamic_peak_upper) / 2, NA_real_)
    out$event_peak_time <- focal$dynamic_observed_peak_time
    out$event_peak_lower <- focal$dynamic_observed_peak_lower
    out$event_peak_upper <- focal$dynamic_observed_peak_upper
  })
  events <- data.table::copy(data.table::as.data.table(fork_contrasts))[, event_row := .I]
  events <- events[time_coordinate == "observed_calendar" | (evaluable %in% TRUE & is.finite(sibling_contrast))]
  switch(as.character(nrow(events) == 0L), "TRUE" = return(out))
  if (time_coordinate == "state_pseudotime" &&
      identical(temporal_lr$observation_scope, "native_event_segment_populations")) {
    panel <- events[, .(next_node = as.integer(unlist(strsplit(sibling_children, ",", fixed = TRUE)))),
      by = event_row]
    panel <- unique(data.table::rbindlist(list(panel, events[, .(event_row, next_node = child)])))
    temporal <- data.table::copy(data.table::as.data.table(temporal_lr$scores))
    temporal[, biological_sample := as.character(sample_labels)]
    score_field <- switch(expression_mean_mode, all_cells = "score", positive_cells = "score_positive")
    temporal[, activity := get(score_field)]
    if (identical(attr(fork_contrasts, "event_units")$activity_transform, "log_exposure")) {
      temporal[, activity := -log1p(-activity)]
    }
    joined <- merge(events[, .(event_row, event_node, child, pair_idx, source, target)],
      temporal[, .(event_node, pair_idx, source, target, next_node, bin, mid,
        activity, structural_valid, weight_mass_cells, biological_sample)],
      by = c("event_node", "pair_idx", "source", "target"), allow.cartesian = TRUE)
    joined <- merge(joined, panel, by = c("event_row", "next_node"))
    if (!is.null(temporal_lr$paired_events)) {
      joined <- merge(joined, temporal_lr$paired_events,
        by = c("pair_idx", "source", "target", "event_node", "child", "biological_sample", "next_node"))
    }
  } else {
  paths <- data.table::data.table(lineage = rep(seq_along(lineage_paths), lengths(lineage_paths)),
    node = unlist(lineage_paths))
  paths[, next_node := data.table::shift(node, type = "lead"), by = lineage]
  members <- merge(events[, .(event_row, event_node, child, pair_idx, source, target)],
    paths, by.x = "event_node", by.y = "node", allow.cartesian = TRUE)
  switch(as.character(time_coordinate=="observed_calendar" &&
    identical(temporal_lr$observation_scope,"observed_native_branch_subtrees")),"TRUE"={
    terminal_nodes <- vapply(lineage_paths, tail, numeric(1L), 1L)
    members <- members[terminal_nodes[lineage] == next_node]
  })
  panel <- events[, .(next_node = as.integer(unlist(strsplit(sibling_children, ",", fixed = TRUE)))),
    by = event_row]
  panel <- unique(data.table::rbindlist(list(panel, events[, .(event_row, next_node = child)])))
  members <- merge(members, panel, by = c("event_row", "next_node"))
  bounds <- data.table::data.table(lineage = rep(seq_along(temporal_lr$intermediate$bin_bounds),
    lengths(temporal_lr$intermediate$bin_bounds)), edge = unlist(temporal_lr$intermediate$bin_bounds))
  bounds[, bin := seq_len(.N) - 1L, by = lineage]
  bounds[, mid := (edge + data.table::shift(edge)) / 2, by = lineage]
  score_field <- switch(expression_mean_mode, all_cells = "score", positive_cells = "score_positive")
  temporal <- data.table::copy(data.table::as.data.table(temporal_lr$scores))
  temporal[, biological_sample := as.character(sample_labels)]
  temporal[, activity := get(score_field)]
  temporal[, observed_calendar_time := switch(time_coordinate,
    state_pseudotime = NA_real_, observed_calendar = observed_time)]
  switch(time_coordinate, observed_calendar = {
    stopifnot(all(is.finite(temporal$observed_calendar_time)),
      all(temporal[, .(n_times = data.table::uniqueN(observed_calendar_time)),
        by = biological_sample]$n_times == 1L))
  })
  switch(as.character(identical(attr(fork_contrasts, "event_units")$activity_transform,
    "log_exposure")), "TRUE" = temporal[, activity := -log1p(-activity)])
  temporal <- temporal[structural_valid %in% TRUE & is.finite(activity) &
    is.finite(weight_mass_cells) & weight_mass_cells > 0]
  joined <- merge(members, temporal[, .(pair_idx, source, target, lineage, bin,
    activity, structural_valid, weight_mass_cells, biological_sample, observed_calendar_time)],
    by = c("pair_idx", "source", "target", "lineage"), allow.cartesian = TRUE)
  if (time_coordinate == "state_pseudotime" && !is.null(temporal_lr$paired_events)) {
    joined <- merge(joined,temporal_lr$paired_events,
      by=c("pair_idx","source","target","event_node","child","biological_sample","next_node"))
  }
  joined <- merge(joined, bounds[bin > 0L, .(lineage, bin, mid)], by = c("lineage", "bin"))
  switch(time_coordinate, observed_calendar = joined[, mid := observed_calendar_time])
  }
  curve <- joined[structural_valid %in% TRUE & is.finite(activity) &
    is.finite(weight_mass_cells) & weight_mass_cells > 0,
    .(activity = sum(activity * weight_mass_cells) / sum(weight_mass_cells),
      mass = sum(weight_mass_cells), focal = next_node[1L] == child[1L]),
    by = .(event_row, biological_sample, mid, next_node)]
  shared_parent <- !identical(attr(fork_contrasts,"event_units")$prefork_weighting,"descendant")
  anchors <- events[,.(event_row,event_anchor=if(shared_parent && "fork_time"%in%names(events)) fork_time else NA_real_)]
  curve[anchors,on="event_row",event_anchor:=i.event_anchor]
  curve[, matched := any(focal) & any(!focal), by = .(event_row, biological_sample, mid)]
  if (isTRUE(temporal_lr$equal_sister_reference) && time_coordinate=="state_pseudotime") {
    if (is.null(temporal_lr$paired_events)) {
      sizes <- panel[,.(expected=.N),by=event_row]
      curve[sizes,on="event_row",expected:=i.expected]
    } else {
      sizes <- merge(events[,.(event_row,pair_idx,source,target,event_node,child)],
        temporal_lr$paired_events,by=c("pair_idx","source","target","event_node","child"),
        allow.cartesian=TRUE)[,.(expected=.N),by=.(event_row,biological_sample)]
      curve[sizes,on=c("event_row","biological_sample"),expected:=i.expected]
    }
    curve[,matched:=any(focal) & any(!focal) & .N==expected[1L],
      by=.(event_row,biological_sample,mid)]
  }
  curve <- switch(background, branch_conditioned = curve[matched %in% TRUE],
    focal_only = curve[focal %in% TRUE])
  switch(as.character(background=="branch_conditioned" && time_coordinate=="state_pseudotime"),"TRUE"={
    curve <- curve[is.na(event_anchor) | mid>event_anchor]
  })
  curve[, n_time := data.table::uniqueN(mid), by = event_row]
  curve <- curve[n_time >= 2L]
  switch(as.character(nrow(curve) == 0L), "TRUE" = return(out))
  # Identical support and weights share one design factorization across all LR
  # responses. No response value or known biological label defines a batch.
  curve[, context := .GRP, by = .(biological_sample, mid, next_node, mass, focal,event_anchor)]
  data.table::setorder(curve, event_row, context)
  layouts <- curve[, .(support = paste(context, collapse = ",")), by = event_row]
  layouts[, layout := .GRP, by = support]
  curve[layouts, on = "event_row", layout := i.layout]
  fits <- curve[, {
    support <- unique(.SD[, .(context, biological_sample, mid, next_node, mass, focal,event_anchor)])[order(context)]
    wide <- data.table::dcast(.SD, context ~ event_row, value.var = "activity")[order(context)]
    response <- as.matrix(wide[, -1L])
    event_rows <- as.integer(colnames(response))
    weight <- support$mass / mean(support$mass)
    if (isTRUE(temporal_lr$equal_sister_reference) && time_coordinate=="state_pseudotime") {
      support[,reference_weight:=data.table::fifelse(focal,1,1/max(1L,sum(!focal))),
        by=.(biological_sample,mid)]
      weight <- support$reference_weight/ave(support$reference_weight,support$biological_sample,FUN=sum)
      weight <- weight/mean(weight)
    }
    switch(time_coordinate, observed_calendar = {
      # Every mouse has one calendar observation; cell counts do not replicate it.
      weight <- 1 / ave(rep.int(1, nrow(support)), support$biological_sample, FUN = length)
      weight <- weight / mean(weight)
    })
    time_start <- min(support$mid)
    time_end <- max(support$mid)
    shared_anchor <- background=="branch_conditioned" && time_coordinate=="state_pseudotime" &&
      all(is.finite(support$event_anchor))
    switch(as.character(shared_anchor),"TRUE"={
      time_start <- min(time_end,unique(support$event_anchor))
    })
    time <- (support$mid - time_start) / (time_end - time_start)
    switch(as.character(shared_anchor),"TRUE"={
      time <- pmax(0,time)
    })
    focal <- support$focal
    support[, sample_child := .GRP, by = .(biological_sample, next_node)]
    support[, sample_time := .GRP, by = .(biological_sample, mid)]
    baseline <- switch(background,
      branch_conditioned = switch(as.character(condition_on_sample),
        "TRUE" = switch(time_coordinate,
          state_pseudotime=switch(as.character(shared_anchor),
            "TRUE"=outer(support$sample_time,unique(support$sample_time),"==")*1,
            "FALSE"=cbind(outer(support$sample_child,unique(support$sample_child),"=="),
              outer(support$sample_time,unique(support$sample_time),"=="))*1),
          observed_calendar=cbind(
            outer(support$sample_child,unique(support$sample_child),"=="),
            outer(support$sample_time,unique(support$sample_time),"=="))*1),
        "FALSE" = switch(as.character(shared_anchor),
          "TRUE"=stats::model.matrix(~factor(mid),support),
          "FALSE"=stats::model.matrix(~factor(next_node)+factor(mid),support))),
      focal_only = switch(as.character(condition_on_sample),
        "TRUE" = outer(support$biological_sample, unique(support$biological_sample), "==") * 1,
        "FALSE" = matrix(1, nrow = nrow(support), ncol = 1L)))
    switch(time_coordinate, observed_calendar = {
      # A sample-by-child intercept would absorb the entire calendar interaction.
      baseline <- switch(background,
        branch_conditioned = stats::model.matrix(
          ~ factor(biological_sample) + factor(next_node), support),
        focal_only = matrix(1, nrow = nrow(support), ncol = 1L))
    })
    native_shared <- background == "branch_conditioned" &&
      time_coordinate == "state_pseudotime" && shared_anchor
    level_effect <- native_shared && include_level
    sample_levels <- unique(support$biological_sample)
    if (native_shared && condition_on_sample && length(sample_levels) > 1L) {
      sample_deviation <- stats::contr.sum(length(sample_levels))[
        match(support$biological_sample, sample_levels), , drop = FALSE]
      baseline <- cbind(baseline, sample_deviation * focal)
    }
    profiled_samples <- switch(time_coordinate,
      state_pseudotime = 0L, observed_calendar = switch(background,
        branch_conditioned = data.table::uniqueN(support$biological_sample),
        focal_only = 0L))
    independent_observations <- nrow(support) - profiled_samples
    # Weighted partial regression profiles out the shared baseline once per layout.
    # The residual objective equals the full weighted least-squares objective.
    fit_tolerance <- eval(formals(stats::lm.wfit)$tol)
    sqrt_weight <- sqrt(weight)
    baseline_qr <- qr(baseline * sqrt_weight, tol = fit_tolerance, LAPACK = FALSE)
    response_residual <- qr.resid(baseline_qr, response * sqrt_weight)
    response_energy <- colSums(weight * response^2)
    observed_peak <- observed_lower <- observed_upper <- rep(NA_real_,ncol(response))
    switch(as.character(background=="focal_only"),"TRUE"={
      observed_time <- sort(unique(support$mid))
      time_basis <- outer(support$mid,observed_time[-1L],"==")*1
      time_projected <- qr.resid(baseline_qr,time_basis*sqrt_weight)
      time_fit <- stats::lm.fit(time_projected,response_residual,tol=fit_tolerance)
      switch(as.character(time_fit$rank==ncol(time_basis)),"TRUE"={
        time_effect <- rbind(0,matrix(time_fit$coefficients,nrow=ncol(time_basis)))
        peak_index <- max.col(t(time_effect),ties.method="first")
        lower_index <- pmax(1L,peak_index-1L)
        upper_index <- pmin(length(observed_time),peak_index+1L)
        observed_peak <- observed_time[peak_index]
        observed_lower <- observed_time[lower_index]
        observed_upper <- observed_time[upper_index]
        center_value <- time_effect[cbind(peak_index,seq_len(ncol(response)))]
        lower_value <- time_effect[cbind(lower_index,seq_len(ncol(response)))]
        upper_value <- time_effect[cbind(upper_index,seq_len(ncol(response)))]
        slope_left <- (center_value-lower_value)/(observed_peak-observed_lower)
        slope_right <- (upper_value-center_value)/(observed_upper-observed_peak)
        curvature <- (slope_right-slope_left)/(observed_upper-observed_lower)
        linear_term <- slope_left-curvature*(observed_lower+observed_peak)
        vertex <- -linear_term/(2*curvature)
        interior <- lower_index<peak_index & peak_index<upper_index & is.finite(vertex) &
          curvature<0 & vertex>=observed_lower & vertex<=observed_upper & time_coordinate=="state_pseudotime"
        observed_peak[interior] <- vertex[interior]
        empty_curve <- colSums(response^2)==0
        observed_peak[empty_curve] <- observed_lower[empty_curve] <- observed_upper[empty_curve] <- NA_real_
      })
    })
    knots <- data.table::data.table(knot = sort(unique(time[time > 0 & time < 1])))
    models <- data.table::rbindlist(list(data.table::data.table(
      model = c("constant", "linear", if (level_effect) "step"), knot = NA_real_),
      knots[, .(model = c("plateau", "delayed", "hinge")), by = knot]), use.names = TRUE)
    models[,knot_end:=knot]
    plateau_knots <- switch(time_coordinate,
      state_pseudotime=data.table::CJ(knot=knots$knot,knot_end=knots$knot)[knot_end>knot],
      observed_calendar=data.table::data.table(knot=numeric(),knot_end=numeric()))
    plateau_knots[,model:="trapezoid"]
    models <- data.table::rbindlist(list(models,plateau_knots),use.names=TRUE,fill=TRUE)
    models[,c("knot_lower","knot_upper","end_lower","end_upper"):=rep(list(NA_real_),4L)]
    # Breakpoint linearization (Muggeo, 2003; doi:10.1002/sim.1545).
    # Within one observed-time gap, the fitted line intersections give the
    # continuous knots without an iterative optimizer or a finer search grid.
    continuous_models <- switch(time_coordinate,observed_calendar=NULL,state_pseudotime={
      observed_times <- sort(unique(time))
      gaps <- data.table::data.table(knot_lower=head(observed_times,-1L),knot_upper=tail(observed_times,-1L))
      gaps[,knot:=(knot_lower+knot_upper)/2]
      hinged <- gaps[,.(model="hinge_continuous",knot,knot_end=knot,knot_lower,knot_upper,
        end_lower=knot_lower,end_upper=knot_upper)]
      flat_top <- merge(gaps[,.(knot,knot_lower,knot_upper,join_key=1L)],
        gaps[,.(knot_end=knot,end_lower=knot_lower,end_upper=knot_upper,join_key=1L)],
        by="join_key",allow.cartesian=TRUE)[knot_end>knot]
      flat_top[,`:=`(model="trapezoid_continuous",join_key=NULL)]
      data.table::rbindlist(list(hinged,flat_top),use.names=TRUE,fill=TRUE)
    })
    models <- data.table::rbindlist(list(models,continuous_models),use.names=TRUE,fill=TRUE)
    models[, c("transition_lo", "off_lo", "off_hi") := rep(list(NA_real_), 3L)]
    interval_models <- switch(time_coordinate, state_pseudotime = NULL,
      observed_calendar = {
        ordered_times <- sort(unique(time))
        intervals <- data.table::data.table(
          transition_lo = head(ordered_times, -1L), knot = tail(ordered_times, -1L))
        pulses <- data.table::CJ(knot = intervals$knot, off_hi = intervals$knot)[off_hi > knot]
        pulses <- merge(pulses, intervals, by = "knot")
        pulses[, off_lo := ordered_times[match(off_hi, ordered_times) - 1L]]
        data.table::rbindlist(list(
          intervals[, .(model = "interval_step", knot, transition_lo,
            off_lo = NA_real_, off_hi = NA_real_)],
          pulses[, .(model = "interval_pulse", knot, transition_lo, off_lo, off_hi)]))
      })
    models <- data.table::rbindlist(list(models, interval_models), use.names = TRUE, fill = TRUE)
    models[, {
      interaction <- switch(.BY$model,
        interval_step = matrix((time >= .BY$knot) * focal, ncol = 1L),
        interval_pulse = matrix((time >= .BY$knot & time < .BY$off_hi) * focal, ncol = 1L),
        constant = matrix(numeric(), nrow = nrow(support), ncol = 0L),
        step = matrix(numeric(), nrow = nrow(support), ncol = 0L),
        linear = matrix(time * focal, ncol = 1L),
        plateau = matrix(pmin(time, .BY$knot) * focal, ncol = 1L),
        delayed = matrix(pmax(time - .BY$knot, 0) * focal, ncol = 1L),
        trapezoid = cbind(pmin(time,.BY$knot)*focal,pmax(time-.BY$knot_end,0)*focal),
        hinge_continuous=cbind(pmin(time,.BY$knot)*focal,pmax(time-.BY$knot,0)*focal,
          (time>.BY$knot)*focal),
        trapezoid_continuous=cbind(pmin(time,.BY$knot)*focal,pmax(time-.BY$knot_end,0)*focal,
          (time>.BY$knot)*focal,(time>.BY$knot_end)*focal),
        hinge = cbind(pmin(time, .BY$knot) * focal, pmax(time - .BY$knot, 0) * focal))
      estimate_level <- level_effect && .BY$model != "constant"
      if (estimate_level) interaction <- cbind(focal, interaction)
      projected <- switch(.BY$model, constant = interaction,
        qr.resid(baseline_qr, interaction * sqrt_weight))
      fitted <- switch(.BY$model,
        constant = list(residuals = response_residual,
          coefficients = matrix(numeric(), nrow = 0L, ncol = ncol(response)), rank = 0L),
        stats::lm.fit(projected, response_residual, tol = fit_tolerance))
      residual <- matrix(fitted$residuals, nrow = nrow(support))
      coefficients <- matrix(fitted$coefficients, nrow = ncol(interaction), ncol = ncol(response))
      base_effect <- rep(0, ncol(response))
      slopes <- coefficients
      if (estimate_level) {
        base_effect <- coefficients[1L, ]
        slopes <- coefficients[-1L, , drop = FALSE]
      }
      model_rank <- baseline_qr$rank + fitted$rank - profiled_samples
      identified <- switch(.BY$model, constant = TRUE,
        all(colSums(projected^2) > fit_tolerance^2 * colSums(weight * interaction^2)))
      left <- switch(.BY$model, constant =, step =, delayed = rep(0, ncol(response)), slopes[1L, ])
      right <- switch(.BY$model, constant =, step =, plateau = rep(0, ncol(response)),
        hinge =, trapezoid =, hinge_continuous=,trapezoid_continuous=slopes[2L, ],slopes[1L, ])
      fitted_knot <- rep(.BY$knot,ncol(response))
      fitted_knot_end <- rep(.BY$knot_end,ncol(response))
      admissible <- rep(TRUE,ncol(response))
      switch(.BY$model,hinge_continuous={
        fitted_knot <- .BY$knot+slopes[3L,]/(left-right)
        fitted_knot_end <- fitted_knot
        admissible <- is.finite(fitted_knot) & fitted_knot>=.BY$knot_lower & fitted_knot<=.BY$knot_upper
      },trapezoid_continuous={
        fitted_knot <- .BY$knot+slopes[3L,]/left
        fitted_knot_end <- .BY$knot_end-slopes[4L,]/right
        admissible <- is.finite(fitted_knot) & is.finite(fitted_knot_end) &
          fitted_knot>=.BY$knot_lower & fitted_knot<=.BY$knot_upper &
          fitted_knot_end>=.BY$end_lower & fitted_knot_end<=.BY$end_upper & fitted_knot_end>fitted_knot
      })
      rss <- pmax(colSums(residual^2), .Machine$double.eps *
        pmax(response_energy, .Machine$double.eps))
      parameters <- model_rank +
        as.integer(.BY$model %in% c("plateau", "delayed", "hinge", "interval_step","trapezoid")) +
        as.integer(.BY$model=="trapezoid") +
        2L * as.integer(.BY$model == "interval_pulse")
      data.table::data.table(event_row = event_rows,
        fitted_knot=fitted_knot,fitted_knot_end=fitted_knot_end,
        dynamic_observed_peak_time=observed_peak,dynamic_observed_peak_lower=observed_lower,
        dynamic_observed_peak_upper=observed_upper,
        bic = independent_observations * log(rss / independent_observations) +
          parameters * log(independent_observations),
        left = left, right = right, base_effect = base_effect,
        n_observations = nrow(support),
        n_bins = data.table::uniqueN(time), n_samples = data.table::uniqueN(support$biological_sample),
        time_start = time_start, time_end = time_end,
        observed_time_start = min(support$mid), observed_start = min(time),
        dynamic_sample_ids=paste(unique(support$biological_sample),collapse="\034"),
        full_rank = identified & admissible &
          fitted$rank == ncol(interaction) & independent_observations > model_rank)
    }, by = .(model, knot, knot_end, knot_lower,knot_upper,end_lower,end_upper,transition_lo, off_lo, off_hi)]
  }, by = layout]
  fits <- fits[full_rank %in% TRUE]
  fits[,`:=`(knot=fitted_knot,knot_end=fitted_knot_end)]
  fits[,model:=data.table::fcase(model=="hinge_continuous","hinge",model=="trapezoid_continuous","trapezoid",default=model)]
  transitions <- merge(fits[model %in% c("interval_step", "interval_pulse")],
    fits[model == "constant", .(event_row, constant_bic = bic)], by = "event_row")
  out$dynamic_transition_status[unique(transitions$event_row)] <- "no_supported_transition"
  data.table::setorder(transitions, event_row, bic)
  transition_best <- transitions[, .SD[1L], by = event_row][constant_bic > bic]
  transition_best[, dynamic_transition_lower := time_start + transition_lo * (time_end - time_start)]
  transition_best[, dynamic_transition_upper := time_start + knot * (time_end - time_start)]
  transition_best[, dynamic_transition_effect := left]
  transition_best[, dynamic_transition_bic_gain := constant_bic - bic]
  transition_best[, dynamic_offset_lower := time_start + off_lo * (time_end - time_start)]
  transition_best[, dynamic_offset_upper := time_start + off_hi * (time_end - time_start)]
  out[transition_best$event_row, transition_fields] <- as.data.frame(
    transition_best[, transition_fields, with = FALSE])
  out$dynamic_transition_model[transition_best$event_row] <- transition_best$model
  out$dynamic_transition_status[transition_best$event_row] <- "observed_interval"
  fits <- fits[!model %in% c("interval_step", "interval_pulse")]
  switch(as.character(nrow(fits) == 0L), "TRUE" = return(out))
  fits[, pulse := model %in% c("hinge","trapezoid") & left > 0 & right < 0]
  comparison <- fits[, .(nonpulse_bic = min(bic[!pulse]),
    pulse_bic = sort(bic[pulse])[1L], constant_bic = bic[model == "constant"][1L]), by = event_row]
  comparison[, pulse_bic_gain := data.table::fcoalesce(pmax(0, nonpulse_bic - pulse_bic), 0)]
  data.table::setorder(fits, event_row, bic)
  best <- merge(fits[, .SD[1L], by = event_row], comparison, by = "event_row")
  best[, dynamic_pattern := data.table::fcase(model == "constant", "constant",
    model == "step" & base_effect > 0, "sustained_preference",
    model == "step" & base_effect < 0, "sustained_negative",
    model == "step", "constant",
    model == "plateau" & left > 0, "activation_plateau", model == "plateau", "deactivation_plateau",
    model == "delayed" & right > 0, "late_increase", model == "delayed", "late_decrease",
    left > 0 & right < 0, "pulse", left < 0 & right > 0, "trough",
    left >= 0 & right >= 0, "increasing", default = "decreasing")]
  best[, `:=`(dynamic_time_start = observed_time_start, dynamic_time_end = time_end,
    dynamic_change_time = time_start + knot * (time_end - time_start),
    dynamic_slope_left = left / (time_end - time_start),
    dynamic_slope_right = right / (time_end - time_start), dynamic_bic_gain = constant_bic - bic,
    dynamic_n_bins = n_bins, dynamic_n_observations = n_observations,
    dynamic_n_samples = n_samples)]
  best[,`:=`(dynamic_peak_lower=data.table::fifelse(dynamic_pattern=="pulse",time_start+knot*(time_end-time_start),NA_real_),
    dynamic_peak_upper=data.table::fifelse(dynamic_pattern=="pulse",time_start+knot_end*(time_end-time_start),NA_real_))]
  best[,effect_knot:=data.table::fcoalesce(knot,1)]
  best[,effect_end:=data.table::fcoalesce(knot_end,effect_knot)]
  best[,effect_at_knot:=base_effect+left*effect_knot]
  best[,effect_at_end:=effect_at_knot+right*(1-effect_end)]
  best[, effect_at_start := base_effect + left * observed_start]
  best[, effect_choice := max.col(abs(cbind(effect_at_start, effect_at_knot,
    effect_at_end)), ties.method = "first")]
  best[, dynamic_effect := cbind(effect_at_start, effect_at_knot,
    effect_at_end)[cbind(seq_len(.N), effect_choice)]]
  best[, dynamic_effect_time := cbind(observed_time_start,
    time_start + effect_knot * (time_end - time_start),
    time_end)[cbind(seq_len(.N), effect_choice)]]
  best[, dynamic_level_effect := base_effect]
  out[best$event_row, fields] <- as.data.frame(best[, fields, with = FALSE])
  out$dynamic_sample_ids[best$event_row] <- best$dynamic_sample_ids
  out$dynamic_model[best$event_row] <- best$model
  out$dynamic_pattern[best$event_row] <- best$dynamic_pattern
  out$dynamic_status[best$event_row] <- "observed"
  if (is.null(sample_contrasts)) out$cohort_effect <- out$dynamic_effect
  if (!is.null(sample_contrasts)) {
    ids <- intersect(c("unit_id","pair_idx","source","target","event_node","child"),names(out))
    local_field <- switch(background,branch_conditioned="sibling_contrast",focal_only="post_fork_change")
    native <- data.table::as.data.table(sample_contrasts)[evaluable%in%TRUE&is.finite(get(local_field)),
      c(ids,"block_label",local_field),with=FALSE]
    data.table::setnames(native,c("block_label",local_field),c("biological_sample","local_effect"))
    predicted <- data.table::copy(data.table::as.data.table(out))[,event_row:=.I]
    native <- merge(native,predicted[,c(ids,"event_row","dynamic_effect"),with=FALSE],by=ids)
    observed <- predicted[is.finite(dynamic_effect)&!is.na(dynamic_sample_ids),
      .(biological_sample=unlist(strsplit(dynamic_sample_ids,"\034",fixed=TRUE))),by=event_row]
    observed[,curve_observed:=TRUE]
    native <- merge(native,observed,by=c("event_row","biological_sample"),all.x=TRUE)
    native[,completed:=data.table::fifelse(curve_observed%in%TRUE,dynamic_effect,local_effect)]
    completed <- native[,.(effect=stats::median(completed),n_native=.N),by=event_row]
    out$cohort_effect[completed$event_row] <- completed$effect
    if (transfer_pulse) {
      native_n <- rep(0,nrow(out))
      native_n[completed$event_row] <- completed$n_native
      pulse_gain <- switch(background,branch_conditioned=out$focal_pulse_bic_gain,
        focal_only=out$pulse_bic_gain)
      pulse_n <- switch(background,branch_conditioned=out$focal_dynamic_n_samples,
        focal_only=out$dynamic_n_samples)
      weight <- pmin(1,pulse_n/pmax(1,native_n))*(-expm1(-data.table::fcoalesce(pulse_gain,0)/2))
      weight[!is.finite(weight)|!is.finite(out$dynamic_effect)|native_n==0] <- 0
      out$cohort_pulse_weight <- weight
      curve_effect <- data.table::fcoalesce(out$dynamic_effect,out$cohort_effect)
      out$cohort_effect <- out$cohort_effect+weight*(curve_effect-out$cohort_effect)
    }
  }
  out
}

# Typed node profiles retain observed sampling contexts; graph membership identifies
# resident senders separately from cells passing through an internal tree node.
ComputeResidentSenderEventsFromCells <- function(
    expr, celltypes, lineage_clusters, lineage_paths, event_units,
    knn_conn, pair_gene_map, hill_scale, stage = rep("snapshot", ncol(expr)),
    context_order = unique(stage), embedding = NULL, observed_time = NULL,
    condition = rep("all", ncol(expr)), min_cells = 5, min_pct = 0.05,
    sender_scope = "parent_and_ancestor", structural_reference = "parent_graph") {
  pair_map <- data.table::as.data.table(pair_gene_map)
  genes <- unique(unlist(strsplit(c(pair_map$ligand_genes,
    pair_map$receptor_genes), "_", fixed = TRUE)))
  cells <- data.table::data.table(cell = seq_len(ncol(expr)), node = lineage_clusters,
    type = celltypes, stage = as.character(stage))
  cells[, group := .GRP, by = .(node, type, stage)]
  membership <- Matrix::sparseMatrix(i = cells$cell, j = cells$group, x = 1)
  size <- Matrix::colSums(membership)
  selected <- expr[genes, , drop = FALSE]
  means <- sweep(as.matrix(selected %*% membership), 2L, size, "/")
  detected <- sweep(as.matrix((selected > 0) %*% membership), 2L, size, "/")
  groups <- unique(cells[, .(group, node, type, stage)])
  groups[, cells := size[group]]
  profile <- data.table::data.table(group = rep(seq_along(size), each = length(genes)),
    gene = rep(genes, length(size)), expression = as.numeric(means),
    detected_fraction = as.numeric(detected))
  profile <- merge(profile, groups, by = "group")
  sides <- data.table::rbindlist(lapply(event_units$units, "[",
    c("event_node", "children", "desc")), idcol = "unit_id")
  data.table::setnames(sides, "children", "child")
  tree <- data.table::data.table(lineage = rep(seq_along(lineage_paths), lengths(lineage_paths)),
    node = unlist(lineage_paths), position = sequence(lengths(lineage_paths)))
  side_lineages <- data.table::data.table(unit_id = rep(sides$unit_id, lengths(sides$desc)),
    child = rep(sides$child, lengths(sides$desc)), lineage = unlist(sides$desc))
  side_lineages <- merge(side_lineages, tree[, .(lineage, child = node, child_position = position)],
    by = c("lineage", "child"))
  side_nodes <- unique(merge(side_lineages, tree, by = "lineage", allow.cartesian = TRUE)[
    position >= child_position, .(unit_id, child, node, scope = "branch")])
  side_nodes <- data.table::rbindlist(list(side_nodes,
    unique(sides[, .(unit_id, child = event_node, node = event_node, scope = "parent")])))
  members <- merge(side_nodes, cells[, .(cell, node, target = type)], by = "node", allow.cartesian = TRUE)
  members[, support_group := .GRP, by = .(scope, unit_id, child, target)]
  support_groups <- unique(members[, .(support_group, scope, unit_id, child, target)])
  support_membership <- Matrix::sparseMatrix(i = members$cell, j = members$support_group,
    x = 1, dims = c(ncol(expr), nrow(support_groups)))
  type_membership <- Matrix::sparseMatrix(i = cells$cell, j = cells$type, x = 1)
  same_type_neighbors <- as.matrix(knn_conn %*% type_membership)
  members[, same_type_edges := same_type_neighbors[cbind(cell, target)]]
  support <- members[, .(context_cells = .N, same_type_edges = sum(same_type_edges)), by = support_group]
  data.table::setorder(support, support_group)
  support[, internal_edges := as.numeric(Matrix::colSums(
    support_membership * (knn_conn %*% support_membership)))]
  support <- merge(support_groups, support, by = "support_group")
  structural_history <- NULL
  switch(as.character(is.null(observed_time) &&
    structural_reference %in% c("parent_graph", "parent_nearest")),
    "TRUE" = {
      branches <- data.table::rbindlist(lapply(event_units$units, "[",
        c("event_node", "segment_post_nodes")))
      parent_nodes <- unique(data.table::data.table(
        receiver_event = rep(branches$event_node, lengths(branches$segment_post_nodes)),
        receiver_node = unlist(branches$segment_post_nodes)))
      structural_history <- merge(parent_nodes, groups[cells >= min_cells,
        .(receiver_node = node, query_group = group, stage)], by = "receiver_node")
      structural_history <- merge(structural_history, groups[cells >= min_cells,
        .(receiver_event = node, pre_group = group, stage)],
        by = c("receiver_event", "stage"), allow.cartesian = TRUE)
      switch(structural_reference,
        parent_graph = {
          connections <- as.matrix(Matrix::crossprod(membership,
            ((knn_conn + Matrix::t(knn_conn)) / 2) %*% membership))
          structural_history[, weight := connections[cbind(query_group, pre_group)]]
          structural_history <- structural_history[weight > 0]
          structural_history[, weight := weight / sum(weight), by = .(receiver_event, query_group)]
        },
        parent_nearest = {
          centers <- as.matrix(Matrix::crossprod(membership, embedding)) / groups$cells
          residual <- embedding - centers[cells$group, , drop = FALSE]
          precision <- solve(crossprod(residual) / (nrow(embedding) - nrow(groups)))
          delta <- centers[structural_history$query_group, , drop = FALSE] -
            centers[structural_history$pre_group, , drop = FALSE]
          structural_history[, distance := rowSums((delta %*% precision) * delta)]
          data.table::setorder(structural_history, receiver_event, query_group, distance, pre_group)
          structural_history <- structural_history[, .SD[1L], by = .(receiver_event, query_group)]
          structural_history[, weight := 1]
        })
    })
  history <- NULL
  switch(as.character(!is.null(observed_time)), "TRUE" = {
    context <- unique(data.table::data.table(stage = as.character(stage),
      time = observed_time, condition = condition))
    groups <- merge(groups, context, by = "stage")
    data.table::setorder(groups, group)
    centers <- as.matrix(Matrix::crossprod(membership, embedding)) / groups$cells
    residual <- embedding - centers[cells$group, , drop = FALSE]
    precision <- solve(crossprod(residual) / (nrow(embedding) - nrow(groups)))
    query <- groups[cells >= min_cells,
      .(query_group = group, query_time = time, query_type = type, condition)]
    previous <- groups[cells >= min_cells,
      .(pre_group = group, pre_time = time, pre_type = type, condition)]
    history <- merge(query, previous, by = "condition", allow.cartesian = TRUE)[pre_time < query_time]
    history[, same_type := query_type == pre_type]
    history[, identity_observed := any(same_type), by = query_group]
    history <- history[!identity_observed | same_type]
    history <- history[, .SD[pre_time == max(pre_time)], by = query_group]
    delta <- centers[history$query_group, , drop = FALSE] - centers[history$pre_group, , drop = FALSE]
    history[, distance := rowSums((delta %*% precision) * delta)]
    groups[, population := .GRP, by = .(type, time, condition)]
    populations <- groups[, .(cells = sum(cells)), by = population]
    population_membership <- Matrix::sparseMatrix(i = groups$group,
      j = groups$population, x = groups$cells)
    population_centers <- as.matrix(Matrix::crossprod(population_membership, centers)) / populations$cells
    population_delta <- population_centers[groups$population[history$query_group], , drop = FALSE] -
      population_centers[groups$population[history$pre_group], , drop = FALSE]
    history[, population_distance := rowSums((population_delta %*% precision) * population_delta)]
    data.table::setorder(history, query_group, -same_type, population_distance, distance, pre_group)
    history <- history[, .SD[1L], by = query_group]
  })
  events <- ComputeResidentSenderEvents(profile, pair_map, event_units,
    support[scope == "branch"], support[scope == "parent"], hill_scale, as.character(context_order),
    min_cells = min_cells, min_pct = min_pct, lineage_paths = lineage_paths,
    sender_scope = sender_scope, counterfactual = "fixed_composition", history_matches = history,
    structural_matches = structural_history)
  attr(events, "history_matches") <- history
  attr(events, "structural_matches") <- structural_history
  events
}

ComputeResidentSenderEvents <- function(
    node_profiles, pair_gene_map, event_units, branch_residency,
    parent_residency, hill_scale, context_order,
    min_cells = 5, min_pct = 0.05, lineage_paths = NULL,
    sender_scope = "parent", counterfactual = "original",
    history_matches = NULL, structural_matches = NULL) {
  profile <- data.table::as.data.table(node_profiles)
  pair_map <- data.table::as.data.table(pair_gene_map)
  branch_support <- data.table::as.data.table(branch_residency)
  parent_support <- data.table::as.data.table(parent_residency)
  sides <- data.table::rbindlist(lapply(event_units$units, "[",
    c("event_node", "children", "segment_post_nodes")), idcol = "unit_id")
  data.table::setnames(sides, "children", "child")
  sides[, side_id := seq_len(.N)]
  post_nodes <- data.table::data.table(
    unit_id = rep(sides$unit_id, lengths(sides$segment_post_nodes)),
    child = rep(sides$child, lengths(sides$segment_post_nodes)),
    node = unlist(sides$segment_post_nodes))
  candidates <- merge(sides[, .(unit_id, sender_unit = unit_id,
    receiver_event = event_node, sender_child = event_node,
    receiver_child = child)], post_nodes, by.x = c("unit_id", "receiver_child"),
    by.y = c("unit_id", "child"))
  data.table::setnames(candidates, "node", "receiver_node")
  switch(sender_scope, ancestor =, parent_and_ancestor = {
    paths <- data.table::data.table(
      lineage = rep(seq_along(lineage_paths), lengths(lineage_paths)),
      node = unlist(lineage_paths), step = sequence(lengths(lineage_paths)))
    ancestry <- unique(merge(paths, paths, by = "lineage", allow.cartesian = TRUE)[
      step.x < step.y | (sender_scope == "parent_and_ancestor" & step.x == step.y),
      .(sender_child = node.x, receiver_event = node.y)])
    senders <- unique(sides[, .(sender_unit = unit_id, sender_child = event_node)])
    candidates[, c("sender_unit", "sender_child") := NULL]
    candidates <- merge(candidates, merge(senders, ancestry, by = "sender_child",
      allow.cartesian = TRUE), by = "receiver_event", allow.cartesian = TRUE)
  })
  gene_map <- data.table::rbindlist(list(
    pair_map[, .(gene = unlist(strsplit(ligand_genes, "_", fixed = TRUE)),
      role = "L"), by = pair_idx],
    pair_map[, .(gene = unlist(strsplit(receptor_genes, "_", fixed = TRUE)),
      role = "R"), by = pair_idx]))
  lr <- merge(profile[cells >= min_cells], gene_map, by = "gene",
    allow.cartesian = TRUE)[, .(expression = min(expression),
      pct = min(detected_fraction)), by = .(node, stage, type, cells, pair_idx, role)]
  lr <- data.table::dcast(lr, node + stage + type + cells + pair_idx ~ role,
    value.var = c("expression", "pct"))
  events <- merge(candidates, lr[, .(sender_child = node, stage, source = type,
    source_cells = cells, pair_idx, L_sender = expression_L, L_sender_pct = pct_L)],
    by = "sender_child", allow.cartesian = TRUE)
  events <- merge(events, lr[, .(receiver_node = node, stage, target = type,
    receiver_cells = cells, pair_idx, R_receiver = expression_R, R_receiver_pct = pct_R)],
    by = c("receiver_node", "stage", "pair_idx"), allow.cartesian = TRUE)
  events[, `:=`(event_id = .I, post_stage_index = match(stage, context_order),
    source_supported = FALSE, receiver_supported = FALSE)]
  events[parent_support, on = .(sender_unit = unit_id, sender_child = child, source = target),
    source_supported := 2 * i.internal_edges > i.same_type_edges]
  events[branch_support, on = .(unit_id, receiver_child = child, target),
    receiver_supported := 2 * i.internal_edges > i.same_type_edges]
  pre_context <- unique(profile[cells >= min_cells, .(node, stage, type, cells)])[
    , .(pre_cells = sum(cells)), by = .(receiver_event = node, pre_stage = stage)]
  pre_context[, pre_stage_index := match(pre_stage, context_order)]
  pre <- merge(events[, .(event_id, pair_idx, sender_child, receiver_event, source, target,
    post_stage_index)], pre_context, by = "receiver_event", allow.cartesian = TRUE)[
      pre_stage_index <= post_stage_index]
  switch(sender_scope, ancestor =, parent_and_ancestor = {
    source_contexts <- unique(profile[cells >= min_cells,
      .(sender_child = node, pre_stage = stage)])
    pre <- merge(pre, source_contexts, by = c("sender_child", "pre_stage"))
  })
  pre[, fraction := pre_cells / sum(pre_cells), by = event_id]
  requests <- data.table::rbindlist(list(
    pre[, .(event_id, pair_idx, phase = "pre", stage = pre_stage, fraction,
      role = "pre_L", node = sender_child, query_type = source)],
    pre[, .(event_id, pair_idx, phase = "pre", stage = pre_stage, fraction,
      role = "pre_R", node = receiver_event, query_type = target)],
    events[, .(event_id, pair_idx, phase = "post", stage, fraction = 1,
      role = "receiver_L", node = receiver_node, query_type = source)]))
  controls <- merge(requests, lr, by = c("node", "stage", "pair_idx"),
    allow.cartesian = TRUE)
  controls[, expression := data.table::fifelse(role == "pre_R",
    expression_R, expression_L)]
  alternatives <- merge(post_nodes,
    branch_support[2 * internal_edges > same_type_edges,
      .(unit_id, child, type = target)], by = c("unit_id", "child"),
    allow.cartesian = TRUE)
  alternative_requests <- merge(events[, .(event_id, unit_id, receiver_child,
    stage, pair_idx, target)], alternatives, by = "unit_id", allow.cartesian = TRUE)[
      child != receiver_child]
  alternative_controls <- merge(alternative_requests, lr,
    by = c("node", "stage", "pair_idx", "type"), allow.cartesian = TRUE)
  controls <- data.table::rbindlist(list(
    controls[, .(event_id, pair_idx, phase, stage, fraction, role, node,
      query_type, type, cells, expression)],
    alternative_controls[, .(event_id, pair_idx, phase = "post", stage,
      fraction = 1, role = "receiver_R", node, query_type = target, type, cells,
      expression = expression_R)]))
  controls[, same_type := query_type == type]
  data.table::setorder(controls, event_id, phase, stage, role, -same_type,
    -expression, node, type)
  controls <- controls[, .SD[1L], by = .(event_id, phase, stage, role)]
  pre_values <- data.table::dcast(controls[phase == "pre"],
    event_id + stage + fraction ~ role, value.var = "expression")
  pre_values[, geo := sqrt(pre_L * pre_R)]
  pre_values[, pre_score := geo / (hill_scale + geo)]
  pre_summary <- pre_values[, .(pre_score = sum(fraction * pre_score),
    pre_coverage = sum(fraction[is.finite(pre_score)])), by = event_id]
  post_values <- data.table::dcast(controls[phase == "post"], event_id ~ role,
    value.var = "expression")
  events <- merge(events, pre_summary, by = "event_id", all.x = TRUE)
  events <- merge(events, post_values, by = "event_id", all.x = TRUE)
  events[, `:=`(post_geo = sqrt(L_sender * R_receiver),
    source_control_geo = sqrt(receiver_L * R_receiver),
    receiver_control_geo = sqrt(L_sender * receiver_R))]
  events[, `:=`(post_score = post_geo / (hill_scale + post_geo),
    source_control_score = source_control_geo / (hill_scale + source_control_geo),
    receiver_control_score = receiver_control_geo / (hill_scale + receiver_control_geo))]
  switch(counterfactual, fixed_composition = {
    # Compare molecular changes at the focal population composition.
    groups <- unique(profile[, .(group, node, stage, type, cells)])
    groups[, fraction := cells / sum(cells), by = .(node, stage)]
    events[groups, on = .(sender_child = node, stage, source = type),
      `:=`(source_group = i.group, source_fraction = i.fraction)]
    events[groups, on = .(receiver_node = node, stage, target = type),
      `:=`(target_group = i.group, target_fraction = i.fraction)]
    fields <- merge(profile, gene_map, by = "gene", allow.cartesian = TRUE)[,
      .(expression = min(expression)), by = .(group, node, stage, type, pair_idx, role)]
    fields <- data.table::dcast(fields, group + node + stage + type + pair_idx ~ role,
      value.var = "expression")
    prior <- switch(as.character(is.null(history_matches)),
      "TRUE" = merge(events[, .(event_id, receiver_event, stage, target,
        pair_idx, history_L = L_sender)],
        fields[, .(receiver_event = node, stage, target = type, pair_idx, history_R = R)],
        by = c("receiver_event", "stage", "target", "pair_idx"), all.x = TRUE),
      "FALSE" = {
        history <- merge(data.table::as.data.table(history_matches), fields,
          by.x = "pre_group", by.y = "group", allow.cartesian = TRUE)
        prior <- merge(events[, .(event_id, pair_idx, source_group, target_group)],
          history[, .(source_group = query_group, pair_idx, history_L = L)],
          by = c("source_group", "pair_idx"), all.x = TRUE)
        merge(prior, history[, .(target_group = query_group, pair_idx, history_R = R)],
          by = c("target_group", "pair_idx"), all.x = TRUE)
      })
    events[prior, on = "event_id", `:=`(history_L = i.history_L, history_R = i.history_R)]
    switch(as.character(!is.null(structural_matches)), "TRUE" = {
      structural_prior <- merge(data.table::as.data.table(structural_matches),
        fields[, .(pre_group = group, pair_idx, reference_R = R)],
        by = "pre_group", allow.cartesian = TRUE)[,
          .(reference_R = sum(weight * reference_R)),
          by = .(receiver_event, target_group = query_group, pair_idx)]
      missing_prior <- merge(events[!is.finite(history_R),
        .(event_id, receiver_event, target_group, pair_idx)], structural_prior,
        by = c("receiver_event", "target_group", "pair_idx"))
      events[missing_prior, on = "event_id", history_R := i.reference_R]
    })
    source_reference <- merge(events[, .(event_id, receiver_node, stage, source, pair_idx)],
      fields[, .(receiver_node = node, stage, source = type, pair_idx, reference_L = L)],
      by = c("receiver_node", "stage", "source", "pair_idx"), all.x = TRUE)
    source_reference[is.na(reference_L), reference_L := 0]
    receptor_reference <- merge(events[, .(event_id, unit_id, receiver_child, stage,
      target, pair_idx)], alternatives[, .(unit_id, child, node, target = type)],
      by = c("unit_id", "target"), allow.cartesian = TRUE)[
        child != receiver_child]
    receptor_reference <- merge(receptor_reference,
      fields[, .(node, stage, target = type, pair_idx, reference_R = R)],
      by = c("node", "stage", "target", "pair_idx"))
    receptor_reference <- receptor_reference[, .(reference_R = max(reference_R)), by = event_id]
    events[, `:=`(receiver_L = 0, receiver_R = 0)]
    events[source_reference, on = "event_id", receiver_L := i.reference_L]
    events[receptor_reference, on = "event_id", receiver_R := i.reference_R]
    events[, composition_scale := source_fraction * target_fraction]
    events[, `:=`(post_geo = sqrt(L_sender * R_receiver * composition_scale),
      pre_geo = sqrt(history_L * history_R * composition_scale),
      source_control_geo = sqrt(receiver_L * R_receiver * composition_scale),
      receiver_control_geo = sqrt(L_sender * receiver_R * composition_scale))]
    events[, `:=`(post_score = post_geo / (hill_scale + post_geo),
      pre_score = pre_geo / (hill_scale + pre_geo),
      pre_coverage = as.numeric(is.finite(pre_geo)),
      source_control_score = source_control_geo / (hill_scale + source_control_geo),
      receiver_control_score = receiver_control_geo / (hill_scale + receiver_control_geo))]
  })
  events[, `:=`(induction = post_score - pre_score,
    sender_specificity = post_score - source_control_score,
    receiver_specificity = post_score - receiver_control_score)]
  events <- merge(events, pair_map, by = "pair_idx")
  # DLL3 is not a trans activator of Notch: doi:10.1083/jcb.200503113.
  events[, trans_eligible := !(toupper(ligand_id) == "DLL3" &
    toupper(receptor_id) %in% paste0("NOTCH", 1:4))]
  events[, evaluable := source_supported & receiver_supported & trans_eligible &
    pre_coverage >= 1 - 1e-10 & L_sender_pct >= min_pct & R_receiver_pct >= min_pct &
    is.finite(induction) & is.finite(sender_specificity) & is.finite(receiver_specificity)]
  events[is.na(evaluable), evaluable := FALSE]
  events[, cross_positive_score := data.table::fifelse(evaluable,
    pmax(0, pmin(induction, sender_specificity, receiver_specificity)), NA_real_)]
  events[evaluable %in% TRUE, resource_rank := data.table::frank(-cross_positive_score,
    ties.method = "min"), by = .(unit_id, sender_child, receiver_child, receiver_node,
      stage, source, target)]
  branch_mass <- branch_support[, .(side_cells = sum(context_cells)), by = .(unit_id, child)]
  branch_mass[, event_cells := sum(side_cells), by = unit_id]
  events <- merge(events, branch_mass, by.x = c("unit_id", "receiver_child"),
    by.y = c("unit_id", "child"))
  events[, structural_multiplier := sqrt(pmin(side_cells,
    event_cells - side_cells) / event_cells)]
  events[, `:=`(event_ranking_score = cross_positive_score * structural_multiplier,
    channel = paste0("resident_", sender_scope))]
  data.table::setorder(events, unit_id, receiver_child, receiver_node, stage,
    source, target, pair_idx)
  events
}

# Learn receiver response programmes from RNA-only event activation. Event
# prototypes retain separate expression states instead of averaging their
# response directions. The caller supplies normalized RNA with validation
# genes excluded. State-continuous programmes retain activity magnitudes and
# estimate within-context response directions separately for each RNA state.
ComputeEventResponseProgramme <- function(expr, event_activity, context,
                                           training_cells = seq_len(ncol(expr)),
                                           event_state = NULL,
                                           representation = c("linear", "event_prototypes", "state_continuous")) {
  expr <- .prepare_lr_expr(expr)
  representation <- match.arg(representation)
  reference <- expr[, training_cells, drop = FALSE]
  active <- as.numeric(event_activity[training_cells] > 0)
  context_id <- match(context[training_cells], unique(context[training_cells]))
  model <- switch(representation,
    linear = {
      design <- Matrix::sparseMatrix(i = seq_along(training_cells), j = context_id, x = 1)
      context_mean <- as.numeric(Matrix::crossprod(design, active)) / Matrix::colSums(design)
      contrast <- active - context_mean[context_id]
      contrast_energy <- sum(contrast^2)
      loading <- as.numeric(reference %*% contrast) *
        data.table::fifelse(contrast_energy > 0, 1 / contrast_energy, 0)
      list(loading = stats::setNames(loading, rownames(reference)),
        gene_mean = stats::setNames(Matrix::rowMeans(reference), rownames(reference)))
    },
    state_continuous = {
      stopifnot(length(event_state) == ncol(expr))
      states <- unique(as.character(event_state[training_cells]))
      members <- data.table::data.table(cell = seq_along(training_cells),
        state = match(as.character(event_state[training_cells]), states),
        context = context_id, activity = as.numeric(event_activity[training_cells]))
      members[, population := .GRP, by = .(state, context)]
      members[, contrast := activity - mean(activity), by = population]
      groups <- members[, .(state = state[1L], cells = .N), by = population][order(population)]
      state_design <- Matrix::sparseMatrix(i = members$cell, j = members$state, x = 1,
        dims = c(ncol(reference), length(states)))
      context_design <- Matrix::sparseMatrix(i = members$cell, j = members$population, x = 1)
      context_to_state <- Matrix::sparseMatrix(i = groups$population, j = groups$state, x = 1,
        dims = c(nrow(groups), length(states)))
      contrasts <- Matrix::Diagonal(x = members$contrast) %*% state_design
      energy <- as.numeric(Matrix::colSums(contrasts^2))
      slopes <- sweep(as.matrix(reference %*% contrasts), 2L,
        data.table::fifelse(energy > 0, energy, 1), "/")
      context_sum <- as.matrix(reference %*% context_design)
      context_ss <- as.matrix((reference^2) %*% context_design) -
        sweep(context_sum^2, 2L, groups$cells, "/")
      state_cells <- as.numeric(Matrix::colSums(state_design))
      variance <- sweep(as.matrix(pmax(context_ss, 0) %*% context_to_state), 2L, state_cells, "/")
      weights <- slopes / pmax(variance, .Machine$double.eps)
      weights[variance <= .Machine$double.eps] <- 0
      normalization <- colSums(weights * slopes)
      weights <- sweep(weights, 2L, data.table::fifelse(normalization > 0, normalization, 1), "/")
      centroids <- sweep(as.matrix(reference %*% state_design), 2L, state_cells, "/")
      colnames(weights) <- colnames(centroids) <- states
      list(state_weights = weights, state_centroids = centroids,
        state_info = data.table::data.table(state = states, cells = state_cells, contrast_energy = energy))
    },
    event_prototypes = {
      stopifnot(length(event_state) == ncol(expr))
      members <- data.table::data.table(cell = seq_along(training_cells),
        context = as.character(context[training_cells]),
        event_state = as.character(event_state[training_cells]), active = active)
      members[, population := .GRP, by = .(context, event_state, active)]
      groups <- members[, .(cells = .N), by = .(population, context, event_state, active)]
      supported <- groups[, .(positive = sum(cells * active),
        negative = sum(cells * (1 - active))), by = context][positive > 0 & negative > 0, context]
      positive <- groups[active > 0 & context %in% supported]
      negative <- groups[active == 0 & context %in% supported]
      membership <- Matrix::sparseMatrix(i = members$cell, j = members$population, x = 1,
        dims = c(ncol(reference), nrow(groups)))
      centroids <- sweep(as.matrix(reference %*% membership), 2L, groups$cells, "/")
      background_context <- unique(positive$context)
      negative[, weight := cells / sum(cells), by = context]
      background <- Matrix::sparseMatrix(i = negative$population,
        j = match(negative$context, background_context), x = negative$weight,
        dims = c(nrow(groups), length(background_context)))
      mu_positive <- centroids[, positive$population, drop = FALSE]
      mu_negative <- as.matrix(centroids %*% background)
      prototype_names <- sprintf("event_%d", seq_len(ncol(mu_positive)))
      colnames(mu_positive) <- prototype_names
      colnames(mu_negative) <- background_context
      positive[, prototype := prototype_names]
      list(positive_centroids = mu_positive, negative_centroids = mu_negative,
        prototype_info = positive[, .(prototype, context, event_state, cells)],
        background_context = background_context)
    })
  model$representation <- representation
  model$genes <- rownames(reference)
  response <- ProjectEventResponseProgramme(model, expr, event_state = event_state)
  c(list(response = response), model)
}

# Apply a fitted programme to independent cells, aligning shared genes and
# retaining the fitted expression centroids. Query labels are not required.
ProjectEventResponseProgramme <- function(model, expr, event_state = NULL) {
  expr <- .prepare_lr_expr(expr)
  shared <- intersect(model$genes, rownames(expr))
  stopifnot(length(shared) > 0L)
  query <- expr[shared, , drop = FALSE]
  response <- switch(model$representation,
    linear = {
      loading <- model$loading[shared]
      loading_energy <- sum(loading^2)
      (as.numeric(Matrix::crossprod(loading, query)) -
        sum(loading * model$gene_mean[shared])) *
        data.table::fifelse(loading_energy > 0, 1 / loading_energy, 0)
    },
    state_continuous = {
      stopifnot(length(event_state) == ncol(query))
      weights <- model$state_weights[shared, , drop = FALSE]
      scores <- sweep(as.matrix(Matrix::crossprod(weights, query)), 1L,
        colSums(weights * model$state_centroids[shared, , drop = FALSE]), "-")
      state_index <- match(as.character(event_state), colnames(weights))
      scores[cbind(state_index, seq_len(ncol(query)))]
    },
    event_prototypes = switch(as.character(ncol(model$positive_centroids) > 0L &&
      ncol(model$negative_centroids) > 0L),
      "FALSE" = rep(NA_real_, ncol(query)),
      "TRUE" = {
        mu_positive <- model$positive_centroids[shared, , drop = FALSE]
        mu_negative <- model$negative_centroids[shared, , drop = FALSE]
        positive_logits <- sweep(as.matrix(Matrix::crossprod(mu_positive, query)),
          1L, colSums(mu_positive^2) / 2, "-")
        negative_logits <- sweep(as.matrix(Matrix::crossprod(mu_negative, query)),
          1L, colSums(mu_negative^2) / 2, "-")
        positive_index <- max.col(t(positive_logits), ties.method = "first")
        negative_index <- max.col(t(negative_logits), ties.method = "first")
        positive_logits[cbind(positive_index, seq_len(ncol(query)))] -
          negative_logits[cbind(negative_index, seq_len(ncol(query)))]
      }))
  stats::setNames(response, colnames(expr))
}

#' Observed calendar LR profiles for native branch subtrees
#'
#' Each native child is represented by its observed descendant cells. Membership
#' overlaps only between nested subtrees; it is not a future-fate assignment.
#' Expression is pooled within a subtree before the nonlinear LR transform.
#' @param expr,celltypes,lr_db,min_pct,min_n_eff See ComputeTemporalLR().
#' @param observed_time Numeric sampling time for every cell.
#' @param sample_labels Observation or biological-sample labels, one time per label.
#' @param lineage_clusters Native tree-node label for every cell.
#' @param lineage_paths Root-to-leaf native tree paths.
#' @param hill_scale Fixed pooled Hill scale.
#' @param requested_contexts Optional pair/source/target requests.
#' @param condition_labels Experimental-condition labels for every cell.
#' @return A temporal-LR list with native subtree paths and observed support.
#' @export
ComputeObservedCalendarLR <- function(
    expr, celltypes, lr_db, observed_time, sample_labels,
    lineage_clusters, lineage_paths, hill_scale, min_pct = 0.05,
    min_n_eff = 5, requested_contexts = NULL,
    condition_labels = rep.int("all", ncol(expr)),
    pooled_celltypes = integer(),activity_mode=c("thresholded","continuous")) {
  activity_mode <- match.arg(activity_mode)
  profile_min_pct <- switch(activity_mode,thresholded=min_pct,continuous=0)
  samples <- unique(data.table::data.table(
    biological_sample = as.character(sample_labels),
    observed_time = as.numeric(observed_time), condition = as.character(condition_labels)))
  stopifnot(nrow(samples) == data.table::uniqueN(sample_labels),
    all(is.finite(samples$observed_time)))
  times <- sort(unique(samples$observed_time))
  stopifnot(length(times) >= 2L)
  grid <- (times - min(times)) / diff(range(times))
  coordinates <- (observed_time - min(times)) / diff(range(times))
  bounds <- c(0, (head(grid, -1L) + tail(grid, -1L)) / 2, 1)
  tree <- data.table::data.table(
    lineage = rep(seq_along(lineage_paths), lengths(lineage_paths)),
    node = unlist(lineage_paths), position = sequence(lengths(lineage_paths)))
  nodes <- sort(unique(tree$node))
  prefix_end <- unique(tree[, .(state = node, lineage, end_position = position)], by = "state")
  prefix <- merge(prefix_end, tree, by = "lineage", allow.cartesian = TRUE)[position <= end_position]
  data.table::setorder(prefix, state, position)
  native_paths <- prefix[, .(path = list(node)), by = state][match(nodes, state), path]
  membership <- merge(data.table::data.table(cell = seq_len(ncol(expr)), state = lineage_clusters),
    prefix[, .(state, node)], by = "state", allow.cartesian = TRUE)
  subtree_weights <- Matrix::sparseMatrix(i = membership$cell, j = match(membership$node, nodes),
    x = 1, dims = c(ncol(expr), length(nodes)))
  scores <- samples[, {
    keep <- as.character(sample_labels) == .BY$biological_sample
    sample_bin <- match(.BY$observed_time, times)
    observed_types <- unique(celltypes[keep])
    requested <- switch(as.character(is.null(requested_contexts)), "TRUE" = NULL,
      "FALSE" = data.table::as.data.table(requested_contexts)[
        source %in% observed_types & target %in% observed_types])
    weights <- subtree_weights[keep, , drop = FALSE]
    pooled <- which(celltypes[keep]%in%pooled_celltypes)
    weights[pooled,] <- rep.int(1,length(pooled)*ncol(weights))
    data.table::as.data.table(ComputeTemporalLRFixedBounds(
      expr[, keep, drop = FALSE], celltypes[keep], lr_db, coordinates[keep],
      weights, bin_bounds = rep(list(bounds), length(nodes)),
      hill_scale = hill_scale, min_pct = profile_min_pct, min_n_eff = min_n_eff,
      requested_contexts = requested, return_complete_grid = TRUE,
      normalize_lineage_weights = FALSE))[bin == sample_bin]
  }, by = .(biological_sample, observed_time, condition)]
  list(scores = scores, observations = samples, requested_contexts = requested_contexts,
    intermediate = list(
    bin_bounds = rep(list(bounds), length(nodes)), hill_scale = hill_scale,
    n_bins = length(times)), lineage_paths = native_paths, profile_nodes = nodes,
    observation_scope = "observed_native_branch_subtrees",pooled_celltypes=pooled_celltypes,
    activity_mode=activity_mode,expression_threshold=profile_min_pct)
}

# Average complete within-sample event estimates. Each sample uses the same
# locally selected sister identities, with one vote per sister and per sample.
.paired_event_amplitudes <- function(fork_contrasts, sample_contrasts, sample_levels,
                                     signed = TRUE, measurement_error = TRUE) {
  out <- fork_contrasts
  ids <- c("unit_id", "pair_idx", "source", "target", "event_node", "child")
  local <- c("focal_pre", "focal_post", "other_pre", "other_post",
    "post_fork_change", "sibling_change", "sibling_contrast", "persistent_positive_score")
  late <- c("late_contrast", "persistence_focal_activity", "persistence_other_activity",
    "persistence_change_contrast", "late_focal_activity_native", "late_other_activity_native",
    "late_change_focal", "late_change_other", "late_change_contrast")
  clear <- c(local,late,"theta","focal_late","focal_late_change","theta_measurement_se",
    "persistence_phase_start","persistence_phase_end","late_phase_start","late_phase_end")
  out[,clear] <- rep(list(rep(NA_real_,nrow(out))), length(clear))
  out$n_paired_samples <- out$n_paired_late_samples <- out$n_late_bins <-
    out$n_late_phase_bins <- integer(nrow(out))
  out$n_observed_samples <- length(sample_levels)
  out$planned_sibling_children <- out$sibling_children
  out$paired_sibling_children <- rep(NA_character_,nrow(out))
  out$n_full_panel_samples <- integer(nrow(out))
  out$paired_sample_fraction <- 0
  out$amplitude_unit <- "equal_weight_paired_samples"
  out$paired_sample_status <- out$reason <- "no_paired_sample"
  out$evaluable <- FALSE
  out$persistence_support <- out$phase_pattern <- "not_evaluable"
  if (!nrow(sample_contrasts)) return(out)
  s <- data.table::copy(data.table::as.data.table(sample_contrasts))[
    evaluable %in% TRUE & is.finite(sibling_contrast)]
  if (!nrow(s)) return(out)
  ref <- data.table::as.data.table(fork_contrasts)[,c(ids,"sibling_children"),with=FALSE]
  s <- merge(s,ref,by=ids,suffixes=c("","_reference"),sort=FALSE)
  # Freeze the observable focal/sister panel within each donor. A different
  # donor can have a different supported panel; both phases use the same
  # donor-specific identities and weights, and availability is reported.
  s <- s[!is.na(sibling_children) & nzchar(sibling_children)]
  if (!nrow(s)) return(out)
  stopifnot(!anyDuplicated(s[,c(ids,"block_label"),with=FALSE]))
  s[, observed_late := is.finite(persistence_change_contrast) & n_late_bins >= 1L]
  s[, focal_late_change := persistence_focal_activity - focal_pre]
  s[, focal_late := data.table::fifelse(observed_late,
    sign(post_fork_change)*pmin(abs(post_fork_change),
      pmax(0,sign(post_fork_change)*focal_late_change)),post_fork_change)]
  if (!signed) s[, focal_late := data.table::fifelse(observed_late & post_fork_change>0,
    pmin(post_fork_change,persistence_focal_activity),post_fork_change)]
  effects <- s[,lapply(.SD,mean),by=ids,.SDcols=c(local,"focal_late")]
  counts <- s[,.(n_paired_samples=.N,n_paired_late_samples=sum(observed_late),
    n_full_panel_samples=sum(sibling_children==sibling_children_reference,na.rm=TRUE),
    paired_sibling_children=paste(sort(unique(as.integer(unlist(strsplit(
      sibling_children,",",fixed=TRUE))))),collapse=","),
    n_late_bins=sum(n_late_bins[observed_late]),
    n_late_phase_bins=sum(n_late_phase_bins[observed_late],na.rm=TRUE)),by=ids]
  effects <- merge(effects,counts,by=ids,sort=FALSE)
  observed <- s[observed_late %in% TRUE]
  if (nrow(observed)) {
    late_means <- observed[,lapply(.SD,function(a) .mean_or_na(a[is.finite(a)])),
      by=ids,.SDcols=c(late,"focal_late_change")]
    horizons <- observed[,.(persistence_phase_start=min(persistence_phase_start),
      persistence_phase_end=max(persistence_phase_end)),by=ids]
    effects <- merge(effects,merge(late_means,horizons,by=ids),by=ids,all.x=TRUE,sort=FALSE)
    effects[, 'late_phase_start' := persistence_phase_start]
    effects[, 'late_phase_end' := persistence_phase_end]
  }
  if (measurement_error) {
    errors <- s[,.(theta_measurement_se=sqrt(sum(theta_measurement_se^2))/.N),by=ids]
    effects <- merge(effects,errors,by=ids,sort=FALSE)
  }
  locations <- ref[effects,on=ids,which=TRUE]
  fields <- setdiff(names(effects),ids)
  out[locations,fields] <- as.data.frame(effects[,..fields])
  out$sibling_children[locations] <- out$paired_sibling_children[locations]
  out$theta <- out$sibling_contrast
  out$evaluable[locations] <- TRUE
  out$reason[locations] <- "evaluable"
  out$paired_sample_fraction[locations] <- out$n_paired_samples[locations]/length(sample_levels)
  out$paired_sample_status[locations] <- ifelse(out$n_paired_samples[locations] == 1L,
    "one_paired_sample","multiple_paired_samples")
  out$persistence_support[locations] <- ifelse(out$n_paired_late_samples[locations] == 0L,
    "local_focal_only", ifelse(out$n_paired_late_samples[locations] == out$n_paired_samples[locations],
      "paired_late_observed","paired_late_and_local"))
  out$phase_pattern[locations] <- ifelse(out$n_paired_late_samples[locations] == 0L,
    "late_unobserved",ifelse(out$sibling_contrast[locations]>0,
      ifelse(out$persistence_change_contrast[locations]>0,"early_and_late_positive","early_positive_only"),
      ifelse(out$persistence_change_contrast[locations]>0,"late_positive_only","no_positive_phase")))
  out
}


#' Rank transient molecular events within declared branch contexts.
#' @param fork_dynamics Output of SummarizeForkDynamics.
#' @param component Focal activity or branch-relative activity.
#' @param n Maximum number of molecular queries per context.
#' @param group_by Columns defining one query context.
#' @param aggregation Rank individual events, or mean bounded pulse evidence
#'   times amplitude across the planned children of each molecular query.
#' @return A data frame of pulse-model queries ranked by their BIC advantage.
#' @export
RankTransientForkEvents <- function(fork_dynamics,component=c("relative","focal"),n=1L,
                                    group_by=c("event_node","child","source","target"),
                                    aggregation=c("event","mean_children")) {
  component <- match.arg(component)
  n <- .validate_positive_int(n,"n",min=1L)
  prefix <- switch(component,relative="",focal="focal_")
  x <- data.table::copy(data.table::as.data.table(fork_dynamics))
  x[,temporal_score:=get(paste0(prefix,"pulse_bic_gain"))]
  aggregation <- match.arg(aggregation)
  x <- switch(aggregation,event={
    x[is.finite(temporal_score)&temporal_score>0&
      get(paste0(prefix,"dynamic_pattern"))%in%"pulse"]
  },mean_children={
    observed <- x[[paste0(prefix,"dynamic_status")]]%in%"observed"
    pulse <- observed & x[[paste0(prefix,"dynamic_pattern")]]%in%"pulse"
    evidence <- -expm1(-pmax(data.table::fcoalesce(x$temporal_score,0),0)/2)
    amplitude <- abs(data.table::fcoalesce(x[[paste0(prefix,"dynamic_effect")]],0))
    x[,`:=`(temporal_score=ifelse(pulse,evidence*amplitude,0),
      pulse_observed=pulse,profile_observed=observed)]
    keys <- unique(c(group_by,"pair_idx","source","target"))
    x[,.(temporal_score=mean(temporal_score),n_planned_children=.N,
      n_observed_children=sum(profile_observed),n_pulse_children=sum(pulse_observed)),
      by=keys][temporal_score>0]
  })
  ordering <- unique(c(group_by,"temporal_score","pair_idx","source","target"))
  data.table::setorderv(x,ordering,ifelse(ordering=="temporal_score",-1L,1L))
  as.data.frame(x[,head(.SD,n),by=group_by])
}


#' Sample-paired levels across observed receiver subtrees
#' @export
SummarizeObservedSiblingLevels <- function(temporal_lr, event_units,
    activity_transform = c("log_exposure","hill")) {
  activity_transform <- match.arg(activity_transform)
  raw <- data.table::copy(data.table::as.data.table(temporal_lr$scores))
  observations <- data.table::as.data.table(temporal_lr$observations)
  sides <- data.table::rbindlist(lapply(event_units$units,"[",
    c("event_node","children")),idcol="unit_id")
  data.table::setnames(sides,"children","child")
  contexts <- switch(as.character(is.null(temporal_lr$requested_contexts)),
    "TRUE"=unique(raw[,.(pair_idx,source,target)]),
    "FALSE"=unique(data.table::as.data.table(temporal_lr$requested_contexts)[,
      .(pair_idx,source,target)]))
  conditions <- unique(observations[,.(condition)])
  contexts[,join_key:=1L]
  conditions[,join_key:=1L]
  contexts <- merge(contexts,conditions,by="join_key",allow.cartesian=TRUE)
  sides[,join_key:=1L]
  frame <- merge(sides,contexts,by="join_key",allow.cartesian=TRUE)[,join_key:=NULL]
  planned <- observations[,.(n_planned_samples=data.table::uniqueN(biological_sample)),
    by=condition]
  frame <- merge(frame,planned,by="condition",all.x=TRUE)
  raw[,child:=temporal_lr$profile_nodes[lineage]]
  # Finite zero opportunity is retained; absent cell support remains missing.
  raw <- raw[structural_valid %in% TRUE & is.finite(score)]
  raw[,activity:=switch(activity_transform,log_exposure=-log1p(-score),hill=score)]
  raw <- raw[is.finite(activity)]
  measured <- merge(raw,sides[,.(unit_id,event_node,child)],by="child",
    allow.cartesian=TRUE)
  sample_keys <- c("condition","unit_id","event_node","pair_idx","source","target",
    "biological_sample","observed_time")
  totals <- measured[,.(activity_sum=sum(activity),n_observed_sides=.N,
    observed_children=paste(sort(child),collapse=",")),by=sample_keys]
  samples <- merge(measured,totals,by=sample_keys)[n_observed_sides>1L]
  samples[,sister_activity:=(activity_sum-activity)/(n_observed_sides-1L)]
  samples[,level_difference:=activity-sister_activity]
  keys <- c("condition","unit_id","event_node","child","pair_idx","source","target")
  effects <- samples[,.(level_effect=mean(level_difference),
    focal_activity=mean(activity),sister_activity=mean(sister_activity),
    n_samples=data.table::uniqueN(biological_sample),
    n_dates=data.table::uniqueN(observed_time),first_observed_day=min(observed_time),
    last_observed_day=max(observed_time)),by=keys]
  effects <- merge(frame,effects,by=keys,all.x=TRUE)
  effects[is.na(n_samples),n_samples:=0L]
  effects[is.na(n_dates),n_dates:=0L]
  effects[,evaluable:=n_samples>0L & is.finite(level_effect)]
  effects[,paired_sample_fraction:=n_samples/n_planned_samples]
  effects[,status:=data.table::fifelse(evaluable,"observed_paired_level",
    "no_paired_observed_context")]
  list(effects=effects,sample_effects=samples,
    activity_transform=activity_transform,sampling_unit=temporal_lr$sampling_unit,
    observation_scope="sample_paired_observed_sibling_levels")
}


# Collapse degree-two RNA states into the branch segments actually compared.
# Terminal column order and the fine-state forest remain unchanged.
.build_branch_segment_tree <- function(lineage_paths, lineage_clusters, pseudotime, knn_conn) {
  paths <- data.table::rbindlist(lapply(seq_along(lineage_paths), function(i) {
    p <- as.integer(lineage_paths[[i]])
    data.table::data.table(lineage = i, position = seq_along(p), node = p)
  }))
  paths[, parent := data.table::shift(node), by = lineage]
  edges <- unique(paths[!is.na(parent), .(parent, child = node)])
  forks <- edges[, .N, by = parent][N >= 2L, parent]
  terminals <- vapply(lineage_paths, tail, integer(1L), n = 1L)
  boundaries <- unique(c(forks, terminals))
  paths[, segment_position := rev(cummin(rev(ifelse(node %in% boundaries, position, Inf)))), by = lineage]
  paths[, segment := node[segment_position], by = lineage]
  mapping <- unique(paths[, .(node, segment)])
  stopifnot(!anyDuplicated(mapping$node))
  segment_paths <- lapply(lineage_paths, function(p) as.integer(p[p %in% boundaries]))
  segment_labels <- mapping$segment[match(lineage_clusters, mapping$node)]
  stopifnot(!anyNA(segment_labels), identical(vapply(segment_paths, tail, integer(1L), n = 1L), terminals))
  tree <- data.table::rbindlist(lapply(segment_paths, function(p) {
    data.table::data.table(event_node = head(p, -1L), child = tail(p, -1L))
  }))
  tree <- unique(tree)
  graph <- data.table::as.data.table(summary(as(knn_conn, "CsparseMatrix")))
  graph[, `:=`(event_node = segment_labels[i], child = segment_labels[j],
    time = (pseudotime[i] + pseudotime[j]) / 2)]
  crossings <- merge(graph[x > 0 & event_node != child], tree,
    by = c("event_node", "child"))
  anchors <- crossings[, .(boundary_time = stats::weighted.mean(time, x),
    boundary_edges = .N, boundary_mass = sum(x)), by = .(event_node, child)]
  anchors <- anchors[, .(fork_time = mean(boundary_time), n_observed_children = .N,
    boundary_edges = sum(boundary_edges), boundary_mass = sum(boundary_mass)), by = event_node]
  list(lineage_paths = segment_paths, lineage_clusters = as.integer(segment_labels),
    node_map = as.data.frame(mapping), anchors = as.data.frame(anchors),
    observation_scope = "observed_branch_segments")
}

# One raw population defines both segment contrasts and temporal profiles.
.event_segment_weights <- function(unit, lineage_weights, lineage_clusters, prefork_weighting) {
  weights <- .event_child_weights(unit, as.matrix(lineage_weights))
  mass <- rowSums(weights)
  weights <- weights / ifelse(mass == 0, 1, mass)
  child_columns <- seq_along(unit$children)
  if (prefork_weighting == "shared") {
    parent_cells <- lineage_clusters %in% unit$segment_pre_nodes
    parent_mass <- rowSums(weights[parent_cells, child_columns, drop = FALSE])
    weights[parent_cells, child_columns] <- parent_mass / length(child_columns)
  }
  allowed <- vapply(child_columns, function(k) {
    as.numeric(lineage_clusters %in% c(unit$segment_pre_nodes, unit$segment_post_nodes[[k]]))
  }, numeric(length(lineage_clusters)))
  child_weights <- weights[, child_columns, drop = FALSE] * allowed
  cbind(child_weights, pmax(0, 1 - rowSums(child_weights)))
}

.event_population_profiles <- function(expr, celltypes, lr_db, pseudotime,
    lineage_weights, lineage_clusters, event_units, bin_bounds, hill_scale,
    min_pct = 0.05, min_n_eff = 5, requested_contexts = NULL) {
  bounds <- as.numeric(bin_bounds[[1L]])
  effective_min_pct <- if (identical(event_units$activity_mode, "continuous") &&
    !identical(event_units$expression_mean_mode, "positive_cells")) 0 else min_pct
  scores <- data.table::rbindlist(lapply(event_units$units, function(unit) {
    weights <- .event_segment_weights(unit, lineage_weights, lineage_clusters,
      event_units$prefork_weighting)
    value <- data.table::as.data.table(ComputeTemporalLRFixedBounds(
      expr, celltypes, lr_db, pseudotime, weights,
      bin_bounds = replicate(ncol(weights), bounds, simplify = FALSE),
      hill_scale = hill_scale, min_pct = effective_min_pct, min_n_eff = min_n_eff,
      requested_contexts = requested_contexts, return_complete_grid = TRUE))
    value <- value[lineage <= length(unit$children)]
    value[, `:=`(event_node = unit$event_node,
      next_node = unit$children[lineage], mid = (bounds[bin] + bounds[bin + 1L]) / 2)]
    value
  }), use.names = TRUE, fill = TRUE)
  list(scores = scores, intermediate = list(bin_bounds = bin_bounds,
    hill_scale = hill_scale, n_bins = length(bounds) - 1L),
    observation_scope = "native_event_segment_populations",
    activity_mode = event_units$activity_mode, min_pct = effective_min_pct)
}


#' RNA availability lost on removal of one source population.
#'
#' All-source and source-removed ligand means use the same total-cell
#' denominator. Receptor expression and assayed composition are held fixed.
#' Complex assembly follows averaging, using the limiting subunit.
#' @param expr Genes by cells transformed RNA expression used for LR scoring.
#' @param celltypes Integer population code per cell.
#' @param lr_db Ligand and receptor subunit rows.
#' @param sample_labels Biological sample per cell.
#' @param source,target Source and receiving population codes.
#' @param hill_scale Positive common LR saturation scale.
#' @param min_n_eff Minimum observed source and receiving cells per sample.
#' @return Per-sample, per-pair source-removal and availability estimates.
#' @export
ComputeSourceRemovalLR <- function(expr, celltypes, lr_db, sample_labels,
    source, target, hill_scale, min_n_eff = 5) {
  samples <- unique(as.character(sample_labels))
  sample_index <- match(as.character(sample_labels), samples)
  membership <- Matrix::sparseMatrix(i = seq_len(ncol(expr)), j = sample_index,
    x = 1, dims = c(ncol(expr), length(samples)))
  source_membership <- Matrix::Diagonal(x = as.numeric(celltypes == source)) %*% membership
  target_membership <- Matrix::Diagonal(x = as.numeric(celltypes == target)) %*% membership
  n_total <- as.numeric(Matrix::colSums(membership))
  n_source <- as.numeric(Matrix::colSums(source_membership))
  n_target <- as.numeric(Matrix::colSums(target_membership))
  total_mean <- sweep(as.matrix(expr %*% membership), 2L, n_total, "/")
  source_contribution <- sweep(as.matrix(expr %*% source_membership), 2L, n_total, "/")
  remaining_mean <- pmax(total_mean - source_contribution, 0)
  source_mean <- sweep(as.matrix(expr %*% source_membership), 2L, pmax(n_source, 1), "/")
  target_mean <- sweep(as.matrix(expr %*% target_membership), 2L, pmax(n_target, 1), "/")
  ligands <- data.table::data.table(
    pair_idx = rep(seq_along(lr_db$ligand_rows), lengths(lr_db$ligand_rows)),
    row = unlist(lr_db$ligand_rows, use.names = FALSE))
  receptors <- data.table::data.table(
    pair_idx = rep(seq_along(lr_db$receptor_rows), lengths(lr_db$receptor_rows)),
    row = unlist(lr_db$receptor_rows, use.names = FALSE))
  ligand_values <- data.table::data.table(
    pair_idx = rep(ligands$pair_idx, times = length(samples)),
    sample_index = rep(seq_along(samples), each = nrow(ligands)),
    total = as.vector(total_mean[ligands$row, , drop = FALSE]),
    remaining = as.vector(remaining_mean[ligands$row, , drop = FALSE]),
    source_mean = as.vector(source_mean[ligands$row, , drop = FALSE]))[,
      .(total_ligand = min(total), remaining_ligand = min(remaining),
        source_ligand = min(source_mean)), by = .(pair_idx, sample_index)]
  receptor_values <- data.table::data.table(
    pair_idx = rep(receptors$pair_idx, times = length(samples)),
    sample_index = rep(seq_along(samples), each = nrow(receptors)),
    value = as.vector(target_mean[receptors$row, , drop = FALSE]))[,
      .(receptor = min(value)), by = .(pair_idx, sample_index)]
  values <- merge(ligand_values, receptor_values, by = c("pair_idx", "sample_index"))
  values[, `:=`(biological_sample = samples[sample_index], source = source, target = target,
    n_total = n_total[sample_index], n_source = n_source[sample_index],
    n_target = n_target[sample_index], source_fraction = n_source[sample_index] / n_total[sample_index])]
  values[, `:=`(total_activity = log1p(sqrt(total_ligand * receptor) / hill_scale),
    remaining_activity = log1p(sqrt(remaining_ligand * receptor) / hill_scale),
    source_activity = log1p(sqrt(source_ligand * receptor) / hill_scale))]
  values[, source_loss := pmax(0, total_activity - remaining_activity)]
  values[, source_loss_fraction := source_loss / pmax(total_activity, .Machine$double.eps)]
  values[, evaluable := n_source >= min_n_eff & n_target >= min_n_eff]
  values[evaluable %in% FALSE, c("total_activity", "remaining_activity", "source_activity",
    "source_loss", "source_loss_fraction") := NA_real_]
  as.data.frame(values)
}
