# Shared state-clock stage for full and dependency-aware runs.
  state_clock <- NULL
  switch(as.character(main_candidate), "TRUE" = {
    clock_anchor_order <- as.numeric(pseudotime)
    clock_anchor_start <- as.integer(start_cell)
    clock_input_clusters <- as.integer(lineage_cl)
    clock_path_nodes <- data.table(
      lineage = rep(seq_along(lineage_paths), lengths(lineage_paths)),
      state = unlist(lineage_paths, use.names = FALSE))
    clock_path_nodes[, parent := shift(state), by = lineage]
    clock_tree_edges <- unique(clock_path_nodes[!is.na(parent), .(from = parent, to = state)])
    stopifnot(!anyDuplicated(clock_tree_edges$to))
    clock_parent <- integer(max(clock_input_clusters))
    clock_parent[clock_tree_edges$to] <- clock_tree_edges$from
    state_clock_conn <- as(as(state_clock_conn, "generalMatrix"), "CsparseMatrix")
    clock_geometry_anchor <- clock_anchor_order
    clock_anchor_signal <- list(weight = 0, scope = "geometry_anchor")
    switch(as.character(startsWith(selected_root_potency_mode, "detected_genes_high") &&
      identical(root_evidence_mode, "anchored_bidirectional_potency_ordering") &&
      CONFIG$llr$component_ablation != "no_coherent_clock_signal"), "TRUE" = {
      # Estimate the depth relation within assay x cell state. Standardizing
      # to a common total retains the fitted state/assay intercepts.
      clock_strata_columns <- intersect(CONFIG$llr$clock_depth_strata_columns, names(meta))
      clock_assay <- switch(as.character(length(clock_strata_columns) > 0L),
        "TRUE" = as.character(meta[[clock_strata_columns[1L]]]),
        "FALSE" = rep.int("single_assay", n_cells))
      clock_fit <- data.table(assay = clock_assay, state = celltypes,
        log_detection = log1p(as.numeric(Matrix::colSums(expr > 0))),
        log_total = log1p(as.numeric(Matrix::colSums(expr))))
      clock_fit[, c("centered_detection", "centered_total") := .(
        log_detection - mean(log_detection), log_total - mean(log_total)),
        by = .(assay, state)]
      clock_depth_variance <- sum(clock_fit$centered_total^2)
      clock_depth_identified <- clock_depth_variance >
        .Machine$double.eps * sum(clock_fit$log_total^2)
      clock_unconstrained_depth_slope <- ifelse(clock_depth_identified,
        sum(clock_fit$centered_total * clock_fit$centered_detection) /
          clock_depth_variance, 0)
      clock_depth_slope <- pmax(0, clock_unconstrained_depth_slope)
      clock_standardized_detection <- clock_fit$log_detection -
        clock_depth_slope * (clock_fit$log_total - median(clock_fit$log_total))
      clock_fit[, standardized_detection := clock_standardized_detection]
      clock_assay_rank <- grepl("within_sample_high_root", root_input_source, fixed = TRUE)
      switch(as.character(clock_assay_rank), "TRUE" = {
        clock_fit[, standardized_detection :=
          data.table::frank(standardized_detection, ties.method = "average") / .N,
          by = assay]
      })
      clock_maturity <- (data.table::frank(
        -clock_fit$standardized_detection, ties.method = "average") - 1) / (n_cells - 1)
      clock_neighbor_maturity <- as.numeric(state_clock_conn %*% clock_maturity) /
        as.numeric(Matrix::rowSums(state_clock_conn))
      clock_signal <- data.table(state = clock_input_clusters,
        maturity = clock_maturity, neighbor = clock_neighbor_maturity)
      clock_signal[, c("maturity_local", "neighbor_local") := .(
        maturity - mean(maturity), neighbor - mean(neighbor)), by = state]
      clock_within_fraction <- sum(clock_signal$maturity_local^2) /
        sum((clock_maturity - mean(clock_maturity))^2)
      clock_local_coherence <- cor(clock_signal$maturity_local,
        clock_signal$neighbor_local, method = "spearman")
      clock_coherent_variation <- clock_within_fraction * pmax(clock_local_coherence, 0)^2
      clock_anchor_weight <- clock_coherent_variation /
        (clock_coherent_variation + 1 - clock_within_fraction)
      clock_anchor_order <- (1 - clock_anchor_weight) * clock_geometry_anchor +
        clock_anchor_weight * (clock_maturity + clock_neighbor_maturity) / 2
      clock_anchor_signal <- list(weight = clock_anchor_weight,
        within_fraction = clock_within_fraction, local_coherence = clock_local_coherence,
        depth_slope = clock_depth_slope,
        unconstrained_depth_slope = clock_unconstrained_depth_slope,
        assay_rank_from_root_policy = clock_assay_rank,
        depth_variance = clock_depth_variance,
        depth_identified = clock_depth_identified,
        depth_stratum_column = paste(clock_strata_columns[1L], collapse = ""),
        scope = "depth_standardized_coherent_gene_detection")
    })

    native_clock <- .Call(LineaLinkR:::wrap__rs_compute_state_clock,
      as.numeric(geom_emb), state_clock_conn, clock_input_clusters,
      clock_anchor_order, clock_parent, as.integer(ncol(geom_emb)), as.integer(CONFIG$n_threads))
    pseudotime <- as.numeric(native_clock$pseudotime)
    clock_columns <- c("n_owned", "n_members", "local_rank", "variance_fraction",
      "eigengap_fraction", "anchor_mean", "mean_time", "local_residual")
    clock_charts <- as.data.table(matrix(native_clock$chart_info,
      ncol = length(clock_columns), byrow = TRUE))
    setnames(clock_charts, clock_columns)
    clock_state <- clock_input_clusters[match(seq_len(nrow(clock_charts)), native_clock$calibration_chart)]
    clock_charts[, `:=`(chart = seq_len(.N), component = native_clock$chart_components,
                       state = clock_state)]
    clock_diagnostic_names <- c("anchor_rss", "local_rss", "normal_residual",
      "order_violation", "complementarity", "minimum_eigenvalue", "maximum_eigenvalue",
      "raw_minimum", "raw_maximum", "active_order_constraints")
    clock_diagnostics <- setNames(as.numeric(native_clock$diagnostics), clock_diagnostic_names)
    clock_state_time <- data.table(state = lineage_cl, time = pseudotime,
      anchor = clock_anchor_order)[,
      .(mean_time = mean(time), anchor_mean = mean(anchor)), by = state]
    clock_tree_edges[, mean_increment := clock_state_time$mean_time[match(to, clock_state_time$state)] -
      clock_state_time$mean_time[match(from, clock_state_time$state)]]
    clock_tree_edges[, anchor_increment := clock_state_time$anchor_mean[match(to, clock_state_time$state)] -
      clock_state_time$anchor_mean[match(from, clock_state_time$state)]]
    clock_tree_edges[, `:=`(
      raw_mean_increment = mean_increment * (clock_diagnostics[["raw_maximum"]] -
        clock_diagnostics[["raw_minimum"]]),
      allowed_reversal = pmax(-anchor_increment, 0))]
    state_clock <- list(charts = as.data.frame(clock_charts),
      diagnostics = clock_diagnostics, anchor_signal = clock_anchor_signal,
      geometry_anchor_order = clock_geometry_anchor,
      reference_edges = as.data.frame(clock_tree_edges),
      reference_parent = clock_parent, calibration_chart = native_clock$calibration_chart,
      input_clusters = clock_input_clusters, anchor_order = clock_anchor_order,
      anchor_start_cell = clock_anchor_start,
      cell_component = native_clock$chart_components[native_clock$calibration_chart],
      minimum_cells = which(pseudotime == 0), anchor_diagnostics = direction_diagnostics,
      coordinate_scope = "global_geometry_anchored_ordinal_progress",
      order_constraint = "Bt<=max(Bg,0)", reference_constraint_scale = "before_global_minmax")
    lineage_cluster_source <- paste0(lineage_cluster_source, "+feasible_anchored_expression_clock")
    direction_status <- "conditional_global_geometry_and_reference_forest"
    direction_diagnostics <- list(direction_status = direction_status,
      clock_shape_uses_anchor = TRUE, reference_roots = which(clock_parent == 0L),
      observed_component_count = uniqueN(native_clock$chart_components))
    root_diagnostics$direction_diagnostics <- direction_diagnostics
    cat(sprintf("  [state clock: %d local charts, %d observed components, %d reference paths; normal residual %.3g]\n",
      nrow(clock_charts), uniqueN(native_clock$chart_components), n_lineages,
      clock_diagnostics[["normal_residual"]]))
  })

  branch_clock_calibration <- NULL
  if (!nzchar(Sys.getenv("LLR_NO_BRANCH_CALIBRATE"))) {
    pt_calibrated <- LlrCalibrateBranchClock(as.numeric(pseudotime),
      as.integer(lineage_cl), lineage_paths)
    branch_clock_calibration <- list(
      n_subtrees_calibrated = attr(pt_calibrated, "n_subtrees_calibrated"),
      pseudotime_raw = as.numeric(pseudotime))
    pseudotime <- as.numeric(pt_calibrated)
    # calibration renormalized to [0,1] (rank-preserving; downstream binning contract)
    pseudotime <- (pseudotime - min(pseudotime)) / diff(range(pseudotime))
    cat(sprintf("  [branch clock calibration: %d subtrees rescaled]\n",
      branch_clock_calibration$n_subtrees_calibrated))
  }

