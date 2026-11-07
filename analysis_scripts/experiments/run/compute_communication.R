# Shared communication stage; numerical methods come from the installed package.
  lineage_curve_weights <- switch(as.character(use_assignment_layout),
    "TRUE" = lineage_assignment_weights,
    switch(component_ablation, local_curve_reference = lineage_assignment_weights,
      local_point_reference = lineage_reference_weights, lineage_weights))
  lineage_shape_weights <- switch(component_ablation,
    assigned_shape_reference = lineage_assignment_weights, lineage_curve_weights)
  reference_temporal_lr <- NULL
  lineage_persistence_weights <- switch(as.character(use_assignment_layout),
    "TRUE" = lineage_weights, lineage_curve_weights)
  shape_temporal_lr <- NULL
  sample_temporal_lr <- NULL
  temporal_lr <- NULL
  if(!no_temporal && !is.null(lr_db)) {
    cat("  Step 6/7: Temporal LR (binned) ... ")
    temporal_pseudotime <- if (!is.null(pseudotime_per_lineage)) pseudotime_per_lineage else pseudotime
    temporal_lr <- LineaLinkR::ComputeTemporalLRWithIntermediate(
      expr_lr, celltypes, lr_db, temporal_pseudotime, lineage_curve_weights,
      n_bins = llr_cfg$n_bins %||% 12L,
      min_n_eff = llr_cfg$min_n_eff %||% 5.0,
      common_grid = TRUE
    )
    shape_temporal_lr <- temporal_lr

    cat(sprintf("done [%d rows]\n", nrow(temporal_lr$scores)))
  } else if (is.null(lr_db)) { cat("  Step 6/7: SKIP (no LR DB for this species)\n")
  } else { cat("  Step 6/7: ABLATION (no temporal LR)\n") }

  # ── Step 7: Entropy + bifurcation + role scores ──────────────────
  entropy_raw <- h_tilde <- bif_df <- NULL
  fork_contrasts <- data.frame(status = "not_evaluable")
  fork_evidence <- data.frame(status = "not_evaluable")
  fork_bins <- NULL
  fork_sample_contrasts <- NULL
  fork_sample_summary <- NULL
  fork_sample_column <- FindBiologicalSampleCol(meta, ds_name)
  paired_amplitude <- main_candidate && !is.na(fork_sample_column) &&
    component_ablation != "no_donor_local_conditioning"
  fork_inference_unit <- "pooled_cells_descriptive"
  if(!no_entropy && !no_temporal && !is.null(temporal_lr)) {
    cat("  Step 7/7: Lineage entropy & bifurcation ... ")
    entropy_raw <- LineaLinkR::ComputeLineageEntropyRaw(temporal_lr$intermediate)
    h_tilde <- LineaLinkR::SmoothMA3(entropy_raw)
    bif_df <- LineaLinkR::FindBifurcationBins(h_tilde)
    cat(sprintf("b*: %s\n",
      paste(sprintf("L%d=%s", bif_df$lineage,
        ifelse(is.na(bif_df$b_star), "NA", bif_df$b_star)),
        collapse = ", ")))
  } else if (is.null(temporal_lr)) { cat("  Step 7/7: SKIP (no temporal LR)\n")
  } else { cat("  Step 7/7: ABLATION (no entropy)\n") }

  if (!is.null(temporal_lr) && !is.null(bif_df) &&
      nrow(temporal_lr$scores) > 0L && nrow(sibling_map) > 0L) {
    fork_bins <- LineaLinkR::LocateForkBinsOnGrid(
      sibling_map, lineage_cl, temporal_pseudotime,
      attr(temporal_lr$scores, "bin_bounds")
    )
    if (enable_event_mode) {
      event_expression_mean_mode <- "all_cells"
      switch(as.character(enable_dropout_aware_event_mean), "TRUE" = {
        lr_rows <- unique(unlist(c(lr_db$ligand_rows,
                                   lr_db$receptor_rows), use.names = FALSE))
        lr_expression <- expr[lr_rows, , drop = FALSE]
        lr_detection <- Matrix::rowMeans(lr_expression > 0)
        lr_positive_mean <- Matrix::rowMeans(lr_expression) / lr_detection
        lr_detection_intensity_rho <- stats::cor(
          lr_detection, lr_positive_mean,
          method = "spearman", use = "complete.obs"
        )
        event_expression_mean_mode <- ifelse(
          is.finite(lr_detection_intensity_rho) &&
            lr_detection_intensity_rho < 0,
          "positive_cells", "all_cells"
        )
      })
      # Event-based fork contrast: enumerate all tree-level fork events
      # (L1) or nearest-sibling units (L0), with exact (D1) or fixed-bin
      # (D0) boundaries. Scores computed directly from raw expression.
      fork_contrasts <- event_contrast_method(
        expr = expr_lr,
        measurement_error = llr_cfg$context_error_weighting && !paired_amplitude,
        prefork_weighting = llr_cfg$prefork_weighting,
        activity_transform = llr_cfg$event_activity_transform,
        celltypes = celltypes,
        lr_db = lr_db,
        pseudotime = temporal_pseudotime,
        lineage_weights = lineage_weights,
        lineage_paths = lineage_paths,
        sibling_map = sibling_map,
        lineage_clusters = lineage_cl,
        fork_bins = fork_bins,
        bin_bounds = attr(temporal_lr$scores, "bin_bounds"),
        event_scope = event_scope_val,
        boundary = event_boundary_val,
        window_bins = llr_cfg$fork_window_bins %||% 2L,
        hill_scale = attr(temporal_lr$scores, "hill_scale"),
        min_n_eff = llr_cfg$min_n_eff %||% 5.0,
        min_pct = llr_cfg$min_pct %||% 0.05,
        activity_floor = llr_cfg$fork_activity_floor %||% 1e-4,
        bif_df = bif_df,
        descendant_window_children = transition_fragment_labels,
        anchor_rule = event_anchor_val,
        scoring_mode = event_scoring_val,
        expression_mean_mode = event_expression_mean_mode,
        activity_mode = event_activity_val,
        sibling_selection = switch(as.character(program_lineage), "TRUE" = "program_matched",
          ifelse(run_tag == "L1D0contMaturityMatchedInterval", "time_matched", "all")),
        sibling_program_profiles = program_profiles
      )
      if (nrow(fork_contrasts) > 0L) {
        persistence_temporal_lr <- temporal_lr
        switch(as.character(use_assignment_layout), "TRUE" = {
          # Assignment determines the observation layout. Geometric weights
          # define both the local baseline and the late population amplitude.
          reference_scores <- LineaLinkR::ComputeTemporalLRFixedBounds(
            expr_lr, celltypes, lr_db, temporal_pseudotime, lineage_weights,
            bin_bounds = temporal_lr$intermediate$bin_bounds,
            hill_scale = temporal_lr$intermediate$hill_scale,
            min_pct = llr_cfg$min_pct %||% 0.05, min_n_eff = llr_cfg$min_n_eff %||% 5.0,
            requested_contexts = unique(as.data.table(fork_contrasts)[, .(pair_idx, source, target)]),
            return_complete_grid = TRUE)
          persistence_temporal_lr <- list(scores = reference_scores,
            intermediate = temporal_lr$intermediate[c("bin_bounds","hill_scale","n_bins")])
        })
        switch(component_ablation, assigned_shape_reference = {
          # Hold the geometric population, grid and event definitions fixed.
          # Assignment weights determine conditional temporal shape only.
          shape_scores <- LineaLinkR::ComputeTemporalLRFixedBounds(
            expr_lr, celltypes, lr_db, temporal_pseudotime, lineage_shape_weights,
            bin_bounds = temporal_lr$intermediate$bin_bounds,
            hill_scale = temporal_lr$intermediate$hill_scale,
            min_pct = llr_cfg$min_pct %||% 0.05, min_n_eff = llr_cfg$min_n_eff %||% 5.0,
            requested_contexts = unique(as.data.table(fork_contrasts)[, .(pair_idx, source, target)]),
            return_complete_grid = TRUE)
          shape_temporal_lr <- list(scores = shape_scores,
            intermediate = temporal_lr$intermediate[c("bin_bounds","hill_scale","n_bins")])
        })
        switch(as.character(enable_late_bin_merging), "TRUE" = {
          # Merge adjacent fine bins before scoring; never average missing LR
          # scores or lower the effective-cell threshold.
          bounds <- temporal_lr$intermediate$bin_bounds
          edge_indices <- sort(unique(c(seq.int(1L, length(bounds[[1L]]),
            by = llr_cfg$late_bin_merge_factor), length(bounds[[1L]]))))
          coarse_matrix <- do.call(rbind, bounds)[, edge_indices, drop = FALSE]
          coarse_bounds <- split(coarse_matrix, row(coarse_matrix))
          coarse_scores <- LineaLinkR::ComputeTemporalLRFixedBounds(
            expr_lr, celltypes, lr_db, pseudotime, lineage_persistence_weights,
            bin_bounds = coarse_bounds, hill_scale = temporal_lr$intermediate$hill_scale,
            min_pct = llr_cfg$min_pct %||% 0.05, min_n_eff = llr_cfg$min_n_eff %||% 5.0,
            requested_contexts = unique(as.data.table(fork_contrasts)[, .(pair_idx, source, target)]),
            return_complete_grid = TRUE)
          temporal_lr$late_fallback <- list(scores = coarse_scores,
            intermediate = list(bin_bounds = coarse_bounds,
              n_bins = ncol(coarse_matrix) - 1L,
              hill_scale = temporal_lr$intermediate$hill_scale))
        })

        persistence_temporal_lr$late_fallback <- temporal_lr$late_fallback
        switch(as.character(use_assignment_layout), "TRUE" = {
          reference_temporal_lr <- persistence_temporal_lr
          temporal_lr$late_fallback <- NULL
        })
        fork_contrasts <- llr_event_methods$AddLateForkPersistence(
          fork_contrasts, persistence_temporal_lr, lineage_paths,
          expression_mean_mode = event_expression_mean_mode,
          persistence_mode = ifelse(signed_persistence, "signed", "positive"),
          sibling_scope = llr_cfg$persistence_sibling_scope
        )
        # The main path fits sample-conditioned dynamics below, before any read.
        switch(as.character(!main_candidate), "TRUE" = {
          fork_contrasts <- LineaLinkR::SummarizeForkDynamics(
            fork_contrasts, shape_temporal_lr, lineage_paths,
            expression_mean_mode = event_expression_mean_mode)
        })
        switch(as.character(run_tag %in% c("L1D0contBranchBalanceInterval",
                                              "L1D0contMaturityMatchedInterval",
                                              "L1D0contDonorSplitGeometryInterval",
                                              "L1D0contDonorBranchBalanceInterval",
                                              "L1D0contDonorBalanceNoSpurInterval")),
        "TRUE" = {
          balance_paths <- as.data.table(stack(setNames(
            lineage_paths, seq_along(lineage_paths))))
          setnames(balance_paths, c("values", "ind"), c("node", "lineage"))
          balance_paths[, position := seq_len(.N), by = lineage]
          balance_paths[, event_node := shift(node), by = lineage]
          balance_edges <- unique(balance_paths[!is.na(event_node),
            .(event_node, child = node, lineage, edge_position = position)])
          balance_descendants <- unique(merge(
            balance_edges, balance_paths[, .(lineage, position, node)],
            by = "lineage", allow.cartesian = TRUE)[
              position >= edge_position, .(event_node, child, node)])
          balance_cluster_size <- data.table(node = lineage_cl)[
            , .(n_cells = .N), by = node]
          branch_balance <- merge(balance_descendants,
            balance_cluster_size, by = "node")[
              , .(branch_cells = sum(n_cells)), by = .(event_node, child)]
          branch_balance[, event_cells := sum(branch_cells), by = event_node]
          branch_balance[, balance_fraction :=
            pmin(branch_cells, event_cells - branch_cells) / event_cells]
          event_key <- paste(fork_contrasts$event_node,
                             fork_contrasts$child, sep = ":")
          balance_key <- paste(branch_balance$event_node,
                               branch_balance$child, sep = ":")
          balance_fraction <- branch_balance$balance_fraction[
            match(event_key, balance_key)]
          fork_contrasts$branch_balance_fraction <- balance_fraction
          fork_contrasts$branch_balance_score <-
            fork_contrasts$persistent_positive_score * sqrt(balance_fraction)
          switch(as.character(run_tag == "L1D0contMaturityMatchedInterval"),
            "TRUE" = {
              units <- attr(fork_contrasts, "event_units")$units
              unit_children <- lapply(units, `[[`, "children")
              unit_post_times <- lapply(units, `[[`, "child_post_time")
              child_times <- data.table(
                unit_id = rep(seq_along(units), lengths(unit_children)),
                child = as.integer(unlist(unit_children)),
                post_time = as.numeric(unlist(unit_post_times))
              )
              time_key <- paste(child_times$unit_id, child_times$child, sep = ":")
              focal_time <- child_times$post_time[match(paste(
                fork_contrasts$unit_id, fork_contrasts$child, sep = ":"
              ), time_key)]
              sibling_time <- child_times$post_time[match(paste(
                fork_contrasts$unit_id, fork_contrasts$sibling_children, sep = ":"
              ), time_key)]
              fork_contrasts$terminal_maturity_factor <-
                pmin(focal_time, sibling_time)
              fork_contrasts$branch_balance_score <-
                fork_contrasts$branch_balance_score *
                  fork_contrasts$terminal_maturity_factor
            },
            "FALSE" = NULL)
        })
      }
    } else if (enable_merged_windows) {
      # Variant B: thin sides are re-estimated on merged adjacent bins with
      # boundaries from cell time and support only (see function docs).
      fork_contrasts <- LineaLinkR::ComputeMergedIntervalForkContrasts(
        temporal_lr$scores,
        bif_df,
        sibling_map,
        expr_lr,
        celltypes,
        lr_db,
        temporal_pseudotime,
        lineage_weights,
        fork_bins,
        window_bins = llr_cfg$fork_window_bins %||% 2L,
        min_n_eff = llr_cfg$min_n_eff %||% 5.0,
        activity_floor = llr_cfg$fork_activity_floor %||% 1e-4,
        hill_scale = attr(temporal_lr$scores, "hill_scale")
      )
    } else {
      fork_contrasts <- LineaLinkR::ComputeForkContrasts(
        temporal_lr$scores,
        bif_df,
        sibling_map,
        window_bins = llr_cfg$fork_window_bins %||% 2L,
        activity_floor = llr_cfg$fork_activity_floor %||% 1e-4,
        min_n_eff = llr_cfg$min_n_eff %||% 5.0,
        fork_bins = fork_bins,
        common_grid = TRUE
      )
    }
    if (!is.na(fork_sample_column) &&
        variant_id != "L1D0contDonorStructuralInterval") {
      sample_labels <- as.character(meta[[fork_sample_column]])
      valid_sample_labels <- unique(sample_labels[
        !is.na(sample_labels) & nzchar(trimws(sample_labels))])
      pooled_bounds <- attr(temporal_lr$scores, "bin_bounds")
      pooled_hill_scale <- attr(temporal_lr$scores, "hill_scale")
      # Sample-level support must be scored on the same objects as the pooled
      # estimate: same pseudotime coordinate, same Hill scale, and - for event
      # mode - the pooled event/child set and windows (units_override). For
      # merged-window (B) mode the same estimator is applied per sample.
      sample_contrast_tables <- lapply(valid_sample_labels, function(block) {
        keep <- which(sample_labels == block)
        sample_pt <- if (is.matrix(temporal_pseudotime)) temporal_pseudotime[keep, , drop = FALSE] else temporal_pseudotime[keep]
        if (enable_event_mode && nrow(fork_contrasts) == 0L) {
          contrast_table <- fork_contrasts
        } else if (enable_event_mode) {
          contrast_table <- event_contrast_method(
            expr = expr_lr[, keep, drop = FALSE],
            prefork_weighting = llr_cfg$prefork_weighting,
            activity_transform = llr_cfg$event_activity_transform,
            measurement_error = paired_amplitude && llr_cfg$context_error_weighting,
            celltypes = celltypes[keep],
            lr_db = lr_db,
            pseudotime = sample_pt,
            lineage_weights = lineage_weights[keep, , drop = FALSE],
            lineage_paths = lineage_paths,
            lineage_clusters = lineage_cl[keep],
            units_override = attr(fork_contrasts, "event_units"),
            hill_scale = pooled_hill_scale,
            min_n_eff = llr_cfg$min_n_eff %||% 5.0,
            min_pct = llr_cfg$min_pct %||% 0.05
          )
        } else {
          score_table <- LineaLinkR::ComputeTemporalLRFixedBounds(
            expr_lr[, keep, drop = FALSE],
            celltypes[keep],
            lr_db,
            sample_pt,
            lineage_weights[keep, , drop = FALSE],
            bin_bounds = pooled_bounds,
            hill_scale = pooled_hill_scale,
            n_bins = llr_cfg$n_bins %||% 12L,
            min_n_eff = llr_cfg$min_n_eff %||% 5.0,
            return_complete_grid = TRUE
          )
          score_table$block_label <- block
          contrast_table <- if (enable_merged_windows) {
            LineaLinkR::ComputeMergedIntervalForkContrasts(
              score_table,
              bif_df,
              sibling_map,
              expr_lr[, keep, drop = FALSE],
              celltypes[keep],
              lr_db,
              sample_pt,
              lineage_weights[keep, , drop = FALSE],
              fork_bins,
              window_bins = llr_cfg$fork_window_bins %||% 2L,
              min_n_eff = llr_cfg$min_n_eff %||% 5.0,
              activity_floor = llr_cfg$fork_activity_floor %||% 1e-4,
              hill_scale = pooled_hill_scale,
              sample_col = "block_label"
            )
          } else {
            LineaLinkR::ComputeForkContrasts(
              score_table,
              bif_df,
              sibling_map,
              sample_col = "block_label",
              window_bins = llr_cfg$fork_window_bins %||% 2L,
              activity_floor = llr_cfg$fork_activity_floor %||% 1e-4,
              min_n_eff = llr_cfg$min_n_eff %||% 5.0,
              fork_bins = fork_bins,
              common_grid = TRUE
            )
          }
          rm(score_table)
          invisible(gc(verbose = FALSE))
        }
        contrast_table$block_label <- rep(block,nrow(contrast_table))
        switch(as.character(main_candidate && nrow(contrast_table) > 0L), "TRUE" = {
          requested_curves <- unique(as.data.table(contrast_table)[evaluable %in% TRUE,
            .(pair_idx,source,target)])
          sample_fine <- list()
          sample_shape <- NULL
          if (nrow(requested_curves)) {
            sample_fine <- list(scores=LineaLinkR::ComputeTemporalLRFixedBounds(
              expr_lr[,keep,drop=FALSE],celltypes[keep],lr_db,sample_pt,
              lineage_persistence_weights[keep,,drop=FALSE],bin_bounds=pooled_bounds,
              hill_scale=pooled_hill_scale,min_pct=llr_cfg$min_pct %||% 0.05,
              min_n_eff=llr_cfg$min_n_eff %||% 5.0,return_complete_grid=TRUE,
              requested_contexts=requested_curves),
              intermediate=list(bin_bounds=pooled_bounds,hill_scale=pooled_hill_scale))
            if (enable_late_bin_merging) {
              sample_fine$late_fallback <- list(scores=LineaLinkR::ComputeTemporalLRFixedBounds(
                expr_lr[,keep,drop=FALSE],celltypes[keep],lr_db,sample_pt,
                lineage_persistence_weights[keep,,drop=FALSE],bin_bounds=coarse_bounds,
                hill_scale=pooled_hill_scale,min_pct=llr_cfg$min_pct %||% 0.05,
                min_n_eff=llr_cfg$min_n_eff %||% 5.0,return_complete_grid=TRUE,
                requested_contexts=requested_curves),intermediate=list(bin_bounds=coarse_bounds))
            }
            sample_shape <- if (identical(lineage_shape_weights,lineage_persistence_weights)) {
              sample_fine$scores
            } else LineaLinkR::ComputeTemporalLRFixedBounds(
              expr_lr[,keep,drop=FALSE],celltypes[keep],lr_db,sample_pt,
              lineage_shape_weights[keep,,drop=FALSE],bin_bounds=pooled_bounds,
              hill_scale=pooled_hill_scale,min_pct=llr_cfg$min_pct %||% 0.05,
              min_n_eff=llr_cfg$min_n_eff %||% 5.0,return_complete_grid=TRUE,
              requested_contexts=requested_curves)
            sample_shape <- as.data.table(sample_shape)[structural_valid %in% TRUE &
              is.finite(weight_mass_cells) & weight_mass_cells>0]
            sample_shape[,biological_sample:=block]
          }
          contrast_table <- llr_event_methods$AddLateForkPersistence(contrast_table,
            sample_fine,lineage_paths,expression_mean_mode=event_expression_mean_mode,
            sibling_scope=llr_cfg$persistence_sibling_scope,
            persistence_mode=ifelse(signed_persistence,"signed","positive"))
          attr(contrast_table,"temporal_curve") <- sample_shape
        })
        contrast_table
      })
      fork_sample_contrasts <- rbindlist(sample_contrast_tables, fill = TRUE)
      switch(as.character(main_candidate && nrow(fork_contrasts) > 0L), "TRUE" = {
        sample_temporal_lr <- shape_temporal_lr
        sample_temporal_lr$scores <- rbindlist(
          lapply(sample_contrast_tables, attr, which = "temporal_curve"), use.names = TRUE)
      })
      if (any(fork_sample_contrasts$evaluable %in% TRUE &
              is.finite(fork_sample_contrasts$sibling_contrast))) {
        fork_sample_summary <- LineaLinkR::SummarizeForkContrastsBySample(
          fork_sample_contrasts,
          n_boot = llr_cfg$fork_sample_boot %||% 999L,
          seed = CONFIG$seed
        )
        fork_inference_unit <- "biological_sample_bootstrap"
      } else {
        fork_inference_unit <- "biological_samples_no_evaluable_contrast"
      }
    }
    if (main_candidate && nrow(fork_contrasts)>0L) {
      if (paired_amplitude) {
        fork_contrasts <- llr_event_methods$.paired_event_amplitudes(
          fork_contrasts,fork_sample_contrasts,valid_sample_labels,
          signed=signed_persistence,measurement_error=llr_cfg$context_error_weighting)
        if (nrow(fork_sample_contrasts)) {
          panels <- unique(as.data.table(fork_sample_contrasts)[
            evaluable %in% TRUE & is.finite(sibling_contrast),
            .(pair_idx,source,target,event_node,child,biological_sample=block_label,sibling_children)])
          sample_temporal_lr$paired_events <- unique(panels[,.(next_node=c(child[1L],
            as.integer(strsplit(sibling_children[1L],",",fixed=TRUE)[[1L]]))),
            by=.(pair_idx,source,target,event_node,child,biological_sample)])
        }
        fork_inference_unit <- "equal_weight_paired_samples"
      } else {
        fork_contrasts$amplitude_unit <- "pooled_cells"
        fork_contrasts$n_observed_samples <- fork_contrasts$n_paired_samples <-
          fork_contrasts$n_paired_late_samples <- NA_integer_
        fork_contrasts$paired_sample_fraction <- NA_real_
        fork_contrasts$paired_sample_status <- ifelse(is.na(fork_sample_column),
          "sample_labels_unavailable","pooled_amplitude_control")
      }
    }
    switch(as.character(main_candidate && nrow(fork_contrasts) > 0L), "TRUE" = {
      switch(as.character(is.null(sample_temporal_lr)), "TRUE" = {
        sample_temporal_lr <- shape_temporal_lr
        sample_temporal_lr$scores <- as.data.table(sample_temporal_lr$scores)
        sample_temporal_lr$scores[, biological_sample := "pooled"]
      })
      sample_temporal_lr$equal_sister_reference <- llr_cfg$persistence_sibling_scope=="native"
      fork_contrasts <- llr_event_methods$SummarizeForkDynamics(
        fork_contrasts, sample_temporal_lr, lineage_paths,
        expression_mean_mode = event_expression_mean_mode,
        sample_labels = sample_temporal_lr$scores$biological_sample,
        condition_on_sample = component_ablation != "no_sample_temporal_conditioning")
      switch(as.character(shared_parent_signed), "TRUE" = {
        unadjusted_dynamics <- llr_event_methods$SummarizeForkDynamics(
          fork_contrasts, sample_temporal_lr, lineage_paths,
          expression_mean_mode = event_expression_mean_mode,
          sample_labels = sample_temporal_lr$scores$biological_sample,
          condition_on_sample = FALSE)
      })
    })
    switch(as.character(enable_structural_event_rank &&
                        nrow(fork_contrasts) > 0L), "TRUE" = {
      tree_nodes <- as.data.table(stack(setNames(
        lineage_paths, seq_along(lineage_paths))))
      setnames(tree_nodes, c("values", "ind"), c("node", "lineage"))
      tree_nodes[, position := seq_len(.N), by = lineage]
      tree_nodes[, event_node := shift(node), by = lineage]
      tree_edges <- unique(tree_nodes[!is.na(event_node),
        .(event_node, child = node, lineage, edge_position = position)])
      descendants <- unique(merge(tree_edges,
        tree_nodes[, .(lineage, position, node)],
        by = "lineage", allow.cartesian = TRUE)[
          position >= edge_position, .(event_node, child, node)])
      cluster_size <- data.table(node = lineage_cl)[
        , .(n_cells = .N), by = node]
      branch_sides <- merge(descendants, cluster_size, by = "node")[
        , .(side_cells = sum(n_cells)), by = .(event_node, child)]
      branch_sides[, event_cells := sum(side_cells), by = event_node]
      branch_sides[, balanced_cells := pmin(
        side_cells, event_cells - side_cells)]
      contrast_edge <- paste(fork_contrasts$event_node,
        fork_contrasts$child, sep = ":")
      support_edge <- paste(branch_sides$event_node,
        branch_sides$child, sep = ":")
      support_fraction <- branch_sides$balanced_cells[
        match(contrast_edge, support_edge)] /
        branch_sides$event_cells[match(contrast_edge, support_edge)]
      stopifnot(all(is.finite(support_fraction)))
      base_score <- fork_contrasts$persistent_positive_score
      rank_factor <- sqrt(support_fraction)
      fork_contrasts$legacy_event_ranking_score <- base_score * rank_factor
      fork_contrasts$pulse_weight <- -expm1(
        -data.table::fcoalesce(fork_contrasts$pulse_bic_gain, 0) / 2)
      fork_contrasts$focal_pulse_weight <- -expm1(
        -data.table::fcoalesce(fork_contrasts$focal_pulse_bic_gain, 0) / 2)
      branch_activation_weight <- -expm1(-data.table::fifelse(
        fork_contrasts$dynamic_slope_left > 0 & fork_contrasts$sibling_contrast > 0,
        data.table::fcoalesce(fork_contrasts$dynamic_bic_gain, 0), 0) / 2)
      fork_contrasts$joint_pulse_weight <- pmin(fork_contrasts$focal_pulse_weight,
        pmax(fork_contrasts$pulse_weight, data.table::fcoalesce(branch_activation_weight, 0)))
      switch(component_ablation, relative_pulse_only = {
        fork_contrasts$joint_pulse_weight <- pmin(
          fork_contrasts$pulse_weight, fork_contrasts$focal_pulse_weight)
      })
      switch(component_ablation, no_pulse_evidence = {
        fork_contrasts$joint_pulse_weight <- rep.int(0, nrow(fork_contrasts))
      })
      transient_or_persistent <- fork_contrasts$persistent_positive_score +
        fork_contrasts$joint_pulse_weight *
        (fork_contrasts$sibling_contrast - fork_contrasts$persistent_positive_score)
      fork_contrasts$event_ranking_score <- transient_or_persistent * rank_factor
      switch(as.character(shared_parent_signed), "TRUE" = {
        unadjusted_weight <- pmin(-expm1(-data.table::fcoalesce(
          unadjusted_dynamics$pulse_bic_gain, 0) / 2),
          -expm1(-data.table::fcoalesce(unadjusted_dynamics$focal_pulse_bic_gain, 0) / 2))
        unadjusted_activation <- -expm1(-data.table::fifelse(
          unadjusted_dynamics$dynamic_slope_left > 0 & fork_contrasts$sibling_contrast > 0,
          data.table::fcoalesce(unadjusted_dynamics$dynamic_bic_gain, 0), 0) / 2)
        switch(as.character(component_ablation != "relative_pulse_only"), "TRUE" = {
          unadjusted_weight <- pmax(unadjusted_weight, pmin(
            -expm1(-data.table::fcoalesce(unadjusted_dynamics$focal_pulse_bic_gain, 0) / 2),
            data.table::fcoalesce(unadjusted_activation, 0)))
        })
        switch(component_ablation, no_pulse_evidence = {
          unadjusted_weight <- rep.int(0, length(unadjusted_weight))
        })
        unadjusted_point <- fork_contrasts$persistent_positive_score +
          unadjusted_weight * (fork_contrasts$sibling_contrast - fork_contrasts$persistent_positive_score)
        fork_contrasts$sample_unadjusted_ranking_score <- unadjusted_point * rank_factor
      })
      focal_change <- fork_contrasts$post_fork_change
      if (paired_amplitude) {
        focal_late <- fork_contrasts$focal_late
      } else {
        focal_late <- data.table::fifelse(focal_change > 0 &
          is.finite(fork_contrasts$persistence_focal_activity) &
          fork_contrasts$n_late_bins >= 1L,
          pmin(focal_change, fork_contrasts$persistence_focal_activity), focal_change)
        switch(as.character(signed_persistence), "TRUE" = {
          focal_late_change <- fork_contrasts$persistence_focal_activity - fork_contrasts$focal_pre
          focal_late <- data.table::fifelse(is.finite(focal_late_change) &
            fork_contrasts$n_late_bins >= 1L, sign(focal_change) *
            pmin(abs(focal_change), pmax(0, sign(focal_change) * focal_late_change)), focal_change)
          fork_contrasts$focal_late_change <- focal_late_change
          fork_contrasts$focal_late <- focal_late
        })
      }
      focal_point <- focal_late +
        fork_contrasts$joint_pulse_weight * (focal_change - focal_late)
      fork_contrasts$no_fork_ranking_score <- focal_point * rank_factor
      # Keep the signed effect; shrink by native context error relative to
      # the median conditional error on the same amplitude estimand.
      switch(as.character(llr_cfg$context_error_weighting), "TRUE" = {
        rank_error_floor <- median(fork_contrasts$theta_measurement_se[
          fork_contrasts$evaluable %in% TRUE], na.rm = TRUE)
        rank_error_weight <- data.table::fifelse(fork_contrasts$theta_measurement_se == 0, 1,
          rank_error_floor / (rank_error_floor + fork_contrasts$theta_measurement_se))
        fork_contrasts$rank_error_floor <- rank_error_floor
        fork_contrasts$rank_error_weight <- rank_error_weight
        ranking_fields <- c("event_ranking_score", "no_fork_ranking_score",
          "sample_unadjusted_ranking_score", "legacy_event_ranking_score")
        fork_contrasts[, ranking_fields] <- as.matrix(fork_contrasts[, ranking_fields]) * rank_error_weight
      })
      fork_contrasts$ranking_estimand <-
        paste(fork_contrasts$amplitude_unit,"signed_branch_change",sep="_")
    })
    receiver_samples <- switch(as.character(!is.na(fork_sample_column)),
      "TRUE" = as.character(meta[[fork_sample_column]]),
      "FALSE" = rep("pooled", length(celltypes)))
    receiver_counts <- data.table::data.table(target = celltypes, sample = receiver_samples)[
      !is.na(sample) & nzchar(sample), .(n = .N), by = .(target, sample)]
    receiver_totals <- receiver_counts[, .(total = sum(n)), by = sample]
    receiver_grid <- merge(data.table::CJ(target = unique(celltypes), sample = receiver_totals$sample),
      receiver_totals, by = "sample")
    receiver_grid <- merge(receiver_grid, receiver_counts, by = c("target", "sample"), all.x = TRUE)
    default_receiver_groups <- receiver_grid[, .(
      weight = mean(data.table::fcoalesce(n, 0L) / total)), by = target][, context := "All receivers"]
    fork_evidence <- llr_event_methods$ComposeForkEvidence(
      fork_contrasts, sample_support = fork_sample_summary,
      nomination_groups = switch(as.character(main_candidate),
        "TRUE" = default_receiver_groups, "FALSE" = NULL),
      pair_gene_map = lr_db$pair_gene_map
    )
  }

  # Calendar dynamics use observed experimental time and independent samples.
  # They do not replace the native state contrast or its LR nomination score.
  calendar_fork_dynamics <- calendar_temporal_lr <- NULL
  calendar_inference_status <- "no_registered_calendar_time"
  calendar_time_column <- CONFIG$dataset_calendar_time_columns[ds_name]
  calendar_time_unit <- CONFIG$dataset_calendar_time_units[ds_name]
  calendar_sample_column <- fork_sample_column
  calendar_sampling_unit <- "biological_sample"
  calendar_observation_column <- CONFIG$dataset_calendar_observation_columns[ds_name]
  switch(as.character(length(calendar_observation_column) == 1L &&
    !is.na(calendar_observation_column)), "TRUE" = {
    calendar_sample_column <- unname(calendar_observation_column)
    calendar_sampling_unit <- "experimental_block"
  })
  switch(as.character(main_candidate && length(calendar_time_column) == 1L &&
    !is.na(calendar_time_column) && !is.na(calendar_sample_column) &&
    !is.null(temporal_lr) && !is.null(fork_contrasts)), "TRUE" = {
    calendar_samples <- unique(data.table(
      biological_sample = as.character(meta[[calendar_sample_column]]),
      observed_time = as.numeric(meta[[calendar_time_column]])))
    stopifnot(all(is.finite(calendar_samples$observed_time)),
      calendar_samples[, all(.N == 1L), by = biological_sample]$V1)
    calendar_times <- sort(unique(calendar_samples$observed_time))
    calendar_inference_status <- "insufficient_observed_times"
    switch(as.character(length(calendar_times) >= 2L), "TRUE" = {
      calendar_condition_column <- CONFIG$dataset_calendar_condition_columns[ds_name]
      calendar_conditions <- rep.int("all", length(celltypes))
      switch(as.character(length(calendar_condition_column) == 1L &&
        !is.na(calendar_condition_column)), "TRUE" = {
        calendar_conditions <- as.character(meta[[calendar_condition_column]])
      })
      calendar_temporal_lr <- LineaLinkR::ComputeObservedCalendarLR(
        expr_lr, celltypes, lr_db, as.numeric(meta[[calendar_time_column]]),
        as.character(meta[[calendar_sample_column]]), lineage_cl, lineage_paths,
        hill_scale = temporal_lr$intermediate$hill_scale,
        min_pct = llr_cfg$min_pct %||% 0.05,
        min_n_eff = llr_cfg$min_n_eff %||% 5.0,
        requested_contexts = unique(as.data.table(fork_contrasts)[, .(pair_idx, source, target)]),
        condition_labels = calendar_conditions)
      calendar_scores <- as.data.table(calendar_temporal_lr$scores)
      calendar_fork_dynamics <- llr_event_methods$SummarizeForkDynamics(
        fork_contrasts, calendar_temporal_lr, lineage_paths,
        expression_mean_mode = event_expression_mean_mode,
        sample_labels = calendar_scores$biological_sample,
        time_coordinate = "observed_calendar",
        condition_labels = calendar_scores$condition,
        reference_conditions = unname(na.omit(CONFIG$dataset_calendar_reference_conditions[ds_name])))
      attr(calendar_fork_dynamics, "calendar_time_unit") <- unname(calendar_time_unit)
      attr(calendar_fork_dynamics, "sampling_unit") <- calendar_sampling_unit
      calendar_temporal_lr$sampling_unit <- calendar_sampling_unit
      calendar_inference_status <- ifelse(calendar_sampling_unit == "biological_sample",
        "observed_sample_profiles", "observed_experimental_profiles")
    })
  })

