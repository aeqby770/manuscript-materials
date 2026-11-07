# Auto-style extendr .Call shim for LineaLinkR.
#
# Each `rs_*` function below is a thin R wrapper around the native symbol
# exported by `extendr_module! { mod LineaLinkR; fn rs_*; }` in
# src/extendr_bindings.rs.  The argument *order* and *count* must match
# the corresponding #[extendr] signature exactly.
#
# `wrap__rs_*` symbols are native objects made available by
#   useDynLib(LineaLinkR, .registration = TRUE)
# in NAMESPACE; .Call() takes the symbol object directly (no quotes,
# no PACKAGE = ...).
#
# Hand-written .Call shims (no rextendr required at install time).
# See BUILD_NOTES.md for the build matrix.
#
# nolint start

rs_compute_pseudotime_fuzzy <- function(nn_idx, nn_dist, n_cells, k, start_cell) {
  .Call(wrap__rs_compute_pseudotime_fuzzy,
        nn_idx, nn_dist, n_cells, k, start_cell)
}

rs_compute_pseudotime_fuzzy_multi <- function(nn_idx, nn_dist, n_cells, k, start_cells) {
  .Call(wrap__rs_compute_pseudotime_fuzzy_multi,
        nn_idx, nn_dist, n_cells, k, start_cells)
}

rs_hnsw_cost_edges <- function(nn_idx, nn_dist, n_cells, k) {
  .Call(wrap__rs_hnsw_cost_edges, nn_idx, nn_dist, n_cells, k)
}

rs_knn_component_count <- function(nn_idx, n_cells, k) {
  .Call(wrap__rs_knn_component_count, nn_idx, n_cells, k)
}

rs_knn_component_labels <- function(nn_idx, n_cells, k) {
  .Call(wrap__rs_knn_component_labels, nn_idx, n_cells, k)
}

rs_attractor_unreachable_cells <- function(nn_idx, nn_dist, terminal_cells, n_cells, k) {
  .Call(wrap__rs_attractor_unreachable_cells, nn_idx, nn_dist, terminal_cells, n_cells, k)
}

rs_compute_attractor_fate_f64 <- function(nn_idx, nn_dist, terminal_cells,
                                          n_cells, k, max_iter, tol) {
  .Call(wrap__rs_compute_attractor_fate_f64,
        nn_idx, nn_dist, terminal_cells, n_cells, k, max_iter, tol)
}

rs_sparse_absorption_lu64 <- function(system, rhs) {
  matrix(.Call(wrap__rs_sparse_absorption_lu64,
               system, as.numeric(rhs), as.integer(ncol(rhs))),
         nrow = nrow(rhs), ncol = ncol(rhs))
}

rs_extract_lineages <- function(embedding, cluster_labels, n_cells, n_dims, n_clusters, start_cluster) {
  .Call(wrap__rs_extract_lineages,
        embedding, cluster_labels, n_cells, n_dims, n_clusters, start_cluster)
}

rs_compute_lineage_weights <- function(embedding, cluster_labels, n_cells, n_dims, n_clusters, start_cluster, sigma) {
  .Call(wrap__rs_compute_lineage_weights,
        embedding, cluster_labels, n_cells, n_dims, n_clusters, start_cluster, sigma)
}

rs_compute_lineage_weights_trunk_aware <- function(embedding, cluster_labels, n_cells, n_dims, n_clusters, start_cluster, sigma) {
  .Call(wrap__rs_compute_lineage_weights_trunk_aware,
        embedding, cluster_labels, n_cells, n_dims, n_clusters, start_cluster, sigma)
}

rs_extract_lineages_from_edges <- function(edges_from, edges_to, n_clusters, start_cluster) {
  .Call(wrap__rs_extract_lineages_from_edges,
        edges_from, edges_to, n_clusters, start_cluster)
}

rs_compute_lineage_weights_trunk_aware_from_edges <- function(embedding, cluster_labels, n_cells, n_dims, n_clusters, edges_from, edges_to, start_cluster, sigma) {
  .Call(wrap__rs_compute_lineage_weights_trunk_aware_from_edges,
        embedding, cluster_labels, n_cells, n_dims, n_clusters, edges_from, edges_to,
        start_cluster, sigma)
}

rs_compute_lineage_weights_local_branch_residual <- function(embedding, cluster_labels, n_cells, n_dims, n_clusters, lineage_paths, priming) {
  .Call(wrap__rs_compute_lineage_weights_local_branch_residual,
        embedding, cluster_labels, n_cells, n_dims, n_clusters, lineage_paths, priming)
}

rs_compute_fate_probabilities <- function(nn_idx, nn_dist, pseudotime, terminal_cells, n_cells, k, bandwidth, max_iter, tol) {
  .Call(wrap__rs_compute_fate_probabilities,
        nn_idx, nn_dist, pseudotime, terminal_cells, n_cells, k, bandwidth, max_iter, tol)
}

rs_compute_lr_scores <- function(expr, cell_types, n_types, ligand_rows_list, receptor_rows_list, n_perms, seed, min_pct) {
  .Call(wrap__rs_compute_lr_scores,
        expr, cell_types, n_types, ligand_rows_list, receptor_rows_list, n_perms, seed, min_pct)
}

rs_compute_temporal_lr <- function(expr, cell_types, n_types, ligand_rows_list, receptor_rows_list, pseudotime, lineage_weights, n_lineages, n_bins, min_pct, min_n_eff, return_complete_grid, req_pair_idx, req_source, req_target, common_grid) {
  .Call(wrap__rs_compute_temporal_lr,
        expr, cell_types, n_types, ligand_rows_list, receptor_rows_list,
        pseudotime, lineage_weights, n_lineages, n_bins, min_pct,
        min_n_eff, return_complete_grid, req_pair_idx, req_source, req_target,
        common_grid)
}

rs_compute_temporal_lr_with_inter <- function(expr, cell_types, n_types, ligand_rows_list, receptor_rows_list, pseudotime, lineage_weights, n_lineages, n_bins, min_pct, min_n_eff, return_complete_grid, req_pair_idx, req_source, req_target, common_grid) {
  .Call(wrap__rs_compute_temporal_lr_with_inter,
        expr, cell_types, n_types, ligand_rows_list, receptor_rows_list,
        pseudotime, lineage_weights, n_lineages, n_bins, min_pct,
        min_n_eff, return_complete_grid, req_pair_idx, req_source, req_target,
        common_grid)
}

rs_compute_bin_bounds <- function(pseudotime, lineage_weights, n_lineages, n_bins) {
  .Call(wrap__rs_compute_bin_bounds,
        pseudotime, lineage_weights, n_lineages, n_bins)
}

rs_compute_common_bin_bounds <- function(pseudotime, lineage_weights, n_lineages, n_bins) {
  .Call(wrap__rs_compute_common_bin_bounds,
        pseudotime, lineage_weights, n_lineages, n_bins)
}

rs_compute_temporal_lr_fixed_bounds <- function(expr, cell_types, n_types, ligand_rows_list, receptor_rows_list, pseudotime, lineage_weights, n_lineages, n_bins, bin_bounds, hill_scale, min_pct, min_n_eff, return_complete_grid, req_pair_idx, req_source, req_target) {
  .Call(wrap__rs_compute_temporal_lr_fixed_bounds,
        expr, cell_types, n_types, ligand_rows_list, receptor_rows_list,
        pseudotime, lineage_weights, n_lineages, n_bins, bin_bounds,
        hill_scale, min_pct, min_n_eff, return_complete_grid,
        req_pair_idx, req_source, req_target)
}

rs_compute_temporal_permutation <- function(expr, cell_types, n_types,
                                            ligand_rows_list, receptor_rows_list,
                                            pseudotime, lineage_weights, block_labels,
                                            n_lineages, n_bins, min_pct,
                                            min_n_eff, min_valid_bins,
                                            min_valid_fraction,
                                            n_perms, seed, stat_code, top_k,
                                            permutation_stratum,
                                            req_pair_idx, req_source, req_target,
                                            req_lineage) {
  .Call(wrap__rs_compute_temporal_permutation,
        expr, cell_types, n_types,
        ligand_rows_list, receptor_rows_list,
        pseudotime, lineage_weights, block_labels,
        n_lineages, n_bins, min_pct,
        min_n_eff, min_valid_bins, min_valid_fraction,
        n_perms, seed, stat_code, top_k, permutation_stratum,
        req_pair_idx, req_source, req_target, req_lineage)
}

rs_compute_lineage_entropy_raw <- function(hwsums, wsums, n_types, n_lineages, n_bins) {
  .Call(wrap__rs_compute_lineage_entropy_raw,
        hwsums, wsums, n_types, n_lineages, n_bins)
}

rs_smooth_ma3 <- function(h_raw, n_lineages, n_bins) {
  .Call(wrap__rs_smooth_ma3,
        h_raw, n_lineages, n_bins)
}

rs_find_bifurcation_bins <- function(h_tilde, n_lineages, n_bins, h_fork) {
  .Call(wrap__rs_find_bifurcation_bins,
        h_tilde, n_lineages, n_bins, h_fork)
}

rs_compute_role_scores <- function(row_bin, row_lineage, spec_self, spec_max_other,
                                    delta_tau, row_n_eff_cells,
                                    lineage_b_star, lineage_reason_code,
                                    tau_0, h_fork, min_n_eff_cells,
                                    transition_window = 0L) {
  .Call(wrap__rs_compute_role_scores,
        row_bin, row_lineage, spec_self, spec_max_other,
        delta_tau, row_n_eff_cells,
        lineage_b_star, lineage_reason_code,
        tau_0, h_fork, min_n_eff_cells, transition_window)
}

rs_trajectory_de <- function(expr_flat, pseudotime, n_genes, n_cells, min_cells) {
  .Call(wrap__rs_trajectory_de,
        expr_flat, pseudotime, n_genes, n_cells, min_cells)
}

rs_compute_fake_lr_null <- function(expr, cell_types, n_types,
                                    candidate_fake_genes,
                                    real_ligand_rows_list, real_receptor_rows_list,
                                    pseudotime, lineage_weights,
                                    n_lineages, n_bins, min_pct,
                                    n_fake_per_real, n_strata, caliper_k,
                                    max_match_distance, seed,
                                    return_bin_scores, stat_code, top_k,
                                    min_n_eff, min_valid_bins, min_valid_fraction,
                                    fake_null_mode, complex_null_mode, block_labels) {
  .Call(wrap__rs_compute_fake_lr_null,
        expr, cell_types, n_types,
        candidate_fake_genes,
        real_ligand_rows_list, real_receptor_rows_list,
        pseudotime, lineage_weights,
        n_lineages, n_bins, min_pct,
        n_fake_per_real, n_strata, caliper_k, max_match_distance, seed,
        return_bin_scores, stat_code, top_k,
        min_n_eff, min_valid_bins, min_valid_fraction, fake_null_mode,
        complex_null_mode, block_labels)
}

rs_evidence_decomp <- function(real_scores_pair_idx, real_scores_source,
                                real_scores_target, real_scores_lineage,
                                real_scores_bin, real_scores_n_eff,
                                expr, cell_types, n_types,
                                ligand_rows_list, receptor_rows_list,
                                pseudotime, lineage_weights,
                                n_lineages, n_bins, min_pct) {
  .Call(wrap__rs_evidence_decomp,
        real_scores_pair_idx, real_scores_source,
        real_scores_target, real_scores_lineage,
        real_scores_bin, real_scores_n_eff,
        expr, cell_types, n_types,
        ligand_rows_list, receptor_rows_list,
        pseudotime, lineage_weights,
        n_lineages, n_bins, min_pct)
}

# nolint end
