//! LineaLinkR: trajectory-linked ligand-receptor communication analysis.

#[path = "../ccc/mod.rs"]
pub mod ccc;
#[path = "../core/mod.rs"]
pub mod core;
#[path = "../trajectory/mod.rs"]
pub mod trajectory;

// R (extendr) bridge. Compiled only when the `extendr` feature is enabled
// (R package build via Makevars). Plain `cargo test --lib` does not need extendr.
#[cfg(feature = "extendr")]
mod extendr_bindings;

pub use ccc::{
    calibrate_fake_trajectories, cluster_stats_csc, compute_gene_bin_stats, compute_gene_features,
    compute_hill_s_from_stats, compute_lineage_entropy_raw, compute_lr_scores,
    compute_lr_scores_bins, compute_p_fake_by_position, compute_role_scores_one,
    compute_temporal_permutation, construct_fake_pairs_stratified, evidence_decomp_batch,
    find_bifurcation_bins, score_fake_pairs_bins, smooth_ma3, BifurcationReason, BifurcationResult,
    EvidenceDecomp, FakeLrBinResult, FakeLrPair, FakeTrajectoryCalibration, GeneBinStats,
    GeneFeatures, LineageEntropyTensor, LrBinIntermediate, LrBinResult, LrPair, LrResult, Role,
    RoleScoreParams, RoleScores, TemporalPermutationBin, TemporalPermutationTrajectory,
};
pub use core::SparseMatrix;
pub use trajectory::{
    cluster_centroids, cluster_lineage_count, compute_lineage_weights,
    compute_lineage_weights_trunk_aware, compute_lineage_weights_local_branch_residual,
    compute_pseudotime, compute_pseudotime_from_hnsw_fuzzy,
    compute_pseudotime_from_hnsw_fuzzy_multi, compute_pseudotime_fuzzy, count_knn_components,
    directed_markov, dtw_similarity, extract_lineages, fate_probabilities, knn_component_labels,
    lineage_entropy, mst_prim, pseudotime_de,
};
