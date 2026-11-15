pub mod evidence;
pub mod fake_lr;
pub mod lr_scoring;
pub mod rank;
pub mod taxonomy;
pub mod temporal_null;
pub mod trajectory_curve;

pub use lr_scoring::{
    cluster_stats_csc, compute_lr_scores, compute_lr_scores_bins, LrBinResult, LrPair, LrResult,
};

pub use taxonomy::{
    compute_lineage_entropy_raw, compute_role_scores_one, find_bifurcation_bins, smooth_ma3,
    BifurcationReason, BifurcationResult, LineageEntropyTensor, LrBinIntermediate, Role,
    RoleScoreParams, RoleScores,
};

pub use fake_lr::{
    calibrate_fake_trajectories, compute_gene_bin_stats, compute_gene_features,
    compute_hill_s_from_stats, compute_p_fake_by_position, construct_fake_pairs_context_aware,
    construct_fake_pairs_stratified, score_fake_pairs_bins, ComplexNullMode, FakeLrBinResult,
    FakeLrPair, FakeMatchingBlocks, FakeNullMode, FakeTrajectoryCalibration, GeneBinStats,
    GeneFeatures, FAKE_NULL_GENERIC, FAKE_NULL_NETWORK_REWIRE_LIGAND,
    FAKE_NULL_NETWORK_REWIRE_RECEPTOR, FAKE_NULL_REWIRE_LIGAND, FAKE_NULL_REWIRE_RECEPTOR,
};

pub use evidence::{evidence_decomp_batch, EvidenceDecomp};

pub use temporal_null::{
    compute_temporal_permutation, TemporalPermutationBin, TemporalPermutationTrajectory,
};
