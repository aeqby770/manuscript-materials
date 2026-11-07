//! Stratified fake-LR null model for score calibration.
//!
//! ## Purpose
//!
//! Constructs a null distribution of LR scores from gene pairs that are
//! **not** known to form ligand-receptor interactions, matched to real
//! LR genes on marginal expression features.  The resulting p_fake
//! values test whether LineaLinkR's scoring responds to LR **identity**
//! (beyond marginal co-expression), complementary to the temporal permutation
//! null which tests whether trajectory structure exceeds random time ordering.
//!
//! ## Terminology
//!
//! - `G_fake` — candidate genes eligible for fake pairing (expressed but
//!   absent from all LR databases + STRING/BioGRID). Supplied by the caller.
//! - **stratum** — a cell in the 4-D discretisation grid
//!   (mean × pct × Tau × trend). Each axis is min–max normalised and split into
//!   `n_strata` equal-width bins; within-stratum nearest-neighbour sampling
//!   follows the caliper rule below.
//! - **caliper** — within-stratum nearest-neighbour selection: for each
//!   real LR gene, pick the `caliper_k` closest candidate genes (by
//!   Euclidean distance in the 4-D feature space) in its own stratum.
//!
//! ## Design decisions
//!
//! - Gene-pair swap (a) is the **primary null** — structure-guaranteed
//!   invalid, testing LR identity specificity head-on.
//! - Stratified sampling + within-stratum caliper (a2), NOT propensity
//!   matching (a3) — semantic overlap would damage interpretability.
//! - `|D*| = n_fake_mult × |D|` target.
//! - Row-level p_fake uses the actually observed matched fake rows.
//! - Publication inference should aggregate bins to fake trajectories and use
//!   enough fake pairs for the requested p-value resolution.
//!
//! ## Relationship with `lr_scoring.rs`
//!
//! This module **does not** re-run the full CSC pass for fake genes.
//! Instead the caller pre-computes per-gene bin stats (mean/pct per
//! (gene, ct, lin, bin)) via the `compute_gene_bin_stats` entry point,
//! which runs the CSC walk once for the candidate fake gene list.
//! Subsequent fake-pair scoring walks the pre-computed lookup tables —
//! O(|D*| × n_types² × n_lineages × n_bins), no cell walk.
//!
//! The publication path streams trajectory statistics and does not materialise
//! the full fake-bin table.

use crate::ccc::rank::spearman_f32;
use crate::core::SparseMatrix;
use rand::{rngs::StdRng, Rng, SeedableRng};
use std::collections::HashSet;

// ─── Per-gene bin statistics ───────────────────────────────────────────

/// Per-(gene, cell_type, lineage, bin) statistics kept for fast fake-pair
/// scoring.  Layout: `gene * (n_types * n_lineages * n_bins) + ct * slb + lin * sb + bin`.
#[derive(Debug, Clone)]
pub struct GeneBinStats {
    /// Weighted means: E[x_g | ct, lin, bin] under lineage weight Ψ_l(i).
    pub means: Vec<f32>,
    /// Weighted pct expressed: fraction of weighted cells with x > 0.
    pub pcts: Vec<f32>,
    /// Total weight per (ct, lin, bin).
    pub wsums: Vec<f64>,
    /// Sum of squared weights per (ct, lin, bin), used for Kish ESS.
    pub w2sums: Vec<f64>,
    /// Global 0-based row indices; local gene index is the vector position.
    pub gene_rows: Vec<usize>,
    pub n_genes: usize,
    pub n_types: usize,
    pub n_lineages: usize,
    pub n_bins: usize,
}

impl GeneBinStats {
    #[inline]
    pub fn local_gene(&self, global_gene: usize) -> Option<usize> {
        self.gene_rows.iter().position(|&g| g == global_gene)
    }

    #[inline]
    pub fn get_mean_global(&self, global_gene: usize, ct: usize, lin: usize, bin: usize) -> f32 {
        self.local_gene(global_gene)
            .map(|local| self.get_mean(local, ct, lin, bin))
            .unwrap_or(0.0)
    }

    #[inline]
    pub fn get_pct_global(&self, global_gene: usize, ct: usize, lin: usize, bin: usize) -> f32 {
        self.local_gene(global_gene)
            .map(|local| self.get_pct(local, ct, lin, bin))
            .unwrap_or(0.0)
    }

    /// Resolve gene expression from `primary` when present, else `secondary`.
    #[inline]
    pub fn get_mean_resolved(
        primary: &GeneBinStats,
        secondary: &GeneBinStats,
        global_gene: usize,
        ct: usize,
        lin: usize,
        bin: usize,
    ) -> f32 {
        if primary.local_gene(global_gene).is_some() {
            primary.get_mean_global(global_gene, ct, lin, bin)
        } else {
            secondary.get_mean_global(global_gene, ct, lin, bin)
        }
    }

    /// Resolve detection fraction from `primary` when present, else `secondary`.
    #[inline]
    pub fn get_pct_resolved(
        primary: &GeneBinStats,
        secondary: &GeneBinStats,
        global_gene: usize,
        ct: usize,
        lin: usize,
        bin: usize,
    ) -> f32 {
        if primary.local_gene(global_gene).is_some() {
            primary.get_pct_global(global_gene, ct, lin, bin)
        } else {
            secondary.get_pct_global(global_gene, ct, lin, bin)
        }
    }

    /// Weighted mean across cell types for one (lineage, bin), using `wsums`.
    pub fn get_mean_lin_bin_weighted_across_types(
        &self,
        global_gene: usize,
        lin: usize,
        bin: usize,
    ) -> f32 {
        let mut weighted_sum = 0.0f64;
        let mut weight_mass = 0.0f64;
        for ct in 0..self.n_types {
            let w = self.wsums[ct * (self.n_lineages * self.n_bins) + lin * self.n_bins + bin];
            if w > 1e-8 {
                weighted_sum += self.get_mean_global(global_gene, ct, lin, bin) as f64 * w;
                weight_mass += w;
            }
        }
        if weight_mass > 1e-8 {
            (weighted_sum / weight_mass) as f32
        } else {
            f32::NAN
        }
    }

    #[inline]
    pub fn idx(&self, gene: usize, ct: usize, lin: usize, bin: usize) -> usize {
        gene * (self.n_types * self.n_lineages * self.n_bins)
            + ct * (self.n_lineages * self.n_bins)
            + lin * self.n_bins
            + bin
    }

    #[inline]
    pub fn get_mean(&self, gene: usize, ct: usize, lin: usize, bin: usize) -> f32 {
        self.means[self.idx(gene, ct, lin, bin)]
    }

    #[inline]
    pub fn get_pct(&self, gene: usize, ct: usize, lin: usize, bin: usize) -> f32 {
        self.pcts[self.idx(gene, ct, lin, bin)]
    }
}

/// Compute per-gene weighted binned statistics for a specified set of
/// gene row indices (0-based).  One CSC walk — same algorithm as
/// `lr_scoring.rs:compute_lr_scores_bins` but outputs raw per-gene
/// arrays instead of LR pair scores.
///
/// `lineage_weights`: flat [n_cells × n_lineages], row-major.
/// When `block_labels` and `focus_block` are set, only cells in that block contribute.
pub fn compute_gene_bin_stats(
    expr: &SparseMatrix,
    cell_types: &[usize],
    n_types: usize,
    genes: &[usize],
    pseudotime: &[f32],
    lineage_weights: &[f32],
    n_lineages: usize,
    n_bins: usize,
    block_labels: Option<&[usize]>,
    focus_block: Option<usize>,
) -> GeneBinStats {
    let n_cells = expr.ncols;
    let mut seen = HashSet::with_capacity(genes.len());
    let gene_rows: Vec<usize> = genes
        .iter()
        .copied()
        .filter(|gene| seen.insert(*gene))
        .collect();
    let n_local = gene_rows.len();

    if n_cells == 0 || n_lineages == 0 || n_bins == 0 || n_local == 0 {
        return GeneBinStats {
            means: Vec::new(),
            pcts: Vec::new(),
            wsums: Vec::new(),
            w2sums: Vec::new(),
            gene_rows,
            n_genes: n_local,
            n_types,
            n_lineages,
            n_bins,
        };
    }

    // ── Equal-frequency bin boundaries (reuses lr_scoring::compute_bin_bounds) ─
    let bin_bounds = crate::ccc::lr_scoring::compute_bin_bounds(
        n_cells,
        pseudotime,
        lineage_weights,
        n_lineages,
        n_bins,
    );

    // O(1) lookup: global row → local index
    let mut g2local = vec![usize::MAX; expr.nrows];
    for (local, &global) in gene_rows.iter().enumerate() {
        if global < expr.nrows {
            g2local[global] = local;
        }
    }

    // ── Accumulators ────────────────────────────────────────────────────
    let sl = n_types * n_lineages * n_bins;
    let slb = n_lineages * n_bins;
    let sb = n_bins;
    let total = n_local * sl;
    let mut sums = vec![0.0f64; total];
    let mut xwsums = vec![0.0f64; total];
    let mut wsums_all = vec![0.0f64; n_types * n_lineages * n_bins];
    let mut w2sums_all = vec![0.0f64; n_types * n_lineages * n_bins];

    for cell in 0..n_cells {
        if let (Some(labels), Some(block)) = (block_labels, focus_block) {
            if labels.get(cell) != Some(&block) {
                continue;
            }
        }
        let ct = cell_types[cell];
        if ct >= n_types {
            continue;
        }
        let pt_raw = pseudotime[cell];
        if pt_raw.is_nan() {
            continue;
        }
        let pt = pt_raw.clamp(0.0, 1.0);

        for lin in 0..n_lineages {
            let w = lineage_weights[cell * n_lineages + lin] as f64;
            if w < 1e-8 {
                continue;
            }
            let bin = bin_bounds[lin]
                .partition_point(|&b| b <= pt)
                .saturating_sub(1)
                .min(n_bins - 1);
            let wsum_idx = ct * slb + lin * sb + bin;
            wsums_all[wsum_idx] += w;
            w2sums_all[wsum_idx] += w * w;
            for nnz in expr.col_ptrs[cell]..expr.col_ptrs[cell + 1] {
                let local = g2local[expr.rows[nnz]];
                if local == usize::MAX {
                    continue;
                }
                let base = local * sl + ct * slb + lin * sb + bin;
                let v = expr.values[nnz] as f64;
                sums[base] += v * w;
                if v > 0.0 {
                    xwsums[base] += w;
                }
            }
        }
    }

    // Normalise
    let means: Vec<f32> = (0..total)
        .map(|i| {
            let ct = (i % sl) / slb;
            let lin = (i % slb) / sb;
            let bin = i % sb;
            let ws = wsums_all[ct * slb + lin * sb + bin];
            if ws > 1e-8 {
                (sums[i] / ws) as f32
            } else {
                0.0
            }
        })
        .collect();
    let pcts: Vec<f32> = (0..total)
        .map(|i| {
            let ct = (i % sl) / slb;
            let lin = (i % slb) / sb;
            let bin = i % sb;
            let ws = wsums_all[ct * slb + lin * sb + bin];
            if ws > 1e-8 {
                (xwsums[i] / ws) as f32
            } else {
                0.0
            }
        })
        .collect();

    GeneBinStats {
        means,
        pcts,
        wsums: wsums_all,
        w2sums: w2sums_all,
        gene_rows,
        n_genes: n_local,
        n_types,
        n_lineages,
        n_bins,
    }
}

// ─── Per-gene features for stratification ─────────────────────────────

/// Scalar features per gene, used for stratified matching of real LR
/// genes against candidate fake genes.  All values are computed from
/// the expression matrix so the caller need not supply them.
#[derive(Debug, Clone)]
pub struct GeneFeatures {
    /// 0-based gene row index in the global expression matrix.
    pub gene_idx: usize,
    /// Mean expression per cell type [n_types].
    pub ct_mean: Vec<f32>,
    /// Pct expressed per cell type [n_types].
    pub ct_pct: Vec<f32>,
    /// Tau tissue specificity ∈ [0, 1] (higher = more specific).
    pub tau: f32,
    /// Signed pseudotime trend — mean Spearman ρ of mean expression against
    /// pseudotime-bin rank across lineages.
    pub trend: f32,
}

/// Compute per-gene scalar features from gene bin statistics.
/// `gene_indices` are 0-based global row indices.
///
/// Tau formula (Yanai et al. 2005):
///   tau = Σ_ct (1 - x_ct / max_ct) / (n_types - 1)
/// where x_ct is the mean expression per cell type.
///
/// Trend: computes mean per-bin expression per lineage, then signed Spearman ρ
/// with bin rank, averaged across lineages.
pub fn compute_gene_features(stats: &GeneBinStats, gene_indices: &[usize]) -> Vec<GeneFeatures> {
    let n_types = stats.n_types;
    let n_lineages = stats.n_lineages;
    let n_bins = stats.n_bins;

    gene_indices
        .iter()
        .map(|&g| {
            // Aggregate per-gene CT-level features from bin stats.
            // Weighted average across (lin, bin) to get per-CT summary.
            let mut ct_mean = vec![0.0f32; n_types];
            let mut ct_pct = vec![0.0f32; n_types];

            let mut ct_w = vec![0.0f64; n_types];
            for ct in 0..n_types {
                for lin in 0..n_lineages {
                    for bin in 0..n_bins {
                        let w = stats.wsums[ct * (n_lineages * n_bins) + lin * n_bins + bin];
                        if w > 1e-8 {
                            let wf = w as f32;
                            ct_mean[ct] += stats.get_mean_global(g, ct, lin, bin) * wf;
                            ct_pct[ct] += stats.get_pct_global(g, ct, lin, bin) * wf;
                            ct_w[ct] += w;
                        }
                    }
                }
                if ct_w[ct] > 1e-8 {
                    ct_mean[ct] /= ct_w[ct] as f32;
                    ct_pct[ct] /= ct_w[ct] as f32;
                }
            }

            // Tau specificity
            let max_mean = ct_mean.iter().cloned().fold(0.0f32, f32::max);
            let tau: f32 = if max_mean > 0.0 && n_types > 1 {
                let s: f32 = ct_mean.iter().map(|&m| 1.0 - m / max_mean).sum();
                s / (n_types as f32 - 1.0)
            } else {
                0.0
            };

            // Signed trend: mean Spearman ρ of weighted (bin, expr) per lineage.
            let mut trend_sum = 0.0f32;
            let mut trend_n = 0usize;
            for lin in 0..n_lineages {
                let bin_vals: Vec<(usize, f32)> = (0..n_bins)
                    .filter_map(|bin| {
                        let m = stats.get_mean_lin_bin_weighted_across_types(g, lin, bin);
                        if m.is_finite() {
                            Some((bin, m))
                        } else {
                            None
                        }
                    })
                    .collect();
                if bin_vals.len() < 3 {
                    continue;
                }
                let bin_indices: Vec<f32> = bin_vals.iter().map(|(b, _)| *b as f32).collect();
                let values: Vec<f32> = bin_vals.iter().map(|(_, v)| *v).collect();
                let rho = spearman_f32(&bin_indices, &values);
                if rho.is_finite() {
                    trend_sum += rho;
                    trend_n += 1;
                }
            }
            let trend = if trend_n > 0 {
                trend_sum / trend_n as f32
            } else {
                0.0
            };

            GeneFeatures {
                gene_idx: g,
                ct_mean,
                ct_pct,
                tau,
                trend,
            }
        })
        .collect()
}

// ─── Stratified fake-pair construction ────────────────────────────────

/// Kind of fake-null draw used for calibration.
pub const FAKE_NULL_GENERIC: u8 = 0;
pub const FAKE_NULL_REWIRE_RECEPTOR: u8 = 1;
pub const FAKE_NULL_REWIRE_LIGAND: u8 = 2;
pub const FAKE_NULL_NETWORK_REWIRE_RECEPTOR: u8 = 3;
pub const FAKE_NULL_NETWORK_REWIRE_LIGAND: u8 = 4;
/// Sentinel coordinates for pooled fake pairs that are not context-scoped.
pub const FAKE_CTX_WILDCARD: usize = usize::MAX;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FakeNullMode {
    ContextGeneric,
    EdgeRewire,
    Both,
    /// LR-network rewiring: replace one endpoint from lr_db ligand/receptor pools.
    LrNetworkRewire,
    /// All three null layers: generic + edge (ligand/receptor pooled) + network.
    All,
}

fn split_null_budget_layers(total: usize, mode: FakeNullMode) -> (usize, usize, usize) {
    match mode {
        FakeNullMode::Both => {
            let generic = total / 2;
            let edge = total.saturating_sub(generic);
            (generic, edge, 0)
        }
        FakeNullMode::All => {
            let generic = total / 3;
            let remaining = total.saturating_sub(generic);
            let edge = remaining / 2;
            let network = remaining.saturating_sub(edge);
            (generic, edge, network)
        }
        FakeNullMode::ContextGeneric => (total, 0, 0),
        FakeNullMode::EdgeRewire => (0, total, 0),
        FakeNullMode::LrNetworkRewire => (0, 0, total),
    }
}

/// Unique ligand/receptor genes and canonical single-gene LR pairs from lr_db.
fn build_lr_network_pools(
    lr_pairs: &[crate::ccc::LrPair],
) -> (
    Vec<usize>,
    Vec<usize>,
    std::collections::HashSet<(usize, usize)>,
) {
    let mut ligand_set = std::collections::HashSet::new();
    let mut receptor_set = std::collections::HashSet::new();
    let mut real_lr_pairs = std::collections::HashSet::new();
    for pair in lr_pairs {
        for &lig in &pair.ligand_rows {
            ligand_set.insert(lig);
        }
        for &rec in &pair.receptor_rows {
            receptor_set.insert(rec);
        }
        if pair.ligand_rows.len() == 1 && pair.receptor_rows.len() == 1 {
            let lig = pair.ligand_rows[0];
            let rec = pair.receptor_rows[0];
            real_lr_pairs.insert((lig, rec));
        }
    }
    let ligand_pool: Vec<usize> = ligand_set.into_iter().collect();
    let receptor_pool: Vec<usize> = receptor_set.into_iter().collect();
    (ligand_pool, receptor_pool, real_lr_pairs)
}

/// Split endpoint-replacement budget across ligand vs receptor rewires.
/// Even draw indices replace the receptor; odd indices replace the ligand.
fn split_edge_budget(total: usize) -> (usize, usize) {
    (total / 2, (total + 1) / 2)
}

/// How multi-subunit LR pairs participate in FakeLR null construction.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ComplexNullMode {
    /// Skip complex pairs (default).
    Skip,
    /// Draw cardinality-matched fake subunit sets across generic, edge, and
    /// network null layers (max 8 subunits/side).
    CardinalityMatched,
}

pub fn is_complex_lr_pair(pair: &crate::ccc::LrPair) -> bool {
    pair.ligand_rows.len() > 1 || pair.receptor_rows.len() > 1
}

/// Maximum subunits per ligand/receptor side for cardinality-matched complex null.
pub const MAX_COMPLEX_CARDINALITY: usize = 8;

fn complex_cardinality_supported(pair: &crate::ccc::LrPair) -> bool {
    !pair.ligand_rows.is_empty()
        && !pair.receptor_rows.is_empty()
        && pair.ligand_rows.len() <= MAX_COMPLEX_CARDINALITY
        && pair.receptor_rows.len() <= MAX_COMPLEX_CARDINALITY
}

#[cfg(test)]
fn fake_null_mode_supports_complex(mode: FakeNullMode) -> bool {
    matches!(
        mode,
        FakeNullMode::ContextGeneric
            | FakeNullMode::EdgeRewire
            | FakeNullMode::Both
            | FakeNullMode::LrNetworkRewire
            | FakeNullMode::All
    )
}

/// A sampled fake ligand-receptor pair, annotated with which real
/// LR pair it shadows (for downstream pair-grouped p_fake computation).
#[derive(Debug, Clone)]
pub struct FakeLrPair {
    pub ligand_gene: usize,
    pub receptor_gene: usize,
    /// When set, score using min-aggregation across all fake ligand subunits.
    pub complex_ligand_rows: Option<Vec<usize>>,
    /// When set, score using min-aggregation across all fake receptor subunits.
    pub complex_receptor_rows: Option<Vec<usize>>,
    pub shadowing_pair_idx: usize,
    pub source: usize,
    pub target: usize,
    pub lineage: usize,
    pub null_kind: u8,
    pub ligand_match_distance: f32,
    pub receptor_match_distance: f32,
}

impl FakeLrPair {
    pub fn is_complex(&self) -> bool {
        self.complex_ligand_rows.is_some()
    }

    pub fn ligand_rows(&self) -> &[usize] {
        self.complex_ligand_rows
            .as_deref()
            .unwrap_or(std::slice::from_ref(&self.ligand_gene))
    }

    pub fn receptor_rows(&self) -> &[usize] {
        self.complex_receptor_rows
            .as_deref()
            .unwrap_or(std::slice::from_ref(&self.receptor_gene))
    }

    pub fn draw_key(&self) -> (u64, u64) {
        if let (Some(lig), Some(rec)) = (&self.complex_ligand_rows, &self.complex_receptor_rows) {
            (hash_gene_set(lig), hash_gene_set(rec))
        } else {
            (self.ligand_gene as u64, self.receptor_gene as u64)
        }
    }
}

fn hash_gene_set(genes: &[usize]) -> u64 {
    let mut h = 0u64;
    for &g in genes {
        h = h.wrapping_mul(31).wrapping_add(g as u64);
    }
    h
}

/// A scored fake LR row.  `fake_idx` identifies which sampled fake pair
/// generated the row, so downstream code can aggregate bins into a fake
/// trajectory instead of pooling all fake rows for the same real LR pair.
#[derive(Debug, Clone)]
pub struct FakeLrBinResult {
    pub fake_idx: usize,
    pub pair_idx: usize,
    pub source: usize,
    pub target: usize,
    pub lineage: usize,
    pub bin: usize,
    pub score: f32,
    pub n_eff_cells: f32,
    pub weight_mass_cells: f32,
}

#[derive(Debug, Clone)]
pub struct FakeTrajectoryCalibration {
    pub pair_idx: usize,
    pub source: usize,
    pub target: usize,
    pub lineage: usize,
    pub statistic: f32,
    pub pvalue: f64,
    pub n_fake_expected: usize,
    pub n_fake_observed: usize,
    pub fake_coverage: f64,
    pub n_valid_bins: usize,
    pub valid_bin_fraction: f32,
    pub n_expression_supported_bins: usize,
    pub n_eff_at_max: f32,
    pub n_eff_at_min: f32,
    pub valid_at_max: bool,
    pub valid_at_min: bool,
    pub trajectory_evaluable: bool,
    pub pvalue_generic: f64,
    pub pvalue_edge: f64,
    pub pvalue_replace_ligand: f64,
    pub pvalue_replace_receptor: f64,
    pub mean_ligand_match_distance: f32,
    pub mean_receptor_match_distance: f32,
    pub n_unique_fake_draws: usize,
    pub n_unique_fake_draws_generic: usize,
    pub n_unique_fake_draws_edge: usize,
    pub n_fake_expected_generic: usize,
    pub n_fake_expected_edge: usize,
    pub n_fake_observed_generic: usize,
    pub n_fake_observed_edge: usize,
    pub fake_coverage_generic: f64,
    pub fake_coverage_edge: f64,
    pub n_fake_expected_replace_ligand: usize,
    pub n_fake_expected_replace_receptor: usize,
    pub n_fake_observed_replace_ligand: usize,
    pub n_fake_observed_replace_receptor: usize,
    pub fake_coverage_replace_ligand: f64,
    pub fake_coverage_replace_receptor: f64,
    pub n_unique_fake_draws_replace_ligand: usize,
    pub n_unique_fake_draws_replace_receptor: usize,
    pub pvalue_network: f64,
    pub pvalue_network_rewire_ligand: f64,
    pub pvalue_network_rewire_receptor: f64,
    pub n_fake_expected_network: usize,
    pub n_fake_observed_network: usize,
    pub fake_coverage_network: f64,
    pub n_unique_fake_draws_network: usize,
    pub n_fake_expected_network_rewire_ligand: usize,
    pub n_fake_expected_network_rewire_receptor: usize,
    pub n_fake_observed_network_rewire_ligand: usize,
    pub n_fake_observed_network_rewire_receptor: usize,
    pub fake_coverage_network_rewire_ligand: f64,
    pub fake_coverage_network_rewire_receptor: f64,
    pub n_unique_fake_draws_network_rewire_ligand: usize,
    pub n_unique_fake_draws_network_rewire_receptor: usize,
    /// Cardinality-matched generic null for multi-subunit LR pairs (NA for simple pairs).
    pub pvalue_complex: f64,
    /// Median fake-draw trajectory statistic across evaluable null layers.
    pub fake_null_stat_median: f32,
}

/// Context-aware 4-D feature vector for matching in `(cell_type, lineage)`.
pub fn context_feature_vec4(
    stats: &GeneBinStats,
    gene: usize,
    ct: usize,
    lin: usize,
    tau: f32,
) -> [f32; 4] {
    let n_bins = stats.n_bins;
    let slb = stats.n_lineages * n_bins;
    let sb = n_bins;
    let mut w_sum = 0.0f64;
    let mut mean_acc = 0.0f64;
    let mut pct_acc = 0.0f64;
    let mut bin_indices = Vec::new();
    let mut values = Vec::new();
    for bin in 0..n_bins {
        let w = stats.wsums[ct * slb + lin * sb + bin];
        if w > 1e-8 {
            let m = stats.get_mean_global(gene, ct, lin, bin);
            let p = stats.get_pct_global(gene, ct, lin, bin);
            mean_acc += m as f64 * w;
            pct_acc += p as f64 * w;
            w_sum += w;
            if m.is_finite() {
                bin_indices.push(bin as f32);
                values.push(m);
            }
        }
    }
    let (mean, pct) = if w_sum > 1e-8 {
        ((mean_acc / w_sum) as f32, (pct_acc / w_sum) as f32)
    } else {
        (0.0, 0.0)
    };
    let trend = if values.len() >= 3 {
        spearman_f32(&bin_indices, &values)
    } else {
        0.0
    };
    [mean, pct, tau, if trend.is_finite() { trend } else { 0.0 }]
}

fn bottleneck_gene_in_context(
    stats: &GeneBinStats,
    genes: &[usize],
    ct: usize,
    lin: usize,
    tau_by_gene: &std::collections::HashMap<usize, f32>,
) -> usize {
    genes
        .iter()
        .copied()
        .min_by(|&a, &b| {
            let ta = *tau_by_gene.get(&a).unwrap_or(&0.0);
            let tb = *tau_by_gene.get(&b).unwrap_or(&0.0);
            let va = context_feature_vec4(stats, a, ct, lin, ta)[0];
            let vb = context_feature_vec4(stats, b, ct, lin, tb)[0];
            va.partial_cmp(&vb).unwrap_or(std::cmp::Ordering::Equal)
        })
        .unwrap_or_else(|| genes[0])
}

fn top_k_nearest_context(
    target: &[f32; 4],
    pool: &[usize],
    candidate_genes: &[usize],
    stats: &GeneBinStats,
    ct: usize,
    lin: usize,
    tau_by_gene: &std::collections::HashMap<usize, f32>,
    k: usize,
    mins: &[f32],
    maxs: &[f32],
    max_match_distance: f32,
) -> Vec<(usize, f32)> {
    let mut dists: Vec<(f32, usize)> = pool
        .iter()
        .copied()
        .map(|idx| {
            let gene = candidate_genes[idx];
            let tau = *tau_by_gene.get(&gene).unwrap_or(&0.0);
            let cv = context_feature_vec4(stats, gene, ct, lin, tau);
            (norm_euclid_sq(target, &cv, mins, maxs), idx)
        })
        .filter(|(d2, _)| *d2 <= max_match_distance * max_match_distance)
        .collect();
    dists.sort_unstable_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    dists.truncate(k);
    dists
        .into_iter()
        .map(|(d2, idx)| (idx, d2.sqrt()))
        .collect()
}

fn collect_context_ranges(
    real_stats: &GeneBinStats,
    fake_stats: &GeneBinStats,
    fake_genes: &[usize],
    keys: &[crate::ccc::trajectory_curve::TrajKey],
    lr_pairs: &[crate::ccc::LrPair],
    contexts: &[(usize, usize)],
    tau_by_gene: &std::collections::HashMap<usize, f32>,
) -> (Vec<f32>, Vec<f32>) {
    let mut mins = [f32::INFINITY; 4];
    let mut maxs = [f32::NEG_INFINITY; 4];
    let mut consider = |gene: usize, ct: usize, lin: usize| {
        let stats = if real_stats.local_gene(gene).is_some() {
            real_stats
        } else {
            fake_stats
        };
        let v = context_feature_vec4(
            stats,
            gene,
            ct,
            lin,
            *tau_by_gene.get(&gene).unwrap_or(&0.0),
        );
        for (i, &x) in v.iter().enumerate() {
            if x < mins[i] {
                mins[i] = x;
            }
            if x > maxs[i] {
                maxs[i] = x;
            }
        }
    };
    for &(_pair_idx, source, target, lineage) in keys {
        let pair = &lr_pairs[_pair_idx];
        for &gene in pair.ligand_rows.iter().chain(pair.receptor_rows.iter()) {
            consider(gene, source, lineage);
            consider(gene, target, lineage);
        }
    }
    for &gene in fake_genes {
        for &(ct, lin) in contexts {
            consider(gene, ct, lin);
        }
    }
    for i in 0..4 {
        if !mins[i].is_finite() {
            mins[i] = 0.0;
        }
        if !maxs[i].is_finite() {
            maxs[i] = 1.0;
        }
    }
    (mins.to_vec(), maxs.to_vec())
}

fn feature_to_stratum_from_vec(
    vec: &[f32; 4],
    n_strata: usize,
    mins: &[f32],
    maxs: &[f32],
) -> usize {
    let ns = n_strata.max(1);
    let mut idx = 0usize;
    for i in 0..4 {
        let range = maxs[i] - mins[i];
        let s = if range > 1e-8 {
            ((vec[i] - mins[i]) / range * ns as f32) as usize
        } else {
            0
        };
        idx = idx * ns + s.min(ns - 1);
    }
    idx
}

/// Per-block context matching for donor/sample-aware FakeLR draws.
pub struct FakeMatchingBlocks<'a> {
    pub block_labels: &'a [usize],
    pub cell_types: &'a [usize],
    pub lineage_weights: &'a [f32],
    pub n_lineages: usize,
    pub n_blocks: usize,
    pub real_stats: &'a [GeneBinStats],
    pub fake_stats: &'a [GeneBinStats],
}

fn trajectory_block_weights(
    blocks: &FakeMatchingBlocks<'_>,
    source: usize,
    target: usize,
    lineage: usize,
) -> Vec<f64> {
    let mut source_weights = vec![0.0f64; blocks.n_blocks];
    let mut target_weights = vec![0.0f64; blocks.n_blocks];
    let n_cells = blocks
        .block_labels
        .len()
        .min(blocks.cell_types.len())
        .min(blocks.lineage_weights.len() / blocks.n_lineages.max(1));
    for cell in 0..n_cells {
        let ct = blocks.cell_types[cell];
        if ct != source && ct != target {
            continue;
        }
        let block = blocks.block_labels[cell];
        if block >= blocks.n_blocks {
            continue;
        }
        let w = blocks.lineage_weights[cell * blocks.n_lineages + lineage] as f64;
        if w >= 1e-8 {
            if ct == source {
                source_weights[block] += w;
            }
            if ct == target {
                target_weights[block] += w;
            }
        }
    }
    source_weights
        .into_iter()
        .zip(target_weights)
        .map(|(source_weight, target_weight)| source_weight.min(target_weight))
        .collect()
}

fn allocate_draws_across_blocks(total: usize, weights: &[f64]) -> Vec<usize> {
    if total == 0 {
        return vec![0; weights.len()];
    }
    let sum: f64 = weights.iter().sum();
    if sum <= 1e-8 {
        return vec![0; weights.len()];
    }
    let mut out = vec![0usize; weights.len()];
    let mut assigned = 0usize;
    let mut remainders: Vec<(usize, f64)> = weights
        .iter()
        .enumerate()
        .map(|(i, &w)| {
            let exact = (total as f64) * w / sum;
            let base = exact.floor() as usize;
            out[i] = base;
            assigned += base;
            (i, exact - base as f64)
        })
        .collect();
    remainders.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
    let mut left = total.saturating_sub(assigned);
    for (i, _) in remainders {
        if left == 0 {
            break;
        }
        out[i] += 1;
        left -= 1;
    }
    out
}

fn draw_one_complex_generic(
    rng: &mut StdRng,
    ligand_subunits: &[usize],
    receptor_subunits: &[usize],
    source: usize,
    target: usize,
    lineage: usize,
    real_stats: &GeneBinStats,
    fake_stats: &GeneBinStats,
    fake_gene_rows: &[usize],
    all_fake_idx: &[usize],
    tau_by_gene: &std::collections::HashMap<usize, f32>,
    n_strata: usize,
    use_strata: bool,
    caliper_k: usize,
    max_match_distance: f32,
    mins: &[f32],
    maxs: &[f32],
    stratum_pools: Option<&std::collections::HashMap<(usize, usize, usize), Vec<usize>>>,
) -> Option<(Vec<usize>, Vec<usize>, f32, f32)> {
    let mut lig_topks: Vec<Vec<(usize, f32)>> = Vec::with_capacity(ligand_subunits.len());
    for &sub in ligand_subunits {
        let tau = *tau_by_gene.get(&sub).unwrap_or(&0.0);
        let ctx = context_feature_vec4(real_stats, sub, source, lineage, tau);
        let topk = context_topk_pool(
            &ctx,
            all_fake_idx,
            fake_gene_rows,
            fake_stats,
            source,
            lineage,
            tau_by_gene,
            n_strata,
            use_strata,
            mins,
            maxs,
            caliper_k,
            max_match_distance,
            stratum_pools,
        );
        if topk.is_empty() {
            return None;
        }
        lig_topks.push(topk);
    }
    let mut rec_topks: Vec<Vec<(usize, f32)>> = Vec::with_capacity(receptor_subunits.len());
    for &sub in receptor_subunits {
        let tau = *tau_by_gene.get(&sub).unwrap_or(&0.0);
        let ctx = context_feature_vec4(real_stats, sub, target, lineage, tau);
        let topk = context_topk_pool(
            &ctx,
            all_fake_idx,
            fake_gene_rows,
            fake_stats,
            target,
            lineage,
            tau_by_gene,
            n_strata,
            use_strata,
            mins,
            maxs,
            caliper_k,
            max_match_distance,
            stratum_pools,
        );
        if topk.is_empty() {
            return None;
        }
        rec_topks.push(topk);
    }

    for _ in 0..50 {
        let mut used = HashSet::new();
        let mut fake_lig = Vec::with_capacity(ligand_subunits.len());
        let mut fake_rec = Vec::with_capacity(receptor_subunits.len());
        let mut lig_dist_sum = 0.0f32;
        let mut rec_dist_sum = 0.0f32;
        let mut ok = true;
        for topk in &lig_topks {
            let (idx, dist) = topk[rng.gen_range(0..topk.len())];
            let gene = fake_gene_rows[idx];
            if used.contains(&gene) {
                ok = false;
                break;
            }
            used.insert(gene);
            fake_lig.push(gene);
            lig_dist_sum += dist;
        }
        if !ok {
            continue;
        }
        for topk in &rec_topks {
            let (idx, dist) = topk[rng.gen_range(0..topk.len())];
            let gene = fake_gene_rows[idx];
            if used.contains(&gene) {
                ok = false;
                break;
            }
            used.insert(gene);
            fake_rec.push(gene);
            rec_dist_sum += dist;
        }
        if !ok {
            continue;
        }
        let lig_mean = lig_dist_sum / ligand_subunits.len().max(1) as f32;
        let rec_mean = rec_dist_sum / receptor_subunits.len().max(1) as f32;
        return Some((fake_lig, fake_rec, lig_mean, rec_mean));
    }
    None
}

fn draw_one_complex_side(
    rng: &mut StdRng,
    subunits: &[usize],
    ct: usize,
    lineage: usize,
    real_stats: &GeneBinStats,
    fake_stats: &GeneBinStats,
    fake_gene_rows: &[usize],
    all_fake_idx: &[usize],
    tau_by_gene: &std::collections::HashMap<usize, f32>,
    n_strata: usize,
    use_strata: bool,
    caliper_k: usize,
    max_match_distance: f32,
    mins: &[f32],
    maxs: &[f32],
    stratum_pools: Option<&std::collections::HashMap<(usize, usize, usize), Vec<usize>>>,
) -> Option<(Vec<usize>, f32)> {
    let mut topks: Vec<Vec<(usize, f32)>> = Vec::with_capacity(subunits.len());
    for &sub in subunits {
        let tau = *tau_by_gene.get(&sub).unwrap_or(&0.0);
        let ctx = context_feature_vec4(real_stats, sub, ct, lineage, tau);
        let topk = context_topk_pool(
            &ctx,
            all_fake_idx,
            fake_gene_rows,
            fake_stats,
            ct,
            lineage,
            tau_by_gene,
            n_strata,
            use_strata,
            mins,
            maxs,
            caliper_k,
            max_match_distance,
            stratum_pools,
        );
        if topk.is_empty() {
            return None;
        }
        topks.push(topk);
    }
    for _ in 0..50 {
        let mut used = HashSet::new();
        let mut fake_subs = Vec::with_capacity(subunits.len());
        let mut dist_sum = 0.0f32;
        let mut ok = true;
        for topk in &topks {
            let (idx, dist) = topk[rng.gen_range(0..topk.len())];
            let gene = fake_gene_rows[idx];
            if used.contains(&gene) {
                ok = false;
                break;
            }
            used.insert(gene);
            fake_subs.push(gene);
            dist_sum += dist;
        }
        if ok {
            return Some((fake_subs, dist_sum / subunits.len().max(1) as f32));
        }
    }
    None
}

fn draw_complex_edge_for_key(
    out: &mut Vec<FakeLrPair>,
    rng: &mut StdRng,
    pair_idx: usize,
    source: usize,
    target: usize,
    lineage: usize,
    pair: &crate::ccc::LrPair,
    real_stats: &GeneBinStats,
    fake_stats: &GeneBinStats,
    fake_gene_rows: &[usize],
    all_fake_idx: &[usize],
    tau_by_gene: &std::collections::HashMap<usize, f32>,
    n_edge_budget: usize,
    n_strata: usize,
    use_strata: bool,
    caliper_k: usize,
    max_match_distance: f32,
    mins: &[f32],
    maxs: &[f32],
    stratum_pools: Option<&std::collections::HashMap<(usize, usize, usize), Vec<usize>>>,
) {
    if n_edge_budget == 0 {
        return;
    }
    let (n_lig, n_rec) = split_edge_budget(n_edge_budget);
    let real_rec = pair.receptor_rows.clone();
    let real_lig = pair.ligand_rows.clone();
    for _ in 0..n_lig {
        let Some((fake_lig, lig_dist)) = draw_one_complex_side(
            rng,
            &pair.ligand_rows,
            source,
            lineage,
            real_stats,
            fake_stats,
            fake_gene_rows,
            all_fake_idx,
            tau_by_gene,
            n_strata,
            use_strata,
            caliper_k,
            max_match_distance,
            mins,
            maxs,
            stratum_pools,
        ) else {
            continue;
        };
        let lig_export =
            bottleneck_gene_in_context(real_stats, &fake_lig, source, lineage, tau_by_gene);
        let rec_export =
            bottleneck_gene_in_context(real_stats, &real_rec, target, lineage, tau_by_gene);
        out.push(FakeLrPair {
            ligand_gene: lig_export,
            receptor_gene: rec_export,
            complex_ligand_rows: Some(fake_lig),
            complex_receptor_rows: Some(real_rec.clone()),
            shadowing_pair_idx: pair_idx,
            source,
            target,
            lineage,
            null_kind: FAKE_NULL_REWIRE_LIGAND,
            ligand_match_distance: lig_dist,
            receptor_match_distance: 0.0,
        });
    }
    for _ in 0..n_rec {
        let Some((fake_rec, rec_dist)) = draw_one_complex_side(
            rng,
            &pair.receptor_rows,
            target,
            lineage,
            real_stats,
            fake_stats,
            fake_gene_rows,
            all_fake_idx,
            tau_by_gene,
            n_strata,
            use_strata,
            caliper_k,
            max_match_distance,
            mins,
            maxs,
            stratum_pools,
        ) else {
            continue;
        };
        let lig_export =
            bottleneck_gene_in_context(real_stats, &real_lig, source, lineage, tau_by_gene);
        let rec_export =
            bottleneck_gene_in_context(real_stats, &fake_rec, target, lineage, tau_by_gene);
        out.push(FakeLrPair {
            ligand_gene: lig_export,
            receptor_gene: rec_export,
            complex_ligand_rows: Some(real_lig.clone()),
            complex_receptor_rows: Some(fake_rec),
            shadowing_pair_idx: pair_idx,
            source,
            target,
            lineage,
            null_kind: FAKE_NULL_REWIRE_RECEPTOR,
            ligand_match_distance: 0.0,
            receptor_match_distance: rec_dist,
        });
    }
}

fn draw_complex_network_for_key(
    out: &mut Vec<FakeLrPair>,
    rng: &mut StdRng,
    pair_idx: usize,
    source: usize,
    target: usize,
    lineage: usize,
    pair: &crate::ccc::LrPair,
    stats: &GeneBinStats,
    _ligand_pool: &[usize],
    _receptor_pool: &[usize],
    _ligand_pool_idx: &[usize],
    _receptor_pool_idx: &[usize],
    fake_gene_rows: &[usize],
    all_fake_idx: &[usize],
    real_lr_pairs: &std::collections::HashSet<(usize, usize)>,
    tau_by_gene: &std::collections::HashMap<usize, f32>,
    n_network_budget: usize,
    n_strata: usize,
    use_strata_lig: bool,
    use_strata_rec: bool,
    caliper_k: usize,
    max_match_distance: f32,
    mins: &[f32],
    maxs: &[f32],
    lig_stratum_pools: Option<&std::collections::HashMap<(usize, usize, usize), Vec<usize>>>,
    rec_stratum_pools: Option<&std::collections::HashMap<(usize, usize, usize), Vec<usize>>>,
) {
    if n_network_budget == 0 {
        return;
    }
    let real_lig = pair.ligand_rows.clone();
    let real_rec = pair.receptor_rows.clone();
    let (n_lig, n_rec) = split_edge_budget(n_network_budget);
    for _ in 0..n_lig {
        let Some((fake_lig, lig_dist)) = draw_one_complex_side(
            rng,
            &pair.ligand_rows,
            source,
            lineage,
            stats,
            stats,
            fake_gene_rows,
            all_fake_idx,
            tau_by_gene,
            n_strata,
            use_strata_lig,
            caliper_k,
            max_match_distance,
            mins,
            maxs,
            lig_stratum_pools,
        ) else {
            continue;
        };
        let lig_export = bottleneck_gene_in_context(stats, &fake_lig, source, lineage, tau_by_gene);
        let rec_export = bottleneck_gene_in_context(stats, &real_rec, target, lineage, tau_by_gene);
        if !real_lr_pairs.contains(&(lig_export, rec_export)) {
            out.push(FakeLrPair {
                ligand_gene: lig_export,
                receptor_gene: rec_export,
                complex_ligand_rows: Some(fake_lig),
                complex_receptor_rows: Some(real_rec.clone()),
                shadowing_pair_idx: pair_idx,
                source,
                target,
                lineage,
                null_kind: FAKE_NULL_NETWORK_REWIRE_LIGAND,
                ligand_match_distance: lig_dist,
                receptor_match_distance: 0.0,
            });
        }
    }
    for _ in 0..n_rec {
        let Some((fake_rec, rec_dist)) = draw_one_complex_side(
            rng,
            &pair.receptor_rows,
            target,
            lineage,
            stats,
            stats,
            fake_gene_rows,
            all_fake_idx,
            tau_by_gene,
            n_strata,
            use_strata_rec,
            caliper_k,
            max_match_distance,
            mins,
            maxs,
            rec_stratum_pools,
        ) else {
            continue;
        };
        let lig_export = bottleneck_gene_in_context(stats, &real_lig, source, lineage, tau_by_gene);
        let rec_export = bottleneck_gene_in_context(stats, &fake_rec, target, lineage, tau_by_gene);
        if !real_lr_pairs.contains(&(lig_export, rec_export)) {
            out.push(FakeLrPair {
                ligand_gene: lig_export,
                receptor_gene: rec_export,
                complex_ligand_rows: Some(real_lig.clone()),
                complex_receptor_rows: Some(fake_rec),
                shadowing_pair_idx: pair_idx,
                source,
                target,
                lineage,
                null_kind: FAKE_NULL_NETWORK_REWIRE_RECEPTOR,
                ligand_match_distance: 0.0,
                receptor_match_distance: rec_dist,
            });
        }
    }
}

fn draw_complex_null_layers_for_key(
    out: &mut Vec<FakeLrPair>,
    rng: &mut StdRng,
    pair_idx: usize,
    source: usize,
    target: usize,
    lineage: usize,
    pair: &crate::ccc::LrPair,
    real_stats: &GeneBinStats,
    fake_stats: &GeneBinStats,
    fake_gene_rows: &[usize],
    all_fake_idx: &[usize],
    tau_by_gene: &std::collections::HashMap<usize, f32>,
    n_generic_budget: usize,
    n_edge_budget: usize,
    n_network_budget: usize,
    n_strata: usize,
    use_strata: bool,
    network_use_strata_lig: bool,
    network_use_strata_rec: bool,
    caliper_k: usize,
    max_match_distance: f32,
    mins: &[f32],
    maxs: &[f32],
    stratum_pools: Option<&std::collections::HashMap<(usize, usize, usize), Vec<usize>>>,
    ligand_pool: &[usize],
    receptor_pool: &[usize],
    ligand_pool_idx: &[usize],
    receptor_pool_idx: &[usize],
    real_lr_pairs: &std::collections::HashSet<(usize, usize)>,
    lig_stratum_pools: Option<&std::collections::HashMap<(usize, usize, usize), Vec<usize>>>,
    rec_stratum_pools: Option<&std::collections::HashMap<(usize, usize, usize), Vec<usize>>>,
) {
    draw_complex_generic_for_key(
        out,
        rng,
        pair_idx,
        source,
        target,
        lineage,
        pair,
        real_stats,
        fake_stats,
        fake_gene_rows,
        all_fake_idx,
        tau_by_gene,
        n_generic_budget,
        n_strata,
        use_strata,
        caliper_k,
        max_match_distance,
        mins,
        maxs,
        stratum_pools,
    );
    draw_complex_edge_for_key(
        out,
        rng,
        pair_idx,
        source,
        target,
        lineage,
        pair,
        real_stats,
        fake_stats,
        fake_gene_rows,
        all_fake_idx,
        tau_by_gene,
        n_edge_budget,
        n_strata,
        use_strata,
        caliper_k,
        max_match_distance,
        mins,
        maxs,
        stratum_pools,
    );
    draw_complex_network_for_key(
        out,
        rng,
        pair_idx,
        source,
        target,
        lineage,
        pair,
        real_stats,
        ligand_pool,
        receptor_pool,
        ligand_pool_idx,
        receptor_pool_idx,
        fake_gene_rows,
        all_fake_idx,
        real_lr_pairs,
        tau_by_gene,
        n_network_budget,
        n_strata,
        network_use_strata_lig,
        network_use_strata_rec,
        caliper_k,
        max_match_distance,
        mins,
        maxs,
        lig_stratum_pools,
        rec_stratum_pools,
    );
}

fn draw_complex_generic_for_key(
    out: &mut Vec<FakeLrPair>,
    rng: &mut StdRng,
    pair_idx: usize,
    source: usize,
    target: usize,
    lineage: usize,
    pair: &crate::ccc::LrPair,
    real_stats: &GeneBinStats,
    fake_stats: &GeneBinStats,
    fake_gene_rows: &[usize],
    all_fake_idx: &[usize],
    tau_by_gene: &std::collections::HashMap<usize, f32>,
    n_generic_budget: usize,
    n_strata: usize,
    use_strata: bool,
    caliper_k: usize,
    max_match_distance: f32,
    mins: &[f32],
    maxs: &[f32],
    stratum_pools: Option<&std::collections::HashMap<(usize, usize, usize), Vec<usize>>>,
) {
    if n_generic_budget == 0 {
        return;
    }
    for _ in 0..n_generic_budget {
        let Some((fake_lig, fake_rec, lig_dist, rec_dist)) = draw_one_complex_generic(
            rng,
            &pair.ligand_rows,
            &pair.receptor_rows,
            source,
            target,
            lineage,
            real_stats,
            fake_stats,
            fake_gene_rows,
            all_fake_idx,
            tau_by_gene,
            n_strata,
            use_strata,
            caliper_k,
            max_match_distance,
            mins,
            maxs,
            stratum_pools,
        ) else {
            continue;
        };
        let lig_export =
            bottleneck_gene_in_context(real_stats, &fake_lig, source, lineage, tau_by_gene);
        let rec_export =
            bottleneck_gene_in_context(real_stats, &fake_rec, target, lineage, tau_by_gene);
        out.push(FakeLrPair {
            ligand_gene: lig_export,
            receptor_gene: rec_export,
            complex_ligand_rows: Some(fake_lig),
            complex_receptor_rows: Some(fake_rec),
            shadowing_pair_idx: pair_idx,
            source,
            target,
            lineage,
            null_kind: FAKE_NULL_GENERIC,
            ligand_match_distance: lig_dist,
            receptor_match_distance: rec_dist,
        });
    }
}

fn draw_fake_pairs_for_key(
    out: &mut Vec<FakeLrPair>,
    rng: &mut StdRng,
    pair_idx: usize,
    source: usize,
    target: usize,
    lineage: usize,
    real_lig: usize,
    real_rec: usize,
    real_stats: &GeneBinStats,
    fake_stats: &GeneBinStats,
    fake_gene_rows: &[usize],
    all_fake_idx: &[usize],
    tau_by_gene: &std::collections::HashMap<usize, f32>,
    n_generic_budget: usize,
    n_edge_budget: usize,
    n_strata: usize,
    use_strata: bool,
    caliper_k: usize,
    max_match_distance: f32,
    mins: &[f32],
    maxs: &[f32],
    stratum_pools: Option<&std::collections::HashMap<(usize, usize, usize), Vec<usize>>>,
) {
    if n_generic_budget == 0 && n_edge_budget == 0 {
        return;
    }
    let lig_tau = *tau_by_gene.get(&real_lig).unwrap_or(&0.0);
    let rec_tau = *tau_by_gene.get(&real_rec).unwrap_or(&0.0);
    let lig_ctx = context_feature_vec4(real_stats, real_lig, source, lineage, lig_tau);
    let rec_ctx = context_feature_vec4(real_stats, real_rec, target, lineage, rec_tau);
    let lig_topk = context_topk_pool(
        &lig_ctx,
        all_fake_idx,
        fake_gene_rows,
        fake_stats,
        source,
        lineage,
        tau_by_gene,
        n_strata,
        use_strata,
        &mins,
        &maxs,
        caliper_k,
        max_match_distance,
        stratum_pools,
    );
    let rec_topk = context_topk_pool(
        &rec_ctx,
        all_fake_idx,
        fake_gene_rows,
        fake_stats,
        target,
        lineage,
        tau_by_gene,
        n_strata,
        use_strata,
        &mins,
        &maxs,
        caliper_k,
        max_match_distance,
        stratum_pools,
    );
    for _ in 0..n_generic_budget {
        if lig_topk.is_empty() || rec_topk.is_empty() {
            break;
        }
        for _ in 0..50 {
            let (li, ld) = lig_topk[rng.gen_range(0..lig_topk.len())];
            let (ri, rd) = rec_topk[rng.gen_range(0..rec_topk.len())];
            let lg = fake_gene_rows[li];
            let rg = fake_gene_rows[ri];
            if lg != rg {
                out.push(FakeLrPair {
                    ligand_gene: lg,
                    receptor_gene: rg,
                    complex_ligand_rows: None,
                    complex_receptor_rows: None,
                    shadowing_pair_idx: pair_idx,
                    source,
                    target,
                    lineage,
                    null_kind: FAKE_NULL_GENERIC,
                    ligand_match_distance: ld,
                    receptor_match_distance: rd,
                });
                break;
            }
        }
    }
    for j in 0..n_edge_budget {
        if j % 2 == 0 {
            if !rec_topk.is_empty() {
                for _ in 0..50 {
                    let (ri, rd) = rec_topk[rng.gen_range(0..rec_topk.len())];
                    let rg = fake_gene_rows[ri];
                    if real_lig != rg {
                        out.push(FakeLrPair {
                            ligand_gene: real_lig,
                            receptor_gene: rg,
                            complex_ligand_rows: None,
                            complex_receptor_rows: None,
                            shadowing_pair_idx: pair_idx,
                            source,
                            target,
                            lineage,
                            null_kind: FAKE_NULL_REWIRE_RECEPTOR,
                            ligand_match_distance: 0.0,
                            receptor_match_distance: rd,
                        });
                        break;
                    }
                }
            }
        } else if !lig_topk.is_empty() {
            for _ in 0..50 {
                let (li, ld) = lig_topk[rng.gen_range(0..lig_topk.len())];
                let lg = fake_gene_rows[li];
                if lg != real_rec {
                    out.push(FakeLrPair {
                        ligand_gene: lg,
                        receptor_gene: real_rec,
                        complex_ligand_rows: None,
                        complex_receptor_rows: None,
                        shadowing_pair_idx: pair_idx,
                        source,
                        target,
                        lineage,
                        null_kind: FAKE_NULL_REWIRE_LIGAND,
                        ligand_match_distance: ld,
                        receptor_match_distance: 0.0,
                    });
                    break;
                }
            }
        }
    }
}

fn draw_network_rewire_for_key(
    out: &mut Vec<FakeLrPair>,
    rng: &mut StdRng,
    pair_idx: usize,
    source: usize,
    target: usize,
    lineage: usize,
    real_lig: usize,
    real_rec: usize,
    stats: &GeneBinStats,
    ligand_pool: &[usize],
    receptor_pool: &[usize],
    ligand_pool_idx: &[usize],
    receptor_pool_idx: &[usize],
    real_lr_pairs: &std::collections::HashSet<(usize, usize)>,
    tau_by_gene: &std::collections::HashMap<usize, f32>,
    n_network_budget: usize,
    n_strata: usize,
    use_strata: bool,
    caliper_k: usize,
    max_match_distance: f32,
    mins: &[f32],
    maxs: &[f32],
    lig_stratum_pools: Option<&std::collections::HashMap<(usize, usize, usize), Vec<usize>>>,
    rec_stratum_pools: Option<&std::collections::HashMap<(usize, usize, usize), Vec<usize>>>,
) {
    if n_network_budget == 0 || ligand_pool.is_empty() || receptor_pool.is_empty() {
        return;
    }
    let lig_tau = *tau_by_gene.get(&real_lig).unwrap_or(&0.0);
    let rec_tau = *tau_by_gene.get(&real_rec).unwrap_or(&0.0);
    let lig_ctx = context_feature_vec4(stats, real_lig, source, lineage, lig_tau);
    let rec_ctx = context_feature_vec4(stats, real_rec, target, lineage, rec_tau);
    let lig_topk = context_topk_pool(
        &lig_ctx,
        ligand_pool_idx,
        ligand_pool,
        stats,
        source,
        lineage,
        tau_by_gene,
        n_strata,
        use_strata,
        mins,
        maxs,
        caliper_k,
        max_match_distance,
        lig_stratum_pools,
    );
    let rec_topk = context_topk_pool(
        &rec_ctx,
        receptor_pool_idx,
        receptor_pool,
        stats,
        target,
        lineage,
        tau_by_gene,
        n_strata,
        use_strata,
        mins,
        maxs,
        caliper_k,
        max_match_distance,
        rec_stratum_pools,
    );
    for j in 0..n_network_budget {
        if j % 2 == 0 {
            if !rec_topk.is_empty() {
                for _ in 0..50 {
                    let (ri, rd) = rec_topk[rng.gen_range(0..rec_topk.len())];
                    let rg = receptor_pool[ri];
                    if real_lig != rg && !real_lr_pairs.contains(&(real_lig, rg)) {
                        out.push(FakeLrPair {
                            ligand_gene: real_lig,
                            receptor_gene: rg,
                            complex_ligand_rows: None,
                            complex_receptor_rows: None,
                            shadowing_pair_idx: pair_idx,
                            source,
                            target,
                            lineage,
                            null_kind: FAKE_NULL_NETWORK_REWIRE_RECEPTOR,
                            ligand_match_distance: 0.0,
                            receptor_match_distance: rd,
                        });
                        break;
                    }
                }
            }
        } else if !lig_topk.is_empty() {
            for _ in 0..50 {
                let (li, ld) = lig_topk[rng.gen_range(0..lig_topk.len())];
                let lg = ligand_pool[li];
                if lg != real_rec && lg != real_lig && !real_lr_pairs.contains(&(lg, real_rec)) {
                    out.push(FakeLrPair {
                        ligand_gene: lg,
                        receptor_gene: real_rec,
                        complex_ligand_rows: None,
                        complex_receptor_rows: None,
                        shadowing_pair_idx: pair_idx,
                        source,
                        target,
                        lineage,
                        null_kind: FAKE_NULL_NETWORK_REWIRE_LIGAND,
                        ligand_match_distance: ld,
                        receptor_match_distance: 0.0,
                    });
                    break;
                }
            }
        }
    }
}

/// Build trajectory-specific fake draws with source/target/lineage context matching.
pub fn construct_fake_pairs_context_aware(
    keys: &[crate::ccc::trajectory_curve::TrajKey],
    lr_pairs: &[crate::ccc::LrPair],
    real_stats: &GeneBinStats,
    fake_stats: &GeneBinStats,
    fake_gene_rows: &[usize],
    tau_by_gene: &std::collections::HashMap<usize, f32>,
    n_fake_per_key: usize,
    n_strata: usize,
    caliper_k: usize,
    max_match_distance: f32,
    mode: FakeNullMode,
    complex_null_mode: ComplexNullMode,
    seed: u64,
    block_matching: Option<&FakeMatchingBlocks<'_>>,
) -> Vec<FakeLrPair> {
    if keys.is_empty() || n_fake_per_key == 0 {
        return Vec::new();
    }
    if fake_gene_rows.is_empty()
        && mode != FakeNullMode::LrNetworkRewire
        && mode != FakeNullMode::All
    {
        return Vec::new();
    }
    let mut rng = StdRng::seed_from_u64(seed);
    let all_fake_idx: Vec<usize> = (0..fake_gene_rows.len()).collect();
    let use_strata = n_strata > 0 && fake_gene_rows.len() >= 100;
    let use_blocks = block_matching.map(|b| b.n_blocks > 1).unwrap_or(false);
    let (ligand_pool, receptor_pool, real_lr_pairs) = build_lr_network_pools(lr_pairs);
    let ligand_pool_idx: Vec<usize> = (0..ligand_pool.len()).collect();
    let receptor_pool_idx: Vec<usize> = (0..receptor_pool.len()).collect();
    let network_use_strata_lig = n_strata > 0 && ligand_pool.len() >= 100;
    let network_use_strata_rec = n_strata > 0 && receptor_pool.len() >= 100;
    let mut context_pairs: std::collections::HashSet<(usize, usize)> =
        std::collections::HashSet::new();
    for &(_pair_idx, source, target, lineage) in keys {
        context_pairs.insert((source, lineage));
        context_pairs.insert((target, lineage));
    }
    let context_list: Vec<(usize, usize)> = context_pairs.into_iter().collect();
    let global_ranges = if !use_blocks {
        Some(collect_context_ranges(
            real_stats,
            fake_stats,
            fake_gene_rows,
            keys,
            lr_pairs,
            &context_list,
            tau_by_gene,
        ))
    } else {
        None
    };
    let block_ranges: Vec<(Vec<f32>, Vec<f32>)> = if use_blocks {
        if let Some(bp) = block_matching {
            bp.real_stats
                .iter()
                .zip(bp.fake_stats.iter())
                .map(|(rs, fs)| {
                    collect_context_ranges(
                        rs,
                        fs,
                        fake_gene_rows,
                        keys,
                        lr_pairs,
                        &context_list,
                        tau_by_gene,
                    )
                })
                .collect()
        } else {
            Vec::new()
        }
    } else {
        Vec::new()
    };
    let global_stratum_pools = if !use_blocks && use_strata {
        let (ref mins, ref maxs) = global_ranges.as_ref().unwrap();
        Some(build_fake_stratum_pools(
            &all_fake_idx,
            fake_gene_rows,
            fake_stats,
            tau_by_gene,
            n_strata,
            &mins,
            &maxs,
            &context_list,
        ))
    } else {
        None
    };
    let block_stratum_pools: Vec<std::collections::HashMap<(usize, usize, usize), Vec<usize>>> =
        if use_blocks && use_strata {
            if let Some(bp) = block_matching {
                bp.real_stats
                    .iter()
                    .zip(bp.fake_stats.iter())
                    .zip(block_ranges.iter())
                    .map(|((_rs, fs), (mins, maxs))| {
                        build_fake_stratum_pools(
                            &all_fake_idx,
                            fake_gene_rows,
                            fs,
                            tau_by_gene,
                            n_strata,
                            mins,
                            maxs,
                            &context_list,
                        )
                    })
                    .collect()
            } else {
                Vec::new()
            }
        } else {
            Vec::new()
        };
    let global_lig_stratum_pools = if !use_blocks && network_use_strata_lig {
        let (ref mins, ref maxs) = global_ranges.as_ref().unwrap();
        Some(build_fake_stratum_pools(
            &ligand_pool_idx,
            &ligand_pool,
            real_stats,
            tau_by_gene,
            n_strata,
            mins,
            maxs,
            &context_list,
        ))
    } else {
        None
    };
    let global_rec_stratum_pools = if !use_blocks && network_use_strata_rec {
        let (ref mins, ref maxs) = global_ranges.as_ref().unwrap();
        Some(build_fake_stratum_pools(
            &receptor_pool_idx,
            &receptor_pool,
            real_stats,
            tau_by_gene,
            n_strata,
            mins,
            maxs,
            &context_list,
        ))
    } else {
        None
    };
    let block_lig_stratum_pools: Vec<std::collections::HashMap<(usize, usize, usize), Vec<usize>>> =
        if use_blocks && network_use_strata_lig {
            if let Some(bp) = block_matching {
                bp.real_stats
                    .iter()
                    .zip(block_ranges.iter())
                    .map(|(rs, (mins, maxs))| {
                        build_fake_stratum_pools(
                            &ligand_pool_idx,
                            &ligand_pool,
                            rs,
                            tau_by_gene,
                            n_strata,
                            mins,
                            maxs,
                            &context_list,
                        )
                    })
                    .collect()
            } else {
                Vec::new()
            }
        } else {
            Vec::new()
        };
    let block_rec_stratum_pools: Vec<std::collections::HashMap<(usize, usize, usize), Vec<usize>>> =
        if use_blocks && network_use_strata_rec {
            if let Some(bp) = block_matching {
                bp.real_stats
                    .iter()
                    .zip(block_ranges.iter())
                    .map(|(rs, (mins, maxs))| {
                        build_fake_stratum_pools(
                            &receptor_pool_idx,
                            &receptor_pool,
                            rs,
                            tau_by_gene,
                            n_strata,
                            mins,
                            maxs,
                            &context_list,
                        )
                    })
                    .collect()
            } else {
                Vec::new()
            }
        } else {
            Vec::new()
        };
    let mut out = Vec::with_capacity(keys.len() * n_fake_per_key);

    for &(pair_idx, source, target, lineage) in keys {
        let pair = &lr_pairs[pair_idx];
        if pair.ligand_rows.is_empty() || pair.receptor_rows.is_empty() {
            continue;
        }
        if is_complex_lr_pair(pair) {
            if complex_null_mode != ComplexNullMode::CardinalityMatched
                || !complex_cardinality_supported(pair)
            {
                continue;
            }
            let (n_generic_total, n_edge_total, n_network_total) =
                split_null_budget_layers(n_fake_per_key, mode);
            if n_generic_total + n_edge_total + n_network_total == 0 {
                continue;
            }
            if use_blocks {
                let bp = block_matching.unwrap();
                let weights = trajectory_block_weights(bp, source, target, lineage);
                let generic_budgets = allocate_draws_across_blocks(n_generic_total, &weights);
                let edge_budgets = allocate_draws_across_blocks(n_edge_total, &weights);
                let network_budgets = allocate_draws_across_blocks(n_network_total, &weights);
                for block_id in 0..bp.n_blocks {
                    if generic_budgets[block_id] == 0
                        && edge_budgets[block_id] == 0
                        && network_budgets[block_id] == 0
                    {
                        continue;
                    }
                    draw_complex_null_layers_for_key(
                        &mut out,
                        &mut rng,
                        pair_idx,
                        source,
                        target,
                        lineage,
                        pair,
                        &bp.real_stats[block_id],
                        &bp.fake_stats[block_id],
                        fake_gene_rows,
                        &all_fake_idx,
                        tau_by_gene,
                        generic_budgets[block_id],
                        edge_budgets[block_id],
                        network_budgets[block_id],
                        n_strata,
                        use_strata,
                        network_use_strata_lig,
                        network_use_strata_rec,
                        caliper_k,
                        max_match_distance,
                        &block_ranges[block_id].0,
                        &block_ranges[block_id].1,
                        if use_strata {
                            Some(&block_stratum_pools[block_id])
                        } else {
                            None
                        },
                        &ligand_pool,
                        &receptor_pool,
                        &ligand_pool_idx,
                        &receptor_pool_idx,
                        &real_lr_pairs,
                        if network_use_strata_lig {
                            Some(&block_lig_stratum_pools[block_id])
                        } else {
                            None
                        },
                        if network_use_strata_rec {
                            Some(&block_rec_stratum_pools[block_id])
                        } else {
                            None
                        },
                    );
                }
            } else {
                let (ref mins, ref maxs) = global_ranges.as_ref().unwrap();
                draw_complex_null_layers_for_key(
                    &mut out,
                    &mut rng,
                    pair_idx,
                    source,
                    target,
                    lineage,
                    pair,
                    real_stats,
                    fake_stats,
                    fake_gene_rows,
                    &all_fake_idx,
                    tau_by_gene,
                    n_generic_total,
                    n_edge_total,
                    n_network_total,
                    n_strata,
                    use_strata,
                    network_use_strata_lig,
                    network_use_strata_rec,
                    caliper_k,
                    max_match_distance,
                    mins,
                    maxs,
                    global_stratum_pools.as_ref(),
                    &ligand_pool,
                    &receptor_pool,
                    &ligand_pool_idx,
                    &receptor_pool_idx,
                    &real_lr_pairs,
                    global_lig_stratum_pools.as_ref(),
                    global_rec_stratum_pools.as_ref(),
                );
            }
            continue;
        }
        let real_lig =
            bottleneck_gene_in_context(real_stats, &pair.ligand_rows, source, lineage, tau_by_gene);
        let real_rec = bottleneck_gene_in_context(
            real_stats,
            &pair.receptor_rows,
            target,
            lineage,
            tau_by_gene,
        );
        let (n_generic_total, n_edge_total, n_network_total) =
            split_null_budget_layers(n_fake_per_key, mode);
        if use_blocks {
            let bp = block_matching.unwrap();
            let weights = trajectory_block_weights(bp, source, target, lineage);
            let generic_budgets = allocate_draws_across_blocks(n_generic_total, &weights);
            let edge_budgets = allocate_draws_across_blocks(n_edge_total, &weights);
            let network_budgets = allocate_draws_across_blocks(n_network_total, &weights);
            for block_id in 0..bp.n_blocks {
                if generic_budgets[block_id] == 0
                    && edge_budgets[block_id] == 0
                    && network_budgets[block_id] == 0
                {
                    continue;
                }
                draw_fake_pairs_for_key(
                    &mut out,
                    &mut rng,
                    pair_idx,
                    source,
                    target,
                    lineage,
                    real_lig,
                    real_rec,
                    &bp.real_stats[block_id],
                    &bp.fake_stats[block_id],
                    fake_gene_rows,
                    &all_fake_idx,
                    tau_by_gene,
                    generic_budgets[block_id],
                    edge_budgets[block_id],
                    n_strata,
                    use_strata,
                    caliper_k,
                    max_match_distance,
                    &block_ranges[block_id].0,
                    &block_ranges[block_id].1,
                    if use_strata {
                        Some(&block_stratum_pools[block_id])
                    } else {
                        None
                    },
                );
                draw_network_rewire_for_key(
                    &mut out,
                    &mut rng,
                    pair_idx,
                    source,
                    target,
                    lineage,
                    real_lig,
                    real_rec,
                    &bp.real_stats[block_id],
                    &ligand_pool,
                    &receptor_pool,
                    &ligand_pool_idx,
                    &receptor_pool_idx,
                    &real_lr_pairs,
                    tau_by_gene,
                    network_budgets[block_id],
                    n_strata,
                    network_use_strata_lig || network_use_strata_rec,
                    caliper_k,
                    max_match_distance,
                    &block_ranges[block_id].0,
                    &block_ranges[block_id].1,
                    if network_use_strata_lig {
                        Some(&block_lig_stratum_pools[block_id])
                    } else {
                        None
                    },
                    if network_use_strata_rec {
                        Some(&block_rec_stratum_pools[block_id])
                    } else {
                        None
                    },
                );
            }
        } else {
            let (ref mins, ref maxs) = global_ranges.as_ref().unwrap();
            draw_fake_pairs_for_key(
                &mut out,
                &mut rng,
                pair_idx,
                source,
                target,
                lineage,
                real_lig,
                real_rec,
                real_stats,
                fake_stats,
                fake_gene_rows,
                &all_fake_idx,
                tau_by_gene,
                n_generic_total,
                n_edge_total,
                n_strata,
                use_strata,
                caliper_k,
                max_match_distance,
                mins,
                maxs,
                global_stratum_pools.as_ref(),
            );
            draw_network_rewire_for_key(
                &mut out,
                &mut rng,
                pair_idx,
                source,
                target,
                lineage,
                real_lig,
                real_rec,
                real_stats,
                &ligand_pool,
                &receptor_pool,
                &ligand_pool_idx,
                &receptor_pool_idx,
                &real_lr_pairs,
                tau_by_gene,
                n_network_total,
                n_strata,
                network_use_strata_lig || network_use_strata_rec,
                caliper_k,
                max_match_distance,
                mins,
                maxs,
                global_lig_stratum_pools.as_ref(),
                global_rec_stratum_pools.as_ref(),
            );
        }
    }
    out
}

fn build_fake_stratum_pools(
    all_fake_idx: &[usize],
    fake_gene_rows: &[usize],
    fake_stats: &GeneBinStats,
    tau_by_gene: &std::collections::HashMap<usize, f32>,
    n_strata: usize,
    mins: &[f32],
    maxs: &[f32],
    contexts: &[(usize, usize)],
) -> std::collections::HashMap<(usize, usize, usize), Vec<usize>> {
    let mut pools: std::collections::HashMap<(usize, usize, usize), Vec<usize>> =
        std::collections::HashMap::new();
    for &(ct, lin) in contexts {
        for &idx in all_fake_idx {
            let gene = fake_gene_rows[idx];
            let tau = *tau_by_gene.get(&gene).unwrap_or(&0.0);
            let cv = context_feature_vec4(fake_stats, gene, ct, lin, tau);
            let sid = feature_to_stratum_from_vec(&cv, n_strata, mins, maxs);
            pools.entry((ct, lin, sid)).or_default().push(idx);
        }
    }
    pools
}

fn context_topk_pool(
    ctx: &[f32; 4],
    all_fake_idx: &[usize],
    fake_gene_rows: &[usize],
    fake_stats: &GeneBinStats,
    ct: usize,
    lin: usize,
    tau_by_gene: &std::collections::HashMap<usize, f32>,
    n_strata: usize,
    use_strata: bool,
    mins: &[f32],
    maxs: &[f32],
    caliper_k: usize,
    max_match_distance: f32,
    stratum_pools: Option<&std::collections::HashMap<(usize, usize, usize), Vec<usize>>>,
) -> Vec<(usize, f32)> {
    if use_strata {
        let sid = feature_to_stratum_from_vec(ctx, n_strata, mins, maxs);
        let pool: Vec<usize> = if let Some(pools) = stratum_pools {
            pools.get(&(ct, lin, sid)).cloned().unwrap_or_default()
        } else {
            all_fake_idx
                .iter()
                .copied()
                .filter(|&idx| {
                    let g = fake_gene_rows[idx];
                    let tau = *tau_by_gene.get(&g).unwrap_or(&0.0);
                    let cv = context_feature_vec4(fake_stats, g, ct, lin, tau);
                    feature_to_stratum_from_vec(&cv, n_strata, mins, maxs) == sid
                })
                .collect()
        };
        let pool_ref: &[usize] = if pool.is_empty() { all_fake_idx } else { &pool };
        top_k_nearest_context(
            ctx,
            pool_ref,
            fake_gene_rows,
            fake_stats,
            ct,
            lin,
            tau_by_gene,
            caliper_k.min(pool_ref.len().max(1)),
            mins,
            maxs,
            max_match_distance,
        )
    } else {
        top_k_nearest_context(
            ctx,
            all_fake_idx,
            fake_gene_rows,
            fake_stats,
            ct,
            lin,
            tau_by_gene,
            caliper_k.min(all_fake_idx.len()),
            mins,
            maxs,
            max_match_distance,
        )
    }
}

/// Construct fake ligand-receptor pairs via stratified sampling +
/// within-stratum caliper matching.
///
/// For each real LR pair `(L_i, R_i)`:
/// 1. Determine L_i's 4-D feature stratum.
/// 2. Select top `caliper_k` candidate fake genes in the same stratum
///    (closest by Euclidean distance in 4-D normalised feature space).
/// 3. Same for R_i.
/// 4. Sample one fake L' and one fake R' uniformly with replacement from
///    their respective shortlists, rejecting self-pairs.
///
/// Sampling with replacement keeps every draw inside the matched calipers.
/// If a shortlist pair contains no valid non-self combination, that real pair
/// receives fewer draws and downstream coverage reports the deficit.
///
/// `features_grid` = list of (stratum_id, gene_features) for all
/// candidate fake genes, pre-computed by the caller.
pub fn construct_fake_pairs_stratified(
    real_lig_features: &[GeneFeatures],
    real_rec_features: &[GeneFeatures],
    n_pairs: usize,
    candidate_features: &[GeneFeatures],
    n_fake_per_real: usize,
    n_strata: usize, // per-axis (4 axes → n_strata^4 cells)
    caliper_k: usize,
    max_match_distance: f32,
    seed: u64,
) -> Vec<FakeLrPair> {
    if candidate_features.is_empty() || n_pairs == 0 {
        return Vec::new();
    }

    let mut rng = StdRng::seed_from_u64(seed);

    // ── Compute global feature ranges for normalisation ────────────────
    let normalization_features: Vec<GeneFeatures> = candidate_features
        .iter()
        .chain(real_lig_features.iter())
        .chain(real_rec_features.iter())
        .cloned()
        .collect();
    let (mins, maxs) = feature_ranges(&normalization_features);

    // ── No stratum collapse — for small datasets just pool all ─────────
    if n_strata == 0 || candidate_features.len() < 100 {
        return construct_fake_pairs_pooled(
            real_lig_features,
            real_rec_features,
            n_pairs,
            candidate_features,
            n_fake_per_real,
            caliper_k,
            max_match_distance,
            &mins,
            &maxs,
            &mut rng,
        );
    }

    // ── Assign candidate genes to strata ───────────────────────────────
    let stratum_idx = assign_strata(candidate_features, n_strata, &mins, &maxs);

    // ── Build stratum → candidate gene index map ──────────────────────
    let n_cells = n_strata * n_strata * n_strata * n_strata;
    let mut stratum_map: Vec<Vec<usize>> = vec![Vec::new(); n_cells];
    for (i, &sid) in stratum_idx.iter().enumerate() {
        stratum_map[sid].push(i);
    }

    let mut out = Vec::with_capacity(n_pairs * n_fake_per_real);

    for p in 0..n_pairs {
        // Stratum for real L
        let ls = feature_to_stratum(&real_lig_features[p], n_strata, &mins, &maxs);
        // Stratum for real R
        let rs = feature_to_stratum(&real_rec_features[p], n_strata, &mins, &maxs);

        let lig_cands = &stratum_map[ls];
        let rec_cands = &stratum_map[rs];

        if lig_cands.is_empty() || rec_cands.is_empty() {
            continue; // stratum empty → skip this pair
        }

        // Top-k nearest in 4-D Euclidean distance within stratum
        let lig_topk = top_k_nearest(
            &real_lig_features[p],
            lig_cands,
            candidate_features,
            caliper_k.min(lig_cands.len()),
            &mins,
            &maxs,
            max_match_distance,
        );
        let rec_topk = top_k_nearest(
            &real_rec_features[p],
            rec_cands,
            candidate_features,
            caliper_k.min(rec_cands.len()),
            &mins,
            &maxs,
            max_match_distance,
        );

        if lig_topk.is_empty() || rec_topk.is_empty() {
            continue;
        }

        for _ in 0..n_fake_per_real {
            for _ in 0..50 {
                let li = lig_topk[rng.gen_range(0..lig_topk.len())];
                let ri = rec_topk[rng.gen_range(0..rec_topk.len())];
                let lg = candidate_features[li].gene_idx;
                let rg = candidate_features[ri].gene_idx;
                if lg != rg {
                    out.push(FakeLrPair {
                        ligand_gene: lg,
                        receptor_gene: rg,
                        complex_ligand_rows: None,
                        complex_receptor_rows: None,
                        shadowing_pair_idx: p,
                        source: FAKE_CTX_WILDCARD,
                        target: FAKE_CTX_WILDCARD,
                        lineage: FAKE_CTX_WILDCARD,
                        null_kind: FAKE_NULL_GENERIC,
                        ligand_match_distance: norm_euclid_sq(
                            &feature_vec4(&real_lig_features[p]),
                            &feature_vec4(&candidate_features[li]),
                            &mins,
                            &maxs,
                        )
                        .sqrt(),
                        receptor_match_distance: norm_euclid_sq(
                            &feature_vec4(&real_rec_features[p]),
                            &feature_vec4(&candidate_features[ri]),
                            &mins,
                            &maxs,
                        )
                        .sqrt(),
                    });
                    break;
                }
            }
        }
    }
    out
}

/// Pooled (no-stratum) version — fallback when stratification is not
/// feasible (few genes, or `n_strata = 0`).
fn construct_fake_pairs_pooled(
    real_lig_features: &[GeneFeatures],
    real_rec_features: &[GeneFeatures],
    n_pairs: usize,
    candidate_features: &[GeneFeatures],
    n_fake_per_real: usize,
    caliper_k: usize,
    max_match_distance: f32,
    mins: &[f32],
    maxs: &[f32],
    rng: &mut StdRng,
) -> Vec<FakeLrPair> {
    let all_idx: Vec<usize> = (0..candidate_features.len()).collect();
    let mut out = Vec::with_capacity(n_pairs * n_fake_per_real);

    for p in 0..n_pairs {
        let lig_topk = top_k_nearest(
            &real_lig_features[p],
            &all_idx,
            candidate_features,
            caliper_k.min(all_idx.len()),
            mins,
            maxs,
            max_match_distance,
        );
        let rec_topk = top_k_nearest(
            &real_rec_features[p],
            &all_idx,
            candidate_features,
            caliper_k.min(all_idx.len()),
            mins,
            maxs,
            max_match_distance,
        );
        if lig_topk.is_empty() || rec_topk.is_empty() {
            continue;
        }
        for _ in 0..n_fake_per_real {
            for _ in 0..50 {
                let li = lig_topk[rng.gen_range(0..lig_topk.len())];
                let ri = rec_topk[rng.gen_range(0..rec_topk.len())];
                let lg = candidate_features[li].gene_idx;
                let rg = candidate_features[ri].gene_idx;
                if lg != rg {
                    out.push(FakeLrPair {
                        ligand_gene: lg,
                        receptor_gene: rg,
                        complex_ligand_rows: None,
                        complex_receptor_rows: None,
                        shadowing_pair_idx: p,
                        source: FAKE_CTX_WILDCARD,
                        target: FAKE_CTX_WILDCARD,
                        lineage: FAKE_CTX_WILDCARD,
                        null_kind: FAKE_NULL_GENERIC,
                        ligand_match_distance: norm_euclid_sq(
                            &feature_vec4(&real_lig_features[p]),
                            &feature_vec4(&candidate_features[li]),
                            mins,
                            maxs,
                        )
                        .sqrt(),
                        receptor_match_distance: norm_euclid_sq(
                            &feature_vec4(&real_rec_features[p]),
                            &feature_vec4(&candidate_features[ri]),
                            mins,
                            maxs,
                        )
                        .sqrt(),
                    });
                    break;
                }
            }
        }
    }
    out
}

// ─── Scoring ───────────────────────────────────────────────────────────

/// Score fake LR pairs using pre-computed gene bin statistics.
///
/// For each fake pair `(L', R')`, looks up `mean_L'[ct, lin, bin]` and
/// `mean_R'[ct, lin, bin]` from `gene_stats`, computes Michaelis-Menten
/// scores `√(L'·R') / (hill_s + √(L'·R'))` per (ct_s, ct_t, lin, bin),
/// and returns `LrBinResult` rows.
///
/// This is O(|D*| × n_types² × n_lineages × n_bins) — NO cell walk.
pub fn score_fake_pairs_bins(
    fake_pairs: &[FakeLrPair],
    real_stats: &GeneBinStats,
    fake_stats: &GeneBinStats,
    hill_s: f32,
    n_types: usize,
    n_lineages: usize,
    n_bins: usize,
    min_pct: f32,
) -> Vec<FakeLrBinResult> {
    let mut out = Vec::new();
    let slb = n_lineages * n_bins;
    let sb = n_bins;

    for (fi, fp) in fake_pairs.iter().enumerate() {
        let (sources, targets, lineages): (Vec<usize>, Vec<usize>, Vec<usize>) =
            if fp.source == FAKE_CTX_WILDCARD {
                (
                    (0..n_types).collect(),
                    (0..n_types).collect(),
                    (0..n_lineages).collect(),
                )
            } else {
                (vec![fp.source], vec![fp.target], vec![fp.lineage])
            };

        for &lin in &lineages {
            for bin in 0..n_bins {
                for &s in &sources {
                    let lm = fp
                        .ligand_rows()
                        .iter()
                        .map(|&g| {
                            GeneBinStats::get_mean_resolved(real_stats, fake_stats, g, s, lin, bin)
                        })
                        .fold(f32::INFINITY, f32::min);
                    let lp = fp
                        .ligand_rows()
                        .iter()
                        .map(|&g| {
                            GeneBinStats::get_pct_resolved(real_stats, fake_stats, g, s, lin, bin)
                        })
                        .fold(f32::INFINITY, f32::min);
                    if lp < min_pct {
                        continue;
                    }
                    for &t in &targets {
                        let rm = fp
                            .receptor_rows()
                            .iter()
                            .map(|&g| {
                                GeneBinStats::get_mean_resolved(
                                    real_stats, fake_stats, g, t, lin, bin,
                                )
                            })
                            .fold(f32::INFINITY, f32::min);
                        let rp = fp
                            .receptor_rows()
                            .iter()
                            .map(|&g| {
                                GeneBinStats::get_pct_resolved(
                                    real_stats, fake_stats, g, t, lin, bin,
                                )
                            })
                            .fold(f32::INFINITY, f32::min);
                        if rp < min_pct {
                            continue;
                        }
                        if lm > 0.0 && rm > 0.0 {
                            let geo = (lm * rm).sqrt();
                            let ws = real_stats.wsums[s * slb + lin * sb + bin] as f32;
                            let wt = real_stats.wsums[t * slb + lin * sb + bin] as f32;
                            let w2s = real_stats.w2sums[s * slb + lin * sb + bin];
                            let w2t = real_stats.w2sums[t * slb + lin * sb + bin];
                            let ess_s = if w2s > 1e-12 {
                                (ws as f64 * ws as f64 / w2s) as f32
                            } else {
                                0.0
                            };
                            let ess_t = if w2t > 1e-12 {
                                (wt as f64 * wt as f64 / w2t) as f32
                            } else {
                                0.0
                            };
                            out.push(FakeLrBinResult {
                                fake_idx: fi,
                                pair_idx: fp.shadowing_pair_idx,
                                source: s,
                                target: t,
                                lineage: lin,
                                bin,
                                score: geo / (hill_s + geo),
                                n_eff_cells: ess_s.min(ess_t),
                                weight_mass_cells: ws.min(wt),
                            });
                        }
                    }
                }
            }
        }
    }
    out
}

/// Stream trajectory-level FakeLR calibration without materialising every
/// fake `(source, target, lineage, bin)` row.
pub fn calibrate_fake_trajectories(
    lr_pairs: &[crate::ccc::LrPair],
    real_scores: &[crate::ccc::LrBinResult],
    fake_pairs: &[FakeLrPair],
    real_grid: &crate::ccc::lr_scoring::LrBinGrid,
    real_stats: &GeneBinStats,
    fake_stats: &GeneBinStats,
    min_pct: f32,
    min_n_eff: f32,
    min_valid_bins: usize,
    min_valid_fraction: f32,
    stat_code: usize,
    top_k: usize,
    n_fake_expected: usize,
    null_mode: FakeNullMode,
) -> Vec<FakeTrajectoryCalibration> {
    use crate::ccc::trajectory_curve::{
        build_fake_complex_trajectory_curve, build_fake_trajectory_curve,
        build_lr_trajectory_curve_from_grid, evaluate_trajectory_curve,
        trajectory_keys_from_scores, TrajectoryEvalConfig,
    };

    let cfg = TrajectoryEvalConfig {
        min_pct,
        min_n_eff,
        min_valid_bins,
        min_valid_fraction,
    };

    let keys = trajectory_keys_from_scores(real_scores);
    let key_to_idx: std::collections::HashMap<crate::ccc::trajectory_curve::TrajKey, usize> =
        keys.iter().enumerate().map(|(i, key)| (*key, i)).collect();

    let real_eval: Vec<_> = keys
        .iter()
        .map(|key| {
            let pair = &lr_pairs[key.0];
            let curve = build_lr_trajectory_curve_from_grid(
                pair, real_grid, key.1, key.2, key.3, min_pct, min_n_eff,
            );
            evaluate_trajectory_curve(&curve, stat_code, top_k, &cfg)
        })
        .collect();

    let (n_generic_exp, n_edge_exp, n_network_exp) =
        split_null_budget_layers(n_fake_expected, null_mode);
    let (n_replace_ligand_exp, n_replace_receptor_exp) = split_edge_budget(n_edge_exp);
    let (n_network_ligand_exp, n_network_receptor_exp) = split_edge_budget(n_network_exp);

    let mut observed_all = vec![0usize; keys.len()];
    let mut count_geq_generic = vec![0usize; keys.len()];
    let mut observed_generic = vec![0usize; keys.len()];
    let mut count_geq_edge = vec![0usize; keys.len()];
    let mut observed_edge = vec![0usize; keys.len()];
    let mut count_geq_replace_ligand = vec![0usize; keys.len()];
    let mut observed_replace_ligand = vec![0usize; keys.len()];
    let mut count_geq_replace_receptor = vec![0usize; keys.len()];
    let mut observed_replace_receptor = vec![0usize; keys.len()];
    let mut count_geq_network = vec![0usize; keys.len()];
    let mut observed_network = vec![0usize; keys.len()];
    let mut count_geq_network_ligand = vec![0usize; keys.len()];
    let mut observed_network_ligand = vec![0usize; keys.len()];
    let mut count_geq_network_receptor = vec![0usize; keys.len()];
    let mut observed_network_receptor = vec![0usize; keys.len()];
    let mut lig_dist_sum = vec![0f32; keys.len()];
    let mut rec_dist_sum = vec![0f32; keys.len()];
    let mut dist_n = vec![0usize; keys.len()];
    let mut unique_draws: Vec<std::collections::HashSet<(u64, u64)>> = (0..keys.len())
        .map(|_| std::collections::HashSet::new())
        .collect();
    let mut unique_draws_generic: Vec<std::collections::HashSet<(u64, u64)>> = (0..keys.len())
        .map(|_| std::collections::HashSet::new())
        .collect();
    let mut unique_draws_edge: Vec<std::collections::HashSet<(u64, u64)>> = (0..keys.len())
        .map(|_| std::collections::HashSet::new())
        .collect();
    let mut unique_draws_replace_ligand: Vec<std::collections::HashSet<(u64, u64)>> = (0..keys
        .len())
        .map(|_| std::collections::HashSet::new())
        .collect();
    let mut unique_draws_replace_receptor: Vec<std::collections::HashSet<(u64, u64)>> = (0..keys
        .len())
        .map(|_| std::collections::HashSet::new())
        .collect();
    let mut unique_draws_network: Vec<std::collections::HashSet<(u64, u64)>> = (0..keys.len())
        .map(|_| std::collections::HashSet::new())
        .collect();
    let mut unique_draws_network_ligand: Vec<std::collections::HashSet<(u64, u64)>> = (0..keys
        .len())
        .map(|_| std::collections::HashSet::new())
        .collect();
    let mut unique_draws_network_receptor: Vec<std::collections::HashSet<(u64, u64)>> = (0..keys
        .len())
        .map(|_| std::collections::HashSet::new())
        .collect();
    let mut traj_fake_stats: Vec<Vec<f32>> = (0..keys.len()).map(|_| Vec::new()).collect();

    let mut evaluate_fake =
        |fp: &FakeLrPair, idx: usize, source: usize, target: usize, lin: usize| {
            if !real_eval[idx].evaluable {
                return;
            }
            let curve = if fp.is_complex() {
                build_fake_complex_trajectory_curve(
                    fp.ligand_rows(),
                    fp.receptor_rows(),
                    source,
                    target,
                    lin,
                    real_stats,
                    fake_stats,
                    real_grid,
                    min_pct,
                    min_n_eff,
                )
            } else {
                build_fake_trajectory_curve(
                    fp.ligand_gene,
                    fp.receptor_gene,
                    source,
                    target,
                    lin,
                    real_stats,
                    fake_stats,
                    real_grid,
                    min_pct,
                    min_n_eff,
                )
            };
            let fake_eval = evaluate_trajectory_curve(&curve, stat_code, top_k, &cfg);
            if !fake_eval.evaluable {
                return;
            }
            traj_fake_stats[idx].push(fake_eval.statistic);
            let draw_key = fp.draw_key();
            observed_all[idx] += 1;
            unique_draws[idx].insert(draw_key);
            lig_dist_sum[idx] += fp.ligand_match_distance;
            rec_dist_sum[idx] += fp.receptor_match_distance;
            dist_n[idx] += 1;
            if fp.null_kind == FAKE_NULL_GENERIC {
                observed_generic[idx] += 1;
                unique_draws_generic[idx].insert(draw_key);
                if fake_eval.statistic >= real_eval[idx].statistic {
                    count_geq_generic[idx] += 1;
                }
            } else if fp.null_kind == FAKE_NULL_REWIRE_LIGAND {
                observed_edge[idx] += 1;
                observed_replace_ligand[idx] += 1;
                unique_draws_edge[idx].insert(draw_key);
                unique_draws_replace_ligand[idx].insert(draw_key);
                if fake_eval.statistic >= real_eval[idx].statistic {
                    count_geq_edge[idx] += 1;
                    count_geq_replace_ligand[idx] += 1;
                }
            } else if fp.null_kind == FAKE_NULL_REWIRE_RECEPTOR {
                observed_edge[idx] += 1;
                observed_replace_receptor[idx] += 1;
                unique_draws_edge[idx].insert(draw_key);
                unique_draws_replace_receptor[idx].insert(draw_key);
                if fake_eval.statistic >= real_eval[idx].statistic {
                    count_geq_edge[idx] += 1;
                    count_geq_replace_receptor[idx] += 1;
                }
            } else if fp.null_kind == FAKE_NULL_NETWORK_REWIRE_LIGAND {
                observed_network[idx] += 1;
                observed_network_ligand[idx] += 1;
                unique_draws_network[idx].insert(draw_key);
                unique_draws_network_ligand[idx].insert(draw_key);
                if fake_eval.statistic >= real_eval[idx].statistic {
                    count_geq_network[idx] += 1;
                    count_geq_network_ligand[idx] += 1;
                }
            } else if fp.null_kind == FAKE_NULL_NETWORK_REWIRE_RECEPTOR {
                observed_network[idx] += 1;
                observed_network_receptor[idx] += 1;
                unique_draws_network[idx].insert(draw_key);
                unique_draws_network_receptor[idx].insert(draw_key);
                if fake_eval.statistic >= real_eval[idx].statistic {
                    count_geq_network[idx] += 1;
                    count_geq_network_receptor[idx] += 1;
                }
            }
        };

    for fp in fake_pairs {
        if fp.source == FAKE_CTX_WILDCARD {
            continue;
        }
        let key = (fp.shadowing_pair_idx, fp.source, fp.target, fp.lineage);
        let Some(&idx) = key_to_idx.get(&key) else {
            continue;
        };
        evaluate_fake(fp, idx, fp.source, fp.target, fp.lineage);
    }

    let laplace = |geq: usize, obs: usize| -> f64 {
        if obs > 0 {
            (geq as f64 + 1.0) / (obs as f64 + 1.0)
        } else {
            f64::NAN
        }
    };

    keys.iter()
        .enumerate()
        .map(|(i, key)| {
            let ev = &real_eval[i];
            let coverage = if n_fake_expected > 0 {
                observed_all[i] as f64 / n_fake_expected as f64
            } else {
                f64::NAN
            };
            let pvalue_generic = laplace(count_geq_generic[i], observed_generic[i]);
            let pvalue_edge = laplace(count_geq_edge[i], observed_edge[i]);
            let pvalue_replace_ligand =
                laplace(count_geq_replace_ligand[i], observed_replace_ligand[i]);
            let pvalue_replace_receptor =
                laplace(count_geq_replace_receptor[i], observed_replace_receptor[i]);
            let pvalue_network = laplace(count_geq_network[i], observed_network[i]);
            let pvalue_network_rewire_ligand =
                laplace(count_geq_network_ligand[i], observed_network_ligand[i]);
            let pvalue_network_rewire_receptor =
                laplace(count_geq_network_receptor[i], observed_network_receptor[i]);
            let pair = &lr_pairs[key.0];
            let pvalue_complex = if is_complex_lr_pair(pair) {
                pvalue_generic
            } else {
                f64::NAN
            };
            let fake_null_stat_median = if traj_fake_stats[i].is_empty() {
                f32::NAN
            } else {
                let mut stats = traj_fake_stats[i].clone();
                stats.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
                let m = stats.len() / 2;
                if stats.len() % 2 == 1 {
                    stats[m]
                } else {
                    (stats[m - 1] + stats[m]) / 2.0
                }
            };
            // Mixed-layer modes (`Both`, `All`) leave `pvalue` unset here; R resolves
            // the primary p_fake hierarchy (All: p_network_iu > p_edge_iu > p_fake_generic).
            let pvalue = match null_mode {
                FakeNullMode::ContextGeneric => pvalue_generic,
                FakeNullMode::EdgeRewire => pvalue_edge,
                FakeNullMode::Both | FakeNullMode::All => f64::NAN,
                FakeNullMode::LrNetworkRewire => pvalue_network,
            };
            let statistic = if ev.evaluable { ev.statistic } else { f32::NAN };
            let (mean_lig, mean_rec) = if dist_n[i] > 0 {
                (
                    lig_dist_sum[i] / dist_n[i] as f32,
                    rec_dist_sum[i] / dist_n[i] as f32,
                )
            } else {
                (f32::NAN, f32::NAN)
            };
            FakeTrajectoryCalibration {
                pair_idx: key.0,
                source: key.1,
                target: key.2,
                lineage: key.3,
                statistic,
                pvalue,
                n_fake_expected,
                n_fake_observed: observed_all[i],
                fake_coverage: coverage,
                n_valid_bins: ev.n_valid_bins,
                valid_bin_fraction: ev.valid_bin_fraction,
                n_expression_supported_bins: ev.n_expression_supported_bins,
                n_eff_at_max: ev.n_eff_at_max,
                n_eff_at_min: ev.n_eff_at_min,
                valid_at_max: ev.valid_at_max,
                valid_at_min: ev.valid_at_min,
                trajectory_evaluable: ev.evaluable,
                pvalue_generic,
                pvalue_edge,
                pvalue_replace_ligand,
                pvalue_replace_receptor,
                mean_ligand_match_distance: mean_lig,
                mean_receptor_match_distance: mean_rec,
                n_unique_fake_draws: unique_draws[i].len(),
                n_unique_fake_draws_generic: unique_draws_generic[i].len(),
                n_unique_fake_draws_edge: unique_draws_edge[i].len(),
                n_fake_expected_generic: n_generic_exp,
                n_fake_expected_edge: n_edge_exp,
                n_fake_observed_generic: observed_generic[i],
                n_fake_observed_edge: observed_edge[i],
                fake_coverage_generic: if n_generic_exp > 0 {
                    observed_generic[i] as f64 / n_generic_exp as f64
                } else {
                    f64::NAN
                },
                fake_coverage_edge: if n_edge_exp > 0 {
                    observed_edge[i] as f64 / n_edge_exp as f64
                } else {
                    f64::NAN
                },
                n_fake_expected_replace_ligand: n_replace_ligand_exp,
                n_fake_expected_replace_receptor: n_replace_receptor_exp,
                n_fake_observed_replace_ligand: observed_replace_ligand[i],
                n_fake_observed_replace_receptor: observed_replace_receptor[i],
                fake_coverage_replace_ligand: if n_replace_ligand_exp > 0 {
                    observed_replace_ligand[i] as f64 / n_replace_ligand_exp as f64
                } else {
                    f64::NAN
                },
                fake_coverage_replace_receptor: if n_replace_receptor_exp > 0 {
                    observed_replace_receptor[i] as f64 / n_replace_receptor_exp as f64
                } else {
                    f64::NAN
                },
                n_unique_fake_draws_replace_ligand: unique_draws_replace_ligand[i].len(),
                n_unique_fake_draws_replace_receptor: unique_draws_replace_receptor[i].len(),
                pvalue_network,
                pvalue_network_rewire_ligand,
                pvalue_network_rewire_receptor,
                n_fake_expected_network: n_network_exp,
                n_fake_observed_network: observed_network[i],
                fake_coverage_network: if n_network_exp > 0 {
                    observed_network[i] as f64 / n_network_exp as f64
                } else {
                    f64::NAN
                },
                n_unique_fake_draws_network: unique_draws_network[i].len(),
                n_fake_expected_network_rewire_ligand: n_network_ligand_exp,
                n_fake_expected_network_rewire_receptor: n_network_receptor_exp,
                n_fake_observed_network_rewire_ligand: observed_network_ligand[i],
                n_fake_observed_network_rewire_receptor: observed_network_receptor[i],
                fake_coverage_network_rewire_ligand: if n_network_ligand_exp > 0 {
                    observed_network_ligand[i] as f64 / n_network_ligand_exp as f64
                } else {
                    f64::NAN
                },
                fake_coverage_network_rewire_receptor: if n_network_receptor_exp > 0 {
                    observed_network_receptor[i] as f64 / n_network_receptor_exp as f64
                } else {
                    f64::NAN
                },
                n_unique_fake_draws_network_rewire_ligand: unique_draws_network_ligand[i].len(),
                n_unique_fake_draws_network_rewire_receptor: unique_draws_network_receptor[i].len(),
                pvalue_complex,
                fake_null_stat_median,
            }
        })
        .collect()
}

/// Compute global Hill saturation scale S from a GeneBinStats mean
/// array (arithmetic mean of all non-zero entries).
pub fn compute_hill_s_from_stats(stats: &GeneBinStats) -> f32 {
    let mut sum = 0.0f64;
    let mut cnt = 0u64;
    for &v in &stats.means {
        if v > 0.0 {
            sum += v as f64;
            cnt += 1;
        }
    }
    if cnt == 0 {
        1.0
    } else {
        (sum / cnt as f64) as f32
    }
}

// ─── p_fake computation ────────────────────────────────────────────────

/// Compute per-real-LR-row p_fake values.
///
/// For each real `LrBinResult` row, count how many fake rows (across the
/// same `(source, target, lineage, bin, pair_idx)` key, or pooled across
/// all fake rows shadowing the same `pair_idx`) score ≥ the real score,
/// then apply Laplace smoothing:
///
///   p_fake = (count_geq + 1) / (n_fake_observed + 1)
///
/// Returns a Vec<f64> of length `n_real_rows`, ordered the same as `real_scores`.
///
/// Matching strategy: for each real row, only compare against
/// fake rows that share the same `(source, target, lineage, bin)` AND
/// are from fake pairs that shadow the SAME real pair_idx.  This is the
/// strictest pairwise comparison — "this real LR at this position vs.
/// its fake shadows at the same position."
pub fn compute_p_fake_by_position(
    real_scores: &[crate::ccc::LrBinResult],
    fake_scores: &[crate::ccc::LrBinResult],
    _n_fake_pairs: usize,
) -> Vec<f64> {
    // Key: (pair_idx, source, target, lineage, bin)
    // For each fake row, accumulate a score list keyed by position.
    let mut fake_by_pos: std::collections::HashMap<(usize, usize, usize, usize, usize), Vec<f32>> =
        std::collections::HashMap::new();

    for f in fake_scores {
        let key = (f.pair_idx, f.source, f.target, f.lineage, f.bin);
        fake_by_pos.entry(key).or_default().push(f.score);
    }

    real_scores
        .iter()
        .map(|r| {
            let key = (r.pair_idx, r.source, r.target, r.lineage, r.bin);
            let fakes = fake_by_pos.get(&key);
            let Some(scores) = fakes else {
                return f64::NAN;
            };
            let count_geq = scores.iter().filter(|&&s| s >= r.score).count() as f64;
            (count_geq + 1.0) / (scores.len() as f64 + 1.0)
        })
        .collect()
}

// ─── Helpers ───────────────────────────────────────────────────────────

/// Compute per-axis [min, max] from a set of GeneFeatures.
fn feature_ranges(features: &[GeneFeatures]) -> (Vec<f32>, Vec<f32>) {
    if features.is_empty() {
        return (vec![0.0; 4], vec![1.0; 4]);
    }
    // 4 Dims: mean_global, pct_global, tau, trend
    let mut mins = [f32::INFINITY; 4];
    let mut maxs = [f32::NEG_INFINITY; 4];
    for f in features {
        let m: f32 = if f.ct_mean.is_empty() {
            0.0
        } else {
            f.ct_mean.iter().sum::<f32>() / f.ct_mean.len() as f32
        };
        let p: f32 = if f.ct_pct.is_empty() {
            0.0
        } else {
            f.ct_pct.iter().sum::<f32>() / f.ct_pct.len() as f32
        };
        for (i, &v) in [m, p, f.tau, f.trend].iter().enumerate() {
            if v < mins[i] {
                mins[i] = v;
            }
            if v > maxs[i] {
                maxs[i] = v;
            }
        }
    }
    (mins.to_vec(), maxs.to_vec())
}

/// 4-D feature vector for a GeneFeatures.
#[inline]
fn feature_vec4(f: &GeneFeatures) -> [f32; 4] {
    let m: f32 = if f.ct_mean.is_empty() {
        0.0
    } else {
        f.ct_mean.iter().sum::<f32>() / f.ct_mean.len() as f32
    };
    let p: f32 = if f.ct_pct.is_empty() {
        0.0
    } else {
        f.ct_pct.iter().sum::<f32>() / f.ct_pct.len() as f32
    };
    [m, p, f.tau, f.trend]
}

/// Compute normalised squared Euclidean distance between two feature vecs.
fn norm_euclid_sq(a: &[f32; 4], b: &[f32; 4], mins: &[f32], maxs: &[f32]) -> f32 {
    let mut d2 = 0.0f32;
    for i in 0..4 {
        let range = maxs[i] - mins[i];
        if range > 1e-8 {
            let da = (a[i] - mins[i]) / range;
            let db = (b[i] - mins[i]) / range;
            d2 += (da - db) * (da - db);
        }
    }
    d2
}

/// Assign each feature set to a discrete stratum.
fn assign_strata(
    features: &[GeneFeatures],
    n_strata: usize,
    mins: &[f32],
    maxs: &[f32],
) -> Vec<usize> {
    features
        .iter()
        .map(|f| feature_to_stratum(f, n_strata, mins, maxs))
        .collect()
}

fn feature_to_stratum(f: &GeneFeatures, n_strata: usize, mins: &[f32], maxs: &[f32]) -> usize {
    let vec = feature_vec4(f);
    let ns = n_strata.max(1);
    let mut idx = 0usize;
    for i in 0..4 {
        let range = maxs[i] - mins[i];
        let s = if range > 1e-8 {
            ((vec[i] - mins[i]) / range * ns as f32) as usize
        } else {
            0
        };
        idx = idx * ns + s.min(ns - 1);
    }
    idx
}

/// Return indices of top-k nearest candidate genes (by normalised
/// Euclidean distance to `target`) from `pool`.
fn top_k_nearest(
    target: &GeneFeatures,
    pool: &[usize],
    features: &[GeneFeatures],
    k: usize,
    mins: &[f32],
    maxs: &[f32],
    max_match_distance: f32,
) -> Vec<usize> {
    let tv = feature_vec4(target);
    let mut dists: Vec<(f32, usize)> = pool
        .iter()
        .map(|&idx| {
            (
                norm_euclid_sq(&tv, &feature_vec4(&features[idx]), mins, maxs),
                idx,
            )
        })
        .filter(|(distance_sq, _)| *distance_sq <= max_match_distance * max_match_distance)
        .collect();
    dists.sort_unstable_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    dists.truncate(k);
    dists.into_iter().map(|(_, idx)| idx).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn network_mode_allocates_full_budget_to_network_layer() {
        let (generic, edge, network) = split_null_budget_layers(5, FakeNullMode::LrNetworkRewire);
        assert_eq!((generic, edge, network), (0, 0, 5));
    }

    #[test]
    fn both_mode_keeps_generic_edge_split_without_network() {
        let (generic, edge, network) = split_null_budget_layers(3, FakeNullMode::Both);
        assert_eq!((generic, edge, network), (1, 2, 0));
    }

    #[test]
    fn all_mode_splits_budget_across_three_layers() {
        for total in 1..=12 {
            let (generic, edge, network) = split_null_budget_layers(total, FakeNullMode::All);
            assert_eq!(
                generic + edge + network,
                total,
                "budget not conserved for total={total}"
            );
            assert!(generic <= edge + 1 || total <= 2);
            assert!(edge <= network + 1 || total <= 2);
        }
        assert_eq!(split_null_budget_layers(6, FakeNullMode::All), (2, 2, 2));
        assert_eq!(split_null_budget_layers(7, FakeNullMode::All), (2, 2, 3));
    }

    #[test]
    fn both_mode_splits_budget_once_before_block_allocation() {
        let weights = vec![1.0, 1.0, 1.0];
        let (generic_total, edge_total, network_total) =
            split_null_budget_layers(3, FakeNullMode::Both);
        let generic_by_block = allocate_draws_across_blocks(generic_total, &weights);
        let edge_by_block = allocate_draws_across_blocks(edge_total, &weights);

        assert_eq!((generic_total, edge_total, network_total), (1, 2, 0));
        assert_eq!(generic_by_block.iter().sum::<usize>(), 1);
        assert_eq!(edge_by_block.iter().sum::<usize>(), 2);
    }

    #[test]
    fn edge_budget_splits_ligand_and_receptor_draws() {
        assert_eq!(split_edge_budget(0), (0, 0));
        assert_eq!(split_edge_budget(1), (0, 1));
        assert_eq!(split_edge_budget(2), (1, 1));
        assert_eq!(split_edge_budget(3), (1, 2));
        assert_eq!(split_edge_budget(4), (2, 2));
    }

    #[test]
    fn block_matching_requires_sender_and_receiver_support() {
        let block_labels = vec![0usize, 1, 1];
        let cell_types = vec![0usize, 0, 1];
        let lineage_weights = vec![1.0f32, 0.7, 0.3];
        let empty_stats: Vec<GeneBinStats> = Vec::new();
        let blocks = FakeMatchingBlocks {
            block_labels: &block_labels,
            cell_types: &cell_types,
            lineage_weights: &lineage_weights,
            n_lineages: 1,
            n_blocks: 2,
            real_stats: &empty_stats,
            fake_stats: &empty_stats,
        };

        let weights = trajectory_block_weights(&blocks, 0, 1, 0);
        assert_eq!(weights[0], 0.0);
        assert!((weights[1] - 0.3).abs() < 1e-6);
    }

    #[test]
    fn block_allocation_with_no_shared_support_requests_no_draws() {
        assert_eq!(allocate_draws_across_blocks(7, &[0.0, 0.0]), vec![0, 0]);
    }

    // ─── GeneBinStats layout ─────────────────────────────────────────

    #[test]
    fn gene_bin_stats_deduplicates_genes_and_counts_positive_values() {
        let expr = SparseMatrix::new_csc(
            vec![0, 1, 0, 0],
            vec![0, 2, 3, 4],
            vec![0.0, -2.0, 1.0, 2.0],
            2,
            3,
        );
        let stats = compute_gene_bin_stats(
            &expr,
            &[0, 0, 0],
            1,
            &[0, 0, 1],
            &[0.0, 0.5, 1.0],
            &[1.0, 1.0, 1.0],
            1,
            1,
            None,
            None,
        );
        assert_eq!(stats.gene_rows, vec![0, 1]);
        assert_eq!(stats.n_genes, 2);
        assert!((stats.get_pct_global(0, 0, 0, 0) - 2.0 / 3.0).abs() < 1e-6);
        assert_eq!(stats.get_pct_global(1, 0, 0, 0), 0.0);
    }

    #[test]
    fn gene_bin_stats_idx_matches_documented() {
        let stats = GeneBinStats {
            means: vec![0.0; 3 * 2 * 2 * 4], // 3 genes × 2 ct × 2 lin × 4 bins
            pcts: vec![0.0; 3 * 2 * 2 * 4],
            wsums: vec![0.0; 2 * 2 * 4],
            w2sums: vec![0.0; 2 * 2 * 4],
            gene_rows: vec![0, 1, 2],
            n_genes: 3,
            n_types: 2,
            n_lineages: 2,
            n_bins: 4,
        };
        // gene=0, ct=0, lin=0, bin=0
        assert_eq!(stats.idx(0, 0, 0, 0), 0);
        // gene=2, ct=1, lin=1, bin=3
        // = 2 * (2*2*4) + 1 * (2*4) + 1 * 4 + 3 = 32 + 8 + 4 + 3 = 47
        assert_eq!(stats.idx(2, 1, 1, 3), 47);
    }

    // ─── Feature computation ─────────────────────────────────────────

    #[test]
    fn compute_gene_features_uses_global_gene_rows() {
        let stats = GeneBinStats {
            means: vec![1.0, 2.0, 3.0, 4.0],
            pcts: vec![0.5; 4],
            wsums: vec![1.0; 4],
            w2sums: vec![1.0; 4],
            gene_rows: vec![10, 11],
            n_genes: 2,
            n_types: 1,
            n_lineages: 1,
            n_bins: 2,
        };
        let feats = compute_gene_features(&stats, &[10, 11]);
        assert_eq!(feats.len(), 2);
        assert_eq!(feats[0].gene_idx, 10);
        assert!((feats[0].ct_mean[0] - 1.5).abs() < 1e-6);
        assert!((feats[1].ct_mean[0] - 3.5).abs() < 1e-6);
    }

    #[test]
    fn compute_gene_features_single_gene_uniform() {
        // 1 type, 1 lineage, 1 bin — trivial case.
        let stats = GeneBinStats {
            means: vec![5.0],
            pcts: vec![0.5],
            wsums: vec![1.0],
            w2sums: vec![1.0],
            gene_rows: vec![0],
            n_genes: 1,
            n_types: 1,
            n_lineages: 1,
            n_bins: 1,
        };
        let feats = compute_gene_features(&stats, &[0]);
        assert_eq!(feats.len(), 1);
        assert_eq!(feats[0].gene_idx, 0);
        assert_eq!(feats[0].ct_mean, vec![5.0]);
        assert_eq!(feats[0].ct_pct, vec![0.5]);
        assert_eq!(feats[0].tau, 0.0); // 1 type ⇒ tau = 0
        assert_eq!(feats[0].trend, 0.0); // 1 bin ⇒ can't compute trend
    }

    #[test]
    fn gene_trend_uses_weighted_expression_across_cell_types() {
        let n_types = 2usize;
        let n_lineages = 1usize;
        let n_bins = 5usize;
        let mut means = Vec::new();
        let pcts = vec![1.0f32; n_types * n_lineages * n_bins];
        let wsums = vec![1.0f64; n_types * n_lineages * n_bins];
        for ct in 0..n_types {
            for bin in 0..n_bins {
                means.push(if ct == 1 { (bin + 1) as f32 } else { 0.0 });
            }
        }
        let stats = GeneBinStats {
            means,
            pcts,
            w2sums: wsums.clone(),
            wsums,
            gene_rows: vec![0],
            n_genes: 1,
            n_types,
            n_lineages,
            n_bins,
        };
        let feats = compute_gene_features(&stats, &[0]);
        assert!(
            feats[0].trend > 0.5,
            "trend should reflect non-CT0 signal, got {}",
            feats[0].trend
        );
    }

    #[test]
    fn gene_trend_invariant_to_cell_type_label_permutation() {
        fn make_stats(ct_first: &[f32; 5], ct_second: &[f32; 5]) -> GeneBinStats {
            let n_types = 2usize;
            let n_lineages = 1usize;
            let n_bins = 5usize;
            let mut means = Vec::new();
            for vals in [ct_first, ct_second] {
                for &v in vals {
                    means.push(v);
                }
            }
            GeneBinStats {
                means,
                pcts: vec![1.0f32; n_types * n_lineages * n_bins],
                wsums: vec![1.0f64; n_types * n_lineages * n_bins],
                w2sums: vec![1.0f64; n_types * n_lineages * n_bins],
                gene_rows: vec![0],
                n_genes: 1,
                n_types,
                n_lineages,
                n_bins,
            }
        }
        let increasing = [1.0f32, 2.0, 3.0, 4.0, 5.0];
        let flat = [0.0f32; 5];
        let a = compute_gene_features(&make_stats(&flat, &increasing), &[0]);
        let b = compute_gene_features(&make_stats(&increasing, &flat), &[0]);
        assert!(
            (a[0].trend - b[0].trend).abs() < 1e-5,
            "trend must not depend on cell-type label order"
        );
    }

    #[test]
    fn compute_gene_features_tau_specificity() {
        // 3 types: means [1.0, 2.0, 7.0], max=7
        // tau = sum(1 - m/max) / (n-1) = ((1-1/7)+(1-2/7)+(1-7/7))/2
        //   = (6/7 + 5/7 + 0) / 2 = (11/7)/2 = 11/14 ≈ 0.7857
        let stats = GeneBinStats {
            means: vec![1.0, 2.0, 7.0],
            pcts: vec![0.1, 0.2, 0.7],
            wsums: vec![1.0, 1.0, 1.0],
            w2sums: vec![1.0, 1.0, 1.0],
            gene_rows: vec![0],
            n_genes: 1,
            n_types: 3,
            n_lineages: 1,
            n_bins: 1,
        };
        let feats = compute_gene_features(&stats, &[0]);
        let expected_tau = (1.0 - 1.0 / 7.0 + 1.0 - 2.0 / 7.0 + 0.0) / 2.0; // ≈ 0.857...
        assert!(
            (feats[0].tau - expected_tau).abs() < 1e-4,
            "tau = {}, expected ≈ {}",
            feats[0].tau,
            expected_tau
        );
    }

    // ─── Stratified matching ─────────────────────────────────────────

    fn dummy_features(idx: usize, vals: [f32; 4]) -> GeneFeatures {
        GeneFeatures {
            gene_idx: idx,
            ct_mean: vec![vals[0]],
            ct_pct: vec![vals[1]],
            tau: vals[2],
            trend: vals[3],
        }
    }

    #[test]
    fn construct_fake_pairs_simple() {
        // 1 real pair, 10 candidate genes, no strata.
        let real_lig = vec![dummy_features(0, [0.5, 0.3, 0.5, 0.2])];
        let real_rec = vec![dummy_features(1, [0.6, 0.4, 0.4, 0.3])];
        let candidates: Vec<GeneFeatures> = (10..20)
            .map(|i| dummy_features(i, [0.5 + i as f32 * 0.01, 0.3, 0.5, 0.2]))
            .collect();

        let pairs = construct_fake_pairs_stratified(
            &real_lig,
            &real_rec,
            1,
            &candidates,
            3,
            5,
            5,
            f32::INFINITY,
            42,
        );
        assert_eq!(pairs.len(), 3);
        for p in &pairs {
            assert_eq!(p.shadowing_pair_idx, 0);
            assert_ne!(p.ligand_gene, p.receptor_gene);
        }
    }

    #[test]
    fn construct_fake_pairs_empty_candidates_returns_empty() {
        let real_lig = vec![dummy_features(0, [0.5, 0.5, 0.5, 0.5])];
        let real_rec = vec![dummy_features(0, [0.5, 0.5, 0.5, 0.5])];
        let pairs = construct_fake_pairs_stratified(
            &real_lig,
            &real_rec,
            1,
            &[],
            3,
            5,
            5,
            f32::INFINITY,
            42,
        );
        assert!(pairs.is_empty());
    }

    #[test]
    fn construct_fake_pairs_never_falls_back_outside_caliper() {
        let real_lig = vec![dummy_features(0, [0.1, 0.1, 0.1, 0.1])];
        let real_rec = vec![dummy_features(1, [0.9, 0.9, 0.9, 0.9])];
        let candidates = vec![
            dummy_features(10, [0.1, 0.1, 0.1, 0.1]),
            dummy_features(11, [0.9, 0.9, 0.9, 0.9]),
            dummy_features(99, [0.5, 0.5, 0.5, 0.5]),
        ];

        let pairs = construct_fake_pairs_stratified(
            &real_lig,
            &real_rec,
            1,
            &candidates,
            20,
            0,
            1,
            f32::INFINITY,
            42,
        );
        assert_eq!(pairs.len(), 20);
        assert!(pairs.iter().all(|p| p.ligand_gene == 10));
        assert!(pairs.iter().all(|p| p.receptor_gene == 11));
    }

    #[test]
    fn construct_fake_pairs_enforces_numeric_distance_caliper() {
        let real_lig = vec![dummy_features(0, [0.0, 0.0, 0.0, 0.0])];
        let real_rec = vec![dummy_features(1, [1.0, 1.0, 1.0, 1.0])];
        let candidates = vec![
            dummy_features(10, [0.5, 0.5, 0.5, 0.5]),
            dummy_features(11, [0.5, 0.5, 0.5, 0.5]),
        ];
        let pairs = construct_fake_pairs_stratified(
            &real_lig,
            &real_rec,
            1,
            &candidates,
            10,
            0,
            2,
            0.25,
            42,
        );
        assert!(pairs.is_empty());
    }

    // ─── Hill S from stats ───────────────────────────────────────────

    #[test]
    fn compute_hill_s_from_stats_basic() {
        let stats = GeneBinStats {
            means: vec![1.0, 2.0, 0.0, 4.0, 3.0, 0.0],
            pcts: vec![0.0; 6],
            wsums: vec![0.0; 6],
            w2sums: vec![0.0; 6],
            gene_rows: vec![0],
            n_genes: 1,
            n_types: 6,
            n_lineages: 1,
            n_bins: 1,
        };
        assert!((compute_hill_s_from_stats(&stats) - 2.5).abs() < 1e-6);
    }

    #[test]
    fn compute_hill_s_from_stats_all_zero() {
        let stats = GeneBinStats {
            means: vec![0.0; 4],
            pcts: vec![0.0; 4],
            wsums: vec![0.0; 4],
            w2sums: vec![0.0; 4],
            gene_rows: vec![0],
            n_genes: 1,
            n_types: 4,
            n_lineages: 1,
            n_bins: 1,
        };
        assert_eq!(compute_hill_s_from_stats(&stats), 1.0);
    }

    // ─── p_fake computation ──────────────────────────────────────────

    #[test]
    fn compute_p_fake_null_case() {
        // No matched fake scores → p_fake is undefined.
        let real = vec![crate::ccc::LrBinResult {
            pair_idx: 0,
            source: 0,
            target: 1,
            lineage: 0,
            bin: 2,
            score: 0.8,
            n_eff_cells: 10.0,
            weight_mass_cells: 10.0,
            lineage_entropy: 0.3,
        }];
        let fake: Vec<crate::ccc::LrBinResult> = vec![];
        let pf = compute_p_fake_by_position(&real, &fake, 1);
        assert_eq!(pf.len(), 1);
        assert!(pf[0].is_nan());
    }

    #[test]
    fn compute_p_fake_all_lower() {
        // Real score = 0.9; all 9 fake scores < 0.9 → count_geq = 0
        // p_fake = (0 + 1) / (9 + 1) = 0.1
        let real = vec![crate::ccc::LrBinResult {
            pair_idx: 0,
            source: 0,
            target: 1,
            lineage: 0,
            bin: 2,
            score: 0.9,
            n_eff_cells: 10.0,
            weight_mass_cells: 10.0,
            lineage_entropy: 0.3,
        }];
        let mut fake = Vec::new();
        for i in 0..9 {
            fake.push(crate::ccc::LrBinResult {
                pair_idx: 0,
                source: 0,
                target: 1,
                lineage: 0,
                bin: 2,
                score: 0.1 + i as f32 * 0.05,
                n_eff_cells: 10.0,
                weight_mass_cells: 10.0,
                lineage_entropy: 0.3,
            });
        }
        let pf = compute_p_fake_by_position(&real, &fake, 9);
        assert!((pf[0] - 0.1).abs() < 1e-6, "p_fake = {}", pf[0]);
    }

    #[test]
    fn compute_p_fake_all_higher() {
        // Real score = 0.1; all 9 fake scores >= 0.1 → count_geq = 9
        // p_fake = (9 + 1) / (9 + 1) = 1.0
        let real = vec![crate::ccc::LrBinResult {
            pair_idx: 0,
            source: 0,
            target: 1,
            lineage: 0,
            bin: 2,
            score: 0.1,
            n_eff_cells: 10.0,
            weight_mass_cells: 10.0,
            lineage_entropy: 0.3,
        }];
        let mut fake = Vec::new();
        for i in 0..9 {
            fake.push(crate::ccc::LrBinResult {
                pair_idx: 0,
                source: 0,
                target: 1,
                lineage: 0,
                bin: 2,
                score: 0.2 + i as f32 * 0.05,
                n_eff_cells: 10.0,
                weight_mass_cells: 10.0,
                lineage_entropy: 0.3,
            });
        }
        let pf = compute_p_fake_by_position(&real, &fake, 9);
        assert!((pf[0] - 1.0).abs() < 1e-6, "p_fake = {}", pf[0]);
    }

    #[test]
    fn compute_p_fake_ignores_different_position() {
        // Real score at (source=0, target=1, bin=2), but fake scores
        // are at different positions → no matched fake support.
        let real = vec![crate::ccc::LrBinResult {
            pair_idx: 0,
            source: 0,
            target: 1,
            lineage: 0,
            bin: 2,
            score: 0.5,
            n_eff_cells: 10.0,
            weight_mass_cells: 10.0,
            lineage_entropy: 0.3,
        }];
        let fake = vec![crate::ccc::LrBinResult {
            pair_idx: 0,
            source: 1,
            target: 2, // different source/target
            lineage: 0,
            bin: 2,
            score: 0.99,
            n_eff_cells: 10.0,
            weight_mass_cells: 10.0,
            lineage_entropy: 0.3,
        }];
        let pf = compute_p_fake_by_position(&real, &fake, 1);
        assert!(
            pf[0].is_nan(),
            "different position should not define p_fake"
        );
    }

    // ─── Feature normalisation ───────────────────────────────────────

    #[test]
    fn normalized_distance_same_features() {
        let a = [0.5, 0.3, 0.2, 0.1];
        let mins = &[0.0, 0.0, 0.0, 0.0];
        let maxs = &[1.0, 1.0, 1.0, 1.0];
        assert!((norm_euclid_sq(&a, &a, mins, maxs)).abs() < 1e-6);
    }

    #[test]
    fn normalized_distance_opposite_corners() {
        let a = [0.0, 0.0, 0.0, 0.0];
        let b = [1.0, 1.0, 1.0, 1.0];
        let mins = &[0.0, 0.0, 0.0, 0.0];
        let maxs = &[1.0, 1.0, 1.0, 1.0];
        let d2 = norm_euclid_sq(&a, &b, mins, maxs);
        assert!((d2 - 4.0).abs() < 1e-6, "expected 4.0, got {}", d2);
    }

    #[test]
    fn complex_lr_pair_cardinality_gate() {
        use crate::ccc::LrPair;
        let ok = LrPair {
            ligand_rows: vec![0, 1],
            receptor_rows: vec![2, 3],
        };
        assert!(is_complex_lr_pair(&ok));
        assert!(complex_cardinality_supported(&ok));
        let three_subunits = LrPair {
            ligand_rows: vec![0, 1, 2],
            receptor_rows: vec![3],
        };
        assert!(is_complex_lr_pair(&three_subunits));
        assert!(complex_cardinality_supported(&three_subunits));
        let too_many = LrPair {
            ligand_rows: (0..=MAX_COMPLEX_CARDINALITY).collect(),
            receptor_rows: vec![3],
        };
        assert!(is_complex_lr_pair(&too_many));
        assert!(!complex_cardinality_supported(&too_many));
    }

    #[test]
    fn fake_null_mode_supports_complex_all_layers() {
        assert!(fake_null_mode_supports_complex(
            FakeNullMode::ContextGeneric
        ));
        assert!(fake_null_mode_supports_complex(FakeNullMode::Both));
        assert!(fake_null_mode_supports_complex(FakeNullMode::All));
        assert!(fake_null_mode_supports_complex(FakeNullMode::EdgeRewire));
        assert!(fake_null_mode_supports_complex(
            FakeNullMode::LrNetworkRewire
        ));
    }

    #[test]
    fn complex_fake_pair_draw_key_hashes_subunit_sets() {
        let fp = FakeLrPair {
            ligand_gene: 1,
            receptor_gene: 4,
            complex_ligand_rows: Some(vec![1, 2]),
            complex_receptor_rows: Some(vec![3, 4]),
            shadowing_pair_idx: 0,
            source: 0,
            target: 1,
            lineage: 0,
            null_kind: FAKE_NULL_GENERIC,
            ligand_match_distance: 0.1,
            receptor_match_distance: 0.2,
        };
        assert!(fp.is_complex());
        assert_eq!(fp.ligand_rows(), &[1, 2]);
        assert_eq!(fp.receptor_rows(), &[3, 4]);
        let k1 = fp.draw_key();
        let k2 = fp.draw_key();
        assert_eq!(k1, k2);
        assert_ne!(k1, (1, 4));
    }
}
