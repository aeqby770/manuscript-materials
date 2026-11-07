//! Evidence decomposition via aligned profile correlation.
//!
//! Within each `(source, target, lineage, bin)` stratum, source and target
//! cells are split into shared pseudotime micro-bins. Ligand and receptor
//! complex means are computed per micro-bin and correlated across the aligned
//! profiles. This is descriptive stratum-level evidence, not partial
//! correlation or a formal conditional-independence test.  The profile
//! correlation is based on at most five micro-bins, so no formal Fisher
//! z significance is reported.

use crate::core::SparseMatrix;
use rayon::prelude::*;

/// Result of the evidence decomposition for a single (ligand, receptor)
/// within one (source, target, lineage, bin) stratum.
#[derive(Debug, Clone, Copy)]
pub struct EvidenceDecomp {
    /// Aligned profile-correlation coefficient r in [-1, 1].
    pub aligned_r: f32,
    /// Deprecated placeholder.  Always NaN; kept only for API compatibility.
    pub z_score: f32,
    /// Kish effective sample size for the source/target stratum.
    pub n_eff: f32,
    /// Number of aligned micro-bins used by `aligned_r`.
    pub n_microbins: u32,
    /// Whether a finite profile concordance could be computed.
    pub evidence_sufficient: bool,
}

/// Build `k` equal-frequency sub-bin means from `(pseudotime, value)` pairs.
#[inline]
fn subbin_means(pairs: &[(f32, f32)], k: usize) -> Vec<f32> {
    let kk = k.min(pairs.len());
    let n = pairs.len();
    let mut means = vec![0.0f32; kk];
    let mut counts = vec![0u32; kk];
    for (i, &(_pt, val)) in pairs.iter().enumerate() {
        let sub = (i * kk / n).min(kk - 1);
        means[sub] += val;
        counts[sub] += 1;
    }
    for i in 0..kk {
        if counts[i] > 0 {
            means[i] /= counts[i] as f32;
        }
    }
    means
}

#[inline]
fn shared_microbin_count(source_n: usize, target_n: usize) -> usize {
    5usize.min(source_n).min(target_n)
}

/// Compute Pearson correlation coefficient from two
/// equal-length vectors.
/// Returns `(r, n_microbins, sufficient)` where:
/// - `r` ∈ [-1, 1] is the Pearson correlation
/// - `sufficient` is true when r is finite and the profile has >= 3 bins
#[inline]
fn pearson_profile_r(a: &[f32], b: &[f32]) -> (f32, u32, bool) {
    let k = a.len().min(b.len());
    if k < 3 {
        return (f32::NAN, k as u32, false);
    }
    let mean_a = a.iter().sum::<f32>() / k as f32;
    let mean_b = b.iter().sum::<f32>() / k as f32;
    let mut cov = 0.0f64;
    let mut var_a = 0.0f64;
    let mut var_b = 0.0f64;
    for i in 0..k {
        let da = a[i] as f64 - mean_a as f64;
        let db = b[i] as f64 - mean_b as f64;
        cov += da * db;
        var_a += da * da;
        var_b += db * db;
    }
    let denom = (var_a * var_b).sqrt();
    if denom <= 1e-12 {
        return (f32::NAN, k as u32, false);
    }
    let r = (cov / denom).clamp(-1.0, 1.0) as f32;
    (r, k as u32, r.is_finite())
}

/// Compute per-stratum aligned profile correlations for LR score rows.
///
/// For each real `LrBinResult` row in `real_scores`, this function:
/// 1. Identifies cells in the relevant `(source, target, lineage, bin)`
///    strata (using `pseudotime`, `lineage_weights`, `cell_types`).
/// 2. Extracts the ligand and receptor expression values for those cells.
/// 3. Uses one shared micro-bin count for source and target profiles.
/// 4. Computes a descriptive Pearson profile correlation; no p-value or
///    Fisher z significance is assigned.
///
/// Returns a vector of `EvidenceDecomp` parallel to `real_scores`.
///
/// # Arguments
/// * `real_scores` — output rows from `compute_lr_scores_bins`.
/// * `expr` — the CSC expression matrix.
/// * `cell_types` — per-cell cluster labels (0-based).
/// * `pseudotime` — per-cell pseudotime ∈ [0,1] or NaN.
/// * `lineage_weights` — flat [n_cells × n_lineages], row-major.
/// * `n_lineages`, `n_bins` — tensor dimensions.
/// * `ligand_rows_per_pair` — global row indices (0-based) for each LR
///   pair's ligand subunit(s).
/// * `receptor_rows_per_pair` — same for receptors.
/// * `n_threads` — rayon thread count (pass 0 for default).
pub fn evidence_decomp_batch(
    real_scores: &[crate::ccc::LrBinResult],
    expr: &SparseMatrix,
    cell_types: &[usize],
    n_types: usize,
    pseudotime: &[f32],
    lineage_weights: &[f32],
    n_lineages: usize,
    n_bins: usize,
    ligand_rows_per_pair: &[Vec<usize>],
    receptor_rows_per_pair: &[Vec<usize>],
    min_pct: f32,
) -> Vec<EvidenceDecomp> {
    if real_scores.is_empty() {
        return Vec::new();
    }

    let n_cells = expr.ncols;

    // ── Pre-compute bin boundaries (reuses lr_scoring::compute_bin_bounds) ─
    let bin_bounds = crate::ccc::lr_scoring::compute_bin_bounds(
        n_cells,
        pseudotime,
        lineage_weights,
        n_lineages,
        n_bins,
    );

    // Row-parallel: each score row is independent
    real_scores
        .par_iter()
        .map(|row| {
            let pair_idx = row.pair_idx;
            let lin = row.lineage;
            let b = row.bin;
            let s_ct = row.source;
            let t_ct = row.target;

            // ── Collect (pseudotime, expression) pairs for source and target
            //     cells in this (lineage, bin) ───────────────────────────────
            let lig_genes = &ligand_rows_per_pair[pair_idx];
            let rec_genes = &receptor_rows_per_pair[pair_idx];

            let mut src_pairs: Vec<(f32, f32)> = Vec::new();
            let mut tgt_pairs: Vec<(f32, f32)> = Vec::new();

            for cell in 0..n_cells {
                let ct = cell_types[cell];
                if ct >= n_types {
                    continue;
                }
                let pt_raw = pseudotime[cell];
                if pt_raw.is_nan() {
                    continue;
                }
                let pt = pt_raw.clamp(0.0, 1.0);
                let w = lineage_weights[cell * n_lineages + lin];
                if w < 1e-8 {
                    continue;
                }
                let cell_bin = bin_bounds[lin]
                    .partition_point(|&bnd| bnd <= pt)
                    .saturating_sub(1)
                    .min(n_bins - 1);
                if cell_bin != b {
                    continue;
                }

                if ct == s_ct {
                    let mut lig_sum = 0.0f32;
                    let mut lig_n = 0usize;
                    for nnz in expr.col_ptrs[cell]..expr.col_ptrs[cell + 1] {
                        if lig_genes.contains(&expr.rows[nnz]) {
                            lig_sum += expr.values[nnz];
                            lig_n += 1;
                        }
                    }
                    if lig_n > 0 && lig_n as f32 / lig_genes.len().max(1) as f32 >= min_pct {
                        src_pairs.push((pt, lig_sum / lig_n as f32));
                    }
                }

                if ct == t_ct {
                    let mut rec_sum = 0.0f32;
                    let mut rec_n = 0usize;
                    for nnz in expr.col_ptrs[cell]..expr.col_ptrs[cell + 1] {
                        if rec_genes.contains(&expr.rows[nnz]) {
                            rec_sum += expr.values[nnz];
                            rec_n += 1;
                        }
                    }
                    if rec_n > 0 && rec_n as f32 / rec_genes.len().max(1) as f32 >= min_pct {
                        tgt_pairs.push((pt, rec_sum / rec_n as f32));
                    }
                }
            }

            let n_eff = row.n_eff_cells;

            if src_pairs.len() < 3 || tgt_pairs.len() < 3 || n_eff < 5.0 {
                return EvidenceDecomp {
                    aligned_r: f32::NAN,
                    z_score: f32::NAN,
                    n_eff,
                    n_microbins: 0,
                    evidence_sufficient: false,
                };
            }

            // Sort by pseudotime for temporal sub-bin grouping
            src_pairs.sort_unstable_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
            tgt_pairs.sort_unstable_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

            // Build sub-bin means and correlate aligned profiles
            let kk = shared_microbin_count(src_pairs.len(), tgt_pairs.len());
            let src_means = subbin_means(&src_pairs, kk);
            let tgt_means = subbin_means(&tgt_pairs, kk);

            let (r, n_microbins, sufficient) = pearson_profile_r(&src_means, &tgt_means);

            EvidenceDecomp {
                aligned_r: r,
                z_score: f32::NAN,
                n_eff,
                n_microbins,
                evidence_sufficient: sufficient,
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn evidence_decomp_empty_input() {
        let empty: Vec<crate::ccc::LrBinResult> = vec![];
        let result = evidence_decomp_batch(
            &empty,
            &SparseMatrix::new_csc(vec![], vec![0], vec![], 0, 0),
            &[],
            0,
            &[],
            &[],
            0,
            0,
            &[],
            &[],
            0.0,
        );
        assert!(result.is_empty());
    }

    #[test]
    fn evidence_decomp_default_structure() {
        let d = EvidenceDecomp {
            aligned_r: 0.5,
            z_score: f32::NAN,
            n_eff: 25.0,
            n_microbins: 5,
            evidence_sufficient: true,
        };
        assert!(d.evidence_sufficient);
        assert!((d.aligned_r - 0.5).abs() < 1e-6);
        assert!(d.z_score.is_nan());
        assert_eq!(d.n_microbins, 5);
    }

    #[test]
    fn shared_microbins_use_the_smaller_profile() {
        assert_eq!(shared_microbin_count(8, 3), 3);
        assert_eq!(shared_microbin_count(4, 9), 4);
        assert_eq!(shared_microbin_count(8, 9), 5);
    }

    #[test]
    fn constant_profiles_have_no_defined_concordance() {
        let (r, k, sufficient) = pearson_profile_r(&[1.0, 1.0, 1.0], &[2.0, 3.0, 4.0]);
        assert!(r.is_nan());
        assert_eq!(k, 3);
        assert!(!sufficient);
    }
}
