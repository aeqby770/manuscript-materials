//! Temporal permutation null for pseudotime-binned LR trajectories.
//!
//! The null fixes expression, cell type labels, and lineage weights, then
//! permutes pseudotime within each `(cell type, dominant lineage)` stratum.
//! Trajectory statistics use the same structural-valid bin rules as real
//! scores and FakeLR null curves.

use rand::seq::SliceRandom;
use rand::{rngs::StdRng, SeedableRng};
use std::collections::HashMap;

use crate::ccc::lr_scoring::{
    compute_bin_bounds, compute_lr_scores_bins_fixed_scale_grid_with_bounds,
    compute_lr_scores_bins_with_scale_grid, LrBinResult, LrPair,
};
use crate::ccc::trajectory_curve::{
    build_lr_trajectory_curve_from_grid, evaluate_trajectory_curve, trajectory_keys_from_scores,
    TrajKey, TrajectoryCurve, TrajectoryEvalConfig,
};
use crate::core::SparseMatrix;

#[derive(Debug, Clone)]
pub struct TemporalPermutationTrajectory {
    pub pair_idx: usize,
    pub source: usize,
    pub target: usize,
    pub lineage: usize,
    pub statistic: f32,
    pub pvalue: f64,
    pub n_eff_cells: f32,
    pub lineage_entropy: f32,
    pub n_valid_bins: usize,
    pub valid_bin_fraction: f32,
    pub n_expression_supported_bins: usize,
    pub n_eff_at_max: f32,
    pub n_eff_at_min: f32,
    pub valid_at_max: bool,
    pub valid_at_min: bool,
    pub trajectory_evaluable: bool,
    pub n_perm_expected: usize,
    pub n_perm_observed: usize,
    pub perm_coverage: f64,
    /// Median trajectory statistic across evaluable permutation draws.
    pub null_stat_median: f32,
}

#[derive(Debug, Clone)]
pub struct TemporalPermutationBin {
    pub pair_idx: usize,
    pub source: usize,
    pub target: usize,
    pub lineage: usize,
    pub bin: usize,
    pub score: f32,
    pub deviation: f32,
    pub pvalue: f64,
    pub structural_valid: bool,
    pub expression_supported: bool,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PermutationStratumMode {
    /// `(block, cell_type, dominant_lineage)` — default.
    BlockCelltypeDominantLineage,
    /// `(block, cell_type)` — use when soft lineage weights assign mass to
    /// multiple lineages per cell; one pseudotime shuffle per cell.
    BlockCelltype,
}

fn dominant_lineage(lineage_weights: &[f32], cell: usize, n_lineages: usize) -> Option<usize> {
    let row = &lineage_weights[cell * n_lineages..cell * n_lineages + n_lineages];
    let mut best: Option<(usize, f32)> = None;
    for (lineage, &weight) in row.iter().enumerate() {
        if weight.is_finite() && best.map_or(true, |(_, current)| weight > current) {
            best = Some((lineage, weight));
        }
    }
    best.filter(|(_, weight)| *weight > 1e-8)
        .map(|(lineage, _)| lineage)
}

fn shuffle_pseudotime_within_strata(
    pseudotime: &[f32],
    cell_types: &[usize],
    n_types: usize,
    lineage_weights: &[f32],
    n_lineages: usize,
    block_labels: &[usize],
    mode: PermutationStratumMode,
    seed: u64,
) -> Vec<f32> {
    let mut groups: HashMap<(usize, usize, usize), Vec<usize>> = HashMap::new();
    for cell in 0..pseudotime.len() {
        if !pseudotime[cell].is_finite() {
            continue;
        }
        let ct = cell_types[cell];
        if ct >= n_types {
            continue;
        }
        let stratum = match mode {
            PermutationStratumMode::BlockCelltypeDominantLineage => {
                let Some(lin) = dominant_lineage(lineage_weights, cell, n_lineages) else {
                    continue;
                };
                (block_labels[cell], ct, lin)
            }
            PermutationStratumMode::BlockCelltype => (block_labels[cell], ct, 0usize),
        };
        groups.entry(stratum).or_default().push(cell);
    }

    let mut out = pseudotime.to_vec();
    let mut rng = StdRng::seed_from_u64(seed);
    let mut strata: Vec<(usize, usize, usize)> = groups.keys().copied().collect();
    strata.sort_unstable();
    for stratum in strata {
        let cells = &groups[&stratum];
        let mut values: Vec<f32> = cells.iter().map(|&i| pseudotime[i]).collect();
        values.shuffle(&mut rng);
        for (cell, value) in cells.iter().zip(values.iter()) {
            out[*cell] = *value;
        }
    }
    out
}

fn deviations_valid(curve: &TrajectoryCurve) -> Vec<f32> {
    let valid_scores: Vec<f32> = curve
        .bins
        .iter()
        .filter(|b| b.structural_valid)
        .map(|b| b.score)
        .collect();
    if valid_scores.is_empty() {
        return vec![f32::NAN; curve.bins.len()];
    }
    let mean = valid_scores.iter().sum::<f32>() / valid_scores.len() as f32;
    curve
        .bins
        .iter()
        .map(|b| {
            if b.structural_valid {
                (b.score - mean).abs()
            } else {
                f32::NAN
            }
        })
        .collect()
}

fn update_max_t_bin_counts(
    real_curve: &TrajectoryCurve,
    real_deviation: &[f32],
    perm_deviation: &[f32],
    counts: &mut [u32],
    observed: &mut [u32],
) {
    let perm_max = perm_deviation
        .iter()
        .copied()
        .filter(|value| value.is_finite())
        .fold(f32::NEG_INFINITY, f32::max);
    if !perm_max.is_finite() {
        return;
    }
    for bin in 0..real_curve.bins.len() {
        if real_curve.bins[bin].structural_valid && real_deviation[bin].is_finite() {
            observed[bin] += 1;
            if perm_max >= real_deviation[bin] {
                counts[bin] += 1;
            }
        }
    }
}

fn peak_quality(scores: &[LrBinResult], key: TrajKey) -> (f32, f32) {
    scores
        .iter()
        .filter(|r| (r.pair_idx, r.source, r.target, r.lineage) == key)
        .max_by(|a, b| {
            a.score
                .partial_cmp(&b.score)
                .unwrap_or(std::cmp::Ordering::Equal)
        })
        .map(|r| (r.n_eff_cells, r.lineage_entropy))
        .unwrap_or((0.0, 1.0))
}

pub fn compute_temporal_permutation(
    expr: &SparseMatrix,
    cell_types: &[usize],
    n_types: usize,
    lr_pairs: &[LrPair],
    pseudotime: &[f32],
    lineage_weights: &[f32],
    n_lineages: usize,
    n_bins: usize,
    min_pct: f32,
    min_n_eff: f32,
    min_valid_bins: usize,
    min_valid_fraction: f32,
    n_perms: usize,
    seed: u64,
    stat_code: usize,
    top_k: usize,
) -> (
    Vec<TemporalPermutationTrajectory>,
    Vec<TemporalPermutationBin>,
) {
    let block_labels = vec![0usize; pseudotime.len()];
    let (trajectories, bins, _) = compute_temporal_permutation_with_scale(
        expr,
        cell_types,
        n_types,
        lr_pairs,
        pseudotime,
        lineage_weights,
        n_lineages,
        n_bins,
        min_pct,
        min_n_eff,
        min_valid_bins,
        min_valid_fraction,
        n_perms,
        seed,
        stat_code,
        top_k,
        &block_labels,
        PermutationStratumMode::BlockCelltypeDominantLineage,
        None,
    );
    (trajectories, bins)
}

pub(crate) fn compute_temporal_permutation_with_scale(
    expr: &SparseMatrix,
    cell_types: &[usize],
    n_types: usize,
    lr_pairs: &[LrPair],
    pseudotime: &[f32],
    lineage_weights: &[f32],
    n_lineages: usize,
    n_bins: usize,
    min_pct: f32,
    min_n_eff: f32,
    min_valid_bins: usize,
    min_valid_fraction: f32,
    n_perms: usize,
    seed: u64,
    stat_code: usize,
    top_k: usize,
    block_labels: &[usize],
    stratum_mode: PermutationStratumMode,
    requested_keys: Option<Vec<TrajKey>>,
) -> (
    Vec<TemporalPermutationTrajectory>,
    Vec<TemporalPermutationBin>,
    f32,
) {
    let cfg = TrajectoryEvalConfig {
        min_pct,
        min_n_eff,
        min_valid_bins,
        min_valid_fraction,
    };

    let (real_scores, _, hill_s, real_grid) = compute_lr_scores_bins_with_scale_grid(
        expr,
        cell_types,
        n_types,
        lr_pairs,
        pseudotime,
        lineage_weights,
        n_lineages,
        n_bins,
        min_pct,
    );
    let bin_bounds = compute_bin_bounds(
        pseudotime.len(),
        pseudotime,
        lineage_weights,
        n_lineages,
        n_bins,
    );

    let keys = if let Some(req) = requested_keys {
        if req.is_empty() {
            trajectory_keys_from_scores(&real_scores)
        } else {
            req
        }
    } else {
        trajectory_keys_from_scores(&real_scores)
    };
    let real_curves: Vec<TrajectoryCurve> = keys
        .iter()
        .map(|key| {
            build_lr_trajectory_curve_from_grid(
                &lr_pairs[key.0],
                &real_grid,
                key.1,
                key.2,
                key.3,
                min_pct,
                min_n_eff,
            )
        })
        .collect();
    let real_eval: Vec<_> = real_curves
        .iter()
        .map(|curve| evaluate_trajectory_curve(curve, stat_code, top_k, &cfg))
        .collect();
    let real_dev: Vec<Vec<f32>> = real_curves.iter().map(deviations_valid).collect();
    let mut traj_counts = vec![0u32; keys.len()];
    let mut traj_perm_observed = vec![0u32; keys.len()];
    let mut bin_counts = vec![vec![0u32; n_bins]; keys.len()];
    let mut bin_perm_observed = vec![vec![0u32; n_bins]; keys.len()];
    let mut traj_perm_stats: Vec<Vec<f32>> = (0..keys.len()).map(|_| Vec::new()).collect();

    for perm in 0..n_perms {
        let pt_perm = shuffle_pseudotime_within_strata(
            pseudotime,
            cell_types,
            n_types,
            lineage_weights,
            n_lineages,
            block_labels,
            stratum_mode,
            seed.wrapping_add(perm as u64 + 1),
        );
        let (_, _, perm_grid) = compute_lr_scores_bins_fixed_scale_grid_with_bounds(
            expr,
            cell_types,
            n_types,
            lr_pairs,
            &pt_perm,
            lineage_weights,
            n_lineages,
            n_bins,
            min_pct,
            hill_s,
            Some(&bin_bounds),
        );
        for (i, key) in keys.iter().enumerate() {
            if !real_eval[i].evaluable {
                continue;
            }
            let perm_curve = build_lr_trajectory_curve_from_grid(
                &lr_pairs[key.0],
                &perm_grid,
                key.1,
                key.2,
                key.3,
                min_pct,
                min_n_eff,
            );
            let perm_eval = evaluate_trajectory_curve(&perm_curve, stat_code, top_k, &cfg);
            if perm_eval.evaluable {
                traj_perm_observed[i] += 1;
                traj_perm_stats[i].push(perm_eval.statistic);
                if perm_eval.statistic >= real_eval[i].statistic {
                    traj_counts[i] += 1;
                }
                let perm_dev = deviations_valid(&perm_curve);
                if stat_code == 0 {
                    update_max_t_bin_counts(
                        &real_curves[i],
                        &real_dev[i],
                        &perm_dev,
                        &mut bin_counts[i],
                        &mut bin_perm_observed[i],
                    );
                }
            }
        }
    }

    let trajectories: Vec<TemporalPermutationTrajectory> = keys
        .iter()
        .enumerate()
        .map(|(i, key)| {
            let ev = &real_eval[i];
            let (n_eff_cells, lineage_entropy) = peak_quality(&real_scores, *key);
            let statistic = if ev.evaluable { ev.statistic } else { f32::NAN };
            let coverage = if n_perms > 0 {
                traj_perm_observed[i] as f64 / n_perms as f64
            } else {
                f64::NAN
            };
            let pvalue = if ev.evaluable && traj_perm_observed[i] > 0 {
                (traj_counts[i] as f64 + 1.0) / (traj_perm_observed[i] as f64 + 1.0)
            } else {
                f64::NAN
            };
            let null_stat_median = if traj_perm_stats[i].is_empty() {
                f32::NAN
            } else {
                let mut stats = traj_perm_stats[i].clone();
                stats.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
                let m = stats.len() / 2;
                if stats.len() % 2 == 1 {
                    stats[m]
                } else {
                    (stats[m - 1] + stats[m]) / 2.0
                }
            };
            TemporalPermutationTrajectory {
                pair_idx: key.0,
                source: key.1,
                target: key.2,
                lineage: key.3,
                statistic,
                pvalue,
                n_eff_cells,
                lineage_entropy,
                n_valid_bins: ev.n_valid_bins,
                valid_bin_fraction: ev.valid_bin_fraction,
                n_expression_supported_bins: ev.n_expression_supported_bins,
                n_eff_at_max: ev.n_eff_at_max,
                n_eff_at_min: ev.n_eff_at_min,
                valid_at_max: ev.valid_at_max,
                valid_at_min: ev.valid_at_min,
                trajectory_evaluable: ev.evaluable,
                n_perm_expected: n_perms,
                n_perm_observed: traj_perm_observed[i] as usize,
                perm_coverage: coverage,
                null_stat_median,
            }
        })
        .collect();

    let mut bins: Vec<TemporalPermutationBin> = Vec::with_capacity(keys.len() * n_bins);
    for (i, key) in keys.iter().enumerate() {
        let curve = &real_curves[i];
        for bin in 0..n_bins {
            let cell = &curve.bins[bin];
            let pvalue = if stat_code == 0
                && cell.structural_valid
                && real_dev[i][bin].is_finite()
                && bin_perm_observed[i][bin] > 0
            {
                (bin_counts[i][bin] as f64 + 1.0) / (bin_perm_observed[i][bin] as f64 + 1.0)
            } else {
                f64::NAN
            };
            bins.push(TemporalPermutationBin {
                pair_idx: key.0,
                source: key.1,
                target: key.2,
                lineage: key.3,
                bin,
                score: if cell.structural_valid {
                    cell.score
                } else {
                    f32::NAN
                },
                deviation: real_dev[i][bin],
                pvalue,
                structural_valid: cell.structural_valid,
                expression_supported: cell.expression_supported,
            });
        }
    }

    (trajectories, bins, hill_s)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ccc::trajectory_curve::{BinCell, TrajectoryCurve};

    #[test]
    fn deviations_respect_structural_valid_mask() {
        let curve = TrajectoryCurve {
            bins: vec![
                BinCell::invalid(),
                BinCell {
                    score: 1.0,
                    structural_valid: true,
                    expression_supported: true,
                    n_eff_source: 10.0,
                    n_eff_target: 10.0,
                },
                BinCell {
                    score: 3.0,
                    structural_valid: true,
                    expression_supported: true,
                    n_eff_source: 10.0,
                    n_eff_target: 10.0,
                },
            ],
        };
        let d = deviations_valid(&curve);
        assert!(d[0].is_nan());
        assert!((d[1] - 1.0).abs() < 1e-6);
        assert!((d[2] - 1.0).abs() < 1e-6);
    }

    #[test]
    fn max_t_uses_largest_permutation_deviation_for_every_real_bin() {
        let real_curve = TrajectoryCurve {
            bins: vec![
                BinCell {
                    score: 0.0,
                    structural_valid: true,
                    expression_supported: true,
                    n_eff_source: 10.0,
                    n_eff_target: 10.0,
                },
                BinCell {
                    score: 0.0,
                    structural_valid: true,
                    expression_supported: true,
                    n_eff_source: 10.0,
                    n_eff_target: 10.0,
                },
            ],
        };
        let real_deviation = vec![0.7, 0.2];
        let perm_deviation = vec![0.1, 0.8];
        let mut counts = vec![0u32; 2];
        let mut observed = vec![0u32; 2];

        update_max_t_bin_counts(
            &real_curve,
            &real_deviation,
            &perm_deviation,
            &mut counts,
            &mut observed,
        );

        assert_eq!(observed, vec![1, 1]);
        assert_eq!(counts, vec![1, 1]);
    }

    #[test]
    fn temporal_shuffle_is_reproducible_and_excludes_zero_weight_cells() {
        let pseudotime = vec![0.1, 0.2, 0.9];
        let cell_types = vec![0, 0, 0];
        let lineage_weights = vec![1.0, 0.0, 1.0, 0.0, 0.0, 0.0];
        let blocks = vec![0, 0, 0];
        let a = shuffle_pseudotime_within_strata(
            &pseudotime,
            &cell_types,
            1,
            &lineage_weights,
            2,
            &blocks,
            PermutationStratumMode::BlockCelltypeDominantLineage,
            17,
        );
        let b = shuffle_pseudotime_within_strata(
            &pseudotime,
            &cell_types,
            1,
            &lineage_weights,
            2,
            &blocks,
            PermutationStratumMode::BlockCelltypeDominantLineage,
            17,
        );
        assert_eq!(a, b);
        assert_eq!(a[2], 0.9);
        let mut active = a[..2].to_vec();
        active.sort_by(|x, y| x.partial_cmp(y).unwrap());
        assert_eq!(active, vec![0.1, 0.2]);
    }

    #[test]
    fn temporal_shuffle_never_crosses_blocks() {
        let pseudotime = vec![0.1, 0.2, 0.8, 0.9];
        let cell_types = vec![0, 0, 0, 0];
        let lineage_weights = vec![1.0; 4];
        let blocks = vec![0, 0, 1, 1];
        let shuffled = shuffle_pseudotime_within_strata(
            &pseudotime,
            &cell_types,
            1,
            &lineage_weights,
            1,
            &blocks,
            PermutationStratumMode::BlockCelltypeDominantLineage,
            29,
        );
        let mut first = shuffled[..2].to_vec();
        let mut second = shuffled[2..].to_vec();
        first.sort_by(|a, b| a.partial_cmp(b).unwrap());
        second.sort_by(|a, b| a.partial_cmp(b).unwrap());
        assert_eq!(first, vec![0.1, 0.2]);
        assert_eq!(second, vec![0.8, 0.9]);
    }
}
