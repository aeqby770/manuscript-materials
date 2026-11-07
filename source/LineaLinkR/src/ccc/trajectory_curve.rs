//! Shared trajectory-curve evaluation with structural vs expression bin semantics.
//!
//! Rules (applied identically to real, FakeLR, and temporal-permutation curves):
//! - `structural_valid = FALSE` → score is NaN, excluded from trajectory statistics.
//! - `structural_valid = TRUE` and expression below `min_pct` → score = 0 (detectable absence).
//! - `structural_valid = TRUE` and expression supported → normal Hill score.
//! - Amplitude (default stat) requires ≥2 structural-valid bins; max/min must both be structural-valid.

use crate::ccc::fake_lr::GeneBinStats;
use crate::ccc::lr_scoring::{LrBinGrid, LrPair};

pub type TrajKey = (usize, usize, usize, usize);

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BinCell {
    pub score: f32,
    pub structural_valid: bool,
    pub expression_supported: bool,
    pub n_eff_source: f32,
    pub n_eff_target: f32,
}

impl BinCell {
    #[inline]
    pub fn invalid() -> Self {
        Self {
            score: f32::NAN,
            structural_valid: false,
            expression_supported: false,
            n_eff_source: 0.0,
            n_eff_target: 0.0,
        }
    }
}

#[derive(Debug, Clone)]
pub struct TrajectoryCurve {
    pub bins: Vec<BinCell>,
}

#[derive(Debug, Clone, Copy)]
pub struct TrajectoryEvalConfig {
    pub min_pct: f32,
    pub min_n_eff: f32,
    pub min_valid_bins: usize,
    pub min_valid_fraction: f32,
}

impl Default for TrajectoryEvalConfig {
    fn default() -> Self {
        Self {
            min_pct: 0.05,
            min_n_eff: 5.0,
            min_valid_bins: 3,
            min_valid_fraction: 0.25,
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct TrajectoryEvalSummary {
    pub statistic: f32,
    pub evaluable: bool,
    pub n_valid_bins: usize,
    pub valid_bin_fraction: f32,
    pub n_expression_supported_bins: usize,
    pub n_eff_at_max: f32,
    pub n_eff_at_min: f32,
    pub valid_at_max: bool,
    pub valid_at_min: bool,
}

#[inline]
pub fn hill_score(lm: f32, rm: f32, hill_s: f32) -> f32 {
    if lm > 0.0 && rm > 0.0 {
        let geo = (lm * rm).sqrt();
        geo / (hill_s + geo)
    } else {
        0.0
    }
}

#[inline]
fn kish_ess(weight_sum: f64, squared_weight_sum: f64) -> f32 {
    if squared_weight_sum > 1e-12 {
        (weight_sum * weight_sum / squared_weight_sum) as f32
    } else {
        0.0
    }
}

#[inline]
pub fn bin_cell_from_support(
    ws: f32,
    wt: f32,
    lp: f32,
    rp: f32,
    lm: f32,
    rm: f32,
    min_pct: f32,
    min_n_eff: f32,
    hill_s: f32,
) -> BinCell {
    let structural_valid = ws >= min_n_eff && wt >= min_n_eff;
    if !structural_valid {
        return BinCell::invalid();
    }
    let expression_supported = lp >= min_pct && rp >= min_pct;
    let score = if expression_supported {
        hill_score(lm, rm, hill_s)
    } else {
        0.0
    };
    BinCell {
        score,
        structural_valid: true,
        expression_supported,
        n_eff_source: ws,
        n_eff_target: wt,
    }
}

pub fn build_fake_trajectory_curve(
    lig_gene: usize,
    rec_gene: usize,
    source: usize,
    target: usize,
    lineage: usize,
    real_stats: &GeneBinStats,
    fake_stats: &GeneBinStats,
    structural: &crate::ccc::lr_scoring::LrBinGrid,
    min_pct: f32,
    min_n_eff: f32,
) -> TrajectoryCurve {
    let n_bins = structural.n_bins;
    let slb = structural.slb();
    let sb = structural.sb();
    let hill_s = structural.hill_s;
    let mut bins = Vec::with_capacity(n_bins);
    for bin in 0..n_bins {
        let source_idx = source * slb + lineage * sb + bin;
        let target_idx = target * slb + lineage * sb + bin;
        let ws = kish_ess(structural.wsums[source_idx], structural.w2sums[source_idx]);
        let wt = kish_ess(structural.wsums[target_idx], structural.w2sums[target_idx]);
        let lp =
            GeneBinStats::get_pct_resolved(real_stats, fake_stats, lig_gene, source, lineage, bin);
        let rp =
            GeneBinStats::get_pct_resolved(real_stats, fake_stats, rec_gene, target, lineage, bin);
        let lm =
            GeneBinStats::get_mean_resolved(real_stats, fake_stats, lig_gene, source, lineage, bin);
        let rm =
            GeneBinStats::get_mean_resolved(real_stats, fake_stats, rec_gene, target, lineage, bin);
        bins.push(bin_cell_from_support(
            ws, wt, lp, rp, lm, rm, min_pct, min_n_eff, hill_s,
        ));
    }
    TrajectoryCurve { bins }
}

/// Fake complex LR curve: min-aggregate ligand/receptor subunits (bottleneck).
pub fn build_fake_complex_trajectory_curve(
    ligand_rows: &[usize],
    receptor_rows: &[usize],
    source: usize,
    target: usize,
    lineage: usize,
    real_stats: &GeneBinStats,
    fake_stats: &GeneBinStats,
    structural: &crate::ccc::lr_scoring::LrBinGrid,
    min_pct: f32,
    min_n_eff: f32,
) -> TrajectoryCurve {
    let n_bins = structural.n_bins;
    let slb = structural.slb();
    let sb = structural.sb();
    let hill_s = structural.hill_s;
    let mut bins = Vec::with_capacity(n_bins);
    for bin in 0..n_bins {
        let source_idx = source * slb + lineage * sb + bin;
        let target_idx = target * slb + lineage * sb + bin;
        let ws = kish_ess(structural.wsums[source_idx], structural.w2sums[source_idx]);
        let wt = kish_ess(structural.wsums[target_idx], structural.w2sums[target_idx]);
        let lm = ligand_rows
            .iter()
            .map(|&g| {
                GeneBinStats::get_mean_resolved(real_stats, fake_stats, g, source, lineage, bin)
            })
            .fold(f32::INFINITY, f32::min);
        let lp = ligand_rows
            .iter()
            .map(|&g| {
                GeneBinStats::get_pct_resolved(real_stats, fake_stats, g, source, lineage, bin)
            })
            .fold(f32::INFINITY, f32::min);
        let rm = receptor_rows
            .iter()
            .map(|&g| {
                GeneBinStats::get_mean_resolved(real_stats, fake_stats, g, target, lineage, bin)
            })
            .fold(f32::INFINITY, f32::min);
        let rp = receptor_rows
            .iter()
            .map(|&g| {
                GeneBinStats::get_pct_resolved(real_stats, fake_stats, g, target, lineage, bin)
            })
            .fold(f32::INFINITY, f32::min);
        bins.push(bin_cell_from_support(
            ws, wt, lp, rp, lm, rm, min_pct, min_n_eff, hill_s,
        ));
    }
    TrajectoryCurve { bins }
}

fn lig_rec_at_bin(
    lig_local: &[usize],
    rec_local: &[usize],
    grid: &LrBinGrid,
    source: usize,
    target: usize,
    lineage: usize,
    bin: usize,
) -> (f32, f32, f32, f32) {
    let sl = grid.sl();
    let slb = grid.slb();
    let sb = grid.sb();
    let get = |local: usize, ct: usize, lin: usize, b: usize| -> (f32, f32) {
        let i = local * sl + ct * slb + lin * sb + b;
        (grid.means[i], grid.pcts[i])
    };
    let lm = lig_local
        .iter()
        .map(|&r| get(r, source, lineage, bin).0)
        .fold(f32::INFINITY, f32::min);
    let lp = lig_local
        .iter()
        .map(|&r| get(r, source, lineage, bin).1)
        .fold(f32::INFINITY, f32::min);
    let rm = rec_local
        .iter()
        .map(|&r| get(r, target, lineage, bin).0)
        .fold(f32::INFINITY, f32::min);
    let rp = rec_local
        .iter()
        .map(|&r| get(r, target, lineage, bin).1)
        .fold(f32::INFINITY, f32::min);
    (lm, lp, rm, rp)
}

pub fn positive_score_from_grid(
    pair: &LrPair,
    grid: &LrBinGrid,
    source: usize,
    target: usize,
    lineage: usize,
    bin: usize,
) -> f32 {
    let lig_local = local_indices(&pair.ligand_rows, &grid.g2local);
    let rec_local = local_indices(&pair.receptor_rows, &grid.g2local);
    if lig_local.is_empty() || rec_local.is_empty() {
        return 0.0;
    }
    let sl = grid.sl();
    let slb = grid.slb();
    let sb = grid.sb();
    let conditional_mean = |rows: &[usize], cell_type: usize| -> f32 {
        rows.iter()
            .map(|&local| {
                let i = local * sl + cell_type * slb + lineage * sb + bin;
                if grid.pcts[i] > 0.0 {
                    grid.means[i] / grid.pcts[i]
                } else {
                    0.0
                }
            })
            .fold(f32::INFINITY, f32::min)
    };
    let ligand = conditional_mean(&lig_local, source);
    let receptor = conditional_mean(&rec_local, target);
    let geo = (ligand * receptor).sqrt();
    if geo > 0.0 { geo / (grid.hill_s + geo) } else { 0.0 }
}

fn local_indices(rows: &[usize], g2local: &[usize]) -> Vec<usize> {
    rows.iter()
        .filter_map(|&g| {
            if g < g2local.len() {
                let local = g2local[g];
                if local != usize::MAX {
                    Some(local)
                } else {
                    None
                }
            } else {
                None
            }
        })
        .collect()
}

pub fn build_lr_trajectory_curve_from_grid(
    pair: &LrPair,
    grid: &LrBinGrid,
    source: usize,
    target: usize,
    lineage: usize,
    min_pct: f32,
    min_n_eff: f32,
) -> TrajectoryCurve {
    if source >= grid.n_types || target >= grid.n_types || lineage >= grid.n_lineages {
        return TrajectoryCurve {
            bins: vec![BinCell::invalid(); grid.n_bins],
        };
    }
    let lig_local = local_indices(&pair.ligand_rows, &grid.g2local);
    let rec_local = local_indices(&pair.receptor_rows, &grid.g2local);
    let slb = grid.slb();
    let sb = grid.sb();
    let mut bins = Vec::with_capacity(grid.n_bins);
    for bin in 0..grid.n_bins {
        let source_idx = source * slb + lineage * sb + bin;
        let target_idx = target * slb + lineage * sb + bin;
        let ws = kish_ess(grid.wsums[source_idx], grid.w2sums[source_idx]);
        let wt = kish_ess(grid.wsums[target_idx], grid.w2sums[target_idx]);
        if lig_local.is_empty() || rec_local.is_empty() {
            bins.push(BinCell::invalid());
            continue;
        }
        let (lm, lp, rm, rp) =
            lig_rec_at_bin(&lig_local, &rec_local, grid, source, target, lineage, bin);
        bins.push(bin_cell_from_support(
            ws,
            wt,
            lp,
            rp,
            lm,
            rm,
            min_pct,
            min_n_eff,
            grid.hill_s,
        ));
    }
    TrajectoryCurve { bins }
}

fn lig_rec_at_bin_gene_stats(
    ligand_rows: &[usize],
    receptor_rows: &[usize],
    stats: &GeneBinStats,
    source: usize,
    target: usize,
    lineage: usize,
    bin: usize,
) -> (f32, f32, f32, f32) {
    let lm = ligand_rows
        .iter()
        .map(|&g| stats.get_mean_global(g, source, lineage, bin))
        .fold(f32::INFINITY, f32::min);
    let lp = ligand_rows
        .iter()
        .map(|&g| stats.get_pct_global(g, source, lineage, bin))
        .fold(f32::INFINITY, f32::min);
    let rm = receptor_rows
        .iter()
        .map(|&g| stats.get_mean_global(g, target, lineage, bin))
        .fold(f32::INFINITY, f32::min);
    let rp = receptor_rows
        .iter()
        .map(|&g| stats.get_pct_global(g, target, lineage, bin))
        .fold(f32::INFINITY, f32::min);
    (lm, lp, rm, rp)
}

pub fn build_lr_trajectory_curve_from_gene_stats(
    pair: &LrPair,
    stats: &GeneBinStats,
    source: usize,
    target: usize,
    lineage: usize,
    hill_s: f32,
    n_bins: usize,
    min_pct: f32,
    min_n_eff: f32,
) -> TrajectoryCurve {
    let slb = stats.n_lineages * n_bins;
    let sb = n_bins;
    let mut bins = Vec::with_capacity(n_bins);
    for bin in 0..n_bins {
        let source_idx = source * slb + lineage * sb + bin;
        let target_idx = target * slb + lineage * sb + bin;
        let ws = kish_ess(stats.wsums[source_idx], stats.w2sums[source_idx]);
        let wt = kish_ess(stats.wsums[target_idx], stats.w2sums[target_idx]);
        let (lm, lp, rm, rp) = lig_rec_at_bin_gene_stats(
            &pair.ligand_rows,
            &pair.receptor_rows,
            stats,
            source,
            target,
            lineage,
            bin,
        );
        bins.push(bin_cell_from_support(
            ws, wt, lp, rp, lm, rm, min_pct, min_n_eff, hill_s,
        ));
    }
    TrajectoryCurve { bins }
}

pub fn trajectory_keys_from_scores(scores: &[crate::ccc::lr_scoring::LrBinResult]) -> Vec<TrajKey> {
    let mut keys: Vec<TrajKey> = scores
        .iter()
        .map(|r| (r.pair_idx, r.source, r.target, r.lineage))
        .collect();
    keys.sort_unstable();
    keys.dedup();
    keys
}

fn min_n_eff_cell(cell: &BinCell) -> f32 {
    cell.n_eff_source.min(cell.n_eff_target)
}

pub fn evaluate_trajectory_curve(
    curve: &TrajectoryCurve,
    stat_code: usize,
    top_k: usize,
    cfg: &TrajectoryEvalConfig,
) -> TrajectoryEvalSummary {
    let n_bins = curve.bins.len();
    let n_valid = curve.bins.iter().filter(|b| b.structural_valid).count();
    let n_expression_supported_bins = curve
        .bins
        .iter()
        .filter(|b| b.structural_valid && b.expression_supported)
        .count();
    let valid_bin_fraction = if n_bins > 0 {
        n_valid as f32 / n_bins as f32
    } else {
        0.0
    };

    let mut evaluable =
        n_valid >= cfg.min_valid_bins && valid_bin_fraction >= cfg.min_valid_fraction;
    let mut statistic = f32::NAN;
    let mut n_eff_at_max = f32::NAN;
    let mut n_eff_at_min = f32::NAN;
    let mut valid_at_max = false;
    let mut valid_at_min = false;

    if evaluable {
        let valid_indices: Vec<usize> = curve
            .bins
            .iter()
            .enumerate()
            .filter(|(_, b)| b.structural_valid)
            .map(|(i, _)| i)
            .collect();

        match stat_code {
            1 => {
                statistic = valid_indices
                    .iter()
                    .map(|&i| curve.bins[i].score)
                    .fold(f32::NEG_INFINITY, f32::max);
                if let Some(&imax) = valid_indices.iter().max_by(|&&a, &&b| {
                    curve.bins[a]
                        .score
                        .partial_cmp(&curve.bins[b].score)
                        .unwrap_or(std::cmp::Ordering::Equal)
                }) {
                    n_eff_at_max = min_n_eff_cell(&curve.bins[imax]);
                    valid_at_max = true;
                }
            }
            2 => {
                // mean_score: equal-weight average over structural-valid bins (not time integral).
                let sum: f32 = valid_indices.iter().map(|&i| curve.bins[i].score).sum();
                statistic = sum / valid_indices.len() as f32;
                valid_at_max = true;
            }
            3 => {
                let mut scores: Vec<f32> =
                    valid_indices.iter().map(|&i| curve.bins[i].score).collect();
                scores.sort_by(|a, b| b.partial_cmp(a).unwrap_or(std::cmp::Ordering::Equal));
                let kk = top_k.max(1).min(scores.len());
                statistic = scores[..kk].iter().sum::<f32>() / kk as f32;
                valid_at_max = true;
            }
            _ => {
                // amplitude: max - min among structural-valid bins; need ≥2 valid bins
                if valid_indices.len() < 2 {
                    evaluable = false;
                } else {
                    let (imax, vmax) = valid_indices
                        .iter()
                        .map(|&i| (i, curve.bins[i].score))
                        .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal))
                        .unwrap();
                    let (imin, vmin) = valid_indices
                        .iter()
                        .map(|&i| (i, curve.bins[i].score))
                        .min_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal))
                        .unwrap();
                    statistic = vmax - vmin;
                    n_eff_at_max = min_n_eff_cell(&curve.bins[imax]);
                    n_eff_at_min = min_n_eff_cell(&curve.bins[imin]);
                    valid_at_max = curve.bins[imax].structural_valid;
                    valid_at_min = curve.bins[imin].structural_valid;
                }
            }
        }
        if !statistic.is_finite() {
            evaluable = false;
        }
    }

    TrajectoryEvalSummary {
        statistic,
        evaluable,
        n_valid_bins: n_valid,
        valid_bin_fraction,
        n_expression_supported_bins,
        n_eff_at_max,
        n_eff_at_min,
        valid_at_max,
        valid_at_min,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn kish_ess_distinguishes_weight_mass_from_independent_support() {
        assert!((kish_ess(1.0, 1.0) - 1.0).abs() < 1e-6);
        assert!((kish_ess(1.0, 0.25) - 4.0).abs() < 1e-6);
        assert!((kish_ess(2.0, 1.0) - 4.0).abs() < 1e-6);
    }

    #[test]
    fn missing_structural_bin_not_treated_as_zero_amplitude() {
        let curve = TrajectoryCurve {
            bins: vec![
                BinCell::invalid(),
                BinCell {
                    score: 0.8,
                    structural_valid: true,
                    expression_supported: true,
                    n_eff_source: 20.0,
                    n_eff_target: 20.0,
                },
                BinCell::invalid(),
            ],
        };
        let cfg = TrajectoryEvalConfig {
            min_valid_bins: 2,
            min_valid_fraction: 0.25,
            ..Default::default()
        };
        let ev = evaluate_trajectory_curve(&curve, 0, 3, &cfg);
        assert!(
            !ev.evaluable,
            "single structural bin cannot define amplitude"
        );
    }

    #[test]
    fn expression_zero_is_structural_valid_score_zero() {
        let cell = bin_cell_from_support(10.0, 10.0, 0.01, 0.5, 0.0, 0.5, 0.05, 5.0, 1.0);
        assert!(cell.structural_valid);
        assert!(!cell.expression_supported);
        assert_eq!(cell.score, 0.0);
    }

    #[test]
    fn amplitude_uses_two_structural_bins() {
        let curve = TrajectoryCurve {
            bins: vec![
                BinCell {
                    score: 0.2,
                    structural_valid: true,
                    expression_supported: true,
                    n_eff_source: 20.0,
                    n_eff_target: 20.0,
                },
                BinCell::invalid(),
                BinCell {
                    score: 0.9,
                    structural_valid: true,
                    expression_supported: true,
                    n_eff_source: 20.0,
                    n_eff_target: 20.0,
                },
            ],
        };
        let cfg = TrajectoryEvalConfig {
            min_valid_bins: 2,
            min_valid_fraction: 0.25,
            ..Default::default()
        };
        let ev = evaluate_trajectory_curve(&curve, 0, 3, &cfg);
        assert!(ev.evaluable);
        assert!((ev.statistic - 0.7).abs() < 1e-6);
        assert!(ev.valid_at_max && ev.valid_at_min);
    }

    #[test]
    fn mean_score_does_not_mark_valid_at_min() {
        let curve = TrajectoryCurve {
            bins: vec![
                BinCell {
                    score: 0.2,
                    structural_valid: true,
                    expression_supported: true,
                    n_eff_source: 20.0,
                    n_eff_target: 20.0,
                },
                BinCell {
                    score: 0.8,
                    structural_valid: true,
                    expression_supported: true,
                    n_eff_source: 20.0,
                    n_eff_target: 20.0,
                },
            ],
        };
        let cfg = TrajectoryEvalConfig {
            min_valid_bins: 2,
            min_valid_fraction: 0.25,
            ..Default::default()
        };
        let ev = evaluate_trajectory_curve(&curve, 2, 3, &cfg);
        assert!(ev.valid_at_max);
        assert!(!ev.valid_at_min);
    }

    #[test]
    fn edge_rewire_resolves_real_gene_from_real_stats() {
        use crate::ccc::fake_lr::GeneBinStats;
        let real_stats = GeneBinStats {
            means: vec![3.0, 0.0],
            pcts: vec![0.8, 0.0],
            wsums: vec![10.0, 10.0],
            w2sums: vec![10.0, 10.0],
            gene_rows: vec![5],
            n_genes: 1,
            n_types: 1,
            n_lineages: 1,
            n_bins: 2,
        };
        let fake_stats = GeneBinStats {
            means: vec![1.0, 2.0],
            pcts: vec![0.9, 0.9],
            wsums: vec![10.0, 10.0],
            w2sums: vec![10.0, 10.0],
            gene_rows: vec![10],
            n_genes: 1,
            n_types: 1,
            n_lineages: 1,
            n_bins: 2,
        };
        let structural = crate::ccc::lr_scoring::LrBinGrid {
            means: vec![],
            pcts: vec![],
            wsums: vec![10.0, 10.0],
            w2sums: vec![10.0, 10.0],
            g2local: vec![],
            bin_bounds: vec![],
            hill_s: 1.0,
            n_types: 1,
            n_lineages: 1,
            n_bins: 2,
        };
        let curve = build_fake_trajectory_curve(
            5,
            10,
            0,
            0,
            0,
            &real_stats,
            &fake_stats,
            &structural,
            0.05,
            5.0,
        );
        assert!(
            curve.bins[0].score > 0.0,
            "real ligand expression must resolve from real_stats, not zero"
        );
    }

    #[test]
    fn grid_and_gene_stats_real_curves_agree_for_simple_pair() {
        use crate::ccc::lr_scoring::LrBinGrid;
        use crate::ccc::LrPair;
        let pair = LrPair {
            ligand_rows: vec![0],
            receptor_rows: vec![1],
        };
        let stats = GeneBinStats {
            means: vec![2.0, 0.0, 0.0, 3.0],
            pcts: vec![0.8, 0.0, 0.0, 0.9],
            wsums: vec![10.0, 10.0],
            w2sums: vec![10.0, 10.0],
            gene_rows: vec![0, 1],
            n_genes: 2,
            n_types: 1,
            n_lineages: 1,
            n_bins: 2,
        };
        let grid = LrBinGrid {
            means: stats.means.clone(),
            pcts: stats.pcts.clone(),
            wsums: stats.wsums.clone(),
            w2sums: stats.w2sums.clone(),
            g2local: {
                let mut g2 = vec![usize::MAX; 2];
                g2[0] = 0;
                g2[1] = 1;
                g2
            },
            bin_bounds: vec![],
            hill_s: 1.5,
            n_types: 1,
            n_lineages: 1,
            n_bins: 2,
        };
        let from_stats =
            build_lr_trajectory_curve_from_gene_stats(&pair, &stats, 0, 0, 0, 1.5, 2, 0.05, 5.0);
        let from_grid = build_lr_trajectory_curve_from_grid(&pair, &grid, 0, 0, 0, 0.05, 5.0);
        for (a, b) in from_stats.bins.iter().zip(from_grid.bins.iter()) {
            assert_eq!(a.structural_valid, b.structural_valid);
            assert!((a.score - b.score).abs() < 1e-6);
        }
    }
}
