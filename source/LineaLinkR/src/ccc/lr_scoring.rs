//! Pseudotime-binned and cluster-level ligand-receptor scoring.
//!
//! Cluster-level scores use Michaelis-Menten saturation with cell-type label
//! permutation p-values. Temporal scores bin expression along lineage-weighted
//! equal-frequency pseudotime and report score, effective cell count, and lineage
//! entropy per `(pair, source, target, lineage, bin)`.

use crate::core::SparseMatrix;
use rand::{rngs::StdRng, seq::SliceRandom, SeedableRng};
use rayon::prelude::*;
use std::collections::HashSet;

#[derive(Debug, Clone)]
pub struct LrPair {
    pub ligand_rows: Vec<usize>,
    pub receptor_rows: Vec<usize>,
}

/// Per-gene binned means/pcts and cluster weights from one CSC walk.
/// Used to build trajectory curves with structural-valid bin semantics.
#[derive(Debug, Clone)]
pub struct LrBinGrid {
    pub means: Vec<f32>,
    pub pcts: Vec<f32>,
    pub wsums: Vec<f64>,
    pub w2sums: Vec<f64>,
    pub g2local: Vec<usize>,
    pub bin_bounds: Vec<Vec<f32>>,
    pub hill_s: f32,
    pub n_types: usize,
    pub n_lineages: usize,
    pub n_bins: usize,
}

impl LrBinGrid {
    #[inline]
    pub fn sl(&self) -> usize {
        self.n_types * self.n_lineages * self.n_bins
    }

    #[inline]
    pub fn slb(&self) -> usize {
        self.n_lineages * self.n_bins
    }

    #[inline]
    pub fn sb(&self) -> usize {
        self.n_bins
    }
}

#[derive(Debug, Clone)]
pub struct LrResult {
    pub source: usize,
    pub target: usize,
    pub pair_idx: usize,
    pub score: f32,
    pub pvalue: f32,
    pub pct_ligand: f32,
    pub pct_receptor: f32,
}

#[derive(Debug, Clone)]
pub struct LrBinResult {
    pub pair_idx: usize,
    pub source: usize,
    pub target: usize,
    pub lineage: usize,
    pub bin: usize,
    pub score: f32,
    /// Minimum Kish effective sample size between source and target in this
    /// (lineage, bin): `(sum w)^2 / sum(w^2)`.
    pub n_eff_cells: f32,
    /// Minimum raw lineage-weight mass between source and target.  Kept
    /// separate from Kish ESS so support thresholds have an unambiguous scale.
    pub weight_mass_cells: f32,
    /// Max(mean entropy(source), mean entropy(target)) in this bin.
    /// Entropy is normalised H(Ψ)/log(L) in [0,1]; high values indicate
    /// pre-bifurcation lineage ambiguity.
    pub lineage_entropy: f32,
}

/// One CSC pass → (means, pct_expressed) per gene per cluster.  O(nnz)
///
/// `means[g][ct]` = arithmetic mean over ALL cells (including dropout zeros).
/// `pcts[g][ct]`  = fraction of cells in cluster ct with observed value > 0.
///
/// The score uses the full-population mean, matching CellPhoneDB v2 and avoiding
/// inflation from rare expressors. `pct_expressed` provides the population-level gate.
pub fn cluster_stats_csc(
    expr: &SparseMatrix,
    cell_types: &[usize],
    n_types: usize,
) -> (Vec<Vec<f32>>, Vec<Vec<f32>>) {
    let n_cells = expr.ncols;
    let n_genes = expr.nrows;
    let mut counts = vec![0u32; n_types];
    for &ct in cell_types.iter().take(n_cells) {
        if ct < n_types {
            counts[ct] += 1;
        }
    }

    // Parallel fold over cells (columns). Write conflict on [g*n_types+ct] is
    // resolved by per-thread accumulator Vecs; reduce merges them in O(n_genes×n_types).
    let flat = n_genes * n_types;
    let (sums, n_expr) = (0..n_cells)
        .into_par_iter()
        .fold(
            || (vec![0.0f64; flat], vec![0u32; flat]),
            |(mut s, mut ne), cell| {
                let ct = cell_types[cell];
                if ct >= n_types {
                    return (s, ne);
                } // skip NA / out-of-range
                for nnz in expr.col_ptrs[cell]..expr.col_ptrs[cell + 1] {
                    let g = expr.rows[nnz];
                    let value = expr.values[nnz];
                    s[g * n_types + ct] += value as f64;
                    if value > 0.0 {
                        ne[g * n_types + ct] += 1;
                    }
                }
                (s, ne)
            },
        )
        .reduce(
            || (vec![0.0f64; flat], vec![0u32; flat]),
            |(mut a, mut an), (b, bn)| {
                for i in 0..flat {
                    a[i] += b[i];
                    an[i] += bn[i];
                }
                (a, an)
            },
        );

    let means: Vec<Vec<f32>> = (0..n_genes)
        .map(|g| {
            (0..n_types)
                .map(|ct| {
                    let c = counts[ct];
                    if c > 0 {
                        (sums[g * n_types + ct] / c as f64) as f32
                    } else {
                        0.0
                    }
                })
                .collect()
        })
        .collect();
    let pcts: Vec<Vec<f32>> = (0..n_genes)
        .map(|g| {
            (0..n_types)
                .map(|ct| {
                    let c = counts[ct];
                    if c > 0 {
                        n_expr[g * n_types + ct] as f32 / c as f32
                    } else {
                        0.0
                    }
                })
                .collect()
        })
        .collect();
    (means, pcts)
}

/// Cluster statistics restricted to the genes referenced by the LR database.
///
/// The permutation null only consumes LR genes. Keeping the global-to-local
/// map lets each CSC column be scanned once while avoiding full-transcriptome
/// thread-local accumulators and reductions on every permutation.
fn cluster_stats_csc_selected(
    expr: &SparseMatrix,
    cell_types: &[usize],
    n_types: usize,
    gene_list: &[usize],
    g2local: &[usize],
) -> (Vec<Vec<f32>>, Vec<Vec<f32>>) {
    let n_cells = expr.ncols;
    let n_local = gene_list.len();
    let mut counts = vec![0u32; n_types];
    for &ct in cell_types.iter().take(n_cells) {
        if ct < n_types {
            counts[ct] += 1;
        }
    }

    let flat = n_local * n_types;
    let (sums, n_expr) = (0..n_cells)
        .into_par_iter()
        .fold(
            || (vec![0.0f64; flat], vec![0u32; flat]),
            |(mut s, mut ne), cell| {
                let ct = cell_types[cell];
                if ct >= n_types {
                    return (s, ne);
                }
                for nnz in expr.col_ptrs[cell]..expr.col_ptrs[cell + 1] {
                    let local = g2local[expr.rows[nnz]];
                    if local == usize::MAX {
                        continue;
                    }
                    let value = expr.values[nnz];
                    let idx = local * n_types + ct;
                    s[idx] += value as f64;
                    if value > 0.0 {
                        ne[idx] += 1;
                    }
                }
                (s, ne)
            },
        )
        .reduce(
            || (vec![0.0f64; flat], vec![0u32; flat]),
            |(mut a, mut an), (b, bn)| {
                for i in 0..flat {
                    a[i] += b[i];
                    an[i] += bn[i];
                }
                (a, an)
            },
        );

    let means = (0..n_local)
        .map(|g| {
            (0..n_types)
                .map(|ct| {
                    let c = counts[ct];
                    if c > 0 {
                        (sums[g * n_types + ct] / c as f64) as f32
                    } else {
                        0.0
                    }
                })
                .collect()
        })
        .collect();
    let pcts = (0..n_local)
        .map(|g| {
            (0..n_types)
                .map(|ct| {
                    let c = counts[ct];
                    if c > 0 {
                        n_expr[g * n_types + ct] as f32 / c as f32
                    } else {
                        0.0
                    }
                })
                .collect()
        })
        .collect();
    (means, pcts)
}

/// Fast means-only variant for the permutation inner loop.  O(nnz)
/// Skips pct tracking to minimise allocation in the hot path (1000+ calls).
fn cluster_means_only(expr: &SparseMatrix, cell_types: &[usize], n_types: usize) -> Vec<Vec<f32>> {
    let n_cells = expr.ncols;
    let n_genes = expr.nrows;
    let mut counts = vec![0u32; n_types];
    for &ct in cell_types.iter().take(n_cells) {
        if ct < n_types {
            counts[ct] += 1;
        }
    }
    let mut sums = vec![0.0f64; n_genes * n_types];
    for cell in 0..n_cells {
        let ct = cell_types[cell];
        if ct >= n_types {
            continue;
        } // skip NA / out-of-range
        for nnz in expr.col_ptrs[cell]..expr.col_ptrs[cell + 1] {
            sums[expr.rows[nnz] * n_types + ct] += expr.values[nnz] as f64;
        }
    }
    (0..n_genes)
        .map(|g| {
            (0..n_types)
                .map(|ct| {
                    let c = counts[ct];
                    if c > 0 {
                        (sums[g * n_types + ct] / c as f64) as f32
                    } else {
                        0.0
                    }
                })
                .collect()
        })
        .collect()
}

/// Returns (complex_mean, complex_min_pct) per cluster.
/// Biological: a hetero-oligomeric complex is rate-limited by its least-expressed
/// subunit — both mean activity AND percent expressed are bounded by the minimum.
fn complex_stats(
    rows: &[usize],
    means: &[Vec<f32>],
    pcts: &[Vec<f32>],
    n: usize,
) -> (Vec<f32>, Vec<f32>) {
    let cm: Vec<f32> = (0..n)
        .map(|ct| {
            rows.iter()
                .map(|&r| means[r][ct])
                .fold(f32::INFINITY, f32::min)
        })
        .collect();
    let cp: Vec<f32> = (0..n)
        .map(|ct| {
            rows.iter()
                .map(|&r| pcts[r][ct])
                .fold(f32::INFINITY, f32::min)
        })
        .collect();
    (cm, cp)
}

/// Global Hill saturation scale S for the Michaelis-Menten-style LR score:
///   score(L,R) = sqrt(L·R) / (S + sqrt(L·R))
///
/// S is the **arithmetic mean of cluster-level expression across all non-zero
/// (gene, cluster) entries** — i.e. mean of `means[g][ct]` over `means[g][ct] > 0`.
/// This is the SingleCellSignalR scale (Cabello-Aguilar 2020,
/// Bioinformatics), where μ is the whole-dataset mean expression. It is not the
/// mean of `sqrt(L·R)` over LR pairs.
///
/// At this scale the LR score evaluates to ~0.5 for an "average" interaction,
/// gracefully saturating at 1.0 for highly abundant pairs and approaching 0 for
/// low-abundance pairs. The saturation prevents high-abundance housekeeping
/// genes (GAPDH, ACTB, etc.) from inflating LR scores when they appear in the
/// LR database.
///
/// **Do not change the implementation without updating the lock test in
/// `tests::hill_s_lock_simple`.** Any deviation from the documented arithmetic
/// mean alters the saturation scale across all downstream LR scores and
/// permutation null distributions.
fn compute_hill_s(means: &[Vec<f32>], n_types: usize) -> f32 {
    let mut sum = 0.0f64;
    let mut cnt = 0u64;
    for row in means {
        for &v in row.iter().take(n_types) {
            if v > 0.0 {
                sum += v as f64;
                cnt += 1;
            }
        }
    }
    if cnt == 0 {
        1.0
    } else {
        (sum / cnt as f64) as f32
    }
}

/// Dataset-level Hill scale shared by every LR subset and temporal null.
///
/// The scale is computed from all genes and cell-type means, so adding or
/// removing LR pairs cannot silently rescale scores for pairs already present.
pub(crate) fn compute_dataset_hill_s(
    expr: &SparseMatrix,
    cell_types: &[usize],
    n_types: usize,
    pseudotime: &[f32],
) -> f32 {
    let temporal_cell_types: Vec<usize> = cell_types
        .iter()
        .zip(pseudotime)
        .map(|(&cell_type, &pt)| if pt.is_nan() { n_types } else { cell_type })
        .collect();
    let means = cluster_means_only(expr, &temporal_cell_types, n_types);
    compute_hill_s(&means, n_types)
}

/// Equal-frequency bin boundaries for one (pseudotime, weight) series.
///
/// `wpt` is (pseudotime, weight) per contributing cell.  Returns
/// `bounds[k]` = pseudotime cutoff where bin k ends; length `n_bins+1`
/// with `bounds[0]=0`, `bounds[n_bins]=1`.
fn equal_frequency_bound_row(wpt: Vec<(f32, f32)>, n_bins: usize) -> Vec<f32> {
    let mut wpt = wpt;
    wpt.sort_unstable_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    let mut grouped: Vec<(f32, f32)> = Vec::new();
    for (pt, weight) in wpt {
        if let Some((last_pt, last_weight)) = grouped.last_mut() {
            if *last_pt == pt {
                *last_weight += weight;
                continue;
            }
        }
        grouped.push((pt, weight));
    }
    let total_w: f32 = grouped.iter().map(|(_, w)| *w).sum();
    if total_w < 1e-8 {
        return (0..=n_bins).map(|k| k as f32 / n_bins as f32).collect();
    }
    if grouped.len() < n_bins {
        return (0..=n_bins).map(|k| k as f32 / n_bins as f32).collect();
    }
    let target = total_w / n_bins as f32;
    let mut bounds = vec![0.0f32; n_bins + 1];
    bounds[n_bins] = 1.0;
    let cumulative: Vec<f32> = grouped
        .iter()
        .scan(0.0f32, |sum, (_, weight)| {
            *sum += *weight;
            Some(*sum)
        })
        .collect();
    let mut previous_cut: Option<usize> = None;
    for k in 1..n_bins {
        let ideal = cumulative
            .partition_point(|&value| value < target * k as f32)
            .min(grouped.len() - 2);
        let lower = previous_cut.map_or(0, |cut| cut + 1);
        let upper = grouped.len() - (n_bins - k) - 1;
        let cut = ideal.clamp(lower, upper);
        bounds[k] = (grouped[cut].0 + grouped[cut + 1].0) * 0.5;
        previous_cut = Some(cut);
    }
    bounds
}

/// Equal-frequency bin boundaries per lineage.
///
/// For each lineage, collects (pseudotime, weight) pairs, sorts by pseudotime,
/// and computes boundary cut points so that each bin receives approximately
/// equal total lineage weight regardless of pseudotime distribution shape.
///
/// Returns `bin_bounds[lin][k]` = pseudotime cutoff where bin k ends.
/// Length `n_bins+1`; `bounds[0]=0`, `bounds[n_bins]=1`.
pub fn compute_bin_bounds(
    n_cells: usize,
    pseudotime: &[f32],
    lineage_weights: &[f32],
    n_lineages: usize,
    n_bins: usize,
) -> Vec<Vec<f32>> {
    (0..n_lineages)
        .map(|lin| {
            let wpt: Vec<(f32, f32)> = (0..n_cells)
                .filter_map(|c| {
                    let pt = pseudotime[c];
                    if pt.is_nan() {
                        return None;
                    }
                    let w = lineage_weights[c * n_lineages + lin];
                    if w < 1e-8 {
                        return None;
                    }
                    Some((pt, w))
                })
                .collect();
            equal_frequency_bound_row(wpt, n_bins)
        })
        .collect()
}

/// Common temporal grid shared by every lineage.
///
/// Per-lineage equal-frequency bounds make a bin index lineage-relative: the
/// same bin number can cover different pseudotime ranges on sibling branches,
/// so differences of per-lineage pre/post changes need not cancel a shared
/// smooth temporal trend.  Fork-local contrasts compare sister branches at
/// one shared topology event on one shared graph-time coordinate, which
/// requires every lineage to be scored on a single boundary vector.
///
/// Bounds are equal-frequency over the pooled on-trajectory mass; a cell's
/// pooled mass is its maximum lineage weight.  The returned matrix repeats
/// the same boundary row for every lineage and can be passed directly to the
/// binned scorer's `bin_bounds_override`.
pub fn compute_common_bin_bounds(
    n_cells: usize,
    pseudotime: &[f32],
    lineage_weights: &[f32],
    n_lineages: usize,
    n_bins: usize,
) -> Vec<Vec<f32>> {
    let wpt: Vec<(f32, f32)> = (0..n_cells)
        .filter_map(|c| {
            let pt = pseudotime[c];
            if pt.is_nan() {
                return None;
            }
            let w = lineage_weights[c * n_lineages..c * n_lineages + n_lineages]
                .iter()
                .copied()
                .filter(|x| x.is_finite() && *x > 0.0)
                .fold(0.0f32, f32::max);
            if w < 1e-8 {
                return None;
            }
            Some((pt, w))
        })
        .collect();
    let shared = equal_frequency_bound_row(wpt, n_bins);
    (0..n_lineages).map(|_| shared.clone()).collect()
}

/// Collect unique gene indices referenced by LR pairs and build O(1) lookup.
fn collect_lr_genes(lr_pairs: &[LrPair], nrows: usize) -> (Vec<usize>, Vec<usize>) {
    let mut gene_set: HashSet<usize> = HashSet::new();
    for pair in lr_pairs {
        for &g in &pair.ligand_rows {
            gene_set.insert(g);
        }
        for &g in &pair.receptor_rows {
            gene_set.insert(g);
        }
    }
    let mut gene_list: Vec<usize> = gene_set.into_iter().collect();
    gene_list.sort_unstable();
    let mut g2local = vec![usize::MAX; nrows];
    for (local, &global) in gene_list.iter().enumerate() {
        g2local[global] = local;
    }
    (gene_list, g2local)
}

/// LR score matrix with population-level pct_expressed gate.
///
/// Score function: sqrt(L·R) / (S + sqrt(L·R)) — Michaelis-Menten saturation.
/// This is the exact formula from SingleCellSignalR (Cabello-Aguilar 2020):
///   LRscore = sqrt(L)·sqrt(R) / (S + sqrt(L)·sqrt(R))
/// where S is the global mean expression scale (see `compute_hill_s`).
/// The saturation term prevents high-abundance housekeeping genes (e.g., GAPDH,
/// ACTB) from generating inflated scores when they happen to appear in the
/// LR database.  Without S, the score grows linearly with abundance — with S,
/// it saturates toward 1.0 and gives a fair comparison across expression levels.
fn lr_score_matrix(
    lm: &[f32],
    rm: &[f32],
    lp: &[f32],
    rp: &[f32],
    n: usize,
    min_pct: f32,
    hill_s: f32,
) -> Vec<f32> {
    let mut out = vec![0.0f32; n * n];
    for s in 0..n {
        if lp[s] < min_pct {
            continue;
        }
        for t in 0..n {
            if rp[t] < min_pct {
                continue;
            }
            let l = lm[s];
            let r = rm[t];
            if l > 0.0 && r > 0.0 {
                let x = (l * r).sqrt();
                out[s * n + t] = x / (hill_s + x);
            }
        }
    }
    out
}

/// LR communication scores with permutation p-values.
/// `min_pct`: minimum fraction of cells expressing a subunit (default 0.05 in R).
pub fn compute_lr_scores(
    expr: &SparseMatrix,
    cell_types: &[usize],
    n_types: usize,
    lr_pairs: &[LrPair],
    n_perms: usize,
    seed: u64,
    min_pct: f32,
) -> Vec<LrResult> {
    let (means, pcts) = cluster_stats_csc(expr, cell_types, n_types);
    let n_pairs = lr_pairs.len();
    let nt2 = n_types * n_types;

    let hill_s = compute_hill_s(&means, n_types);
    let (gene_list, g2local) = collect_lr_genes(lr_pairs, expr.nrows);
    let local_pairs: Vec<LrPair> = lr_pairs
        .iter()
        .map(|pair| LrPair {
            ligand_rows: pair.ligand_rows.iter().map(|&g| g2local[g]).collect(),
            receptor_rows: pair.receptor_rows.iter().map(|&g| g2local[g]).collect(),
        })
        .collect();
    let selected_means: Vec<Vec<f32>> = gene_list.iter().map(|&g| means[g].clone()).collect();
    let selected_pcts: Vec<Vec<f32>> = gene_list.iter().map(|&g| pcts[g].clone()).collect();

    // Real scores with pct gate; also store per-pair pct for output
    let real_data: Vec<(Vec<f32>, Vec<f32>, Vec<f32>)> = local_pairs
        .par_iter()
        .map(|pair| {
            let (lm, lp) =
                complex_stats(&pair.ligand_rows, &selected_means, &selected_pcts, n_types);
            let (rm, rp) = complex_stats(
                &pair.receptor_rows,
                &selected_means,
                &selected_pcts,
                n_types,
            );
            let score = lr_score_matrix(&lm, &rm, &lp, &rp, n_types, min_pct, hill_s);
            (score, lp, rp)
        })
        .collect();

    let real: Vec<f32> = real_data
        .iter()
        .flat_map(|(s, _, _)| s.iter().copied())
        .collect();

    // Permutation counts: recompute both means and detection fractions so the
    // null applies the identical min_pct gate as the observed statistic.
    let perm_counts: Vec<u32> = (0..n_perms)
        .into_par_iter()
        .map(|pi| {
            let mut rng = StdRng::seed_from_u64(seed.wrapping_add(pi as u64));
            let mut shuf = cell_types.to_vec();
            shuf.shuffle(&mut rng);
            let (pm, pp) = cluster_stats_csc_selected(expr, &shuf, n_types, &gene_list, &g2local);
            let mut local = vec![0u32; n_pairs * nt2];
            for (p, pair) in local_pairs.iter().enumerate() {
                let (lm, lp) = complex_stats(&pair.ligand_rows, &pm, &pp, n_types);
                let (rm, rp) = complex_stats(&pair.receptor_rows, &pm, &pp, n_types);
                let ps = lr_score_matrix(&lm, &rm, &lp, &rp, n_types, min_pct, hill_s);
                for (k, (&ps_k, &rs_k)) in ps.iter().zip(real[p * nt2..].iter()).enumerate() {
                    if ps_k >= rs_k {
                        local[p * nt2 + k] += 1;
                    }
                }
            }
            local
        })
        .reduce(
            || vec![0u32; n_pairs * nt2],
            |mut a, b| {
                for i in 0..a.len() {
                    a[i] += b[i];
                }
                a
            },
        );

    let mut results = Vec::new();
    for p in 0..n_pairs {
        let (_, ref lp, ref rp) = real_data[p];
        for s in 0..n_types {
            for t in 0..n_types {
                let k = p * nt2 + s * n_types + t;
                let score = real[k];
                if score > 0.0 {
                    results.push(LrResult {
                        source: s,
                        target: t,
                        pair_idx: p,
                        score,
                        pvalue: if n_perms > 0 {
                            (perm_counts[k] as f32 + 1.0) / (n_perms as f32 + 1.0)
                        } else {
                            1.0
                        },
                        pct_ligand: lp[s],
                        pct_receptor: rp[t],
                    });
                }
            }
        }
    }
    results
}

/// Pseudotime-binned LR scoring (LineaLinkR Module 1) — memory-efficient.
///
/// The accumulator is sized from the unique gene indices used by the LR
/// database rather than from the full expression matrix, so memory scales with
/// the scored gene set.
///
/// lineage_weights: flat [n_cells × n_lineages], row-major.
/// min_pct: minimum weighted fraction of cells expressing a subunit per bin (default 0.05).
pub fn compute_lr_scores_bins(
    expr: &SparseMatrix,
    cell_types: &[usize],
    n_types: usize,
    lr_pairs: &[LrPair],
    pseudotime: &[f32],
    lineage_weights: &[f32],
    n_lineages: usize,
    n_bins: usize,
    min_pct: f32,
) -> (Vec<LrBinResult>, crate::ccc::taxonomy::LrBinIntermediate) {
    let (scores, intermediate, _, _) = compute_lr_scores_bins_internal(
        expr,
        cell_types,
        n_types,
        lr_pairs,
        pseudotime,
        lineage_weights,
        n_lineages,
        n_bins,
        min_pct,
        None,
        None,
    );
    (scores, intermediate)
}

#[cfg(test)]
pub(crate) fn compute_lr_scores_bins_with_scale(
    expr: &SparseMatrix,
    cell_types: &[usize],
    n_types: usize,
    lr_pairs: &[LrPair],
    pseudotime: &[f32],
    lineage_weights: &[f32],
    n_lineages: usize,
    n_bins: usize,
    min_pct: f32,
) -> (
    Vec<LrBinResult>,
    crate::ccc::taxonomy::LrBinIntermediate,
    f32,
) {
    let (scores, intermediate, hill_s, _) = compute_lr_scores_bins_internal(
        expr,
        cell_types,
        n_types,
        lr_pairs,
        pseudotime,
        lineage_weights,
        n_lineages,
        n_bins,
        min_pct,
        None,
        None,
    );
    (scores, intermediate, hill_s)
}

pub(crate) fn compute_lr_scores_bins_with_scale_grid(
    expr: &SparseMatrix,
    cell_types: &[usize],
    n_types: usize,
    lr_pairs: &[LrPair],
    pseudotime: &[f32],
    lineage_weights: &[f32],
    n_lineages: usize,
    n_bins: usize,
    min_pct: f32,
) -> (
    Vec<LrBinResult>,
    crate::ccc::taxonomy::LrBinIntermediate,
    f32,
    LrBinGrid,
) {
    compute_lr_scores_bins_internal(
        expr,
        cell_types,
        n_types,
        lr_pairs,
        pseudotime,
        lineage_weights,
        n_lineages,
        n_bins,
        min_pct,
        None,
        None,
    )
}

/// Reuses precomputed equal-frequency bin boundaries from the real pseudotime,
/// as required for valid temporal permutation nulls.
pub fn compute_lr_scores_bins_fixed_scale_grid_with_bounds(
    expr: &SparseMatrix,
    cell_types: &[usize],
    n_types: usize,
    lr_pairs: &[LrPair],
    pseudotime: &[f32],
    lineage_weights: &[f32],
    n_lineages: usize,
    n_bins: usize,
    min_pct: f32,
    hill_s: f32,
    bin_bounds: Option<&[Vec<f32>]>,
) -> (
    Vec<LrBinResult>,
    crate::ccc::taxonomy::LrBinIntermediate,
    LrBinGrid,
) {
    let (scores, intermediate, _, grid) = compute_lr_scores_bins_internal(
        expr,
        cell_types,
        n_types,
        lr_pairs,
        pseudotime,
        lineage_weights,
        n_lineages,
        n_bins,
        min_pct,
        Some(hill_s),
        bin_bounds,
    );
    (scores, intermediate, grid)
}

/// Binned scoring on a caller-selected temporal grid.
///
/// `common_grid = true` scores every lineage on one shared boundary vector
/// (see `compute_common_bin_bounds`); `false` reproduces the legacy
/// per-lineage equal-frequency grid.  The Hill scale is still estimated from
/// the data, as in `compute_lr_scores_bins_with_scale_grid`.
pub(crate) fn compute_lr_scores_bins_with_scale_grid_opts(
    expr: &SparseMatrix,
    cell_types: &[usize],
    n_types: usize,
    lr_pairs: &[LrPair],
    pseudotime: &[f32],
    lineage_weights: &[f32],
    n_lineages: usize,
    n_bins: usize,
    min_pct: f32,
    common_grid: bool,
) -> (
    Vec<LrBinResult>,
    crate::ccc::taxonomy::LrBinIntermediate,
    f32,
    LrBinGrid,
) {
    let bin_bounds = if common_grid {
        compute_common_bin_bounds(expr.ncols, pseudotime, lineage_weights, n_lineages, n_bins)
    } else {
        compute_bin_bounds(expr.ncols, pseudotime, lineage_weights, n_lineages, n_bins)
    };
    compute_lr_scores_bins_internal(
        expr,
        cell_types,
        n_types,
        lr_pairs,
        pseudotime,
        lineage_weights,
        n_lineages,
        n_bins,
        min_pct,
        None,
        Some(&bin_bounds),
    )
}

pub(crate) fn compute_lr_scores_bins_internal(
    expr: &SparseMatrix,
    cell_types: &[usize],
    n_types: usize,
    lr_pairs: &[LrPair],
    pseudotime: &[f32],
    lineage_weights: &[f32],
    n_lineages: usize,
    n_bins: usize,
    min_pct: f32,
    hill_s_override: Option<f32>,
    bin_bounds_override: Option<&[Vec<f32>]>,
) -> (
    Vec<LrBinResult>,
    crate::ccc::taxonomy::LrBinIntermediate,
    f32,
    LrBinGrid,
) {
    // Return signature includes LrBinIntermediate for downstream entropy.
    // from `Vec<LrBinResult>` to `(Vec<LrBinResult>, LrBinIntermediate)` so
    // that the `hwsums` / `wsums` accumulators live past the rayon region
    // and can feed `compute_lineage_entropy_raw` downstream.  `LrBinResult`
    // itself stays at 8 fields, one byte unchanged.
    use crate::ccc::taxonomy::LrBinIntermediate;
    let n_cells = expr.ncols;

    // ── Equal-frequency bin boundaries per lineage ───────────────────────────
    //
    // Weighted equal-frequency binning gives each bin approximately equal total
    // lineage weight regardless of the pseudotime distribution.
    //
    // bin_bounds[lin][k] = pseudotime cutoff above which bin > k for lineage lin.
    // Length n_bins+1; bounds[0]=0, bounds[n_bins]=1.
    if n_cells == 0 || n_lineages == 0 || n_bins == 0 {
        return (
            Vec::new(),
            LrBinIntermediate {
                hwsums: Vec::new(),
                wsums: Vec::new(),
                n_types,
                n_lineages,
                n_bins,
            },
            hill_s_override.unwrap_or(1.0),
            LrBinGrid {
                means: Vec::new(),
                pcts: Vec::new(),
                wsums: Vec::new(),
                w2sums: Vec::new(),
                g2local: Vec::new(),
                bin_bounds: Vec::new(),
                hill_s: hill_s_override.unwrap_or(1.0),
                n_types,
                n_lineages,
                n_bins,
            },
        );
    }

    let bin_bounds = if let Some(bounds) = bin_bounds_override {
        bounds.to_vec()
    } else {
        compute_bin_bounds(n_cells, pseudotime, lineage_weights, n_lineages, n_bins)
    };

    let (gene_list, g2local) = collect_lr_genes(lr_pairs, expr.nrows);
    let n_local = gene_list.len();

    // ── Allocate only for n_local genes (not all n_genes) ────────────────────
    let sl = n_types * n_lineages * n_bins; // stride per local gene
    let slb = n_lineages * n_bins; // stride per cluster
    let sb = n_bins; // stride per lineage
    let total = n_local * sl;
    let mut sums = vec![0.0f64; total]; // weighted expression sums
    let mut wsums = vec![0.0f64; n_types * n_lineages * n_bins]; // weight sums per (ct,lin,bin)
    let mut w2sums = vec![0.0f64; n_types * n_lineages * n_bins]; // squared-weight sums for Kish ESS
    let mut xwsums = vec![0.0f64; total]; // weight sums for expressing cells (pct numerator)
    let mut hwsums = vec![0.0f64; n_types * n_lineages * n_bins]; // weighted lineage entropy sums

    let ln_l = if n_lineages > 1 {
        (n_lineages as f32).ln().max(1e-6)
    } else {
        1.0
    };

    for cell in 0..n_cells {
        let ct = cell_types[cell];
        if ct >= n_types {
            continue;
        } // skip NA / out-of-range cell type
        let pt_raw = pseudotime[cell];
        if pt_raw.is_nan() {
            continue;
        } // disconnected cell: skip
        let pt = pt_raw.clamp(0.0, 1.0);

        // Per-cell lineage entropy H(Ψ)/log(L), computed once and reused across bins.
        let row0 = cell * n_lineages;
        let z: f32 = lineage_weights[row0..row0 + n_lineages]
            .iter()
            .copied()
            .filter(|w| w.is_finite() && *w > 0.0)
            .sum();
        let h_cell: f64 = if z > 1e-8 && n_lineages > 1 {
            let mut h = 0.0f32;
            for &w in &lineage_weights[row0..row0 + n_lineages] {
                if w.is_finite() && w > 0.0 {
                    let p = (w / z).max(1e-12);
                    h -= p * p.ln();
                }
            }
            (h / ln_l).clamp(0.0, 1.0) as f64
        } else {
            0.0f64
        };

        for lin in 0..n_lineages {
            let w = lineage_weights[cell * n_lineages + lin] as f64;
            if w < 1e-8 {
                continue;
            }
            // equal-frequency bin: binary search in pre-computed boundaries
            let bin = bin_bounds[lin]
                .partition_point(|&b| b <= pt)
                .saturating_sub(1)
                .min(n_bins - 1);
            let support_idx = ct * slb + lin * sb + bin;
            wsums[support_idx] += w;
            w2sums[support_idx] += w * w;
            hwsums[support_idx] += w * h_cell;
            for nnz in expr.col_ptrs[cell]..expr.col_ptrs[cell + 1] {
                let local = g2local[expr.rows[nnz]];
                if local == usize::MAX {
                    continue;
                }
                let base = local * sl + ct * slb + lin * sb + bin;
                let v = expr.values[nnz];
                sums[base] += v as f64 * w;
                if v > 0.0 {
                    xwsums[base] += w; // weighted pct numerator
                }
            }
        }
    }

    // Normalise: weighted_mean and weighted_pct per (local_gene, ct, lin, bin)
    let means: Vec<f32> = (0..total)
        .map(|i| {
            let ct = (i % sl) / slb;
            let lin = (i % slb) / sb;
            let bin = i % sb;
            let w = wsums[ct * slb + lin * sb + bin];
            if w > 1e-8 {
                (sums[i] / w) as f32
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
            let w = wsums[ct * slb + lin * sb + bin];
            if w > 1e-8 {
                (xwsums[i] / w) as f32
            } else {
                0.0
            }
        })
        .collect();

    // Compute S at dataset level. It must not depend on which LR pairs happen to be
    // supplied, otherwise adding an unrelated pair changes existing scores.
    // Permutation and FakeLR callers pass this same value back as an override.
    let hill_s_bins = hill_s_override
        .unwrap_or_else(|| compute_dataset_hill_s(expr, cell_types, n_types, pseudotime));

    let get = |local: usize, ct: usize, lin: usize, bin: usize| -> (f32, f32) {
        let i = local * sl + ct * slb + lin * sb + bin;
        (means[i], pcts[i])
    };

    let scores: Vec<LrBinResult> = lr_pairs
        .par_iter()
        .enumerate()
        .flat_map(|(p, pair)| {
            let lig_local: Vec<usize> = pair
                .ligand_rows
                .iter()
                .filter_map(|&g| {
                    if g < g2local.len() {
                        Some(g2local[g])
                    } else {
                        None
                    }
                })
                .filter(|&x| x != usize::MAX)
                .collect();
            let rec_local: Vec<usize> = pair
                .receptor_rows
                .iter()
                .filter_map(|&g| {
                    if g < g2local.len() {
                        Some(g2local[g])
                    } else {
                        None
                    }
                })
                .filter(|&x| x != usize::MAX)
                .collect();
            if lig_local.is_empty() || rec_local.is_empty() {
                return Vec::new();
            }
            let mut out = Vec::new();
            for lin in 0..n_lineages {
                for bin in 0..n_bins {
                    let lm: Vec<f32> = (0..n_types)
                        .map(|c| {
                            lig_local
                                .iter()
                                .map(|&r| get(r, c, lin, bin).0)
                                .fold(f32::INFINITY, f32::min)
                        })
                        .collect();
                    let lp: Vec<f32> = (0..n_types)
                        .map(|c| {
                            lig_local
                                .iter()
                                .map(|&r| get(r, c, lin, bin).1)
                                .fold(f32::INFINITY, f32::min)
                        })
                        .collect();
                    let rm: Vec<f32> = (0..n_types)
                        .map(|c| {
                            rec_local
                                .iter()
                                .map(|&r| get(r, c, lin, bin).0)
                                .fold(f32::INFINITY, f32::min)
                        })
                        .collect();
                    let rp: Vec<f32> = (0..n_types)
                        .map(|c| {
                            rec_local
                                .iter()
                                .map(|&r| get(r, c, lin, bin).1)
                                .fold(f32::INFINITY, f32::min)
                        })
                        .collect();
                    for s in 0..n_types {
                        if lp[s] < min_pct {
                            continue;
                        }
                        for t in 0..n_types {
                            if rp[t] < min_pct {
                                continue;
                            }
                            let x = lm[s];
                            let y = rm[t];
                            if x > 0.0 && y > 0.0 {
                                let geo = (x * y).sqrt();
                                let s_idx = s * slb + lin * sb + bin;
                                let t_idx = t * slb + lin * sb + bin;
                                let ws = wsums[s_idx] as f32;
                                let wt = wsums[t_idx] as f32;
                                let ess_s = if w2sums[s_idx] > 1e-12 {
                                    (wsums[s_idx] * wsums[s_idx] / w2sums[s_idx]) as f32
                                } else {
                                    0.0
                                };
                                let ess_t = if w2sums[t_idx] > 1e-12 {
                                    (wsums[t_idx] * wsums[t_idx] / w2sums[t_idx]) as f32
                                } else {
                                    0.0
                                };
                                let hs = if wsums[s_idx] > 1e-8 {
                                    (hwsums[s_idx] / wsums[s_idx]) as f32
                                } else {
                                    0.0
                                };
                                let ht = if wsums[t_idx] > 1e-8 {
                                    (hwsums[t_idx] / wsums[t_idx]) as f32
                                } else {
                                    0.0
                                };
                                out.push(LrBinResult {
                                    pair_idx: p,
                                    source: s,
                                    target: t,
                                    lineage: lin,
                                    bin,
                                    score: geo / (hill_s_bins + geo), // Michaelis-Menten saturation
                                    n_eff_cells: ess_s.min(ess_t),
                                    weight_mass_cells: ws.min(wt),
                                    lineage_entropy: hs.max(ht),
                                });
                            }
                        }
                    }
                }
            }
            out
        })
        .collect();

    // per-lineage entropy: `hwsums` / `wsums` borrow ends here (rayon closure above
    // captured `&hwsums` / `&wsums` for line 614-617's read-only lookup).
    // Move them both into the intermediate now that the parallel region is
    // done — zero-copy, layout unchanged (`ct * slb + lin * sb + bin`).
    let intermediate = LrBinIntermediate {
        hwsums: hwsums.clone(),
        wsums: wsums.clone(),
        n_types,
        n_lineages,
        n_bins,
    };
    let grid = LrBinGrid {
        means,
        pcts,
        wsums,
        w2sums,
        g2local,
        bin_bounds,
        hill_s: hill_s_bins,
        n_types,
        n_lineages,
        n_bins,
    };
    (scores, intermediate, hill_s_bins, grid)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn weighted_bin_bounds_do_not_split_tied_pseudotime() {
        let pseudotime = vec![0.0, 0.2, 0.2, 0.2, 0.8, 1.0];
        let weights = vec![1.0; pseudotime.len()];
        let bounds = compute_bin_bounds(pseudotime.len(), &pseudotime, &weights, 1, 3);

        assert_eq!(bounds.len(), 1);
        assert!(bounds[0].windows(2).all(|pair| pair[0] < pair[1]));
        assert!(bounds[0].iter().all(|&bound| bound != 0.2));
    }

    #[test]
    fn weighted_bin_bounds_are_strict_when_unique_values_equal_bins() {
        let pseudotime = vec![0.0, 0.25, 0.5, 0.75];
        let weights = vec![100.0, 1.0, 1.0, 1.0];
        let bounds = compute_bin_bounds(pseudotime.len(), &pseudotime, &weights, 1, 4);

        assert!(bounds[0].windows(2).all(|pair| pair[0] < pair[1]));
    }

    #[test]
    fn common_bin_bounds_share_one_boundary_row_across_lineages() {
        // Disjoint lineage support: lineage 0 occupies early pseudotime,
        // lineage 1 occupies late pseudotime.  Per-lineage grids disagree,
        // but the common grid must return one identical boundary row.
        let pseudotime = vec![0.1, 0.2, 0.3, 0.5, 0.7, 0.9];
        let weights = vec![
            1.0, 0.0, 1.0, 0.0, 1.0, 0.0, // lineage 0 cells
            0.0, 1.0, 0.0, 1.0, 0.0, 1.0, // lineage 1 cells
        ];
        let common = compute_common_bin_bounds(6, &pseudotime, &weights, 2, 2);
        let per = compute_bin_bounds(6, &pseudotime, &weights, 2, 2);

        assert_eq!(common.len(), 2);
        assert_eq!(common[0], common[1]);
        assert_eq!(common[0].first(), Some(&0.0));
        assert_eq!(common[0].last(), Some(&1.0));
        assert_ne!(per[0], per[1]);
        // Pooled equal-frequency cut: 6 unit-mass cells over 2 bins → the
        // shared boundary sits between cells 2 and 3.
        assert!((common[0][1] - 0.4).abs() < 1e-6);
    }

    #[test]
    fn common_bin_bounds_pools_max_lineage_weight_and_skips_nan() {
        // A cell's pooled mass is its maximum lineage weight; NaN cells are
        // excluded entirely.
        let pseudotime = vec![0.1, 0.9, f32::NAN, 0.5];
        let weights = vec![
            1.0, 0.2, // cell 0: max = 1.0
            0.1, 1.0, // cell 1: max = 1.0
            1.0, 1.0, // cell 2: NaN pseudotime → skipped
            0.4, 0.4, // cell 3: max = 0.4
        ];
        let common = compute_common_bin_bounds(4, &pseudotime, &weights, 2, 2);
        assert_eq!(common[0], common[1]);
        assert!(common[0][0] == 0.0 && common[0][2] == 1.0);
        // total mass = 1.0 + 0.4 + 1.0 = 2.4; target per bin = 1.2.  The
        // first cut index is 1 (cumulative 1.0 < 1.2, next 1.4 >= 1.2), so
        // the boundary is the midpoint of the 0.5 and 0.9 cells.
        assert!((common[0][1] - 0.7).abs() < 1e-6);
    }

    #[test]
    fn common_grid_scoring_reuses_identical_bounds_for_every_lineage() {
        // Two genes × four cells, one LR pair, one cell type, two lineages
        // with disjoint pseudotime support.
        let expr = SparseMatrix::new_csc(
            vec![0, 1, 0, 1, 0, 1, 0, 1],
            vec![0, 2, 4, 6, 8],
            vec![1.0; 8],
            2,
            4,
        );
        let cell_types = vec![0usize; 4];
        let pseudotime = vec![0.1f32, 0.3, 0.7, 0.9];
        let weights = vec![
            1.0, 0.0, 1.0, 0.0, // lineage 0: cells 0,1
            0.0, 1.0, 0.0, 1.0, // lineage 1: cells 2,3
        ];
        let pairs = vec![LrPair {
            ligand_rows: vec![0],
            receptor_rows: vec![1],
        }];

        let (_, _, _, grid_common) = compute_lr_scores_bins_with_scale_grid_opts(
            &expr, &cell_types, 1, &pairs, &pseudotime, &weights, 2, 2, 0.0, true,
        );
        assert_eq!(grid_common.bin_bounds[0], grid_common.bin_bounds[1]);

        let (_, _, _, grid_local) = compute_lr_scores_bins_with_scale_grid_opts(
            &expr, &cell_types, 1, &pairs, &pseudotime, &weights, 2, 2, 0.0, false,
        );
        assert_ne!(grid_local.bin_bounds[0], grid_local.bin_bounds[1]);
    }

    #[test]
    fn detection_fraction_counts_only_strictly_positive_values() {
        let expr = SparseMatrix::new_csc(
            vec![0, 1, 0, 0],
            vec![0, 2, 3, 4],
            vec![0.0, -2.0, 1.0, 2.0],
            2,
            3,
        );
        let (means, pcts) = cluster_stats_csc(&expr, &[0, 0, 0], 1);
        assert!((means[0][0] - 1.0).abs() < 1e-6);
        assert!((pcts[0][0] - 2.0 / 3.0).abs() < 1e-6);
        assert!((means[1][0] + 2.0 / 3.0).abs() < 1e-6);
        assert_eq!(pcts[1][0], 0.0);
    }

    #[test]
    fn selected_cluster_stats_match_full_rows_exactly() {
        let expr = SparseMatrix::new_csc(
            vec![0, 1, 3, 0, 2, 3, 1, 2],
            vec![0, 3, 6, 8],
            vec![1.0, 2.0, 4.0, 3.0, 5.0, 6.0, 7.0, 8.0],
            4,
            3,
        );
        let cell_types = vec![0, 1, 0];
        let selected = vec![0usize, 2usize];
        let mut g2local = vec![usize::MAX; 4];
        g2local[0] = 0;
        g2local[2] = 1;

        let (full_means, full_pcts) = cluster_stats_csc(&expr, &cell_types, 2);
        let (selected_means, selected_pcts) =
            cluster_stats_csc_selected(&expr, &cell_types, 2, &selected, &g2local);

        assert_eq!(selected_means[0], full_means[0]);
        assert_eq!(selected_means[1], full_means[2]);
        assert_eq!(selected_pcts[0], full_pcts[0]);
        assert_eq!(selected_pcts[1], full_pcts[2]);
    }

    #[test]
    fn selected_permutation_path_matches_full_reference_scores_and_pvalues() {
        let expr = SparseMatrix::new_csc(
            vec![0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3],
            vec![0, 4, 8, 12, 16],
            vec![
                1.0, 2.0, 3.0, 4.0, 2.0, 1.0, 4.0, 3.0, 4.0, 3.0, 2.0, 1.0, 3.0, 4.0, 1.0, 2.0,
            ],
            4,
            4,
        );
        let cell_types = vec![0, 0, 1, 1];
        let lr_pairs = vec![
            LrPair {
                ligand_rows: vec![0],
                receptor_rows: vec![1],
            },
            LrPair {
                ligand_rows: vec![2],
                receptor_rows: vec![3],
            },
        ];
        let n_types = 2usize;
        let n_perms = 7usize;
        let seed = 42u64;
        let min_pct = 0.0f32;
        let observed = compute_lr_scores(
            &expr,
            &cell_types,
            n_types,
            &lr_pairs,
            n_perms,
            seed,
            min_pct,
        );

        let (means, pcts) = cluster_stats_csc(&expr, &cell_types, n_types);
        let hill_s = compute_hill_s(&means, n_types);
        let nt2 = n_types * n_types;
        let real_data: Vec<(Vec<f32>, Vec<f32>, Vec<f32>)> = lr_pairs
            .iter()
            .map(|pair| {
                let (lm, lp) = complex_stats(&pair.ligand_rows, &means, &pcts, n_types);
                let (rm, rp) = complex_stats(&pair.receptor_rows, &means, &pcts, n_types);
                (
                    lr_score_matrix(&lm, &rm, &lp, &rp, n_types, min_pct, hill_s),
                    lp,
                    rp,
                )
            })
            .collect();
        let real: Vec<f32> = real_data
            .iter()
            .flat_map(|(score, _, _)| score.iter().copied())
            .collect();
        let mut perm_counts = vec![0u32; lr_pairs.len() * nt2];
        for pi in 0..n_perms {
            let mut rng = StdRng::seed_from_u64(seed.wrapping_add(pi as u64));
            let mut shuf = cell_types.clone();
            shuf.shuffle(&mut rng);
            let (pm, pp) = cluster_stats_csc(&expr, &shuf, n_types);
            for (p, pair) in lr_pairs.iter().enumerate() {
                let (lm, lp) = complex_stats(&pair.ligand_rows, &pm, &pp, n_types);
                let (rm, rp) = complex_stats(&pair.receptor_rows, &pm, &pp, n_types);
                let scores = lr_score_matrix(&lm, &rm, &lp, &rp, n_types, min_pct, hill_s);
                for (k, (&score, &reference)) in
                    scores.iter().zip(real[p * nt2..].iter()).enumerate()
                {
                    if score >= reference {
                        perm_counts[p * nt2 + k] += 1;
                    }
                }
            }
        }

        let mut expected = Vec::new();
        for p in 0..lr_pairs.len() {
            let (_, ref lp, ref rp) = real_data[p];
            for source in 0..n_types {
                for target in 0..n_types {
                    let k = p * nt2 + source * n_types + target;
                    if real[k] > 0.0 {
                        expected.push((
                            p,
                            source,
                            target,
                            real[k],
                            (perm_counts[k] as f32 + 1.0) / (n_perms as f32 + 1.0),
                            lp[source],
                            rp[target],
                        ));
                    }
                }
            }
        }

        assert_eq!(observed.len(), expected.len());
        for (row, reference) in observed.iter().zip(expected.iter()) {
            assert_eq!(row.pair_idx, reference.0);
            assert_eq!(row.source, reference.1);
            assert_eq!(row.target, reference.2);
            assert_eq!(row.score, reference.3);
            assert_eq!(row.pvalue, reference.4);
            assert_eq!(row.pct_ligand, reference.5);
            assert_eq!(row.pct_receptor, reference.6);
        }
    }

    /// `compute_hill_s` uses the arithmetic mean over non-zero `means[g][ct]`
    /// entries with `ct < n_types`, matching the SingleCellSignalR scale
    /// (Cabello-Aguilar 2020, Bioinformatics).
    #[test]
    fn hill_s_lock_simple() {
        // means[g][ct]; non-zero values across (g, ct) with ct < 2:
        //   (0,0)=1.0, (0,1)=2.0, (1,1)=4.0, (2,0)=3.0
        // sum = 10.0, count = 4 → mean = 2.5
        let means = vec![vec![1.0, 2.0], vec![0.0, 4.0], vec![3.0, 0.0]];
        let s = compute_hill_s(&means, 2);
        assert!((s - 2.5).abs() < 1e-6, "expected 2.5, got {}", s);
    }

    /// Boundary case: all-zero means must fall back to 1.0 to avoid division
    /// by zero in `score = sqrt(L·R) / (S + sqrt(L·R))`.
    #[test]
    fn hill_s_lock_all_zero() {
        let means = vec![vec![0.0, 0.0], vec![0.0, 0.0]];
        let s = compute_hill_s(&means, 2);
        assert_eq!(s, 1.0);
    }

    /// `n_types` truncation: values beyond `means[g][..n_types]` must be
    /// excluded. Locks the `.take(n_types)` slicing behaviour.
    #[test]
    fn hill_s_lock_partial_n_types() {
        // n_types = 2; means[g][2] = 999.0 must be ignored.
        let means = vec![vec![1.0, 2.0, 999.0]];
        let s = compute_hill_s(&means, 2);
        // Only 1.0 and 2.0 count → mean = 1.5
        assert!(
            (s - 1.5).abs() < 1e-6,
            "expected 1.5 (must ignore col >= n_types), got {}",
            s
        );
    }

    #[test]
    fn temporal_hill_scale_is_invariant_to_lr_subset() {
        let expr = SparseMatrix::new_csc(
            vec![0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3],
            vec![0, 4, 8, 12, 16],
            vec![
                1.0, 2.0, 3.0, 4.0, 1.0, 2.0, 3.0, 4.0, 2.0, 1.0, 4.0, 3.0, 2.0, 1.0, 4.0, 3.0,
            ],
            4,
            4,
        );
        let cell_types = vec![0, 0, 1, 1];
        let pseudotime = vec![0.0, 0.3, 0.7, 1.0];
        let lineage_weights = vec![1.0; 4];
        let primary = LrPair {
            ligand_rows: vec![0],
            receptor_rows: vec![1],
        };
        let unrelated = LrPair {
            ligand_rows: vec![2],
            receptor_rows: vec![3],
        };

        let (scores_one, _, s_one) = compute_lr_scores_bins_with_scale(
            &expr,
            &cell_types,
            2,
            std::slice::from_ref(&primary),
            &pseudotime,
            &lineage_weights,
            1,
            2,
            0.0,
        );
        let (scores_two, _, s_two) = compute_lr_scores_bins_with_scale(
            &expr,
            &cell_types,
            2,
            &[primary, unrelated],
            &pseudotime,
            &lineage_weights,
            1,
            2,
            0.0,
        );
        let first_pair: Vec<&LrBinResult> =
            scores_two.iter().filter(|row| row.pair_idx == 0).collect();

        assert!((s_one - s_two).abs() < 1e-6);
        assert_eq!(scores_one.len(), first_pair.len());
        assert!(scores_one
            .iter()
            .zip(first_pair)
            .all(|(left, right)| (left.score - right.score).abs() < 1e-6));
    }
}
