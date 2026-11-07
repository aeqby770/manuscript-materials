//! Pseudotime, lineage weights, and fate probabilities.
//!
//! ## Comparison with Slingshot (Street et al., 2018)
//!
//! Slingshot constructs trajectories via:
//!   1. MST on cluster centroids (same as this impl.)
//!   2. **Simultaneous principal curves** fitted iteratively in the ambient
//!      embedding space (`princurve::principal_curve`, smoother = "lowess")
//!   3. `shrink.common`: pre-bifurcation cells are pulled toward the shared
//!      portion of diverging lineages, preventing double-counting of early cells
//!   4. Cell pseudotime = arc-length projection onto each principal curve
//!
//! This implementation simplifies to:
//!   1. MST on cluster centroids (step 1 identical)
//!   2. **Piecewise-linear path through centroids** — no smooth curve fitting
//!   3. Gaussian soft-assignment to MST branches (no shrink.common)
//!   4. Dijkstra pseudotime on cell-level kNN graph (fuzzy weights option)
//!
//! Practical impact: for well-separated, linear trajectories the difference is
//! minor.  For curved or branching trajectories in high dimensions, principal
//! curves fit the data manifold more accurately.  Arc-length pseudotime from
//! Slingshot is smoother than Dijkstra pseudotime near branch points.
//!
//! ## Comparison with Palantir (Setty et al., 2019)
//!
//! Palantir uses:
//!   1. **Diffusion maps** to construct a low-dimensional Markov chain that
//!      respects the data manifold (not a kNN graph with Euclidean weights)
//!   2. Multi-scale space from diffusion eigenvectors (entropy-based terminal
//!      state identification)
//!   3. **Probabilistic pseudotime** via harmonic function on Markov chain
//!   4. Fate probabilities via absorbing Markov chain (same concept, different
//!      graph construction)
//!
//! This implementation:
//!   1. Uses Dijkstra on a kNN L2-distance graph, optionally after converting
//!      distances to locally scaled fuzzy memberships
//!   2. Identifies terminal states by user-supplied cluster labels
//!   3. Pseudotime = shortest path from root cell
//!   4. Fate probabilities via sigmoid-gated directed Markov chain (power iter.)
//!
//! Locally scaled fuzzy memberships reduce the influence of low-confidence
//! kNN edges before shortest-path distances are computed.

use rayon::prelude::*;
use std::cmp::Ordering;
use std::collections::{BinaryHeap, HashMap, HashSet, VecDeque};

// ─── Dijkstra helpers ────────────────────────────────────────────────────────

#[derive(Clone)]
struct HeapNode {
    cost: f32,
    node: usize,
}
impl PartialEq for HeapNode {
    fn eq(&self, o: &Self) -> bool {
        self.cost == o.cost
    }
}
impl Eq for HeapNode {}
impl PartialOrd for HeapNode {
    fn partial_cmp(&self, o: &Self) -> Option<Ordering> {
        Some(self.cmp(o))
    }
}
impl Ord for HeapNode {
    // min-heap via reversed compare
    fn cmp(&self, o: &Self) -> Ordering {
        o.cost.partial_cmp(&self.cost).unwrap_or(Ordering::Equal)
    }
}

/// Dijkstra shortest-path distances from `start` on an undirected kNN graph.
///
/// `neighbors`: for each cell, list of (neighbor_idx, distance)
/// Returns dist[cell], f32::INFINITY if unreachable.
pub fn dijkstra(n_cells: usize, neighbors: &[Vec<(usize, f32)>], start: usize) -> Vec<f32> {
    let mut dist = vec![f32::INFINITY; n_cells];
    dist[start] = 0.0;
    let mut heap = BinaryHeap::new();
    heap.push(HeapNode {
        cost: 0.0,
        node: start,
    });

    while let Some(HeapNode { cost, node }) = heap.pop() {
        if cost > dist[node] {
            continue;
        }
        for &(nb, w) in &neighbors[node] {
            let next = cost + w;
            if next < dist[nb] {
                dist[nb] = next;
                heap.push(HeapNode {
                    cost: next,
                    node: nb,
                });
            }
        }
    }
    dist
}

#[derive(Clone)]
struct PreciseHeapNode {
    cost: f64,
    node: usize,
}
impl PartialEq for PreciseHeapNode {
    fn eq(&self, other: &Self) -> bool {
        self.cost == other.cost
    }
}
impl Eq for PreciseHeapNode {}
impl PartialOrd for PreciseHeapNode {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for PreciseHeapNode {
    fn cmp(&self, other: &Self) -> Ordering {
        other.cost.partial_cmp(&self.cost).unwrap_or(Ordering::Equal)
    }
}

fn dijkstra_precise(
    n_cells: usize,
    neighbors: &[Vec<(usize, f32)>],
    start: usize,
) -> Vec<f64> {
    let mut dist = vec![f64::INFINITY; n_cells];
    dist[start] = 0.0;
    let mut heap = BinaryHeap::new();
    heap.push(PreciseHeapNode { cost: 0.0, node: start });

    while let Some(PreciseHeapNode { cost, node }) = heap.pop() {
        if cost > dist[node] {
            continue;
        }
        for &(neighbor, edge_cost) in &neighbors[node] {
            let next = cost + f64::from(edge_cost);
            if next < dist[neighbor] {
                dist[neighbor] = next;
                heap.push(PreciseHeapNode { cost: next, node: neighbor });
            }
        }
    }
    dist
}

/// Compute per-cell normalized lineage entropy H(Ψ)/log(L), in [0,1].
///
/// - 0   => lineage-certain cells (weight concentrated in one lineage)
/// - 1   => maximally ambiguous cells (uniform weights across lineages)
///
/// This provides a direct, quantitative flag for bifurcation ambiguity and is
/// intended to operationalize the `H(Ψ)` caveat discussed in LineaLinkR Module 1.
pub fn lineage_entropy(lineage_weights: &[f32], n_cells: usize, n_lineages: usize) -> Vec<f32> {
    if n_lineages <= 1 {
        return vec![0.0; n_cells];
    }
    let ln_l = (n_lineages as f32).ln().max(1e-6);
    (0..n_cells)
        .into_par_iter()
        .map(|i| {
            let row = &lineage_weights[i * n_lineages..(i + 1) * n_lineages];
            let z: f32 = row
                .iter()
                .copied()
                .filter(|w| w.is_finite() && *w > 0.0)
                .sum();
            if z <= 0.0 {
                return 0.0;
            }
            let mut h = 0.0f32;
            for &w in row {
                if w.is_finite() && w > 0.0 {
                    let p = (w / z).max(1e-12);
                    h -= p * p.ln();
                }
            }
            (h / ln_l).clamp(0.0, 1.0)
        })
        .collect()
}

/// Lightweight DTW-based trajectory profile similarity in [0,1].
///
/// This absorbs TraSig's key advantage (cross-trajectory temporal alignment)
/// while remaining computationally cheap for per-pair profile comparison.
///
/// Similarity is computed as:
///   sim = 1 / (1 + dist / max(n, m))
/// where `dist` is absolute-difference DTW cost with Sakoe-Chiba band `window`.
pub fn dtw_similarity(seq_a: &[f32], seq_b: &[f32], window: usize) -> f32 {
    let n = seq_a.len();
    let m = seq_b.len();
    if n == 0 || m == 0 {
        return 0.0;
    }

    let band = window.max(n.abs_diff(m));
    let mut prev = vec![f32::INFINITY; m + 1];
    let mut curr = vec![f32::INFINITY; m + 1];
    prev[0] = 0.0;

    for i in 1..=n {
        curr.fill(f32::INFINITY);
        let j_lo = 1usize.max(i.saturating_sub(band));
        let j_hi = m.min(i + band);
        for j in j_lo..=j_hi {
            let cost = (seq_a[i - 1] - seq_b[j - 1]).abs();
            let best = prev[j].min(curr[j - 1]).min(prev[j - 1]);
            curr[j] = cost + best;
        }
        std::mem::swap(&mut prev, &mut curr);
    }

    let dist = prev[m];
    if !dist.is_finite() {
        return 0.0;
    }
    let norm = dist / (n.max(m) as f32);
    (1.0 / (1.0 + norm)).clamp(0.0, 1.0)
}

// ─── kNN graph from dense embedding ─────────────────────────────────────────

/// Build approximate kNN graph from embedding matrix (n_cells × n_dims).
///
/// **HARD LIMIT**: n_cells must not exceed 5 000.  This function is O(n² × d)
/// and is intended only for unit tests and small exploratory datasets.
/// For production use, supply HNSW-computed indices via `knn_graph_from_hnsw`.
///
/// Panics immediately with a descriptive message if the limit is exceeded,
/// rather than silently hanging for hours on atlas-scale data.
///
/// Returns neighbors[cell] = Vec<(neighbor_idx, l2_dist)>
pub fn knn_graph_brute(
    embedding: &[f32], // row-major [n_cells × n_dims]
    n_cells: usize,
    n_dims: usize,
    k: usize,
) -> Vec<Vec<(usize, f32)>> {
    assert!(
        n_cells <= 5_000,
        "knn_graph_brute is O(n²×d) and must not be called with n_cells={} > 5000. \
         At this scale it would take ~{:.1} minutes on a 16-core machine. \
         Pass pre-computed HNSW indices via knn_graph_from_hnsw instead.",
        n_cells,
        (n_cells as f64).powi(2) * n_dims as f64 / (16.0 * 2e9) / 60.0
    );
    (0..n_cells)
        .into_par_iter()
        .map(|i| {
            let xi = &embedding[i * n_dims..(i + 1) * n_dims];
            let mut dists: Vec<(usize, f32)> = (0..n_cells)
                .filter(|&j| j != i)
                .map(|j| {
                    let xj = &embedding[j * n_dims..(j + 1) * n_dims];
                    let d = xi
                        .iter()
                        .zip(xj)
                        .map(|(a, b)| (a - b).powi(2))
                        .sum::<f32>()
                        .sqrt();
                    (j, d)
                })
                .collect();
            dists.sort_unstable_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
            dists.truncate(k);
            dists
        })
        .collect()
}

/// Build kNN graph from pre-computed HNSW neighbor indices + distances.
/// `nn_idx`: flat [n_cells × k], 0-based
/// `nn_dist`: flat [n_cells × k]
pub fn knn_graph_from_hnsw(
    nn_idx: &[usize],
    nn_dist: &[f32],
    n_cells: usize,
    k: usize,
) -> Vec<Vec<(usize, f32)>> {
    let directed: Vec<Vec<(usize, f32)>> = (0..n_cells)
        .map(|i| {
            (0..k)
                .map(|j| (nn_idx[i * k + j], nn_dist[i * k + j]))
                .collect()
        })
        .collect();
    symmetrize_knn_graph(directed, n_cells)
}

/// Inject TICCI-style cell-cell communication (CCI) confidence into kNN edge lengths.
///
/// For each directed edge i -> j in the HNSW neighbour list, if cell pair (i,j)
/// appears in the CCI table with probability p, the edge distance becomes:
///
///   d'_ij = d_ij / (1 + cci_k * p_ij)
///
/// where p_ij is the summed probability across duplicated entries for the same
/// unordered cell pair.
///
/// This is equivalent to TICCI's idea of strengthening communication-supported
/// neighbourhood links before trajectory inference, while preserving the original
/// sparse HNSW graph structure and O(n_cells * k) complexity.
pub fn inject_cci_weights_into_knn(
    nn_idx: &[usize],
    nn_dist: &[f32],
    n_cells: usize,
    k: usize,
    sender_cells: &[usize],
    receiver_cells: &[usize],
    probs: &[f32],
    cci_k: f32,
) -> Vec<f32> {
    let mut adjusted = nn_dist.to_vec();
    if cci_k <= 0.0 {
        return adjusted;
    }

    let n_pairs = sender_cells
        .len()
        .min(receiver_cells.len())
        .min(probs.len());
    if n_pairs == 0 {
        return adjusted;
    }

    // Hash unordered pair (min, max) into a compact u64 key.
    let mut pair_bonus: HashMap<u64, f32> = HashMap::with_capacity(n_pairs * 2 + 1);
    for x in 0..n_pairs {
        let s = sender_cells[x];
        let r = receiver_cells[x];
        let p = probs[x];
        if s >= n_cells || r >= n_cells || s == r || !p.is_finite() || p <= 0.0 {
            continue;
        }
        let (a, b) = if s < r { (s, r) } else { (r, s) };
        let key = ((a as u64) << 32) | (b as u64);
        *pair_bonus.entry(key).or_insert(0.0) += p;
    }

    if pair_bonus.is_empty() {
        return adjusted;
    }

    let cci_scale = cci_k.max(0.0);
    for i in 0..n_cells {
        let base = i * k;
        for j in 0..k {
            let nb = nn_idx[base + j];
            if nb >= n_cells || nb == i {
                continue;
            }
            let (a, b) = if i < nb { (i, nb) } else { (nb, i) };
            let key = ((a as u64) << 32) | (b as u64);
            if let Some(&p) = pair_bonus.get(&key) {
                let factor = 1.0 + cci_scale * p.max(0.0);
                adjusted[base + j] = (nn_dist[base + j] / factor).max(1e-6);
            }
        }
    }

    adjusted
}

fn symmetrize_knn_graph(
    directed: Vec<Vec<(usize, f32)>>,
    n_cells: usize,
) -> Vec<Vec<(usize, f32)>> {
    let mut undirected: Vec<Vec<(usize, f32)>> = vec![Vec::new(); n_cells];

    for i in 0..n_cells {
        for &(j, w) in &directed[i] {
            if j >= n_cells || j == i || !w.is_finite() || w < 0.0 {
                continue;
            }
            undirected[i].push((j, w));
            undirected[j].push((i, w));
        }
    }

    for row in undirected.iter_mut() {
        if row.len() <= 1 {
            continue;
        }
        row.sort_unstable_by_key(|&(j, _)| j);
        let mut write = 0usize;
        for read in 0..row.len() {
            if write == 0 || row[read].0 != row[write - 1].0 {
                row[write] = row[read];
                write += 1;
            } else if row[read].1 < row[write - 1].1 {
                row[write - 1].1 = row[read].1;
            }
        }
        row.truncate(write);
    }

    undirected
}

/// Build kNN graph from UMAP fuzzy simplicial set edge weights.
///
/// Edge cost = −log(max(w, ε_local)), mapping high-membership edges to low cost
/// and low-membership edges to higher cost.
///
/// ## Effect relative to raw Euclidean distances
///
/// The transformation gives low-membership edges a larger cost and
/// high-membership edges a smaller cost. It depends only on the supplied fuzzy
/// weights and is not equivalent to diffusion-map or Markov-chain pseudotime.
///
/// ## Adaptive local floor
///
/// The floor prevents −log(w) from diverging when a supplied membership is zero.
///
/// The per-cell floor is ε_local = max(1e-6, 1e-3 × mean_w_cell):
///   - Scales with the local edge-weight distribution
///   - Bounds the largest finite edge cost at −log(ε_local)
///   - Negligible effect when all weights are large (normal dense regions)
///
/// `nn_idx`:        flat [n_cells × k], 0-based neighbor indices
/// `fuzzy_weights`: flat [n_cells × k], UMAP fuzzy membership strengths ∈ (0,1]
pub fn knn_graph_from_fuzzy(
    nn_idx: &[usize],
    fuzzy_weights: &[f32],
    n_cells: usize,
    k: usize,
) -> Vec<Vec<(usize, f32)>> {
    let directed: Vec<Vec<(usize, f32)>> = (0..n_cells)
        .map(|i| {
            let base = i * k;
            // Per-cell adaptive floor: 1e-3 × local mean weight, hard minimum 1e-6
            let mean_w: f32 = fuzzy_weights[base..base + k].iter().sum::<f32>() / k as f32;
            let eps_local: f32 = (1e-3 * mean_w).max(1e-6);
            (0..k)
                .map(|j| {
                    let nb = nn_idx[base + j];
                    let w = fuzzy_weights[base + j].max(eps_local);
                    let cost = -(w.ln()); // ∈ (0, +∞); higher-confidence edges have lower cost
                    (nb, cost)
                })
                .collect()
        })
        .collect();
    symmetrize_knn_graph(directed, n_cells)
}

/// Per-cell adaptive sigma for UMAP fuzzy membership (smooth_knn_dist, McInnes 2018 §2.2).
///
/// For cell i with sorted kNN distances `dists[0..k]`:
///   ρ_i = dists[0]  (nearest-neighbour distance)
///   σ_i is solved by binary search so that:
///     Σ_{j=1..k} exp(-(dists[j] - ρ_i).max(0) / σ_i) ≈ log2(k)
///
/// Returns (rho, sigma).  Minimum sigma is clamped to 1e-6.
fn smooth_knn_dist(dists: &[f32], k: usize) -> (f32, f32) {
    if dists.is_empty() {
        return (0.0, 1e-6);
    }

    let target = (k.max(2) as f32).log2().max(1.0);
    let rho = dists[0];
    let mut lo = 0.0f32;
    let mut hi = f32::INFINITY;
    let mut mid = 1.0f32;
    let mut sigma_best = mid;
    let mut adiff_min = f32::INFINITY;
    let mut converged = false;

    for _ in 0..64 {
        let sum: f32 = dists[1..]
            .iter()
            .map(|&d| (-(d - rho).max(0.0) / mid).exp())
            .sum();
        let adiff = (sum - target).abs();
        if adiff < 1e-5 {
            converged = true;
            break;
        }
        if adiff < adiff_min {
            adiff_min = adiff;
            sigma_best = mid;
        }
        if sum > target {
            hi = mid;
            mid = (lo + hi) / 2.0;
        } else {
            lo = mid;
            mid = if hi.is_infinite() {
                mid * 2.0
            } else {
                (lo + hi) / 2.0
            };
        }
    }

    if !converged {
        mid = sigma_best;
    }

    let mean_i = dists.iter().sum::<f32>() / dists.len() as f32;
    (rho, mid.max(1e-3 * mean_i).max(1e-6))
}

/// Convert raw HNSW kNN distances to UMAP fuzzy membership weights, then
/// compute pseudotime via Dijkstra on the fuzzy graph.
///
/// The local scaling penalizes low-membership kNN edges before shortest-path
/// distances are computed. This is not a diffusion-map pseudotime.
///
/// `nn_idx`:   flat [n_cells × k], 0-based
/// `nn_dist`:  flat [n_cells × k], L2 distances from HNSW
pub fn compute_pseudotime_from_hnsw_fuzzy(
    nn_idx: &[usize],
    nn_dist: &[f32],
    n_cells: usize,
    k: usize,
    start_cell: usize,
) -> Vec<f64> {
    // Edge cost = d / w  (distance weighted by inverse fuzzy confidence)
    //
    // Use d / w rather than -ln(w):
    //   UMAP smooth_knn_dist sets w[j=0] = exp(0) = 1.0 for the nearest neighbour
    //   (by definition: rho = dists[0], so (d - rho).max(0) = 0 for j=0).
    //   This gives -ln(1) = 0 for all nearest-neighbour edges.
    //   When every cell's nearest-neighbour edge has cost 0, Dijkstra propagates
    //   with total cost 0 from the root to every cell, and pseudotime ≡ 0.
    //
    // d / w semantics:
    //   - Nearest neighbour (w=1.0): cost = d  (raw L2 distance, no penalty)
    //   - Subsequent neighbours (w<1): cost = d/w > d  (slightly penalised)
    //   - Low-membership edges (w near 0): cost becomes large
    //
    // Dijkstra therefore runs on locally rescaled kNN edge costs.
    let neighbors = hnsw_fuzzy_graph(nn_idx, nn_dist, n_cells, k);
    normalized_dijkstra_precise(n_cells, &neighbors, start_cell)
}

/// Compute several root-specific pseudotime vectors while constructing the
/// UMAP-fuzzy graph only once. Results are flattened root-major.
pub fn compute_pseudotime_from_hnsw_fuzzy_multi(
    nn_idx: &[usize],
    nn_dist: &[f32],
    n_cells: usize,
    k: usize,
    start_cells: &[usize],
) -> Vec<f32> {
    let neighbors = hnsw_fuzzy_graph(nn_idx, nn_dist, n_cells, k);
    start_cells
        .iter()
        .flat_map(|&start_cell| normalized_dijkstra(n_cells, &neighbors, start_cell))
        .collect()
}

fn hnsw_fuzzy_graph(
    nn_idx: &[usize],
    nn_dist: &[f32],
    n_cells: usize,
    k: usize,
) -> Vec<Vec<(usize, f32)>> {
    let directed: Vec<Vec<(usize, f32)>> = (0..n_cells)
        .map(|i| {
            let mut pairs: Vec<(usize, f32)> = (0..k)
                .filter_map(|j| {
                    let nb = nn_idx[i * k + j];
                    let d = nn_dist[i * k + j];
                    if nb < n_cells && nb != i && d.is_finite() && d >= 0.0 {
                        Some((nb, d))
                    } else {
                        None
                    }
                })
                .collect();
            pairs.sort_unstable_by(|a, b| {
                a.1.partial_cmp(&b.1)
                    .unwrap_or(Ordering::Equal)
                    .then_with(|| a.0.cmp(&b.0))
            });
            let mut seen = HashSet::new();
            pairs.retain(|(nb, _)| seen.insert(*nb));
            let valid_distances: Vec<f32> = pairs.iter().map(|(_, d)| *d).collect();
            let (rho, sigma) = smooth_knn_dist(&valid_distances, valid_distances.len());
            pairs
                .into_iter()
                .map(|(nb, d)| {
                    let w = (-(d - rho).max(0.0) / sigma).exp().max(1e-6);
                    let cost = d / w;
                    (nb, cost)
                })
                .collect()
        })
        .collect();
    symmetrize_knn_graph(directed, n_cells)
}

/// One undirected edge from the exact graph used by fuzzy Dijkstra together
/// with the corresponding symmetrized raw-distance cost.
#[derive(Clone, Debug)]
pub struct KnnCostEdge {
    pub from: usize,
    pub to: usize,
    pub fuzzy_cost: f32,
    pub raw_distance: f32,
}

/// Expose the realized fuzzy and raw edge costs for numerical diagnostics.
///
/// The fuzzy column is produced by `hnsw_fuzzy_graph`, including its local
/// scale, `d / w` transformation, duplicate resolution and undirected
/// symmetrization. The raw column uses the same supplied kNN graph and the same
/// symmetrization, changing only the edge cost back to `d`. Each undirected
/// pair is returned once with zero-based endpoints.
pub fn hnsw_cost_edges(
    nn_idx: &[usize],
    nn_dist: &[f32],
    n_cells: usize,
    k: usize,
) -> Vec<KnnCostEdge> {
    let fuzzy = hnsw_fuzzy_graph(nn_idx, nn_dist, n_cells, k);
    let raw = knn_graph_from_hnsw(nn_idx, nn_dist, n_cells, k);
    let mut edges = Vec::new();

    for from in 0..n_cells {
        for &(to, fuzzy_cost) in &fuzzy[from] {
            if from >= to {
                continue;
            }
            let raw_at = raw[from]
                .binary_search_by_key(&to, |&(neighbor, _)| neighbor)
                .expect("fuzzy and raw graphs must share the same topology");
            edges.push(KnnCostEdge {
                from,
                to,
                fuzzy_cost,
                raw_distance: raw[from][raw_at].1,
            });
        }
    }
    edges
}

fn normalized_dijkstra(
    n_cells: usize,
    neighbors: &[Vec<(usize, f32)>],
    start_cell: usize,
) -> Vec<f32> {
    let mut dist = dijkstra(n_cells, neighbors, start_cell);

    // Disconnected cells remain NaN; finite pseudotime normalised to [0, 1].
    let max_d = dist
        .iter()
        .cloned()
        .filter(|d| d.is_finite())
        .fold(0.0f32, f32::max);
    dist.iter_mut().for_each(|d| {
        if !d.is_finite() {
            *d = f32::NAN;
        }
    });
    if max_d > 0.0 {
        dist.iter_mut().for_each(|d| {
            if d.is_finite() {
                *d /= max_d;
            }
        });
    }
    dist
}

fn normalized_dijkstra_precise(
    n_cells: usize,
    neighbors: &[Vec<(usize, f32)>],
    start_cell: usize,
) -> Vec<f64> {
    let mut dist = dijkstra_precise(n_cells, neighbors, start_cell);
    let max_d = dist
        .iter()
        .copied()
        .filter(|value| value.is_finite())
        .fold(0.0f64, f64::max);
    dist.iter_mut().for_each(|value| {
        if !value.is_finite() {
            *value = f64::NAN;
        } else if max_d > 0.0 {
            *value /= max_d;
        }
    });
    dist
}

/// Count connected components in the undirected form of a row-major kNN graph.
pub fn count_knn_components(nn_idx: &[usize], n_cells: usize, k: usize) -> usize {
    knn_component_labels(nn_idx, n_cells, k)
        .into_iter()
        .max()
        .map_or(0, |label| label + 1)
}

/// Return zero-based connected-component labels for a row-major kNN graph.
pub fn knn_component_labels(nn_idx: &[usize], n_cells: usize, k: usize) -> Vec<usize> {
    let mut adjacency = vec![Vec::<usize>::new(); n_cells];
    for cell in 0..n_cells {
        for neighbor_rank in 0..k {
            let neighbor = nn_idx[cell * k + neighbor_rank];
            if neighbor < n_cells && neighbor != cell {
                adjacency[cell].push(neighbor);
                adjacency[neighbor].push(cell);
            }
        }
    }
    let mut seen = vec![false; n_cells];
    let mut labels = vec![usize::MAX; n_cells];
    let mut component = 0usize;
    let mut queue = VecDeque::new();
    for start in 0..n_cells {
        if seen[start] {
            continue;
        }
        seen[start] = true;
        labels[start] = component;
        queue.push_back(start);
        while let Some(cell) = queue.pop_front() {
            for &neighbor in &adjacency[cell] {
                if !seen[neighbor] {
                    seen[neighbor] = true;
                    labels[neighbor] = component;
                    queue.push_back(neighbor);
                }
            }
        }
        component += 1;
    }
    labels
}

/// Compute pseudotime using UMAP fuzzy simplicial set weights as Dijkstra costs.
///
/// Low-membership edges receive larger costs than high-membership edges. This
/// provides a fuzzy-weight alternative to the raw-distance graph used by
/// `compute_pseudotime`.
///
/// `nn_idx`:        flat [n_cells × k], 0-based
/// `fuzzy_weights`: flat [n_cells × k], UMAP fuzzy membership strengths
/// `start_cell`:    index of the graph origin
pub fn compute_pseudotime_fuzzy(
    nn_idx: &[usize],
    fuzzy_weights: &[f32],
    n_cells: usize,
    k: usize,
    start_cell: usize,
) -> Vec<f32> {
    let neighbors = knn_graph_from_fuzzy(nn_idx, fuzzy_weights, n_cells, k);
    let mut dist = dijkstra(n_cells, &neighbors, start_cell);

    // Disconnected cells → NaN (same policy as compute_pseudotime)
    let max_d = dist
        .iter()
        .cloned()
        .filter(|d| d.is_finite())
        .fold(0.0f32, f32::max);
    dist.iter_mut().for_each(|d| {
        if !d.is_finite() {
            *d = f32::NAN;
        }
    });
    if max_d > 0.0 {
        dist.iter_mut().for_each(|d| {
            if d.is_finite() {
                *d /= max_d;
            }
        });
    }
    dist
}

// ─── Public API ──────────────────────────────────────────────────────────────

/// Compute pseudotime for all cells as Dijkstra distance from start_cell,
/// normalised to [0, 1].
///
/// `nn_idx`, `nn_dist`: HNSW output (flat, row-major [n_cells × k])
pub fn compute_pseudotime(
    nn_idx: &[usize],
    nn_dist: &[f32],
    n_cells: usize,
    k: usize,
    start_cell: usize,
) -> Vec<f32> {
    let neighbors = knn_graph_from_hnsw(nn_idx, nn_dist, n_cells, k);
    let mut dist = dijkstra(n_cells, &neighbors, start_cell);

    // Disconnected cells (unreachable from start_cell) → NaN so R callers can
    // detect them via is.na(). Do NOT replace with max_d: silently treating
    // disconnected cells as "late" cells misrepresents the trajectory topology.
    // (Palantir uses connected_components to separate such cells explicitly.)
    let max_d = dist
        .iter()
        .cloned()
        .filter(|d| d.is_finite())
        .fold(0.0f32, f32::max);
    dist.iter_mut().for_each(|d| {
        if !d.is_finite() {
            *d = f32::NAN;
        }
    });

    // normalise finite values to [0, 1]; NaN cells remain NaN
    if max_d > 0.0 {
        dist.iter_mut().for_each(|d| {
            if d.is_finite() {
                *d /= max_d;
            }
        });
    }
    dist
}

// ─── Lineage weights (soft assignment via cluster centroid polylines) ─────────

/// Compute cluster centroids from embedding.
/// `cluster_labels`: per-cell cluster id (0-based)
pub fn cluster_centroids(
    embedding: &[f32],
    cluster_labels: &[usize],
    n_cells: usize,
    n_dims: usize,
    n_clusters: usize,
) -> Vec<Vec<f32>> {
    let mut sums = vec![vec![0.0f64; n_dims]; n_clusters];
    let mut counts = vec![0u32; n_clusters];
    for i in 0..n_cells {
        if i >= cluster_labels.len() {
            break;
        }
        let c = cluster_labels[i];
        if c >= n_clusters {
            continue;
        }
        counts[c] += 1;
        for d in 0..n_dims {
            if i * n_dims + d >= embedding.len() {
                break;
            }
            sums[c][d] += embedding[i * n_dims + d] as f64;
        }
    }
    sums.into_iter()
        .zip(counts.iter())
        .map(|(row, &cnt)| {
            row.into_iter()
                .map(|s| {
                    if cnt > 0 {
                        (s / cnt as f64) as f32
                    } else {
                        0.0
                    }
                })
                .collect()
        })
        .collect()
}

/// Minimal Spanning Tree via Prim's algorithm on cluster centroids.
/// Returns edges as (from, to, dist).
pub fn mst_prim(centroids: &[Vec<f32>]) -> Vec<(usize, usize, f32)> {
    let n = centroids.len();
    if n <= 1 {
        return vec![];
    }
    let dist2 = |a: &[f32], b: &[f32]| -> f32 {
        a.iter()
            .zip(b)
            .map(|(x, y)| (x - y).powi(2))
            .sum::<f32>()
            .sqrt()
    };

    let mut in_tree = vec![false; n];
    let mut min_dist = vec![f32::INFINITY; n];
    let mut parent = vec![0usize; n];
    min_dist[0] = 0.0;
    let mut edges = Vec::with_capacity(n - 1);

    for _ in 0..n {
        // pick cheapest node not in tree
        let u = (0..n)
            .filter(|&i| !in_tree[i])
            .min_by(|&a, &b| min_dist[a].partial_cmp(&min_dist[b]).unwrap())
            .unwrap();
        in_tree[u] = true;
        if u != 0 {
            edges.push((parent[u], u, min_dist[u]));
        }
        // relax neighbours
        for v in 0..n {
            if !in_tree[v] {
                let d = dist2(&centroids[u], &centroids[v]);
                if d < min_dist[v] {
                    min_dist[v] = d;
                    parent[v] = u;
                }
            }
        }
    }
    edges
}

/// Extract lineages (root→leaf paths) from MST.
/// `start_cluster`: index of root cluster.
/// Returns lineages: Vec<Vec<usize>> (each a path of cluster ids)
pub fn extract_lineages(
    mst_edges: &[(usize, usize, f32)],
    n_clusters: usize,
    start_cluster: usize,
) -> Vec<Vec<usize>> {
    // build adjacency list
    let mut adj = vec![vec![]; n_clusters];
    for &(u, v, _) in mst_edges {
        adj[u].push(v);
        adj[v].push(u);
    }

    // find leaf clusters (degree 1, not root)
    let leaves: Vec<usize> = (0..n_clusters)
        .filter(|&c| c != start_cluster && adj[c].len() == 1)
        .collect();

    // BFS once from the root, then trace each leaf back through the same parent map.
    let mut parent = vec![usize::MAX; n_clusters];
    let mut visited = vec![false; n_clusters];
    let mut queue = std::collections::VecDeque::new();
    queue.push_back(start_cluster);
    visited[start_cluster] = true;
    while let Some(cur) = queue.pop_front() {
        for &nb in &adj[cur] {
            if !visited[nb] {
                visited[nb] = true;
                parent[nb] = cur;
                queue.push_back(nb);
            }
        }
    }

    leaves
        .into_iter()
        .map(|leaf| {
            // trace path back from leaf
            let mut path = vec![];
            let mut cur = leaf;
            while cur != usize::MAX {
                path.push(cur);
                cur = parent[cur];
            }
            path.reverse();
            path
        })
        .collect()
}

#[inline]
fn point_to_lineage_sq(
    xi: &[f32],
    lineage: &[usize],
    centroids: &[Vec<f32>],
    n_dims: usize,
) -> f32 {
    match lineage.len() {
        0 => f32::INFINITY,
        1 => centroids[lineage[0]]
            .iter()
            .zip(xi)
            .map(|(a, b)| (a - b).powi(2))
            .sum::<f32>(),
        _ => lineage
            .windows(2)
            .map(|seg| {
                let a = &centroids[seg[0]];
                let b = &centroids[seg[1]];
                let mut ab_norm2 = 0.0f32;
                let mut ap_dot = 0.0f32;
                for d in 0..n_dims {
                    let ab = b[d] - a[d];
                    let ap = xi[d] - a[d];
                    ab_norm2 += ab * ab;
                    ap_dot += ap * ab;
                }
                if ab_norm2 <= 1e-12 {
                    a.iter().zip(xi).map(|(u, v)| (u - v).powi(2)).sum::<f32>()
                } else {
                    let t = (ap_dot / ab_norm2).clamp(0.0, 1.0);
                    let mut dist2 = 0.0f32;
                    for d in 0..n_dims {
                        let proj = a[d] + t * (b[d] - a[d]);
                        let diff = xi[d] - proj;
                        dist2 += diff * diff;
                    }
                    dist2
                }
            })
            .fold(f32::INFINITY, f32::min),
    }
}

/// Assign soft lineage weights to cells based on distance to the lineage polyline.
/// For each cell, weight[lineage] ∝ exp(-dist_to_lineage_path² / (2 sigma²)),
/// except that membership in a cluster on the lineage path receives raw weight 1.
/// Returns flat [n_cells × n_lineages], row-major.
pub fn compute_lineage_weights(
    embedding: &[f32],
    cluster_labels: &[usize],
    n_cells: usize,
    n_dims: usize,
    centroids: &[Vec<f32>],
    lineages: &[Vec<usize>],
    sigma: f32,
) -> Vec<f32> {
    let n_lineages = lineages.len();
    // Pre-build per-lineage cluster HashSets once: O(total path length) instead of
    // O(n_cells × n_lineages × path_length) with Vec::contains() inside the hot loop.
    let lineage_sets: Vec<HashSet<usize>> = lineages
        .iter()
        .map(|lin| lin.iter().copied().collect())
        .collect();
    let inv_2s2 = 1.0 / (2.0 * sigma * sigma);

    // Each cell writes to its own row [i*n_lineages..(i+1)*n_lineages]: no write conflict.
    // flat_map_iter: outer loop parallel over cells; inner row.into_iter() sequential
    // (n_lineages is typically 3-10, nested parallelism would add more overhead than gain).
    (0..n_cells)
        .into_par_iter()
        .flat_map_iter(|i| {
            let xi = &embedding[i * n_dims..(i + 1) * n_dims];
            let ci = cluster_labels[i];
            let mut row = vec![0.0f32; n_lineages];
            let mut row_sum = 0.0f32;
            for (lin, lineage) in lineages.iter().enumerate() {
                let min_d2 = point_to_lineage_sq(xi, lineage, centroids, n_dims);
                let w = if lineage_sets[lin].contains(&ci) {
                    1.0f32
                } else {
                    (-min_d2 * inv_2s2).exp()
                };
                row[lin] = w;
                row_sum += w;
            }
            if row_sum > 0.0 {
                for w in &mut row {
                    *w /= row_sum;
                }
            }
            row.into_iter()
        })
        .collect()
}

// ─── Fate probability (absorbing Markov chain, power iteration) ──────────────

/// Build a directed Markov chain with a sigmoid soft gate. A hard
/// forward-only filter would make each local pseudotime maximum absorbing;
/// the soft gate preserves transitions while weighting them by the
/// pseudotime direction.
///
/// ## Soft-gate formula
///
/// For each neighbor j of cell i:
///   gate(i→j) = σ((t_j − t_i) / bw) + ε
///   w(i→j)   = edge_weight(d_ij) × gate(i→j)
///   T(i→j)   = w(i→j) / Σ_k w(i→k)     (row-stochastic)
///
/// where:
///   σ(x) = 1/(1+e^{-x})  (sigmoid, maps Δt to (0,1))
///   bw    = bandwidth     (controls sharpness of directional preference)
///   ε     = DIFFUSION_LEAK (small constant guaranteeing ergodicity)
///
/// Properties:
///   - Forward edges (t_j > t_i): gate → 1+ε, strongly preferred
///   - Backward edges (t_j < t_i): gate → ε, allowed but heavily suppressed
///   - No cell can ever become a 100% absorbing state unless it has ZERO neighbours
///     (truly isolated node, not a noise artifact)
///   - True biological terminals are defined externally via `terminal_cells` in
///     `fate_probabilities`, not by local pseudotime maxima
///
/// NaN pseudotime cells (disconnected from root) receive a self-loop and do NOT
/// participate in fate probability diffusion — this is the correct behaviour since
/// Palantir also excludes disconnected components explicitly.
pub fn directed_markov(
    neighbors: &[Vec<(usize, f32)>],
    pseudotime: &[f32],
    bandwidth: f32,
) -> Vec<Vec<(usize, f32)>> {
    const DIFFUSION_LEAK: f32 = 1e-3;
    let bw = bandwidth.max(1e-6);
    let n = neighbors.len();
    (0..n)
        .into_par_iter()
        .map(|i| {
            // NaN pseudotime = disconnected from root → self-loop, excluded from diffusion
            if pseudotime[i].is_nan() {
                return vec![(i, 1.0f32)];
            }
            let ti = pseudotime[i];

            // Soft-gate: ALL neighbours get a transition weight, forward ones dominate.
            //
            // Edge weight = 1/(d+ε): distance-inverse, independent of bandwidth.
            // Using a Gaussian kernel exp(-d²/2bw²) here is WRONG for PCA-space distances
            // (d ~ 0.3-5.0) combined with pseudotime bandwidth (bw ~ 0.05): the exponent
            // reaches -200, causing f32 underflow to exactly 0.  When all w_edge=0 the
            // row sum s=0 and fate_probabilities converges to all-zeros.
            // The kNN graph structure already encodes proximity; no secondary Gaussian
            // is needed.  distance-inverse keeps near neighbours more influential without
            // scale-dependent underflow.
            let mut edges: Vec<(usize, f32)> = neighbors[i]
                .iter()
                .filter_map(|&(j, d)| {
                    if j >= n || pseudotime[j].is_nan() {
                        return None;
                    }
                    // Distance-inverse edge weight: near neighbours dominate, no underflow
                    let w_edge = 1.0f32 / (d + 1e-6);
                    // Sigmoid directional gate: > 0.5 for forward, < 0.5 for backward
                    let dt = (pseudotime[j] - ti) / bw;
                    let sigma = 1.0 / (1.0 + (-dt).exp());
                    // Leak term ε ensures no zero-weight rows (prevents fake absorbing states)
                    Some((j, w_edge * (sigma + DIFFUSION_LEAK)))
                })
                .collect();

            // Truly isolated cell (no valid neighbours): self-loop is the only option
            if edges.is_empty() {
                return vec![(i, 1.0f32)];
            }

            let s: f32 = edges.iter().map(|(_, w)| w).sum();
            // s > 0 guaranteed by DIFFUSION_LEAK > 0, but guard against floating-point underflow
            if s > 0.0 {
                edges.iter_mut().for_each(|(_, w)| *w /= s);
            }
            edges
        })
        .collect()
}

/// Power iteration to compute fate probabilities for absorbing states.
/// `terminal_mask[i]` = true → cell is a terminal/absorbing state.
/// Returns fate_prob[i * n_terminals + k] = probability cell i reaches terminal k.
pub fn fate_probabilities(
    n_cells: usize,
    trans: &[Vec<(usize, f32)>],
    terminal_cells: &[usize],
    max_iter: usize,
    tol: f32,
) -> Vec<f32> {
    let n_t = terminal_cells.len();
    if n_t == 0 {
        return vec![0.0f32; n_cells];
    }

    let mut terminal_index = vec![usize::MAX; n_cells];
    for (k, &tc) in terminal_cells.iter().enumerate() {
        terminal_index[tc] = k;
    }

    let mut prob = vec![1.0f32 / n_t as f32; n_cells * n_t];
    for (k, &tc) in terminal_cells.iter().enumerate() {
        for kk in 0..n_t {
            prob[tc * n_t + kk] = if kk == k { 1.0 } else { 0.0 };
        }
    }

    let mut new_prob = vec![0.0f32; n_cells * n_t];
    for _iter in 0..max_iter {
        let delta = new_prob
            .par_chunks_mut(n_t)
            .enumerate()
            .map(|(i, row)| {
                let terminal_k = terminal_index[i];
                if terminal_k != usize::MAX {
                    row.fill(0.0);
                    row[terminal_k] = 1.0;
                    0.0
                } else {
                    let mut max_change = 0.0f32;
                    for k in 0..n_t {
                        let nv: f32 = trans[i].iter().map(|&(j, w)| w * prob[j * n_t + k]).sum();
                        let ch = (nv - prob[i * n_t + k]).abs();
                        row[k] = nv;
                        if ch > max_change {
                            max_change = ch;
                        }
                    }
                    max_change
                }
            })
            .reduce(|| 0.0f32, f32::max);

        std::mem::swap(&mut prob, &mut new_prob);
        if delta < tol {
            break;
        }
    }
    prob
}

/// Cells without a directed path to any supplied absorbing cell.
pub fn attractor_unreachable_cells(
    nn_idx: &[usize],
    nn_dist: &[f64],
    n_cells: usize,
    k: usize,
    terminal_cells: &[usize],
) -> Vec<usize> {
    let mut incoming = vec![Vec::new(); n_cells];
    for cell in 0..n_cells {
        for rank in 0..k {
            let offset = cell * k + rank;
            let neighbor = nn_idx[offset];
            let distance = nn_dist[offset];
            if neighbor < n_cells && distance.is_finite() && 1.0 / (distance + 1e-6) > 0.0 {
                incoming[neighbor].push(cell);
            }
        }
    }
    let mut reachable = vec![false; n_cells];
    let mut pending = terminal_cells.to_vec();
    for &cell in terminal_cells {
        reachable[cell] = true;
    }
    while let Some(cell) = pending.pop() {
        for &predecessor in &incoming[cell] {
            if !reachable[predecessor] {
                reachable[predecessor] = true;
                pending.push(predecessor);
            }
        }
    }
    reachable.iter().enumerate()
        .filter_map(|(cell, &value)| (!value).then_some(cell))
        .collect()
}

/// Double-precision absorption probabilities on the supplied directed kNN
/// transition. Solve the harmonic equations rather than stopping a diffusion
/// when its per-step change is small: rare absorption can make that change tiny
/// while the probabilities are still close to their uniform initialization.
pub fn attractor_fate_probabilities_f64(
    nn_idx: &[usize],
    nn_dist: &[f64],
    n_cells: usize,
    k: usize,
    terminal_cells: &[usize],
    max_iter: usize,
    tol: f64,
) -> Result<Vec<f64>, String> {
    let n_terminals = terminal_cells.len();
    if n_terminals == 0 {
        return Ok(Vec::new());
    }

    let transition: Vec<Vec<(usize, f64)>> = (0..n_cells)
        .into_par_iter()
        .map(|cell| {
            let mut edges: Vec<(usize, f64)> = (0..k)
                .filter_map(|rank| {
                    let offset = cell * k + rank;
                    let neighbor = nn_idx[offset];
                    let distance = nn_dist[offset];
                    (neighbor < n_cells && distance.is_finite())
                        .then_some((neighbor, 1.0 / (distance + 1e-6)))
                })
                .collect();
            if edges.is_empty() {
                return vec![(cell, 1.0)];
            }
            let total: f64 = edges.iter().map(|(_, weight)| weight).sum();
            edges.iter_mut().for_each(|(_, weight)| *weight /= total);
            edges
        })
        .collect();

    let mut terminal_index = vec![usize::MAX; n_cells];
    for (terminal, &cell) in terminal_cells.iter().enumerate() {
        terminal_index[cell] = terminal;
    }

    // A closed communicating class without any absorbing cell has no defined
    // terminal distribution. Its harmonic residual can be zero for any prior.
    let unreachable = attractor_unreachable_cells(nn_idx, nn_dist, n_cells, k, terminal_cells);
    if !unreachable.is_empty() {
        return Err(format!(
            "Absorption is undefined: {} cells cannot reach any terminal", unreachable.len()
        ));
    }

    if n_terminals == 1 {
        return Ok(vec![1.0; n_cells]);
    }

    // The harmonic operator is an M-matrix. A zero-fill incomplete LU keeps
    // the kNN sparsity and removes the near-closed-class scale from the Krylov
    // iteration without changing the graph or adding artificial absorption.
    let mut factors: Vec<Vec<(usize, f64)>> = transition
        .iter()
        .enumerate()
        .map(|(cell, edges)| {
            let mut row = std::collections::BTreeMap::new();
            row.insert(cell, f64::from(terminal_index[cell] != usize::MAX));
            if terminal_index[cell] == usize::MAX {
                for &(neighbor, weight) in edges {
                    if neighbor != cell {
                        *row.entry(cell).or_insert(0.0) += weight;
                        *row.entry(neighbor).or_insert(0.0) -= weight;
                    }
                }
            }
            row.into_iter().collect()
        })
        .collect();
    let diagonal: Vec<usize> = factors.iter().enumerate()
        .map(|(cell, row)| row.binary_search_by_key(&cell, |entry| entry.0).unwrap())
        .collect();
    for cell in 0..n_cells {
        let mut row = std::mem::take(&mut factors[cell]);
        for entry in 0..diagonal[cell] {
            let previous = row[entry].0;
            let multiplier = row[entry].1 / factors[previous][diagonal[previous]].1;
            row[entry].1 = multiplier;
            for &(column, value) in &factors[previous][diagonal[previous] + 1..] {
                if let Ok(slot) = row.binary_search_by_key(&column, |item| item.0) {
                    row[slot].1 -= multiplier * value;
                }
            }
        }
        if !row[diagonal[cell]].1.is_finite() || row[diagonal[cell]].1 <= 0.0 {
            return Err(format!("Absorption preconditioner has a nonpositive pivot at cell {}", cell + 1));
        }
        factors[cell] = row;
    }
    let precondition = |rhs: &[f64]| -> Vec<f64> {
        let mut answer = rhs.to_vec();
        for cell in 0..n_cells {
            for &(column, value) in &factors[cell][..diagonal[cell]] {
                answer[cell] -= value * answer[column];
            }
        }
        for cell in (0..n_cells).rev() {
            for &(column, value) in &factors[cell][diagonal[cell] + 1..] {
                answer[cell] -= value * answer[column];
            }
            answer[cell] /= factors[cell][diagonal[cell]].1;
        }
        answer
    };

    // Right-preconditioned GMRES minimizes the actual harmonic residual.
    // Restarting bounds memory; reorthogonalization and residual recomputation
    // avoid the unstable recursive residual seen with unpreconditioned BiCGSTAB.
    // Absorbing rows are identities and their search directions stay zero.
    // A tighter backward tolerance is needed because rare terminal visits make
    // the harmonic system ill-conditioned. The original tol also bounds any
    // probability-range violation after the independent columns are assembled.
    let columns: Result<Vec<Vec<f64>>, String> = (0..n_terminals.saturating_sub(1))
        .into_par_iter()
        .map(|terminal| {
            let apply = |values: &[f64]| -> Vec<f64> {
                transition
                    .par_iter()
                    .enumerate()
                    .map(|(cell, edges)| {
                        if terminal_index[cell] != usize::MAX {
                            values[cell]
                        } else {
                            edges
                                .iter()
                                .map(|&(neighbor, weight)| {
                                    weight * (values[cell] - values[neighbor])
                                })
                                .sum()
                        }
                    })
                    .collect()
            };
            let dot = |left: &[f64], right: &[f64]| -> f64 {
                left.par_iter().zip(right.par_iter()).map(|(a, b)| a * b).sum()
            };
            let mut solution = vec![1.0 / n_terminals as f64; n_cells];
            for (index, &cell) in terminal_cells.iter().enumerate() {
                solution[cell] = f64::from(index == terminal);
            }
            let mut residual = apply(&solution);
            residual.iter_mut().for_each(|value| *value = -*value);
            residual[terminal_cells[terminal]] += 1.0;
            let initial_norm = dot(&residual, &residual).sqrt();
            let threshold = (tol * f64::EPSILON.sqrt() * initial_norm)
                .max(8.0 * f64::EPSILON * (n_cells as f64).sqrt());
            let mut converged = initial_norm <= threshold;
            let mut iterations = 0;
            const RESTART: usize = 100;
            while !converged && iterations < max_iter {
                let cycle_length = RESTART.min(n_cells).min(max_iter - iterations);
                let residual_norm = dot(&residual, &residual).sqrt();
                let mut basis = vec![residual.iter().map(|value| value / residual_norm).collect::<Vec<_>>()];
                let mut hessenberg = vec![vec![0.0; cycle_length]; cycle_length + 1];
                let mut cosine = vec![0.0; cycle_length];
                let mut sine = vec![0.0; cycle_length];
                let mut projected_rhs = vec![0.0; cycle_length + 1];
                projected_rhs[0] = residual_norm;
                let mut used = 0;
                for column in 0..cycle_length {
                    let mut work = apply(&precondition(&basis[column]));
                    // Two modified Gram-Schmidt passes preserve the basis for
                    // almost-dependent directions in rare-absorption systems.
                    for _ in 0..2 {
                        for row in 0..=column {
                            let projection = dot(&basis[row], &work);
                            hessenberg[row][column] += projection;
                            work.par_iter_mut().zip(basis[row].par_iter())
                                .for_each(|(value, coordinate)| *value -= projection * *coordinate);
                        }
                    }
                    let next_norm = dot(&work, &work).sqrt();
                    hessenberg[column + 1][column] = next_norm;
                    if next_norm > 0.0 {
                        basis.push(work.iter().map(|value| value / next_norm).collect());
                    }
                    for row in 0..column {
                        let upper = hessenberg[row][column];
                        let lower = hessenberg[row + 1][column];
                        hessenberg[row][column] = cosine[row] * upper + sine[row] * lower;
                        hessenberg[row + 1][column] = -sine[row] * upper + cosine[row] * lower;
                    }
                    let rotation_norm = hessenberg[column][column].hypot(hessenberg[column + 1][column]);
                    if !rotation_norm.is_finite() || rotation_norm == 0.0 {
                        return Err(format!("Absorption solve has a singular search direction for terminal {}", terminal + 1));
                    }
                    cosine[column] = hessenberg[column][column] / rotation_norm;
                    sine[column] = hessenberg[column + 1][column] / rotation_norm;
                    hessenberg[column][column] = rotation_norm;
                    hessenberg[column + 1][column] = 0.0;
                    projected_rhs[column + 1] = -sine[column] * projected_rhs[column];
                    projected_rhs[column] *= cosine[column];
                    used = column + 1;
                    iterations += 1;
                    if projected_rhs[column + 1].abs() <= threshold || next_norm == 0.0 {
                        break;
                    }
                }
                let mut coefficients = projected_rhs[..used].to_vec();
                for row in (0..used).rev() {
                    for column in row + 1..used {
                        coefficients[row] -= hessenberg[row][column] * coefficients[column];
                    }
                    coefficients[row] /= hessenberg[row][row];
                }
                let mut correction = vec![0.0; n_cells];
                for column in 0..used {
                    for cell in 0..n_cells {
                        correction[cell] += coefficients[column] * basis[column][cell];
                    }
                }
                let correction = precondition(&correction);
                for cell in 0..n_cells {
                    solution[cell] += correction[cell];
                }
                residual = apply(&solution);
                residual.iter_mut().for_each(|value| *value = -*value);
                residual[terminal_cells[terminal]] += 1.0;
                let actual_norm = dot(&residual, &residual).sqrt();
                if !actual_norm.is_finite() {
                    return Err(format!("Absorption solve has a nonfinite residual for terminal {}", terminal + 1));
                }
                converged = actual_norm <= threshold;
            }
            let mut actual_residual = apply(&solution);
            actual_residual[terminal_cells[terminal]] -= 1.0;
            let residual_norm = dot(&actual_residual, &actual_residual).sqrt();
            // The iteration target is stricter than the caller's relative tolerance.
            // Accept a finite solution that meets the requested tolerance at max_iter.
            if !residual_norm.is_finite()
                || residual_norm > threshold.max(tol * initial_norm)
            {
                return Err(format!(
                    "Absorption solve did not converge for terminal {} in {} iterations (residual {:e})",
                    terminal + 1, max_iter, residual_norm
                ));
            }
            Ok(solution)
        })
        .collect();
    let columns = columns?;
    let mut probability = vec![0.0; n_cells * n_terminals];
    for (cell, row) in probability.chunks_mut(n_terminals).enumerate() {
        for (terminal, column) in columns.iter().enumerate() {
            row[terminal] = column[cell];
        }
        row[n_terminals - 1] = 1.0 - row[..n_terminals - 1].iter().sum::<f64>();
    }
    if probability.iter().any(|value| !value.is_finite() || *value < -tol || *value > 1.0 + tol) {
        return Err("Absorption solve produced probabilities outside [0, 1]".into());
    }
    Ok(probability)
}

// ─── Trajectory differential expression ──────────────────────────────────────

/// Rank a slice of f32 values (1-based, average-rank ties).
/// NaN values receive rank 0.0 and are excluded from tie resolution.
fn rank_vec(v: &[f32]) -> Vec<f64> {
    crate::ccc::rank::average_rank_f32(v)
}

/// Pearson correlation of two f64 slices (pre-ranked for Spearman).
fn pearson_f64(x: &[f64], y: &[f64]) -> f64 {
    let n = x.len() as f64;
    if n < 3.0 {
        return 0.0;
    }
    let mx = x.iter().sum::<f64>() / n;
    let my = y.iter().sum::<f64>() / n;
    let num: f64 = x.iter().zip(y).map(|(a, b)| (a - mx) * (b - my)).sum();
    let dx = x.iter().map(|a| (a - mx).powi(2)).sum::<f64>().sqrt();
    let dy = y.iter().map(|b| (b - my).powi(2)).sum::<f64>().sqrt();
    if dx > 1e-12 && dy > 1e-12 {
        num / (dx * dy)
    } else {
        0.0
    }
}

/// Complementary error function approximation (Abramowitz & Stegun 7.1.26).
/// Max absolute error: 1.5e-7.
fn erfc_approx(x: f64) -> f64 {
    let t = 1.0 / (1.0 + 0.3275911 * x);
    let poly = t
        * (0.254829592
            + t * (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429))));
    poly * (-x * x).exp()
}

/// Two-tailed p-value from t-statistic using normal approximation.
/// Valid for df >= 10; adequate for the n >= min_cells constraint used below.
fn t_pval(t: f64, df: f64) -> f64 {
    if df < 1.0 {
        return 1.0;
    }
    // Apply a df-correction to the t-statistic before the normal approximation.
    // The Cornish-Fisher correction z ≈ t*(1 - 1/(4*df)) shrinks small-df tails.
    let z = t.abs() * (1.0 - 1.0 / (4.0 * df)).max(0.0).sqrt();
    erfc_approx(z / std::f64::consts::SQRT_2).min(1.0)
}

/// Benjamini-Hochberg FDR adjustment.
///
/// Input: p-values in original gene order.
/// Output: q-values (adjusted p-values) in the same order, clamped to [0, 1].
fn bh_adjust(pvals: &[f64]) -> Vec<f64> {
    let n = pvals.len();
    if n == 0 {
        return vec![];
    }
    let mut indexed: Vec<(usize, f64)> = pvals.iter().enumerate().map(|(i, &p)| (i, p)).collect();
    indexed.sort_unstable_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));

    let mut qvals = vec![1.0f64; n];
    let mut cummin = 1.0f64;
    for rank0 in (0..n).rev() {
        let (orig, p) = indexed[rank0];
        let q = (p * n as f64 / (rank0 + 1) as f64).min(1.0);
        cummin = cummin.min(q);
        qvals[orig] = cummin;
    }
    qvals
}

/// Compute per-gene Spearman correlation with pseudotime (trajectory DE).
///
/// For each gene the function:
///   1. Identifies cells with BOTH finite pseudotime AND expression > 0
///      (cells with NaN pseudotime are disconnected from the root and excluded).
///   2. Skips the gene if fewer than `min_cells` cells qualify.
///   3. Computes Spearman ρ between expression and pseudotime for qualifying cells.
///   4. Converts ρ to a t-statistic: t = ρ√((n-2)/(1-ρ²))
///   5. Computes two-tailed p-value via normal approximation (adequate for n ≥ 10).
/// After all genes, applies Benjamini-Hochberg FDR correction.
///
/// ## Returns
/// `(rho, pval, qval)` — three `Vec<f64>` of length `n_genes`.
/// Genes with fewer than `min_cells` qualifying cells receive rho=0, pval=1, qval=1.
///
/// ## Input layout
/// `expr_flat`: row-major [n_genes × n_cells] — raw or normalised counts.
///   Row i = gene i, column j = cell j.
/// `pseudotime`: length n_cells; NaN for cells disconnected from root.
/// `min_cells`:  minimum qualifying cells (default 10).
pub fn pseudotime_de(
    expr_flat: &[f32],
    pseudotime: &[f32],
    n_genes: usize,
    n_cells: usize,
    min_cells: usize,
) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    // Pre-compute which cells have valid (finite) pseudotime — avoids per-gene filtering
    let valid_cells: Vec<usize> = (0..n_cells)
        .filter(|&c| !pseudotime[c].is_nan() && pseudotime[c].is_finite())
        .collect();

    let results: Vec<(f64, f64)> = (0..n_genes)
        .into_par_iter()
        .map(|g| {
            let row = &expr_flat[g * n_cells..(g + 1) * n_cells];

            // Filter to cells with valid pseudotime AND expression > 0
            let sel: Vec<usize> = valid_cells
                .iter()
                .filter(|&&c| row[c] > 0.0)
                .copied()
                .collect();
            let n = sel.len();
            if n < min_cells {
                return (0.0f64, 1.0f64);
            }

            let expr_sel: Vec<f32> = sel.iter().map(|&c| row[c]).collect();
            let pt_sel: Vec<f32> = sel.iter().map(|&c| pseudotime[c]).collect();

            let r_expr = rank_vec(&expr_sel);
            let r_pt = rank_vec(&pt_sel);

            let rho = pearson_f64(&r_expr, &r_pt);
            let df = (n as f64) - 2.0;
            let denom = (1.0 - rho * rho).max(1e-12);
            let t = rho * (df / denom).sqrt();
            let p = t_pval(t, df);
            (rho, p)
        })
        .collect();

    let rho: Vec<f64> = results.iter().map(|(r, _)| *r).collect();
    let pval: Vec<f64> = results.iter().map(|(_, p)| *p).collect();
    let qval = bh_adjust(&pval);
    (rho, pval, qval)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A local pseudotime maximum must retain transitions to real neighbours
    /// rather than becoming a single-entry self-loop.
    #[test]
    fn test_directed_markov_no_fake_absorbing_states() {
        // Cell 2 has the highest pseudotime among its neighbours → was a fake sink.
        // Topology: 0↔1↔2↔3, cell 2 is a local maximum (t[2] > t[1] and t[2] > t[3]).
        let neighbors = vec![
            vec![(1, 0.5f32)],
            vec![(0, 0.5), (2, 0.5)],
            vec![(1, 0.5), (3, 0.5)], // local max: both neighbours have lower t
            vec![(2, 0.5)],
        ];
        let pseudotime = vec![0.0f32, 0.3, 0.6, 0.4]; // cell 2 is local max

        let trans = directed_markov(&neighbors, &pseudotime, 0.2);

        // Every non-isolated cell must have ≥ 2 entries (no pure self-loop)
        // Cell 2: neighbours are 1 (t=0.3 < 0.6) and 3 (t=0.4 < 0.6)
        // The soft gate diffuses from this local maximum to both neighbours.
        assert!(
            trans[2].len() >= 2,
            "cell 2 (local pseudotime max) must have transitions to both neighbours, got {:?}",
            trans[2]
        );

        // All rows must be row-stochastic (sum ≈ 1.0)
        for (i, row) in trans.iter().enumerate() {
            let s: f32 = row.iter().map(|(_, w)| w).sum();
            assert!(
                (s - 1.0).abs() < 1e-5,
                "row {} is not stochastic: sum = {}",
                i,
                s
            );
        }

        // Forward direction must dominate: cell 0→1 (forward) > cell 1→0 (backward)
        // Cell 1 transitions to both 0 (backward) and 2 (forward).
        let w_fwd = trans[1]
            .iter()
            .find(|&&(j, _)| j == 2)
            .map(|(_, w)| *w)
            .unwrap_or(0.0);
        let w_bwd = trans[1]
            .iter()
            .find(|&&(j, _)| j == 0)
            .map(|(_, w)| *w)
            .unwrap_or(0.0);
        assert!(
            w_fwd > w_bwd,
            "forward edge (1→2) should outweigh backward edge (1→0): {:.4} vs {:.4}",
            w_fwd,
            w_bwd
        );
    }

    /// Cells with non-finite pseudotime are isolated with self-loops and do not
    /// affect neighbouring transition weights.
    #[test]
    fn test_directed_markov_nan_pseudotime_isolated() {
        let neighbors = vec![
            vec![(1, 0.5f32), (2, 0.5)],
            vec![(0, 0.5)],
            vec![(0, 0.5)], // cell 2 has NaN pseudotime
        ];
        let pseudotime = vec![0.0f32, 0.5, f32::NAN];

        let trans = directed_markov(&neighbors, &pseudotime, 0.2);

        // Cell 2 (NaN) → self-loop
        assert_eq!(trans[2].len(), 1);
        assert_eq!(trans[2][0].0, 2);
        assert!((trans[2][0].1 - 1.0).abs() < 1e-6);

        // Cell 0: neighbour 2 has NaN → filtered out, only 1 remains
        let nb_indices: Vec<usize> = trans[0].iter().map(|(j, _)| *j).collect();
        assert!(
            !nb_indices.contains(&2),
            "NaN cell should not appear in transitions"
        );
    }

    #[test]
    fn test_knn_graph_from_hnsw_is_symmetrized_and_reachable() {
        // Directed 1-NN input:
        // 0 -> 1, 1 -> 0, 2 -> 1
        // Without symmetrization, 0 cannot reach 2.
        let nn_idx = vec![1usize, 0usize, 1usize];
        let nn_dist = vec![1.0f32, 1.0f32, 1.0f32];
        let g = knn_graph_from_hnsw(&nn_idx, &nn_dist, 3, 1);

        // Graph must contain both directions for every edge.
        for i in 0..g.len() {
            for &(j, _) in &g[i] {
                assert!(
                    g[j].iter().any(|&(jj, _)| jj == i),
                    "edge {} -> {} missing reverse edge",
                    i,
                    j
                );
            }
        }

        // Dijkstra from 0 should now reach all nodes.
        let dist = dijkstra(3, &g, 0);
        assert!(dist.iter().all(|d| d.is_finite()));
    }

    #[test]
    fn fuzzy_pseudotime_sorts_distances_with_their_indices() {
        let unsorted_idx = vec![2usize, 1, 2, 0, 0, 1];
        let unsorted_dist = vec![2.0f32, 1.0, 1.0, 1.0, 2.0, 1.0];
        let sorted_idx = vec![1usize, 2, 0, 2, 1, 0];
        let sorted_dist = vec![1.0f32, 2.0, 1.0, 1.0, 1.0, 2.0];
        let left = compute_pseudotime_from_hnsw_fuzzy(&unsorted_idx, &unsorted_dist, 3, 2, 0);
        let right = compute_pseudotime_from_hnsw_fuzzy(&sorted_idx, &sorted_dist, 3, 2, 0);
        assert_eq!(left, right);
    }

    #[test]
    fn fuzzy_pseudotime_ignores_invalid_padding_before_local_scale() {
        let compact_idx = vec![1usize, 0, 1];
        let compact_dist = vec![1.0f32, 1.0, 1.0];
        let padded_idx = vec![1usize, usize::MAX, 0, usize::MAX, 1, usize::MAX];
        let padded_dist = vec![1.0f32, 0.0, 1.0, 0.0, 1.0, 0.0];
        let compact = compute_pseudotime_from_hnsw_fuzzy(&compact_idx, &compact_dist, 3, 1, 0);
        let padded = compute_pseudotime_from_hnsw_fuzzy(&padded_idx, &padded_dist, 3, 2, 0);
        assert_eq!(compact, padded);
    }

    #[test]
    fn test_lineage_entropy_expected_bounds() {
        // cell0: deterministic lineage [1,0] => entropy 0
        // cell1: uniform lineage [0.5,0.5] => entropy 1
        let lw = vec![1.0f32, 0.0, 0.5, 0.5];
        let h = lineage_entropy(&lw, 2, 2);
        assert!((h[0] - 0.0).abs() < 1e-6);
        assert!((h[1] - 1.0).abs() < 1e-6);
    }

    #[test]
    fn attractor_fate_f64_preserves_symmetric_absorption() {
        let nn_idx = vec![1usize, 2, 0, 2, 0, 1];
        let nn_dist = vec![1.0f64; 6];
        let probability =
            attractor_fate_probabilities_f64(&nn_idx, &nn_dist, 3, 2, &[0, 2], 300, 1e-12).unwrap();
        assert_eq!(&probability[0..2], &[1.0, 0.0]);
        assert!((probability[2] - 0.5).abs() < 1e-12);
        assert!((probability[3] - 0.5).abs() < 1e-12);
        assert_eq!(&probability[4..6], &[0.0, 1.0]);
    }

    #[test]
    fn attractor_fate_resolves_rare_absorption_without_uniform_prior_bias() {
        let nn_idx = vec![0usize, 1, 2, 0, 1, 2, 0, 1, 2];
        let nn_dist = vec![1.0f64, 1.0, 1.0, 100.0, 0.0, 10.0, 1.0, 1.0, 1.0];
        let probability =
            attractor_fate_probabilities_f64(&nn_idx, &nn_dist, 3, 3, &[0, 2], 300, 1e-5).unwrap();
        let first_weight = 1.0 / (100.0 + 1e-6);
        let last_weight = 1.0 / (10.0 + 1e-6);
        let expected = first_weight / (first_weight + last_weight);
        assert!((probability[2] - expected).abs() < 1e-8);
        assert!((probability[3] - (1.0 - expected)).abs() < 1e-8);
    }

    #[test]
    fn attractor_fate_resolves_a_nearly_closed_directed_cycle() {
        let nn_idx = vec![0usize, 0, 0, 2, 0, 5, 3, 0, 5, 4, 0, 5, 1, 0, 5, 5, 5, 5];
        let nn_dist = vec![
            1.0, 1.0, 1.0, 1.0, 1e9, 3e9, 1.0, 1e9, 3e9,
            1.0, 1e9, 3e9, 1.0, 1e9, 3e9, 1.0, 1.0, 1.0,
        ];
        let probability = attractor_fate_probabilities_f64(
            &nn_idx, &nn_dist, 6, 3, &[0, 5], 300, 1e-5,
        ).unwrap();
        for cell in 1..5 {
            assert!((probability[cell * 2] - 0.75).abs() < 1e-7);
            assert!((probability[cell * 2 + 1] - 0.25).abs() < 1e-7);
        }
    }

    #[test]
    fn attractor_fate_reports_exhausted_solver_instead_of_returning_prior() {
        let nn_idx = vec![1usize, 2, 0, 2, 0, 1];
        let nn_dist = vec![1.0f64, 1.0, 1.0, 2.0, 1.0, 1.0];
        let result = attractor_fate_probabilities_f64(&nn_idx, &nn_dist, 3, 2, &[0, 2], 0, 1e-5);
        assert!(result.unwrap_err().contains("did not converge"));
    }

    #[test]
    fn attractor_fate_rejects_closed_classes_without_absorbing_cells() {
        let result = attractor_fate_probabilities_f64(
            &[0usize, 1, 2], &[1.0f64; 3], 3, 1, &[0, 2], 300, 1e-5,
        );
        assert!(result.unwrap_err().contains("1 cells cannot reach any terminal"));
    }

    #[test]
    fn test_dtw_similarity_identical_higher_than_shifted() {
        let a = vec![0.1f32, 0.2, 0.4, 0.8, 1.0];
        let b = vec![0.1f32, 0.2, 0.4, 0.8, 1.0];
        let c = vec![1.0f32, 0.8, 0.4, 0.2, 0.1];
        let sim_same = dtw_similarity(&a, &b, 2);
        let sim_diff = dtw_similarity(&a, &c, 2);
        assert!(
            sim_same > sim_diff,
            "expected identical profiles to have higher DTW similarity"
        );
    }
}

// Trunk-aware geometric lineage weighting
//
// `compute_lineage_weights` assigns a hard weight to clusters on a lineage
// path. `compute_lineage_weights_trunk_aware` instead uses geometric weights
// for clusters shared by multiple lineage paths.

/// Count how many lineages each cluster id appears on.
///
/// Used by trunk-aware weighting to classify each cluster as "trunk"
/// (`count >= 2`, i.e. the cluster sits on `>=2` lineage root→leaf paths)
/// vs "single-lineage".  O(sum_l |lineage_l|) = O(N_clusters · L) in the
/// worst case.  Any cluster id outside `[0, n_clusters)` is silently
/// ignored at the container boundary.
pub fn cluster_lineage_count(lineages: &[Vec<usize>], n_clusters: usize) -> Vec<u32> {
    let mut counts = vec![0u32; n_clusters];
    for path in lineages {
        // A lineage path is a simple root→leaf path through the cluster
        // MST; it cannot visit the same cluster twice, so no dedup set
        // is needed here.
        for &c in path {
            if c < n_clusters {
                counts[c] += 1;
            }
        }
    }
    counts
}

/// Trunk-aware geometric lineage weighting.
///
/// Three cases are evaluated:
/// - **trunk cluster** (`cluster_lineage_count[ci] >= 2`): pure Gaussian
///   `exp(-d² / (2σ²))`, without a hard membership weight.
/// - **single-lineage cluster, on `lineages[lin]`'s path**: hard `1.0`
///   membership weight.
/// - **otherwise** (single-lineage cluster off-path, or any cluster
///   not appearing on `lineages[lin]`): pure Gaussian `exp(-d² / (2σ²))`.
///
/// `is_trunk(ci)` is a function of `ci` ONLY, not of `lin`, so the
/// branch decision is hoisted out of the per-lineage inner loop for
/// `O(n_cells)` instead of `O(n_cells × n_lineages)` cluster-membership
/// checks.
///
/// The signature mirrors `compute_lineage_weights`. The return value is a flat
/// row-major `[n_cells × n_lineages]` matrix; each nonzero row is normalized to
/// sum to one.
pub fn compute_lineage_weights_trunk_aware(
    embedding: &[f32],
    cluster_labels: &[usize],
    n_cells: usize,
    n_dims: usize,
    centroids: &[Vec<f32>],
    lineages: &[Vec<usize>],
    sigma: f32,
) -> Vec<f32> {
    let n_lineages = lineages.len();
    let n_clusters = centroids.len();

    // A cluster is a trunk cluster when it lies on at least two lineage paths.
    // Counting path membership adds O(N_clusters * L) preprocessing.
    let counts = cluster_lineage_count(lineages, n_clusters);

    // Per-lineage HashSets avoid O(n_cells * n_lineages * path_length)
    // Vec::contains() calls inside the hot loop.
    let lineage_sets: Vec<HashSet<usize>> = lineages
        .iter()
        .map(|lin| lin.iter().copied().collect())
        .collect();
    let inv_2s2 = 1.0 / (2.0 * sigma * sigma);

    (0..n_cells)
        .into_par_iter()
        .flat_map_iter(|i| {
            let xi = &embedding[i * n_dims..(i + 1) * n_dims];
            let ci = cluster_labels[i];
            // Hoisted-out trunk decision: depends on ci alone, not on lin.
            // Bounds checking keeps malformed cluster labels from indexing
            // outside the path-membership vector.
            let is_trunk_ci = ci < n_clusters && counts[ci] >= 2;

            let mut row = vec![0.0f32; n_lineages];
            let mut max_log_weight = f32::NEG_INFINITY;
            for (lin, lineage) in lineages.iter().enumerate() {
                let min_d2 = point_to_lineage_sq(xi, lineage, centroids, n_dims);
                let log_weight = if is_trunk_ci {
                    // Shared trunk clusters use geometry only.
                    -min_d2 * inv_2s2
                } else if lineage_sets[lin].contains(&ci) {
                    // A single-lineage cluster on its own path receives weight 1.0.
                    0.0f32
                } else {
                    // Off-path geometric fall-back.
                    -min_d2 * inv_2s2
                };
                row[lin] = log_weight;
                max_log_weight = max_log_weight.max(log_weight);
            }
            let mut row_sum = 0.0f32;
            if max_log_weight.is_finite() {
                for w in &mut row {
                    *w = (*w - max_log_weight).exp();
                    row_sum += *w;
                }
                for w in &mut row {
                    *w /= row_sum;
                }
            }
            row.into_iter()
        })
        .collect()
}

/// Assign lineage weights from local expression deviations at each tree fork.
/// In assignment mode, a cell downstream of a child keeps that branch; cells at or above
/// the fork use their displacement from their own cluster center to choose
/// between the outgoing child directions. Each fork's contrast is measured
/// in units of its within-parent expression spread, so embedding scale does
/// not force early cells into nearly deterministic branches. Priming mode
/// keeps continuous fork-centered expression differences in descendant cells.
pub fn compute_lineage_weights_local_branch_residual(
    embedding: &[f32],
    cluster_labels: &[usize],
    n_cells: usize,
    n_dims: usize,
    centroids: &[Vec<f32>],
    lineages: &[Vec<usize>],
    priming: bool,
) -> Vec<f32> {
    let n_clusters = centroids.len();
    let mut children = vec![Vec::<usize>::new(); n_clusters];
    let mut parent = vec![usize::MAX; n_clusters];
    for path in lineages {
        for edge in path.windows(2) {
            let (from, to) = (edge[0], edge[1]);
            if !children[from].contains(&to) {
                children[from].push(to);
                parent[to] = from;
            }
        }
    }

    let mut descendant = vec![vec![false; n_clusters]; n_clusters];
    for node in 0..n_clusters {
        let mut current = node;
        while current != usize::MAX {
            descendant[current][node] = true;
            current = parent[current];
        }
    }

    let mut cluster_cells = vec![Vec::<usize>::new(); n_clusters];
    for (cell, &cluster) in cluster_labels.iter().enumerate() {
        cluster_cells[cluster].push(cell);
    }
    let mut fork_scale = vec![1.0f64; n_clusters];
    for (fork, branches) in children.iter().enumerate() {
        if branches.len() < 2 {
            continue;
        }
        let mut squared_contrast = 0.0f64;
        let mut n_pairs = 0usize;
        for left in 0..branches.len() {
            for right in left + 1..branches.len() {
                n_pairs += 1;
                for &cell in &cluster_cells[fork] {
                    let expression = &embedding[cell * n_dims..(cell + 1) * n_dims];
                    let projection: f64 = (0..n_dims)
                        .map(|dim| {
                            (expression[dim] - centroids[fork][dim]) as f64
                                * (centroids[branches[left]][dim]
                                    - centroids[branches[right]][dim]) as f64
                        })
                        .sum();
                    squared_contrast += projection * projection;
                }
            }
        }
        fork_scale[fork] = (squared_contrast
            / (cluster_cells[fork].len() * n_pairs) as f64)
            .sqrt();
    }

    (0..n_cells)
        .into_par_iter()
        .flat_map_iter(|cell| {
            let cell_cluster = cluster_labels[cell];
            let expression = &embedding[cell * n_dims..(cell + 1) * n_dims];
            let mut child_weight = vec![0.0f64; n_clusters];

            for (fork, branches) in children.iter().enumerate() {
                if branches.is_empty() {
                    continue;
                }
                if !priming {
                    if let Some(&chosen) = branches
                        .iter()
                        .find(|&&branch| descendant[branch][cell_cluster])
                    {
                        child_weight[chosen] = 1.0;
                        continue;
                    }
                }
                let center = if priming { &centroids[fork] } else { &centroids[cell_cluster] };

                let logits: Vec<f64> = branches
                    .iter()
                    .map(|&branch| {
                        (0..n_dims)
                            .map(|dim| {
                                (expression[dim] - center[dim]) as f64
                                    * (centroids[branch][dim] - centroids[fork][dim]) as f64
                            })
                            .sum::<f64>()
                            / fork_scale[fork]
                    })
                    .collect();
                let max_logit = logits.iter().copied().fold(f64::NEG_INFINITY, f64::max);
                let normalizer: f64 = logits.iter().map(|v| (v - max_logit).exp()).sum();
                for (&branch, &logit) in branches.iter().zip(logits.iter()) {
                    child_weight[branch] = (logit - max_logit).exp() / normalizer;
                }
            }

            let mut row: Vec<f64> = lineages
                .iter()
                .map(|path| path.iter().skip(1).map(|&node| child_weight[node]).product())
                .collect();
            let total: f64 = row.iter().sum();
            for weight in &mut row {
                *weight /= total;
            }
            row.into_iter().map(|weight| weight as f32)
        })
        .collect()
}

#[cfg(test)]
mod trunk_aware_tests {
    //! Unit tests for trunk-aware geometric lineage weighting.
    //!
    //! Topology used by most tests below:
    //! ```text
    //!   cluster 0 (root) ── 1 (trunk) ── 2 (A end)
    //!                                ╲── 3 (B end)
    //!
    //!   lineages: A = [0, 1, 2]   B = [0, 1, 3]
    //!   centroids (2D):
    //!     0 = (0, 0)
    //!     1 = (1, 0)   ← trunk
    //!     2 = (2, 1)   ← A-only end
    //!     3 = (2,-1)   ← B-only end
    //! ```
    //! Trunk clusters: {0, 1}.  Single-lineage clusters: {2, 3}.

    use super::{
        cluster_lineage_count, compute_lineage_weights, compute_lineage_weights_trunk_aware,
        compute_lineage_weights_local_branch_residual,
    };

    /// 2-lineage bifurcation fixture used by the trunk-aware tests.
    /// Returns `(centroids_flat_2d, lineages, n_clusters, n_dims)`.
    fn bifurcation_fixture() -> (Vec<Vec<f32>>, Vec<Vec<usize>>, usize, usize) {
        let centroids: Vec<Vec<f32>> = vec![
            vec![0.0, 0.0],  // 0  root  (trunk)
            vec![1.0, 0.0],  // 1  trunk
            vec![2.0, 1.0],  // 2  A end (single-lineage A)
            vec![2.0, -1.0], // 3  B end (single-lineage B)
        ];
        let lineages = vec![vec![0, 1, 2], vec![0, 1, 3]];
        (centroids, lineages, 4, 2)
    }

    #[test]
    fn local_branch_residual_separates_trunk_priming_and_keeps_descendants() {
        let (centroids, lineages, _, _) = bifurcation_fixture();
        let embedding = vec![1.0, 0.4, 1.0, -0.4, 2.0, 1.0, 2.0, -1.0];
        let clusters = vec![1, 1, 2, 3];
        let weights = compute_lineage_weights_local_branch_residual(
            &embedding, &clusters, 4, 2, &centroids, &lineages, false,
        );
        assert!(weights[0] > weights[1]);
        assert!(weights[2] < weights[3]);
        assert_eq!(&weights[4..6], &[1.0, 0.0]);
        assert_eq!(&weights[6..8], &[0.0, 1.0]);
    }

    #[test]
    fn local_branch_residual_priming_keeps_within_branch_expression() {
        let (centroids, lineages, _, _) = bifurcation_fixture();
        let embedding = vec![1.0, 0.4, 1.0, -0.4, 2.0, 0.6, 2.0, 1.4, 2.0, -0.6, 2.0, -1.4];
        let clusters = vec![1, 1, 2, 2, 3, 3];
        let weights = compute_lineage_weights_local_branch_residual(
            &embedding, &clusters, 6, 2, &centroids, &lineages, true,
        );
        assert!((weights[4] - 0.8175745).abs() < 1e-6);
        assert!((weights[6] - 0.9706878).abs() < 1e-6);
        assert!((weights[8] - 0.1824255).abs() < 1e-6);
        assert!((weights[10] - 0.0293122).abs() < 1e-6);
        assert!(weights.chunks_exact(2).all(|row| {
            row.iter().all(|value| value.is_finite() && *value >= 0.0)
                && (row.iter().sum::<f32>() - 1.0).abs() < 1e-6
        }));
    }

    #[test]
    fn cluster_lineage_count_trunk_detection() {
        // Topology (MST branch):
        //   cluster 0 ─ 1 ─ 2          (lineage A)
        //               ╲
        //                3 ─ 4          (lineage B)
        //
        // Clusters 0 and 1 are on BOTH lineages (count = 2, trunk).
        // Clusters 2, 3, 4 are on exactly one lineage (single-lineage).
        let lineages = vec![vec![0, 1, 2], vec![0, 1, 3, 4]];
        let counts = cluster_lineage_count(&lineages, 5);
        assert_eq!(counts, vec![2, 2, 1, 1, 1]);

        //  binary trunk definition: is_trunk(c) = count(c) >= 2.
        let is_trunk: Vec<bool> = counts.iter().map(|&c| c >= 2).collect();
        assert_eq!(is_trunk, vec![true, true, false, false, false]);
    }

    #[test]
    fn cluster_lineage_count_isolated_cluster_is_zero() {
        // Cluster 4 is on no lineage (e.g. dropped from all root→leaf
        // paths).  Must report count = 0, not panic.
        let lineages = vec![vec![0, 1, 2], vec![0, 3]];
        let counts = cluster_lineage_count(&lineages, 5);
        assert_eq!(counts, vec![2, 1, 1, 1, 0]);
    }

    #[test]
    fn cluster_lineage_count_ignores_out_of_range_ids() {
        // Ignore cluster ids outside the declared range while retaining counts
        // for valid segments of each path.
        let lineages = vec![
            vec![0, 1, 99], // 99 out of range
            vec![0, 2],
        ];
        let counts = cluster_lineage_count(&lineages, 3);
        assert_eq!(counts, vec![2, 1, 1]);
    }

    /// A shared-trunk cell displaced toward lineage A must receive more weight
    /// on A than on B. The path-membership weighting function returns equal
    /// weights because both lineages contain the trunk cluster.
    #[test]
    fn trunk_aware_uses_geometry_on_shared_trunk() {
        let (centroids, lineages, _n_clusters, n_dims) = bifurcation_fixture();
        // Cell at (1.2, 0.4): cluster_label = 1 (trunk).  Slightly toward
        // the A end (positive y).  σ = 1.0, default.
        let embedding = vec![1.2f32, 0.4];
        let cluster_labels = vec![1usize];
        let sigma = 1.0f32;

        let new = compute_lineage_weights_trunk_aware(
            &embedding,
            &cluster_labels,
            1,
            n_dims,
            &centroids,
            &lineages,
            sigma,
        );
        assert_eq!(new.len(), 2);

        // Geometry should break the equal path-membership tie on the trunk.
        let (w_a, w_b) = (new[0], new[1]);
        assert!((w_a + w_b - 1.0).abs() < 1e-5, "row must normalise to 1");
        assert!(
            (w_a - 0.5).abs() > 1e-3,
            "trunk-aware must not produce uniform (0.5, 0.5) at off-centroid trunk cell; got w_a = {w_a}"
        );

        // Specifically: positive-y cell should lean toward A (which
        // ends at (2, 1)), so w_a > w_b.
        assert!(
            w_a > w_b,
            "off-trunk-toward-A cell should weight A higher; got ({w_a}, {w_b})"
        );

        // Path-membership weighting returns (0.5, 0.5) because both lineages
        // contain cluster 1.
        let path_membership = compute_lineage_weights(
            &embedding,
            &cluster_labels,
            1,
            n_dims,
            &centroids,
            &lineages,
            sigma,
        );
        assert!(
            (path_membership[0] - 0.5).abs() < 1e-6 && (path_membership[1] - 0.5).abs() < 1e-6,
            "path-membership weighting must return (0.5, 0.5); got ({}, {})",
            path_membership[0],
            path_membership[1]
        );
    }

    /// A single-lineage cluster on its own path retains a hard weight of 1.0.
    #[test]
    fn trunk_aware_keeps_hard_one_on_single_lineage_path() {
        let (centroids, lineages, _, n_dims) = bifurcation_fixture();
        // Cell sitting exactly on cluster 2's centroid; cluster 2 is
        // single-lineage (only on A's path [0, 1, 2]).
        let embedding = vec![2.0f32, 1.0];
        let cluster_labels = vec![2usize];
        let sigma = 1.0f32;

        let w = compute_lineage_weights_trunk_aware(
            &embedding,
            &cluster_labels,
            1,
            n_dims,
            &centroids,
            &lineages,
            sigma,
        );

        // unnormalised:  w_A = 1.0 (hard boost, cluster 2 is in A path),
        //                w_B = exp(-min_d²_B / 2σ²).  After normalisation
        //                w_A > w_B and both > 0.
        assert!(
            w[0] > w[1],
            "A-only cluster on A path must out-weigh B (got A={}, B={})",
            w[0],
            w[1]
        );
        assert!((w[0] + w[1] - 1.0).abs() < 1e-5);

        // Both weighting functions use the same single-lineage branch.
        let path_membership = compute_lineage_weights(
            &embedding,
            &cluster_labels,
            1,
            n_dims,
            &centroids,
            &lineages,
            sigma,
        );
        assert!(
            (w[0] - path_membership[0]).abs() < 1e-6 && (w[1] - path_membership[1]).abs() < 1e-6,
            "single-lineage on-path weights must match"
        );
    }

    /// Row-sum invariant: every row of the
    /// returned matrix sums to exactly 1.0 within float tolerance.
    #[test]
    fn trunk_aware_row_normalisation() {
        let (centroids, lineages, _, n_dims) = bifurcation_fixture();
        // Mixture: cells at every cluster centroid + two off-centroid
        // points covering trunk, A-end, B-end, plus a deliberately
        // off-path point.
        let embedding: Vec<f32> = vec![
            0.0, 0.0, // on cluster 0 (trunk)
            1.0, 0.0, // on cluster 1 (trunk)
            2.0, 1.0, // on cluster 2 (A end)
            2.0, -1.0, // on cluster 3 (B end)
            1.5, 0.5, // off-trunk-toward-A
            5.0, 0.0, // far off-path
        ];
        let cluster_labels = vec![0usize, 1, 2, 3, 1, 0];
        let n_cells = 6;
        let sigma = 1.0f32;

        let w = compute_lineage_weights_trunk_aware(
            &embedding,
            &cluster_labels,
            n_cells,
            n_dims,
            &centroids,
            &lineages,
            sigma,
        );
        assert_eq!(w.len(), n_cells * 2);
        for i in 0..n_cells {
            let row_sum = w[i * 2] + w[i * 2 + 1];
            assert!(
                (row_sum - 1.0).abs() < 1e-5,
                "row {i} must sum to 1; got {row_sum}"
            );
            assert!(w[i * 2].is_finite() && w[i * 2] >= 0.0);
            assert!(w[i * 2 + 1].is_finite() && w[i * 2 + 1] >= 0.0);
        }
    }

    #[test]
    fn trunk_aware_normalises_far_trunk_cells_without_underflow() {
        let (centroids, lineages, _, n_dims) = bifurcation_fixture();
        let embedding = vec![1000.0f32, 1000.0];
        let cluster_labels = vec![1usize];

        let weights = compute_lineage_weights_trunk_aware(
            &embedding,
            &cluster_labels,
            1,
            n_dims,
            &centroids,
            &lineages,
            1.0,
        );

        assert!(weights.iter().all(|weight| weight.is_finite()));
        assert!(weights.iter().all(|weight| *weight >= 0.0));
        assert!(
            (weights.iter().sum::<f32>() - 1.0).abs() < 1e-5,
            "far trunk cells must retain a normalised geometric assignment"
        );
    }

    /// σ sensitivity direction: smaller σ steepens the Gaussian decay
    /// on the trunk segment, so an off-trunk-toward-A cell's
    /// `w_A / w_B` ratio MUST grow as σ decreases.  This is the
    /// expected monotone behavior under bandwidth sensitivity.
    #[test]
    fn trunk_aware_smaller_sigma_steepens_decay() {
        let (centroids, lineages, _, n_dims) = bifurcation_fixture();
        let embedding = vec![1.2f32, 0.4];
        let cluster_labels = vec![1usize];

        let ratio_at_sigma = |sigma: f32| -> f32 {
            let w = compute_lineage_weights_trunk_aware(
                &embedding,
                &cluster_labels,
                1,
                n_dims,
                &centroids,
                &lineages,
                sigma,
            );
            w[0] / w[1]
        };

        let r_sharp = ratio_at_sigma(0.5);
        let r_default = ratio_at_sigma(1.0);
        let r_wide = ratio_at_sigma(2.0);

        // All three ratios exceed 1 because the positive-y cell leans toward A.
        assert!(r_sharp > 1.0 && r_default > 1.0 && r_wide > 1.0);

        // Monotone decreasing as σ widens.
        assert!(
            r_sharp > r_default,
            "σ=0.5 should steepen decay vs σ=1.0; got ratios {} vs {}",
            r_sharp,
            r_default
        );
        assert!(
            r_default > r_wide,
            "σ=1.0 should steepen decay vs σ=2.0; got ratios {} vs {}",
            r_default,
            r_wide
        );
    }

    /// When no trunk clusters exist, trunk-aware and path-membership weighting
    /// follow the same branches and return identical values.
    #[test]
    fn trunk_aware_matches_path_membership_weights_without_trunk() {
        // 2 disjoint lineages on disjoint clusters — no trunk anywhere.
        let centroids: Vec<Vec<f32>> = vec![
            vec![0.0, 0.0],
            vec![1.0, 0.0],
            vec![2.0, 1.0],
            vec![2.0, -1.0],
        ];
        let lineages = vec![vec![0, 1], vec![2, 3]];
        // counts = [1, 1, 1, 1], so every cluster is non-trunk.
        let embedding: Vec<f32> = vec![
            0.5, 0.0, // near lineage A
            2.0, 0.5, // near lineage B
            1.0, 0.5, // between
        ];
        let cluster_labels = vec![0usize, 2, 0];
        let n_cells = 3;
        let sigma = 1.0f32;

        let new = compute_lineage_weights_trunk_aware(
            &embedding,
            &cluster_labels,
            n_cells,
            2,
            &centroids,
            &lineages,
            sigma,
        );
        let path_membership = compute_lineage_weights(
            &embedding,
            &cluster_labels,
            n_cells,
            2,
            &centroids,
            &lineages,
            sigma,
        );

        assert_eq!(new.len(), path_membership.len());
        for (i, (a, b)) in new.iter().zip(path_membership.iter()).enumerate() {
            assert!((a - b).abs() < 1e-7,
                "no-trunk topology must give equal output at index {i}; trunk_aware={a}, path_membership={b}");
        }
    }
}
