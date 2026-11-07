//! R (extendr) bindings for LineaLinkR.
//!
//! Gated behind the `extendr` Cargo feature so that `cargo test --lib` does
//! not need R / extendr-api installed.  The R package build sets
//! `--features extendr` from `src/Makevars[.win]`.
//!
//! Each `#[extendr]` here corresponds 1:1 to one of the 7 public R functions
//! exported by `R/linealinkr.R`.  See `R/wrap__LineaLinkR.R` for the
//! `.Call` shim (regeneratable via `rextendr::document()`).
//!
//! ## Index / matrix conventions (boundary contract)
//!
//! - R passes **1-based** indices (cells, clusters, genes, terminal cells).
//!   Rust core operates on **0-based**.  This file does the −1 translation at
//!   the boundary, in exactly one place per index variable.  Downstream Rust
//!   code never sees 1-based indices.
//! - R matrices are **column-major**; LineaLinkR Rust core expects
//!   **row-major** flat slices.  The R wrappers in `R/linealinkr.R`
//!   transpose (`t()`) before passing in.  This file does NOT re-transpose.
//!
//! ## Boundary behavior
//!
//! The bridge forwards R arguments to the Rust core without adding wrapper
//! types or changing core semantics.

#[cfg(feature = "extendr")]
mod extendr_impl {
    use extendr_api::prelude::*;
    use std::hash::Hasher;

    use crate::ccc::lr_scoring::{
        compute_bin_bounds, compute_common_bin_bounds,
        compute_lr_scores_bins_fixed_scale_grid_with_bounds,
        compute_lr_scores_bins_with_scale_grid_opts, LrBinGrid,
    };
    use crate::ccc::taxonomy::LrBinIntermediate;
    use crate::ccc::trajectory_curve::{
        build_lr_trajectory_curve_from_grid, positive_score_from_grid,
        trajectory_keys_from_scores, TrajKey,
    };
    use crate::ccc::{compute_lr_scores, LrBinResult, LrPair};
    use crate::core::SparseMatrix;
    use crate::trajectory::{
        attractor_fate_probabilities_f64, cluster_centroids, compute_lineage_weights,
        compute_lineage_weights_trunk_aware, compute_lineage_weights_local_branch_residual,
        compute_pseudotime_from_hnsw_fuzzy,
        compute_pseudotime_from_hnsw_fuzzy_multi, count_knn_components, directed_markov,
        extract_lineages, fate_probabilities, knn_component_labels, mst_prim, pseudotime_de,
    };
    // `knn_graph_from_hnsw` is public in `trajectory::pseudotime` but is not
    // re-exported by `trajectory::mod.rs`; the full path preserves the public API.
    use crate::trajectory::pseudotime::{attractor_unreachable_cells, hnsw_cost_edges, knn_graph_from_hnsw};

    // ── helpers ─────────────────────────────────────────────────────────────

    /// Convert R dgCMatrix (Matrix package) → Rust SparseMatrix (CSC).
    /// `i` and `p` slots of dgCMatrix are 0-based, so no offset translation
    /// is needed here — only widening from `i32` to `usize` and from
    /// `f64` (R numeric) to `f32` (LineaLinkR core).
    fn r_dgcmatrix_to_sparse(mat: &S4) -> SparseMatrix {
        let dim_r = mat.get_slot("Dim").unwrap();
        let p_r = mat.get_slot("p").unwrap();
        let x_r = mat.get_slot("x").unwrap();
        let i_r = mat.get_slot("i").unwrap();
        let dim = dim_r.as_integer_slice().unwrap();
        let p = p_r.as_integer_slice().unwrap();
        let x = x_r.as_real_slice().unwrap();
        let i = i_r.as_integer_slice().unwrap();
        SparseMatrix::new_csc(
            i.iter().map(|&v| v as usize).collect(),
            p.iter().map(|&v| v as usize).collect(),
            x.iter().map(|&v| v as f32).collect(),
            dim[0] as usize,
            dim[1] as usize,
        )
    }

    /// Convert an R list of integer vectors (1-based gene row indices, one
    /// vector per ligand or receptor subunit set) into `Vec<Vec<usize>>`
    /// (0-based, Rust-native).  Used by both `rs_compute_lr_scores` and
    /// `rs_compute_temporal_lr`.
    fn list_of_int_to_rows0(list: &List) -> Vec<Vec<usize>> {
        list.iter()
            .map(|(_, rv)| {
                Integers::try_from(rv)
                    .unwrap()
                    .iter()
                    .map(|v| (v.inner() - 1) as usize)
                    .collect()
            })
            .collect()
    }

    #[inline]
    fn output_kish_ess(weight_sum: f64, squared_weight_sum: f64) -> f64 {
        if squared_weight_sum > 1e-12 {
            weight_sum * weight_sum / squared_weight_sum
        } else {
            0.0
        }
    }

    fn requested_lr_contexts_from_r(
        req_pair: &Integers,
        req_source: &Integers,
        req_target: &Integers,
    ) -> Option<Vec<(usize, usize, usize)>> {
        let n = req_pair.len();
        if n == 0 || req_source.len() != n || req_target.len() != n {
            return None;
        }
        let mut contexts: Vec<(usize, usize, usize)> = req_pair
            .iter()
            .zip(req_source.iter())
            .zip(req_target.iter())
            .map(|((p, s), t)| {
                (
                    (p.inner() - 1) as usize,
                    (s.inner() - 1) as usize,
                    (t.inner() - 1) as usize,
                )
            })
            .collect();
        contexts.sort_unstable();
        contexts.dedup();
        Some(contexts)
    }

    fn requested_trajectory_keys_from_r(
        req_pair: &Integers,
        req_source: &Integers,
        req_target: &Integers,
        req_lineage: &Integers,
    ) -> Option<Vec<TrajKey>> {
        let n = req_pair.len();
        if n == 0 || req_source.len() != n || req_target.len() != n || req_lineage.len() != n {
            return None;
        }
        let mut keys: Vec<TrajKey> = req_pair
            .iter()
            .zip(req_source.iter())
            .zip(req_target.iter())
            .zip(req_lineage.iter())
            .map(|(((p, s), t), l)| {
                (
                    (p.inner() - 1) as usize,
                    (s.inner() - 1) as usize,
                    (t.inner() - 1) as usize,
                    (l.inner() - 1) as usize,
                )
            })
            .collect();
        keys.sort_unstable();
        keys.dedup();
        Some(keys)
    }

    /// Convert the temporal scorer into one explicit bin representation.
    /// Complete mode expands every lineage and bin for LR/source/target contexts
    /// discovered by the prescreening scorer:
    /// structural invalid -> NA, structural valid/expression unsupported -> 0.
    fn temporal_scores_to_list(
        results: &[LrBinResult],
        lr_pairs: &[LrPair],
        grid: &LrBinGrid,
        inter: &LrBinIntermediate,
        min_pct: f32,
        min_n_eff: f32,
        return_complete_grid: bool,
        requested_contexts: Option<Vec<(usize, usize, usize)>>,
    ) -> List {
        let mut pair = Vec::<i32>::new();
        let mut source = Vec::<i32>::new();
        let mut target = Vec::<i32>::new();
        let mut lineage = Vec::<i32>::new();
        let mut bin = Vec::<i32>::new();
        let mut score = Vec::<f64>::new();
        let mut score_positive = Vec::<f64>::new();
        let mut n_eff = Vec::<f64>::new();
        let mut weight_mass = Vec::<f64>::new();
        let mut lineage_entropy = Vec::<f64>::new();
        let mut structural_valid = Vec::<i32>::new();
        let mut expression_supported = Vec::<i32>::new();

        if return_complete_grid {
            let mut contexts = if let Some(req) = requested_contexts {
                req
            } else {
                let prescreen_keys = trajectory_keys_from_scores(results);
                let mut derived: Vec<(usize, usize, usize)> = prescreen_keys
                    .iter()
                    .map(|key| (key.0, key.1, key.2))
                    .collect();
                derived.sort_unstable();
                derived.dedup();
                derived
            };
            contexts.sort_unstable();
            contexts.dedup();
            let slb = grid.slb();
            let sb = grid.sb();
            for context in contexts {
                if context.0 >= lr_pairs.len()
                    || context.1 >= grid.n_types
                    || context.2 >= grid.n_types
                {
                    continue;
                }
                for lineage_index in 0..grid.n_lineages {
                    let curve = build_lr_trajectory_curve_from_grid(
                        &lr_pairs[context.0],
                        grid,
                        context.1,
                        context.2,
                        lineage_index,
                        min_pct,
                        min_n_eff,
                    );
                    for (bin_index, cell) in curve.bins.iter().enumerate() {
                        let source_index = context.1 * slb + lineage_index * sb + bin_index;
                        let target_index = context.2 * slb + lineage_index * sb + bin_index;
                        let ess_source =
                            output_kish_ess(grid.wsums[source_index], grid.w2sums[source_index]);
                        let ess_target =
                            output_kish_ess(grid.wsums[target_index], grid.w2sums[target_index]);
                        let entropy_source = if inter.wsums[source_index] > 1e-8 {
                            inter.hwsums[source_index] / inter.wsums[source_index]
                        } else {
                            0.0
                        };
                        let entropy_target = if inter.wsums[target_index] > 1e-8 {
                            inter.hwsums[target_index] / inter.wsums[target_index]
                        } else {
                            0.0
                        };
                        pair.push(context.0 as i32 + 1);
                        source.push(context.1 as i32 + 1);
                        target.push(context.2 as i32 + 1);
                        lineage.push(lineage_index as i32 + 1);
                        bin.push(bin_index as i32 + 1);
                        score.push(if cell.structural_valid {
                            cell.score as f64
                        } else {
                            f64::NAN
                        });
                        score_positive.push(if cell.structural_valid {
                            if cell.expression_supported {
                                positive_score_from_grid(
                                    &lr_pairs[context.0], grid, context.1, context.2,
                                    lineage_index, bin_index,
                                ) as f64
                            } else {
                                0.0
                            }
                        } else {
                            f64::NAN
                        });
                        n_eff.push(ess_source.min(ess_target));
                        weight_mass.push(grid.wsums[source_index].min(grid.wsums[target_index]));
                        lineage_entropy.push(entropy_source.max(entropy_target));
                        structural_valid.push(cell.structural_valid as i32);
                        expression_supported.push(cell.expression_supported as i32);
                    }
                }
            }
        } else {
            for row in results {
                let valid = row.n_eff_cells >= min_n_eff;
                pair.push(row.pair_idx as i32 + 1);
                source.push(row.source as i32 + 1);
                target.push(row.target as i32 + 1);
                lineage.push(row.lineage as i32 + 1);
                bin.push(row.bin as i32 + 1);
                score.push(if valid { row.score as f64 } else { f64::NAN });
                score_positive.push(if valid {
                    positive_score_from_grid(
                        &lr_pairs[row.pair_idx], grid, row.source, row.target,
                        row.lineage, row.bin,
                    ) as f64
                } else {
                    f64::NAN
                });
                n_eff.push(row.n_eff_cells as f64);
                weight_mass.push(row.weight_mass_cells as f64);
                lineage_entropy.push(row.lineage_entropy as f64);
                structural_valid.push(valid as i32);
                expression_supported.push(valid as i32);
            }
        }

        List::from_names_and_values(
            [
                "pair_idx",
                "source",
                "target",
                "lineage",
                "bin",
                "score",
                "score_positive",
                "n_eff_cells",
                "weight_mass_cells",
                "lineage_entropy",
                "structural_valid",
                "expression_supported",
            ],
            [
                Integers::from_values(pair).into_robj(),
                Integers::from_values(source).into_robj(),
                Integers::from_values(target).into_robj(),
                Integers::from_values(lineage).into_robj(),
                Integers::from_values(bin).into_robj(),
                Doubles::from_values(score).into_robj(),
                Doubles::from_values(score_positive).into_robj(),
                Doubles::from_values(n_eff).into_robj(),
                Doubles::from_values(weight_mass).into_robj(),
                Doubles::from_values(lineage_entropy).into_robj(),
                Integers::from_values(structural_valid).into_robj(),
                Integers::from_values(expression_supported).into_robj(),
            ],
        )
        .unwrap()
    }

    /// Hash the shared scoring/evaluation configuration used by FakeLR and
    /// temporal permutation. Mechanism-specific controls (number of null draws,
    /// candidate pool, seed) are intentionally excluded so independently run
    /// calibration layers can still be checked for scoring compatibility.
    fn calibration_config_hash(
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
        stat_code: usize,
        top_k: usize,
    ) -> String {
        let mut hasher = crate::core::fxhash::FxHasher::default();
        hasher.write_u64(0x4c4c_4346_4730_3031);
        hasher.write_u64(expr.nrows as u64);
        hasher.write_u64(expr.ncols as u64);
        hasher.write_u64(expr.rows.len() as u64);
        for &value in &expr.rows {
            hasher.write_u64(value as u64);
        }
        hasher.write_u64(expr.col_ptrs.len() as u64);
        for &value in &expr.col_ptrs {
            hasher.write_u64(value as u64);
        }
        hasher.write_u64(expr.values.len() as u64);
        for &value in &expr.values {
            hasher.write_u32(value.to_bits());
        }
        hasher.write_u64(cell_types.len() as u64);
        for &value in cell_types {
            hasher.write_u64(value as u64);
        }
        hasher.write_u64(n_types as u64);
        hasher.write_u64(lr_pairs.len() as u64);
        for pair in lr_pairs {
            hasher.write_u64(pair.ligand_rows.len() as u64);
            for &value in &pair.ligand_rows {
                hasher.write_u64(value as u64);
            }
            hasher.write_u64(pair.receptor_rows.len() as u64);
            for &value in &pair.receptor_rows {
                hasher.write_u64(value as u64);
            }
        }
        hasher.write_u64(pseudotime.len() as u64);
        for &value in pseudotime {
            hasher.write_u32(value.to_bits());
        }
        hasher.write_u64(lineage_weights.len() as u64);
        for &value in lineage_weights {
            hasher.write_u32(value.to_bits());
        }
        hasher.write_u64(n_lineages as u64);
        hasher.write_u64(n_bins as u64);
        hasher.write_u32(min_pct.to_bits());
        hasher.write_u32(min_n_eff.to_bits());
        hasher.write_u64(min_valid_bins as u64);
        hasher.write_u32(min_valid_fraction.to_bits());
        hasher.write_u64(stat_code as u64);
        hasher.write_u64(top_k as u64);
        format!("{:016x}", hasher.finish())
    }

    // ── 1. ComputePseudotimeFuzzy ───────────────────────────────────────────

    /// Pseudotime via UMAP-fuzzy-weighted Dijkstra on a pre-computed kNN graph.
    ///
    /// @param nn_idx integer vector, row-major flat \[n_cells × k\], 1-based
    /// @param nn_dist numeric vector, row-major flat \[n_cells × k\], L2 distances
    /// @param n_cells integer
    /// @param k integer
    /// @param start_cell integer, 1-based root cell index
    /// @return numeric vector of length n_cells, pseudotime in \[0,1\]; NaN for
    ///   cells disconnected from the root.
    #[extendr]
    pub fn rs_compute_pseudotime_fuzzy(
        nn_idx: Integers,
        nn_dist: Doubles,
        n_cells: i32,
        k: i32,
        start_cell: i32,
    ) -> Doubles {
        let nc = n_cells as usize;
        let kk = k as usize;
        let idx: Vec<usize> = nn_idx.iter().map(|v| (v.inner() - 1) as usize).collect();
        let dist: Vec<f32> = nn_dist.iter().map(|v| v.inner() as f32).collect();
        let pt = compute_pseudotime_from_hnsw_fuzzy(&idx, &dist, nc, kk, (start_cell - 1) as usize);
        Doubles::from_values(pt.into_iter())
    }

    #[extendr]
    pub fn rs_compute_pseudotime_fuzzy_multi(
        nn_idx: Integers,
        nn_dist: Doubles,
        n_cells: i32,
        k: i32,
        start_cells: Integers,
    ) -> Doubles {
        let nc = n_cells as usize;
        let kk = k as usize;
        let idx: Vec<usize> = nn_idx.iter().map(|v| (v.inner() - 1) as usize).collect();
        let dist: Vec<f32> = nn_dist.iter().map(|v| v.inner() as f32).collect();
        let starts: Vec<usize> = start_cells
            .iter()
            .map(|v| (v.inner() - 1) as usize)
            .collect();
        let pt = compute_pseudotime_from_hnsw_fuzzy_multi(&idx, &dist, nc, kk, &starts);
        Doubles::from_values(pt.into_iter().map(|v| v as f64))
    }

    #[extendr]
    pub fn rs_hnsw_cost_edges(nn_idx: Integers, nn_dist: Doubles, n_cells: i32, k: i32) -> List {
        let nc = n_cells as usize;
        let kk = k as usize;
        let idx: Vec<usize> = nn_idx.iter().map(|v| (v.inner() - 1) as usize).collect();
        let dist: Vec<f32> = nn_dist.iter().map(|v| v.inner() as f32).collect();
        let edges = hnsw_cost_edges(&idx, &dist, nc, kk);
        List::from_names_and_values(
            ["from", "to", "fuzzy_cost", "raw_distance"],
            [
                Integers::from_values(edges.iter().map(|edge| (edge.from + 1) as i32)).into_robj(),
                Integers::from_values(edges.iter().map(|edge| (edge.to + 1) as i32)).into_robj(),
                Doubles::from_values(edges.iter().map(|edge| edge.fuzzy_cost as f64)).into_robj(),
                Doubles::from_values(edges.iter().map(|edge| edge.raw_distance as f64)).into_robj(),
            ],
        )
        .unwrap()
    }

    #[extendr]
    pub fn rs_knn_component_count(nn_idx: Integers, n_cells: i32, k: i32) -> i32 {
        let idx: Vec<usize> = nn_idx.iter().map(|v| (v.inner() - 1) as usize).collect();
        count_knn_components(&idx, n_cells as usize, k as usize) as i32
    }

    #[extendr]
    pub fn rs_knn_component_labels(nn_idx: Integers, n_cells: i32, k: i32) -> Integers {
        let idx: Vec<usize> = nn_idx.iter().map(|v| (v.inner() - 1) as usize).collect();
        let labels = knn_component_labels(&idx, n_cells as usize, k as usize);
        Integers::from_values(labels.into_iter().map(|label| (label + 1) as i32))
    }

    #[extendr]
    pub fn rs_attractor_unreachable_cells(
        nn_idx: Integers,
        nn_dist: Doubles,
        terminal_cells: Integers,
        n_cells: i32,
        k: i32,
    ) -> Integers {
        let indices: Vec<usize> = nn_idx.iter().map(|value| (value.inner() - 1) as usize).collect();
        let distances: Vec<f64> = nn_dist.iter().map(|value| value.inner()).collect();
        let terminals: Vec<usize> = terminal_cells.iter().map(|value| (value.inner() - 1) as usize).collect();
        Integers::from_values(attractor_unreachable_cells(
            &indices, &distances, n_cells as usize, k as usize, &terminals,
        ).into_iter().map(|cell| (cell + 1) as i32))
    }

    #[extendr]
    pub fn rs_compute_attractor_fate_f64(
        nn_idx: Integers,
        nn_dist: Doubles,
        terminal_cells: Integers,
        n_cells: i32,
        k: i32,
        max_iter: i32,
        tol: f64,
    ) -> Doubles {
        let indices: Vec<usize> = nn_idx
            .iter()
            .map(|value| (value.inner() - 1) as usize)
            .collect();
        let distances: Vec<f64> = nn_dist.iter().map(|value| value.inner()).collect();
        let terminals: Vec<usize> = terminal_cells
            .iter()
            .map(|value| (value.inner() - 1) as usize)
            .collect();
        let probability = attractor_fate_probabilities_f64(
            &indices,
            &distances,
            n_cells as usize,
            k as usize,
            &terminals,
            max_iter as usize,
            tol,
        );
        match probability {
            Ok(values) => Doubles::from_values(values),
            Err(message) => {
                drop(indices);
                drop(distances);
                drop(terminals);
                throw_r_error(message)
            }
        }
    }

    extern "C" {
        fn llr_absorption_lu64(
            size: i32, outer: *const i32, inner: *const i32, values: *const f64,
            rhs: *const f64, n_rhs: i32, output: *mut f64,
            error: *mut std::os::raw::c_char, error_size: i32,
        ) -> i32;
        fn llr_state_clock(
            n_cells: i32, n_dims: i32, expression: *const f64,
            outer: *const i32, inner: *const i32, edge_weight: *const f64,
            input_state: *const i32, anchor_order: *const f64,
            n_states: i32, reference_parent: *const i32, threads: i32,
            output_time: *mut f64, output_chart: *mut i32, output_component: *mut i32,
            output_info: *mut f64, output_diagnostics: *mut f64, output_ncharts: *mut i32,
            error: *mut std::os::raw::c_char, error_size: i32,
        ) -> i32;
    }

    /// Unlike the trajectory kernels, this solver keeps R's column-major layout.
    /// dgCMatrix CSC indices are already zero-based and must not be shifted.
    #[extendr]
    pub fn rs_sparse_absorption_lu64(system: S4, rhs: Doubles, n_rhs: i32) -> Doubles {
        let dimensions = system.get_slot("Dim").unwrap();
        let outer = system.get_slot("p").unwrap();
        let inner = system.get_slot("i").unwrap();
        let values = system.get_slot("x").unwrap();
        let dimensions = dimensions.as_integer_slice().unwrap();
        let outer = outer.as_integer_slice().unwrap();
        let inner = inner.as_integer_slice().unwrap();
        let values = values.as_real_slice().unwrap();
        let rhs: Vec<f64> = rhs.iter().map(|value| value.inner()).collect();
        if dimensions[0] != dimensions[1] || n_rhs < 1
            || rhs.len() != dimensions[0] as usize * n_rhs as usize {
            throw_r_error("Invalid sparse absorption system dimensions");
        }
        let mut output = vec![0.0; rhs.len()];
        let mut error = [0 as std::os::raw::c_char; 1024];
        let status = unsafe {
            llr_absorption_lu64(
                dimensions[0], outer.as_ptr(), inner.as_ptr(), values.as_ptr(),
                rhs.as_ptr(), n_rhs, output.as_mut_ptr(), error.as_mut_ptr(), error.len() as i32,
            )
        };
        if status != 0 {
            let message = unsafe { std::ffi::CStr::from_ptr(error.as_ptr()) }
                .to_string_lossy().into_owned();
            drop(rhs);
            drop(output);
            throw_r_error(message);
        }
        Doubles::from_values(output)
    }

    /// Expression-constrained correction of a globally anchored coordinate.
    /// Calibration patches and the supplied reference states are distinct objects.
    #[extendr]
    pub fn rs_compute_state_clock(
        expression: Doubles, graph: S4, input_state: Integers,
        anchor_order: Doubles, reference_parent: Integers, n_dims: i32, threads: i32,
    ) -> Result<List> {
        let dimensions = graph.get_slot("Dim").unwrap();
        let outer = graph.get_slot("p").unwrap();
        let inner = graph.get_slot("i").unwrap();
        let edge_weight = graph.get_slot("x").unwrap();
        let dimensions = dimensions.as_integer_slice().unwrap();
        let outer = outer.as_integer_slice().unwrap();
        let inner = inner.as_integer_slice().unwrap();
        let edge_weight = edge_weight.as_real_slice().unwrap();
        let n_cells = dimensions[0] as usize;
        let expression: Vec<f64> = expression.iter().map(|x| x.inner()).collect();
        let input_state: Vec<i32> = input_state.iter().map(|x| x.inner()).collect();
        let anchor_order: Vec<f64> = anchor_order.iter().map(|x| x.inner()).collect();
        let parent: Vec<i32> = reference_parent.iter().map(|x| x.inner()).collect();
        if dimensions[0] != dimensions[1] || n_dims < 1 || n_cells < 2 || parent.is_empty()
            || expression.len() != n_cells * n_dims as usize
            || input_state.len() != n_cells || anchor_order.len() != n_cells
            || input_state.iter().any(|&x| x < 1 || x as usize > parent.len())
            || parent.iter().any(|&x| x < 0 || x as usize > parent.len()) {
            return Err(Error::Other("Invalid anchored clock or reference forest dimensions".into()));
        }
        let mut time = vec![0.0; n_cells];
        let mut labels = vec![0i32; n_cells];
        let mut components = vec![0i32; n_cells];
        let mut chart_info = vec![0.0; 8 * n_cells];
        let mut diagnostics = vec![0.0; 10];
        let mut n_charts = 0i32;
        let mut error = [0 as std::os::raw::c_char; 1024];
        let status = unsafe {
            llr_state_clock(
                n_cells as i32, n_dims, expression.as_ptr(), outer.as_ptr(), inner.as_ptr(),
                edge_weight.as_ptr(), input_state.as_ptr(), anchor_order.as_ptr(),
                parent.len() as i32, parent.as_ptr(), threads,
                time.as_mut_ptr(), labels.as_mut_ptr(), components.as_mut_ptr(),
                chart_info.as_mut_ptr(), diagnostics.as_mut_ptr(), &mut n_charts,
                error.as_mut_ptr(), error.len() as i32,
            )
        };
        if status != 0 {
            let message = unsafe { std::ffi::CStr::from_ptr(error.as_ptr()) }
                .to_string_lossy().into_owned();
            return Err(Error::Other(message));
        }
        components.truncate(n_charts as usize);
        chart_info.truncate(8 * n_charts as usize);
        List::from_names_and_values(
            ["pseudotime", "calibration_chart", "chart_components", "chart_info", "diagnostics"],
            [Doubles::from_values(time).into_robj(), Integers::from_values(labels).into_robj(),
             Integers::from_values(components).into_robj(), Doubles::from_values(chart_info).into_robj(),
             Doubles::from_values(diagnostics).into_robj()],
        )
    }
    // ── 2. ExtractLineages ──────────────────────────────────────────────────

    /// Extract root→leaf lineage paths from cluster centroids via MST.
    /// Internally: cluster_centroids → mst_prim → extract_lineages.
    ///
    /// @param embedding numeric vector, row-major flat \[n_cells × n_dims\]
    /// @param cluster_labels integer vector of length n_cells, 1-based
    /// @param n_cells, n_dims, n_clusters integers
    /// @param start_cluster integer, 1-based root cluster index
    /// @return list of integer vectors; each element is one lineage path
    ///   (sequence of 1-based cluster ids from root to leaf).
    #[extendr]
    pub fn rs_extract_lineages(
        embedding: Doubles,
        cluster_labels: Integers,
        n_cells: i32,
        n_dims: i32,
        n_clusters: i32,
        start_cluster: i32,
    ) -> List {
        let nc = n_cells as usize;
        let nd = n_dims as usize;
        let nclust = n_clusters as usize;
        let emb: Vec<f32> = embedding.iter().map(|v| v.inner() as f32).collect();
        let cl: Vec<usize> = cluster_labels
            .iter()
            .map(|v| (v.inner() - 1) as usize)
            .collect();
        let cents = cluster_centroids(&emb, &cl, nc, nd, nclust);
        let edges = mst_prim(&cents);
        let lins = extract_lineages(&edges, nclust, (start_cluster - 1) as usize);
        List::from_values(
            lins.into_iter().map(|path| {
                Integers::from_values(path.into_iter().map(|c| c as i32 + 1)).into_robj()
            }),
        )
    }

    // ── 3. ComputeLineageWeights ────────────────────────────────────────────

    /// Soft per-cell lineage weights along each lineage's polyline.
    /// Internally re-runs cluster_centroids → mst_prim → extract_lineages
    /// to keep the R-facing signature (start_cluster) free of upstream
    /// lineage objects.  Cost: O(n_clusters · n_dims + |MST|); negligible
    /// vs. the cell-level kernel pass below.
    ///
    /// Clusters on a lineage path receive a hard weight of 1.0. For geometric
    /// weighting on shared trunk clusters, use `ComputeLineageWeightsTrunkAware`.
    ///
    /// @return numeric vector, row-major flat \[n_cells × n_lineages\];
    ///   the R wrapper reshapes to matrix.  n_lineages can be recovered by
    ///   the caller from a separate `rs_extract_lineages` call.
    #[extendr]
    pub fn rs_compute_lineage_weights(
        embedding: Doubles,
        cluster_labels: Integers,
        n_cells: i32,
        n_dims: i32,
        n_clusters: i32,
        start_cluster: i32,
        sigma: f64,
    ) -> Doubles {
        let nc = n_cells as usize;
        let nd = n_dims as usize;
        let nclust = n_clusters as usize;
        let emb: Vec<f32> = embedding.iter().map(|v| v.inner() as f32).collect();
        let cl: Vec<usize> = cluster_labels
            .iter()
            .map(|v| (v.inner() - 1) as usize)
            .collect();
        let cents = cluster_centroids(&emb, &cl, nc, nd, nclust);
        let edges = mst_prim(&cents);
        let lins = extract_lineages(&edges, nclust, (start_cluster - 1) as usize);
        let w = compute_lineage_weights(&emb, &cl, nc, nd, &cents, &lins, sigma as f32);
        Doubles::from_values(w.into_iter().map(|v| v as f64))
    }

    // ── 3b. ComputeLineageWeightsTrunkAware ─────────────────────────────────

    /// Trunk-aware geometric lineage weighting. Same
    /// signature as `rs_compute_lineage_weights` so R callers switch
    /// by changing only the function name.
    ///
    /// Behaviour difference: trunk clusters (appearing on >= 2 lineage
    /// paths) use pure Gaussian `exp(-d^2/2σ^2)` rather than a hard membership
    /// weight, so `H̃_l(b)` reflects real geometry on trunk
    /// instead of being pinned near `log K / log L`.  Single-lineage
    /// on-path cells keep hard 1.0; off-path cells are geometric.
    ///
    /// `rs_compute_lineage_weights` remains available for callers that require
    /// path-membership weighting.
    #[extendr]
    pub fn rs_compute_lineage_weights_trunk_aware(
        embedding: Doubles,
        cluster_labels: Integers,
        n_cells: i32,
        n_dims: i32,
        n_clusters: i32,
        start_cluster: i32,
        sigma: f64,
    ) -> Doubles {
        let nc = n_cells as usize;
        let nd = n_dims as usize;
        let nclust = n_clusters as usize;
        let emb: Vec<f32> = embedding.iter().map(|v| v.inner() as f32).collect();
        let cl: Vec<usize> = cluster_labels
            .iter()
            .map(|v| (v.inner() - 1) as usize)
            .collect();
        let cents = cluster_centroids(&emb, &cl, nc, nd, nclust);
        let edges = mst_prim(&cents);
        let lins = extract_lineages(&edges, nclust, (start_cluster - 1) as usize);
        let w = compute_lineage_weights_trunk_aware(&emb, &cl, nc, nd, &cents, &lins, sigma as f32);
        Doubles::from_values(w.into_iter().map(|v| v as f64))
    }

    /// Extract root->leaf lineage paths from an explicit cluster tree.
    ///
    /// Variant A uses this with a tree built from the native cell kNN graph
    /// instead of the centroid MST, so paths, weights and fork placement all
    /// come from one observed-support structure.
    ///
    /// @param edges_from,edges_to 1-based cluster ids of the tree edges.
    /// @param n_clusters integer.
    /// @param start_cluster 1-based root cluster.
    /// @return List of integer vectors (1-based cluster ids, root first).
    #[extendr]
    pub fn rs_extract_lineages_from_edges(
        edges_from: Integers,
        edges_to: Integers,
        n_clusters: i32,
        start_cluster: i32,
    ) -> List {
        let nclust = n_clusters as usize;
        let edges: Vec<(usize, usize, f32)> = edges_from
            .iter()
            .zip(edges_to.iter())
            .map(|(a, b)| {
                ((a.inner() - 1) as usize, (b.inner() - 1) as usize, 1.0f32)
            })
            .collect();
        let lins = extract_lineages(&edges, nclust, (start_cluster - 1) as usize);
        List::from_values(lins.into_iter().map(|path| {
            Integers::from_values(
                path.into_iter().map(|c| (c + 1) as i32).collect::<Vec<_>>(),
            )
        }))
    }

    /// Trunk-aware lineage weights over an explicit cluster tree.
    ///
    /// Same weighting arithmetic as `rs_compute_lineage_weights_trunk_aware`;
    /// only the tree source changes (caller-provided edges).
    #[extendr]
    pub fn rs_compute_lineage_weights_trunk_aware_from_edges(
        embedding: Doubles,
        cluster_labels: Integers,
        n_cells: i32,
        n_dims: i32,
        n_clusters: i32,
        edges_from: Integers,
        edges_to: Integers,
        start_cluster: i32,
        sigma: f64,
    ) -> Doubles {
        let nc = n_cells as usize;
        let nd = n_dims as usize;
        let nclust = n_clusters as usize;
        let emb: Vec<f32> = embedding.iter().map(|v| v.inner() as f32).collect();
        let cl: Vec<usize> = cluster_labels
            .iter()
            .map(|v| (v.inner() - 1) as usize)
            .collect();
        let edges: Vec<(usize, usize, f32)> = edges_from
            .iter()
            .zip(edges_to.iter())
            .map(|(a, b)| {
                ((a.inner() - 1) as usize, (b.inner() - 1) as usize, 1.0f32)
            })
            .collect();
        let cents = cluster_centroids(&emb, &cl, nc, nd, nclust);
        let lins = extract_lineages(&edges, nclust, (start_cluster - 1) as usize);
        let w = compute_lineage_weights_trunk_aware(&emb, &cl, nc, nd, &cents, &lins, sigma as f32);
        Doubles::from_values(w.into_iter().map(|v| v as f64))
    }

    #[extendr]
    pub fn rs_compute_lineage_weights_local_branch_residual(
        embedding: Doubles,
        cluster_labels: Integers,
        n_cells: i32,
        n_dims: i32,
        n_clusters: i32,
        lineage_paths: List,
        priming: bool,
    ) -> Doubles {
        let nc = n_cells as usize;
        let nd = n_dims as usize;
        let nclust = n_clusters as usize;
        let emb: Vec<f32> = embedding.iter().map(|v| v.inner() as f32).collect();
        let cl: Vec<usize> = cluster_labels
            .iter()
            .map(|v| (v.inner() - 1) as usize)
            .collect();
        let paths = list_of_int_to_rows0(&lineage_paths);
        let cents = cluster_centroids(&emb, &cl, nc, nd, nclust);
        let weights = compute_lineage_weights_local_branch_residual(
            &emb, &cl, nc, nd, &cents, &paths, priming,
        );
        Doubles::from_values(weights.into_iter().map(|v| v as f64))
    }

    // ── 4. ComputeFateProbabilities ─────────────────────────────────────────

    /// Fate probabilities via absorbing-state Markov power iteration on a
    /// pseudotime-soft-gated kNN graph.
    /// Internally: knn_graph_from_hnsw → directed_markov → fate_probabilities.
    ///
    /// @param nn_idx integer vector, row-major flat \[n_cells × k\], 1-based
    /// @param nn_dist numeric vector, row-major flat \[n_cells × k\]
    /// @param pseudotime numeric vector of length n_cells; NaN cells receive
    ///   self-loops and are excluded from diffusion (matches Palantir).
    /// @param terminal_cells integer vector of 1-based terminal cell indices
    /// @param n_cells, k integers
    /// @param bandwidth numeric, sigmoid soft-gate bandwidth on pseudotime Δ
    /// @param max_iter, tol convergence params
    /// @return numeric vector, row-major flat \[n_cells × n_terminals\]
    #[extendr]
    pub fn rs_compute_fate_probabilities(
        nn_idx: Integers,
        nn_dist: Doubles,
        pseudotime: Doubles,
        terminal_cells: Integers,
        n_cells: i32,
        k: i32,
        bandwidth: f64,
        max_iter: i32,
        tol: f64,
    ) -> Doubles {
        let nc = n_cells as usize;
        let kk = k as usize;
        let idx: Vec<usize> = nn_idx.iter().map(|v| (v.inner() - 1) as usize).collect();
        let dist: Vec<f32> = nn_dist.iter().map(|v| v.inner() as f32).collect();
        let pt: Vec<f32> = pseudotime.iter().map(|v| v.inner() as f32).collect();
        let tc: Vec<usize> = terminal_cells
            .iter()
            .map(|v| (v.inner() - 1) as usize)
            .collect();
        let neighbors = knn_graph_from_hnsw(&idx, &dist, nc, kk);
        let trans = directed_markov(&neighbors, &pt, bandwidth as f32);
        let fp = fate_probabilities(nc, &trans, &tc, max_iter as usize, tol as f32);
        Doubles::from_values(fp.into_iter().map(|v| v as f64))
    }

    // ── 5. ComputeLRScores ──────────────────────────────────────────────────

    /// Cluster-level LR communication scores with permutation p-values
    /// (CellPhoneDB-style; SingleCellSignalR-style Michaelis-Menten saturation).
    ///
    /// @param expr dgCMatrix (genes × cells), log-normalised
    /// @param cell_types integer vector of length n_cells, 1-based cluster id
    /// @param n_types integer
    /// @param ligand_rows_list list of integer vectors (1-based gene row idx
    ///   per pair, ligand subunits — usually length 1 but can be a complex)
    /// @param receptor_rows_list same shape as `ligand_rows_list` for receptors
    /// @param n_perms integer
    /// @param seed numeric (R 32-bit integers cannot hold full u64 seeds)
    /// @param min_pct numeric, minimum per-cluster expressing fraction
    /// @return list with columns source/target/pair_idx (1-based),
    ///   score/pvalue/pct_ligand/pct_receptor; an R-side wrapper assembles
    ///   the data.frame.
    #[extendr]
    pub fn rs_compute_lr_scores(
        expr: S4,
        cell_types: Integers,
        n_types: i32,
        ligand_rows_list: List,
        receptor_rows_list: List,
        n_perms: i32,
        seed: f64,
        min_pct: f64,
    ) -> List {
        let sparse = r_dgcmatrix_to_sparse(&expr);
        let ct: Vec<usize> = cell_types
            .iter()
            .map(|v| (v.inner() - 1) as usize)
            .collect();
        let nt = n_types as usize;

        let lig = list_of_int_to_rows0(&ligand_rows_list);
        let rec = list_of_int_to_rows0(&receptor_rows_list);
        let lr_pairs: Vec<LrPair> = lig
            .into_iter()
            .zip(rec.into_iter())
            .map(|(l, r)| LrPair {
                ligand_rows: l,
                receptor_rows: r,
            })
            .collect();

        let results = compute_lr_scores(
            &sparse,
            &ct,
            nt,
            &lr_pairs,
            n_perms as usize,
            seed as u64,
            min_pct as f32,
        );

        let source: Vec<i32> = results.iter().map(|r| r.source as i32 + 1).collect();
        let target: Vec<i32> = results.iter().map(|r| r.target as i32 + 1).collect();
        let pair: Vec<i32> = results.iter().map(|r| r.pair_idx as i32 + 1).collect();
        let score: Vec<f64> = results.iter().map(|r| r.score as f64).collect();
        let pvalue: Vec<f64> = results.iter().map(|r| r.pvalue as f64).collect();
        let pct_lig: Vec<f64> = results.iter().map(|r| r.pct_ligand as f64).collect();
        let pct_rec: Vec<f64> = results.iter().map(|r| r.pct_receptor as f64).collect();

        List::from_names_and_values(
            [
                "source",
                "target",
                "pair_idx",
                "score",
                "pvalue",
                "pct_ligand",
                "pct_receptor",
            ],
            [
                Integers::from_values(source).into_robj(),
                Integers::from_values(target).into_robj(),
                Integers::from_values(pair).into_robj(),
                Doubles::from_values(score).into_robj(),
                Doubles::from_values(pvalue).into_robj(),
                Doubles::from_values(pct_lig).into_robj(),
                Doubles::from_values(pct_rec).into_robj(),
            ],
        )
        .unwrap()
    }

    // ── 6. ComputeTemporalLR ────────────────────────────────────────────────

    /// Pseudotime-binned LR scores per (lineage, bin) — the core LineaLinkR Module 1
    /// output used downstream by taxonomy and lineage-weighting analyses.
    ///
    /// Returns `n_eff_cells` and `lineage_entropy` per row so downstream R-
    /// side reliability gates (`min_n_eff_cells`) can be applied without
    /// re-touching Rust.
    ///
    /// This returns `LrBinResult.lineage_entropy`, defined as the maximum of
    /// source and target mean entropy. The per-lineage tensor `H̃_l(b)` is
    /// available through `ComputeTemporalLRWithIntermediate`.
    #[extendr]
    pub fn rs_compute_temporal_lr(
        expr: S4,
        cell_types: Integers,
        n_types: i32,
        ligand_rows_list: List,
        receptor_rows_list: List,
        pseudotime: Doubles,
        lineage_weights: Doubles,
        n_lineages: i32,
        n_bins: i32,
        min_pct: f64,
        min_n_eff: f64,
        return_complete_grid: bool,
        req_pair_idx: Integers,
        req_source: Integers,
        req_target: Integers,
        common_grid: bool,
    ) -> List {
        let sparse = r_dgcmatrix_to_sparse(&expr);
        let ct: Vec<usize> = cell_types
            .iter()
            .map(|v| (v.inner() - 1) as usize)
            .collect();
        let nt = n_types as usize;
        let pt: Vec<f32> = pseudotime.iter().map(|v| v.inner() as f32).collect();
        let lw: Vec<f32> = lineage_weights.iter().map(|v| v.inner() as f32).collect();
        let lig = list_of_int_to_rows0(&ligand_rows_list);
        let rec = list_of_int_to_rows0(&receptor_rows_list);
        let lr_pairs: Vec<LrPair> = lig
            .into_iter()
            .zip(rec.into_iter())
            .map(|(l, r)| LrPair {
                ligand_rows: l,
                receptor_rows: r,
            })
            .collect();

        // This binding returns scores only. The intermediate tensors are
        // exposed by `rs_compute_temporal_lr_with_inter`.
        let (results, inter, _hill_s, grid) = compute_lr_scores_bins_with_scale_grid_opts(
            &sparse,
            &ct,
            nt,
            &lr_pairs,
            &pt,
            &lw,
            n_lineages as usize,
            n_bins as usize,
            min_pct as f32,
            common_grid,
        );
        let requested_contexts =
            requested_lr_contexts_from_r(&req_pair_idx, &req_source, &req_target);
        temporal_scores_to_list(
            &results,
            &lr_pairs,
            &grid,
            &inter,
            min_pct as f32,
            min_n_eff as f32,
            return_complete_grid,
            requested_contexts,
        )
    }

    // ── 6b. ComputeTemporalLRWithIntermediate ──────────────────────────────

    /// Same computation as `rs_compute_temporal_lr` but also exposes the
    /// `LrBinIntermediate` (`hwsums` / `wsums`) needed by per-lineage entropy to build
    /// per-lineage entropy `H_l(b)`.
    ///
    /// Layout: `hwsums` and `wsums` are row-major flat of length
    /// `n_types * n_lineages * n_bins` with `ct` outermost,
    /// index = `ct * n_lineages * n_bins + l * n_bins + b`.  This matches
    /// the `LrBinIntermediate` struct (`src/ccc/taxonomy.rs`).
    ///
    /// @return list with two named elements:
    ///   - `scores`: same 11-column list `rs_compute_temporal_lr` returns,
    ///     including `structural_valid` and `expression_supported`
    ///   - `intermediate`: list(hwsums, wsums, n_types, n_lineages, n_bins);
    ///     consumed by `rs_compute_lineage_entropy_raw`.
    #[extendr]
    pub fn rs_compute_temporal_lr_with_inter(
        expr: S4,
        cell_types: Integers,
        n_types: i32,
        ligand_rows_list: List,
        receptor_rows_list: List,
        pseudotime: Doubles,
        lineage_weights: Doubles,
        n_lineages: i32,
        n_bins: i32,
        min_pct: f64,
        min_n_eff: f64,
        return_complete_grid: bool,
        req_pair_idx: Integers,
        req_source: Integers,
        req_target: Integers,
        common_grid: bool,
    ) -> List {
        let sparse = r_dgcmatrix_to_sparse(&expr);
        let ct: Vec<usize> = cell_types
            .iter()
            .map(|v| (v.inner() - 1) as usize)
            .collect();
        let nt = n_types as usize;
        let pt: Vec<f32> = pseudotime.iter().map(|v| v.inner() as f32).collect();
        let lw: Vec<f32> = lineage_weights.iter().map(|v| v.inner() as f32).collect();

        let lig = list_of_int_to_rows0(&ligand_rows_list);
        let rec = list_of_int_to_rows0(&receptor_rows_list);
        let lr_pairs: Vec<LrPair> = lig
            .into_iter()
            .zip(rec.into_iter())
            .map(|(l, r)| LrPair {
                ligand_rows: l,
                receptor_rows: r,
            })
            .collect();

        let (results, inter, hill_s, grid) = compute_lr_scores_bins_with_scale_grid_opts(
            &sparse,
            &ct,
            nt,
            &lr_pairs,
            &pt,
            &lw,
            n_lineages as usize,
            n_bins as usize,
            min_pct as f32,
            common_grid,
        );
        let requested_contexts =
            requested_lr_contexts_from_r(&req_pair_idx, &req_source, &req_target);
        let scores = temporal_scores_to_list(
            &results,
            &lr_pairs,
            &grid,
            &inter,
            min_pct as f32,
            min_n_eff as f32,
            return_complete_grid,
            requested_contexts,
        );

        let flat_bin_bounds: Vec<f64> = grid
            .bin_bounds
            .iter()
            .flat_map(|row| row.iter().map(|&value| value as f64))
            .collect();
        let intermediate = List::from_names_and_values(
            [
                "hwsums",
                "wsums",
                "n_types",
                "n_lineages",
                "n_bins",
                "hill_scale",
                "bin_bounds",
            ],
            [
                Doubles::from_values(inter.hwsums.into_iter()).into_robj(),
                Doubles::from_values(inter.wsums.into_iter()).into_robj(),
                Integers::from_values([inter.n_types as i32]).into_robj(),
                Integers::from_values([inter.n_lineages as i32]).into_robj(),
                Integers::from_values([inter.n_bins as i32]).into_robj(),
                Doubles::from_values([hill_s as f64]).into_robj(),
                Doubles::from_values(flat_bin_bounds.into_iter()).into_robj(),
            ],
        )
        .unwrap();

        List::from_names_and_values(
            ["scores", "intermediate", "hill_scale"],
            [
                scores.into_robj(),
                intermediate.into_robj(),
                Doubles::from_values([hill_s as f64]).into_robj(),
            ],
        )
        .unwrap()
    }

    /// Production equal-frequency temporal boundaries used by Rust scoring.
    #[extendr]
    pub fn rs_compute_bin_bounds(
        pseudotime: Doubles,
        lineage_weights: Doubles,
        n_lineages: i32,
        n_bins: i32,
    ) -> List {
        let pt: Vec<f32> = pseudotime.iter().map(|v| v.inner() as f32).collect();
        let lw: Vec<f32> = lineage_weights.iter().map(|v| v.inner() as f32).collect();
        let bounds = compute_bin_bounds(pt.len(), &pt, &lw, n_lineages as usize, n_bins as usize);
        List::from_values(
            bounds.into_iter().map(|row| {
                Doubles::from_values(row.into_iter().map(|value| value as f64)).into_robj()
            }),
        )
    }

    /// Shared temporal grid: one boundary vector common to every lineage.
    ///
    /// Returns the single `n_bins + 1` boundary row used for all lineages by
    /// `common_grid = TRUE` scoring.  Fork-local contrasts require this
    /// shared coordinate so a bin index denotes the same pseudotime interval
    /// on every branch.
    #[extendr]
    pub fn rs_compute_common_bin_bounds(
        pseudotime: Doubles,
        lineage_weights: Doubles,
        n_lineages: i32,
        n_bins: i32,
    ) -> Doubles {
        let pt: Vec<f32> = pseudotime.iter().map(|v| v.inner() as f32).collect();
        let lw: Vec<f32> = lineage_weights.iter().map(|v| v.inner() as f32).collect();
        let bounds = compute_common_bin_bounds(
            pt.len(),
            &pt,
            &lw,
            n_lineages as usize,
            n_bins as usize,
        );
        Doubles::from_values(bounds[0].iter().map(|&value| value as f64))
    }

    /// Temporal LR scoring on externally supplied pooled bin boundaries and
    /// one fixed positive Hill scale.
    #[extendr]
    pub fn rs_compute_temporal_lr_fixed_bounds(
        expr: S4,
        cell_types: Integers,
        n_types: i32,
        ligand_rows_list: List,
        receptor_rows_list: List,
        pseudotime: Doubles,
        lineage_weights: Doubles,
        n_lineages: i32,
        n_bins: i32,
        bin_bounds: Doubles,
        hill_scale: f64,
        min_pct: f64,
        min_n_eff: f64,
        return_complete_grid: bool,
        req_pair_idx: Integers,
        req_source: Integers,
        req_target: Integers,
    ) -> List {
        let sparse = r_dgcmatrix_to_sparse(&expr);
        let ct: Vec<usize> = cell_types
            .iter()
            .map(|v| (v.inner() - 1) as usize)
            .collect();
        let nt = n_types as usize;
        let pt: Vec<f32> = pseudotime.iter().map(|v| v.inner() as f32).collect();
        let lw: Vec<f32> = lineage_weights.iter().map(|v| v.inner() as f32).collect();
        let flat_bounds: Vec<f32> = bin_bounds.iter().map(|v| v.inner() as f32).collect();
        let bounds: Vec<Vec<f32>> = flat_bounds
            .chunks((n_bins as usize) + 1)
            .map(|row| row.to_vec())
            .collect();
        let lig = list_of_int_to_rows0(&ligand_rows_list);
        let rec = list_of_int_to_rows0(&receptor_rows_list);
        let lr_pairs: Vec<LrPair> = lig
            .into_iter()
            .zip(rec.into_iter())
            .map(|(l, r)| LrPair {
                ligand_rows: l,
                receptor_rows: r,
            })
            .collect();
        let (results, inter, grid) = compute_lr_scores_bins_fixed_scale_grid_with_bounds(
            &sparse,
            &ct,
            nt,
            &lr_pairs,
            &pt,
            &lw,
            n_lineages as usize,
            n_bins as usize,
            min_pct as f32,
            hill_scale as f32,
            Some(&bounds),
        );
        let requested_contexts =
            requested_lr_contexts_from_r(&req_pair_idx, &req_source, &req_target);
        temporal_scores_to_list(
            &results,
            &lr_pairs,
            &grid,
            &inter,
            min_pct as f32,
            min_n_eff as f32,
            return_complete_grid,
            requested_contexts,
        )
    }

    // ── 6e. ComputeTemporalPermutation ──────────────────────────────────────

    /// Temporal permutation null for trajectory-level LR dynamics.
    #[extendr]
    pub fn rs_compute_temporal_permutation(
        expr: S4,
        cell_types: Integers,
        n_types: i32,
        ligand_rows_list: List,
        receptor_rows_list: List,
        pseudotime: Doubles,
        lineage_weights: Doubles,
        block_labels: Integers,
        n_lineages: i32,
        n_bins: i32,
        min_pct: f64,
        min_n_eff: f64,
        min_valid_bins: i32,
        min_valid_fraction: f64,
        n_perms: i32,
        seed: f64,
        stat_code: i32,
        top_k: i32,
        permutation_stratum: i32,
        req_pair_idx: Integers,
        req_source: Integers,
        req_target: Integers,
        req_lineage: Integers,
    ) -> List {
        use crate::ccc::temporal_null::{
            compute_temporal_permutation_with_scale, PermutationStratumMode,
        };

        let sparse = r_dgcmatrix_to_sparse(&expr);
        let ct: Vec<usize> = cell_types
            .iter()
            .map(|v| (v.inner() - 1) as usize)
            .collect();
        let nt = n_types as usize;
        let lig = list_of_int_to_rows0(&ligand_rows_list);
        let rec = list_of_int_to_rows0(&receptor_rows_list);
        let lr_pairs: Vec<LrPair> = lig
            .into_iter()
            .zip(rec.into_iter())
            .map(|(l, r)| LrPair {
                ligand_rows: l,
                receptor_rows: r,
            })
            .collect();
        let pt: Vec<f32> = pseudotime.iter().map(|v| v.inner() as f32).collect();
        let lw: Vec<f32> = lineage_weights.iter().map(|v| v.inner() as f32).collect();
        let blocks: Vec<usize> = block_labels
            .iter()
            .map(|v| (v.inner() - 1) as usize)
            .collect();
        let config_hash = calibration_config_hash(
            &sparse,
            &ct,
            nt,
            &lr_pairs,
            &pt,
            &lw,
            n_lineages as usize,
            n_bins as usize,
            min_pct as f32,
            min_n_eff as f32,
            min_valid_bins.max(0) as usize,
            min_valid_fraction as f32,
            stat_code.max(0) as usize,
            top_k.max(1) as usize,
        );

        let stratum_mode = if permutation_stratum == 1 {
            PermutationStratumMode::BlockCelltype
        } else {
            PermutationStratumMode::BlockCelltypeDominantLineage
        };
        let requested_keys =
            requested_trajectory_keys_from_r(&req_pair_idx, &req_source, &req_target, &req_lineage);
        let (traj, bins, hill_s) = compute_temporal_permutation_with_scale(
            &sparse,
            &ct,
            nt,
            &lr_pairs,
            &pt,
            &lw,
            n_lineages as usize,
            n_bins as usize,
            min_pct as f32,
            min_n_eff as f32,
            min_valid_bins.max(0) as usize,
            min_valid_fraction as f32,
            n_perms.max(0) as usize,
            seed as u64,
            stat_code.max(0) as usize,
            top_k.max(1) as usize,
            &blocks,
            stratum_mode,
            requested_keys,
        );

        let tr_pair: Vec<i32> = traj.iter().map(|r| r.pair_idx as i32 + 1).collect();
        let tr_source: Vec<i32> = traj.iter().map(|r| r.source as i32 + 1).collect();
        let tr_target: Vec<i32> = traj.iter().map(|r| r.target as i32 + 1).collect();
        let tr_lineage: Vec<i32> = traj.iter().map(|r| r.lineage as i32 + 1).collect();
        let tr_stat: Vec<f64> = traj.iter().map(|r| r.statistic as f64).collect();
        let tr_p: Vec<f64> = traj.iter().map(|r| r.pvalue).collect();
        let tr_n_eff: Vec<f64> = traj.iter().map(|r| r.n_eff_cells as f64).collect();
        let tr_entropy: Vec<f64> = traj.iter().map(|r| r.lineage_entropy as f64).collect();
        let tr_n_valid: Vec<i32> = traj.iter().map(|r| r.n_valid_bins as i32).collect();
        let tr_valid_frac: Vec<f64> = traj.iter().map(|r| r.valid_bin_fraction as f64).collect();
        let tr_n_expr: Vec<i32> = traj
            .iter()
            .map(|r| r.n_expression_supported_bins as i32)
            .collect();
        let tr_n_eff_max: Vec<f64> = traj.iter().map(|r| r.n_eff_at_max as f64).collect();
        let tr_n_eff_min: Vec<f64> = traj.iter().map(|r| r.n_eff_at_min as f64).collect();
        let tr_valid_max: Vec<i32> = traj.iter().map(|r| r.valid_at_max as i32).collect();
        let tr_valid_min: Vec<i32> = traj.iter().map(|r| r.valid_at_min as i32).collect();
        let tr_evaluable: Vec<i32> = traj.iter().map(|r| r.trajectory_evaluable as i32).collect();
        let tr_perm_expected: Vec<i32> = traj.iter().map(|r| r.n_perm_expected as i32).collect();
        let tr_perm_observed: Vec<i32> = traj.iter().map(|r| r.n_perm_observed as i32).collect();
        let tr_perm_coverage: Vec<f64> = traj.iter().map(|r| r.perm_coverage).collect();
        let tr_null_median: Vec<f64> = traj.iter().map(|r| r.null_stat_median as f64).collect();
        let trajectories = List::from_names_and_values(
            [
                "pair_idx",
                "source",
                "target",
                "lineage",
                "temporal_stat",
                "p_temporal",
                "n_eff_cells",
                "lineage_entropy",
                "n_valid_bins",
                "valid_bin_fraction",
                "n_expression_supported_bins",
                "n_eff_at_max",
                "n_eff_at_min",
                "valid_at_max",
                "valid_at_min",
                "trajectory_evaluable",
                "n_perm_expected",
                "n_perm_observed",
                "perm_coverage",
                "temporal_null_median",
            ],
            [
                Integers::from_values(tr_pair).into_robj(),
                Integers::from_values(tr_source).into_robj(),
                Integers::from_values(tr_target).into_robj(),
                Integers::from_values(tr_lineage).into_robj(),
                Doubles::from_values(tr_stat).into_robj(),
                Doubles::from_values(tr_p).into_robj(),
                Doubles::from_values(tr_n_eff).into_robj(),
                Doubles::from_values(tr_entropy).into_robj(),
                Integers::from_values(tr_n_valid).into_robj(),
                Doubles::from_values(tr_valid_frac).into_robj(),
                Integers::from_values(tr_n_expr).into_robj(),
                Doubles::from_values(tr_n_eff_max).into_robj(),
                Doubles::from_values(tr_n_eff_min).into_robj(),
                Integers::from_values(tr_valid_max).into_robj(),
                Integers::from_values(tr_valid_min).into_robj(),
                Integers::from_values(tr_evaluable).into_robj(),
                Integers::from_values(tr_perm_expected).into_robj(),
                Integers::from_values(tr_perm_observed).into_robj(),
                Doubles::from_values(tr_perm_coverage).into_robj(),
                Doubles::from_values(tr_null_median).into_robj(),
            ],
        )
        .unwrap();

        let b_pair: Vec<i32> = bins.iter().map(|r| r.pair_idx as i32 + 1).collect();
        let b_source: Vec<i32> = bins.iter().map(|r| r.source as i32 + 1).collect();
        let b_target: Vec<i32> = bins.iter().map(|r| r.target as i32 + 1).collect();
        let b_lineage: Vec<i32> = bins.iter().map(|r| r.lineage as i32 + 1).collect();
        let b_bin: Vec<i32> = bins.iter().map(|r| r.bin as i32 + 1).collect();
        let b_score: Vec<f64> = bins.iter().map(|r| r.score as f64).collect();
        let b_dev: Vec<f64> = bins.iter().map(|r| r.deviation as f64).collect();
        let b_p: Vec<f64> = bins.iter().map(|r| r.pvalue).collect();
        let b_struct: Vec<i32> = bins.iter().map(|r| r.structural_valid as i32).collect();
        let b_expr: Vec<i32> = bins.iter().map(|r| r.expression_supported as i32).collect();
        let bin_table = List::from_names_and_values(
            [
                "pair_idx",
                "source",
                "target",
                "lineage",
                "bin",
                "score",
                "temporal_deviation",
                "p_temporal_bin",
                "structural_valid",
                "expression_supported",
            ],
            [
                Integers::from_values(b_pair).into_robj(),
                Integers::from_values(b_source).into_robj(),
                Integers::from_values(b_target).into_robj(),
                Integers::from_values(b_lineage).into_robj(),
                Integers::from_values(b_bin).into_robj(),
                Doubles::from_values(b_score).into_robj(),
                Doubles::from_values(b_dev).into_robj(),
                Doubles::from_values(b_p).into_robj(),
                Integers::from_values(b_struct).into_robj(),
                Integers::from_values(b_expr).into_robj(),
            ],
        )
        .unwrap();

        List::from_names_and_values(
            [
                "trajectories",
                "bins",
                "calibration_config_hash",
                "hill_scale",
            ],
            [
                trajectories.into_robj(),
                bin_table.into_robj(),
                Strings::from_values([config_hash]).into_robj(),
                Doubles::from_values([hill_s as f64]).into_robj(),
            ],
        )
        .unwrap()
    }

    // ── 6c. ComputeLineageEntropyRaw ───────────────────────────────────────

    /// Build the raw per-lineage entropy tensor `H_l(b)` using
    /// `H_l(b) = Σ_ct hwsums[ct,l,b] / Σ_ct wsums[ct,l,b]`.
    ///
    /// Bins where the denominator is zero carry NaN so the C2
    /// check can count invalid bins with `is.finite()`.
    ///
    /// @param hwsums numeric vector length n_types*n_lineages*n_bins,
    ///   row-major with ct outermost — produced by
    ///   `rs_compute_temporal_lr_with_inter`.
    /// @param wsums same shape as `hwsums`.
    /// @param n_types, n_lineages, n_bins integer dims.
    /// @return numeric vector length `n_lineages * n_bins`, row-major
    ///   (lineage outer, bin inner).  R wrapper reshapes to matrix.
    #[extendr]
    pub fn rs_compute_lineage_entropy_raw(
        hwsums: Doubles,
        wsums: Doubles,
        n_types: i32,
        n_lineages: i32,
        n_bins: i32,
    ) -> Doubles {
        use crate::ccc::{compute_lineage_entropy_raw, LrBinIntermediate};
        let inter = LrBinIntermediate {
            hwsums: hwsums.iter().map(|v| v.inner()).collect(),
            wsums: wsums.iter().map(|v| v.inner()).collect(),
            n_types: n_types as usize,
            n_lineages: n_lineages as usize,
            n_bins: n_bins as usize,
        };
        let out = compute_lineage_entropy_raw(&inter);
        Doubles::from_values(out.data.into_iter().map(|v| v as f64))
    }

    // ── 6d. SmoothMA3 ──────────────────────────────────────────────────────

    /// Apply the MA(3) smoother: 2-bin partial window at b=0 and
    /// b=n_bins-1, 3-bin window for interior b.  NaN propagates through
    /// any window it enters; the C1 bypass (n_bins < 6) returns
    /// the input unchanged.
    ///
    /// @param h_raw numeric vector length `n_lineages * n_bins`, row-major
    ///   (output of `rs_compute_lineage_entropy_raw`).
    /// @param n_lineages, n_bins integer dims.
    /// @return numeric vector same shape as `h_raw`, smoothed.
    #[extendr]
    pub fn rs_smooth_ma3(h_raw: Doubles, n_lineages: i32, n_bins: i32) -> Doubles {
        use crate::ccc::{smooth_ma3, LineageEntropyTensor};
        let raw = LineageEntropyTensor {
            data: h_raw.iter().map(|v| v.inner() as f32).collect(),
            n_lineages: n_lineages as usize,
            n_bins: n_bins as usize,
        };
        let out = smooth_ma3(&raw);
        Doubles::from_values(out.data.into_iter().map(|v| v as f64))
    }

    // ── 6e. FindBifurcationBins ────────────────────────────────────────────

    /// Locate the bifurcation bin `b*(l)` per lineage using the
    /// continuous-2-bin strict-less-than rule: first `b ∈ [0, n_bins-3]`
    /// whose own value AND `b+1` value both lie strictly below `h_fork`
    /// (default 0.5, a scale-natural choice).
    ///
    /// Edge-case tagging is folded into the returned `reason`:
    /// - `"crossing"`           — normal case, `b_star = 1-based b + 1`.
    /// - `"c2_not_informative"` — curve has < 4 finite bins OR
    ///   finite max - min < 0.05 (C2). `b_star = NA`.
    /// - `"c3_no_low_entropy"`  — informative curve that never satisfies
    ///   the continuous-2-bin rule (C3). `b_star = NA`.
    /// - `"c4_boundary_onset"`  — crossing happened at the first two
    ///   bins (`b ≤ 1`, 0-based), pre-fork taxonomy terms must be
    ///   gated off downstream (C4). `b_star = 1-based b + 1`.
    ///
    /// @param h_tilde numeric vector length `n_lineages * n_bins`,
    ///   row-major (typically output of `rs_smooth_ma3`).
    /// @param n_lineages, n_bins integer dims.
    /// @param h_fork numeric threshold (default 0.5).
    /// @return list with named elements:
    ///   - `b_star`: numeric vector of length `n_lineages`; NaN values
    ///     are `NA` on the R side after `as.integer(...)`.  1-based.
    ///   - `reason`: character vector of length `n_lineages`.
    #[extendr]
    pub fn rs_find_bifurcation_bins(
        h_tilde: Doubles,
        n_lineages: i32,
        n_bins: i32,
        h_fork: f64,
    ) -> List {
        use crate::ccc::{find_bifurcation_bins, BifurcationReason, LineageEntropyTensor};
        let tensor = LineageEntropyTensor {
            data: h_tilde.iter().map(|v| v.inner() as f32).collect(),
            n_lineages: n_lineages as usize,
            n_bins: n_bins as usize,
        };
        let results = find_bifurcation_bins(&tensor, h_fork as f32);

        // b_star carried as f64 with NaN ↔ None so the R side collapses
        // into NA via `as.integer(suppressWarnings(...))` without
        // needing extendr's (comparatively fussier) Rint-NA path.
        // 1-based (+ 1) per the LineaLinkR R-side index convention.
        let b_star: Vec<f64> = results
            .iter()
            .map(|r| match r.b_star {
                Some(b) => (b + 1) as f64,
                None => f64::NAN,
            })
            .collect();

        let reason: Vec<String> = results
            .iter()
            .map(|r| {
                match r.reason {
                    BifurcationReason::Crossing => "crossing",
                    BifurcationReason::C2NotInformative => "c2_not_informative",
                    BifurcationReason::C3NoLowEntropy => "c3_no_low_entropy",
                    BifurcationReason::C4BoundaryOnset => "c4_boundary_onset",
                }
                .to_string()
            })
            .collect();

        List::from_names_and_values(
            ["b_star", "reason"],
            [
                Doubles::from_values(b_star).into_robj(),
                Strings::from_values(reason).into_robj(),
            ],
        )
        .unwrap()
    }

    /// Medians of consecutive, finite bootstrap samples already drawn in R.
    /// This kernel consumes no random numbers and keeps the sample order.
    #[extendr]
    pub fn rs_bootstrap_medians(draws: Doubles, n_samples: i32) -> Result<Doubles> {
        let mut values: Vec<f64> = draws.iter().map(|x| x.inner()).collect();
        if n_samples < 1 || values.len() % n_samples as usize != 0 {
            return Err(Error::Other("Invalid bootstrap draw dimensions".into()));
        }
        let width = n_samples as usize;
        let middle = width / 2;
        let medians: Vec<f64> = values.chunks_exact_mut(width).map(|column| {
            let (lower, upper, _) = column.select_nth_unstable_by(middle, f64::total_cmp);
            let upper = *upper;
            if width % 2 == 1 {
                upper
            } else {
                let lower = lower.iter().copied().fold(f64::NEG_INFINITY, f64::max);
                let sum = lower + upper;
                if sum.is_finite() {
                    sum * 0.5
                } else {
                    lower * 0.5 + upper * 0.5
                }
            }
        }).collect();
        Ok(Doubles::from_values(medians))
    }

    // ── 6f. ComputeRoleScores ──────────────────────────────────────────────

    /// Compute the five role scores (instructive / reinforcing / feedback /
    /// lateral / bystander) per row. Batched over all
    /// `LrBinResult` rows; each row's own `bin` and `lineage` address
    /// into the per-lineage `BifurcationResult`.
    ///
    /// **Row-level inputs** (all length `n_rows`, same order as the
    /// `scores` data.frame produced by `ComputeTemporalLRWithIntermediate`):
    /// - `row_bin`: 1-based bin index for this row (Rust side decrements).
    /// - `row_lineage`: 1-based lineage index (Rust side decrements).
    /// - `spec_self`: specificity `spec(p, s, t, b, l)` on this row's
    ///   own lineage. Caller must supply it.
    /// - `spec_max_other`: `max_{l' != l} spec(p, s, t, b, l')` — used
    ///   for the `lateral` term only.  Set to 0 when `n_lineages == 1`.
    /// - `delta_tau`: Δτ = source_pt - target_pt. NaN
    ///   propagates through sigmoid and flags the row.
    ///
    /// **Per-lineage inputs** (all length `n_lineages`, matching
    /// `FindBifurcationBins` output row order):
    /// - `lineage_b_star`: 1-based `b*(l)` or `NaN` when C2/C3
    ///   downgraded the lineage to time-only.
    /// - `lineage_reason_code`: 0 = crossing, 1 = c2_not_informative,
    ///   2 = c3_no_low_entropy, 3 = c4_boundary_onset. The C4
    ///   pre-fork gate activates on code 3 even with a valid `b_star`.
    ///
    /// **row-level reliability gate**:
    /// - `row_n_eff_cells`: per-row `n_eff_cells = min(source, target)`
    ///   Kish effective sample size (typically read from the `n_eff_cells`
    ///   column of `ComputeTemporalLRWithIntermediate`'s `scores` df).
    /// - `min_n_eff_cells`: scalar threshold; rows with strictly less
    ///   are forced to `(0,0,0,0,1)` regardless of formula inputs.
    ///   Pass `0.0` to disable the gate.
    ///
    /// @param tau_0 numeric, direction-sigmoid scale (typically
    ///   10% of per-lineage pseudotime std; caller-computed).
    /// @param h_fork numeric threshold retained in the call metadata.
    /// @return list with 6 numeric vectors each of length `n_rows`:
    ///   `instructive`, `transition_like`, `reinforcing`, `feedback`,
    ///   `lateral`, `bystander`.
    /// @param min_n_eff_cells numeric row-level reliability gate (default 5.0).
    /// @param transition_window non-negative integer; bins with
    ///   `abs(bin - b_star) <= transition_window` receive non-zero
    ///   `transition_like` (default 0 = exact bifurcation bin only).
    #[extendr]
    pub fn rs_compute_role_scores(
        row_bin: Integers,
        row_lineage: Integers,
        spec_self: Doubles,
        spec_max_other: Doubles,
        delta_tau: Doubles,
        row_n_eff_cells: Doubles,
        lineage_b_star: Doubles,
        lineage_reason_code: Integers,
        tau_0: f64,
        h_fork: f64,
        min_n_eff_cells: f64,
        transition_window: i32,
    ) -> List {
        use crate::ccc::{
            compute_role_scores_one, BifurcationReason, BifurcationResult, LrBinResult, Role,
            RoleScoreParams,
        };
        let _ = Role::Instructive; // silence unused-import pattern

        let n_rows = row_bin.len();
        let n_lineages = lineage_b_star.len();
        assert_eq!(row_lineage.len(), n_rows, "row_lineage length mismatch");
        assert_eq!(spec_self.len(), n_rows, "spec_self length mismatch");
        assert_eq!(
            spec_max_other.len(),
            n_rows,
            "spec_max_other length mismatch"
        );
        assert_eq!(delta_tau.len(), n_rows, "delta_tau length mismatch");
        assert_eq!(
            row_n_eff_cells.len(),
            n_rows,
            "row_n_eff_cells length mismatch"
        );
        assert_eq!(
            lineage_reason_code.len(),
            n_lineages,
            "lineage_reason_code length mismatch with lineage_b_star"
        );

        // ── Pre-parse per-lineage BifurcationResult vector ─────────────
        // Matches FindBifurcationBins output: NaN b_star = None, else
        // 1-based → 0-based; reason_code 0..3 maps to the 4 enum variants.
        let mut bif: Vec<BifurcationResult> = Vec::with_capacity(n_lineages);
        for i in 0..n_lineages {
            let b_raw = lineage_b_star[i].inner();
            let b_opt = if b_raw.is_nan() {
                None
            } else {
                Some((b_raw as usize).saturating_sub(1))
            };
            let code = lineage_reason_code[i].inner();
            let reason = match code {
                0 => BifurcationReason::Crossing,
                1 => BifurcationReason::C2NotInformative,
                2 => BifurcationReason::C3NoLowEntropy,
                3 => BifurcationReason::C4BoundaryOnset,
                _ => BifurcationReason::C2NotInformative, // default for unknown codes
            };
            bif.push(BifurcationResult {
                b_star: b_opt,
                reason,
            });
        }

        let params = RoleScoreParams {
            tau_0: tau_0 as f32,
            h_fork: h_fork as f32,
            min_n_eff_cells: min_n_eff_cells as f32,
            transition_window: transition_window.max(0) as u32,
        };

        // ── Row-wise score computation ─────────────────────────────────
        // LrBinResult is built as a stub because compute_role_scores_one only
        // reads `.bin`; the full struct is still required by the function
        // signature. The lineage field is populated in the returned table.
        let mut instructive = Vec::with_capacity(n_rows);
        let mut transition_like = Vec::with_capacity(n_rows);
        let mut reinforcing = Vec::with_capacity(n_rows);
        let mut feedback = Vec::with_capacity(n_rows);
        let mut lateral = Vec::with_capacity(n_rows);
        let mut bystander = Vec::with_capacity(n_rows);

        for i in 0..n_rows {
            // 1-based → 0-based for both indices.
            let bin_idx = (row_bin[i].inner() as usize).saturating_sub(1);
            let lineage_idx = (row_lineage[i].inner() as usize).saturating_sub(1);
            assert!(
                lineage_idx < n_lineages,
                "row {i}: lineage index {lineage_idx} out of range (n_lineages = {n_lineages})"
            );

            let stub = LrBinResult {
                pair_idx: 0,
                source: 0,
                target: 0,
                lineage: lineage_idx,
                bin: bin_idx,
                score: 0.0,
                // gate reads this; populated from caller-supplied
                // per-row `n_eff_cells` (typically scores_df$n_eff_cells).
                n_eff_cells: row_n_eff_cells[i].inner() as f32,
                weight_mass_cells: row_n_eff_cells[i].inner() as f32,
                lineage_entropy: 0.0,
            };
            let scores = compute_role_scores_one(
                &stub,
                spec_self[i].inner() as f32,
                spec_max_other[i].inner() as f32,
                delta_tau[i].inner() as f32,
                bif[lineage_idx],
                &params,
            );
            instructive.push(scores.instructive as f64);
            transition_like.push(scores.transition_like as f64);
            reinforcing.push(scores.reinforcing as f64);
            feedback.push(scores.feedback as f64);
            lateral.push(scores.lateral as f64);
            bystander.push(scores.bystander as f64);
        }

        List::from_names_and_values(
            [
                "instructive",
                "transition_like",
                "reinforcing",
                "feedback",
                "lateral",
                "bystander",
            ],
            [
                Doubles::from_values(instructive).into_robj(),
                Doubles::from_values(transition_like).into_robj(),
                Doubles::from_values(reinforcing).into_robj(),
                Doubles::from_values(feedback).into_robj(),
                Doubles::from_values(lateral).into_robj(),
                Doubles::from_values(bystander).into_robj(),
            ],
        )
        .unwrap()
    }

    // ── 7. TrajectoryDE ─────────────────────────────────────────────────────

    /// Per-gene Spearman correlation of expression vs. pseudotime, BH-corrected.
    /// Inputs/outputs are 0-indexed and dimensionless; no R-side index
    /// translation is needed because there are no index references.
    ///
    /// @param expr_flat numeric vector, row-major flat \[n_genes × n_cells\]
    /// @param pseudotime numeric vector of length n_cells; NaN cells excluded
    /// @param n_genes, n_cells, min_cells integers
    /// @return list with `rho`, `pval`, `qval` (each length n_genes); rows
    ///   with < min_cells qualifying cells get rho=0, pval=1, qval=1.
    #[extendr]
    pub fn rs_trajectory_de(
        expr_flat: Doubles,
        pseudotime: Doubles,
        n_genes: i32,
        n_cells: i32,
        min_cells: i32,
    ) -> List {
        let ef: Vec<f32> = expr_flat.iter().map(|v| v.inner() as f32).collect();
        let pt: Vec<f32> = pseudotime.iter().map(|v| v.inner() as f32).collect();
        let (rho, pval, qval) = pseudotime_de(
            &ef,
            &pt,
            n_genes as usize,
            n_cells as usize,
            min_cells as usize,
        );
        List::from_names_and_values(
            ["rho", "pval", "qval"],
            [
                Doubles::from_values(rho).into_robj(),
                Doubles::from_values(pval).into_robj(),
                Doubles::from_values(qval).into_robj(),
            ],
        )
        .unwrap()
    }

    // ── 8. ComputeFakeLRNull ───────────────────────────────────────────────

    ///  Fake LR null model: stratified gene-pair swap matching + scoring.
    ///
    /// Constructs fake ligand-receptor pairs from a candidate gene pool
    /// (genes NOT in any LR database), matched to real LR genes on marginal
    /// expression features, then scores them and computes per-row p_fake.
    ///
    /// @param expr dgCMatrix (genes × cells), log-normalised.
    /// @param cell_types integer vector, 1-based.
    /// @param n_types integer.
    /// @param candidate_fake_genes integer vector, 1-based gene row indices
    ///   eligible for fake pairing (must be pre-filtered in R to exclude
    ///   all known LR databases + STRING/BioGRID).
    /// @param real_ligand_rows_list list of int vecs (1-based), per-LR-pair.
    /// @param real_receptor_rows_list same.
    /// @param pseudotime numeric vector.
    /// @param lineage_weights numeric vector, flat row-major [n_cells × n_lineages].
    /// @param n_lineages, n_bins integers.
    /// @param min_pct numeric.
    /// @param n_fake_per_real integer, target |D*| / |D| multiplier (default 10).
    /// @param n_strata integer, per-dimension equal-width strata count (default 5).
    /// @param caliper_k integer, top-k nearest in stratum (default 10).
    /// @param seed numeric.
    /// @param return_bin_scores logical; skip the large fake bin table when false.
    /// @param stat_code integer trajectory statistic code.
    /// @param top_k integer used by top-k mean trajectory statistic.
    /// @return list with: fake_scores (data.frame columns), p_fake (numeric vector
    ///   parallel to caller's real scores), fake_pairs (data.frame of
    ///   ligand_gene/receptor_gene/shadowing_pair_idx), candidate_pool_size (int).
    #[extendr]
    pub fn rs_compute_fake_lr_null(
        expr: S4,
        cell_types: Integers,
        n_types: i32,
        candidate_fake_genes: Integers,
        real_ligand_rows_list: List,
        real_receptor_rows_list: List,
        pseudotime: Doubles,
        lineage_weights: Doubles,
        n_lineages: i32,
        n_bins: i32,
        min_pct: f64,
        n_fake_per_real: i32,
        n_strata: i32,
        caliper_k: i32,
        max_match_distance: f64,
        seed: f64,
        return_bin_scores: bool,
        stat_code: i32,
        top_k: i32,
        min_n_eff: f64,
        min_valid_bins: i32,
        min_valid_fraction: f64,
        fake_null_mode: i32,
        complex_null_mode: i32,
        block_labels: Integers,
    ) -> List {
        use crate::ccc::lr_scoring::compute_lr_scores_bins_with_scale_grid;
        use crate::ccc::trajectory_curve::trajectory_keys_from_scores;
        use crate::ccc::{
            calibrate_fake_trajectories, compute_gene_bin_stats, compute_gene_features,
            construct_fake_pairs_context_aware, score_fake_pairs_bins, ComplexNullMode,
            FakeMatchingBlocks, FakeNullMode, LrPair,
        };

        let sparse = r_dgcmatrix_to_sparse(&expr);
        let ct: Vec<usize> = cell_types
            .iter()
            .map(|v| (v.inner() - 1) as usize)
            .collect();
        let nt = n_types as usize;
        let pt: Vec<f32> = pseudotime.iter().map(|v| v.inner() as f32).collect();
        let lw: Vec<f32> = lineage_weights.iter().map(|v| v.inner() as f32).collect();
        let nl = n_lineages as usize;
        let nb = n_bins as usize;
        let blocks: Vec<usize> = block_labels
            .iter()
            .map(|v| (v.inner() - 1) as usize)
            .collect();
        let n_matching_blocks = blocks.iter().max().map(|m| m + 1).unwrap_or(1usize);

        let lig = list_of_int_to_rows0(&real_ligand_rows_list);
        let rec = list_of_int_to_rows0(&real_receptor_rows_list);
        let n_pairs = lig.len();

        // ── Step 0: candidate pool construction and exclusion ───────────
        let mut real_gene_set = std::collections::HashSet::new();
        for rows in lig.iter().chain(rec.iter()) {
            real_gene_set.extend(rows.iter().copied());
        }
        let mut seen_fake = std::collections::HashSet::new();
        let fake_genes: Vec<usize> = candidate_fake_genes
            .iter()
            .map(|v| (v.inner() - 1) as usize)
            .filter(|gene| !real_gene_set.contains(gene))
            .filter(|gene| seen_fake.insert(*gene))
            .collect();
        let n_candidate = fake_genes.len();

        // ── Step 1: compute GeneBinStats for candidate fake genes ─────────
        let fake_stats =
            compute_gene_bin_stats(&sparse, &ct, nt, &fake_genes, &pt, &lw, nl, nb, None, None);

        // ── Step 2: compute features for candidate fake genes ─────────────
        let fake_features = compute_gene_features(&fake_stats, &fake_stats.gene_rows);

        // ── Step 3: compute features for real LR genes (unique set) ───────
        // Collect all unique gene indices in the real LR database.
        let mut real_gene_list: Vec<usize> = real_gene_set.into_iter().collect();
        real_gene_list.sort_unstable();

        let real_stats = compute_gene_bin_stats(
            &sparse,
            &ct,
            nt,
            &real_gene_list,
            &pt,
            &lw,
            nl,
            nb,
            None,
            None,
        );
        let real_features = compute_gene_features(&real_stats, &real_gene_list);

        let block_real_stats: Vec<_> = if n_matching_blocks > 1 {
            (0..n_matching_blocks)
                .map(|b| {
                    compute_gene_bin_stats(
                        &sparse,
                        &ct,
                        nt,
                        &real_gene_list,
                        &pt,
                        &lw,
                        nl,
                        nb,
                        Some(&blocks),
                        Some(b),
                    )
                })
                .collect()
        } else {
            Vec::new()
        };
        let block_fake_stats: Vec<_> = if n_matching_blocks > 1 {
            (0..n_matching_blocks)
                .map(|b| {
                    compute_gene_bin_stats(
                        &sparse,
                        &ct,
                        nt,
                        &fake_genes,
                        &pt,
                        &lw,
                        nl,
                        nb,
                        Some(&blocks),
                        Some(b),
                    )
                })
                .collect()
        } else {
            Vec::new()
        };

        // ── Step 4: build LR pairs and real scores (trajectory keys) ─────
        let lr_pairs: Vec<LrPair> = lig
            .iter()
            .zip(rec.iter())
            .map(|(l, r)| LrPair {
                ligand_rows: l.clone(),
                receptor_rows: r.clone(),
            })
            .collect();

        let (real_scores, _, hill_s, real_grid) = compute_lr_scores_bins_with_scale_grid(
            &sparse,
            &ct,
            nt,
            &lr_pairs,
            &pt,
            &lw,
            nl,
            nb,
            min_pct as f32,
        );
        let traj_keys = trajectory_keys_from_scores(&real_scores);
        let n_traj_keys = traj_keys.len();

        let mut tau_by_gene: std::collections::HashMap<usize, f32> =
            std::collections::HashMap::new();
        for f in fake_features.iter() {
            tau_by_gene.insert(f.gene_idx, f.tau);
        }
        for (i, g) in real_gene_list.iter().enumerate() {
            tau_by_gene.entry(*g).or_insert(real_features[i].tau);
        }

        // ── Step 5: context-aware fake pair construction ───────────────
        let n_fpr = n_fake_per_real.max(1) as usize;
        let ns = n_strata.max(1) as usize;
        let ck = caliper_k.max(1) as usize;
        let null_mode = match fake_null_mode {
            1 => FakeNullMode::EdgeRewire,
            2 => FakeNullMode::Both,
            3 => FakeNullMode::LrNetworkRewire,
            4 => FakeNullMode::All,
            _ => FakeNullMode::ContextGeneric,
        };
        let complex_mode = if complex_null_mode == 1 {
            ComplexNullMode::CardinalityMatched
        } else {
            ComplexNullMode::Skip
        };
        let block_matching = if n_matching_blocks > 1 {
            Some(FakeMatchingBlocks {
                block_labels: &blocks,
                cell_types: &ct,
                lineage_weights: &lw,
                n_lineages: nl,
                n_blocks: n_matching_blocks,
                real_stats: &block_real_stats,
                fake_stats: &block_fake_stats,
            })
        } else {
            None
        };
        let fake_pairs = construct_fake_pairs_context_aware(
            &traj_keys,
            &lr_pairs,
            &real_stats,
            &fake_stats,
            &fake_genes,
            &tau_by_gene,
            n_fpr,
            ns,
            ck,
            max_match_distance as f32,
            null_mode,
            complex_mode,
            seed as u64,
            block_matching.as_ref(),
        );

        let config_hash = calibration_config_hash(
            &sparse,
            &ct,
            nt,
            &lr_pairs,
            &pt,
            &lw,
            nl,
            nb,
            min_pct as f32,
            min_n_eff as f32,
            min_valid_bins.max(0) as usize,
            min_valid_fraction as f32,
            stat_code.max(0) as usize,
            top_k.max(1) as usize,
        );

        // ── Step 6: score fake pairs on the observed real-score scale ─────
        // A shared Hill S is required for direct real-vs-fake comparison.
        let fake_scores = if return_bin_scores {
            score_fake_pairs_bins(
                &fake_pairs,
                &real_stats,
                &fake_stats,
                hill_s,
                nt,
                nl,
                nb,
                min_pct as f32,
            )
        } else {
            Vec::new()
        };
        let trajectory_calibration = calibrate_fake_trajectories(
            &lr_pairs,
            &real_scores,
            &fake_pairs,
            &real_grid,
            &real_stats,
            &fake_stats,
            min_pct as f32,
            min_n_eff as f32,
            min_valid_bins.max(0) as usize,
            min_valid_fraction as f32,
            stat_code.max(0) as usize,
            top_k.max(1) as usize,
            n_fpr,
            null_mode,
        );

        // ── Step 7: prepare return values ──────────────────────────────────
        // fake_scores as data.frame columns (1-based like real scores)
        let f_fake_id: Vec<i32> = fake_scores.iter().map(|r| r.fake_idx as i32 + 1).collect();
        let f_pair: Vec<i32> = fake_scores.iter().map(|r| r.pair_idx as i32 + 1).collect();
        let f_source: Vec<i32> = fake_scores.iter().map(|r| r.source as i32 + 1).collect();
        let f_target: Vec<i32> = fake_scores.iter().map(|r| r.target as i32 + 1).collect();
        let f_lineage: Vec<i32> = fake_scores.iter().map(|r| r.lineage as i32 + 1).collect();
        let f_bin: Vec<i32> = fake_scores.iter().map(|r| r.bin as i32 + 1).collect();
        let f_score: Vec<f64> = fake_scores.iter().map(|r| r.score as f64).collect();
        let f_n_eff: Vec<f64> = fake_scores.iter().map(|r| r.n_eff_cells as f64).collect();
        let f_weight_mass: Vec<f64> = fake_scores
            .iter()
            .map(|r| r.weight_mass_cells as f64)
            .collect();

        let fake_scores_list = List::from_names_and_values(
            [
                "fake_id",
                "pair_idx",
                "source",
                "target",
                "lineage",
                "bin",
                "score",
                "n_eff_cells",
                "weight_mass_cells",
            ],
            [
                Integers::from_values(f_fake_id).into_robj(),
                Integers::from_values(f_pair).into_robj(),
                Integers::from_values(f_source).into_robj(),
                Integers::from_values(f_target).into_robj(),
                Integers::from_values(f_lineage).into_robj(),
                Integers::from_values(f_bin).into_robj(),
                Doubles::from_values(f_score).into_robj(),
                Doubles::from_values(f_n_eff).into_robj(),
                Doubles::from_values(f_weight_mass).into_robj(),
            ],
        )
        .unwrap();

        // fake_pairs info
        let fp_lig: Vec<i32> = fake_pairs
            .iter()
            .map(|p| p.ligand_gene as i32 + 1)
            .collect();
        let fp_rec: Vec<i32> = fake_pairs
            .iter()
            .map(|p| p.receptor_gene as i32 + 1)
            .collect();
        let fp_shadow: Vec<i32> = fake_pairs
            .iter()
            .map(|p| p.shadowing_pair_idx as i32 + 1)
            .collect();
        let fp_lig_distance: Vec<f64> = fake_pairs
            .iter()
            .map(|p| p.ligand_match_distance as f64)
            .collect();
        let fp_rec_distance: Vec<f64> = fake_pairs
            .iter()
            .map(|p| p.receptor_match_distance as f64)
            .collect();
        let fp_source: Vec<i32> = fake_pairs.iter().map(|p| p.source as i32 + 1).collect();
        let fp_target: Vec<i32> = fake_pairs.iter().map(|p| p.target as i32 + 1).collect();
        let fp_lineage: Vec<i32> = fake_pairs.iter().map(|p| p.lineage as i32 + 1).collect();
        let fp_null_kind: Vec<i32> = fake_pairs.iter().map(|p| p.null_kind as i32).collect();
        let fake_pairs_info = List::from_names_and_values(
            [
                "ligand_gene",
                "receptor_gene",
                "shadowing_pair",
                "source",
                "target",
                "lineage",
                "null_kind",
                "ligand_match_distance",
                "receptor_match_distance",
            ],
            [
                Integers::from_values(fp_lig).into_robj(),
                Integers::from_values(fp_rec).into_robj(),
                Integers::from_values(fp_shadow).into_robj(),
                Integers::from_values(fp_source).into_robj(),
                Integers::from_values(fp_target).into_robj(),
                Integers::from_values(fp_lineage).into_robj(),
                Integers::from_values(fp_null_kind).into_robj(),
                Doubles::from_values(fp_lig_distance).into_robj(),
                Doubles::from_values(fp_rec_distance).into_robj(),
            ],
        )
        .unwrap();

        let tc_pair: Vec<i32> = trajectory_calibration
            .iter()
            .map(|r| r.pair_idx as i32 + 1)
            .collect();
        let tc_source: Vec<i32> = trajectory_calibration
            .iter()
            .map(|r| r.source as i32 + 1)
            .collect();
        let tc_target: Vec<i32> = trajectory_calibration
            .iter()
            .map(|r| r.target as i32 + 1)
            .collect();
        let tc_lineage: Vec<i32> = trajectory_calibration
            .iter()
            .map(|r| r.lineage as i32 + 1)
            .collect();
        let tc_stat: Vec<f64> = trajectory_calibration
            .iter()
            .map(|r| r.statistic as f64)
            .collect();
        let tc_fake_null_median: Vec<f64> = trajectory_calibration
            .iter()
            .map(|r| r.fake_null_stat_median as f64)
            .collect();
        let tc_p: Vec<f64> = trajectory_calibration.iter().map(|r| r.pvalue).collect();
        let tc_p_generic: Vec<f64> = trajectory_calibration
            .iter()
            .map(|r| r.pvalue_generic)
            .collect();
        let tc_p_edge: Vec<f64> = trajectory_calibration
            .iter()
            .map(|r| r.pvalue_edge)
            .collect();
        let tc_p_replace_ligand: Vec<f64> = trajectory_calibration
            .iter()
            .map(|r| r.pvalue_replace_ligand)
            .collect();
        let tc_p_replace_receptor: Vec<f64> = trajectory_calibration
            .iter()
            .map(|r| r.pvalue_replace_receptor)
            .collect();
        let tc_p_network: Vec<f64> = trajectory_calibration
            .iter()
            .map(|r| r.pvalue_network)
            .collect();
        let tc_p_network_ligand: Vec<f64> = trajectory_calibration
            .iter()
            .map(|r| r.pvalue_network_rewire_ligand)
            .collect();
        let tc_p_network_receptor: Vec<f64> = trajectory_calibration
            .iter()
            .map(|r| r.pvalue_network_rewire_receptor)
            .collect();
        let tc_p_complex: Vec<f64> = trajectory_calibration
            .iter()
            .map(|r| r.pvalue_complex)
            .collect();
        let tc_mean_lig: Vec<f64> = trajectory_calibration
            .iter()
            .map(|r| r.mean_ligand_match_distance as f64)
            .collect();
        let tc_mean_rec: Vec<f64> = trajectory_calibration
            .iter()
            .map(|r| r.mean_receptor_match_distance as f64)
            .collect();
        let tc_unique_draws: Vec<i32> = trajectory_calibration
            .iter()
            .map(|r| r.n_unique_fake_draws as i32)
            .collect();
        let tc_unique_generic: Vec<i32> = trajectory_calibration
            .iter()
            .map(|r| r.n_unique_fake_draws_generic as i32)
            .collect();
        let tc_unique_edge: Vec<i32> = trajectory_calibration
            .iter()
            .map(|r| r.n_unique_fake_draws_edge as i32)
            .collect();
        let tc_exp_generic: Vec<i32> = trajectory_calibration
            .iter()
            .map(|r| r.n_fake_expected_generic as i32)
            .collect();
        let tc_exp_edge: Vec<i32> = trajectory_calibration
            .iter()
            .map(|r| r.n_fake_expected_edge as i32)
            .collect();
        let tc_obs_generic: Vec<i32> = trajectory_calibration
            .iter()
            .map(|r| r.n_fake_observed_generic as i32)
            .collect();
        let tc_obs_edge: Vec<i32> = trajectory_calibration
            .iter()
            .map(|r| r.n_fake_observed_edge as i32)
            .collect();
        let tc_cov_generic: Vec<f64> = trajectory_calibration
            .iter()
            .map(|r| r.fake_coverage_generic)
            .collect();
        let tc_cov_edge: Vec<f64> = trajectory_calibration
            .iter()
            .map(|r| r.fake_coverage_edge)
            .collect();
        let tc_exp_replace_ligand: Vec<i32> = trajectory_calibration
            .iter()
            .map(|r| r.n_fake_expected_replace_ligand as i32)
            .collect();
        let tc_exp_replace_receptor: Vec<i32> = trajectory_calibration
            .iter()
            .map(|r| r.n_fake_expected_replace_receptor as i32)
            .collect();
        let tc_obs_replace_ligand: Vec<i32> = trajectory_calibration
            .iter()
            .map(|r| r.n_fake_observed_replace_ligand as i32)
            .collect();
        let tc_obs_replace_receptor: Vec<i32> = trajectory_calibration
            .iter()
            .map(|r| r.n_fake_observed_replace_receptor as i32)
            .collect();
        let tc_cov_replace_ligand: Vec<f64> = trajectory_calibration
            .iter()
            .map(|r| r.fake_coverage_replace_ligand)
            .collect();
        let tc_cov_replace_receptor: Vec<f64> = trajectory_calibration
            .iter()
            .map(|r| r.fake_coverage_replace_receptor)
            .collect();
        let tc_unique_replace_ligand: Vec<i32> = trajectory_calibration
            .iter()
            .map(|r| r.n_unique_fake_draws_replace_ligand as i32)
            .collect();
        let tc_unique_replace_receptor: Vec<i32> = trajectory_calibration
            .iter()
            .map(|r| r.n_unique_fake_draws_replace_receptor as i32)
            .collect();
        let tc_exp_network: Vec<i32> = trajectory_calibration
            .iter()
            .map(|r| r.n_fake_expected_network as i32)
            .collect();
        let tc_obs_network: Vec<i32> = trajectory_calibration
            .iter()
            .map(|r| r.n_fake_observed_network as i32)
            .collect();
        let tc_cov_network: Vec<f64> = trajectory_calibration
            .iter()
            .map(|r| r.fake_coverage_network)
            .collect();
        let tc_unique_network: Vec<i32> = trajectory_calibration
            .iter()
            .map(|r| r.n_unique_fake_draws_network as i32)
            .collect();
        let tc_exp_network_ligand: Vec<i32> = trajectory_calibration
            .iter()
            .map(|r| r.n_fake_expected_network_rewire_ligand as i32)
            .collect();
        let tc_exp_network_receptor: Vec<i32> = trajectory_calibration
            .iter()
            .map(|r| r.n_fake_expected_network_rewire_receptor as i32)
            .collect();
        let tc_obs_network_ligand: Vec<i32> = trajectory_calibration
            .iter()
            .map(|r| r.n_fake_observed_network_rewire_ligand as i32)
            .collect();
        let tc_obs_network_receptor: Vec<i32> = trajectory_calibration
            .iter()
            .map(|r| r.n_fake_observed_network_rewire_receptor as i32)
            .collect();
        let tc_cov_network_ligand: Vec<f64> = trajectory_calibration
            .iter()
            .map(|r| r.fake_coverage_network_rewire_ligand)
            .collect();
        let tc_cov_network_receptor: Vec<f64> = trajectory_calibration
            .iter()
            .map(|r| r.fake_coverage_network_rewire_receptor)
            .collect();
        let tc_unique_network_ligand: Vec<i32> = trajectory_calibration
            .iter()
            .map(|r| r.n_unique_fake_draws_network_rewire_ligand as i32)
            .collect();
        let tc_unique_network_receptor: Vec<i32> = trajectory_calibration
            .iter()
            .map(|r| r.n_unique_fake_draws_network_rewire_receptor as i32)
            .collect();
        let tc_expected: Vec<i32> = trajectory_calibration
            .iter()
            .map(|r| r.n_fake_expected as i32)
            .collect();
        let tc_observed: Vec<i32> = trajectory_calibration
            .iter()
            .map(|r| r.n_fake_observed as i32)
            .collect();
        let tc_coverage: Vec<f64> = trajectory_calibration
            .iter()
            .map(|r| r.fake_coverage)
            .collect();
        let tc_n_valid: Vec<i32> = trajectory_calibration
            .iter()
            .map(|r| r.n_valid_bins as i32)
            .collect();
        let tc_valid_frac: Vec<f64> = trajectory_calibration
            .iter()
            .map(|r| r.valid_bin_fraction as f64)
            .collect();
        let tc_n_expr: Vec<i32> = trajectory_calibration
            .iter()
            .map(|r| r.n_expression_supported_bins as i32)
            .collect();
        let tc_n_eff_max: Vec<f64> = trajectory_calibration
            .iter()
            .map(|r| r.n_eff_at_max as f64)
            .collect();
        let tc_n_eff_min: Vec<f64> = trajectory_calibration
            .iter()
            .map(|r| r.n_eff_at_min as f64)
            .collect();
        let tc_valid_max: Vec<i32> = trajectory_calibration
            .iter()
            .map(|r| r.valid_at_max as i32)
            .collect();
        let tc_valid_min: Vec<i32> = trajectory_calibration
            .iter()
            .map(|r| r.valid_at_min as i32)
            .collect();
        let tc_evaluable: Vec<i32> = trajectory_calibration
            .iter()
            .map(|r| r.trajectory_evaluable as i32)
            .collect();
        let trajectory_calibration_list = List::from_names_and_values(
            [
                "pair_idx",
                "source",
                "target",
                "lineage",
                "trajectory_stat",
                "fake_null_median",
                "p_fake",
                "p_fake_generic",
                "p_fake_edge",
                "p_replace_ligand",
                "p_replace_receptor",
                "p_fake_network",
                "p_network_rewire_ligand",
                "p_network_rewire_receptor",
                "p_fake_complex",
                "n_fake_expected",
                "n_fake_observed",
                "fake_coverage",
                "mean_ligand_match_distance",
                "mean_receptor_match_distance",
                "n_unique_fake_draws",
                "n_unique_fake_draws_generic",
                "n_unique_fake_draws_edge",
                "n_fake_expected_generic",
                "n_fake_expected_edge",
                "n_fake_observed_generic",
                "n_fake_observed_edge",
                "fake_coverage_generic",
                "fake_coverage_edge",
                "n_fake_expected_replace_ligand",
                "n_fake_expected_replace_receptor",
                "n_fake_observed_replace_ligand",
                "n_fake_observed_replace_receptor",
                "fake_coverage_replace_ligand",
                "fake_coverage_replace_receptor",
                "n_unique_fake_draws_replace_ligand",
                "n_unique_fake_draws_replace_receptor",
                "n_fake_expected_network",
                "n_fake_observed_network",
                "fake_coverage_network",
                "n_unique_fake_draws_network",
                "n_fake_expected_network_rewire_ligand",
                "n_fake_expected_network_rewire_receptor",
                "n_fake_observed_network_rewire_ligand",
                "n_fake_observed_network_rewire_receptor",
                "fake_coverage_network_rewire_ligand",
                "fake_coverage_network_rewire_receptor",
                "n_unique_fake_draws_network_rewire_ligand",
                "n_unique_fake_draws_network_rewire_receptor",
                "n_valid_bins",
                "valid_bin_fraction",
                "n_expression_supported_bins",
                "n_eff_at_max",
                "n_eff_at_min",
                "valid_at_max",
                "valid_at_min",
                "trajectory_evaluable",
            ],
            [
                Integers::from_values(tc_pair).into_robj(),
                Integers::from_values(tc_source).into_robj(),
                Integers::from_values(tc_target).into_robj(),
                Integers::from_values(tc_lineage).into_robj(),
                Doubles::from_values(tc_stat).into_robj(),
                Doubles::from_values(tc_fake_null_median).into_robj(),
                Doubles::from_values(tc_p).into_robj(),
                Doubles::from_values(tc_p_generic).into_robj(),
                Doubles::from_values(tc_p_edge).into_robj(),
                Doubles::from_values(tc_p_replace_ligand).into_robj(),
                Doubles::from_values(tc_p_replace_receptor).into_robj(),
                Doubles::from_values(tc_p_network).into_robj(),
                Doubles::from_values(tc_p_network_ligand).into_robj(),
                Doubles::from_values(tc_p_network_receptor).into_robj(),
                Doubles::from_values(tc_p_complex).into_robj(),
                Integers::from_values(tc_expected).into_robj(),
                Integers::from_values(tc_observed).into_robj(),
                Doubles::from_values(tc_coverage).into_robj(),
                Doubles::from_values(tc_mean_lig).into_robj(),
                Doubles::from_values(tc_mean_rec).into_robj(),
                Integers::from_values(tc_unique_draws).into_robj(),
                Integers::from_values(tc_unique_generic).into_robj(),
                Integers::from_values(tc_unique_edge).into_robj(),
                Integers::from_values(tc_exp_generic).into_robj(),
                Integers::from_values(tc_exp_edge).into_robj(),
                Integers::from_values(tc_obs_generic).into_robj(),
                Integers::from_values(tc_obs_edge).into_robj(),
                Doubles::from_values(tc_cov_generic).into_robj(),
                Doubles::from_values(tc_cov_edge).into_robj(),
                Integers::from_values(tc_exp_replace_ligand).into_robj(),
                Integers::from_values(tc_exp_replace_receptor).into_robj(),
                Integers::from_values(tc_obs_replace_ligand).into_robj(),
                Integers::from_values(tc_obs_replace_receptor).into_robj(),
                Doubles::from_values(tc_cov_replace_ligand).into_robj(),
                Doubles::from_values(tc_cov_replace_receptor).into_robj(),
                Integers::from_values(tc_unique_replace_ligand).into_robj(),
                Integers::from_values(tc_unique_replace_receptor).into_robj(),
                Integers::from_values(tc_exp_network).into_robj(),
                Integers::from_values(tc_obs_network).into_robj(),
                Doubles::from_values(tc_cov_network).into_robj(),
                Integers::from_values(tc_unique_network).into_robj(),
                Integers::from_values(tc_exp_network_ligand).into_robj(),
                Integers::from_values(tc_exp_network_receptor).into_robj(),
                Integers::from_values(tc_obs_network_ligand).into_robj(),
                Integers::from_values(tc_obs_network_receptor).into_robj(),
                Doubles::from_values(tc_cov_network_ligand).into_robj(),
                Doubles::from_values(tc_cov_network_receptor).into_robj(),
                Integers::from_values(tc_unique_network_ligand).into_robj(),
                Integers::from_values(tc_unique_network_receptor).into_robj(),
                Integers::from_values(tc_n_valid).into_robj(),
                Doubles::from_values(tc_valid_frac).into_robj(),
                Integers::from_values(tc_n_expr).into_robj(),
                Doubles::from_values(tc_n_eff_max).into_robj(),
                Doubles::from_values(tc_n_eff_min).into_robj(),
                Integers::from_values(tc_valid_max).into_robj(),
                Integers::from_values(tc_valid_min).into_robj(),
                Integers::from_values(tc_evaluable).into_robj(),
            ],
        )
        .unwrap();

        let n_fake_total = fake_pairs.len();

        List::from_names_and_values(
            [
                "fake_scores",
                "fake_pairs",
                "trajectory_calibration",
                "n_fake_total",
                "candidate_pool_size",
                "n_real_pairs",
                "n_traj_keys",
                "n_matching_blocks",
                "fake_null_mode",
                "calibration_config_hash",
                "hill_scale",
            ],
            [
                fake_scores_list.into_robj(),
                fake_pairs_info.into_robj(),
                trajectory_calibration_list.into_robj(),
                Integers::from_values([n_fake_total as i32]).into_robj(),
                Integers::from_values([n_candidate as i32]).into_robj(),
                Integers::from_values([n_pairs as i32]).into_robj(),
                Integers::from_values([n_traj_keys as i32]).into_robj(),
                Integers::from_values([n_matching_blocks as i32]).into_robj(),
                Integers::from_values([fake_null_mode]).into_robj(),
                Strings::from_values([config_hash]).into_robj(),
                Doubles::from_values([hill_s as f64]).into_robj(),
            ],
        )
        .unwrap()
    }

    // ─── 9. ComputeEvidenceDecomp ────────────

    /// Evidence decomposition: per-stratum aligned profile correlation.
    ///
    /// For each row in `real_scores`, computes the descriptive aligned profile
    /// correlation within the same (source, target, lineage, bin) stratum.
    /// `z_score` is retained as an all-NaN compatibility field.
    ///
    /// @param real_scores data.frame from `ComputeTemporalLR`/`ComputeTemporalLRWithIntermediate`.
    ///   Must have columns: pair_idx, source, target, lineage, bin,
    ///   n_eff_cells (all 1-based except n_eff_cells which is numeric).
    /// @param expr dgCMatrix (genes × cells), log-normalised.
    /// @param cell_types integer vector of length n_cells, 1-based.
    /// @param ligand_rows_list list of int vectors, per-LR-pair, 1-based gene rows.
    /// @param receptor_rows_list same.
    /// @param pseudotime numeric vector.
    /// @param lineage_weights numeric vector, flat row-major [n_cells × n_lineages].
    /// @param n_lineages, n_bins integers.
    /// @param min_pct numeric, minimum within-stratum pct for a gene to be counted.
    /// @return data.frame with aligned_r, deprecated z_score, n_eff,
    ///   n_microbins, and evidence_sufficient.
    #[extendr]
    pub fn rs_evidence_decomp(
        real_scores_pair_idx: Integers,
        real_scores_source: Integers,
        real_scores_target: Integers,
        real_scores_lineage: Integers,
        real_scores_bin: Integers,
        real_scores_n_eff: Doubles,
        expr: S4,
        cell_types: Integers,
        n_types: i32,
        ligand_rows_list: List,
        receptor_rows_list: List,
        pseudotime: Doubles,
        lineage_weights: Doubles,
        n_lineages: i32,
        n_bins: i32,
        min_pct: f64,
    ) -> List {
        use crate::ccc::{evidence_decomp_batch, EvidenceDecomp, LrBinResult};

        let n_rows = real_scores_pair_idx.len();
        assert_eq!(real_scores_source.len(), n_rows);
        assert_eq!(real_scores_target.len(), n_rows);
        assert_eq!(real_scores_lineage.len(), n_rows);
        assert_eq!(real_scores_bin.len(), n_rows);
        assert_eq!(real_scores_n_eff.len(), n_rows);

        let real_scores: Vec<LrBinResult> = (0..n_rows)
            .map(|i| LrBinResult {
                pair_idx: (real_scores_pair_idx[i].inner() - 1) as usize,
                source: (real_scores_source[i].inner() - 1) as usize,
                target: (real_scores_target[i].inner() - 1) as usize,
                lineage: (real_scores_lineage[i].inner() - 1) as usize,
                bin: (real_scores_bin[i].inner() - 1) as usize,
                score: 0.0,
                n_eff_cells: real_scores_n_eff[i].inner() as f32,
                weight_mass_cells: real_scores_n_eff[i].inner() as f32,
                lineage_entropy: 0.0,
            })
            .collect();

        let sparse = r_dgcmatrix_to_sparse(&expr);
        let ct: Vec<usize> = cell_types
            .iter()
            .map(|v| (v.inner() - 1) as usize)
            .collect();
        let nt = n_types as usize;
        let pt: Vec<f32> = pseudotime.iter().map(|v| v.inner() as f32).collect();
        let lw: Vec<f32> = lineage_weights.iter().map(|v| v.inner() as f32).collect();
        let nl = n_lineages as usize;
        let nb = n_bins as usize;

        let lig = list_of_int_to_rows0(&ligand_rows_list);
        let rec = list_of_int_to_rows0(&receptor_rows_list);

        let results: Vec<EvidenceDecomp> = evidence_decomp_batch(
            &real_scores,
            &sparse,
            &ct,
            nt,
            &pt,
            &lw,
            nl,
            nb,
            &lig,
            &rec,
            min_pct as f32,
        );

        let aligned_r: Vec<f64> = results.iter().map(|d| d.aligned_r as f64).collect();
        let z_score: Vec<f64> = results.iter().map(|d| d.z_score as f64).collect();
        let n_eff: Vec<f64> = results.iter().map(|d| d.n_eff as f64).collect();
        let n_micro: Vec<i32> = results.iter().map(|d| d.n_microbins as i32).collect();
        let sufficient: Vec<String> = results
            .iter()
            .map(|d| {
                if d.evidence_sufficient {
                    "TRUE"
                } else {
                    "FALSE"
                }
                .to_string()
            })
            .collect();

        List::from_names_and_values(
            [
                "aligned_r",
                "z_score",
                "n_eff",
                "n_microbins",
                "evidence_sufficient",
            ],
            [
                Doubles::from_values(aligned_r).into_robj(),
                Doubles::from_values(z_score).into_robj(),
                Doubles::from_values(n_eff).into_robj(),
                Integers::from_values(n_micro).into_robj(),
                Strings::from_values(sufficient).into_robj(),
            ],
        )
        .unwrap()
    }

    // ── extendr module registration ─────────────────────────────────────────
    //
    // The `mod LineaLinkR` name MUST match `useDynLib(LineaLinkR, ...)` in
    // R/NAMESPACE and `R_init_LineaLinkR` in `src/entrypoint.c`.

    extendr_module! {
        mod LineaLinkR;
        fn rs_compute_pseudotime_fuzzy;
        fn rs_compute_pseudotime_fuzzy_multi;
        fn rs_hnsw_cost_edges;
        fn rs_knn_component_count;
        fn rs_knn_component_labels;
        fn rs_attractor_unreachable_cells;
        fn rs_compute_attractor_fate_f64;
        fn rs_sparse_absorption_lu64;
        fn rs_compute_state_clock;
        fn rs_extract_lineages;
        fn rs_compute_lineage_weights;
        fn rs_compute_lineage_weights_trunk_aware;
        fn rs_extract_lineages_from_edges;
        fn rs_compute_lineage_weights_trunk_aware_from_edges;
        fn rs_compute_lineage_weights_local_branch_residual;
        fn rs_compute_fate_probabilities;
        fn rs_compute_lr_scores;
        fn rs_compute_temporal_lr;
        fn rs_compute_temporal_lr_with_inter;
        fn rs_compute_bin_bounds;
        fn rs_compute_common_bin_bounds;
        fn rs_compute_temporal_lr_fixed_bounds;
        fn rs_compute_temporal_permutation;
        fn rs_compute_lineage_entropy_raw;
        fn rs_smooth_ma3;
        fn rs_find_bifurcation_bins;
        fn rs_bootstrap_medians;
        fn rs_compute_role_scores;
        fn rs_compute_fake_lr_null;
        fn rs_evidence_decomp;
        fn rs_trajectory_de;
    }
}
