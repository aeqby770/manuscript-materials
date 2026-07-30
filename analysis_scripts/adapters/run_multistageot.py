"""MultistageOT adapter for LineaLinkR benchmark pipeline.

Usage: python run_multistageot.py <canonical_h5ad_path> <output_csv_path>

Reads a canonical h5ad, identifies initial/terminal cells using clusters,
fits MultistageOT, and outputs pseudotime as CSV.
"""

import os
import sys
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse
from pathlib import Path
from _benchmark_contract import (
    benchmark_track,
    common_root_position,
    expression_root_position,
    farthest_positions,
    nearest_positions,
    pick_cluster_column,
)

VENDOR_ROOT = Path(
    os.environ.get(
        "LLR_VENDOR_ROOT",
        r"F:\biotrade\LineaLinkR_benchmark_assets\02_vendor_algorithms",
    )
)
CONFIG = {
    "paths": {
        "assets": Path(__file__).resolve().parents[2],
    },
    "threads": os.cpu_count() or 1,
    "seed": 12345,
    "plot": {},
    "memory_budget_gb": 55.0,
    "vendor_live_float64_square_arrays": 7,
    "expression_float64_copies": 2,
    "python_object_margin_gb": 2.0,
    "allocator_headroom": 1.15,
    "hvg_count": 2000,
    "pca_components": 50,
    "multistageot_epsilon": 1.0,
    "multistageot_tolerance": 1e-4,
    "multistageot_patience": 100,
}

os.environ["OMP_NUM_THREADS"] = str(CONFIG["threads"])
np.random.seed(CONFIG["seed"])
CONFIG["paths"]["vendor"] = VENDOR_ROOT / "MultistageOT"
sys.path.insert(0, str(CONFIG["paths"]["vendor"]))
from packages.models.multistageot import MultistageOT

import _mem_guard
_mem_guard.cap_self(CONFIG["memory_budget_gb"])


def _multistageot_preflight_gb(n_cells, n_genes):
    hvg_features = min(CONFIG["hvg_count"], max(1, n_genes - 1))
    square_gb = _mem_guard.square_gb(
        n_cells,
        bytes_per=np.dtype(np.float64).itemsize,
        copies=CONFIG["vendor_live_float64_square_arrays"],
    )
    expression_gb = _mem_guard.dense_gb(
        n_cells,
        hvg_features,
        bytes_per=np.dtype(np.float64).itemsize,
        copies=CONFIG["expression_float64_copies"],
    )
    return (
        square_gb
        + expression_gb
        + CONFIG["python_object_margin_gb"]
    ) * CONFIG["allocator_headroom"]


def identify_initial_cells(adata, n_init=20):
    """Identify initial cells using cluster with lowest mean pseudotime proxy.
    
    If timepoint column exists, use earliest timepoint.
    Otherwise use cluster with lowest UMAP coordinate 1.
    """
    if benchmark_track() == "common_root":
        anchor = common_root_position(adata.obs_names)
        return nearest_positions(adata.obsm["X_pca"], anchor, n_init)

    if benchmark_track() == "annotation_free":
        anchor = expression_root_position(adata)
        return nearest_positions(adata.obsm["X_pca"], anchor, n_init)

    if benchmark_track() == "legacy_annotated":
        time_columns = [column for column in ("timepoint", "age(days)", "day", "stage") if column in adata.obs.columns]
        if time_columns:
            values = pd.to_numeric(adata.obs[time_columns[0]], errors="coerce")
            valid = values.notna()
            if valid.sum() > 10:
                earliest_idx = values[valid].nsmallest(max(1, n_init)).index
                return adata.obs_names.get_indexer(earliest_idx).astype(int).tolist()
    
    # Unified cluster selection with non-empty validation
    cluster_col = pick_cluster_column(adata.obs)
    if cluster_col is None:
        # Fallback: Leiden clustering
        if "X_pca" not in adata.obsm:
            sc.pp.pca(adata, n_comps=min(30, adata.n_vars - 1, adata.n_obs - 1))
        sc.pp.neighbors(adata, use_rep="X_pca", n_neighbors=min(30, adata.n_obs - 1))
        sc.tl.leiden(adata, resolution=1.0, key_added="benchmark_clusters")
        cluster_col = "benchmark_clusters"
        print("MultistageOT: using fallback Leiden clusters")
    
    clusters = adata.obs[cluster_col].astype(str).values
    umap_key = "X_umap"
    if umap_key not in adata.obsm:
        sc.pp.neighbors(adata)
        sc.tl.umap(adata)
    
    cluster_mean = {}
    for cl in np.unique(clusters):
        mask = clusters == cl
        cluster_mean[cl] = adata.obsm[umap_key][mask, 0].mean()
    
    root_cluster = min(cluster_mean, key=cluster_mean.get)
    root_indices = np.where(clusters == root_cluster)[0]
    root_cells = root_indices[:n_init].tolist()
    return root_cells


def identify_terminal_cells(adata, n_term=20):
    """Identify terminal cells using cluster with highest UMAP coordinate 1."""
    if benchmark_track() == "common_root":
        anchor = common_root_position(adata.obs_names)
        return farthest_positions(adata.obsm["X_pca"], anchor, n_term)

    if benchmark_track() == "annotation_free":
        anchor = expression_root_position(adata)
        return farthest_positions(adata.obsm["X_pca"], anchor, n_term)

    # Reuse unified cluster selection
    cluster_col = pick_cluster_column(adata.obs)
    if cluster_col is None:
        if "benchmark_clusters" in adata.obs.columns:
            cluster_col = "benchmark_clusters"
        else:
            # Without clusters, a storage-order terminal set is not defined.
            return None
    
    clusters = adata.obs[cluster_col].astype(str).values
    umap_key = "X_umap"
    if umap_key not in adata.obsm:
        sc.pp.neighbors(adata)
        sc.tl.umap(adata)
    
    cluster_mean = {}
    for cl in np.unique(clusters):
        mask = clusters == cl
        cluster_mean[cl] = adata.obsm[umap_key][mask, 0].mean()
    
    sorted_clusters = sorted(cluster_mean, key=cluster_mean.get, reverse=True)
    terminal_cells = []
    for cl in sorted_clusters:
        cl_mask = clusters == cl
        cl_indices = np.where(cl_mask)[0].tolist()
        terminal_cells.extend(cl_indices[:max(1, n_term // len(sorted_clusters))])
        if len(terminal_cells) >= n_term:
            break
    
    return terminal_cells[:n_term]


def run_multistageot(h5ad_path, output_csv_path, n_groups=5, epsilon=CONFIG["multistageot_epsilon"]):
    """Run MultistageOT on a canonical h5ad and save pseudotime CSV."""
    
    adata = sc.read(h5ad_path)
    canonical_cell_ids = adata.obs_names.astype(str).to_numpy(copy=True)
    if canonical_cell_ids.size != adata.n_obs or pd.Index(canonical_cell_ids).has_duplicates:
        raise RuntimeError("MultistageOT: canonical obs_names must be complete and unique")
    estimated_peak_gb = _multistageot_preflight_gb(
        adata.n_obs,
        adata.n_vars,
    )
    if estimated_peak_gb > CONFIG["memory_budget_gb"]:
        raise MemoryError(
            "MultistageOT preflight: "
            f"n_cells={adata.n_obs}, n_genes={adata.n_vars}, "
            f"required={estimated_peak_gb:.2f} GiB, "
            f"budget={CONFIG['memory_budget_gb']:.2f} GiB"
        )

    if "lognorm" in adata.layers:
        X = adata.layers["lognorm"].copy()
    else:
        X = adata.X.copy()
    
    # Filter out infinity and NaN (sparse: .data, dense: whole array)
    if scipy.sparse.issparse(X):
        X.data = np.nan_to_num(X.data, nan=0.0, posinf=0.0, neginf=0.0)
    else:
        X = np.nan_to_num(np.asarray(X), nan=0.0, posinf=0.0, neginf=0.0)
    
    # Reduce to HVGs: manual top-variance selection (avoids seurat expm1 overflow
    # on non-standard lognorm, skmisc dependency, and cell_ranger duplicate bins).
    _Xh = adata.X
    if hasattr(_Xh, "power"):
        _mh = np.asarray(_Xh.mean(axis=0)).ravel()
        _vh = np.asarray(_Xh.power(2).mean(axis=0)).ravel() - _mh ** 2
    else:
        _mh = np.asarray(_Xh).mean(axis=0)
        _vh = np.square(np.asarray(_Xh)).mean(axis=0) - _mh ** 2
    _vh = np.nan_to_num(_vh, nan=0.0, posinf=0.0, neginf=0.0)
    _nh = min(CONFIG["hvg_count"], adata.n_vars - 1)
    _hvg = np.argsort(_vh)[-_nh:]
    hvg_idx = np.zeros(adata.n_vars, dtype=bool)
    hvg_idx[_hvg] = True
    X = X[:, hvg_idx]
    
    # Now convert to dense after HVG selection.
    if hasattr(X, "toarray"):
        X = X.toarray()
    
    # Build PCA for input
    from sklearn.decomposition import PCA
    pca = PCA(
        n_components=min(
            CONFIG["pca_components"],
            X.shape[1],
            X.shape[0],
        )
    )
    X_pca = pca.fit_transform(X).astype(np.float32)
    adata.obsm["X_pca"] = X_pca
    
    # Convert to DataFrame (MultistageOT expects pandas DataFrame)
    data_df = pd.DataFrame(
        X_pca,
        index=[f"cell_{i}" for i in range(X_pca.shape[0])],
        columns=[f"PC{i+1}" for i in range(X_pca.shape[1])]
    )
    
    # Identify initial and terminal cells
    init_cells = identify_initial_cells(adata, n_init=min(30, adata.n_obs // 10))
    term_cells = identify_terminal_cells(adata, n_term=min(30, adata.n_obs // 10))
    if init_cells is None or term_cells is None:
        print("FAIL: MultistageOT could not determine initial/terminal cells from "
              "the observed metadata.", file=sys.stderr)
        sys.exit(1)

    n_cells = adata.n_obs
    init_idx = [f"cell_{i}" for i in init_cells if isinstance(i, (int, np.integer)) and i < n_cells]
    term_idx = [f"cell_{i}" for i in term_cells if isinstance(i, (int, np.integer)) and i < n_cells]
    
    if len(init_idx) < 2 or len(term_idx) < 2:
        print(f"FAIL: insufficient initial ({len(init_idx)}) or terminal ({len(term_idx)}) cells",
              file=sys.stderr)
        sys.exit(1)
    
    model = MultistageOT(
        initial_cells=init_idx,
        terminal_cells=term_idx,
        n_groups=n_groups,
        epsilon=epsilon,
    )
    
    model.fit(
        data_df,
        verbose=False,
        patience=CONFIG["multistageot_patience"],
        tolerance=CONFIG["multistageot_tolerance"],
    )
    
    # Get pseudotime (index may have duplicates due to cell overlap in initial/intermediate/terminal)
    pseudotime = model.pseudotemporal_order()
    pseudotime = pseudotime.groupby(pseudotime.index).mean()  # deduplicate
    
    # Map back to cell order
    full_idx = pd.Index([f"cell_{i}" for i in range(n_cells)])
    pt_values = pd.to_numeric(pseudotime.reindex(full_idx), errors="coerce").to_numpy(dtype=float)
    if not np.array_equal(adata.obs_names.astype(str).to_numpy(), canonical_cell_ids):
        raise RuntimeError("MultistageOT: preprocessing changed canonical cell order")
    if pt_values.size != canonical_cell_ids.size:
        raise RuntimeError(
            f"MultistageOT: pseudotime length {pt_values.size} != canonical cells {canonical_cell_ids.size}"
        )
    if not np.isfinite(pt_values).all() or np.ptp(pt_values) <= 0:
        n_missing = int((~np.isfinite(pt_values)).sum())
        raise RuntimeError(
            "MultistageOT: pseudotime must map to every canonical cell with finite, "
            f"non-constant values; invalid cells={n_missing}"
        )
    
    # Save
    out = pd.DataFrame({"cell_id": canonical_cell_ids, "pseudotime": pt_values})
    out.to_csv(output_csv_path, index=False)
    written = pd.read_csv(output_csv_path, dtype={"cell_id": str})
    written_pt = pd.to_numeric(written["pseudotime"], errors="coerce").to_numpy(dtype=float)
    if (
        written.columns.tolist() != ["cell_id", "pseudotime"]
        or not np.array_equal(written["cell_id"].astype(str).to_numpy(), canonical_cell_ids)
        or written["cell_id"].astype(str).duplicated().any()
        or not np.isfinite(written_pt).all()
        or np.ptp(written_pt) <= 0
    ):
        raise RuntimeError("MultistageOT: written CSV failed the full-cell pseudotime contract")
    
    print(f"MultistageOT OK: {n_cells} cells, {len(init_idx)} initial, {len(term_idx)} terminal")
    return pseudotime


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python run_multistageot.py <canonical_h5ad> <output_csv>")
        sys.exit(1)
    
    h5ad_path = sys.argv[1]
    output_csv = sys.argv[2]
    run_multistageot(h5ad_path, output_csv)
