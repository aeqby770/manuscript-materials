"""DPT (Diffusion Pseudotime) adapter using the official Scanpy API.

Usage: python run_dpt.py <canonical_h5ad> <output_csv> [root_strategy]

Root selection strategy:
- "auto": cluster-based root selection (consistent with scanpy tutorial).
  Picks root cluster by label semantics ("0"/"root"/"stem"/min numeric),
  then selects the PCA-centroid cell within that cluster.
  Falls back to Leiden clustering if no annotation column exists.
  Ground-truth columns, including true_pseudotime, are excluded.

Official reference: scanpy.tl.dpt(adata, n_dcs=10)
- https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.dpt.html
"""
import sys, os
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from scipy.sparse.csgraph import connected_components, minimum_spanning_tree
from scipy.spatial.distance import cdist
import _mem_guard
from _benchmark_contract import expression_root_position
_mem_guard.cap_self(55)

CONFIG = {
    "paths": {},
    "seed": 12345,
    "threads": os.cpu_count(),
    "plot": {},
    "n_dcs": 10,
    "n_neighbors": 15,
    "n_pcs": 30,
    "component_bridge_percentile": 10.0,
}

os.environ["OMP_NUM_THREADS"] = str(CONFIG["threads"])
np.random.seed(CONFIG["seed"])


def benchmark_track():
    track = os.environ.get("LLR_BENCHMARK_TRACK", "annotation_free")
    if track not in {"annotation_free", "common_root", "legacy_annotated"}:
        raise ValueError(f"Unknown LLR_BENCHMARK_TRACK: {track}")
    return track


def pick_cluster_col(ad):
    """Unified cluster column selection with non-empty validation (mirrors Palantir)."""
    candidates = (
        ("clusters", "cluster", "cell_type", "celltype", "seurat_clusters")
        if benchmark_track() == "legacy_annotated"
        else ("benchmark_clusters",)
    )
    for cn in candidates:
        if cn not in ad.obs.columns:
            continue
        v = ad.obs[cn].astype(str)
        valid = v.notna() & (v != "")
        if v[valid].nunique() >= 2:
            return cn
    return None


def attach_fallback_clusters(ad):
    """Derive clusters from PCA+Leiden when no annotation available (mirrors Palantir)."""
    if "X_pca" not in ad.obsm:
        sc.pp.pca(ad, n_comps=min(30, ad.n_vars - 1, ad.n_obs - 1))
    sc.pp.neighbors(ad, use_rep="X_pca", n_neighbors=min(30, ad.n_obs - 1))
    sc.tl.leiden(ad, resolution=1.0, key_added="benchmark_clusters")
    return "benchmark_clusters"


ROOT_CANDIDATES = ("0", "root", "Root", "STEM", "stem", "branch_0", "early", "Early")


def determine_root_label(labels):
    """Choose root cluster label consistent with PickMonocle3RootCells semantics."""
    cl = pd.Series(labels.astype(str))
    root = None
    for candid in ROOT_CANDIDATES:
        if (cl == candid).any():
            root = candid
            break
    if root is None:
        num = pd.to_numeric(cl, errors="coerce")
        if np.isfinite(num).any():
            root = str(np.nanmin(num.values))
    if root is None and len(cl) > 0:
        vc = cl.value_counts()
        if len(vc) >= 1:
            root = str(vc.index[0])
    return root


def _dense_finite(x):
    if sparse.issparse(x):
        x = x.toarray()
    x = np.asarray(x, dtype=np.float64)
    x = np.nan_to_num(x, nan=0.0, posinf=0.0, neginf=0.0)
    pos = x[x > 0]
    if pos.size:
        cap = float(np.percentile(pos, 99.9))
        if np.isfinite(cap) and cap > 0:
            x = np.clip(x, 0, cap)
    return x


def _sanitize_x(adata):
    if sparse.issparse(adata.X):
        adata.X.data = np.nan_to_num(
            adata.X.data, nan=0.0, posinf=0.0, neginf=0.0
        )
        adata.X.eliminate_zeros()
        return adata
    adata.X = _dense_finite(adata.X)
    return adata


def _is_lognorm_matrix(x):
    if sparse.issparse(x):
        sample = x.data
        if sample.size == 0:
            sample = x.toarray().ravel()
    else:
        sample = np.asarray(x).ravel()
    sample = sample[np.isfinite(sample)]
    if sample.size == 0:
        return True
    # log-normalized: non-integer and modest dynamic range
    return (
        np.max(sample) < 50
        and not np.allclose(sample, np.round(sample), rtol=0, atol=1e-6)
    )


def _select_hvg(adata, n_top_genes):
    if sparse.issparse(adata.X):
        gene_means = np.asarray(adata.X.mean(axis=0)).ravel()
    else:
        gene_means = np.asarray(adata.X.mean(axis=0)).ravel()
    finite_genes = np.isfinite(gene_means)
    np.testing.assert_equal(finite_genes.all(), True)
    n_top_genes = min(n_top_genes, adata.n_vars - 1)
    sc.pp.highly_variable_genes(
        adata, n_top_genes=n_top_genes, flavor="seurat_v3"
    )
    return adata[:, adata.var.highly_variable].copy()


def select_root_cell(adata, strategy="auto"):
    """Select the root cell and return its zero-based index.

    Strategy "auto" uses cluster-based root selection, following the Scanpy
    tutorial. The true_pseudotime column is excluded from root selection.
    """
    if benchmark_track() == "common_root":
        root_cell = os.environ["LLR_BENCHMARK_ROOT_CELL"]
        if root_cell not in adata.obs_names:
            raise ValueError(f"Common-root cell is absent from input: {root_cell}")
        return int(adata.obs_names.get_loc(root_cell))

    cluster_col = pick_cluster_col(adata)
    if cluster_col is None:
        cluster_col = attach_fallback_clusters(adata)
        print(f"DPT: using fallback Leiden clusters for root selection")

    clusters = adata.obs[cluster_col].astype(str).values
    if benchmark_track() == "annotation_free":
        root_idx = expression_root_position(adata, clusters)
        print(f"Root selection: shared expression-only anchor, cell_idx={root_idx}")
        return root_idx

    root_label = determine_root_label(clusters)
    if root_label is None:
        print("Root selection: fallback to cell_idx=0 (no root cluster determined)")
        return 0

    mask = clusters == root_label
    if not np.any(mask):
        print(f"Root selection: fallback to cell_idx=0 (root label '{root_label}' has no cells)")
        return 0

    root_cells = np.where(mask)[0]
    if "X_pca" in adata.obsm:
        pca = adata.obsm["X_pca"][root_cells, :]
        centroid = pca.mean(axis=0)
        dists = np.linalg.norm(pca - centroid, axis=1)
        root_idx = int(root_cells[np.argmin(dists)])
    else:
        root_idx = int(root_cells[0])
    print(f"Root selection: root_label='{root_label}' center, cell_idx={root_idx}")
    return root_idx


def run_dpt(h5ad_path, output_csv, root_strategy="auto"):
    adata = sc.read(h5ad_path)
    # DPT densifies the expression matrix (toarray) for PCA/diffmap.
    _mem_guard.require_under(_mem_guard.dense_gb(adata.n_obs, adata.n_vars, copies=2), 55, "DPT")

    # Prepare data following the official Scanpy workflow.
    if "lognorm" in adata.layers:
        adata.X = adata.layers["lognorm"].copy()
    elif "counts" in adata.layers:
        adata.X = adata.layers["counts"].copy()
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
    else:
        adata.X = adata.X.copy()
        if not _is_lognorm_matrix(adata.X):
            sc.pp.normalize_total(adata, target_sum=1e4)
            sc.pp.log1p(adata)

    adata = _sanitize_x(adata)
    sc.pp.filter_genes(adata, min_cells=3)
    adata = _sanitize_x(adata)

    n_top = min(2000, adata.n_vars - 1)
    if adata.n_vars > 3000:
        adata = _select_hvg(adata, n_top)
    else:
        print(f"Using all {adata.n_vars} genes (< 3000 threshold)")

    # The DPT path requires a dense representation before PCA.
    adata.X = _dense_finite(adata.X)
    
    # PCA and neighborhood graph (official workflow).
    sc.pp.pca(adata, n_comps=CONFIG["n_pcs"])
    sc.pp.neighbors(adata, n_neighbors=CONFIG["n_neighbors"])

    # Root selection under the documented protocol.
    root_idx = select_root_cell(adata, root_strategy)

    # DPT is infinite outside the root's connected component. Connect the
    # component centroids by their PCA-space MST, adding only C - 1 graph edges.
    # This restores a single diffusion geometry without using any GT columns.
    connectivities = adata.obsp["connectivities"].tocsr()
    n_components, component_labels = connected_components(
        connectivities, directed=False
    )
    component_sizes = np.bincount(component_labels, minlength=n_components)
    component_membership = sparse.csr_matrix(
        (
            np.ones(adata.n_obs, dtype=np.float64),
            (component_labels, np.arange(adata.n_obs)),
        ),
        shape=(n_components, adata.n_obs),
    )
    pca = np.asarray(adata.obsm["X_pca"], dtype=np.float64)
    component_centroids = (
        np.asarray(component_membership @ pca) / component_sizes[:, None]
    )
    component_distances = np.maximum(
        cdist(component_centroids, component_centroids),
        np.finfo(np.float64).eps,
    )
    np.fill_diagonal(component_distances, 0.0)
    component_tree = minimum_spanning_tree(component_distances)
    source_components, target_components = component_tree.nonzero()

    cell_norm = np.einsum("ij,ij->i", pca, pca)
    centroid_norm = np.einsum(
        "ij,ij->i", component_centroids, component_centroids
    )
    cell_to_centroid_sq = np.maximum(
        cell_norm[:, None]
        + centroid_norm[None, :]
        - 2.0 * pca @ component_centroids.T,
        0.0,
    )
    bridge_sources = np.argmin(
        np.where(
            component_labels[:, None] == source_components[None, :],
            cell_to_centroid_sq[:, target_components],
            np.inf,
        ),
        axis=0,
    )
    bridge_targets = np.argmin(
        np.where(
            component_labels[:, None] == target_components[None, :],
            cell_to_centroid_sq[:, source_components],
            np.inf,
        ),
        axis=0,
    )
    bridge_distances = np.linalg.norm(
        pca[bridge_sources] - pca[bridge_targets], axis=1
    )
    bridge_weights = np.full(
        bridge_sources.size,
        np.percentile(
            connectivities.data, CONFIG["component_bridge_percentile"]
        ),
    )
    connectivity_bridges = sparse.coo_matrix(
        (
            np.concatenate((bridge_weights, bridge_weights)),
            (
                np.concatenate((bridge_sources, bridge_targets)),
                np.concatenate((bridge_targets, bridge_sources)),
            ),
        ),
        shape=connectivities.shape,
    ).tocsr()
    adata.obsp["connectivities"] = connectivities.maximum(
        connectivity_bridges
    )

    distances = adata.obsp["distances"].tocsr()
    distance_bridges = sparse.coo_matrix(
        (
            np.concatenate((bridge_distances, bridge_distances)),
            (
                np.concatenate((bridge_sources, bridge_targets)),
                np.concatenate((bridge_targets, bridge_sources)),
            ),
        ),
        shape=distances.shape,
    ).tocsr()
    adata.obsp["distances"] = distances.maximum(distance_bridges)
    n_components_after, _ = connected_components(
        adata.obsp["connectivities"], directed=False
    )
    np.testing.assert_equal(n_components_after, 1)
    print(
        f"DPT graph: {n_components} components, "
        f"{bridge_sources.size} PCA-MST bridges, "
        f"{n_components_after} component after repair"
    )

    sc.tl.diffmap(adata)
    adata.uns["iroot"] = root_idx
    
    # DPT (official API).
    sc.tl.dpt(adata, n_dcs=CONFIG["n_dcs"])
    
    # Output
    out = pd.DataFrame({
        "cell_id": adata.obs_names,
        "pseudotime": adata.obs["dpt_pseudotime"].values
    })
    n_valid = np.sum(np.isfinite(out["pseudotime"]))
    np.testing.assert_equal(n_valid, adata.n_obs)
    out.to_csv(output_csv, index=False)

    print(f"DPT OK: {adata.n_obs} cells, {n_valid} valid pseudotime, root_idx={root_idx}, root_strategy={root_strategy}")
    return out


if __name__ == "__main__":
    if len(sys.argv) not in (3, 4):
        print(__doc__)
        sys.exit(1)
    run_dpt(sys.argv[1], sys.argv[2], sys.argv[3] if len(sys.argv) == 4 else "auto")
