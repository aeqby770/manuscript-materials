#!/usr/bin/env python
"""GeneTrajectory adapter for canonical h5ad input.

Usage: python run_genetrajectory.py <canonical.h5ad> <output_gene_trajectory.csv>
"""
import os
import sys
from pathlib import Path

import numpy as np
import scanpy as sc
from scipy import sparse
import _mem_guard
_mem_guard.cap_self(55)

CONFIG = {
    "paths": {
        "vendor": Path(os.environ["LLR_VENDOR_ROOT"]).resolve() / "GeneTrajectory",
    },
    "threads": min(os.cpu_count(), 8),
    "seed": 12345,
    "n_pcs": 30,
    "graph_k": 10,
    "graph_dims": 5,
    "n_genes": 2000,
    "embedding_k": 10,
    "embedding_eigs": 30,
    "n_trajectories": 2,
    "diffusion_steps": 3,
}

os.environ["OMP_NUM_THREADS"] = str(CONFIG["threads"])
np.random.seed(CONFIG["seed"])
sys.path.insert(0, str(CONFIG["paths"]["vendor"]))

from gene_trajectory.extract_gene_trajectory import extract_gene_trajectory, get_gene_embedding  # noqa: E402
from gene_trajectory.gene_distance_shared import cal_ot_mat  # noqa: E402
from gene_trajectory.get_graph_distance import get_graph_distance  # noqa: E402


def main():
    if len(sys.argv) != 3:
        print("Usage: run_genetrajectory.py <canonical.h5ad> <output_gene_trajectory.csv>", file=sys.stderr)
        sys.exit(2)

    inpath, outpath = sys.argv[1], sys.argv[2]
    adata = sc.read_h5ad(inpath)

    # No subsampling: GeneTrajectory runs on the full cell set. The cell-cell graph
    # distance and the gene-gene EMD (cal_ot_mat) are parallelised across
    # CONFIG["threads"] rather than reduced in size.

    if "lognorm" in adata.layers:
        adata.X = adata.layers["lognorm"].copy()
    elif "counts" in adata.layers:
        adata.X = adata.layers["counts"].copy()
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
    # else: use adata.X as-is

    # NaN filter on sparse data BEFORE any dense conversion
    if sparse.issparse(adata.X):
        adata.X.data = np.nan_to_num(adata.X.data, nan=0.0, posinf=0.0, neginf=0.0)
        adata.X.data = np.clip(adata.X.data, 0.0, None)
    else:
        adata.X = np.nan_to_num(adata.X, nan=0.0, posinf=0.0, neginf=0.0)
        adata.X = np.clip(adata.X, 0.0, None)

    # Rank genes by sparse variance without densifying the full cell-by-gene matrix.
    # This avoids non-unique expression bins in the simulated datasets.
    n_top = min(3000, adata.n_vars - 1)
    if sparse.issparse(adata.X):
        # Var(X) = E[X^2] - (E[X])^2, works on sparse without per-column .data.var()
        x_csc = adata.X.tocsc()
        n = adata.n_obs
        mean_sq = np.array(x_csc.multiply(x_csc).sum(axis=0)).ravel() / n
        mean_x = np.array(x_csc.sum(axis=0)).ravel() / n
        gene_var = mean_sq - mean_x ** 2
        gene_var = np.nan_to_num(gene_var, nan=0.0)
    else:
        gene_var = np.asarray(adata.X.var(axis=0))
    top_var_idx = np.argsort(gene_var)[::-1][:n_top]
    adata = adata[:, top_var_idx].copy()
    adata.X = adata.X.toarray() if sparse.issparse(adata.X) else np.asarray(adata.X)

    n_pcs = min(CONFIG["n_pcs"], adata.n_obs - 1, adata.n_vars - 1)
    sc.pp.pca(adata, n_comps=n_pcs, random_state=CONFIG["seed"])

    x = adata.X
    gene_sum = x.sum(axis=0)
    gene_var = x.var(axis=0)
    valid_idx = np.where(gene_sum > 0)[0]
    if valid_idx.size < 3:
        print("GeneTrajectory requires at least 3 expressed genes", file=sys.stderr)
        sys.exit(1)

    ranked_valid = valid_idx[np.argsort(gene_var[valid_idx])[::-1]]
    feature_idx = ranked_valid[: min(CONFIG["n_genes"], ranked_valid.size)]
    feature_names = adata.var_names.astype(str).to_numpy()[feature_idx]
    gene_expr = x[:, feature_idx].astype(np.float64)

    # Ensure the kNN graph is connected (hard requirement of get_graph_distance):
    # pick the smallest k (from CONFIG k, doubling) that yields a single connected
    # component, decided by an EXPLICIT connected-components check rather than an
    # exception-driven retry.
    from sklearn.neighbors import kneighbors_graph
    from scipy.sparse.csgraph import connected_components
    dims_use = min(CONFIG["graph_dims"], n_pcs)
    pca_emb = adata.obsm["X_pca"][:, :dims_use]
    k = min(CONFIG["graph_k"], adata.n_obs - 1)
    while k < adata.n_obs - 1:
        knn = kneighbors_graph(pca_emb, n_neighbors=k, mode="connectivity")
        if connected_components(knn, directed=False, return_labels=False) == 1:
            break
        k = min(k * 2, adata.n_obs - 1)
    graph_dist = get_graph_distance(
        adata,
        reduction="X_pca",
        k=k,
        dims=dims_use,
    )
    gene_dist = cal_ot_mat(
        ot_cost=graph_dist,
        gene_expr=gene_expr,
        show_progress_bar=False,
        processes=float(CONFIG["threads"]),
    )
    gene_embedding, eigenvalues = get_gene_embedding(
        gene_dist,
        k=min(CONFIG["embedding_k"], gene_dist.shape[0] - 1),
        n_ev=min(CONFIG["embedding_eigs"], gene_dist.shape[0] - 1),
    )

    result = extract_gene_trajectory(
        gene_embedding=gene_embedding,
        dist_mat=gene_dist,
        gene_names=feature_names.tolist(),
        t_list=np.repeat(CONFIG["diffusion_steps"], CONFIG["n_trajectories"]),
        n=CONFIG["n_trajectories"],
        dims=min(5, gene_embedding.shape[1]),
        k=min(CONFIG["embedding_k"], gene_dist.shape[0] - 1),
    )
    result.insert(0, "gene_id", result.index.astype(str))
    result.insert(1, "gene_subset_n", len(feature_names))
    result.to_csv(outpath, index=False)

    print(f"GeneTrajectory OK: {len(result)} genes -> {outpath}")


if __name__ == "__main__":
    main()
