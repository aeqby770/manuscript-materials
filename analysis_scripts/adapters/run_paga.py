"""PAGA adapter: compute PAGA connectivity matrix.
Usage: python run_paga.py <canonical_h5ad> <output_csv> [metadata|leiden]
"""
import sys, os
import numpy as np
import pandas as pd
import scanpy as sc
import _mem_guard
_mem_guard.cap_self(55)

CONFIG = {
    "seed": 12345,
    "threads": os.cpu_count(),
    "default_neighbors": 15,
    "leiden_neighbors": 50,
    "leiden_resolution": 1.0,
}

os.environ["OMP_NUM_THREADS"] = str(CONFIG["threads"])
np.random.seed(CONFIG["seed"])


def run_paga(h5ad_path, output_csv, group_mode="metadata"):
    adata = sc.read(h5ad_path)
    
    if "counts" in adata.layers:
        adata.X = adata.layers["counts"].copy()
    if hasattr(adata.X, "toarray"):
        adata.X = adata.X.toarray()
    adata.X = np.nan_to_num(adata.X, nan=0, posinf=0, neginf=0)
    
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    adata_hvg = adata[:, adata.var.highly_variable].copy()
    sc.pp.pca(adata_hvg, n_comps=30)

    if group_mode == "leiden":
        sc.pp.neighbors(adata_hvg, n_neighbors=CONFIG["leiden_neighbors"])
        sc.tl.leiden(
            adata_hvg,
            resolution=CONFIG["leiden_resolution"],
            key_added="paga_group",
        )
        group_col = "paga_group"
    else:
        sc.pp.neighbors(adata_hvg, n_neighbors=CONFIG["default_neighbors"])
        group_col = None
        for c in ['clusters', 'cluster', 'cell_type', 'celltype', 'seurat_clusters']:
            if c not in adata_hvg.obs.columns:
                continue
            v = adata_hvg.obs[c].astype(str)
            valid = v.notna() & (v != "")
            if v[valid].nunique() >= 2:
                group_col = c
                break
        if group_col is None:
            # Fallback: Leiden clustering
            sc.tl.leiden(adata_hvg, resolution=1.0, key_added="paga_group")
            group_col = "paga_group"
            print("PAGA: using fallback Leiden clusters")
        adata_hvg.obs['paga_group'] = adata_hvg.obs[group_col].astype(str)

    adata_hvg.obs['paga_group'] = adata_hvg.obs['paga_group'].astype('category')
    sc.tl.paga(adata_hvg, groups='paga_group')
    
    connectivities = adata_hvg.uns['paga']['connectivities']
    labels = adata_hvg.obs['paga_group'].cat.categories.tolist()
    
    connectivities = adata_hvg.uns['paga']['connectivities']
    labels = adata_hvg.obs['paga_group'].cat.categories.tolist()
    
    conn_df = pd.DataFrame(connectivities.toarray() if hasattr(connectivities, 'toarray') else connectivities,
                          index=labels, columns=labels)
    conn_df.to_csv(output_csv)
    
    n_groups = len(labels)
    n_edges = int((connectivities.toarray() > 0).sum() // 2) if hasattr(connectivities, 'toarray') else 0
    print(f"PAGA OK: {n_groups} groups, {n_edges} edges, group_col={group_col}, group_mode={group_mode}")

if __name__ == "__main__":
    if len(sys.argv) not in (3, 4):
        print("Usage: python run_paga.py <canonical_h5ad> <output_csv> [metadata|leiden]")
        sys.exit(1)
    run_paga(sys.argv[1], sys.argv[2], sys.argv[3] if len(sys.argv) == 4 else "metadata")
