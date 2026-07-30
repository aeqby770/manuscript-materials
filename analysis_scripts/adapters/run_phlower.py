"""PHLOWER adapter: Hodge Laplacian pseudotime.

Usage: python run_phlower.py <canonical_h5ad> <output_csv>
"""
import os
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.decomposition import PCA

CONFIG = {
    "paths": {
        "vendor": Path(os.environ["LLR_VENDOR_ROOT"]).resolve() / "PHLOWER",
    },
    "threads": os.cpu_count() or 1,
    "seed": 12345,
    "plot": {},
}

os.environ["OMP_NUM_THREADS"] = str(CONFIG["threads"])
np.random.seed(CONFIG["seed"])
warnings.filterwarnings("ignore")
sys.path.insert(0, str(CONFIG["paths"]["vendor"]))
import phlower as ph
import _mem_guard
_mem_guard.cap_self(55)

def run_phlower(h5ad_path, output_csv):
    adata = sc.read(h5ad_path)
    # PHLOWER densifies the full matrix (toarray).
    _mem_guard.require_under(_mem_guard.dense_gb(adata.n_obs, adata.n_vars), 55, "PHLOWER")
    if "lognorm" in adata.layers:
        adata.X = adata.layers["lognorm"].copy()
    if hasattr(adata.X, "toarray"):
        X = np.asarray(adata.X.toarray(), dtype=np.float64)
    else:
        X = np.asarray(adata.X, dtype=np.float64)
    X = np.nan_to_num(X, nan=0.0, posinf=0.0, neginf=0.0)
    X = np.clip(X, 0, None)
    adata = sc.AnnData(X, obs=adata.obs, var=adata.var)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    adata = adata[:, adata.var.highly_variable].copy()
    X = np.asarray(adata.X, dtype=np.float64)
    pca = PCA(n_components=min(40, X.shape[0] - 1, X.shape[1] - 1))
    X_pca = pca.fit_transform(X)
    n_cells = X_pca.shape[0]
    adata.obsm["X_pca"] = X_pca
    # Unified cluster selection (mirrors Palantir pick_cluster_and_labels)
    group_col = None
    for cn in ("clusters", "cluster", "cell_type", "celltype", "seurat_clusters"):
        if cn not in adata.obs.columns:
            continue
        v = adata.obs[cn].astype(str)
        valid = v.notna() & (v != "")
        if v[valid].nunique() >= 2:
            group_col = cn
            break
    if group_col is not None:
        group = adata.obs[group_col].astype(str)
    else:
        # Fallback: PCA + Leiden clustering
        sc.pp.neighbors(adata, use_rep="X_pca", n_neighbors=min(30, n_cells - 1))
        sc.tl.leiden(adata, resolution=1.0, key_added="benchmark_clusters")
        group = adata.obs["benchmark_clusters"].astype(str)
        group_col = "benchmark_clusters"
        print(f"PHLOWER: using fallback Leiden clusters")
    adata.obs["group"] = group.values
    pc1 = pd.Series(X_pca[:, 0], index=adata.obs_names)
    root_group = pc1.groupby(adata.obs["group"]).median().idxmin()
    roots = (adata.obs["group"] == root_group).to_numpy()
    k = min(7, n_cells - 1)
    ph.ext.ddhodge(adata, basis="X_pca", roots=roots, k=k, npc=X_pca.shape[1], ndc=min(30, X_pca.shape[1]), s=1, lstsq_method="cholesky", layout=None, verbose=False)
    pt = adata.obs.get("u", None)
    if pt is None:
        # Do not write an all-NaN table that could be mistaken for a completed output.
        print("PHLOWER FAIL: ddhodge produced no pseudotime ('u')", file=sys.stderr)
        sys.exit(1)
    pt = np.asarray(pt, dtype=float)
    pt = (pt - np.nanmin(pt)) / max(np.nanmax(pt) - np.nanmin(pt), 1e-10)
    out = pd.DataFrame({"pseudotime": pt})
    out.to_csv(output_csv, index=False)
    n_valid = np.sum(np.isfinite(pt))
    print(f"PHLOWER OK: {n_cells} cells, {n_valid} assigned")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python run_phlower.py <canonical_h5ad> <output_csv>")
        sys.exit(1)
    run_phlower(sys.argv[1], sys.argv[2])
