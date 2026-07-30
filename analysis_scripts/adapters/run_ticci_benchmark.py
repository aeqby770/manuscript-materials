#!/usr/bin/env python
# coding: utf-8
"""
TICCI benchmark adapter (full cells, no subsampling).

Runs the vendor `Python_tools/ticci.py` workflow: single-cell entropy,
clustering, PAGA, and entropy-directed Chu-Liu/Edmonds lineage inference.
Exceptions from `ticci.getLineage` propagate to the benchmark runner.
"""
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
import _mem_guard
_mem_guard.cap_self(55)

import matplotlib
matplotlib.use("Agg")

CONFIG = {
    "paths": {
        "vendor": Path(os.environ["LLR_TICCI_ROOT"]).resolve(),
    },
    "threads": os.cpu_count() or 1,
    "seed": 12345,
    "plot": {},
}

os.environ["OMP_NUM_THREADS"] = str(CONFIG["threads"])
np.random.seed(CONFIG["seed"])
TICCI_ROOT = CONFIG["paths"]["vendor"]
if not TICCI_ROOT.is_dir():
    sys.stderr.write(f"FAIL: TICCI vendor not found at {TICCI_ROOT}\n")
    sys.exit(1)


def top_variable_genes(adata, n_top):
    """Top genes by variance; robust when seurat HVG overflows on non-model h5ad."""
    x = adata.X
    if sparse.issparse(x):
        mean = np.asarray(x.mean(axis=0)).ravel()
        sq_mean = np.asarray(x.multiply(x).mean(axis=0)).ravel()
    else:
        mean = np.asarray(x).mean(axis=0)
        sq_mean = np.square(np.asarray(x)).mean(axis=0)
    var = np.nan_to_num(sq_mean - mean ** 2, nan=0.0, posinf=0.0, neginf=0.0)
    n_top = min(int(n_top), adata.n_vars)
    idx = np.argsort(var)[::-1][:n_top]
    return adata.var_names[np.asarray(idx, dtype=int)]


def prepare_expression(adata):
    if "lognorm" in adata.layers:
        x = adata.layers["lognorm"]
        adata.X = x.toarray() if sparse.issparse(x) else np.asarray(x, dtype=np.float64)
        return
    if "counts" in adata.layers:
        x = adata.layers["counts"]
        adata.X = x.toarray() if sparse.issparse(x) else np.asarray(x, dtype=np.float64)
        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)
        return
    x = adata.X
    xmax = float(x.max()) if sparse.issparse(x) else float(np.max(x))
    if xmax > 50.0:
        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)
    elif sparse.issparse(x):
        adata.X = x.toarray().astype(np.float64, copy=False)


def select_work_adata(adata, n_top=2000):
    if adata.n_vars <= n_top:
        return adata.copy()
    n_top = min(n_top, adata.n_vars - 1)
    try:
        if "counts" in adata.layers:
            sc.pp.highly_variable_genes(
                adata, n_top_genes=n_top, flavor="seurat_v3", layer="counts"
            )
        else:
            sc.pp.highly_variable_genes(adata, n_top_genes=n_top, flavor="cell_ranger")
        if int(adata.var["highly_variable"].sum()) >= 10:
            return adata[:, adata.var.highly_variable].copy()
    except (ValueError, FloatingPointError, ZeroDivisionError):
        pass
    keep = top_variable_genes(adata, n_top)
    return adata[:, keep].copy()


def shannon_entropy_per_cell(X):
    """scEntropy-style per-cell Shannon entropy on (cells x genes) expression."""
    if sparse.issparse(X):
        X = X.toarray()
    X = np.clip(np.nan_to_num(np.asarray(X, dtype=np.float64), nan=0.0, posinf=0.0, neginf=0.0), 0.0, None)
    row_sum = X.sum(axis=1, keepdims=True)
    row_sum[row_sum == 0] = 1.0
    P = X / row_sum
    with np.errstate(divide="ignore", invalid="ignore"):
        logP = np.where(P > 0, np.log(P), 0.0)
    return -(P * logP).sum(axis=1)


def main():
    if len(sys.argv) < 3:
        print("Usage: run_ticci_benchmark.py <canonical.h5ad> <output_dir>", file=sys.stderr)
        sys.exit(2)
    inpath, outdir = sys.argv[1], sys.argv[2]
    os.makedirs(outdir, exist_ok=True)
    marker = os.path.join(outdir, "ticci_summary.csv")

    # ---- Load, full cell set (no subsampling) ----
    adata = sc.read_h5ad(inpath)
    prepare_expression(adata)
    adata.X = np.nan_to_num(adata.X, nan=0.0, posinf=0.0, neginf=0.0)
    sc.pp.filter_genes(adata, min_cells=3)

    # Single-cell entropy on the full gene set, before any HVG reduction.
    entropy = shannon_entropy_per_cell(adata.X)

    # Dimensionality reduction + graph (HVG only affects PCA/neighbors, not entropy).
    work = select_work_adata(adata, n_top=2000)
    work.obs["entropy"] = np.asarray(entropy, dtype=float)

    sc.pp.pca(work, n_comps=min(30, work.n_vars - 1, work.n_obs - 1))
    sc.pp.neighbors(work)

    # TICCI expects integer-like categories in obs['louvain']; use Leiden to
    # construct the required grouping without the deprecated louvain backend.
    sc.tl.leiden(work, key_added="louvain")
    work.obs["louvain"] = work.obs["louvain"].astype(int).astype("category")
    n_clusters = int(work.obs["louvain"].nunique())

    sc.tl.paga(work, groups="louvain")
    n_paga_edges = int(work.uns["paga"]["connectivities"].nnz)
    # PAGA layout positions are required by ticci.getLineage.
    sc.pl.paga(work, plot=False, show=False)

    # ---- TICCI lineage inference ----
    parent = TICCI_ROOT.parent
    for p in (parent, TICCI_ROOT):
        if str(p) not in sys.path:
            sys.path.insert(0, str(p))

    from Python_tools.ticci import ticci
    import matplotlib.pyplot as pl

    se_path = os.path.join(outdir, "ticci_stable_entropy.svg")
    lin_path = os.path.join(outdir, "ticci_lineage.svg")

    # Vendor nx.draw + pl.colorbar() fails headless on newer matplotlib; non-fatal for lineage.
    _orig_colorbar = pl.colorbar
    def _safe_colorbar(*args, **kwargs):
        try:
            return _orig_colorbar(*args, **kwargs)
        except ValueError:
            return None
    pl.colorbar = _safe_colorbar
    try:
        se_suffix = "_ticci_stable_entropy"
        ticci.getLineage(work, se_suffix, lin_path)
    finally:
        pl.colorbar = _orig_colorbar

    fig_dir = os.path.join(os.getcwd(), "figures")
    for ext in (".svg", ".pdf"):
        fig_src = os.path.join(fig_dir, f"paga{se_suffix}{ext}")
        if os.path.isfile(fig_src):
            import shutil
            shutil.copy2(fig_src, se_path.replace(".svg", ext))
            break
    tree = work.uns["paga"]["directed_tree"]
    rows, cols = tree.nonzero()
    weights = np.asarray(tree[rows, cols]).ravel()
    edge_df = pd.DataFrame({"from_cluster": rows, "to_cluster": cols, "weight": weights})
    edge_df.to_csv(os.path.join(outdir, "ticci_lineage_edges.csv"), index=False)
    n_tree_edges = int(len(edge_df))

    summary = {
        "n_cells": int(work.n_obs),
        "n_clusters": n_clusters,
        "n_paga_edges": n_paga_edges,
        "n_lineage_tree_edges": n_tree_edges,
        "method": "ticci_lineage",
    }
    pd.DataFrame([summary]).to_csv(marker, index=False)
    print(f"OK (ticci_lineage): {marker}")


if __name__ == "__main__":
    main()
