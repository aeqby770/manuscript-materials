# -*- coding: utf-8 -*-
"""Compute the independent scVelo axis with a Windows-compatible process guard."""
from __future__ import annotations


import sys
from pathlib import Path
_LLR_CFG = Path(__file__).resolve().parents[1] / "config"
if str(_LLR_CFG) not in sys.path:
    sys.path.insert(0, str(_LLR_CFG))
from llr_paths import assets_dir as _llr_assets_dir  # noqa: E402
from llr_paths import input_data_dir as _llr_input_data_dir  # noqa: E402
_LLR_ASSETS = _llr_assets_dir()
_LLR_INPUT_DATA = _llr_input_data_dir()

import os

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
from joblib import parallel_config

CONFIG = {
    "paths": {
        "input": _LLR_INPUT_DATA / "literature_benchmarks/pancreas_scvelo/pancreas_endocrinogenesis.h5ad",
        "output": _LLR_ASSETS / "results/paper_tables/supplement/tables/orthogonal_dynamic_validation/pancreas_velocity_axis.csv",
    },
    "threads": os.cpu_count() or 1,
    "seed": 12345,
    "min_shared_counts": 20,
    "n_top_genes": 2000,
    "n_pcs": 30,
    "n_neighbors": 30,
    "velocity_mode": "stochastic",
    "graph_workers": 8,
}

os.environ["OMP_NUM_THREADS"] = str(CONFIG["threads"])
np.random.seed(CONFIG["seed"])


def main() -> None:
    pancreas = ad.read_h5ad(CONFIG["paths"]["input"])
    scv.pp.filter_genes(pancreas, min_shared_counts=CONFIG["min_shared_counts"])
    scv.pp.normalize_per_cell(pancreas)
    sc.pp.log1p(pancreas)
    sc.pp.highly_variable_genes(
        pancreas,
        n_top_genes=CONFIG["n_top_genes"],
        flavor="seurat",
    )
    pancreas = pancreas[:, pancreas.var["highly_variable"]].copy()
    scv.pp.moments(
        pancreas,
        n_pcs=CONFIG["n_pcs"],
        n_neighbors=CONFIG["n_neighbors"],
    )
    scv.tl.velocity(pancreas, mode=CONFIG["velocity_mode"])
    with parallel_config(backend="loky", inner_max_num_threads=1):
        scv.tl.velocity_graph(pancreas, n_jobs=CONFIG["graph_workers"])
    scv.tl.velocity_pseudotime(pancreas)
    output = pd.DataFrame(
        {
            "cell": pancreas.obs_names.astype(str),
            "velocity_pseudotime": pancreas.obs["velocity_pseudotime"].to_numpy(float),
            "umap1": pancreas.obsm["X_umap"][:, 0],
            "umap2": pancreas.obsm["X_umap"][:, 1],
            "cluster": pancreas.obs["clusters"].astype(str).to_numpy(),
        }
    )
    CONFIG["paths"]["output"].parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(CONFIG["paths"]["output"], index=False)
    print(output.shape)


if __name__ == "__main__":
    main()
