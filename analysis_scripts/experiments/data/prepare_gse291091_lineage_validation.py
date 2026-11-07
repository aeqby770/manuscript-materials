import os
from pathlib import Path

CONFIG = {"paths": {"root": Path(os.environ["LLR_WORK_ROOT"]),
    "counts": "results/validation/gse291091/retained_counts.npz",
    "features": "results/validation/gse291091/features.tsv",
    "metadata": "results/validation/gse291091/cell_metadata.csv.gz",
    "output": "inputs/external_validation/GSE291091/GSE291091_fap_depletion.h5ad"},
    "threads": os.cpu_count(), "plot": {}}
os.environ["OMP_NUM_THREADS"] = str(CONFIG["threads"])

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

root = CONFIG["paths"]["root"]
counts = sparse.load_npz(root / CONFIG["paths"]["counts"])
features = pd.read_csv(root / CONFIG["paths"]["features"], sep="\t",
    header=None, names=["ensembl", "symbol", "feature_type"], encoding="utf-8")
metadata = pd.read_csv(root / CONFIG["paths"]["metadata"], encoding="utf-8")
codes, symbols = pd.factorize(features["symbol"], sort=False)
membership = sparse.csr_matrix((np.ones(len(features)), (codes, np.arange(len(features)))),
    shape=(len(symbols), len(features)))
collapsed = (membership @ counts).T.tocsr()
obs = metadata[["sample_id"]].set_axis(metadata["cell_id"].astype(str))
var = features.groupby("symbol", sort=False)["ensembl"].agg("|".join).to_frame()
var = var.reindex(symbols)
dataset = ad.AnnData(X=collapsed, obs=obs, var=var)
dataset.uns["source"] = "Existing retained counts; all four 3dpi samples; condition and annotation withheld"
assert dataset.n_obs == len(metadata) and dataset.obs_names.is_unique
assert dataset.n_vars == len(symbols) and dataset.var_names.is_unique
assert dataset.X.sum() == counts.sum()
output = root / CONFIG["paths"]["output"]
output.parent.mkdir(parents=True, exist_ok=True)
dataset.write_h5ad(output, compression="gzip")
loaded = ad.read_h5ad(output)
assert loaded.obs_names.equals(dataset.obs_names)
assert loaded.var_names.equals(dataset.var_names)
assert (loaded.X != dataset.X).nnz == 0
print({"cells": loaded.n_obs, "genes": loaded.n_vars,
    "metadata_columns": loaded.obs.columns.tolist(), "samples": loaded.obs["sample_id"].unique().tolist(),
    "count_sum_preserved": True})
