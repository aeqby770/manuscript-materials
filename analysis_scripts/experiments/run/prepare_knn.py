# -*- coding: utf-8 -*-
"""One LineaLinkR analysis preprocessing replicate (through PCA and kNN)."""
from __future__ import annotations


import sys
from pathlib import Path
_LLR_CFG = Path(__file__).resolve().parents[1] / "config"
if str(_LLR_CFG) not in sys.path:
    sys.path.insert(0, str(_LLR_CFG))
from llr_paths import assets_dir as _llr_assets_dir  # noqa: E402
_LLR_ASSETS = _llr_assets_dir()

import os
import sys

import numpy as np
from scipy.io import mmread, mmwrite
from scipy.sparse import csr_matrix, diags
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors

CONFIG = {
    "paths": {
        "workspace": _LLR_ASSETS / "cache/linealinkr_workspace",
    },
    "threads": os.cpu_count() or 1,
    "seed": 12345,
    "neighbor_seed": 0,
    "n_top_genes": 2000,
    "n_pcs": 30,
    "pearson_theta": 100,
    "n_neighbors": 50,
    "brute_threshold": 8192,
}

import pandas as pd

os.environ["OMP_NUM_THREADS"] = str(CONFIG["threads"])
np.random.seed(CONFIG["seed"])

dataset = sys.argv[1]
output_dir = Path(sys.argv[2])
matrix_path = (
    Path(sys.argv[3])
    if len(sys.argv) > 3
    else CONFIG["paths"]["workspace"]
    / dataset
    / "extract"
    / f"{dataset}_matrix.mtx"
)
n_neighbors_requested = (
    int(sys.argv[4]) if len(sys.argv) > 4 else CONFIG["n_neighbors"]
)
write_umap = len(sys.argv) > 5 and sys.argv[5] == "with_umap"
n_pcs_requested = int(sys.argv[6]) if len(sys.argv) > 6 else CONFIG["n_pcs"]
geometry_mode = sys.argv[7] if len(sys.argv) > 7 else "standard"
output_dir.mkdir(parents=True, exist_ok=True)

with matrix_path.open("rb") as stream:
    expression = mmread(stream).tocsr().T.astype(np.float32)
gene_path = matrix_path.parent / f"{dataset}_genes.txt"
genes = pd.read_csv(gene_path, header=None)[0]
exclusions = pd.read_csv(_LLR_CFG / "feature_exclusions.csv", encoding="utf-8")
excluded_genes = exclusions.loc[exclusions["dataset"].eq(dataset), "gene"]
keep_genes = ~genes.isin(excluded_genes).to_numpy()
assert len(genes) == expression.shape[1]
expression = expression[:, keep_genes]
genes = genes.loc[keep_genes].reset_index(drop=True)
library_size = np.asarray(expression.sum(axis=1)).ravel()
target_sum = float(np.median(library_size[library_size > 0]))
scale = np.divide(
    target_sum,
    library_size,
    out=np.zeros_like(library_size, dtype=np.float32),
    where=library_size > 0,
)
normalized_counts = diags(scale) @ expression
statistics_counts = normalized_counts.astype(np.float64)
gene_mean = np.asarray(statistics_counts.mean(axis=0)).ravel()
gene_variance = (
    np.asarray(statistics_counts.power(2).mean(axis=0)).ravel() - gene_mean**2
) * (normalized_counts.shape[0] / max(1, normalized_counts.shape[0] - 1))
gene_mean[gene_mean == 0] = 1e-12
gene_dispersion = gene_variance / gene_mean
gene_dispersion[gene_dispersion == 0] = np.nan
gene_dispersion = np.log(gene_dispersion)
logged_mean = np.log1p(gene_mean)
dispersion_frame = pd.DataFrame(
    {"mean": logged_mean, "dispersion": gene_dispersion}
)
dispersion_frame["mean_bin"] = pd.cut(dispersion_frame["mean"], bins=20)
bin_statistics = (
    dispersion_frame.groupby("mean_bin", observed=True)["dispersion"]
    .agg(["mean", "std"])
)
single_gene_bin = bin_statistics["std"].isna()
bin_statistics.loc[single_gene_bin, "std"] = bin_statistics.loc[
    single_gene_bin, "mean"
]
bin_statistics.loc[single_gene_bin, "mean"] = 0
bin_mean = dispersion_frame["mean_bin"].map(bin_statistics["mean"]).astype(float)
bin_std = dispersion_frame["mean_bin"].map(bin_statistics["std"]).astype(float)
dispersion_residual = gene_dispersion - bin_mean.to_numpy()
normalized_dispersion = np.full_like(dispersion_residual, np.nan, dtype=float)
np.divide(
    dispersion_residual,
    bin_std.to_numpy(),
    out=normalized_dispersion,
    where=bin_std.to_numpy() != 0,
)
n_top_genes = min(CONFIG["n_top_genes"], expression.shape[1])
finite_dispersion = normalized_dispersion[np.isfinite(normalized_dispersion)]
n_top_finite = min(n_top_genes, finite_dispersion.size)
dispersion_cutoff = np.sort(finite_dispersion)[-n_top_finite]
hvg_mask = np.nan_to_num(normalized_dispersion, nan=-np.inf) >= dispersion_cutoff
logged_expression = normalized_counts.copy()
logged_expression.data = np.log1p(logged_expression.data)
expression_hvg = logged_expression[:, hvg_mask]
if geometry_mode == "donor_center_lr":
    donor = pd.Categorical(pd.read_csv(sys.argv[8], header=None, encoding="utf-8")[0])
    lr_genes = pd.read_csv(sys.argv[9], header=None, encoding="utf-8")[0]
    lr_columns = genes.loc[hvg_mask].isin(lr_genes).to_numpy()
    donor_indicator = csr_matrix(
        (np.ones(len(donor.codes), dtype=np.float32),
         (np.arange(len(donor.codes)), donor.codes)),
        shape=(len(donor.codes), len(donor.categories)),
    )
    centered = expression_hvg.toarray()
    lr_expression = centered[:, lr_columns]
    donor_counts = np.asarray(donor_indicator.sum(axis=0)).ravel()
    donor_means = np.asarray(donor_indicator.T @ lr_expression) / donor_counts[:, None]
    centered[:, lr_columns] = (
        lr_expression - donor_means[donor.codes] + lr_expression.mean(axis=0)
    )
    expression_hvg = centered
if geometry_mode == "pearson":
    pearson_counts = expression[:, hvg_mask].toarray().astype(np.float64)
    pearson_expected = (
        pearson_counts.sum(axis=1)[:, None]
        * pearson_counts.sum(axis=0)[None, :]
        / pearson_counts.sum()
    )
    expression_hvg = np.clip(
        (pearson_counts - pearson_expected)
        / np.sqrt(
            pearson_expected
            + pearson_expected**2 / CONFIG["pearson_theta"]
            + 1e-12
        ),
        -np.sqrt(expression.shape[0]),
        np.sqrt(expression.shape[0]),
    )
n_pcs = min(n_pcs_requested, expression_hvg.shape[0] - 1, expression_hvg.shape[1])
pca = PCA(
    n_components=n_pcs,
    svd_solver="arpack",
    random_state=0,
)
pca_embedding = pca.fit_transform(expression_hvg).astype(np.float32, copy=False)

n_neighbors = min(n_neighbors_requested, max(2, expression.shape[0] - 1))
use_brute = expression.shape[0] < CONFIG["brute_threshold"]
if use_brute:
    neighbor_distance, neighbor_index = NearestNeighbors(
        n_neighbors=n_neighbors,
        algorithm="brute",
        metric="euclidean",
        n_jobs=CONFIG["threads"],
    ).fit(pca_embedding).kneighbors(pca_embedding)
    nonself = neighbor_index != np.arange(expression.shape[0])[:, None]
    nonself_rank = np.cumsum(nonself, axis=1)
    keep_neighbor = nonself & (nonself_rank < n_neighbors)
    retained_neighbors = n_neighbors - 1
    neighbor_index = neighbor_index[keep_neighbor].reshape(
        expression.shape[0], retained_neighbors
    )
    neighbor_distance = neighbor_distance[keep_neighbor].reshape(
        expression.shape[0], retained_neighbors
    )
    row_index = np.repeat(np.arange(expression.shape[0]), retained_neighbors)
    knn_distance = csr_matrix(
        (
            neighbor_distance.ravel(),
            (row_index, neighbor_index.ravel()),
        ),
        shape=(expression.shape[0], expression.shape[0]),
    )
    neighbor_backend = "sklearn_brute"
else:
    from pynndescent import PyNNDescentTransformer

    neighbor_transformer = PyNNDescentTransformer(
        n_neighbors=n_neighbors,
        metric="euclidean",
        metric_kwds={},
        random_state=CONFIG["neighbor_seed"],
        n_jobs=1,
        n_trees=min(64, 5 + round(expression.shape[0] ** 0.5 / 20.0)),
        n_iters=max(5, round(np.log2(expression.shape[0]))),
    )
    knn_distance = neighbor_transformer.fit_transform(pca_embedding).tocsr()
    # Preserve explicit zero-distance edges between duplicate PCA profiles.
    # Removing them can leave duplicate cells with no usable neighbours.
    neighbor_backend = "pynndescent_scanpy_defaults"
knn_connectivity = knn_distance.copy()
knn_connectivity.data = np.ones_like(knn_connectivity.data)

with (output_dir / "knn_dist.mtx").open("wb") as stream:
    mmwrite(stream, knn_distance)
with (output_dir / "knn_conn.mtx").open("wb") as stream:
    mmwrite(stream, knn_connectivity)
np.savetxt(output_dir / "pca.csv", pca_embedding, delimiter=",")
if write_umap:
    sys.modules["tensorflow"] = None
    from umap import UMAP

    umap_embedding = UMAP(
        n_neighbors=n_neighbors,
        n_components=2,
        metric="euclidean",
        random_state=CONFIG["seed"],
        n_jobs=1,
        low_memory=True,
    ).fit_transform(pca_embedding)
    np.savetxt(output_dir / "umap.csv", umap_embedding, delimiter=",")
print(
    f"{dataset},{expression.shape[0]},{int(hvg_mask.sum())},{neighbor_backend}"
)
