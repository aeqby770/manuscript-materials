import os

import numpy as np
import pandas as pd
from scipy import sparse


CONFIG = {
    "paths": {},
    "threads": os.cpu_count() or 1,
    "seed": 12345,
    "plot": {},
}

os.environ["OMP_NUM_THREADS"] = str(CONFIG["threads"])
np.random.seed(CONFIG["seed"])

BENCHMARK_TRACKS = ("annotation_free", "common_root", "legacy_annotated")
ANNOTATION_FREE_COLUMNS = ("benchmark_clusters",)
LEGACY_CLUSTER_COLUMNS = (
    "clusters",
    "cluster",
    "cell_type",
    "celltype",
    "seurat_clusters",
)


def benchmark_track():
    track = os.environ.get("LLR_BENCHMARK_TRACK", "annotation_free")
    if track not in BENCHMARK_TRACKS:
        raise ValueError(f"Unsupported LLR_BENCHMARK_TRACK: {track}")
    return track


def pick_cluster_column(obs):
    candidates = (
        LEGACY_CLUSTER_COLUMNS
        if benchmark_track() == "legacy_annotated"
        else ANNOTATION_FREE_COLUMNS
    )
    usable = [
        column
        for column in candidates
        if column in obs.columns
        and obs[column].astype("string").fillna("").loc[lambda values: values != ""].nunique() >= 2
    ]
    return usable[0] if usable else None


def common_root_position(obs_names):
    root_cell = os.environ["LLR_BENCHMARK_ROOT_CELL"]
    position = obs_names.get_indexer([root_cell])[0]
    if position < 0:
        raise ValueError(f"Common-root cell is absent from input: {root_cell}")
    return int(position)


def detected_gene_counts(adata):
    matrix = adata.layers["counts"] if "counts" in adata.layers else adata.X
    if sparse.issparse(matrix):
        return np.asarray((matrix != 0).sum(axis=1)).ravel().astype(float)
    return np.count_nonzero(np.asarray(matrix), axis=1).astype(float)


def expression_root_position(adata, labels=None):
    detected = detected_gene_counts(adata)
    positions = np.arange(adata.n_obs, dtype=int)
    cluster_labels = labels
    if cluster_labels is None:
        cluster_column = pick_cluster_column(adata.obs)
        cluster_labels = (
            adata.obs[cluster_column].astype("string").fillna("").to_numpy()
            if cluster_column is not None
            else None
        )

    if cluster_labels is None:
        candidates = positions[detected == np.nanmax(detected)]
    else:
        frame = pd.DataFrame(
            {
                "label": np.asarray(cluster_labels, dtype=str),
                "detected": detected,
                "position": positions,
            }
        )
        summary = (
            frame.groupby("label", sort=False, observed=True)
            .agg(
                detected_median=("detected", "median"),
                detected_mean=("detected", "mean"),
                cluster_size=("detected", "size"),
                first_position=("position", "min"),
            )
            .sort_values(
                ["detected_median", "detected_mean", "cluster_size", "first_position"],
                ascending=[False, False, False, True],
                kind="stable",
            )
        )
        root_label = str(summary.index[0])
        candidates = positions[frame["label"].to_numpy() == root_label]

    if "X_pca" not in adata.obsm:
        candidate_detected = detected[candidates]
        return int(candidates[np.flatnonzero(candidate_detected == np.nanmax(candidate_detected))[0]])

    embedding = np.asarray(adata.obsm["X_pca"], dtype=float)
    candidate_embedding = embedding[candidates]
    centroid = np.nanmean(candidate_embedding, axis=0)
    distances = np.linalg.norm(candidate_embedding - centroid, axis=1)
    return int(candidates[np.nanargmin(distances)])


def root_principal_point(adata, root_position=None):
    if root_position is None:
        root_position = common_root_position(adata.obs_names)
    membership = np.asarray(adata.obsm["X_R"])[root_position]
    tips = np.asarray(adata.uns["graph"].get("tips", []), dtype=int)
    candidates = tips if tips.size else np.arange(membership.size, dtype=int)
    return int(candidates[np.nanargmax(membership[candidates])])


def nearest_positions(embedding, anchor_position, count):
    distances = np.linalg.norm(embedding - embedding[anchor_position], axis=1)
    return np.argsort(distances, kind="stable")[:count].astype(int).tolist()


def farthest_positions(embedding, anchor_position, count):
    distances = np.linalg.norm(embedding - embedding[anchor_position], axis=1)
    return np.argsort(-distances, kind="stable")[:count].astype(int).tolist()
