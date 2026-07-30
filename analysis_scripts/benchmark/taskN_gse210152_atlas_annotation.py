# -*- coding: utf-8 -*-
"""Map GSE210152 cells to the published 2025 plaque-atlas annotations by barcode."""
from __future__ import annotations


import sys
from pathlib import Path
_LLR_CFG = Path(__file__).resolve().parents[1] / "config"
if str(_LLR_CFG) not in sys.path:
    sys.path.insert(0, str(_LLR_CFG))
from llr_paths import input_data_dir, work_dir  # noqa: E402
_LLR_INPUT = input_data_dir()
_LLR_WORK = work_dir()

import os

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap

CONFIG = {
    "paths": {
        "atlas": _LLR_INPUT / "external_validation/plaque_atlas_2025/plaque_atlas_cellxgene.h5ad",
        "local_meta": _LLR_WORK / "cache/linealinkr_workspace/GSE210152_athero/extract/GSE210152_athero_meta.csv",
        "auto_clusters": _LLR_WORK / "cache/linealinkr_workspace/GSE210152_athero/GSE210152_athero_auto_clusters.csv",
        "out": _LLR_WORK / "results/biology/GSE210152",
        "fig": _LLR_WORK / "manuscript/submission_supplementary_figures",
    },
    "threads": os.cpu_count() or 1,
    "seed": 12345,
    "atlas_dataset": "Dib",
    "plot": {
        "figure.dpi": 150,
        "savefig.dpi": 300,
        "font.family": "Arial",
        "font.size": 9,
        "axes.spines.top": False,
        "axes.spines.right": False,
    },
}

os.environ["OMP_NUM_THREADS"] = str(CONFIG["threads"])
np.random.seed(CONFIG["seed"])
plt.rcParams.update(CONFIG["plot"])
CONFIG["paths"]["out"].mkdir(parents=True, exist_ok=True)

atlas = ad.read_h5ad(CONFIG["paths"]["atlas"], backed="r")
atlas_obs = atlas.obs.loc[
    atlas.obs["dataset"].astype(str) == CONFIG["atlas_dataset"],
    ["donor_id", "cell_type_level1", "cell_type_level2", "cell_type"],
].copy()
atlas_obs.index = atlas_obs.index.astype(str)
atlas_obs.index.name = "barcode_id"

local = pd.read_csv(CONFIG["paths"]["local_meta"])
auto_clusters = pd.read_csv(CONFIG["paths"]["auto_clusters"], index_col=0)
local["clusters"] = auto_clusters["auto_clusters"].to_numpy()
mapped = local.merge(atlas_obs.reset_index(), on="barcode_id", how="left", validate="one_to_one")
mapped["atlas_mapped"] = mapped["cell_type_level1"].notna()
mapped.to_csv(
    CONFIG["paths"]["out"] / "gse210152_atlas_cell_annotations.csv.gz",
    index=False,
    compression="gzip",
)

all_clusters = np.sort(mapped["clusters"].unique())

level1_counts = pd.crosstab(
    mapped.loc[mapped["atlas_mapped"], "clusters"],
    mapped.loc[mapped["atlas_mapped"], "cell_type_level1"],
).reindex(all_clusters, fill_value=0)
level1_fractions = level1_counts.div(level1_counts.sum(axis=1).replace(0, np.nan), axis=0).fillna(0)
level1_counts.to_csv(CONFIG["paths"]["out"] / "leiden_atlas_level1_counts.csv")
level1_fractions.to_csv(CONFIG["paths"]["out"] / "leiden_atlas_level1_fractions.csv")

level2_counts = pd.crosstab(
    mapped.loc[mapped["atlas_mapped"], "clusters"],
    mapped.loc[mapped["atlas_mapped"], "cell_type_level2"],
).reindex(all_clusters, fill_value=0)
level2_fractions = level2_counts.div(level2_counts.sum(axis=1).replace(0, np.nan), axis=0).fillna(0)

cluster_total = mapped.groupby("clusters", as_index=True).size()
cluster_mapped = mapped.groupby("clusters", as_index=True)["atlas_mapped"].sum()
cluster_summary = pd.DataFrame(
    {
        "leiden": all_clusters,
        "n_cells": cluster_total.loc[all_clusters].to_numpy(),
        "n_atlas_mapped": cluster_mapped.loc[all_clusters].to_numpy(),
        "atlas_mapping_fraction": (
            cluster_mapped.loc[all_clusters] / cluster_total.loc[all_clusters]
        ).to_numpy(),
        "major_level1": np.where(
            cluster_mapped.loc[all_clusters].to_numpy() > 0,
            level1_counts.idxmax(axis=1).to_numpy(),
            "unmapped",
        ),
        "major_level1_fraction": np.where(
            cluster_mapped.loc[all_clusters].to_numpy() > 0,
            level1_fractions.max(axis=1).to_numpy(),
            np.nan,
        ),
        "major_level2": np.where(
            cluster_mapped.loc[all_clusters].to_numpy() > 0,
            level2_counts.idxmax(axis=1).to_numpy(),
            "unmapped",
        ),
        "major_level2_fraction": np.where(
            cluster_mapped.loc[all_clusters].to_numpy() > 0,
            level2_fractions.max(axis=1).to_numpy(),
            np.nan,
        ),
    }
)
cluster_summary["annotation_status"] = np.select(
    [
        cluster_summary["atlas_mapping_fraction"] < 0.1,
        (cluster_summary["atlas_mapping_fraction"] >= 0.5)
        & (cluster_summary["major_level1_fraction"] >= 0.5),
    ],
    ["insufficient_mapping", "reference_supported_majority"],
    default="partial_or_mixed",
)
cluster_summary.to_csv(
    CONFIG["paths"]["out"] / "leiden_atlas_annotation_summary.csv",
    index=False,
)
cluster_summary.to_csv(
    CONFIG["paths"]["out"] / "gse210152_annotation_status_064.csv",
    index=False,
)

mapping_summary = pd.DataFrame(
    [
        {
            "local_cells": len(mapped),
            "atlas_dib_cells": len(atlas_obs),
            "mapped_cells": int(mapped["atlas_mapped"].sum()),
            "mapping_fraction": float(mapped["atlas_mapped"].mean()),
            "n_local_leiden": mapped["clusters"].nunique(),
            "n_level1_labels": mapped.loc[mapped["atlas_mapped"], "cell_type_level1"].nunique(),
            "n_level2_labels": mapped.loc[mapped["atlas_mapped"], "cell_type_level2"].nunique(),
            "annotation_source": "Traeuble et al. 2025 plaque atlas, CELLxGENE dataset Dib; exact barcode join",
            "boundary": "unmapped local cells remain unlabeled; atlas labels are reference annotations, not original GEO obs fields",
        }
    ]
)
mapping_summary.to_csv(
    CONFIG["paths"]["out"] / "gse210152_atlas_mapping_summary.csv",
    index=False,
)

heatmap_data = level1_fractions.loc[:, level1_fractions.sum(axis=0).sort_values(ascending=False).index]
composition_cmap = LinearSegmentedColormap.from_list(
    "composition", ["#F7F3E8", "#E2B857", "#23867B"]
)
fig, axes = plt.subplots(1, 3, figsize=(11.2, 4.5), constrained_layout=True)
sns.heatmap(
    heatmap_data,
    vmin=0,
    vmax=1,
    cmap=composition_cmap,
    linewidths=0.4,
    linecolor="white",
    cbar_kws={"label": "Fraction among atlas-mapped cells", "shrink": 0.75},
    ax=axes[0],
)
axes[0].set_xlabel("Plaque-atlas level-1 label")
axes[0].set_ylabel("Local Leiden cluster")
axes[0].set_title("A  Level-1 composition", loc="left", fontweight="bold")
axes[0].set_xticklabels(
    heatmap_data.columns,
    rotation=55,
    ha="right",
    fontsize=7,
)

axes[1].barh(
    cluster_summary["leiden"].astype(str),
    cluster_summary["atlas_mapping_fraction"],
    color="#4E79A7",
)
axes[1].set_xlim(0, 1)
axes[1].set_xlabel("Exact-barcode mapping fraction")
axes[1].set_ylabel("Local Leiden cluster")
axes[1].set_title("B  Mapping coverage", loc="left", fontweight="bold")
axes[1].invert_yaxis()

axes[2].barh(
    cluster_summary["leiden"].astype(str),
    cluster_summary["major_level1_fraction"],
    color="#23867B",
)
axes[2].set_xlim(0, 1)
axes[2].set_xlabel("Major-label fraction (mapped cells)")
axes[2].set_ylabel("Local Leiden cluster")
axes[2].set_title("C  Cluster annotation purity", loc="left", fontweight="bold")
axes[2].invert_yaxis()
axes[2].bar_label(
    axes[2].containers[0],
    labels=cluster_summary["major_level1"],
    padding=3,
    fontsize=7,
)

figure_png = CONFIG["paths"]["fig"] / "Supplementary_Figure_S17_gse210152_atlas_annotation.png"
figure_pdf = figure_png.with_suffix(".pdf")
fig.savefig(figure_png, bbox_inches="tight")
fig.savefig(figure_pdf, bbox_inches="tight")
plt.close(fig)

print(mapping_summary.to_string(index=False))
print(cluster_summary.to_string(index=False))
