# -*- coding: utf-8 -*-
"""Orthogonal dynamic-axis validation using RNA velocity and withheld sampling stage."""
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
import sys

os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import anndata as ad
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
from scipy.stats import kendalltau, spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent))
from figure_style import STYLE

CONFIG = {
    "paths": {
        "pancreas_h5ad": _LLR_INPUT_DATA / "literature_benchmarks/pancreas_scvelo/pancreas_endocrinogenesis.h5ad",
        "pancreas_llr": _LLR_ASSETS / "results/paper_tables/supplement/tables/temporal_only_orthogonal_validation/pancreas_llr_pseudotime.csv",
        "pancreas_velocity": _LLR_ASSETS / "results/paper_tables/supplement/tables/orthogonal_dynamic_validation/pancreas_velocity_axis.csv",
        "gastrulation_h5ad": _LLR_INPUT_DATA / "literature_benchmarks/gastrulation_erythroid/gastrulation_erythroid.h5ad",
        "gastrulation_llr": _LLR_ASSETS / "results/paper_tables/supplement/tables/user_root_contrast/gastrulation_erythroid_blind_vs_user_pseudotime.csv",
        "out": _LLR_ASSETS / "results/paper_tables/supplement/tables/orthogonal_dynamic_validation",
        "fig": _LLR_ASSETS / "manuscript/submission_supplementary_figures",
    },
    "threads": os.cpu_count() or 1,
    "seed": 12345,
    "stage_order": ["E7.0", "E7.25", "E7.5", "E7.75", "E8.0", "E8.25", "E8.5"],
    "velocity": {
        "min_shared_counts": 20,
        "n_top_genes": 2000,
        "n_pcs": 30,
        "n_neighbors": 30,
        "mode": "stochastic",
        "n_jobs": 8,
        "blas_threads": 1,
    },
    "plot": STYLE["rc"],
}

os.environ["OMP_NUM_THREADS"] = str(CONFIG["velocity"]["blas_threads"])
os.environ["OPENBLAS_NUM_THREADS"] = str(CONFIG["velocity"]["blas_threads"])
os.environ["MKL_NUM_THREADS"] = str(CONFIG["velocity"]["blas_threads"])
np.random.seed(CONFIG["seed"])
plt.rcParams.update(CONFIG["plot"])
C = STYLE["colors"]
dynamic_cmap = LinearSegmentedColormap.from_list(
    "dynamic_axis", ["#253B73", "#23867B", "#E2B857", "#D65258"]
)
CONFIG["paths"]["out"].mkdir(parents=True, exist_ok=True)
CONFIG["paths"]["fig"].mkdir(parents=True, exist_ok=True)

pancreas_llr = (
    pd.read_csv(CONFIG["paths"]["pancreas_llr"])
    .rename(columns={"cell_id": "cell"})
    .set_index("cell")
)
pancreas_validation = pd.read_csv(CONFIG["paths"]["pancreas_velocity"]).set_index("cell")
pancreas_validation["llr_blind_pseudotime"] = pancreas_llr.loc[
    pancreas_validation.index, "pseudotime"
].to_numpy()
pancreas_validation.index.name = "cell"
pancreas_validation.to_csv(CONFIG["paths"]["out"] / "pancreas_velocity_vs_llr_cells.csv")

pancreas_rho = spearmanr(
    pancreas_validation["llr_blind_pseudotime"],
    pancreas_validation["velocity_pseudotime"],
).statistic
pancreas_tau = kendalltau(
    pancreas_validation["llr_blind_pseudotime"],
    pancreas_validation["velocity_pseudotime"],
).statistic

gastrulation = ad.read_h5ad(CONFIG["paths"]["gastrulation_h5ad"])
gastrulation_llr = pd.read_csv(CONFIG["paths"]["gastrulation_llr"]).set_index("cell")
gastrulation_validation = gastrulation.obs[["stage", "celltype"]].copy()
gastrulation_validation["stage_numeric"] = (
    gastrulation_validation["stage"].astype(str).str.replace("E", "", regex=False).astype(float)
)
gastrulation_validation[["llr_blind_pseudotime", "llr_user_root_pseudotime"]] = (
    gastrulation_llr.loc[
        gastrulation_validation.index,
        ["pseudotime_blind", "pseudotime_user_root"],
    ].to_numpy()
)
gastrulation_validation.index.name = "cell"
gastrulation_validation.to_csv(
    CONFIG["paths"]["out"] / "gastrulation_stage_vs_llr_cells.csv"
)

gastrulation_blind_rho = spearmanr(
    gastrulation_validation["llr_blind_pseudotime"],
    gastrulation_validation["stage_numeric"],
).statistic
gastrulation_blind_tau = kendalltau(
    gastrulation_validation["llr_blind_pseudotime"],
    gastrulation_validation["stage_numeric"],
).statistic
gastrulation_user_rho = spearmanr(
    gastrulation_validation["llr_user_root_pseudotime"],
    gastrulation_validation["stage_numeric"],
).statistic
gastrulation_user_tau = kendalltau(
    gastrulation_validation["llr_user_root_pseudotime"],
    gastrulation_validation["stage_numeric"],
).statistic

summary = pd.DataFrame(
    [
        {
            "dataset": "pancreas",
            "independent_axis": "scVelo stochastic velocity pseudotime",
            "llr_protocol": "blind root",
            "n_cells": len(pancreas_validation),
            "spearman": pancreas_rho,
            "kendall": pancreas_tau,
            "axis_used_by_llr": False,
            "interpretation": "orthogonal dynamic-axis concordance; not lineage tracing or causal validation",
        },
        {
            "dataset": "gastrulation_erythroid",
            "independent_axis": "withheld embryo sampling stage E7.0-E8.5",
            "llr_protocol": "blind root",
            "n_cells": len(gastrulation_validation),
            "spearman": gastrulation_blind_rho,
            "kendall": gastrulation_blind_tau,
            "axis_used_by_llr": False,
            "interpretation": "signed validation of the preregistered blind-root output",
        },
        {
            "dataset": "gastrulation_erythroid",
            "independent_axis": "withheld embryo sampling stage E7.0-E8.5",
            "llr_protocol": "user root sensitivity",
            "n_cells": len(gastrulation_validation),
            "spearman": gastrulation_user_rho,
            "kendall": gastrulation_user_tau,
            "axis_used_by_llr": False,
            "interpretation": "root-specified sensitivity; not the blind headline protocol",
        },
    ]
)
summary.to_csv(CONFIG["paths"]["out"] / "orthogonal_dynamic_validation_summary.csv", index=False)

fig, axes = plt.subplots(2, 2, figsize=(8.4, 6.4), constrained_layout=True)
umap = pancreas_validation[["umap1", "umap2"]].to_numpy()
p0 = axes[0, 0].scatter(
    umap[:, 0],
    umap[:, 1],
    c=pancreas_validation["velocity_pseudotime"],
    s=4,
    cmap=dynamic_cmap,
    linewidths=0,
    rasterized=True,
)
axes[0, 0].set_title("A  Pancreas: independent velocity axis", loc="left", fontweight="bold", color=C["ink"])
axes[0, 0].set_xlabel("UMAP1")
axes[0, 0].set_ylabel("UMAP2")
fig.colorbar(p0, ax=axes[0, 0], label="scVelo pseudotime", shrink=0.8)

axes[0, 1].hexbin(
    pancreas_validation["velocity_pseudotime"],
    pancreas_validation["llr_blind_pseudotime"],
    gridsize=45,
    mincnt=1,
    cmap="YlGnBu",
)
axes[0, 1].set_title(
    f"B  Blind LLR vs velocity (Spearman {pancreas_rho:.2f})",
    loc="left",
    fontweight="bold",
    color=C["ink"],
)
axes[0, 1].set_xlabel("scVelo stochastic pseudotime")
axes[0, 1].set_ylabel("LineaLinkR blind pseudotime")

sns.boxplot(
    data=gastrulation_validation,
    x="stage",
    y="llr_blind_pseudotime",
    color=C["blind"],
    showfliers=False,
    width=0.62,
    order=CONFIG["stage_order"],
    ax=axes[1, 0],
)
sns.stripplot(
    data=gastrulation_validation.sample(n=min(1200, len(gastrulation_validation)), random_state=CONFIG["seed"]),
    x="stage",
    y="llr_blind_pseudotime",
    color=C["ink"],
    size=1.2,
    alpha=0.12,
    order=CONFIG["stage_order"],
    ax=axes[1, 0],
    jitter=0.18,
)
axes[1, 0].set_title(
    f"C  Withheld stage, blind root (Spearman {gastrulation_blind_rho:.2f})",
    loc="left",
    fontweight="bold",
    color=C["ink"],
)
axes[1, 0].set_xlabel("Embryo sampling stage")
axes[1, 0].set_ylabel("LineaLinkR blind pseudotime")

sns.boxplot(
    data=gastrulation_validation,
    x="stage",
    y="llr_user_root_pseudotime",
    color=C["user_root"],
    showfliers=False,
    width=0.62,
    order=CONFIG["stage_order"],
    ax=axes[1, 1],
)
sns.stripplot(
    data=gastrulation_validation.sample(n=min(1200, len(gastrulation_validation)), random_state=CONFIG["seed"]),
    x="stage",
    y="llr_user_root_pseudotime",
    color=C["ink"],
    size=1.2,
    alpha=0.12,
    order=CONFIG["stage_order"],
    ax=axes[1, 1],
    jitter=0.18,
)
axes[1, 1].set_title(
    f"D  Withheld stage, user root (Spearman {gastrulation_user_rho:.2f})",
    loc="left",
    fontweight="bold",
    color=C["ink"],
)
axes[1, 1].set_xlabel("Embryo sampling stage")
axes[1, 1].set_ylabel("LineaLinkR user-root pseudotime")

fig.suptitle(
    "Independent dynamic-axis validation",
    fontsize=12,
    fontweight="bold",
    color=C["ink"],
)
figure_png = CONFIG["paths"]["fig"] / "Supplementary_Figure_S19_orthogonal_dynamic_validation.png"
figure_pdf = figure_png.with_suffix(".pdf")
fig.savefig(figure_png, bbox_inches="tight")
fig.savefig(figure_pdf, bbox_inches="tight")
plt.close(fig)

print(summary.to_string(index=False))
