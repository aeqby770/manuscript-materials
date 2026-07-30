# -*- coding: utf-8 -*-
"""Trajectory-level descriptive calibration of LineaLinkR role weights."""
from __future__ import annotations

import os
import sys
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap

_LLR_CFG = Path(__file__).resolve().parents[1] / "config"
sys.path.insert(0, str(_LLR_CFG))
from llr_paths import assets_dir as _llr_assets_dir  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from figure_style import STYLE  # noqa: E402

_LLR_ASSETS = _llr_assets_dir()

CONFIG = {
    "paths": {
        "full": _LLR_ASSETS / (
            "results/biology/taxonomy/pancreas_taxonomy_full.csv"
        ),
        "output": _LLR_ASSETS / (
            "results/biology/taxonomy"
        ),
        "figure": _LLR_ASSETS / "manuscript/submission_supplementary_figures",
    },
    "threads": os.cpu_count() or 1,
    "seed": 12345,
    "plot": STYLE["rc"],
    "roles": ["instructive", "reinforcing", "feedback", "lateral", "bystander"],
}

os.environ["OMP_NUM_THREADS"] = str(CONFIG["threads"])
np.random.seed(CONFIG["seed"])
plt.rcParams.update(CONFIG["plot"])
CONFIG["paths"]["output"].mkdir(parents=True, exist_ok=True)
CONFIG["paths"]["figure"].mkdir(parents=True, exist_ok=True)

full = pd.read_csv(CONFIG["paths"]["full"])
ligand = full["ligand_genes"].astype(str)
full["pathway"] = np.select(
    [
        ligand.str.match(r"(?i)^FGF"),
        ligand.str.match(r"(?i)^WNT"),
        ligand.isin(["Dlk1", "Dll1", "Dll3", "Dll4", "Jag1", "Jag2"]),
    ],
    ["FGF", "Wnt", "Notch"],
    default="",
)
structural_bins = full[
    full["pathway"].ne("") & full["structural_valid"].fillna(False)
].copy()
bins = structural_bins[
    structural_bins["calibration_eligible"].fillna(False).astype(bool)
].copy()

keys = ["pathway", "pair_idx", "source", "target", "lineage"]
trajectory_scores = (
    bins.groupby(keys, observed=True)[CONFIG["roles"]]
    .mean()
    .reset_index()
)
trajectory_evidence = (
    structural_bins.groupby(keys, observed=True)
    .agg(
        n_structurally_valid_bins=("structural_valid", "size"),
        n_expression_supported_bins=("expression_supported", "sum"),
        n_calibration_eligible_bins=("calibration_eligible", "sum"),
        any_expression_supported=("expression_supported", "any"),
    )
    .reset_index()
)
trajectories = trajectory_scores.merge(
    trajectory_evidence,
    on=keys,
    how="left",
    validate="one_to_one",
)
trajectories["dominant_role"] = trajectories[CONFIG["roles"]].idxmax(axis=1)
trajectories["expected_role"] = trajectories["pathway"].map(
    {"FGF": "instructive", "Wnt": "instructive", "Notch": "feedback"}
)
trajectories["coarse_prior_match"] = (
    trajectories["dominant_role"] == trajectories["expected_role"]
)
trajectories["expected_role_weight"] = np.select(
    [
        trajectories["pathway"].isin(["FGF", "Wnt"]),
        trajectories["pathway"].eq("Notch"),
    ],
    [trajectories["instructive"], trajectories["feedback"]],
    default=np.nan,
)
trajectories["max_alternative_weight"] = np.where(
    trajectories["pathway"].eq("Notch"),
    trajectories[["instructive", "reinforcing", "lateral", "bystander"]].max(axis=1),
    trajectories[["reinforcing", "feedback", "lateral", "bystander"]].max(axis=1),
)
trajectories["expected_role_margin"] = (
    trajectories["expected_role_weight"] - trajectories["max_alternative_weight"]
)

role_counts = pd.crosstab(trajectories["pathway"], trajectories["dominant_role"])
role_counts = role_counts.reindex(
    index=["FGF", "Wnt", "Notch"],
    columns=CONFIG["roles"],
    fill_value=0,
)
role_counts.to_csv(
    CONFIG["paths"]["output"] / "pancreas_pathway_role_confusion_matrix.csv"
)

calibration = (
    trajectories.groupby("pathway", observed=True)
    .agg(
        n_calibration_eligible_trajectories=("coarse_prior_match", "size"),
        n_match_expected=("coarse_prior_match", "sum"),
        recall_expected=("coarse_prior_match", "mean"),
        n_expression_supported_trajectories=("any_expression_supported", "sum"),
        median_expected_role_weight=("expected_role_weight", "median"),
        median_expected_role_margin=("expected_role_margin", "median"),
    )
    .reset_index()
)
calibration["expected_role"] = calibration["pathway"].map(
    {"FGF": "instructive", "Wnt": "instructive", "Notch": "feedback"}
)
calibration["calibration_unit"] = "pair_source_target_lineage_trajectory"
calibration["calibration_population"] = (
    "expression_supported_finite_feature_trajectory"
)
calibration = calibration[
    [
        "pathway",
        "n_calibration_eligible_trajectories",
        "n_match_expected",
        "recall_expected",
        "expected_role",
        "n_expression_supported_trajectories",
        "median_expected_role_weight",
        "median_expected_role_margin",
        "calibration_unit",
        "calibration_population",
    ]
]
calibration.to_csv(
    CONFIG["paths"]["output"] / "pancreas_pathway_role_precision_recall.csv",
    index=False,
)
calibration.to_csv(
    CONFIG["paths"]["output"] / "taxonomy_prior_alignment_064.csv",
    index=False,
)

expression_trajectories = trajectories[trajectories["any_expression_supported"]]
expression_calibration = (
    expression_trajectories.groupby("pathway", observed=True)
    .agg(
        n_trajectories=("coarse_prior_match", "size"),
        n_matching_coarse_prior=("coarse_prior_match", "sum"),
        coarse_prior_agreement=("coarse_prior_match", "mean"),
        median_expected_role_weight=("expected_role_weight", "median"),
        median_expected_role_margin=("expected_role_margin", "median"),
    )
    .reset_index()
)
expression_calibration["calibration_unit"] = "expression_supported_trajectory"

trajectories.to_csv(
    CONFIG["paths"]["output"] / "pancreas_pathway_role_trajectory_rows.csv",
    index=False,
)
expression_calibration.to_csv(
    CONFIG["paths"]["output"] / "pancreas_pathway_role_expression_supported_summary.csv",
    index=False,
)

role_share = role_counts.div(role_counts.sum(axis=1), axis=0)
taxonomy_cmap = LinearSegmentedColormap.from_list(
    "taxonomy_plum_teal",
    ["#F3EEE8", STYLE["colors"]["purple"], STYLE["colors"]["teal"]],
)
role_annotation = np.char.add(
    np.char.add(np.char.mod("%.0f%%", role_share.to_numpy() * 100), "\n"),
    np.char.add("n=", role_counts.to_numpy().astype(str)),
)
figure_calibration = calibration.set_index("pathway").reindex(["FGF", "Wnt", "Notch"])
figure_trajectories = trajectories.loc[
    trajectories["pathway"].isin(["FGF", "Wnt", "Notch"])
].copy()
agreement_upper = 1.22 * max(
    figure_calibration["recall_expected"].max(),
    figure_calibration["median_expected_role_weight"].max(),
)
fig, axes = plt.subplots(
    1,
    3,
    figsize=(9.4, 3.25),
    constrained_layout=True,
    gridspec_kw={"width_ratios": [1.75, 1.0, 1.25]},
)
sns.heatmap(
    role_share,
    ax=axes[0],
    cmap=taxonomy_cmap,
    vmin=0,
    vmax=0.55,
    annot=role_annotation,
    fmt="",
    linewidths=0.6,
    linecolor="white",
    cbar_kws={"label": "Within-pathway share", "shrink": 0.78},
)
axes[0].set_xlabel("Dominant trajectory role")
axes[0].set_ylabel("Pathway timing class")
axes[0].set_title("A  Pathway-to-role distribution", loc="left", fontweight="bold")
axes[0].tick_params(axis="x", rotation=28)
agreement_bars = axes[1].barh(
    figure_calibration.index,
    figure_calibration["recall_expected"],
    color=STYLE["colors"]["teal"],
    alpha=0.88,
)
axes[1].scatter(
    figure_calibration["median_expected_role_weight"],
    figure_calibration.index,
    marker="D",
    s=42,
    color=STYLE["colors"]["purple"],
    edgecolor="white",
    linewidth=0.6,
)
axes[1].bar_label(
    agreement_bars,
    labels=np.char.mod("%.3f", figure_calibration["recall_expected"].to_numpy()),
    padding=3,
    fontsize=6.5,
    color=STYLE["colors"]["teal"],
)
axes[1].set_xlim(0, agreement_upper)
axes[1].set_xlabel("Prior agreement (bar)\nMedian role weight (diamond)")
axes[1].set_ylabel("")
axes[1].invert_yaxis()
axes[1].set_title("B  Expected-role agreement", loc="left", fontweight="bold")
axes[2].axvline(0, color=STYLE["colors"]["grid"], linewidth=0.8)
sns.boxplot(
    data=figure_trajectories,
    x="expected_role_margin",
    y="pathway",
    order=["FGF", "Wnt", "Notch"],
    color=STYLE["colors"]["gold"],
    width=0.55,
    fliersize=1.4,
    linewidth=0.8,
    ax=axes[2],
)
axes[2].set_xlabel("Expected-role margin")
axes[2].set_ylabel("")
axes[2].set_title("C  Trajectory-level margins", loc="left", fontweight="bold")
fig.savefig(
    CONFIG["paths"]["figure"] / "Figure_S26_taxonomy_trajectory_calibration.png",
    dpi=400,
)
fig.savefig(
    CONFIG["paths"]["figure"] / "Figure_S26_taxonomy_trajectory_calibration.pdf"
)
plt.close(fig)

print(calibration.to_string(index=False))
print(expression_calibration.to_string(index=False))
