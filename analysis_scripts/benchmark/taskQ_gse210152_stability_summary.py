# -*- coding: utf-8 -*-
"""Summarize GSE210152 subset and leave-one-patient-out temporal-LR rank stability."""
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

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from figure_style import STYLE

CONFIG = {
    "paths": {
        "tables": _LLR_ASSETS / "results/biology/GSE210152",
        "fig": _LLR_ASSETS / "manuscript/submission_supplementary_figures",
        "fig_mirror": _LLR_ASSETS / "manuscript/submission_figures",
    },
    "threads": os.cpu_count() or 1,
    "seed": 12345,
    "headline_n": 30,
    "rank_set_n": 15,
    "forest_n": 12,
    "plot": STYLE["rc"],
}

os.environ["OMP_NUM_THREADS"] = str(CONFIG["threads"])
np.random.seed(CONFIG["seed"])
plt.rcParams.update(CONFIG["plot"])
C = STYLE["colors"]

tables = CONFIG["paths"]["tables"]
subset = pd.read_csv(tables / "gse210152_subset_temporal_lr_stability.csv")
full = pd.read_csv(tables / "gse210152_full_temporal_lr_rank.csv")
effects = pd.read_csv(tables / "gse210152_patient_temporal_effects.csv")
meta = pd.read_csv(tables / "gse210152_patient_meta_effects.csv")
meta_random = pd.read_csv(tables / "gse210152_patient_meta_random_effects.csv")
lopo = pd.read_csv(tables / "gse210152_patient_lopo_effects.csv")

full["trajectory_key"] = (
    full["pair_idx"].astype(str)
    + "|"
    + full["source"].astype(str)
    + "|"
    + full["target"].astype(str)
    + "|"
    + full["lineage"].astype(str)
)
effects["trajectory_key"] = (
    effects["pair_idx"].astype(str)
    + "|"
    + effects["source"].astype(str)
    + "|"
    + effects["target"].astype(str)
    + "|"
    + effects["lineage"].astype(str)
)
meta["trajectory_key"] = (
    meta["pair_idx"].astype(str)
    + "|"
    + meta["source"].astype(str)
    + "|"
    + meta["target"].astype(str)
    + "|"
    + meta["lineage"].astype(str)
)
meta_random["trajectory_key"] = (
    meta_random["pair_idx"].astype(str)
    + "|"
    + meta_random["source"].astype(str)
    + "|"
    + meta_random["target"].astype(str)
    + "|"
    + meta_random["lineage"].astype(str)
)
lopo["trajectory_key"] = (
    lopo["pair_idx"].astype(str)
    + "|"
    + lopo["source"].astype(str)
    + "|"
    + lopo["target"].astype(str)
    + "|"
    + lopo["lineage"].astype(str)
)

headline = full.head(CONFIG["headline_n"]).copy()
headline_keys = pd.Index(headline["trajectory_key"])
headline_effect = headline.set_index("trajectory_key")["temporal_effect"]
headline_top = pd.Series(
    np.arange(1, len(headline) + 1) <= CONFIG["rank_set_n"],
    index=headline_keys,
)

effects_headline = effects.loc[effects["trajectory_key"].isin(headline_keys) & effects["evaluable"]].copy()
patient_pivot = effects_headline.pivot(index="trajectory_key", columns="block_label", values="temporal_effect")
patient_rho = patient_pivot.corrwith(headline_effect, axis=0, method="spearman")
patient_rank = patient_pivot.rank(axis=0, ascending=False, method="min")
patient_top = patient_rank.le(CONFIG["rank_set_n"]).fillna(False)
headline_top_patient = headline_top.reindex(patient_top.index).fillna(False).to_numpy()[:, None]
patient_intersection = (patient_top.to_numpy() & headline_top_patient).sum(axis=0)
patient_union = (patient_top.to_numpy() | headline_top_patient).sum(axis=0)
patient_jaccard = pd.Series(patient_intersection / patient_union, index=patient_top.columns)
patient_evaluable = effects_headline.groupby("block_label", observed=True).size()

meta_headline = meta.loc[meta["trajectory_key"].isin(headline_keys)].copy()
meta_effect = meta_headline.set_index("trajectory_key")["meta_effect"]
meta_rank = meta_effect.rank(ascending=False, method="min")
meta_top = meta_rank.le(CONFIG["rank_set_n"])

lopo_headline = lopo.loc[lopo["trajectory_key"].isin(headline_keys)].copy()
lopo_pivot = lopo_headline.pivot(index="trajectory_key", columns="left_out_block_label", values="meta_effect")
lopo_rho = lopo_pivot.corrwith(meta_effect, axis=0, method="spearman")
lopo_rank = lopo_pivot.rank(axis=0, ascending=False, method="min")
lopo_top = lopo_rank.le(CONFIG["rank_set_n"]).fillna(False)
meta_top_lopo = meta_top.reindex(lopo_top.index).fillna(False).to_numpy()[:, None]
lopo_intersection = (lopo_top.to_numpy() & meta_top_lopo).sum(axis=0)
lopo_union = (lopo_top.to_numpy() | meta_top_lopo).sum(axis=0)
lopo_jaccard = pd.Series(lopo_intersection / lopo_union, index=lopo_top.columns)
lopo_evaluable = lopo_headline.groupby("left_out_block_label", observed=True).size()

patient_summary = pd.concat(
    [
        patient_evaluable.rename("n_evaluable_headline"),
        patient_rho.rename("patient_rank_spearman_vs_pooled"),
        patient_jaccard.rename("patient_top15_jaccard_vs_pooled"),
        lopo_evaluable.rename("n_lopo_headline"),
        lopo_rho.rename("lopo_rank_spearman_vs_full_meta"),
        lopo_jaccard.rename("lopo_top15_jaccard_vs_full_meta"),
    ],
    axis=1,
)
patient_summary.index.name = "patient"
patient_summary.reset_index().to_csv(tables / "gse210152_patient_rank_stability.csv", index=False)

global_summary = pd.DataFrame(
    [
        {
            "patient_rescore_pool": effects["trajectory_key"].nunique(),
            "headline_trajectories": len(headline),
            "headline_with_patient_meta": len(meta_headline),
            "median_patient_rank_spearman": float(patient_rho.median()),
            "minimum_patient_rank_spearman": float(patient_rho.min()),
            "median_patient_top15_jaccard": float(patient_jaccard.median()),
            "minimum_patient_top15_jaccard": float(patient_jaccard.min()),
            "median_lopo_rank_spearman": float(lopo_rho.median()),
            "minimum_lopo_rank_spearman": float(lopo_rho.min()),
            "median_lopo_top15_jaccard": float(lopo_jaccard.median()),
            "minimum_lopo_top15_jaccard": float(lopo_jaccard.min()),
            "median_meta_heterogeneity_I2": float(meta_headline["heterogeneity_I2"].median()),
            "interpretation": "Patient and LOO stability are rank-based amplitude checks; they do not provide causal or temporal-permutation significance.",
        }
    ]
)
global_summary.to_csv(tables / "gse210152_patient_stability_summary.csv", index=False)
global_summary.to_csv(
    tables / "gse210152_patient_stability_summary_064.csv", index=False
)

subset_plot = subset.loc[subset["analysis"] != "full"].copy()
subset_labels = np.array(["Exclude\nLeiden 2/3", "Exclude\nLeiden 3", "Exact atlas\nmapped only"])
patients = patient_summary.index.to_numpy()
x_subset = np.arange(len(subset_plot))
x_patient = np.arange(len(patients))

fig, axes = plt.subplots(2, 2, figsize=(10.2, 7.2), constrained_layout=True)
axes[0, 0].bar(
    x_subset - 0.18,
    subset_plot["effect_spearman_vs_full"],
    width=0.36,
    color=C["spearman"],
    label="Effect-rank Spearman",
)
axes[0, 0].bar(
    x_subset + 0.18,
    subset_plot["top30_jaccard_vs_full"],
    width=0.36,
    color=C["jaccard"],
    label="Top-30 Jaccard",
)
axes[0, 0].set_xticks(x_subset, subset_labels)
axes[0, 0].set_ylim(0, 1.05)
axes[0, 0].set_ylabel("Stability versus full analysis")
axes[0, 0].set_title("A  Atlas-coverage sensitivity", loc="left", fontweight="bold", color=C["ink"])
axes[0, 0].legend(
    loc="lower center",
    bbox_to_anchor=(0.5, 1.14),
    ncol=2,
    frameon=False,
    fontsize=7,
    borderaxespad=0.0,
)

axes[0, 1].bar(
    x_patient - 0.18,
    patient_summary["patient_rank_spearman_vs_pooled"],
    width=0.36,
    color=C["spearman"],
    label="Rank Spearman",
)
axes[0, 1].bar(
    x_patient + 0.18,
    patient_summary["patient_top15_jaccard_vs_pooled"],
    width=0.36,
    color=C["jaccard"],
    label="Top-15 Jaccard",
)
axes[0, 1].set_xticks(x_patient, patients)
axes[0, 1].set_ylim(-0.5, 1.05)
axes[0, 1].axhline(0, color=C["muted"], linewidth=0.8)
axes[0, 1].set_ylabel("Patient versus pooled ranking")
axes[0, 1].set_title("B  Patient-rescored headline programs", loc="left", fontweight="bold", color=C["ink"])
axes[0, 1].legend(
    loc="lower center",
    bbox_to_anchor=(0.5, 1.14),
    ncol=2,
    frameon=False,
    fontsize=7,
    borderaxespad=0.0,
)

axes[1, 0].bar(
    x_patient - 0.18,
    patient_summary["lopo_rank_spearman_vs_full_meta"],
    width=0.36,
    color=C["spearman"],
    label="Rank Spearman",
)
axes[1, 0].bar(
    x_patient + 0.18,
    patient_summary["lopo_top15_jaccard_vs_full_meta"],
    width=0.36,
    color=C["jaccard"],
    label="Top-15 Jaccard",
)
axes[1, 0].set_xticks(x_patient, patients)
axes[1, 0].set_ylim(0, 1.05)
axes[1, 0].set_ylabel("Leave-one-patient-out stability")
axes[1, 0].set_title("C  LOO meta-ranking stability", loc="left", fontweight="bold", color=C["ink"])
axes[1, 0].legend(
    loc="lower center",
    bbox_to_anchor=(0.5, 1.14),
    ncol=2,
    frameon=False,
    fontsize=7,
    borderaxespad=0.0,
)

forest = meta_random.loc[
    meta_random["trajectory_key"].isin(headline_keys) & meta_random["n_samples"].ge(3)
].nlargest(CONFIG["forest_n"], "meta_effect").sort_values("meta_effect")
forest["label"] = (
    forest["ligand_genes"].astype(str)
    + "–"
    + forest["receptor_genes"].astype(str)
    + "  "
    + forest["source"].astype(str)
    + "→"
    + forest["target"].astype(str)
    + "/L"
    + forest["lineage"].astype(str)
)
y_forest = np.arange(len(forest))
axes[1, 1].errorbar(
    forest["meta_effect"],
    y_forest,
    xerr=1.96 * forest["meta_se"],
    fmt="o",
    color=C["forest"],
    ecolor=C["muted"],
    capsize=2,
    markersize=4,
)
axes[1, 1].set_yticks(y_forest, forest["label"], fontsize=7)
axes[1, 1].set_xlabel("Random-effects temporal amplitude (95% CI)")
axes[1, 1].set_title("D  Random-effects patient meta-analysis", loc="left", fontweight="bold", color=C["ink"])

fig.suptitle(
    "GSE210152 robustness across atlas coverage and patients",
    fontsize=12,
    fontweight="bold",
    color=C["ink"],
)
figure_png = CONFIG["paths"]["fig"] / "Supplementary_Figure_S21_gse210152_lr_stability.png"
figure_pdf = figure_png.with_suffix(".pdf")
fig.savefig(figure_png, bbox_inches="tight")
fig.savefig(figure_pdf, bbox_inches="tight")
fig.savefig(CONFIG["paths"]["fig_mirror"] / figure_png.name, bbox_inches="tight")
fig.savefig(CONFIG["paths"]["fig_mirror"] / figure_pdf.name, bbox_inches="tight")
plt.close(fig)

print(subset.to_string(index=False))
print(patient_summary.to_string())
print(global_summary.to_string(index=False))
