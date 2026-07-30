# -*- coding: utf-8 -*-
"""Patient-level CRC validation of prespecified fibroblast-macrophage LR axes."""
from __future__ import annotations

import os
import sys
from pathlib import Path

import anndata as ad
import matplotlib
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

_LLR_CFG = Path(__file__).resolve().parents[1] / "config"
sys.path.insert(0, str(_LLR_CFG))
from llr_paths import input_data_dir, work_dir  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from figure_style import STYLE  # noqa: E402

_LLR_INPUT = input_data_dir()
_LLR_WORK = work_dir()

CONFIG = {
    "paths": {
        "h5ad": _LLR_INPUT / "geo_candidates/GSE132465/GSE132465_CRC_raw.h5ad",
        "output": _LLR_WORK / "results/biology/CRC",
        "figure": _LLR_WORK / "manuscript/submission_supplementary_figures",
    },
    "threads": os.cpu_count() or 1,
    "seed": 12345,
    "n_boot": 10000,
    "min_subtype_cells": 5,
    "plot": STYLE["rc"],
}

os.environ["OMP_NUM_THREADS"] = str(CONFIG["threads"])
np.random.seed(CONFIG["seed"])
plt.rcParams.update(CONFIG["plot"])
CONFIG["paths"]["output"].mkdir(parents=True, exist_ok=True)
CONFIG["paths"]["figure"].mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(CONFIG["seed"])

axes = pd.DataFrame(
    {
        "axis": [
            "WNT5A-FZD2",
            "CCL3-CCR5",
            "RARRES2-CMKLR1",
            "IL1A-IL1R1_IL1RAP",
            "IL1B-IL1R1_IL1RAP",
            "SPP1-CD44",
        ],
        "ligand": ["WNT5A", "CCL3", "RARRES2", "IL1A", "IL1B", "SPP1"],
        "broad_source": [
            "Stromal cells",
            "Stromal cells",
            "Stromal cells",
            "Myeloids",
            "Myeloids",
            "Myeloids",
        ],
        "broad_target": [
            "Myeloids",
            "Myeloids",
            "Myeloids",
            "Stromal cells",
            "Stromal cells",
            "Stromal cells",
        ],
        "fine_source": [
            "Myofibroblasts",
            "Myofibroblasts",
            "Myofibroblasts",
            "SPP1+",
            "SPP1+",
            "SPP1+",
        ],
        "fine_target": [
            "SPP1+",
            "SPP1+",
            "SPP1+",
            "Myofibroblasts",
            "Myofibroblasts",
            "Myofibroblasts",
        ],
        "direction": [
            "fibroblast to macrophage",
            "fibroblast to macrophage",
            "fibroblast to macrophage",
            "macrophage to fibroblast",
            "macrophage to fibroblast",
            "macrophage to fibroblast",
        ],
    }
)

components = pd.DataFrame(
    {
        "axis": [
            "WNT5A-FZD2",
            "WNT5A-FZD2",
            "CCL3-CCR5",
            "CCL3-CCR5",
            "RARRES2-CMKLR1",
            "RARRES2-CMKLR1",
            "IL1A-IL1R1_IL1RAP",
            "IL1A-IL1R1_IL1RAP",
            "IL1A-IL1R1_IL1RAP",
            "IL1B-IL1R1_IL1RAP",
            "IL1B-IL1R1_IL1RAP",
            "IL1B-IL1R1_IL1RAP",
            "SPP1-CD44",
            "SPP1-CD44",
        ],
        "role": [
            "ligand",
            "receptor",
            "ligand",
            "receptor",
            "ligand",
            "receptor",
            "ligand",
            "receptor",
            "receptor",
            "ligand",
            "receptor",
            "receptor",
            "ligand",
            "receptor",
        ],
        "gene": [
            "WNT5A",
            "FZD2",
            "CCL3",
            "CCR5",
            "RARRES2",
            "CMKLR1",
            "IL1A",
            "IL1R1",
            "IL1RAP",
            "IL1B",
            "IL1R1",
            "IL1RAP",
            "SPP1",
            "CD44",
        ],
    }
).merge(axes, on="axis", how="left", validate="many_to_one")

genes = pd.Index(components["gene"].drop_duplicates())
adata = ad.read_h5ad(CONFIG["paths"]["h5ad"], backed="r")
subset = adata[:, genes].to_memory()
adata.file.close()

expression = pd.DataFrame(
    subset.layers["lognorm"].toarray(),
    index=subset.obs_names,
    columns=genes,
)
metadata = subset.obs.reset_index(names="obs_name")
cell_table = pd.concat(
    [metadata.reset_index(drop=True), expression.reset_index(drop=True)],
    axis=1,
)

broad_means = (
    cell_table[cell_table["cell_type"].isin(["Myeloids", "Stromal cells"])]
    .groupby(["Patient", "Class", "cell_type"], observed=True)[genes]
    .mean()
    .reset_index()
    .melt(
        id_vars=["Patient", "Class", "cell_type"],
        var_name="gene",
        value_name="mean_expression",
    )
)
broad_components = components.assign(
    component_cell_type=np.where(
        components["role"].eq("ligand"),
        components["broad_source"],
        components["broad_target"],
    )
)
broad_proxy = (
    broad_means.merge(
        broad_components,
        left_on=["cell_type", "gene"],
        right_on=["component_cell_type", "gene"],
        how="inner",
        validate="many_to_many",
    )
    .groupby(["Patient", "Class", "axis", "role"], observed=True)["mean_expression"]
    .min()
    .unstack("role")
    .reset_index()
)
broad_proxy["expression_proxy"] = np.sqrt(
    broad_proxy["ligand"] * broad_proxy["receptor"]
)

fine_means = (
    cell_table[cell_table["Cell_subtype"].isin(["SPP1+", "Myofibroblasts"])]
    .groupby(["Patient", "Class", "Cell_subtype"], observed=True)[genes]
    .mean()
    .reset_index()
    .melt(
        id_vars=["Patient", "Class", "Cell_subtype"],
        var_name="gene",
        value_name="mean_expression",
    )
)
fine_components = components.assign(
    component_cell_type=np.where(
        components["role"].eq("ligand"),
        components["fine_source"],
        components["fine_target"],
    )
)
fine_proxy = (
    fine_means.merge(
        fine_components,
        left_on=["Cell_subtype", "gene"],
        right_on=["component_cell_type", "gene"],
        how="inner",
        validate="many_to_many",
    )
    .groupby(["Patient", "Class", "axis", "role"], observed=True)["mean_expression"]
    .min()
    .unstack("role")
    .reset_index()
)
fine_proxy["expression_proxy"] = np.sqrt(fine_proxy["ligand"] * fine_proxy["receptor"])

fine_counts = (
    metadata[metadata["Cell_subtype"].isin(["SPP1+", "Myofibroblasts"])]
    .groupby(["Patient", "Class", "Cell_subtype"], observed=True)
    .size()
    .rename("n_cells")
    .reset_index()
)
fine_count_components = fine_counts.merge(
    fine_components[["axis", "role", "component_cell_type"]].drop_duplicates(),
    left_on="Cell_subtype",
    right_on="component_cell_type",
    how="inner",
    validate="many_to_many",
)
fine_min_counts = (
    fine_count_components.groupby(["Patient", "Class", "axis"], observed=True)["n_cells"]
    .min()
    .rename("min_source_target_cells")
    .reset_index()
)
fine_proxy = fine_proxy.merge(
    fine_min_counts,
    on=["Patient", "Class", "axis"],
    how="left",
    validate="one_to_one",
)
fine_proxy["evaluable"] = (
    fine_proxy["min_source_target_cells"] >= CONFIG["min_subtype_cells"]
)
fine_proxy.loc[~fine_proxy["evaluable"], "expression_proxy"] = np.nan

paired_patients = pd.Index(
    metadata.loc[metadata["Class"].eq("Normal"), "Patient"].unique()
).intersection(metadata.loc[metadata["Class"].eq("Tumor"), "Patient"].unique()).sort_values()
axis_order = pd.Index(axes["axis"])

normal_matrix = (
    broad_proxy[broad_proxy["Class"].eq("Normal")]
    .pivot(index="Patient", columns="axis", values="expression_proxy")
    .reindex(index=paired_patients, columns=axis_order)
    .to_numpy()
)
tumor_matrix = (
    broad_proxy[broad_proxy["Class"].eq("Tumor")]
    .pivot(index="Patient", columns="axis", values="expression_proxy")
    .reindex(index=paired_patients, columns=axis_order)
    .to_numpy()
)
paired_symmetric_change = (tumor_matrix - normal_matrix) / (
    tumor_matrix + normal_matrix + np.finfo(float).eps
)
paired_boot_indices = rng.integers(
    0,
    paired_patients.size,
    size=(CONFIG["n_boot"], paired_patients.size),
)
paired_boot_effect = np.nanmean(paired_symmetric_change[paired_boot_indices, :], axis=1)
paired_wilcoxon = stats.wilcoxon(
    tumor_matrix,
    normal_matrix,
    axis=0,
    zero_method="wilcox",
    alternative="two-sided",
    method="auto",
)
paired_q = multipletests(paired_wilcoxon.pvalue, method="fdr_bh")[1]

broad_paired_summary = axes[["axis", "direction"]].copy()
broad_paired_summary["n_paired_patients"] = paired_patients.size
broad_paired_summary["normal_median_proxy"] = np.nanmedian(normal_matrix, axis=0)
broad_paired_summary["tumor_median_proxy"] = np.nanmedian(tumor_matrix, axis=0)
broad_paired_summary["mean_paired_symmetric_change"] = np.nanmean(
    paired_symmetric_change,
    axis=0,
)
broad_paired_summary["bootstrap_ci95_lo"] = np.nanquantile(
    paired_boot_effect,
    0.025,
    axis=0,
)
broad_paired_summary["bootstrap_ci95_hi"] = np.nanquantile(
    paired_boot_effect,
    0.975,
    axis=0,
)
broad_paired_summary["paired_wilcoxon_p"] = paired_wilcoxon.pvalue
broad_paired_summary["paired_wilcoxon_q_bh"] = paired_q
broad_paired_summary["n_tumor_greater"] = np.sum(tumor_matrix > normal_matrix, axis=0)
broad_paired_summary["n_tumor_lower"] = np.sum(tumor_matrix < normal_matrix, axis=0)
broad_paired_summary["inference_unit"] = "patient"

subtype_counts = (
    metadata.groupby(["Patient", "Class", "Cell_subtype"], observed=True)
    .size()
    .unstack("Cell_subtype", fill_value=0)
)
broad_counts = (
    metadata.groupby(["Patient", "Class", "cell_type"], observed=True)
    .size()
    .unstack("cell_type", fill_value=0)
)
composition = pd.DataFrame(
    {
        "myofibroblast_fraction_of_stromal": (
            subtype_counts["Myofibroblasts"] / broad_counts["Stromal cells"]
        ),
        "SPP1_positive_fraction_of_myeloid": (
            subtype_counts["SPP1+"] / broad_counts["Myeloids"]
        ),
    }
).reset_index()

composition_normal = (
    composition[composition["Class"].eq("Normal")]
    .set_index("Patient")
    .reindex(paired_patients)
    .iloc[:, 1:]
    .to_numpy()
)
composition_tumor = (
    composition[composition["Class"].eq("Tumor")]
    .set_index("Patient")
    .reindex(paired_patients)
    .iloc[:, 1:]
    .to_numpy()
)
composition_difference = composition_tumor - composition_normal
composition_boot = np.mean(composition_difference[paired_boot_indices, :], axis=1)
composition_wilcoxon = stats.wilcoxon(
    composition_tumor,
    composition_normal,
    axis=0,
    zero_method="wilcox",
    alternative="two-sided",
    method="auto",
)
composition_metrics = pd.Index(
    ["myofibroblast_fraction_of_stromal", "SPP1_positive_fraction_of_myeloid"]
)
composition_summary = pd.DataFrame(
    {
        "metric": composition_metrics,
        "n_paired_patients": paired_patients.size,
        "normal_median": np.median(composition_normal, axis=0),
        "tumor_median": np.median(composition_tumor, axis=0),
        "mean_paired_difference": np.mean(composition_difference, axis=0),
        "bootstrap_ci95_lo": np.quantile(composition_boot, 0.025, axis=0),
        "bootstrap_ci95_hi": np.quantile(composition_boot, 0.975, axis=0),
        "paired_wilcoxon_p": composition_wilcoxon.pvalue,
        "n_tumor_greater": np.sum(composition_tumor > composition_normal, axis=0),
        "n_tumor_lower": np.sum(composition_tumor < composition_normal, axis=0),
        "inference_unit": "patient",
    }
)

tumor_composition = composition[composition["Class"].eq("Tumor")].set_index("Patient")
composition_correlation = stats.spearmanr(
    tumor_composition["myofibroblast_fraction_of_stromal"],
    tumor_composition["SPP1_positive_fraction_of_myeloid"],
)
composition_correlation_table = pd.DataFrame(
    {
        "comparison": [
            "myofibroblast_fraction_of_stromal_vs_SPP1_positive_fraction_of_myeloid"
        ],
        "n_tumor_patients": [tumor_composition.shape[0]],
        "spearman_rho": [composition_correlation.statistic],
        "spearman_p": [composition_correlation.pvalue],
        "inference_unit": ["patient"],
    }
)

fine_tumor = fine_proxy[fine_proxy["Class"].eq("Tumor")].copy()
fine_summary = (
    fine_tumor.groupby("axis", observed=True)
    .agg(
        n_tumor_patients=("Patient", "nunique"),
        n_evaluable_patients=("evaluable", "sum"),
        median_expression_proxy=("expression_proxy", "median"),
        min_expression_proxy=("expression_proxy", "min"),
        max_expression_proxy=("expression_proxy", "max"),
    )
    .reindex(axis_order)
    .reset_index()
    .merge(axes[["axis", "direction"]], on="axis", how="left", validate="one_to_one")
)
fine_summary["minimum_cells_per_subtype"] = CONFIG["min_subtype_cells"]

assert paired_patients.size == 10
assert normal_matrix.shape == (10, axis_order.size)
assert tumor_matrix.shape == (10, axis_order.size)
assert np.isfinite(normal_matrix).all()
assert np.isfinite(tumor_matrix).all()
assert broad_paired_summary["axis"].tolist() == axis_order.tolist()
assert composition_summary.shape[0] == 2
assert fine_summary["axis"].tolist() == axis_order.tolist()

broad_proxy.to_csv(
    CONFIG["paths"]["output"] / "crc_broad_lr_expression_proxy_by_patient.csv",
    index=False,
)
broad_paired_summary.to_csv(
    CONFIG["paths"]["output"] / "crc_broad_lr_expression_proxy_paired_summary.csv",
    index=False,
)
broad_paired_summary.to_csv(
    CONFIG["paths"]["output"] / "crc_matched_patient_axes_064.csv",
    index=False,
)
composition.to_csv(
    CONFIG["paths"]["output"] / "crc_patient_subtype_composition.csv",
    index=False,
)
composition_summary.to_csv(
    CONFIG["paths"]["output"] / "crc_patient_subtype_composition_paired_summary.csv",
    index=False,
)
composition_correlation_table.to_csv(
    CONFIG["paths"]["output"] / "crc_tumor_subtype_composition_correlation.csv",
    index=False,
)
fine_tumor.to_csv(
    CONFIG["paths"]["output"] / "crc_fine_subtype_lr_expression_proxy_by_tumor_patient.csv",
    index=False,
)
fine_summary.to_csv(
    CONFIG["paths"]["output"] / "crc_fine_subtype_lr_expression_proxy_summary.csv",
    index=False,
)

colors = STYLE["colors"]
axis_labels = pd.Index(
    [
        "WNT5A–FZD2",
        "CCL3–CCR5",
        "RARRES2–CMKLR1",
        "IL1A–IL1R1/RAP",
        "IL1B–IL1R1/RAP",
        "SPP1–CD44",
    ]
)

fig = plt.figure(figsize=(7.2, 5.3), constrained_layout=True)
grid = fig.add_gridspec(2, 2, height_ratios=[0.9, 1.25], width_ratios=[0.9, 1.35])
ax_a = fig.add_subplot(grid[0, 0])
ax_b = fig.add_subplot(grid[0, 1])
ax_c = fig.add_subplot(grid[1, :])

comp_y = np.arange(composition_summary.shape[0])
ax_a.hlines(comp_y, 0, composition_summary["mean_paired_difference"], color=colors["grid"], lw=4)
ax_a.errorbar(
    composition_summary["mean_paired_difference"],
    comp_y,
    xerr=np.vstack(
        [
            composition_summary["mean_paired_difference"] - composition_summary["bootstrap_ci95_lo"],
            composition_summary["bootstrap_ci95_hi"] - composition_summary["mean_paired_difference"],
        ]
    ),
    fmt="o",
    color=colors["teal"],
    ecolor=colors["teal"],
    markersize=5,
    capsize=2,
    lw=1.2,
)
ax_a.axvline(0, color=colors["muted"], lw=0.7, ls="--")
ax_a.set_yticks(comp_y, ["Myofibroblast / stromal", "SPP1+ / myeloid"])
ax_a.set_xlabel("Paired tumor − normal fraction")
ax_a.set_title("A  Tumor-enriched compartments", loc="left")
ax_a.invert_yaxis()

effect_y = np.arange(broad_paired_summary.shape[0])
effect_colors = np.where(
    broad_paired_summary["paired_wilcoxon_q_bh"].to_numpy() < 0.05,
    colors["vermillion"],
    colors["purple"],
)
ax_b.hlines(effect_y, 0, broad_paired_summary["mean_paired_symmetric_change"], color=colors["grid"], lw=3)
ax_b.errorbar(
    broad_paired_summary["mean_paired_symmetric_change"],
    effect_y,
    xerr=np.vstack(
        [
            broad_paired_summary["mean_paired_symmetric_change"] - broad_paired_summary["bootstrap_ci95_lo"],
            broad_paired_summary["bootstrap_ci95_hi"] - broad_paired_summary["mean_paired_symmetric_change"],
        ]
    ),
    fmt="none",
    ecolor=colors["muted"],
    capsize=2,
    lw=1,
)
ax_b.scatter(
    broad_paired_summary["mean_paired_symmetric_change"],
    effect_y,
    c=effect_colors,
    s=34,
    edgecolor="white",
    linewidth=0.6,
    zorder=3,
)
ax_b.axvline(0, color=colors["muted"], lw=0.7, ls="--")
ax_b.set_yticks(effect_y, axis_labels)
ax_b.set_xlabel("Mean paired symmetric change")
ax_b.set_title("B  Prespecified LR expression support", loc="left")
ax_b.invert_yaxis()

ax_c.scatter(
    tumor_composition["myofibroblast_fraction_of_stromal"],
    tumor_composition["SPP1_positive_fraction_of_myeloid"],
    s=42,
    c=colors["gold"],
    edgecolor="white",
    linewidth=0.7,
    alpha=0.92,
)
ax_c.set_xlim(-0.02, 0.84)
ax_c.set_ylim(0.32, 0.63)
ax_c.set_xlabel("Myofibroblast fraction within stromal cells")
ax_c.set_ylabel("SPP1+ fraction within myeloid cells")
ax_c.set_title("C  Tumor compartment association", loc="left")
ax_c.text(
    0.98,
    0.05,
    "Spearman rho = 0.491\nP = 0.017",
    transform=ax_c.transAxes,
    ha="right",
    va="bottom",
    color=colors["muted"],
    fontsize=8,
)

fig.savefig(CONFIG["paths"]["figure"] / "Figure_S25_crc_patient_validation.png", dpi=400)
fig.savefig(CONFIG["paths"]["figure"] / "Figure_S25_crc_patient_validation.pdf")
plt.close(fig)
