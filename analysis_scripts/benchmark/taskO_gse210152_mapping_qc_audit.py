# -*- coding: utf-8 -*-
"""Audit whether GSE210152 plaque-atlas non-mapping follows QC, patient, or Leiden structure."""
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
import sys

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from figure_style import STYLE

CONFIG = {
    "paths": {
        "raw": _LLR_INPUT / "geo_candidates/GSE210152/GSE210152_raw.h5ad",
        "mapping": _LLR_WORK / "results/biology/GSE210152/gse210152_atlas_cell_annotations.csv.gz",
        "out": _LLR_WORK / "results/biology/GSE210152",
        "fig": _LLR_WORK / "manuscript/submission_supplementary_figures",
    },
    "threads": os.cpu_count() or 1,
    "seed": 12345,
    "qc": {
        "min_genes_exclusive": 200,
        "min_counts_exclusive": 500,
        "max_counts_exclusive": 40000,
        "max_pct_mt_exclusive": 10,
    },
    "plot": STYLE["rc"],
}

os.environ["OMP_NUM_THREADS"] = str(CONFIG["threads"])
np.random.seed(CONFIG["seed"])
plt.rcParams.update(CONFIG["plot"])
C = STYLE["colors"]
CONFIG["paths"]["out"].mkdir(parents=True, exist_ok=True)

adata = ad.read_h5ad(CONFIG["paths"]["raw"])
counts = np.asarray(adata.X.sum(axis=1)).ravel()
detected = np.asarray((adata.X > 0).sum(axis=1)).ravel()
mt_mask = np.asarray(adata.var_names.astype(str).str.upper().str.startswith("MT-"), dtype=bool)
mt_counts = np.asarray(adata.X[:, mt_mask].sum(axis=1)).ravel()
pct_mt = np.divide(mt_counts * 100, counts, out=np.zeros_like(mt_counts, dtype=float), where=counts > 0)

mapping = pd.read_csv(CONFIG["paths"]["mapping"])
qc = pd.DataFrame(
    {
        "barcode_id": adata.obs_names.astype(str),
        "n_counts": counts,
        "n_genes": detected,
        "pct_mt": pct_mt,
    }
)
audit = mapping.merge(qc, on="barcode_id", how="left", validate="one_to_one")
audit["basic_atlas_qc_pass"] = (
    (audit["n_genes"] > CONFIG["qc"]["min_genes_exclusive"])
    & (audit["n_counts"] > CONFIG["qc"]["min_counts_exclusive"])
    & (audit["n_counts"] < CONFIG["qc"]["max_counts_exclusive"])
    & (audit["pct_mt"] < CONFIG["qc"]["max_pct_mt_exclusive"])
)
audit["mapping_group"] = np.where(audit["atlas_mapped"], "Atlas-mapped", "Not in processed atlas")
audit.to_csv(
    CONFIG["paths"]["out"] / "gse210152_cell_qc_mapping.csv.gz",
    index=False,
    compression="gzip",
)

qc_summary = audit.groupby("mapping_group", observed=True).agg(
    n_cells=("barcode_id", "size"),
    basic_qc_pass_n=("basic_atlas_qc_pass", "sum"),
    basic_qc_pass_fraction=("basic_atlas_qc_pass", "mean"),
    n_counts_median=("n_counts", "median"),
    n_genes_median=("n_genes", "median"),
    pct_mt_median=("pct_mt", "median"),
).reset_index()
qc_summary.to_csv(CONFIG["paths"]["out"] / "gse210152_mapping_qc_summary.csv", index=False)

by_patient_leiden = audit.groupby(["patient", "clusters"], observed=True).agg(
    n_cells=("barcode_id", "size"),
    n_atlas_mapped=("atlas_mapped", "sum"),
    atlas_mapping_fraction=("atlas_mapped", "mean"),
    basic_qc_pass_fraction=("basic_atlas_qc_pass", "mean"),
).reset_index()
by_patient_leiden["n_unmapped"] = by_patient_leiden["n_cells"] - by_patient_leiden["n_atlas_mapped"]
by_patient_leiden.to_csv(
    CONFIG["paths"]["out"] / "gse210152_mapping_by_patient_leiden.csv",
    index=False,
)

cluster_table = pd.crosstab(audit["clusters"], audit["atlas_mapped"])
patient_table = pd.crosstab(audit["patient"], audit["atlas_mapped"])
qc_table = pd.crosstab(audit["basic_atlas_qc_pass"], audit["atlas_mapped"])
cluster_chi2, cluster_p, _, _ = stats.chi2_contingency(cluster_table)
patient_chi2, patient_p, _, _ = stats.chi2_contingency(patient_table)
qc_odds, qc_fisher_p = stats.fisher_exact(qc_table)
cluster_v = np.sqrt(cluster_chi2 / (len(audit) * min(cluster_table.shape[0] - 1, cluster_table.shape[1] - 1)))
patient_v = np.sqrt(patient_chi2 / (len(audit) * min(patient_table.shape[0] - 1, patient_table.shape[1] - 1)))
counts_u = stats.mannwhitneyu(
    audit.loc[audit["atlas_mapped"], "n_counts"],
    audit.loc[~audit["atlas_mapped"], "n_counts"],
    alternative="two-sided",
)
genes_u = stats.mannwhitneyu(
    audit.loc[audit["atlas_mapped"], "n_genes"],
    audit.loc[~audit["atlas_mapped"], "n_genes"],
    alternative="two-sided",
)
mt_u = stats.mannwhitneyu(
    audit.loc[audit["atlas_mapped"], "pct_mt"],
    audit.loc[~audit["atlas_mapped"], "pct_mt"],
    alternative="two-sided",
)
mechanism_summary = pd.DataFrame(
    [
        {
            "local_cells": len(audit),
            "mapped_cells": int(audit["atlas_mapped"].sum()),
            "unmapped_cells": int((~audit["atlas_mapped"]).sum()),
            "mapped_basic_qc_pass_fraction": float(audit.loc[audit["atlas_mapped"], "basic_atlas_qc_pass"].mean()),
            "unmapped_basic_qc_pass_fraction": float(audit.loc[~audit["atlas_mapped"], "basic_atlas_qc_pass"].mean()),
            "cluster_mapping_cramers_v": float(cluster_v),
            "cluster_mapping_p": float(cluster_p),
            "patient_mapping_cramers_v": float(patient_v),
            "patient_mapping_p": float(patient_p),
            "basic_qc_mapping_odds_ratio": float(qc_odds),
            "basic_qc_mapping_fisher_p": float(qc_fisher_p),
            "n_counts_mannwhitney_p": float(counts_u.pvalue),
            "n_genes_mannwhitney_p": float(genes_u.pvalue),
            "pct_mt_mannwhitney_p": float(mt_u.pvalue),
            "interpretation": "Exact mapping is complete for the processed Dib atlas subset; non-membership is audited as structured missingness, not label-transfer failure.",
        }
    ]
)
mechanism_summary.to_csv(
    CONFIG["paths"]["out"] / "gse210152_mapping_missingness_summary.csv",
    index=False,
)

heatmap = by_patient_leiden.pivot(index="clusters", columns="patient", values="atlas_mapping_fraction")
mapping_cmap = LinearSegmentedColormap.from_list(
    "mapping", ["#F7F3E8", "#E2B857", "#23867B"]
)
cluster_plot = audit.groupby("clusters", observed=True).agg(
    mapping_fraction=("atlas_mapped", "mean"),
    qc_pass_fraction=("basic_atlas_qc_pass", "mean"),
    n_cells=("barcode_id", "size"),
).reset_index()

fig, axes = plt.subplots(2, 2, figsize=(10.0, 7.0), constrained_layout=True)
sns.ecdfplot(data=audit, x="n_genes", hue="mapping_group", palette=[C["unmapped"], C["mapped"]], ax=axes[0, 0])
axes[0, 0].get_legend().set_title("Atlas membership")
axes[0, 0].get_legend().get_frame().set_visible(False)
axes[0, 0].set_xlim(0, np.quantile(audit["n_genes"], 0.995))
axes[0, 0].set_title("A  Detected genes", loc="left", fontweight="bold", color=C["ink"])
axes[0, 0].set_xlabel("Detected genes per cell")

sns.ecdfplot(data=audit, x="pct_mt", hue="mapping_group", palette=[C["unmapped"], C["mapped"]], ax=axes[0, 1])
axes[0, 1].get_legend().set_title("Atlas membership")
axes[0, 1].get_legend().get_frame().set_visible(False)
axes[0, 1].get_legend().remove()
axes[0, 1].set_xlim(0, np.quantile(audit["pct_mt"], 0.995))
axes[0, 1].axvline(CONFIG["qc"]["max_pct_mt_exclusive"], color=C["muted"], linestyle="--", linewidth=1)
axes[0, 1].set_title("B  Mitochondrial fraction", loc="left", fontweight="bold", color=C["ink"])
axes[0, 1].set_xlabel("Mitochondrial counts (%)")

sns.heatmap(
    heatmap,
    vmin=0,
    vmax=1,
    cmap=mapping_cmap,
    annot=True,
    fmt=".2f",
    linewidths=0.5,
    cbar_kws={"label": "Exact-mapping fraction"},
    ax=axes[1, 0],
)
axes[1, 0].set_title("C  Mapping coverage by patient and Leiden", loc="left", fontweight="bold", color=C["ink"])
axes[1, 0].set_xlabel("Patient")
axes[1, 0].set_ylabel("Leiden cluster")

x = np.arange(len(cluster_plot))
axes[1, 1].bar(x - 0.18, cluster_plot["mapping_fraction"], width=0.36, color=C["blue"], label="Atlas mapping")
axes[1, 1].bar(x + 0.18, cluster_plot["qc_pass_fraction"], width=0.36, color=C["qc_pass"], label="Basic QC pass")
axes[1, 1].set_xticks(x, cluster_plot["clusters"].astype(str))
axes[1, 1].set_ylim(0, 1.05)
axes[1, 1].set_xlabel("Leiden cluster")
axes[1, 1].set_ylabel("Cell fraction")
axes[1, 1].set_title("D  Atlas retention versus basic QC", loc="left", fontweight="bold", color=C["ink"])
axes[1, 1].legend(loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=2)

figure_png = CONFIG["paths"]["fig"] / "Supplementary_Figure_S20_gse210152_mapping_qc.png"
figure_pdf = figure_png.with_suffix(".pdf")
fig.savefig(figure_png, bbox_inches="tight")
fig.savefig(figure_pdf, bbox_inches="tight")
plt.close(fig)

print(qc_summary.to_string(index=False))
print(mechanism_summary.to_string(index=False))
