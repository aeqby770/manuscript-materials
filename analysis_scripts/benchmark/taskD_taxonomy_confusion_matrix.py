# -*- coding: utf-8 -*-
from __future__ import annotations


import sys
from pathlib import Path
_LLR_CFG = Path(__file__).resolve().parents[1] / "config"
if str(_LLR_CFG) not in sys.path:
    sys.path.insert(0, str(_LLR_CFG))
from llr_paths import assets_dir as _llr_assets_dir  # noqa: E402
_LLR_ASSETS = _llr_assets_dir()

import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

CONFIG = {
    "paths": {
        "full": _LLR_ASSETS / "results/paper_tables/table1_main/taxonomy_validation/pancreas_taxonomy_full.csv",
        "out_dir": _LLR_ASSETS / "results/paper_tables/supplement/tables/taxonomy_pilot_calibration",
        "fig_dir": _LLR_ASSETS / "manuscript/submission_figures",
    },
    "threads": os.cpu_count() or 1,
    "seed": 12345,
    "plot": {"figure.dpi": 300, "savefig.dpi": 300, "font.size": 9},
    "expected": {"FGF": "instructive", "Wnt": "instructive", "Notch": "feedback"},
    "roles": ["instructive", "reinforcing", "feedback", "lateral", "bystander"],
}

os.environ["OMP_NUM_THREADS"] = str(CONFIG["threads"])
np.random.seed(CONFIG["seed"])
plt.rcParams.update(CONFIG["plot"])
CONFIG["paths"]["out_dir"].mkdir(parents=True, exist_ok=True)
CONFIG["paths"]["fig_dir"].mkdir(parents=True, exist_ok=True)

full = pd.read_csv(CONFIG["paths"]["full"])


def pathway(lig: str):
    lig = str(lig)
    if lig.startswith("Fgf") or lig.upper().startswith("FGF"):
        return "FGF"
    if lig.startswith("Wnt") or lig.upper().startswith("WNT"):
        return "Wnt"
    if lig in {"Dlk1", "Dll1", "Dll3", "Dll4", "Jag1", "Jag2"}:
        return "Notch"
    return None


full["pathway"] = full["ligand_genes"].map(pathway)
sub = full[full["pathway"].notna() & (full["structural_valid"] == True)].copy()
sub["expected_role"] = sub["pathway"].map(CONFIG["expected"])
sub["match"] = sub["dominant_role"] == sub["expected_role"]

cm = pd.crosstab(sub["pathway"], sub["dominant_role"])
for r in CONFIG["roles"]:
    if r not in cm.columns:
        cm[r] = 0
cm = cm[CONFIG["roles"]]
cm.to_csv(CONFIG["paths"]["out_dir"] / "pancreas_pathway_role_confusion_matrix.csv")

metrics = (
    sub.groupby("pathway")
    .agg(
        n_structurally_valid=("match", "size"),
        n_match_expected=("match", "sum"),
        recall_expected=("match", "mean"),
    )
    .reset_index()
)
metrics["expected_role"] = metrics["pathway"].map(CONFIG["expected"])
# precision among predicted expected-role rows
prec_rows = []
for pw, exp in CONFIG["expected"].items():
    pred = sub[(sub["pathway"] == pw) & (sub["dominant_role"] == exp)]
    true_pos = int((pred["match"]).sum()) if len(pred) else 0
    n_pred = len(pred)
    prec_rows.append(
        {
            "pathway": pw,
            "precision_of_expected_label": (true_pos / n_pred) if n_pred else np.nan,
            "n_predicted_as_expected": n_pred,
        }
    )
prec = pd.DataFrame(prec_rows)
metrics = metrics.merge(prec, on="pathway", how="left")
metrics["known_failure_case"] = metrics["pathway"].eq("Notch")
metrics.to_csv(CONFIG["paths"]["out_dir"] / "pancreas_pathway_role_precision_recall.csv", index=False)

# figure
fig, axes = plt.subplots(1, 2, figsize=(10, 3.8), constrained_layout=True)
sns.heatmap(cm, annot=True, fmt="d", cmap="viridis", ax=axes[0], cbar_kws={"label": "n rows"})
axes[0].set_title("Pancreas pathway × dominant role\n(structurally valid)")
axes[0].set_xlabel("Predicted dominant role")
axes[0].set_ylabel("Canonical pathway prior")
bar = metrics.set_index("pathway")["recall_expected"].reindex(["FGF", "Wnt", "Notch"])
colors = ["#1b9e77", "#1b9e77", "#d95f02"]
axes[1].bar(bar.index, bar.values, color=colors)
axes[1].axhline(0.33, color="grey", ls="--", lw=0.8)
axes[1].set_ylim(0, 1)
axes[1].set_ylabel("Recall of expected role")
axes[1].set_title("Pilot prior alignment\n(Notch = known failure)")
for i, v in enumerate(bar.values):
    axes[1].text(i, v + 0.03, f"{v:.2f}", ha="center", fontsize=8)
outp = CONFIG["paths"]["fig_dir"] / "Figure_5_taxonomy_confusion_pilot.png"
fig.savefig(outp)
fig.savefig(outp.with_suffix(".pdf"))
print("wrote", outp)
print(cm)
print(metrics.to_string(index=False))
