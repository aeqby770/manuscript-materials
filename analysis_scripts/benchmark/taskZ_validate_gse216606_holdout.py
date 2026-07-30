from pathlib import Path
import os

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import spearmanr
from sklearn.metrics import (
    adjusted_rand_score,
    normalized_mutual_info_score,
    roc_auc_score,
)


CONFIG = {
    "paths": {
        "truth": Path(
            r"F:\biotrade\0.64LLR\results\biology\gse216606\gse216606_holdout_truth.csv"
        ),
        "cells": Path(
            r"F:\biotrade\0.64LLR\results\biology\gse216606\gse216606_linealinkr_cells.csv"
        ),
        "run_summary": Path(
            r"F:\biotrade\0.64LLR\results\biology\gse216606\gse216606_linealinkr_run_summary.csv"
        ),
        "umap": Path(
            r"F:\biotrade\0.64LLR\cache\linealinkr_holdout\GSE216606_celltag_lsk\umap.csv"
        ),
        "metrics": Path(
            r"F:\biotrade\0.64LLR\results\biology\gse216606\gse216606_holdout_metrics.csv"
        ),
        "clone_purity": Path(
            r"F:\biotrade\0.64LLR\results\biology\gse216606\gse216606_clone_lineage_purity.csv"
        ),
        "fate_lineage": Path(
            r"F:\biotrade\0.64LLR\results\biology\gse216606\gse216606_fate_lineage_fraction.csv"
        ),
        "figure_png": Path(
            r"F:\biotrade\0.64LLR\results\biology\gse216606\gse216606_holdout_validation_diagnostic.png"
        ),
        "figure_pdf": Path(
            r"F:\biotrade\0.64LLR\results\biology\gse216606\gse216606_holdout_validation_diagnostic.pdf"
        ),
    },
    "threads": os.cpu_count(),
    "seed": 12345,
    "bootstrap": 2000,
    "expected": {
        "package_version": "0.6.4",
        "benchmark_protocol": "linealinkr064_submission",
        "result_revision": "linealinkr064_terminal_leaf_endpoint",
    },
    "plot": {
        "font.family": "Arial",
        "font.size": 8,
        "axes.titlesize": 9,
        "axes.labelsize": 8,
        "axes.linewidth": 0.7,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "legend.fontsize": 7,
        "figure.dpi": 120,
        "savefig.dpi": 500,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    },
    "colors": {
        "ink": "#172033",
        "muted": "#667085",
        "grid": "#E4E7EC",
        "paper": "#FFFFFF",
        "early": "#5F7FA3",
        "late": "#C28E2C",
        "selected": "#237A73",
        "alternative": "#B85C70",
        "llr": "#6A4C93",
        "neutral": "#A7ADB7",
    },
}

os.environ["OMP_NUM_THREADS"] = str(CONFIG["threads"])
os.environ["MKL_NUM_THREADS"] = str(CONFIG["threads"])
os.environ["OPENBLAS_NUM_THREADS"] = str(CONFIG["threads"])
np.random.seed(CONFIG["seed"])
plt.rcParams.update(CONFIG["plot"])

truth = pd.read_csv(CONFIG["paths"]["truth"])
pred = pd.read_csv(CONFIG["paths"]["cells"])
run = pd.read_csv(CONFIG["paths"]["run_summary"])
umap = pd.read_csv(CONFIG["paths"]["umap"], header=None, names=["UMAP1", "UMAP2"])

np.testing.assert_equal(len(pred), len(truth))
np.testing.assert_equal(len(umap), len(pred))
np.testing.assert_equal(run.loc[0, "celltype_column"], "auto_clusters")
np.testing.assert_equal(run.loc[0, "package_version"], CONFIG["expected"]["package_version"])
np.testing.assert_equal(run.loc[0, "benchmark_protocol"], CONFIG["expected"]["benchmark_protocol"])
np.testing.assert_equal(run.loc[0, "result_revision"], CONFIG["expected"]["result_revision"])
np.testing.assert_array_equal(pred["cell_id"].to_numpy(), truth["cell_id"].to_numpy())
cells = pred.merge(
    truth,
    on="cell_id",
    how="inner",
    sort=False,
    validate="one_to_one",
)
np.testing.assert_array_equal(cells["cell_id"].to_numpy(), pred["cell_id"].to_numpy())
cells[["UMAP1", "UMAP2"]] = umap.to_numpy()
cells["is_day5"] = cells["day"].eq("d5").astype(int)

rho_selected = spearmanr(
    cells["pseudotime"], cells["day_numeric"], nan_policy="raise"
).statistic
auc_selected = roc_auc_score(cells["is_day5"], cells["pseudotime"])
auc_selected_unsigned = max(auc_selected, 1.0 - auc_selected)

state_ari = adjusted_rand_score(
    cells["author_cell_type"], cells["inferred_lineage"]
)
state_nmi = normalized_mutual_info_score(
    cells["author_cell_type"], cells["inferred_lineage"]
)

fate_mask = (
    cells["day"].eq("d5")
    & cells["author_fate"].notna()
    & ~cells["author_fate"].isin(["no_fate_cells", "uns"])
)
fate_cells = cells.loc[fate_mask].copy()
fate_ari = adjusted_rand_score(
    fate_cells["author_fate"], fate_cells["inferred_lineage"]
)
fate_nmi = normalized_mutual_info_score(
    fate_cells["author_fate"], fate_cells["inferred_lineage"]
)

day5 = cells.loc[cells["day"].eq("d5")].copy()
clone_sizes = day5.groupby("clone_id", observed=True).size()
eligible_clones = clone_sizes.index[clone_sizes >= 2]
clone_counts = pd.crosstab(
    day5.loc[day5["clone_id"].isin(eligible_clones), "clone_id"],
    day5.loc[day5["clone_id"].isin(eligible_clones), "inferred_lineage"],
)
clone_purity = pd.DataFrame(
    {
        "clone_id": clone_counts.index,
        "n_day5_cells": clone_counts.sum(axis=1).to_numpy(),
        "dominant_lineage": clone_counts.idxmax(axis=1).to_numpy(),
        "lineage_purity": (
            clone_counts.max(axis=1) / clone_counts.sum(axis=1)
        ).to_numpy(),
    }
)
rng = np.random.default_rng(CONFIG["seed"])
bootstrap_indices = rng.integers(
    0,
    len(clone_purity),
    size=(CONFIG["bootstrap"], len(clone_purity)),
)
bootstrap_medians = np.median(
    clone_purity["lineage_purity"].to_numpy()[bootstrap_indices], axis=1
)
clone_purity_median = clone_purity["lineage_purity"].median()
clone_purity_ci = np.quantile(bootstrap_medians, [0.025, 0.975])

clone_pair_numerator = (
    clone_counts.mul(clone_counts - 1).div(2).sum(axis=1).to_numpy()
)
clone_pair_denominator = (
    clone_counts.sum(axis=1).mul(clone_counts.sum(axis=1) - 1).div(2).to_numpy()
)
within_clone_pair_agreement = clone_pair_numerator.sum() / clone_pair_denominator.sum()
global_lineage_counts = day5["inferred_lineage"].value_counts().to_numpy()
global_pair_agreement = (
    (global_lineage_counts * (global_lineage_counts - 1) / 2).sum()
    / (len(day5) * (len(day5) - 1) / 2)
)
clone_pair_excess = (
    (within_clone_pair_agreement - global_pair_agreement)
    / (1.0 - global_pair_agreement)
)
bootstrap_pair_agreement = (
    clone_pair_numerator[bootstrap_indices].sum(axis=1)
    / clone_pair_denominator[bootstrap_indices].sum(axis=1)
)
bootstrap_pair_excess = (
    (bootstrap_pair_agreement - global_pair_agreement)
    / (1.0 - global_pair_agreement)
)
clone_pair_excess_ci = np.quantile(bootstrap_pair_excess, [0.025, 0.975])

selected_root_day = truth.set_index("cell_id").loc[
    run.loc[0, "start_cell_id"], "day"
]

metrics = pd.DataFrame(
    {
        "endpoint": [
            "signed_spearman_selected",
            "signed_day5_auc_selected",
            "direction_agnostic_auc_selected",
            "state_ari",
            "state_nmi",
            "fate_ari",
            "fate_nmi",
            "median_clone_lineage_purity",
            "clone_purity_ci_low",
            "clone_purity_ci_high",
            "within_clone_pair_agreement",
            "global_pair_agreement",
            "normalized_clone_pair_excess",
            "normalized_clone_pair_excess_ci_low",
            "normalized_clone_pair_excess_ci_high",
            "root_support",
            "direction_margin",
        ],
        "value": [
            rho_selected,
            auc_selected,
            auc_selected_unsigned,
            state_ari,
            state_nmi,
            fate_ari,
            fate_nmi,
            clone_purity_median,
            clone_purity_ci[0],
            clone_purity_ci[1],
            within_clone_pair_agreement,
            global_pair_agreement,
            clone_pair_excess,
            clone_pair_excess_ci[0],
            clone_pair_excess_ci[1],
            run.loc[0, "root_support"],
            run.loc[0, "direction_margin"],
        ],
        "n": [
            len(cells),
            len(cells),
            len(cells),
            len(cells),
            len(cells),
            len(fate_cells),
            len(fate_cells),
            len(clone_purity),
            len(clone_purity),
            len(clone_purity),
            len(clone_purity),
            len(day5),
            len(clone_purity),
            len(clone_purity),
            len(clone_purity),
            len(cells),
            len(cells),
        ],
    }
)
metrics["selected_root_day"] = selected_root_day
metrics["direction_status"] = run.loc[0, "direction_status"]

fate_lineage = pd.crosstab(
    fate_cells["author_fate"],
    fate_cells["inferred_lineage"],
    normalize="index",
).reset_index()

CONFIG["paths"]["metrics"].parent.mkdir(parents=True, exist_ok=True)
CONFIG["paths"]["figure_png"].parent.mkdir(parents=True, exist_ok=True)
metrics.to_csv(CONFIG["paths"]["metrics"], index=False)
clone_purity.to_csv(CONFIG["paths"]["clone_purity"], index=False)
fate_lineage.to_csv(CONFIG["paths"]["fate_lineage"], index=False)

colors = CONFIG["colors"]
time_cmap = LinearSegmentedColormap.from_list(
    "llr_time",
    ["#3B315B", "#507A82", "#D1A65A", "#9A4C61"],
)
fig = plt.figure(figsize=(12.6, 6.7), facecolor=colors["paper"])
grid = fig.add_gridspec(
    2,
    3,
    width_ratios=[1.25, 0.95, 1.05],
    height_ratios=[1.0, 1.0],
    wspace=0.56,
    hspace=0.48,
)
axes = np.array(
    [
        fig.add_subplot(grid[0, 0]),
        fig.add_subplot(grid[0, 1]),
        fig.add_subplot(grid[0, 2]),
        fig.add_subplot(grid[1, 0]),
        fig.add_subplot(grid[1, 1]),
        fig.add_subplot(grid[1, 2]),
    ]
)

scatter = axes[0].scatter(
    cells["UMAP1"],
    cells["UMAP2"],
    c=cells["pseudotime"],
    cmap=time_cmap,
    s=1.1,
    alpha=0.78,
    linewidths=0,
    rasterized=True,
)
axes[0].set_title("Blind LineaLinkR ordering", loc="left", fontweight="bold")
axes[0].set_xlabel("UMAP 1")
axes[0].set_ylabel("UMAP 2")
fig.colorbar(scatter, ax=axes[0], fraction=0.045, pad=0.02, label="Pseudotime")

sns.violinplot(
    data=cells,
    x="day",
    y="pseudotime",
    hue="day",
    order=["d2", "d5"],
    palette=[colors["early"], colors["late"]],
    legend=False,
    inner=None,
    cut=0,
    linewidth=0.7,
    ax=axes[1],
)
sns.boxplot(
    data=cells,
    x="day",
    y="pseudotime",
    order=["d2", "d5"],
    width=0.22,
    showfliers=False,
    boxprops={"facecolor": colors["paper"], "edgecolor": colors["ink"]},
    medianprops={"color": colors["ink"], "linewidth": 1.1},
    whiskerprops={"color": colors["ink"]},
    capprops={"color": colors["ink"]},
    ax=axes[1],
)
axes[1].set_title("Withheld sampling day", loc="left", fontweight="bold")
axes[1].set_xlabel("")
axes[1].set_ylabel("Pseudotime")
axes[1].text(
    0.03,
    0.97,
    f"signed ρ = {rho_selected:.3f}\nAUC = {auc_selected:.3f}",
    transform=axes[1].transAxes,
    va="top",
    color=colors["ink"],
)

direction_y = np.array([1.0, 0.0])
direction_value = np.array([rho_selected, auc_selected])
axes[2].hlines(
    direction_y,
    xmin=np.array([-1.0, 0.0]),
    xmax=np.array([1.0, 1.0]),
    color=colors["grid"],
    linewidth=2.2,
)
axes[2].vlines(0.0, 0.78, 1.22, color=colors["ink"], linewidth=0.8)
axes[2].vlines(0.5, -0.22, 0.22, color=colors["ink"], linewidth=0.8)
axes[2].scatter(
    direction_value,
    direction_y,
    s=78,
    color=colors["selected"],
    edgecolor=colors["paper"],
    linewidth=0.8,
    zorder=3,
)
axes[2].text(rho_selected + 0.05, 1.0, f"{rho_selected:.3f}", va="center")
axes[2].text(auc_selected + 0.05, 0.0, f"{auc_selected:.3f}", va="center")
axes[2].set_yticks(direction_y, ["Signed Spearman", "Signed day-5 AUC"])
axes[2].set_xlim(-0.55, 1.0)
axes[2].set_ylim(-0.5, 1.5)
axes[2].set_title("Direction relative to sampling time", loc="left", fontweight="bold")
axes[2].set_xlabel("Signed performance")
axes[2].set_ylabel("")

state_table = pd.crosstab(
    cells["author_cell_type"], cells["inferred_lineage"], normalize="index"
)
sns.heatmap(
    state_table,
    cmap=sns.light_palette(colors["selected"], as_cmap=True),
    vmin=0,
    vmax=1,
    cbar_kws={"shrink": 0.72},
    linewidths=0.2,
    linecolor=colors["paper"],
    ax=axes[3],
)
axes[3].set_title(
    "Published cell-state alignment",
    loc="left",
    fontweight="bold",
)
axes[3].text(
    0.0,
    1.02,
    f"ARI {state_ari:.2f} · NMI {state_nmi:.2f}",
    transform=axes[3].transAxes,
    color=colors["muted"],
)
axes[3].set_xlabel("Inferred lineage")
axes[3].set_ylabel("Published cell state")

fate_table = pd.crosstab(
    fate_cells["author_fate"], fate_cells["inferred_lineage"], normalize="index"
)
sns.heatmap(
    fate_table,
    cmap=sns.light_palette(colors["alternative"], as_cmap=True),
    vmin=0,
    vmax=1,
    cbar_kws={"shrink": 0.72},
    linewidths=0.2,
    linecolor=colors["paper"],
    ax=axes[4],
)
axes[4].set_title(
    "Published fate alignment",
    loc="left",
    fontweight="bold",
)
axes[4].text(
    0.0,
    1.02,
    f"ARI {fate_ari:.2f} · NMI {fate_nmi:.2f}",
    transform=axes[4].transAxes,
    color=colors["muted"],
)
axes[4].set_xlabel("Inferred lineage")
axes[4].set_ylabel("Published day-5 fate")

sns.histplot(
    clone_purity["lineage_purity"],
    bins=np.linspace(0, 1, 21),
    color=colors["late"],
    edgecolor=colors["paper"],
    linewidth=0.4,
    ax=axes[5],
)
axes[5].axvline(
    clone_purity_median,
    color=colors["ink"],
    linewidth=1.1,
    linestyle="--",
)
axes[5].set_title("Clone-level lineage concentration", loc="left", fontweight="bold")
axes[5].set_xlabel("Dominant-lineage fraction per clone")
axes[5].set_ylabel("Clones")
axes[5].text(
    0.03,
    0.95,
    f"median {clone_purity_median:.2f}\nn={len(clone_purity):,} clones",
    transform=axes[5].transAxes,
    va="top",
    color=colors["ink"],
)

axes[0].text(-0.16, 1.08, "A", transform=axes[0].transAxes, fontsize=11, fontweight="bold", va="top", color=colors["ink"])
axes[1].text(-0.16, 1.08, "B", transform=axes[1].transAxes, fontsize=11, fontweight="bold", va="top", color=colors["ink"])
axes[2].text(-0.16, 1.08, "C", transform=axes[2].transAxes, fontsize=11, fontweight="bold", va="top", color=colors["ink"])
axes[3].text(-0.16, 1.08, "D", transform=axes[3].transAxes, fontsize=11, fontweight="bold", va="top", color=colors["ink"])
axes[4].text(-0.16, 1.08, "E", transform=axes[4].transAxes, fontsize=11, fontweight="bold", va="top", color=colors["ink"])
axes[5].text(-0.16, 1.08, "F", transform=axes[5].transAxes, fontsize=11, fontweight="bold", va="top", color=colors["ink"])
sns.despine(fig=fig)

fig.suptitle(
    "CellTag-multi validation of blind ordering and lineage structure",
    x=0.055,
    y=0.995,
    ha="left",
    fontsize=11,
    fontweight="bold",
    color=colors["ink"],
)
fig.savefig(CONFIG["paths"]["figure_png"], bbox_inches="tight", facecolor=colors["paper"])
fig.savefig(CONFIG["paths"]["figure_pdf"], bbox_inches="tight", facecolor=colors["paper"])
plt.close(fig)

print(metrics.to_string(index=False))
