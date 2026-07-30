from pathlib import Path
import os

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import statsmodels.formula.api as smf
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score


CONFIG = {
    "paths": {
        "h5ad": Path(
            r"F:\biotrade\LineaLinkR_benchmark_assets\01_datasets\canonical\GSE274389_signal_recorder_canonical.h5ad"
        ),
        "llr_cells": Path(
            r"F:\biotrade\0.64LLR\results\biology\gse274389\gse274389_linealinkr_cells.csv"
        ),
        "wnt_temporal": Path(
            r"F:\biotrade\0.64LLR\results\biology\gse274389\gse274389_wnt_temporal_scores.csv"
        ),
        "root_candidates": Path(
            r"F:\biotrade\0.64LLR\results\biology\gse274389\gse274389_root_candidates.csv"
        ),
        "summary": Path(
            r"F:\biotrade\0.64LLR\results\biology\gse274389\gse274389_validation_summary.csv"
        ),
        "bins": Path(
            r"F:\biotrade\0.64LLR\results\biology\gse274389\gse274389_lineage_bin_validation.csv"
        ),
        "root_sensitivity": Path(
            r"F:\biotrade\0.64LLR\results\biology\gse274389\gastruloid_root_sensitivity_064.csv"
        ),
        "time_order": Path(
            r"F:\biotrade\0.64LLR\results\biology\gse274389\gastruloid_time_order_064.csv"
        ),
        "wnt_recorder": Path(
            r"F:\biotrade\0.64LLR\results\biology\gse274389\gastruloid_wnt_recorder_064.csv"
        ),
        "figure_cells": Path(
            r"F:\biotrade\0.64LLR\results\biology\gse274389\gse274389_validation_cells.csv"
        ),
        "figure_png": Path(
            r"F:\biotrade\0.64LLR\results\biology\gse274389\gse274389_validation_diagnostic.png"
        ),
        "figure_pdf": Path(
            r"F:\biotrade\0.64LLR\results\biology\gse274389\gse274389_validation_diagnostic.pdf"
        ),
    },
    "threads": os.cpu_count(),
    "seed": 12345,
    "plot": {
        "font.family": "Arial",
        "font.size": 8.5,
        "axes.titlesize": 9.5,
        "axes.labelsize": 8.5,
        "xtick.labelsize": 7.5,
        "ytick.labelsize": 7.5,
        "legend.fontsize": 7.5,
        "axes.linewidth": 0.65,
        "axes.edgecolor": "#A7ADB7",
        "axes.labelcolor": "#172033",
        "xtick.color": "#667085",
        "ytick.color": "#667085",
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
        "savefig.dpi": 500,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    },
    "colors": {
        "ink": "#172033",
        "muted": "#667085",
        "plum": "#6A4C93",
        "jade": "#237A73",
        "rose": "#B85C70",
        "ochre": "#C28E2C",
        "bluegray": "#5F7FA3",
        "neutral": "#A7ADB7",
        "pale": "#F2F4F7",
    },
    "n_bins": 12,
    "wnt_targets": ["Axin2", "Lef1", "Sp5", "Tcf7", "Tcf7l2", "Dkk1", "Notum"],
}

os.environ["OMP_NUM_THREADS"] = str(CONFIG["threads"])
os.environ["MKL_NUM_THREADS"] = str(CONFIG["threads"])
os.environ["OPENBLAS_NUM_THREADS"] = str(CONFIG["threads"])
np.random.seed(CONFIG["seed"])
plt.rcParams.update(CONFIG["plot"])

adata = ad.read_h5ad(CONFIG["paths"]["h5ad"])
llr_cells = pd.read_csv(CONFIG["paths"]["llr_cells"])
wnt_temporal = pd.read_csv(CONFIG["paths"]["wnt_temporal"])
root_candidates = pd.read_csv(CONFIG["paths"]["root_candidates"])

obs = adata.obs.reset_index(names="cell")
obs[["UMAP1", "UMAP2"]] = adata.obsm["X_umap"]
obs = obs.merge(llr_cells, on="cell", how="inner", validate="one_to_one")
obs["time_120"] = (obs["time_hours"] == 120).astype(int)
obs["recorder_positive"] = (obs["wnt_recorded"] == "positive").astype(int)

wnt_targets = adata.var_names.intersection(pd.Index(CONFIG["wnt_targets"])).tolist()
sc.tl.score_genes(
    adata,
    gene_list=wnt_targets,
    score_name="wnt_target_score",
    ctrl_size=50,
    random_state=CONFIG["seed"],
)
obs = obs.merge(
    adata.obs[["wnt_target_score"]].reset_index(names="cell"),
    on="cell",
    how="left",
    validate="one_to_one",
)

obs["lineage_rank"] = obs.groupby("dominant_lineage", observed=True)["llr_pseudotime"].rank(
    method="first", pct=True
)
obs["lineage_bin"] = np.minimum(
    np.ceil(obs["lineage_rank"] * CONFIG["n_bins"]).astype(int),
    CONFIG["n_bins"],
)
wnt_by_lineage_bin = (
    wnt_temporal.groupby(["lineage", "bin"], observed=True)
    .agg(
        llr_wnt_temporal_score=("score", "mean"),
        n_wnt_pairs=("pair_idx", "nunique"),
        n_wnt_rows=("pair_idx", "size"),
    )
    .reset_index()
    .rename(columns={"lineage": "dominant_lineage", "bin": "lineage_bin"})
)
obs = obs.merge(
    wnt_by_lineage_bin,
    on=["dominant_lineage", "lineage_bin"],
    how="left",
    validate="many_to_one",
)

t120 = obs.loc[obs["time_hours"] == 120].copy()
t120["pt_iqr"] = (
    t120["llr_pseudotime"] - t120["llr_pseudotime"].median()
) / (t120["llr_pseudotime"].quantile(0.75) - t120["llr_pseudotime"].quantile(0.25))
t120["target_iqr"] = (
    t120["wnt_target_score"] - t120["wnt_target_score"].median()
) / (t120["wnt_target_score"].quantile(0.75) - t120["wnt_target_score"].quantile(0.25))
t120["temporal_iqr"] = (
    t120["llr_wnt_temporal_score"] - t120["llr_wnt_temporal_score"].median()
) / (
    t120["llr_wnt_temporal_score"].quantile(0.75)
    - t120["llr_wnt_temporal_score"].quantile(0.25)
)

cluster_recorder_counts = pd.crosstab(t120["clusters"], t120["recorder_positive"])
comparable_clusters = cluster_recorder_counts.index[
    cluster_recorder_counts.reindex(columns=[0, 1], fill_value=0).min(axis=1) >= 20
]
model_data = t120.loc[t120["clusters"].isin(comparable_clusters)].copy()
model_data["clusters_model"] = model_data["clusters"].astype(str)
model_pt = smf.logit("recorder_positive ~ pt_iqr + C(clusters_model)", data=model_data).fit(disp=False)
model_target = smf.logit("recorder_positive ~ target_iqr + C(clusters_model)", data=model_data).fit(disp=False)
model_temporal = smf.logit("recorder_positive ~ temporal_iqr + C(clusters_model)", data=model_data).fit(disp=False)

bin_summary = (
    t120.groupby(["dominant_lineage", "lineage_bin"], observed=True)
    .agg(
        n_cells=("cell", "size"),
        recorder_positive_fraction=("recorder_positive", "mean"),
        mean_llr_pseudotime=("llr_pseudotime", "mean"),
        mean_wnt_target_score=("wnt_target_score", "mean"),
        llr_wnt_temporal_score=("llr_wnt_temporal_score", "mean"),
        n_wnt_pairs=("n_wnt_pairs", "max"),
    )
    .reset_index()
)

time_rho = spearmanr(obs["llr_pseudotime"], obs["time_hours"]).statistic
time_auc = roc_auc_score(obs["time_120"], obs["llr_pseudotime"])
time_abs_rho = np.abs(time_rho)
time_direction_agnostic_auc = np.maximum(time_auc, 1.0 - time_auc)
alternative_mask = obs["llr_pseudotime_alternative"].notna()
alternative_scores = obs.loc[alternative_mask, "llr_pseudotime_alternative"]
alternative_labels = obs.loc[alternative_mask, "time_120"].to_numpy()
alternative_n_positive = np.sum(alternative_labels == 1)
alternative_n_negative = np.sum(alternative_labels == 0)
alternative_ranks = alternative_scores.rank(method="average").to_numpy()
alternative_u = (
    np.sum(alternative_ranks[alternative_labels == 1])
    - alternative_n_positive * (alternative_n_positive + 1) / 2
)
alternative_denominator = alternative_n_positive * alternative_n_negative
time_rho_alternative = obs["llr_pseudotime_alternative"].corr(
    obs["time_hours"], method="spearman"
)
time_auc_alternative = np.divide(
    float(alternative_u),
    float(alternative_denominator),
    out=np.array(np.nan),
    where=alternative_denominator > 0,
).item()
alternative_available = bool(alternative_denominator > 0)
recorder_pt_auc = roc_auc_score(t120["recorder_positive"], t120["llr_pseudotime"])
recorder_target_auc = roc_auc_score(t120["recorder_positive"], t120["wnt_target_score"])
recorder_temporal_auc = roc_auc_score(t120["recorder_positive"], t120["llr_wnt_temporal_score"])
bin_target_rho = spearmanr(
    bin_summary["recorder_positive_fraction"], bin_summary["mean_wnt_target_score"]
).statistic
bin_temporal_rho = spearmanr(
    bin_summary["recorder_positive_fraction"], bin_summary["llr_wnt_temporal_score"]
).statistic

summary = pd.DataFrame(
    {
        "endpoint": [
            "direction_agnostic_pseudotime_vs_sampling_time_abs_spearman",
            "direction_agnostic_pseudotime_discrimination_t96_vs_t120_auc",
            "signed_pseudotime_vs_sampling_time_spearman",
            "pseudotime_discrimination_t96_vs_t120_auc",
            "pseudotime_discrimination_wnt_recorder_auc",
            "wnt_target_score_discrimination_recorder_auc",
            "llr_wnt_temporal_score_discrimination_recorder_auc",
            "cluster_adjusted_recorder_or_per_pseudotime_iqr",
            "cluster_adjusted_recorder_or_per_wnt_target_iqr",
            "cluster_adjusted_recorder_or_per_llr_wnt_temporal_iqr",
            "lineage_bin_recorder_fraction_vs_wnt_target_spearman",
            "lineage_bin_recorder_fraction_vs_llr_wnt_temporal_spearman",
        ],
        "estimate": [
            time_abs_rho,
            time_direction_agnostic_auc,
            time_rho,
            time_auc,
            recorder_pt_auc,
            recorder_target_auc,
            recorder_temporal_auc,
            np.exp(model_pt.params["pt_iqr"]),
            np.exp(model_target.params["target_iqr"]),
            np.exp(model_temporal.params["temporal_iqr"]),
            bin_target_rho,
            bin_temporal_rho,
        ],
        "lower_95": [
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.exp(model_pt.conf_int().loc["pt_iqr", 0]),
            np.exp(model_target.conf_int().loc["target_iqr", 0]),
            np.exp(model_temporal.conf_int().loc["temporal_iqr", 0]),
            np.nan,
            np.nan,
        ],
        "upper_95": [
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.exp(model_pt.conf_int().loc["pt_iqr", 1]),
            np.exp(model_target.conf_int().loc["target_iqr", 1]),
            np.exp(model_temporal.conf_int().loc["temporal_iqr", 1]),
            np.nan,
            np.nan,
        ],
        "n_cells": [
            len(obs),
            len(obs),
            len(obs),
            len(obs),
            len(t120),
            len(t120),
            len(t120),
            len(model_data),
            len(model_data),
            len(model_data),
            len(t120),
            len(t120),
        ],
        "inference_scope": [
            "primary direction-agnostic ordering evaluation",
            "primary direction-agnostic ordering evaluation",
            "secondary signed blind-direction diagnostic",
            "secondary signed blind-direction diagnostic",
            "descriptive",
            "descriptive positive control",
            "descriptive",
            "cluster-adjusted cell-level logistic model",
            "cluster-adjusted cell-level logistic model",
            "cluster-adjusted cell-level logistic model",
            "descriptive lineage-bin concordance",
            "descriptive lineage-bin concordance",
        ],
    }
)
root_sensitivity = root_candidates.copy()
root_sensitivity["time_auc"] = np.where(
    root_sensitivity["root_role"].eq("selected"),
    time_auc,
    time_auc_alternative,
)
root_sensitivity["time_spearman"] = np.where(
    root_sensitivity["root_role"].eq("selected"),
    time_rho,
    time_rho_alternative,
)
root_sensitivity["direction_sensitive"] = (
    np.where(
        alternative_available,
        (time_auc - 0.5) * (time_auc_alternative - 0.5) < 0,
        np.nan,
    )
)
root_sensitivity["inference_scope"] = (
    "post hoc cell-level root sensitivity against recorded sampling time"
)
CONFIG["paths"]["summary"].parent.mkdir(parents=True, exist_ok=True)
summary.to_csv(CONFIG["paths"]["summary"], index=False)
bin_summary.to_csv(CONFIG["paths"]["bins"], index=False)
root_sensitivity.to_csv(CONFIG["paths"]["root_sensitivity"], index=False)
summary.loc[
    summary["endpoint"].str.contains("sampling_time|t96_vs_t120"), :
].to_csv(CONFIG["paths"]["time_order"], index=False)
summary.loc[
    summary["endpoint"].str.contains("recorder|wnt", case=False), :
].to_csv(CONFIG["paths"]["wnt_recorder"], index=False)
obs[
    [
        "cell",
        "UMAP1",
        "UMAP2",
        "sample",
        "time_hours",
        "wnt_recorded",
        "cell_type",
        "clusters",
        "llr_pseudotime",
        "dominant_lineage",
        "dominant_lineage_weight",
        "wnt_target_score",
        "recorder_positive",
        "lineage_bin",
        "llr_wnt_temporal_score",
    ]
].rename(
    columns={"llr_pseudotime": "linealinkr_pseudotime"}
).to_csv(CONFIG["paths"]["figure_cells"], index=False)

alternative_label = np.where(
    alternative_available,
    f"{time_auc_alternative:.2f}",
    "not produced",
).item()
colors = CONFIG["colors"]
fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.7), constrained_layout=True)

umap = adata.obsm["X_umap"]
umap_order = np.argsort(obs["llr_pseudotime"].to_numpy())
scatter = axes[0, 0].scatter(
    umap[umap_order, 0],
    umap[umap_order, 1],
    c=obs["llr_pseudotime"].to_numpy()[umap_order],
    cmap=sns.blend_palette([colors["pale"], colors["plum"]], as_cmap=True),
    s=2.6,
    linewidths=0,
    alpha=0.82,
    rasterized=True,
)
axes[0, 0].set(
    title="Blind LineaLinkR ordering",
    xlabel="UMAP 1",
    ylabel="UMAP 2",
)
fig.colorbar(scatter, ax=axes[0, 0], fraction=0.045, pad=0.02, label="Pseudotime")

sns.violinplot(
    data=obs,
    x="sample",
    y="llr_pseudotime",
    order=["t96", "t120_R", "t120_G"],
    palette=[colors["muted"], colors["jade"], colors["rose"]],
    inner=None,
    linewidth=0.7,
    cut=0,
    ax=axes[0, 1],
)
sns.boxplot(
    data=obs,
    x="sample",
    y="llr_pseudotime",
    order=["t96", "t120_R", "t120_G"],
    width=0.18,
    showfliers=False,
    boxprops={"facecolor": "white", "edgecolor": colors["ink"], "linewidth": 0.7},
    whiskerprops={"color": colors["ink"], "linewidth": 0.7},
    medianprops={"color": colors["ink"], "linewidth": 1.0},
    ax=axes[0, 1],
)
axes[0, 1].set(
    title="Withheld sampling time",
    xlabel="Sample",
    ylabel="LineaLinkR pseudotime",
)
axes[0, 1].text(
    0.03,
    0.97,
    (
        f"ordering AUC {time_direction_agnostic_auc:.2f}\n"
        f"signed direction AUC {time_auc:.2f}"
    ),
    transform=axes[0, 1].transAxes,
    va="top",
    color=colors["ink"],
)

sns.violinplot(
    data=t120,
    x="wnt_recorded",
    y="wnt_target_score",
    order=["negative", "positive"],
    palette=[colors["jade"], colors["rose"]],
    inner=None,
    linewidth=0.7,
    cut=0,
    ax=axes[1, 0],
)
sns.boxplot(
    data=t120,
    x="wnt_recorded",
    y="wnt_target_score",
    order=["negative", "positive"],
    width=0.18,
    showfliers=False,
    boxprops={"facecolor": "white", "edgecolor": colors["ink"], "linewidth": 0.7},
    medianprops={"color": colors["ink"], "linewidth": 1.0},
    whiskerprops={"color": colors["ink"], "linewidth": 0.7},
    ax=axes[1, 0],
)
axes[1, 0].set(
    title="Wnt recorder control",
    xlabel="Wnt signal recorded at 96 h",
    ylabel="Wnt-target module score at 120 h",
)
axes[1, 0].text(
    0.03,
    0.97,
    f"AUC {recorder_target_auc:.2f}",
    transform=axes[1, 0].transAxes,
    va="top",
    color=colors["ink"],
)

point_sizes = 12 + 55 * bin_summary["n_cells"] / bin_summary["n_cells"].max()
lineage_colors = sns.blend_palette(
    [colors["plum"], colors["rose"], colors["ochre"], colors["jade"]],
    n_colors=int(bin_summary["dominant_lineage"].max()),
)
axes[1, 1].scatter(
    bin_summary["llr_wnt_temporal_score"],
    bin_summary["recorder_positive_fraction"],
    s=point_sizes,
    c=np.asarray(lineage_colors)[bin_summary["dominant_lineage"].astype(int).to_numpy() - 1],
    alpha=0.82,
    linewidths=0.35,
    edgecolors="white",
)
sns.regplot(
    data=bin_summary,
    x="llr_wnt_temporal_score",
    y="recorder_positive_fraction",
    scatter=False,
    ci=None,
    color=colors["ink"],
    line_kws={"linewidth": 0.9, "alpha": 0.75},
    ax=axes[1, 1],
)
axes[1, 1].set(
    title="Temporal score and recorded Wnt history",
    xlabel="Mean LineaLinkR Wnt temporal score",
    ylabel="Recorder-positive fraction",
)
axes[1, 1].text(
    0.03,
    0.97,
    f"Spearman {bin_temporal_rho:.2f}",
    transform=axes[1, 1].transAxes,
    va="top",
    color=colors["ink"],
)

axes[0, 0].text(-0.16, 1.08, "A", transform=axes[0, 0].transAxes, fontsize=11, fontweight="bold", color=colors["ink"], va="top")
axes[0, 1].text(-0.16, 1.08, "B", transform=axes[0, 1].transAxes, fontsize=11, fontweight="bold", color=colors["ink"], va="top")
axes[1, 0].text(-0.16, 1.08, "C", transform=axes[1, 0].transAxes, fontsize=11, fontweight="bold", color=colors["ink"], va="top")
axes[1, 1].text(-0.16, 1.08, "D", transform=axes[1, 1].transAxes, fontsize=11, fontweight="bold", color=colors["ink"], va="top")
sns.despine(fig=fig)

fig.suptitle(
    "Blind ordering and Wnt-recorded signal in mouse gastruloids",
    fontsize=11,
    fontweight="bold",
    color=colors["ink"],
)
fig.savefig(CONFIG["paths"]["figure_png"], dpi=CONFIG["plot"]["savefig.dpi"], bbox_inches="tight")
fig.savefig(CONFIG["paths"]["figure_pdf"], bbox_inches="tight")
plt.close(fig)

print(summary.to_string(index=False))
print(f"Wnt targets used: {', '.join(wnt_targets)}")
