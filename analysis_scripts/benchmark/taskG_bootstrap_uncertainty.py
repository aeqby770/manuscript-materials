# -*- coding: utf-8 -*-
"""Task G: bootstrap 95% CI for Table1 means + top-LR score uncertainty."""
from __future__ import annotations


import os
import re
import subprocess
import sys
from pathlib import Path
_LLR_CFG = Path(__file__).resolve().parents[1] / "config"
if str(_LLR_CFG) not in sys.path:
    sys.path.insert(0, str(_LLR_CFG))
from llr_paths import assets_dir as _llr_assets_dir  # noqa: E402
_LLR_ASSETS = _llr_assets_dir()
_LLR_RSCRIPT = Path(
    os.environ.get("LLR064_RSCRIPT", r"D:\R\R-4.5.3\bin\Rscript.exe")
)
_LLR_CONFIG_FILE = _LLR_ASSETS / "pipeline/config/00_config.R"
_LLR_CONFIG_CODE = (
    "args <- commandArgs(trailingOnly = TRUE); "
    "e <- new.env(parent = globalenv()); "
    "sink(tempfile()); sys.source(args[[1L]], envir = e); sink(); "
    "cat(e$CONFIG$llr$benchmark_protocol, e$CONFIG$llr$result_revision, "
    "e$CONFIG$llr$package_version, "
    "paste(e$BENCHMARK_COMPLETE_CORE_DATASETS, collapse = '|'), "
    "paste(e$BENCHMARK_PSEUDOTIME_ALGORITHMS, collapse = '|'), "
    "paste(e$BENCHMARK_BRANCH_DATASETS, collapse = '|'), "
    "paste(e$BENCHMARK_BRANCH_TERMINAL_EXPECTED_N, collapse = '|'), sep = '\\n')"
)
_LLR_CONFIG_VALUES = subprocess.run(
    [
        str(_LLR_RSCRIPT),
        "--encoding=UTF-8",
        "-e",
        _LLR_CONFIG_CODE,
        str(_LLR_CONFIG_FILE),
    ],
    check=True,
    capture_output=True,
    text=True,
    encoding="utf-8",
).stdout.strip().splitlines()
(
    _LLR_PROTOCOL,
    _LLR_REVISION,
    _LLR_PACKAGE_VERSION,
    _LLR_CORE_DATASETS_TEXT,
    _LLR_CORE_ALGORITHMS_TEXT,
    _LLR_BRANCH_DATASETS_TEXT,
    _LLR_BRANCH_TERMINAL_EXPECTED_N_TEXT,
) = _LLR_CONFIG_VALUES
_LLR_CORE_DATASETS = _LLR_CORE_DATASETS_TEXT.split("|")
_LLR_CORE_ALGORITHMS = _LLR_CORE_ALGORITHMS_TEXT.split("|")
_LLR_BRANCH_DATASETS = _LLR_BRANCH_DATASETS_TEXT.split("|")
_LLR_TABLE_TAG = re.sub(r"[^A-Za-z0-9_]+", "_", _LLR_PROTOCOL)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from figure_style import STYLE

_LLR_BRANCH_TERMINAL_EXPECTED_N = np.fromstring(
    _LLR_BRANCH_TERMINAL_EXPECTED_N_TEXT, dtype=int, sep="|"
)

CONFIG = {
    "paths": {
        "pt_core": _LLR_ASSETS / f"results/paper_tables/table1_main/table1_common_core_values_{_LLR_TABLE_TAG}.csv",
        "ari": _LLR_ASSETS / f"results/paper_tables/table1_main/table1_lineage_{_LLR_TABLE_TAG}.csv",
        "bf": _LLR_ASSETS / f"results/paper_tables/table1_main/table1_branch_f1_{_LLR_TABLE_TAG}.csv",
        "pairwise": _LLR_ASSETS / f"results/paper_tables/table1_main/table1_pairwise_values_{_LLR_TABLE_TAG}.csv",
        "pair_table": _LLR_ASSETS / "results/paper_tables/supplement/tables/linealinkr_narrative/pancreas_top_temporal_lr_pairs.csv",
        "lineage_table": _LLR_ASSETS / "results/paper_tables/supplement/tables/linealinkr_narrative/pancreas_temporal_lr_lineage_scores.csv",
        "out_dir": _LLR_ASSETS / "results/paper_tables/supplement/tables/uncertainty",
        "fig_dir": _LLR_ASSETS / "manuscript/submission_figures",
        "supp_fig_dir": _LLR_ASSETS / "manuscript/submission_supplementary_figures",
    },
    "threads": os.cpu_count() or 1,
    "seed": 12345,
    "identity": {
        "package_version": _LLR_PACKAGE_VERSION,
        "benchmark_protocol": _LLR_PROTOCOL,
        "result_revision": _LLR_REVISION,
    },
    "benchmark": {
        "core_datasets": _LLR_CORE_DATASETS,
        "core_algorithms": _LLR_CORE_ALGORITHMS,
        "branch_datasets": _LLR_BRANCH_DATASETS,
        "branch_terminal_expected_n": dict(
            zip(_LLR_BRANCH_DATASETS, _LLR_BRANCH_TERMINAL_EXPECTED_N, strict=True)
        ),
    },
    "n_boot": 2000,
    "plot": STYLE["rc"],
}

os.environ["OMP_NUM_THREADS"] = str(CONFIG["threads"])
np.random.seed(CONFIG["seed"])
plt.rcParams.update(CONFIG["plot"])
CONFIG["paths"]["out_dir"].mkdir(parents=True, exist_ok=True)
CONFIG["paths"]["fig_dir"].mkdir(parents=True, exist_ok=True)
CONFIG["paths"]["supp_fig_dir"].mkdir(parents=True, exist_ok=True)

rng = np.random.default_rng(CONFIG["seed"])


def boot_mean_ci(x: np.ndarray, n_boot: int, rng: np.random.Generator):
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return np.nan, np.nan, np.nan, np.nan
    means = np.empty(n_boot, dtype=float)
    n = x.size
    for i in range(n_boot):
        means[i] = np.mean(x[rng.integers(0, n, size=n)])
    lo, hi = np.quantile(means, [0.025, 0.975])
    return float(np.mean(x)), float(lo), float(hi), float(np.std(means, ddof=1))


pt = pd.read_csv(CONFIG["paths"]["pt_core"])
core_datasets = CONFIG["benchmark"]["core_datasets"]
core_algorithms = CONFIG["benchmark"]["core_algorithms"]
core_keys = pt.groupby(["dataset", "algorithm"], observed=True).size()
if not (
    len(core_datasets) == 19
    and len(core_algorithms) == 12
    and len(pt) == len(core_datasets) * len(core_algorithms)
    and len(core_keys) == len(pt)
    and core_keys.eq(1).all()
    and set(pt["dataset"]) == set(core_datasets)
    and set(pt["algorithm"]) == set(core_algorithms)
    and pt["is_duplicate"].eq(False).all()
    and pt["cell_coverage"].eq(1.0).all()
    and np.isfinite(pt[["Spearman", "Pearson", "Kendall", "RMSE"]].to_numpy(dtype=float)).all()
):
    raise RuntimeError("Common-core input is not the configured complete 19-dataset x 12-algorithm matrix")
sim = pt["dataset"].str.startswith("sim_") | pt["dataset"].str.contains("simulated", case=False)
pt["scope"] = np.where(sim, "sim", "real")
ari = pd.read_csv(CONFIG["paths"]["ari"])
bf = pd.read_csv(CONFIG["paths"]["bf"])
branch_datasets = CONFIG["benchmark"]["branch_datasets"]
branch_expected_n = CONFIG["benchmark"]["branch_terminal_expected_n"]
ari_keys = ari.groupby(["dataset", "algorithm"], observed=True).size()
bf_keys = bf.groupby(["dataset", "algorithm"], observed=True).size()
ari_rows_per_algorithm = ari.groupby("algorithm", observed=True).size()
bf_rows_per_algorithm = bf.groupby("algorithm", observed=True).size()
ari_expected_n = ari["dataset"].map(branch_expected_n)
bf_expected_n = bf["dataset"].map(branch_expected_n)
if not (
    len(branch_datasets) == 4
    and len(branch_expected_n) == 4
    and set(ari["dataset"]) == set(branch_datasets)
    and set(bf["dataset"]) == set(branch_datasets)
    and set(ari["algorithm"]) == set(bf["algorithm"])
    and ari_keys.eq(1).all()
    and bf_keys.eq(1).all()
    and ari_rows_per_algorithm.eq(len(branch_datasets)).all()
    and bf_rows_per_algorithm.eq(len(branch_datasets)).all()
    and ari_expected_n.notna().all()
    and bf_expected_n.notna().all()
    and ari["n_aligned"].eq(ari_expected_n).all()
    and bf["n_aligned"].eq(bf_expected_n).all()
):
    raise RuntimeError(
        "Branch inputs do not match the configured four-dataset terminal-cell contract"
    )

rows = []
for algo, g in pt.groupby("algorithm"):
    for scope, sub in [("all", g), ("real", g[g["scope"] == "real"]), ("sim", g[g["scope"] == "sim"])]:
        mean, lo, hi, se = boot_mean_ci(sub["Spearman"].to_numpy(), CONFIG["n_boot"], rng)
        rows.append(
            {
                "algorithm": algo,
                "metric": "Spearman",
                "scope": scope,
                "n": int(np.isfinite(sub["Spearman"]).sum()),
                "mean": mean,
                "ci95_lo": lo,
                "ci95_hi": hi,
                "boot_se": se,
            }
        )

spearman_ci = pd.DataFrame(rows)
spearman_ci = spearman_ci[spearman_ci["metric"] == "Spearman"]
for scope, expected_datasets in {
    "all": set(core_datasets),
    "real": set(pt.loc[pt["scope"] == "real", "dataset"]),
    "sim": set(pt.loc[pt["scope"] == "sim", "dataset"]),
}.items():
    observed = pt.loc[pt["scope"].eq(scope) | (scope == "all")].groupby("algorithm")["dataset"].agg(set)
    if len(observed) != len(core_algorithms) or not observed.map(lambda x: x == expected_datasets).all():
        raise RuntimeError(f"Common-core {scope} bootstrap does not use one shared dataset set across algorithms")
for algo, g in ari.groupby("algorithm"):
    mean, lo, hi, se = boot_mean_ci(g["ARI"].to_numpy(), CONFIG["n_boot"], rng)
    rows.append(
        {
            "algorithm": algo,
            "metric": "lineage_ARI",
            "scope": "sim",
            "n": int(np.isfinite(g["ARI"]).sum()),
            "mean": mean,
            "ci95_lo": lo,
            "ci95_hi": hi,
            "boot_se": se,
        }
    )
for algo, g in bf.groupby("algorithm"):
    mean, lo, hi, se = boot_mean_ci(g["Branch_F1"].to_numpy(), CONFIG["n_boot"], rng)
    rows.append(
        {
            "algorithm": algo,
            "metric": "Branch_F1_one_to_one",
            "scope": "sim",
            "n": int(np.isfinite(g["Branch_F1"]).sum()),
            "mean": mean,
            "ci95_lo": lo,
            "ci95_hi": hi,
            "boot_se": se,
        }
    )
    mean, lo, hi, se = boot_mean_ci(g["Branch_F1_many_to_one_legacy"].to_numpy(), CONFIG["n_boot"], rng)
    rows.append(
        {
            "algorithm": algo,
            "metric": "Branch_F1_many_to_one_legacy",
            "scope": "sim",
            "n": int(np.isfinite(g["Branch_F1_many_to_one_legacy"]).sum()),
            "mean": mean,
            "ci95_lo": lo,
            "ci95_hi": hi,
            "boot_se": se,
        }
    )

ci_dt = pd.DataFrame(rows)
ci_dt = ci_dt.loc[ci_dt["n"] > 0].copy()

pairwise = pd.read_csv(CONFIG["paths"]["pairwise"])
pairwise_keys = pairwise.groupby(["algorithm", "dataset"], observed=True).size()
if not (
    set(pairwise["algorithm"]) == set(core_algorithms) - {"LineaLinkR"}
    and pairwise.groupby("algorithm", observed=True)["dataset"].nunique().between(19, 24).all()
    and pairwise_keys.eq(1).all()
    and pairwise["algorithm"].ne("LineaLinkR").all()
    and np.isfinite(
        pairwise[[
            "competitor_spearman",
            "linealinkr_spearman",
            "delta_linealinkr_minus_competitor",
        ]].to_numpy(dtype=float)
    ).all()
    and np.allclose(
        pairwise["delta_linealinkr_minus_competitor"],
        pairwise["linealinkr_spearman"] - pairwise["competitor_spearman"],
        rtol=0,
        atol=1e-12,
    )
):
    raise RuntimeError("Pairwise extension input is not one valid LineaLinkR-comparator intersection per dataset")
pairwise_rows = []
for algo, g in pairwise.groupby("algorithm"):
    mean, lo, hi, se = boot_mean_ci(
        g["delta_linealinkr_minus_competitor"].to_numpy(), CONFIG["n_boot"], rng
    )
    pairwise_rows.append(
        {
            "competitor": algo,
            "n_common": int(len(g)),
            "analysis_scope": "pairwise_extension_intersection",
            "mean_delta_linealinkr_minus_competitor": mean,
            "ci95_lo": lo,
            "ci95_hi": hi,
            "boot_se": se,
        }
    )
pairwise_ci = pd.DataFrame(pairwise_rows)
pairwise_order = pd.Categorical(
    pairwise_ci["competitor"],
    categories=[x for x in core_algorithms if x != "LineaLinkR"],
    ordered=True,
)
pairwise_ci = pairwise_ci.assign(_configured_order=pairwise_order).sort_values("_configured_order")
pairwise_ci = pairwise_ci.drop(columns="_configured_order")

# Top-LR uncertainty from current LineaLinkR lineage scores.
pair = pd.read_csv(CONFIG["paths"]["pair_table"])
top = pair.sort_values("temporal_mean_score", ascending=False).head(12).copy()
lineage = pd.read_csv(CONFIG["paths"]["lineage_table"])
required_lineage_columns = {
    "dataset",
    "pair_name",
    "lineage",
    "score",
    "benchmark_protocol",
    "result_revision",
}
np.testing.assert_equal(required_lineage_columns.issubset(lineage.columns), True)
np.testing.assert_array_equal(
    lineage[["benchmark_protocol", "result_revision"]].drop_duplicates().to_numpy(),
    np.array([[CONFIG["identity"]["benchmark_protocol"], CONFIG["identity"]["result_revision"]]]),
)
headline = top["pair_name"].to_numpy()
lineage_wide = (
    lineage.loc[lineage["pair_name"].isin(headline)]
    .pivot(index="pair_name", columns="lineage", values="score")
    .reindex(headline)
)
lineage_values = lineage_wide.to_numpy(dtype=float)
lineage_order = np.argsort(~np.isfinite(lineage_values), axis=1)
lineage_values = np.take_along_axis(lineage_values, lineage_order, axis=1)
n_lineages = np.isfinite(lineage_values).sum(axis=1)
np.testing.assert_equal(lineage_values.shape[0], len(headline))
np.testing.assert_equal(np.all(n_lineages > 0), True)
bootstrap_width = int(n_lineages.max())
bootstrap_index = np.floor(
    rng.random((len(headline), CONFIG["n_boot"], bootstrap_width))
    * n_lineages[:, None, None]
).astype(int)
bootstrap_values = np.take_along_axis(
    lineage_values[:, None, :], bootstrap_index, axis=2
)
bootstrap_values = np.where(
    np.arange(bootstrap_width)[None, None, :] < n_lineages[:, None, None],
    bootstrap_values,
    np.nan,
)
bootstrap_means = np.nanmean(bootstrap_values, axis=2)
lr_ci = pd.DataFrame(
    {
        "dataset": "pancreas",
        "pair_name": headline,
        "n_lineages": n_lineages,
        "mean_score": np.nanmean(lineage_values, axis=1),
        "ci95_lo": np.nanquantile(bootstrap_means, 0.025, axis=1),
        "ci95_hi": np.nanquantile(bootstrap_means, 0.975, axis=1),
        "boot_se": np.nanstd(bootstrap_means, axis=1, ddof=1),
    }
)
ci_dt.to_csv(CONFIG["paths"]["out_dir"] / "table1_bootstrap_ci.csv", index=False)
pairwise_ci.to_csv(
    CONFIG["paths"]["out_dir"] / "table1_pairwise_delta_bootstrap_ci.csv", index=False
)
lr_ci.to_csv(CONFIG["paths"]["out_dir"] / "pancreas_top_lr_bootstrap_ci.csv", index=False)

# Benchmark CI forest for selected methods and metrics.
plot_metrics = ci_dt[
    (ci_dt["algorithm"].isin(["LineaLinkR", "slingshot", "SCORPIUS", "Palantir", "scTour"]))
    & (
        ((ci_dt["metric"] == "Spearman") & (ci_dt["scope"] == "all"))
        | (ci_dt["metric"].isin(["lineage_ARI", "Branch_F1_one_to_one"]))
    )
].copy()
plot_metrics["metric_label"] = plot_metrics["metric"].replace(
    {
        "Spearman": "Spearman ρ",
        "lineage_ARI": "Lineage ARI",
        "Branch_F1_one_to_one": "One-to-one branch F1",
    }
)
plot_metrics["label"] = plot_metrics["algorithm"] + " · " + plot_metrics["metric_label"]
plot_metrics = plot_metrics.sort_values(["metric", "mean"], ascending=[True, True])

ordering_metrics = plot_metrics.loc[plot_metrics["metric"].eq("Spearman")].sort_values(
    "mean", ascending=False
)
lineage_metrics = plot_metrics.loc[plot_metrics["metric"].eq("lineage_ARI")].sort_values(
    "mean", ascending=False
)
branch_metrics = plot_metrics.loc[
    plot_metrics["metric"].eq("Branch_F1_one_to_one")
].sort_values("mean", ascending=False)

fig, axes = plt.subplots(
    1,
    3,
    figsize=(10.2, 3.25),
    sharex=True,
    constrained_layout=True,
    gridspec_kw={"width_ratios": [1.35, 1.0, 1.0]},
)
y_ordering = np.arange(len(ordering_metrics))
ordering_labels = np.char.add(
    ordering_metrics["algorithm"].astype(str).to_numpy(),
    np.char.mod("  %.2f", ordering_metrics["mean"].to_numpy()),
)
ordering_colors = np.where(
    ordering_metrics["algorithm"].eq("LineaLinkR"),
    STYLE["colors"]["purple"],
    STYLE["colors"]["muted"],
)
axes[0].hlines(
    y_ordering,
    ordering_metrics["ci95_lo"],
    ordering_metrics["ci95_hi"],
    color=STYLE["colors"]["grid"],
    lw=2.2,
)
axes[0].scatter(
    ordering_metrics["mean"],
    y_ordering,
    c=ordering_colors,
    s=34,
    zorder=3,
    edgecolors="white",
    linewidths=0.6,
)
axes[0].set_yticks(y_ordering, labels=ordering_labels)
axes[0].invert_yaxis()
axes[0].set_title(
    f"A  Pseudotime ordering\n{int(ordering_metrics['n'].min())} common datasets",
    loc="left",
)

y_lineage = np.arange(len(lineage_metrics))
lineage_labels = np.char.add(
    lineage_metrics["algorithm"].astype(str).to_numpy(),
    np.char.mod("  %.2f", lineage_metrics["mean"].to_numpy()),
)
lineage_colors = np.where(
    lineage_metrics["algorithm"].eq("LineaLinkR"),
    STYLE["colors"]["purple"],
    STYLE["colors"]["muted"],
)
axes[1].hlines(
    y_lineage,
    lineage_metrics["ci95_lo"],
    lineage_metrics["ci95_hi"],
    color=STYLE["colors"]["grid"],
    lw=2.2,
)
axes[1].scatter(
    lineage_metrics["mean"],
    y_lineage,
    c=lineage_colors,
    s=34,
    zorder=3,
    edgecolors="white",
    linewidths=0.6,
)
axes[1].set_yticks(y_lineage, labels=lineage_labels)
axes[1].invert_yaxis()
axes[1].set_title(
    f"B  Lineage recovery\n{int(lineage_metrics['n'].min())} branch datasets",
    loc="left",
)

y_branch = np.arange(len(branch_metrics))
branch_labels = np.char.add(
    branch_metrics["algorithm"].astype(str).to_numpy(),
    np.char.mod("  %.2f", branch_metrics["mean"].to_numpy()),
)
branch_colors = np.where(
    branch_metrics["algorithm"].eq("LineaLinkR"),
    STYLE["colors"]["purple"],
    STYLE["colors"]["muted"],
)
axes[2].hlines(
    y_branch,
    branch_metrics["ci95_lo"],
    branch_metrics["ci95_hi"],
    color=STYLE["colors"]["grid"],
    lw=2.2,
)
axes[2].scatter(
    branch_metrics["mean"],
    y_branch,
    c=branch_colors,
    s=34,
    zorder=3,
    edgecolors="white",
    linewidths=0.6,
)
axes[2].set_yticks(y_branch, labels=branch_labels)
axes[2].invert_yaxis()
axes[2].set_title(
    f"C  Terminal branch F1\n{int(branch_metrics['n'].min())} branch datasets",
    loc="left",
)

axes[0].axvline(0, color=STYLE["colors"]["muted"], ls="--", lw=0.8)
axes[1].axvline(0, color=STYLE["colors"]["muted"], ls="--", lw=0.8)
axes[2].axvline(0, color=STYLE["colors"]["muted"], ls="--", lw=0.8)
axes[0].set_xlim(-0.05, 1.03)
axes[0].set_xlabel("Mean (95% CI)")
axes[1].set_xlabel("Mean (95% CI)")
axes[2].set_xlabel("Mean (95% CI)")
fig.suptitle("Dataset-bootstrap uncertainty (2,000 resamples)")
outp = CONFIG["paths"]["fig_dir"] / "Supplementary_Figure_S12_table1_bootstrap_ci.png"
fig.savefig(outp)
fig.savefig(outp.with_suffix(".pdf"))
fig.savefig(CONFIG["paths"]["supp_fig_dir"] / outp.name)
fig.savefig(CONFIG["paths"]["supp_fig_dir"] / outp.with_suffix(".pdf").name)

# LR CI plot
lr_ci = lr_ci.sort_values("mean_score", ascending=False).reset_index(drop=True)
lr_left = lr_ci.iloc[:6]
lr_right = lr_ci.iloc[6:]
lr_xlim_lower = float(lr_ci["ci95_lo"].min()) - 0.015
lr_xlim_upper = float(lr_ci["ci95_hi"].max()) + 0.015

fig2, axes2 = plt.subplots(
    1,
    2,
    figsize=(9.2, 3.45),
    sharex=True,
    constrained_layout=True,
)
y_left = np.arange(len(lr_left))
left_labels = np.char.add(
    lr_left["pair_name"].astype(str).to_numpy(),
    np.char.mod("  %.2f", lr_left["mean_score"].to_numpy()),
)
axes2[0].hlines(
    y_left,
    lr_left["ci95_lo"],
    lr_left["ci95_hi"],
    color=STYLE["colors"]["grid"],
    lw=2.2,
)
axes2[0].scatter(
    lr_left["mean_score"],
    y_left,
    c=STYLE["colors"]["vermillion"],
    s=34,
    zorder=3,
    edgecolors="white",
    linewidths=0.6,
)
axes2[0].set_yticks(y_left, labels=left_labels)
axes2[0].invert_yaxis()
axes2[0].set_title("A  Higher-score programs", loc="left")

y_right = np.arange(len(lr_right))
right_labels = np.char.add(
    lr_right["pair_name"].astype(str).to_numpy(),
    np.char.mod("  %.2f", lr_right["mean_score"].to_numpy()),
)
axes2[1].hlines(
    y_right,
    lr_right["ci95_lo"],
    lr_right["ci95_hi"],
    color=STYLE["colors"]["grid"],
    lw=2.2,
)
axes2[1].scatter(
    lr_right["mean_score"],
    y_right,
    c=STYLE["colors"]["vermillion"],
    s=34,
    zorder=3,
    edgecolors="white",
    linewidths=0.6,
)
axes2[1].set_yticks(y_right, labels=right_labels)
axes2[1].invert_yaxis()
axes2[1].set_title("B  Additional programs", loc="left")
axes2[0].set_xlim(lr_xlim_lower, lr_xlim_upper)
axes2[0].set_xlabel("Mean score (95% CI)")
axes2[1].set_xlabel("Mean score (95% CI)")
fig2.suptitle("Pancreas temporal LR lineage-bootstrap uncertainty")
outp2 = CONFIG["paths"]["fig_dir"] / "Supplementary_Figure_S12b_top_lr_bootstrap_ci.png"
fig2.savefig(outp2)
fig2.savefig(outp2.with_suffix(".pdf"))
fig2.savefig(CONFIG["paths"]["supp_fig_dir"] / outp2.name)
fig2.savefig(CONFIG["paths"]["supp_fig_dir"] / outp2.with_suffix(".pdf").name)
print("wrote", outp2)

print("wrote", outp)
print(ci_dt[(ci_dt["algorithm"] == "LineaLinkR")].to_string(index=False))
print(lr_ci.head(12).to_string(index=False))
