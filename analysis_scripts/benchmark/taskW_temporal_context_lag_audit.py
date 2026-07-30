# -*- coding: utf-8 -*-
"""Context-matched descriptive lag audit for six temporal-only LR axes."""
from __future__ import annotations

import os
import sys
from pathlib import Path

import anndata as ad
import matplotlib
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap

matplotlib.use("Agg")
import matplotlib.pyplot as plt

_LLR_CFG = Path(__file__).resolve().parents[1] / "config"
sys.path.insert(0, str(_LLR_CFG))
from llr_paths import assets_dir as _llr_assets_dir  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from figure_style import STYLE  # noqa: E402

_LLR_ASSETS = _llr_assets_dir()

CONFIG = {
    "paths": {
        "canonical": Path(
            r"F:\biotrade\LineaLinkR_benchmark_assets\01_datasets\canonical"
        ),
        "output": Path(r"F:\biotrade\0.64LLR\results\biology\temporal_context"),
        "figure": _LLR_ASSETS / "manuscript/submission_supplementary_figures",
    },
    "threads": os.cpu_count() or 1,
    "seed": 12345,
    "n_bins": 12,
    "max_lag": 5,
    "plot": STYLE["rc"],
}

os.environ["OMP_NUM_THREADS"] = str(CONFIG["threads"])
np.random.seed(CONFIG["seed"])
plt.rcParams.update(CONFIG["plot"])
CONFIG["paths"]["output"].mkdir(parents=True, exist_ok=True)
CONFIG["paths"]["figure"].mkdir(parents=True, exist_ok=True)

contexts = pd.read_csv(
    CONFIG["paths"]["output"] / "temporal_only_top_lineage_contexts.csv"
)
pan_weights = pd.read_csv(
    CONFIG["paths"]["output"] / "pancreas_cell_pseudotime_lineage_weights.csv"
)
dent_weights = pd.read_csv(
    CONFIG["paths"]["output"] / "dentategyrus_cell_pseudotime_lineage_weights.csv"
)

pan_genes = pd.Index(
    contexts.loc[contexts["dataset"].eq("pancreas"), ["ligand_genes", "receptor_genes"]]
    .to_numpy()
    .ravel()
).drop_duplicates()
dent_genes = pd.Index(
    contexts.loc[
        contexts["dataset"].eq("dentategyrus"), ["ligand_genes", "receptor_genes"]
    ]
    .to_numpy()
    .ravel()
).drop_duplicates()

pan_backed = ad.read_h5ad(
    CONFIG["paths"]["canonical"] / "pancreas_canonical.h5ad", backed="r"
)
pan = pan_backed[:, pan_genes].to_memory()
pan_backed.file.close()
dent_backed = ad.read_h5ad(
    CONFIG["paths"]["canonical"] / "dentategyrus_canonical.h5ad", backed="r"
)
dent = dent_backed[:, dent_genes].to_memory()
dent_backed.file.close()

pan_expression = pd.DataFrame(
    pan.layers["lognorm"].toarray(), index=pan.obs_names, columns=pan_genes
).reset_index(names="cell_id")
dent_expression = pd.DataFrame(
    dent.layers["lognorm"].toarray(), index=dent.obs_names, columns=dent_genes
).reset_index(names="cell_id")
pan_cells = (
    pan.obs[["clusters"]]
    .reset_index(names="cell_id")
    .merge(pan_weights, on="cell_id", how="inner", validate="one_to_one")
    .merge(pan_expression, on="cell_id", how="inner", validate="one_to_one")
    .assign(dataset="pancreas")
)
dent_cells = (
    dent.obs[["clusters"]]
    .reset_index(names="cell_id")
    .merge(dent_weights, on="cell_id", how="inner", validate="one_to_one")
    .merge(dent_expression, on="cell_id", how="inner", validate="one_to_one")
    .assign(dataset="dentategyrus")
)
cells = pd.concat([pan_cells, dent_cells], ignore_index=True, sort=False)

lineage_columns = pd.Index(cells.columns[cells.columns.str.match(r"^L\d+$")])
weights_long = cells.melt(
    id_vars=["dataset", "cell_id", "clusters", "pseudotime", *pan_genes, *dent_genes],
    value_vars=lineage_columns,
    var_name="lineage_name",
    value_name="lineage_weight",
)
weights_long["lineage"] = weights_long["lineage_name"].str.removeprefix("L").astype(int)

ligand_roles = contexts.assign(
    role="Ligand",
    context_cell_type=contexts["source_label"],
    gene=contexts["ligand_genes"],
)
receptor_roles = contexts.assign(
    role="Receptor",
    context_cell_type=contexts["target_label"],
    gene=contexts["receptor_genes"],
)
roles = pd.concat([ligand_roles, receptor_roles], ignore_index=True)

context_cells = weights_long.merge(
    roles[
        [
            "dataset",
            "pair_name",
            "role",
            "gene",
            "context_cell_type",
            "lineage",
            "source_label",
            "target_label",
            "peak_score_bin",
            "peak_score",
        ]
    ],
    on=["dataset", "lineage"],
    how="inner",
    validate="many_to_many",
)
context_cells = context_cells.loc[
    context_cells["clusters"].astype(str).eq(context_cells["context_cell_type"])
    & context_cells["pseudotime"].notna()
    & context_cells["lineage_weight"].gt(0)
].copy()

expression_long = cells.melt(
    id_vars=["dataset", "cell_id"],
    value_vars=pd.Index([*pan_genes, *dent_genes]).drop_duplicates(),
    var_name="gene",
    value_name="lognorm_expression",
).dropna(subset=["lognorm_expression"])
context_cells = context_cells.drop(columns=pd.Index([*pan_genes, *dent_genes]).drop_duplicates())
context_cells = context_cells.merge(
    expression_long,
    on=["dataset", "cell_id", "gene"],
    how="inner",
    validate="many_to_one",
)

group_columns = ["dataset", "pair_name", "role"]
context_cells = context_cells.sort_values([*group_columns, "pseudotime", "cell_id"])
context_cells["cumulative_weight"] = context_cells.groupby(
    group_columns, observed=True
)["lineage_weight"].cumsum()
context_cells["total_weight"] = context_cells.groupby(
    group_columns, observed=True
)["lineage_weight"].transform("sum")
context_cells["weighted_midpoint"] = (
    context_cells["cumulative_weight"] - 0.5 * context_cells["lineage_weight"]
) / context_cells["total_weight"]
context_cells["bin"] = (
    np.floor(context_cells["weighted_midpoint"] * CONFIG["n_bins"]).astype(int) + 1
).clip(1, CONFIG["n_bins"])
context_cells["weighted_expression"] = (
    context_cells["lognorm_expression"] * context_cells["lineage_weight"]
)

profiles = (
    context_cells.groupby([*group_columns, "bin"], observed=True)
    .agg(
        weighted_expression_sum=("weighted_expression", "sum"),
        lineage_weight_sum=("lineage_weight", "sum"),
        n_cells=("cell_id", "nunique"),
        mean_pseudotime=("pseudotime", "mean"),
    )
    .reset_index()
)
profiles["weighted_mean_lognorm"] = (
    profiles["weighted_expression_sum"] / profiles["lineage_weight_sum"]
)
profile_grid = roles[[*group_columns, "gene", "context_cell_type", "lineage"]].merge(
    pd.DataFrame({"bin": np.arange(1, CONFIG["n_bins"] + 1)}), how="cross"
)
profiles = profile_grid.merge(
    profiles,
    on=[*group_columns, "bin"],
    how="left",
    validate="one_to_one",
)
profiles["profile_max"] = profiles.groupby(group_columns, observed=True)[
    "weighted_mean_lognorm"
].transform("max")

peak_rows = (
    profiles.loc[profiles["profile_max"].gt(0)].dropna(
        subset=["weighted_mean_lognorm"]
    )
    .sort_values(
        [*group_columns, "weighted_mean_lognorm", "bin"],
        ascending=[True, True, True, False, True],
    )
    .drop_duplicates(group_columns, keep="last")
    [[*group_columns, "bin", "weighted_mean_lognorm", "n_cells"]]
    .rename(
        columns={
            "bin": "peak_bin",
            "weighted_mean_lognorm": "peak_weighted_mean_lognorm",
            "n_cells": "n_cells_in_peak_bin",
        }
    )
)
peak_wide = peak_rows.pivot(
    index=["dataset", "pair_name"], columns="role", values="peak_bin"
).reset_index()
peak_wide["context_peak_lag_bins"] = peak_wide["Receptor"] - peak_wide["Ligand"]
peak_wide["context_peak_lag_positive"] = peak_wide["context_peak_lag_bins"].gt(0)

ligand_profiles = profiles.loc[profiles["role"].eq("Ligand")].rename(
    columns={"bin": "ligand_bin", "weighted_mean_lognorm": "ligand_expression"}
)[["dataset", "pair_name", "ligand_bin", "ligand_expression"]]
receptor_profiles = profiles.loc[profiles["role"].eq("Receptor")].rename(
    columns={"bin": "receptor_bin", "weighted_mean_lognorm": "receptor_expression"}
)[["dataset", "pair_name", "receptor_bin", "receptor_expression"]]
lag_grid = contexts[["dataset", "pair_name"]].merge(
    pd.DataFrame({"lag": np.arange(-CONFIG["max_lag"], CONFIG["max_lag"] + 1)}),
    how="cross",
)
xcorr_rows = (
    lag_grid.merge(ligand_profiles, on=["dataset", "pair_name"], how="left")
    .merge(receptor_profiles, on=["dataset", "pair_name"], how="left")
)
xcorr_rows = xcorr_rows.loc[
    xcorr_rows["ligand_bin"].sub(xcorr_rows["receptor_bin"]).eq(xcorr_rows["lag"])
    & xcorr_rows["ligand_expression"].notna()
    & xcorr_rows["receptor_expression"].notna()
].copy()
xcorr_rows["xy"] = xcorr_rows["ligand_expression"] * xcorr_rows["receptor_expression"]
xcorr_rows["x2"] = xcorr_rows["ligand_expression"] ** 2
xcorr_rows["y2"] = xcorr_rows["receptor_expression"] ** 2
xcorr = (
    xcorr_rows.groupby(["dataset", "pair_name", "lag"], observed=True)
    .agg(
        n=("ligand_expression", "size"),
        sx=("ligand_expression", "sum"),
        sy=("receptor_expression", "sum"),
        sxy=("xy", "sum"),
        sx2=("x2", "sum"),
        sy2=("y2", "sum"),
    )
    .reset_index()
)
xcorr["numerator"] = xcorr["n"] * xcorr["sxy"] - xcorr["sx"] * xcorr["sy"]
xcorr["denominator"] = np.sqrt(
    (xcorr["n"] * xcorr["sx2"] - xcorr["sx"] ** 2)
    * (xcorr["n"] * xcorr["sy2"] - xcorr["sy"] ** 2)
)
xcorr["correlation"] = xcorr["numerator"] / xcorr["denominator"]
best_xcorr = (
    xcorr.loc[xcorr["n"].ge(3) & xcorr["correlation"].notna()]
    .sort_values(["dataset", "pair_name", "correlation", "lag"])
    .drop_duplicates(["dataset", "pair_name"], keep="last")
    [["dataset", "pair_name", "lag", "correlation", "n"]]
    .rename(
        columns={
            "lag": "context_best_xcorr_lag",
            "correlation": "context_best_xcorr",
            "n": "context_best_xcorr_n_bins",
        }
    )
)

summary = (
    contexts.merge(peak_wide, on=["dataset", "pair_name"], how="left", validate="one_to_one")
    .merge(best_xcorr, on=["dataset", "pair_name"], how="left", validate="one_to_one")
)
summary["context_receptor_later_xcorr_bins"] = -summary["context_best_xcorr_lag"]
summary["context_receptor_later_xcorr_positive"] = summary[
    "context_receptor_later_xcorr_bins"
].gt(0)
summary["context_profile_status"] = np.where(
    summary[["Ligand", "Receptor"]].notna().all(axis=1),
    "both_profiles_nonzero",
    "zero_profile_in_selected_context",
)
summary["expression_scale"] = "log-normalized"
summary["context_definition"] = "source/target cell type and selected lineage weight"
summary["binning"] = "12 approximate lineage-weighted equal-mass bins"
summary["inference"] = "descriptive; no causal or ligand-leading test"

profiles.to_csv(
    CONFIG["paths"]["output"] / "temporal_only_context_weighted_profiles_064.csv",
    index=False,
)
xcorr.to_csv(
    CONFIG["paths"]["output"] / "temporal_only_context_cross_correlation_064.csv",
    index=False,
)
summary.to_csv(
    CONFIG["paths"]["output"] / "temporal_only_lag_audit_064.csv",
    index=False,
)

pair_order = pd.Index(
    contexts.sort_values(["dataset", "peak_score"], ascending=[True, False])[
        "pair_name"
    ]
)
profile_matrix = (
    profiles.assign(row=profiles["pair_name"] + " | " + profiles["role"])
    .pivot(index="row", columns="bin", values="weighted_mean_lognorm")
    .reindex(
        pd.Index(
            np.ravel(
                np.column_stack(
                    [pair_order + " | Ligand", pair_order + " | Receptor"]
                )
            )
        )
    )
)
row_min = np.nanmin(profile_matrix.to_numpy(), axis=1, keepdims=True)
row_max = np.nanmax(profile_matrix.to_numpy(), axis=1, keepdims=True)
profile_scaled = (profile_matrix.to_numpy() - row_min) / (row_max - row_min + np.finfo(float).eps)

colors = STYLE["colors"]
cmap = LinearSegmentedColormap.from_list(
    "temporal_profile", ["#F7F3E8", "#E2B857", "#D65258", "#7A2838"]
)
fig, (ax_h, ax_l) = plt.subplots(
    1,
    2,
    figsize=(7.2, 5.15),
    gridspec_kw={"width_ratios": [1.55, 0.8]},
    constrained_layout=True,
)
image = ax_h.imshow(profile_scaled, aspect="auto", cmap=cmap, vmin=0, vmax=1)
ax_h.set_xticks(np.arange(CONFIG["n_bins"]), np.arange(1, CONFIG["n_bins"] + 1))
ax_h.set_yticks(
    np.arange(profile_matrix.shape[0]),
    profile_matrix.index.str.replace("_", "–", regex=False).str.replace(" | ", "  ", regex=False),
)
ax_h.set_xlabel("Approximate lineage-weighted pseudotime bin")
ax_h.set_title("A  Context-matched expression profiles", loc="left")
colorbar = fig.colorbar(image, ax=ax_h, fraction=0.026, pad=0.02)
colorbar.set_label("Within-row scaled expression")

lag_plot = summary.set_index("pair_name").reindex(pair_order)
y = np.arange(pair_order.size)
ax_l.axvline(0, color=colors["muted"], lw=0.7, ls="--")
ax_l.scatter(
    lag_plot["context_peak_lag_bins"],
    y - 0.12,
    marker="o",
    s=35,
    color=colors["muted"],
    label="peak-bin difference",
)
ax_l.scatter(
    lag_plot["context_receptor_later_xcorr_bins"],
    y + 0.12,
    marker="D",
    s=32,
    color=colors["teal"],
    label="cross-correlation",
)
ax_l.set_yticks(y, pair_order.str.replace("_", "–", regex=False))
ax_l.set_xticks(np.arange(-CONFIG["max_lag"], CONFIG["max_lag"] + 1, 2))
ax_l.set_xlim(-CONFIG["max_lag"] - 0.5, CONFIG["max_lag"] + 0.5)
ax_l.set_xlabel("Receptor-later shift from cross-correlation (bins)")
ax_l.set_title("B  Lag comparison", loc="left")
ax_l.invert_yaxis()
ax_l.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0))

fig.savefig(
    CONFIG["paths"]["figure"] / "Figure_S27_temporal_context_lag_audit.png",
    dpi=400,
)
fig.savefig(CONFIG["paths"]["figure"] / "Figure_S27_temporal_context_lag_audit.pdf")
plt.close(fig)

print(summary.to_string(index=False))
