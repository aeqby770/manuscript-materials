#!/usr/bin/env python
"""WOT optimal transport on pseudo-timepoint bins (official wot.ot.OTModel API)."""
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

CONFIG = {
    "paths": {
        "vendor": Path(os.environ["LLR_VENDOR_ROOT"]).resolve() / "wot",
    },
    "threads": os.cpu_count() or 1,
    "seed": 12345,
    "plot": {},
}

os.environ["OMP_NUM_THREADS"] = str(CONFIG["threads"])
np.random.seed(CONFIG["seed"])
VENDOR = CONFIG["paths"]["vendor"]
if str(VENDOR) not in sys.path:
    sys.path.insert(0, str(VENDOR))
import wot  # noqa: E402
import _mem_guard
_mem_guard.cap_self(55)


def pick_group_col(adata):
    for c in ("cell_type", "cluster", "clusters", "celltype"):
        if c not in adata.obs.columns:
            continue
        v = adata.obs[c].astype(str)
        valid = v.notna() & (v != "")
        if v[valid].nunique() >= 2:
            return c
    sc.pp.neighbors(adata, n_neighbors=min(15, adata.n_obs - 1))
    sc.tl.leiden(adata, resolution=0.5, key_added="wot_leiden")
    return "wot_leiden"


def main():
    if len(sys.argv) < 3:
        print("Usage: run_wot.py <h5ad> <out_dir>", file=sys.stderr)
        sys.exit(2)
    inpath, outdir = sys.argv[1], sys.argv[2]
    os.makedirs(outdir, exist_ok=True)
    marker = os.path.join(outdir, "wot_done.txt")
    out_csv = os.path.join(outdir, "wot_fates_summary.csv")
    if os.path.exists(marker) and os.path.getsize(marker) > 0:
        print("CACHED")
        return

    adata = sc.read_h5ad(inpath)
    if "lognorm" in adata.layers:
        x = adata.layers["lognorm"]
        adata.X = x.toarray() if sparse.issparse(x) else np.asarray(x)
    elif "counts" in adata.layers:
        adata.X = adata.layers["counts"]
        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)
    else:
        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)

    sc.pp.filter_genes(adata, min_cells=3)
    if adata.n_vars > 2000:
        sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat_v3")
        adata = adata[:, adata.var.highly_variable].copy()

    # WOT requires an observed experimental time axis. Do not substitute an
    # embedding or inferred pseudotime for the timepoint field.
    pt = None
    for tp_col in ("timepoint", "time_info", "stage", "day", "age(days)"):
        if tp_col in adata.obs.columns:
            vals = pd.to_numeric(adata.obs[tp_col], errors="coerce")
            if vals.notna().sum() > 10:
                pt = vals.fillna(0).values
                break
    if pt is None:
        print("FAIL: WOT requires real experimental timepoints; none found in obs "
              "(timepoint/time_info/stage/day/age(days))", file=sys.stderr)
        sys.exit(1)

    # Bin the real timepoints into discrete OT days. No artificial 5-bin cap.
    pt_series = pd.Series(np.asarray(pt, dtype=float), index=adata.obs_names)
    if pt_series.nunique() <= 20:
        adata.obs["day"] = pt_series.values
    else:
        n_bins = int(max(3, min(20, adata.n_obs // 200)))
        adata.obs["day"] = pd.qcut(pt_series, q=n_bins, labels=False, duplicates="drop").astype(float)
    adata.obs["covariate"] = "all"
    adata.obs["cell_growth_rate"] = 0.0

    g = pick_group_col(adata)
    adata.obs["cell_set"] = adata.obs[g].astype(str)

    tmap_dir = os.path.join(outdir, "tmaps")
    model = wot.ot.OTModel(
        adata,
        day_field="day",
        covariate_field="covariate",
        growth_rate_field="cell_growth_rate",
        epsilon=0.05,
    )
    model.compute_all_transport_maps(overwrite=True, output_file_format="h5ad", tmap_out=tmap_dir)

    tmap_model = wot.tmap.TransportMapModel.from_directory(tmap_dir)
    cell_sets = adata.obs.groupby("cell_set").groups
    days = sorted(adata.obs["day"].unique())
    day = float(days[len(days) // 2])
    pop = tmap_model.population_from_cell_sets(cell_sets, at_time=day)
    if len(pop) == 0:
        print("FAIL: empty population at day", day)
        sys.exit(1)
    fates = tmap_model.fates(pop)
    summary = pd.DataFrame(
        {
            "cell_id": fates.obs.index.astype(str),
            "day": fates.obs["day"].values,
            "max_fate": fates.X.max(axis=1).A1 if sparse.issparse(fates.X) else fates.X.max(axis=1),
        }
    )
    summary.to_csv(out_csv, index=False)
    with open(marker, "w", encoding="utf-8") as fh:
        fh.write("OK\n")
    print(f"OK: {len(summary)} cells, tmaps in {tmap_dir}")


if __name__ == "__main__":
    main()
