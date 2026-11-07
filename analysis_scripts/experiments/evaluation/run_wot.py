#!/usr/bin/env python
"""WOT optimal transport on pseudo-timepoint bins (official wot.ot.OTModel API)."""
import os
import sys
import argparse
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

parser=argparse.ArgumentParser()
parser.add_argument("input_h5ad")
parser.add_argument("out_dir")
parser.add_argument("--time-column",default="")
parser.add_argument("--group-column",default="")
parser.add_argument("--metadata",default="")
parser.add_argument("--time-origin",default="observed_experimental_time")
args=parser.parse_args()

CONFIG = {
    "paths": {
        "vendor": Path(os.environ["LLR_WORK_ROOT"]).resolve() / "tools/competitors/original/wot",
    },
    "threads": int(os.environ.get("LLR_NATIVE_THREADS","4")),
    "seed": 12345,
    "plot": {},
}

os.environ["OMP_NUM_THREADS"] = str(CONFIG["threads"])
np.random.seed(CONFIG["seed"])
VENDOR = CONFIG["paths"]["vendor"]
if str(VENDOR) not in sys.path:
    sys.path.insert(0, str(VENDOR))
sys.path.insert(0,str(Path(os.environ["LLR_WORK_ROOT"])/"tools/competitors/dependencies/wot"))
import wot  # noqa: E402


def pick_group_col(adata):
    for c in (args.group_column, "cell_type", "cluster", "clusters", "celltype"):
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
    inpath, outdir = args.input_h5ad, args.out_dir
    os.makedirs(outdir, exist_ok=True)
    marker = os.path.join(outdir, "wot_done.txt")
    out_csv = os.path.join(outdir, "wot_fates_summary.csv")
    if os.path.exists(marker) and os.path.getsize(marker) > 0:
        print("CACHED")
        return

    adata = sc.read_h5ad(inpath)
    n_input_cells,n_input_genes=adata.shape
    if args.metadata:
        supplied=pd.read_csv(args.metadata,index_col=0,dtype={0:str})
        supplied.index=supplied.index.astype(str)
        assert supplied.index.is_unique and set(supplied.index)==set(adata.obs_names)
        adata.obs=adata.obs.join(supplied.reindex(adata.obs_names))
    adata.layers["wot_counts_for_hvg"]=sparse.csr_matrix(adata.layers.get("counts",adata.X),dtype=np.float64,copy=True)
    counts_values=adata.layers["wot_counts_for_hvg"].data
    assert np.isfinite(counts_values).all() and (counts_values>=0).all() and np.equal(counts_values,np.round(counts_values)).all()
    started=time.perf_counter()
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
        sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat_v3",layer="wot_counts_for_hvg")
        adata = adata[:, adata.var.highly_variable].copy()

    pd.DataFrame({"gene_id":adata.var_names}).to_csv(Path(outdir)/"wot_selected_genes.csv",index=False,encoding="utf-8")

    # WOT requires an observed experimental time axis. Do not substitute an
    # embedding or inferred pseudotime for the timepoint field.
    pt = None
    for tp_col in (args.time_column,"timepoint", "time_info", "stage", "day", "age(days)"):
        if tp_col in adata.obs.columns:
            vals = pd.to_numeric(adata.obs[tp_col].astype("string").str.extract(r"([-+]?[0-9]+(?:\.[0-9]+)?)",expand=False),errors="coerce")
            if vals.notna().sum() > 10:
                pt = vals.to_numpy(dtype=float,na_value=np.nan)
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
    adata.obs["cell_growth_rate"] = 1.0

    g = pick_group_col(adata)
    adata.obs["cell_set"] = adata.obs[g].astype(str)

    pd.DataFrame({"cell_id":adata.obs_names,"day":adata.obs["day"].to_numpy(),
       "provided_time_finite":np.isfinite(pt),"cell_set":adata.obs["cell_set"].to_numpy(),
       "input_time_column":tp_col,"time_origin":args.time_origin}).to_csv(Path(outdir)/"wot_input_binding.csv",index=False,encoding="utf-8")
    tmap_dir = Path(outdir)/"tmaps"
    tmap_dir.mkdir(parents=True,exist_ok=True)
    model = wot.ot.OTModel(
        adata,
        day_field="day",
        covariate_field="covariate",
        growth_rate_field="cell_growth_rate",
        epsilon=0.05,
    )
    model.compute_all_transport_maps(overwrite=False, output_file_format="h5ad", tmap_out=str(tmap_dir/"tmaps"))

    tmap_model = wot.tmap.TransportMapModel.from_directory(str(tmap_dir/"tmaps"))
    cell_sets = adata.obs.groupby("cell_set").groups
    days = sorted(adata.obs["day"].unique())
    day = float(np.nanmax(days))
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
    probability=pd.DataFrame(fates.X.toarray() if sparse.issparse(fates.X) else np.asarray(fates.X),
       index=fates.obs_names,columns=fates.var_names)
    assert probability.index.is_unique and np.isfinite(probability.to_numpy()).all()
    np.testing.assert_allclose(probability.sum(axis=1).to_numpy(),1.0,rtol=0,atol=1e-8)
    probability.reindex(adata.obs_names).to_csv(Path(outdir)/"wot_target_probabilities.csv.gz",index_label="cell_id",encoding="utf-8")
    fates.write_h5ad(Path(outdir)/"wot_native_fates.h5ad",compression="gzip")
    adata.obs.loc[adata.obs.day.eq(day),["cell_set","day"]].rename_axis("cell_id").to_csv(Path(outdir)/"wot_target_definitions.csv",encoding="utf-8")
    summary.to_csv(out_csv, index=False,encoding="utf-8")
    maps=list(map(sc.read_h5ad,sorted(tmap_dir.glob("*.h5ad"))))
    diagnostics=pd.DataFrame(map(lambda m:dict(n_source=m.n_obs,n_target=m.n_vars,
       n_entries=m.n_obs*m.n_vars,n_finite=int(np.isfinite(m.X).sum()),minimum=float(np.nanmin(m.X)),
       source_ids_unique=m.obs_names.is_unique,target_ids_unique=m.var_names.is_unique),maps))
    assert diagnostics.n_finite.eq(diagnostics.n_entries).all()
    diagnostics.to_csv(Path(outdir)/"wot_transport_diagnostics.csv",index=False,encoding="utf-8")
    pd.DataFrame([dict(dataset=Path(outdir).name,method="WOT",task="cell_transport_and_target_identified_fate",
       n_input_cells=n_input_cells,n_input_genes=n_input_genes,n_model_genes=adata.n_vars,n_time_observed=int(np.isfinite(pt).sum()),
       n_native_fates=len(probability),n_targets=probability.shape[1],n_timepoints=int(adata.obs.day.nunique()),
       n_transport_maps=len(maps),native_seconds=time.perf_counter()-started,result_status="completed_native",
       time_column=tp_col,time_origin=args.time_origin,target_time=day,unit_growth_multiplier=1.0,
       seed=CONFIG["seed"],native_local_PCA_random_state=58951,epsilon=0.05,gene_preprocessing="native_min_cells3_then_seurat_v3_on_original_counts_top2000",native_tolerance=1e-8,
       native_max_iter=1e7,missing_time_policy="preserved_in_full_input_binding;author_API_excludes_from_transport",
       target_definition="provided_labels_at_last_observed_time",source_file=str(Path(outdir)/"wot_native_fates.h5ad"))]).to_csv(Path(outdir)/"coverage.csv",index=False,encoding="utf-8")
    with open(marker, "w", encoding="utf-8") as fh:
        fh.write("OK\n")
    print(f"OK: {len(summary)} cells, tmaps in {tmap_dir}")


if __name__ == "__main__":
    main()
