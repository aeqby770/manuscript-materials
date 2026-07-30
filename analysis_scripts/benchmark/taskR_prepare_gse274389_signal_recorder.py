from pathlib import Path
import os

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt


CONFIG = {
    "paths": {
        "input": Path(r"F:\biotrade\LineaLinkR_benchmark_assets\01_datasets\external_validation\GSE274389_signal_recorder"),
        "output": Path(
            r"F:\biotrade\LineaLinkR_benchmark_assets\01_datasets\canonical\GSE274389_signal_recorder_canonical.h5ad"
        ),
        "qc_table": Path(
            r"F:\biotrade\0.64LLR\results\biology\gse274389\gse274389_qc_summary.csv"
        ),
    },
    "threads": os.cpu_count(),
    "seed": 12345,
    "plot": {
        "font.family": "Arial",
        "font.size": 8,
        "axes.linewidth": 0.7,
        "figure.dpi": 120,
        "savefig.dpi": 400,
    },
    "qc": {
        "min_genes": 500,
        "min_counts": 3000,
        "max_pct_mt": 8.0,
        "max_pct_ribo": 30.0,
    },
    "analysis": {
        "target_sum": 10000,
        "n_top_genes": 3000,
        "n_pcs": 40,
        "n_neighbors": 20,
        "leiden_resolution": 0.6,
    },
}

os.environ["OMP_NUM_THREADS"] = str(CONFIG["threads"])
os.environ["MKL_NUM_THREADS"] = str(CONFIG["threads"])
os.environ["OPENBLAS_NUM_THREADS"] = str(CONFIG["threads"])
np.random.seed(CONFIG["seed"])
plt.rcParams.update(CONFIG["plot"])

input_dir = CONFIG["paths"]["input"]
output_path = CONFIG["paths"]["output"]
qc_path = CONFIG["paths"]["qc_table"]
output_path.parent.mkdir(parents=True, exist_ok=True)
qc_path.parent.mkdir(parents=True, exist_ok=True)

t96 = sc.read_10x_h5(input_dir / "GSM8448423_t96.h5", gex_only=True)
t120_g = sc.read_10x_h5(input_dir / "GSM8448424_t120_dox96-G.h5", gex_only=True)
t120_r = sc.read_10x_h5(input_dir / "GSM8448425_t120_dox96-R.h5", gex_only=True)
t96.var_names_make_unique()
t120_g.var_names_make_unique()
t120_r.var_names_make_unique()

adata = sc.concat(
    {"t96": t96, "t120_G": t120_g, "t120_R": t120_r},
    label="sample",
    index_unique="-",
    join="inner",
)
adata.obs["time_hours"] = adata.obs["sample"].map(
    {"t96": 96, "t120_G": 120, "t120_R": 120}
).astype(int)
adata.obs["wnt_recorded"] = adata.obs["sample"].map(
    {"t96": "not_measured", "t120_G": "positive", "t120_R": "negative"}
).astype(str)
adata.var["mt"] = adata.var_names.str.startswith("mt-")
adata.var["ribo"] = adata.var_names.str.startswith(("Rps", "Rpl"))
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo"], inplace=True)

qc_before = (
    adata.obs.groupby("sample", observed=True)
    .agg(
        n_cells=("sample", "size"),
        median_counts=("total_counts", "median"),
        median_genes=("n_genes_by_counts", "median"),
        median_pct_mt=("pct_counts_mt", "median"),
    )
    .reset_index()
    .assign(stage="before_qc")
)

keep = (
    (adata.obs["n_genes_by_counts"] >= CONFIG["qc"]["min_genes"])
    & (adata.obs["total_counts"] >= CONFIG["qc"]["min_counts"])
    & (adata.obs["pct_counts_mt"] < CONFIG["qc"]["max_pct_mt"])
    & (adata.obs["pct_counts_ribo"] < CONFIG["qc"]["max_pct_ribo"])
)
adata = adata[keep].copy()
sc.pp.filter_genes(adata, min_cells=3)
adata.layers["counts"] = adata.X.copy()

qc_after = (
    adata.obs.groupby("sample", observed=True)
    .agg(
        n_cells=("sample", "size"),
        median_counts=("total_counts", "median"),
        median_genes=("n_genes_by_counts", "median"),
        median_pct_mt=("pct_counts_mt", "median"),
    )
    .reset_index()
    .assign(stage="after_qc")
)
pd.concat([qc_before, qc_after], ignore_index=True).to_csv(qc_path, index=False)

sc.pp.normalize_total(adata, target_sum=CONFIG["analysis"]["target_sum"])
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(
    adata,
    n_top_genes=CONFIG["analysis"]["n_top_genes"],
    flavor="seurat",
)
synthetic = adata.var_names.str.startswith(("BFP", "IRFP", "rTTA", "Stoplight"))
adata.var.loc[synthetic, "highly_variable"] = False
sc.tl.pca(
    adata,
    n_comps=CONFIG["analysis"]["n_pcs"],
    use_highly_variable=True,
    svd_solver="arpack",
    random_state=CONFIG["seed"],
)
sc.pp.neighbors(
    adata,
    n_neighbors=CONFIG["analysis"]["n_neighbors"],
    n_pcs=CONFIG["analysis"]["n_pcs"],
    random_state=CONFIG["seed"],
)
sc.tl.umap(adata, random_state=CONFIG["seed"])
sc.tl.leiden(
    adata,
    resolution=CONFIG["analysis"]["leiden_resolution"],
    random_state=CONFIG["seed"],
    key_added="leiden",
)
adata.obs["clusters"] = adata.obs["leiden"].astype(str)
adata.obs["cluster"] = adata.obs["clusters"].astype(str)
adata.obs["cell_type"] = adata.obs["clusters"].astype(str)
adata.uns["validation_design"] = (
    "Blind LineaLinkR inference: time_hours and wnt_recorded are retained only "
    "for post hoc validation and are not used for root, graph, lineage, or LR scoring."
)
adata.write_h5ad(output_path, compression="gzip")

print(adata)
print(pd.crosstab(adata.obs["clusters"], adata.obs["sample"]))
