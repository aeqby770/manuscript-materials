from pathlib import Path
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc


CONFIG = {
    "paths": {
        "raw": Path(r"F:\biotrade\0.64LLR\cache\downloads\GSE216606_RAW"),
        "truth_source": Path(
            r"F:\biotrade\0.64LLR\external_methods\CellTag-multi-2023\clone_tables\hsc.rna&atac.r1&2_master_v2.csv"
        ),
        "canonical": Path(
            r"F:\biotrade\0.64LLR\cache\biology\GSE216606_celltag_lsk_canonical.h5ad"
        ),
        "truth_holdout": Path(
            r"F:\biotrade\0.64LLR\results\biology\gse216606\gse216606_holdout_truth.csv"
        ),
        "mapping_qc": Path(
            r"F:\biotrade\0.64LLR\results\biology\gse216606\gse216606_barcode_mapping_qc.csv"
        ),
        "batch_map": Path(
            r"F:\biotrade\0.64LLR\results\biology\gse216606\gse216606_batch_mapping.csv"
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
}

os.environ["OMP_NUM_THREADS"] = str(CONFIG["threads"])
os.environ["MKL_NUM_THREADS"] = str(CONFIG["threads"])
os.environ["OPENBLAS_NUM_THREADS"] = str(CONFIG["threads"])
np.random.seed(CONFIG["seed"])
plt.rcParams.update(CONFIG["plot"])

raw_dir = CONFIG["paths"]["raw"]
canonical_path = CONFIG["paths"]["canonical"]
truth_path = CONFIG["paths"]["truth_holdout"]
mapping_qc_path = CONFIG["paths"]["mapping_qc"]
batch_map_path = CONFIG["paths"]["batch_map"]
canonical_path.parent.mkdir(parents=True, exist_ok=True)
truth_path.parent.mkdir(parents=True, exist_ok=True)

d2 = sc.read_10x_h5(
    raw_dir / "GSM6681127_d2_5_filtered_feature_bc_matrix.h5", gex_only=True
)
d51 = sc.read_10x_h5(
    raw_dir / "GSM6681128_d5_1_filtered_feature_bc_matrix.h5", gex_only=True
)
d52 = sc.read_10x_h5(
    raw_dir / "GSM6681129_d5_2_filtered_feature_bc_matrix.h5", gex_only=True
)
d53 = sc.read_10x_h5(
    raw_dir / "GSM6681130_d5_3_filtered_feature_bc_matrix.h5", gex_only=True
)
d54 = sc.read_10x_h5(
    raw_dir / "GSM6681131_d5_4_filtered_feature_bc_matrix.h5", gex_only=True
)
d55 = sc.read_10x_h5(
    raw_dir / "GSM6681132_d5_5_filtered_feature_bc_matrix.h5", gex_only=True
)
d56 = sc.read_10x_h5(
    raw_dir / "GSM6681133_d5_6_filtered_feature_bc_matrix.h5", gex_only=True
)
d57 = sc.read_10x_h5(
    raw_dir / "GSM6681134_d5_7_filtered_feature_bc_matrix.h5", gex_only=True
)
d58 = sc.read_10x_h5(
    raw_dir / "GSM6681135_d5_8_filtered_feature_bc_matrix.h5", gex_only=True
)

d2.var_names_make_unique()
d51.var_names_make_unique()
d52.var_names_make_unique()
d53.var_names_make_unique()
d54.var_names_make_unique()
d55.var_names_make_unique()
d56.var_names_make_unique()
d57.var_names_make_unique()
d58.var_names_make_unique()

d2.obs_names = "d2_5-rna-" + d2.obs_names
d51.obs_names = "d5r1-rna-" + d51.obs_names
d52.obs_names = "d5r1-rna-" + d52.obs_names.str.replace(
    r"-1$", "-2", regex=True
)
d53.obs_names = "d5r1-rna-" + d53.obs_names.str.replace(
    r"-1$", "-3", regex=True
)
d54.obs_names = "d5r1-rna-" + d54.obs_names.str.replace(
    r"-1$", "-4", regex=True
)
d55.obs_names = "d5r2-rna-" + d55.obs_names
d56.obs_names = "d5r2-rna-" + d56.obs_names.str.replace(
    r"-1$", "-2", regex=True
)
d57.obs_names = "d5r2-rna-" + d57.obs_names.str.replace(
    r"-1$", "-3", regex=True
)
d58.obs_names = "d5r2-rna-" + d58.obs_names.str.replace(
    r"-1$", "-4", regex=True
)

d2.obs["technical_batch"] = "batch_01"
d51.obs["technical_batch"] = "batch_02"
d52.obs["technical_batch"] = "batch_03"
d53.obs["technical_batch"] = "batch_04"
d54.obs["technical_batch"] = "batch_05"
d55.obs["technical_batch"] = "batch_06"
d56.obs["technical_batch"] = "batch_07"
d57.obs["technical_batch"] = "batch_08"
d58.obs["technical_batch"] = "batch_09"

adata_all = sc.concat(
    [d2, d51, d52, d53, d54, d55, d56, d57, d58],
    join="inner",
    merge="same",
)
truth_source = pd.read_csv(CONFIG["paths"]["truth_source"])
truth_rna = truth_source.loc[truth_source["assay"].eq("rna")].copy()
truth_rna = truth_rna.drop_duplicates("cell.bc").set_index("cell.bc", drop=False)

mapped_mask = adata_all.obs_names.isin(truth_rna.index)
mapping_qc = (
    adata_all.obs.assign(mapped_to_author_truth=mapped_mask)
    .groupby("technical_batch", observed=True)
    .agg(
        n_geo_filtered_cells=("technical_batch", "size"),
        n_author_truth_matched=("mapped_to_author_truth", "sum"),
    )
    .reset_index()
)
mapping_qc["mapping_fraction"] = (
    mapping_qc["n_author_truth_matched"] / mapping_qc["n_geo_filtered_cells"]
)
mapping_qc.to_csv(mapping_qc_path, index=False)

adata = adata_all[mapped_mask].copy()
annotation = truth_rna.loc[adata.obs_names].copy()
original_barcodes = adata.obs_names.astype(str)
neutral_ids = pd.Index(
    "cell_" + pd.RangeIndex(1, adata.n_obs + 1).astype(str).str.zfill(8),
    name="cell_id",
)

adata.obs_names = neutral_ids
adata.obs = pd.DataFrame(
    {"data_source": "GSE216606_CellTag_multi_RNA"},
    index=neutral_ids,
)
adata.uns["validation_design"] = (
    "Protocol-frozen holdout. Day, clone identity, fate, author cell state, and "
    "technical batch are absent from the inference object and retained only in "
    "separate audit tables. The runner must derive clusters without annotations."
)
adata.uns["barcode_mapping"] = (
    "GSM6681128-1131 map to d5r1 suffixes 1-4; GSM6681132-1135 map "
    "to d5r2 suffixes 1-4, verified by exact barcode intersections."
)

forbidden_obs = pd.Index(adata.obs.columns).str.contains(
    r"day|time|fate|clone|pseudotime|cell.?type|cluster|batch",
    case=False,
    regex=True,
)
assert not forbidden_obs.any()
assert adata.obs_names.is_unique
assert annotation.index.is_unique
assert adata.n_obs == mapping_qc["n_author_truth_matched"].sum()

truth_holdout = pd.DataFrame(
    {
        "cell_id": neutral_ids,
        "original_barcode": original_barcodes,
        "clone_id": annotation["clone.id"].to_numpy(),
        "day": annotation["day"].to_numpy(),
        "day_numeric": annotation["day"].map({"d2": 2.5, "d5": 5.0}).to_numpy(),
        "author_cell_type": annotation["cell_type2"].to_numpy(),
        "author_fate": annotation["fate"].to_numpy(),
        "author_fate_pct": annotation["fate_pct"].to_numpy(),
        "clone_size_rna": annotation["clone.size (RNA)"].to_numpy(),
        "n_d2_cells_rna": annotation["# of D2 cells (RNA)"].to_numpy(),
        "n_d5_cells_rna": annotation["# of D5 cells (RNA)"].to_numpy(),
    }
)
truth_holdout.to_csv(truth_path, index=False)

batch_mapping = pd.DataFrame(
    {
        "geo_accession": [
            "GSM6681127",
            "GSM6681128",
            "GSM6681129",
            "GSM6681130",
            "GSM6681131",
            "GSM6681132",
            "GSM6681133",
            "GSM6681134",
            "GSM6681135",
        ],
        "technical_batch": [
            "batch_01",
            "batch_02",
            "batch_03",
            "batch_04",
            "batch_05",
            "batch_06",
            "batch_07",
            "batch_08",
            "batch_09",
        ],
        "author_prefix": [
            "d2_5-rna",
            "d5r1-rna",
            "d5r1-rna",
            "d5r1-rna",
            "d5r1-rna",
            "d5r2-rna",
            "d5r2-rna",
            "d5r2-rna",
            "d5r2-rna",
        ],
        "author_suffix": [1, 1, 2, 3, 4, 1, 2, 3, 4],
    }
)
batch_mapping.to_csv(batch_map_path, index=False)

adata.write_h5ad(canonical_path, compression="gzip")

print(adata)
print(mapping_qc.to_string(index=False))
print("Inference obs columns:", adata.obs.columns.tolist())
