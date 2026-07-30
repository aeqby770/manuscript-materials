#!/usr/bin/env python
"""CellPhoneDB statistical_analysis from canonical h5ad."""
import os
import sys
import tempfile
from multiprocessing import freeze_support
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

from cellphonedb.src.core.methods import cpdb_statistical_analysis_method
from cellphonedb.utils import db_utils
import _mem_guard

SCRIPT_DIR = Path(__file__).resolve().parent
WORK_ROOT = Path(
    os.environ.get("LLR064_WORK_ROOT", SCRIPT_DIR.parents[1])
).resolve()

CONFIG = {
    "paths": {
        "cache": WORK_ROOT / "cache" / "cellphonedb",
        "linealinkr": WORK_ROOT / "cache" / "linealinkr_workspace",
        "r_export": WORK_ROOT / "cache" / "r_export",
    },
    "threads": min(8, os.cpu_count() or 1),
    "cellphonedb_threads": min(8, os.cpu_count() or 1),
    "seed": 12345,
    "plot": {
        "figure.dpi": 300,
        "savefig.dpi": 300,
        "font.family": "Arial",
        "font.size": 12,
    },
    "release": "v4.1.0",
    # Panel default: 100 permutations (official allows 100+; 1000 x full gene set = days).
    "iterations": 100,
    "max_genes": 3000,
    "mem_budget_gb": 55,
}

os.environ["OMP_NUM_THREADS"] = str(CONFIG["threads"])
np.random.seed(CONFIG["seed"])
plt.rcParams.update(CONFIG["plot"])
_mem_guard.cap_self(CONFIG["mem_budget_gb"])


def pick_group_col(adata):
    for cn in (
        "cell_type", "clusters", "cluster", "celltype",
        "seurat_clusters", "auto_clusters", "state_info", "cell_ontology_class",
        "ct", "ct_subtype", "celltype_l3", "lineage",
    ):
        if cn not in adata.obs.columns:
            continue
        vals = adata.obs[cn].astype(str).fillna("")
        if len({v for v in vals if v and v != "nan"}) >= 2:
            return cn
    raise ValueError("CellPhoneDB: no cluster column with >= 2 groups")


def merge_extract_clusters(adata, ds_name):
    """GEO raw h5ad often lacks cluster; use LLR extract or r_export meta."""
    roots = [
        CONFIG["paths"]["linealinkr"]
        / ds_name
        / "extract"
        / f"{ds_name}_meta.csv",
        CONFIG["paths"]["r_export"] / ds_name / "meta.csv",
    ]
    for meta_path in roots:
        if not os.path.exists(meta_path):
            continue
        ext = pd.read_csv(meta_path)
        id_col = None
        for cn in ("barcode_id", "cell", "Cell", "cell_id"):
            if cn in ext.columns:
                id_col = cn
                break
        if id_col is None and ext.columns[0] not in (
            "clusters", "cluster", "auto_clusters", "cell_type",
        ):
            id_col = ext.columns[0]
        if id_col is None:
            continue
        ext[id_col] = ext[id_col].astype(str)
        ext = ext.set_index(id_col)
        for cn in ("clusters", "cluster", "auto_clusters", "seurat_clusters", "cell_type"):
            if cn not in ext.columns:
                continue
            vals = ext[cn].astype(str)
            mapped = vals.reindex(adata.obs_names.astype(str))
            if mapped.notna().sum() < 2:
                continue
            nu = mapped.dropna().nunique()
            if nu < 2:
                continue
            adata.obs[cn] = mapped.values
            return adata
    return adata


def counts_layer(adata):
    if "counts" in adata.layers:
        return adata.layers["counts"]
    if adata.raw is not None:
        return adata.raw.X
    return adata.X


def select_genes(adata, max_genes):
    x = counts_layer(adata)
    if sparse.issparse(x):
        nnz = np.asarray((x > 0).sum(axis=0)).ravel()
    else:
        nnz = (np.asarray(x) > 0).sum(axis=0)
    keep = np.where(nnz >= max(10, int(0.01 * adata.n_obs)))[0]
    if keep.size == 0:
        raise ValueError("CellPhoneDB: no expressed genes")
    if keep.size > max_genes:
        if sparse.issparse(x):
            means = np.asarray(x[:, keep].mean(axis=0)).ravel()
        else:
            means = np.asarray(x)[:, keep].mean(axis=0)
        order = np.argsort(means)[::-1][:max_genes]
        keep = keep[order]
    return adata.var_names[np.asarray(keep, dtype=int)].astype(str)


def dense_counts(adata, genes):
    x = counts_layer(adata)
    idx = [adata.var_names.get_loc(g) for g in genes]
    if sparse.issparse(x):
        sub = x[:, idx]
        return sub.toarray()
    return np.asarray(x)[:, idx]


if __name__ == '__main__':
    freeze_support()

    inpath, outdir = sys.argv[1], sys.argv[2]
    os.makedirs(outdir, exist_ok=True)
    marker = os.path.join(outdir, "cpdb_done.txt")

    adata = sc.read_h5ad(inpath)
    ds_name = os.path.basename(os.path.normpath(outdir))
    adata = merge_extract_clusters(adata, ds_name)

    group_col = pick_group_col(adata)
    genes = select_genes(adata, CONFIG["max_genes"])
    est_gb = adata.n_obs * len(genes) * 8 * 3 / (1024**3)
    _mem_guard.require_under(est_gb, CONFIG["mem_budget_gb"], "CellPhoneDB/densify")

    os.makedirs(CONFIG["paths"]["cache"], exist_ok=True)
    db_utils.download_database(str(CONFIG["paths"]["cache"]), CONFIG["release"])
    cpdb_file = str(CONFIG["paths"]["cache"] / "cellphonedb.zip")

    with tempfile.TemporaryDirectory() as tmp:
        mat = dense_counts(adata, genes)
        counts = pd.DataFrame(
            mat.T,
            index=genes,
            columns=adata.obs_names.astype(str),
        )
        meta = pd.DataFrame(
            {"cell_type": adata.obs[group_col].astype(str).values},
            index=adata.obs_names.astype(str),
        )
        meta.index.name = "Cell"
        counts_path = os.path.join(tmp, "counts.txt")
        meta_path = os.path.join(tmp, "meta.txt")
        counts.to_csv(counts_path, sep="\t")
        meta.to_csv(meta_path, sep="\t")

        result = cpdb_statistical_analysis_method.call(
            cpdb_file_path=cpdb_file,
            meta_file_path=meta_path,
            counts_file_path=counts_path,
            counts_data="gene_name",
            output_path=outdir,
            iterations=CONFIG["iterations"],
            threshold=0.1,
            threads=CONFIG["cellphonedb_threads"],
            debug_seed=CONFIG["seed"],
            result_precision=3,
            pvalue=0.05,
            subsampling=False,
            separator="|",
            debug=False,
            output_suffix=None,
            score_interactions=False,
        )

    assert not result["means"].empty
    assert not result["pvalues"].empty

    with open(marker, "w", encoding="utf-8") as fh:
        fh.write("OK\n")
