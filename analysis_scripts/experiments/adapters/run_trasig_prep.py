#!/usr/bin/env python
"""Prepare TraSig inputs from registered Slingshot pseudotime and run TraSig."""
import gc
import os
import subprocess
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pickle
import scanpy as sc
from scipy import sparse
import _mem_guard

SCRIPT_DIR = Path(__file__).resolve().parent
SOURCE_ASSETS_ROOT = Path(os.environ.get("LLR_WORK_ROOT", SCRIPT_DIR.parents[1])).resolve()

CONFIG = {
    "paths": {
        "source_assets": SOURCE_ASSETS_ROOT,
    },
    "threads": os.cpu_count() or 1,
    "seed": 12345,
    "plot": {},
    "mem_budget_gb": 55,
    "permutation_workers": 4,
}

os.environ["OMP_NUM_THREADS"] = str(CONFIG["threads"])
np.random.seed(CONFIG["seed"])
plt.rcParams.update(CONFIG["plot"])
_mem_guard.cap_self(CONFIG["mem_budget_gb"])


def main():
    if len(sys.argv) < 3:
        print(
            "Usage: run_trasig_prep.py <h5ad> <out_dir> [dataset] [source_assets] [lr_resource] [group_column|lineage:column] [input_scale]",
            file=sys.stderr,
        )
        sys.exit(2)
    inpath, outdir = sys.argv[1], sys.argv[2]
    ds = sys.argv[3] if len(sys.argv) > 3 else os.path.basename(inpath).replace("_canonical.h5ad", "").replace(".h5ad", "")
    source_assets = (
        Path(sys.argv[4]).resolve()
        if len(sys.argv) > 4
        else CONFIG["paths"]["source_assets"]
    )
    out_dir = Path(outdir).resolve()
    os.makedirs(out_dir, exist_ok=True)
    lr_resource = Path(sys.argv[5]) if len(sys.argv) > 5 else (
        source_assets / "tools" / "competitors" / "original" / "TraSig" /
        "ligand_receptor_lists" / "ligand_receptor_FANTOM.pickle"
    )
    group_column = (sys.argv[6] or None) if len(sys.argv) > 6 else None
    input_scale = sys.argv[7] if len(sys.argv) > 7 else "counts"

    slingshot_dir = out_dir.parents[1] / "slingshot" / ds
    pt_file = os.path.join(slingshot_dir, "slingshot_lineage_pseudotime.csv")
    weight_file = os.path.join(slingshot_dir, "slingshot_weights.csv")

    adata = sc.read_h5ad(inpath)
    exclusions = pd.read_csv(
        SCRIPT_DIR.parent / "config" / "feature_exclusions.csv", encoding="utf-8"
    )
    excluded_genes = exclusions.loc[exclusions["dataset"].eq(ds), "gene"]
    adata = adata[:, ~adata.var_names.isin(excluded_genes)].copy()
    use_lognorm = input_scale == "counts" and "lognorm" in adata.layers and excluded_genes.empty
    if not use_lognorm and input_scale == "counts":
        adata.X = adata.layers["counts"].copy() if "counts" in adata.layers else adata.X
        sc.pp.normalize_total(adata, target_sum=10000)
        sc.pp.log1p(adata)
    print(f"TraSig expression input scale: {input_scale}", flush=True)
    # TraSig consumes ligand/receptor-pair genes extracted below. Apply no HVG cap so
    # the complete detected pair universe remains available to the method.

    trajectory_ids = pd.read_csv(
        slingshot_dir / "slingshot_pseudotime.csv", usecols=["cell_id"], dtype=str
    )["cell_id"]
    if trajectory_ids.duplicated().any():
        raise ValueError("Slingshot cell IDs must be unique.")
    row_order = pd.Index(trajectory_ids).get_indexer(adata.obs_names.astype(str))
    if (row_order < 0).any() or len(trajectory_ids) != adata.n_obs:
        raise ValueError("Slingshot cell IDs do not match the input H5AD.")
    pt_mat = pd.read_csv(pt_file).to_numpy(dtype=float)[row_order]
    weight_mat = pd.read_csv(weight_file).to_numpy(dtype=float)[row_order]
    if pt_mat.shape != weight_mat.shape or pt_mat.shape[0] != adata.n_obs:
        raise ValueError("Slingshot lineage pseudotime and weights must match the H5AD rows.")
    cell_path = np.nanargmax(np.nan_to_num(weight_mat, nan=-1), axis=1)
    cell_time = pt_mat[np.arange(pt_mat.shape[0]), cell_path]
    if group_column is not None:
        lineage_celltype = group_column.startswith("lineage:")
        group_column = group_column.removeprefix("lineage:")
        group_labels = pd.MultiIndex.from_arrays([
            adata.obs[group_column].astype(str),
            np.where(lineage_celltype, cell_path, -1),
        ], names=["cell_type", "trajectory"])
        cell_path, path_labels = pd.factorize(group_labels, sort=True)
        path_labels = path_labels.set_names(group_labels.names)
        global_time = pd.read_csv(
            slingshot_dir / "slingshot_pseudotime.csv"
        )["pseudotime"].to_numpy(dtype=float)[row_order]
        cell_time = np.where(lineage_celltype, cell_time, global_time)
        pd.DataFrame({"path": np.arange(len(path_labels)),
            "cell_type": path_labels.get_level_values("cell_type"),
            "trajectory": path_labels.get_level_values("trajectory")}).to_csv(
            out_dir / "trasig_path_labels.csv", index=False
        )
    cell_time = np.asarray(cell_time, dtype=float)
    if not np.all(np.isfinite(cell_time)):
        # Non-finite slingshot values identify cells without an assigned position.
        print("FAIL: slingshot pseudotime contains non-finite entries.",
              file=sys.stderr)
        sys.exit(1)
    cell_time = (cell_time - cell_time.min()) / (cell_time.max() - cell_time.min())
    pd.DataFrame({
        "cell_id": adata.obs_names.astype(str),
        "cell_path": cell_path, "cell_time": cell_time,
    }).to_csv(out_dir / "trasig_cell_inputs.csv", index=False)

    hid_var = {"cell_path": cell_path.astype(int), "cell_time": cell_time.astype(float)}
    with open(os.path.join(outdir, "hid_var.pickle"), "wb") as fh:
        pickle.dump(hid_var, fh)
    with open(os.path.join(outdir, f"{ds}_log1p_slingshot_it2_hid_var.pickle"), "wb") as fh:
        pickle.dump(hid_var, fh)

    trasig_root = source_assets / "tools" / "competitors" / "original" / "TraSig"
    if lr_resource.suffix == ".pickle":
        with open(lr_resource, "rb") as fh:
            lr_pairs = pickle.load(fh)
    else:
        lr_table = pd.read_csv(lr_resource, sep="\t")
        lr_pairs = list(map(tuple, lr_table[
            ["ligand_gene_symbol", "receptor_gene_symbol"]
        ].drop_duplicates().to_numpy()))
    genes = pd.Index(adata.var_names.astype(str))
    gene_set = set(genes)
    lr_pairs = [(ligand, receptor) for ligand, receptor in lr_pairs if ligand in gene_set and receptor in gene_set]
    if len(lr_pairs) == 0:
        print("FAIL: no ligand-receptor pairs overlap dataset genes", file=sys.stderr)
        sys.exit(1)
    with open(os.path.join(outdir, f"lr_{ds}_log1p.pickle"), "wb") as fh:
        pickle.dump(lr_pairs, fh)
    expr_genes = sorted(set([gene for pair in lr_pairs for gene in pair]))
    adata_lr = adata[:, expr_genes]
    expr = adata_lr.layers["lognorm"] if use_lognorm else adata_lr.X
    if sparse.issparse(expr):
        expr = expr.toarray()
    pd.DataFrame(expr, index=adata.obs_names, columns=expr_genes).to_csv(os.path.join(outdir, f"{ds}_log1p_lr.txt"))
    with open(os.path.join(outdir, f"sampling_time_per_path_{ds}_log1p_slingshot.pickle"), "wb") as fh:
        pickle.dump({0: list(np.unique(cell_path).astype(int))}, fh)
    del expr, adata_lr, adata
    gc.collect()

    if "--prepare-only" in sys.argv[8:]:
        return

    trasig_main = trasig_root / "main.py"
    if not os.path.isfile(trasig_main):
        print(f"FAIL: TraSig main.py not found at {trasig_main}", file=sys.stderr)
        sys.exit(1)
    res_dir = os.path.join(outdir, "results")
    os.makedirs(res_dir, exist_ok=True)
    # Propagate the TraSig exit status to the benchmark runner.
    subprocess.run(
        [
            sys.executable,
            str(trasig_main),
            "-i",
            outdir,
            "-o",
            res_dir,
            "-d",
            ds,
            "-g",
            "log1p",
            "-b",
            "slingshot",
            "-t",
            "lr",
            "-l",
            "20",
            "-m",
            "dot",
            "-n",
            "10000",
            "-p",
            "true",
            "-c",
            str(CONFIG["permutation_workers"]),
        ],
        check=True,
    )
    permutation_file = os.path.join(
        res_dir,
        f"{ds}_lr_log1p_slingshot_smallerWindow_nlap_20_4_permutation_results.pickle",
    )
    with open(permutation_file, "rb") as fh:
        permutation_results = pickle.load(fh)
    permutation_counts = np.concatenate(
        [
            np.asarray(values["dot"], dtype=np.int64)
            for values in permutation_results.values()
        ]
    )
    if permutation_counts.size == 0 or permutation_counts.min() < 0 or permutation_counts.max() > 10000:
        raise ValueError("TraSig permutation counts must be within 0 and 10000.")
    observed_file = os.path.join(
        res_dir,
        f"{ds}_lr_log1p_slingshot_smallerWindow_nlap_20_metrics_0.pickle",
    )
    with open(observed_file, "rb") as fh:
        observed_results = pickle.load(fh)
    score_tables = [
        pd.DataFrame(
            {
                "path_pair": path_pair,
                "ligand": np.asarray(lr_pairs, dtype=object)[:, 0],
                "receptor": np.asarray(lr_pairs, dtype=object)[:, 1],
                "observed_score": np.asarray(observed_results[path_pair]["dot"], dtype=float),
                "permutation_count": np.asarray(permutation_results[path_pair]["dot"], dtype=np.int64),
                "permutation_p_value": (
                    np.asarray(permutation_results[path_pair]["dot"], dtype=np.int64) + 1
                ) / 10001,
            }
        )
        for path_pair in observed_results
    ]
    if any(len(table) != len(lr_pairs) for table in score_tables):
        raise ValueError("TraSig score dimensions do not match the ligand-receptor list.")
    pd.concat(score_tables, ignore_index=True).to_csv(
        out_dir / "trasig_scores.csv", index=False
    )
    print(f"OK: TraSig results in {res_dir}")


if __name__ == "__main__":
    main()
