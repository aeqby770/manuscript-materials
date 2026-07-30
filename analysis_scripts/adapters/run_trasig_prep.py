#!/usr/bin/env python
"""TraSig input prep via slingshot adapter (no rpy2). Optional TraSig main.py."""
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
PROJECT_ROOT = Path(
    os.environ.get("LLR064_WORK_ROOT", SCRIPT_DIR.parents[1])
).resolve()
SOURCE_ASSETS_ROOT = Path(
    os.environ.get(
        "LLR064_SOURCE_ASSETS",
        PROJECT_ROOT.parent / "LineaLinkR_benchmark_assets",
    )
).resolve()

CONFIG = {
    "paths": {
        "rscript": Path(
            os.environ.get("LLR_RSCRIPT", r"D:\R\R-4.5.3\bin\Rscript.exe")
        ).resolve(),
        "slingshot_adapter": SCRIPT_DIR / "run_slingshot.R",
        "source_assets": SOURCE_ASSETS_ROOT,
    },
    "threads": os.cpu_count() or 1,
    "seed": 12345,
    "plot": {},
    "mem_budget_gb": 55,
}

os.environ["OMP_NUM_THREADS"] = str(CONFIG["threads"])
np.random.seed(CONFIG["seed"])
plt.rcParams.update(CONFIG["plot"])
_mem_guard.cap_self(CONFIG["mem_budget_gb"])


def main():
    if len(sys.argv) < 3:
        print(
            "Usage: run_trasig_prep.py <h5ad> <out_dir> [dataset] [source_assets]",
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

    slingshot_dir = out_dir.parents[1] / "slingshot" / ds
    os.makedirs(slingshot_dir, exist_ok=True)
    subprocess.run(
        [
            str(CONFIG["paths"]["rscript"]),
            "--encoding=UTF-8",
            str(CONFIG["paths"]["slingshot_adapter"]),
            ds,
            str(slingshot_dir),
        ],
        check=False,
    )

    pt_file = os.path.join(slingshot_dir, "slingshot_pseudotime.csv")
    if not os.path.exists(pt_file):
        print("FAIL: slingshot output missing", pt_file)
        sys.exit(1)

    adata = sc.read_h5ad(inpath)
    if "lognorm" in adata.layers:
        x = adata.layers["lognorm"]
        adata.X = x.toarray() if sparse.issparse(x) else np.asarray(x)
    else:
        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)
    # TraSig consumes ligand/receptor-pair genes extracted below. Apply no HVG cap so
    # the complete detected pair universe remains available to the method.

    pt_df = pd.read_csv(pt_file)
    if pt_df.shape[1] < 1 or len(pt_df) == 0:
        print("FAIL: empty slingshot pseudotime file", file=sys.stderr)
        sys.exit(1)
    first_col = str(pt_df.columns[0])
    if first_col.startswith("Lineage") or first_col.lower() in ("pseudotime",):
        if len(pt_df) != adata.n_obs:
            print(f"FAIL: slingshot rows {len(pt_df)} != n_obs {adata.n_obs}", file=sys.stderr)
            sys.exit(1)
        pt_df.index = list(adata.obs_names)
    else:
        pt_df = pd.read_csv(pt_file, index_col=0)
    if "Lineage1" in pt_df.columns:
        pt_mat = pt_df.values
        cell_path = np.nanargmax(np.nan_to_num(pt_mat, nan=-1), axis=1)
        cell_time = pt_mat[np.arange(pt_mat.shape[0]), cell_path]
    else:
        col = [c for c in pt_df.columns if "lineage" in c.lower() or c == "pseudotime"]
        cell_time = pt_df[col[0]].values if col else pt_df.iloc[:, 0].values
        cell_path = np.zeros(len(cell_time), dtype=int)
    cell_time = np.asarray(cell_time, dtype=float)
    if not np.all(np.isfinite(cell_time)):
        # Non-finite slingshot values identify cells without an assigned position.
        print("FAIL: slingshot pseudotime contains non-finite entries.",
              file=sys.stderr)
        sys.exit(1)
    cell_time = (cell_time - cell_time.min()) / (cell_time.max() - cell_time.min())

    gcol = None
    for c in ("cell_type", "cluster", "clusters"):
        if c in adata.obs.columns:
            gcol = c
            break
    if gcol is None:
    # TraSig needs at least one observed cell-type or cluster label.
        print("FAIL: TraSig needs cell-type/cluster labels; none found in obs.", file=sys.stderr)
        sys.exit(1)
    labels = pd.Categorical(adata.obs[gcol].astype(str)).codes

    hid_var = {"cell_path": cell_path.astype(int), "cell_time": cell_time.astype(float), "cell_labels": labels}
    with open(os.path.join(outdir, "hid_var.pickle"), "wb") as fh:
        pickle.dump(hid_var, fh)
    with open(os.path.join(outdir, f"{ds}_log1p_slingshot_it2_hid_var.pickle"), "wb") as fh:
        pickle.dump(hid_var, fh)

    trasig_root = source_assets / "02_vendor_algorithms" / "TraSig"
    lr_path = trasig_root / "ligand_receptor_lists" / "ligand_receptor_FANTOM.pickle"
    with open(lr_path, "rb") as fh:
        lr_pairs = pickle.load(fh)
    genes = pd.Index(adata.var_names.astype(str))
    gene_set = set(genes)
    lr_pairs = [(ligand, receptor) for ligand, receptor in lr_pairs if ligand in gene_set and receptor in gene_set]
    if len(lr_pairs) == 0:
        print("FAIL: no ligand-receptor pairs overlap dataset genes", file=sys.stderr)
        sys.exit(1)
    with open(os.path.join(outdir, f"lr_{ds}_log1p.pickle"), "wb") as fh:
        pickle.dump(lr_pairs, fh)
    expr_genes = sorted(set([gene for pair in lr_pairs for gene in pair]))
    expr = adata[:, expr_genes].X
    if sparse.issparse(expr):
        expr = expr.toarray()
    pd.DataFrame(expr, index=adata.obs_names, columns=expr_genes).to_csv(os.path.join(outdir, f"{ds}_log1p_lr.txt"))
    with open(os.path.join(outdir, f"sampling_time_per_path_{ds}_log1p_slingshot.pickle"), "wb") as fh:
        pickle.dump({0: list(np.unique(cell_path).astype(int))}, fh)

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
            "false",
        ],
        check=True,
    )
    print(f"OK: TraSig results in {res_dir}")


if __name__ == "__main__":
    main()
