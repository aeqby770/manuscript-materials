#!/usr/bin/env python
"""scTour adapter for canonical h5ad input.

Usage: python run_sctour.py <canonical.h5ad> <output_pseudotime.csv>
"""
import os
import sys

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
import _mem_guard
_mem_guard.cap_self(55)

VENDOR_ROOT = os.environ.get(
    "LLR_VENDOR_ROOT",
    r"F:\biotrade\LineaLinkR_benchmark_assets\02_vendor_algorithms",
)
CONFIG = {
    "paths": {
        "vendor": os.path.join(VENDOR_ROOT, "scTour"),
    },
    "threads": os.cpu_count(),
    "seed": 12345,
    "plot": {
        "dpi": 300,
        "font": "Arial",
        "fontsize": 12,
        "palette": "viridis",
    },
    # scTour vendor defaults (model.py / train.py): full model capacity, no down-scaling.
    "nepoch": None,       # None -> vendor adaptive: min(round((10000/n_cells)*400), 400)
    "batch_size": 1024,
    "percent": None,      # None -> vendor adaptive: 0.2 if n_cells>10000 else 0.9
    "n_latent": 5,
    "n_ode_hidden": 25,
    "n_vae_hidden": 128,
}

os.environ["OMP_NUM_THREADS"] = str(CONFIG["threads"])
np.random.seed(CONFIG["seed"])
sys.path.insert(0, CONFIG["paths"]["vendor"])

import torch  # noqa: E402
import sctour as sct  # noqa: E402


def main():
    if len(sys.argv) != 3:
        print("Usage: run_sctour.py <canonical.h5ad> <output_pseudotime.csv>", file=sys.stderr)
        sys.exit(2)

    inpath, outpath = sys.argv[1], sys.argv[2]
    torch.manual_seed(CONFIG["seed"])
    adata = sc.read_h5ad(inpath)
    canonical_cell_ids = adata.obs_names.astype(str).to_numpy(copy=True)
    if canonical_cell_ids.size != adata.n_obs or pd.Index(canonical_cell_ids).has_duplicates:
        raise RuntimeError("scTour: canonical obs_names must be complete and unique")

    loss_mode = "nb"
    if "counts" in adata.layers:
        x = adata.layers["counts"]
        x = x.toarray() if sparse.issparse(x) else np.asarray(x)
        adata.X = np.rint(np.nan_to_num(x, nan=0.0, posinf=0.0, neginf=0.0)).astype(np.float32)
    elif "lognorm" in adata.layers:
        loss_mode = "mse"
        x = adata.layers["lognorm"]
        x = x.toarray() if sparse.issparse(x) else np.asarray(x)
        adata.X = np.nan_to_num(x, nan=0.0, posinf=0.0, neginf=0.0).astype(np.float32)
    else:
        x = adata.X.toarray() if sparse.issparse(adata.X) else np.asarray(adata.X)
        x = np.nan_to_num(x, nan=0.0, posinf=0.0, neginf=0.0).astype(np.float32)
        xd = np.asarray(x).ravel()
        if xd.size and np.all(xd >= 0) and np.allclose(xd, np.round(xd)):
            loss_mode = "nb"
            adata.X = np.rint(x).astype(np.float32)
        else:
            loss_mode = "mse"
            adata.X = x.astype(np.float32)

    sc.pp.calculate_qc_metrics(adata, inplace=True, percent_top=None)
    trainer = sct.train.Trainer(
        adata,
        percent=CONFIG["percent"],
        n_latent=CONFIG["n_latent"],
        n_ode_hidden=CONFIG["n_ode_hidden"],
        n_vae_hidden=CONFIG["n_vae_hidden"],
        nepoch=CONFIG["nepoch"],
        batch_size=min(CONFIG["batch_size"], adata.n_obs),
        loss_mode=loss_mode,
        random_state=CONFIG["seed"],
        # Use CUDA when available; otherwise use the CPU backend.
        use_gpu=torch.cuda.is_available(),
    )
    trainer.train()
    pseudotime = np.asarray(trainer.get_time(), dtype=float).reshape(-1)
    if not np.array_equal(adata.obs_names.astype(str).to_numpy(), canonical_cell_ids):
        raise RuntimeError("scTour: training changed canonical cell order")
    if pseudotime.size != canonical_cell_ids.size:
        raise RuntimeError(
            f"scTour: pseudotime length {pseudotime.size} != canonical cells {canonical_cell_ids.size}"
        )
    if not np.isfinite(pseudotime).all() or np.ptp(pseudotime) <= 0:
        raise RuntimeError("scTour: pseudotime must be finite for every cell and have non-zero range")

    pd.DataFrame({"cell_id": canonical_cell_ids, "pseudotime": pseudotime}).to_csv(outpath, index=False)
    written = pd.read_csv(outpath, dtype={"cell_id": str})
    written_pt = pd.to_numeric(written["pseudotime"], errors="coerce").to_numpy(dtype=float)
    if (
        written.columns.tolist() != ["cell_id", "pseudotime"]
        or not np.array_equal(written["cell_id"].astype(str).to_numpy(), canonical_cell_ids)
        or written["cell_id"].astype(str).duplicated().any()
        or not np.isfinite(written_pt).all()
        or np.ptp(written_pt) <= 0
    ):
        raise RuntimeError("scTour: written CSV failed the full-cell pseudotime contract")

    print(f"scTour OK: {adata.n_obs} cells -> {outpath}")


if __name__ == "__main__":
    main()
