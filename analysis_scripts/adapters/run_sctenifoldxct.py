#!/usr/bin/env python
"""scTenifoldXct L-R candidate scoring across ALL ordered cell-type pairs (full matrix)."""
import os
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

CONFIG = {
    "paths": {
        "vendor": Path(os.environ["LLR_VENDOR_ROOT"]).resolve() / "scTenifoldXct",
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
from scTenifoldXct import scTenifoldXct  # noqa: E402

import _mem_guard
_mem_guard.cap_self(55)


def pick_group(adata):
    for c in ("cell_type", "cluster", "clusters", "celltype"):
        if c in adata.obs.columns and adata.obs[c].nunique() >= 2:
            return c
    sc.pp.neighbors(adata, n_neighbors=min(15, adata.n_obs - 1))
    sc.tl.leiden(adata, resolution=0.5, key_added="sctenifold_leiden")
    return "sctenifold_leiden"


def main():
    if len(sys.argv) < 3:
        print("Usage: run_sctenifoldxct.py <h5ad> <output.csv>", file=sys.stderr)
        sys.exit(2)
    inpath, outpath = sys.argv[1], sys.argv[2]
    if os.path.exists(outpath) and os.path.getsize(outpath) > 100:
        print("CACHED")
        return

    adata = sc.read_h5ad(inpath)
    # scTenifoldXct builds gene-gene GRNs + manifold alignment (~O(n_genes^2)).
    _mem_guard.require_under(_mem_guard.square_gb(adata.n_vars, copies=4), 55, "scTenifoldXct")
    if "lognorm" in adata.layers:
        x = adata.layers["lognorm"]
        adata.X = x.toarray() if sparse.issparse(x) else np.asarray(x)
    else:
        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)

    gcol = pick_group(adata)
    groups = sorted(adata.obs[gcol].astype(str).unique().tolist())
    if len(groups) < 2:
        # Fewer than two groups does not satisfy the method's input contract.
        print("FAIL: fewer than 2 cell groups; scTenifoldXct needs >=2 to compare", file=sys.stderr)
        sys.exit(1)

    cache = os.path.join(os.path.dirname(outpath), "grn_cache")
    pair_dir = os.path.join(os.path.dirname(outpath), "sctenifoldxct_pairs")
    os.makedirs(pair_dir, exist_ok=True)
    all_results = []
    # Full cell-type x cell-type matrix (no top-2 restriction). Each ordered pair
    # runs the official scTenifoldXct end to end: build GRNs, train the manifold-
    # aligned NN (get_embeds), then null_test() for the significance test — no step
    # skipped.
    for src in groups:
        for tgt in groups:
            if src == tgt:
                continue
            tag = re.sub(r"[^0-9A-Za-z_.-]", "_", f"{src}__{tgt}")
            pair_csv = os.path.join(pair_dir, f"{tag}.csv")
            # Per-pair disk cache (resumable, identical result).
            if os.path.exists(pair_csv) and os.path.getsize(pair_csv) > 0:
                all_results.append(pd.read_csv(pair_csv))
                continue
            xct = scTenifoldXct(
                data=adata,
                source_celltype=str(src),
                target_celltype=str(tgt),
                obs_label=gcol,
                GRN_file_dir=cache,
                rebuild_GRN=True,
                query_DB="pairs",
                verbose=False,
            )
            xct.get_embeds(train=True)
            df = xct.null_test()
            df = df.copy()
            df["source_celltype"] = src
            df["target_celltype"] = tgt
            df.to_csv(pair_csv, index=False)
            all_results.append(df)

    if not all_results:
        print("FAIL: no cell-type pair produced results")
        sys.exit(1)
    out = pd.concat(all_results, ignore_index=True)
    out["dataset"] = os.path.basename(inpath)
    out.to_csv(outpath, index=False)
    print(f"OK: {len(out)} significant LR rows across {len(groups)} cell types -> {outpath}")


if __name__ == "__main__":
    main()
