#!/usr/bin/env python
"""LIANA rank_aggregate on full canonical h5ad (all cells, no subsampling)."""
import sys
import os
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import pandas as pd
from scipy import sparse
import _mem_guard

CONFIG = {
    "paths": {},
    "threads": os.cpu_count() or 1,
    "seed": 12345,
    "plot": {},
    "mem_budget_gb": 55,
    "content_contract": "liana_rank_aggregate_native_v1",
}

os.environ["OMP_NUM_THREADS"] = str(CONFIG["threads"])
np.random.seed(CONFIG["seed"])
plt.rcParams.update(CONFIG["plot"])
_mem_guard.cap_self(CONFIG["mem_budget_gb"])

def main():
    if len(sys.argv) < 3:
        print("Usage: run_liana.py <canonical.h5ad> <output.csv>", file=sys.stderr)
        sys.exit(2)
    inpath, outpath = sys.argv[1], sys.argv[2]
    resource = sys.argv[3] if len(sys.argv) > 3 else "consensus"
    input_path = Path(inpath).resolve()
    output_path = Path(outpath).resolve()
    output_prefix = output_path.name.removesuffix("_result.csv")
    dataset = output_prefix.removeprefix("liana_")
    cells_path = output_path.with_name(f"{output_prefix}_cells.csv")
    run_path = output_path.with_name(f"{output_prefix}_run.json")

    adata = sc.read_h5ad(input_path)
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

    def pick_group(adata):
        best = None
        for c in ("cell_type", "cluster", "clusters", "celltype", "louvain"):
            if c not in adata.obs.columns:
                continue
            n = adata.obs[c].nunique()
            if n >= 2 and (best is None or n < best[1]):
                best = (c, n)
        if best is not None and best[1] <= 80:
            return best[0]
        sc.pp.neighbors(adata, n_neighbors=min(15, adata.n_obs - 1))
        sc.tl.leiden(adata, resolution=0.5, key_added="liana_leiden")
        return "liana_leiden"

    g = pick_group(adata)
    adata.obs[g] = adata.obs[g].astype(str)

    import liana as ln

    # Official: liana.method.rank_aggregate (see LIANA tests/test_rank_aggregate.py)
    # Accept both human-style HGNC (`TNF`, `CXCL8`) and mouse-style (`Actb`, `Cxcl12`)
    # symbols; reject simulated/non-symbol IDs before LIANA runs expensive steps.
    genes = pd.Index(adata.var_names).astype(str)
    is_mouse_resource = resource == "mouseconsensus"
    human_like = genes.str.match(r"^[A-Z][A-Z0-9.-]+$")
    mouse_like = genes.str.match(r"^[A-Z][a-z0-9.-]+$")
    looks_symbolic = (
        mouse_like.mean() > 0.2 if is_mouse_resource else human_like.mean() > 0.2
    )
    if not looks_symbolic:
        print(
            "FAIL: gene symbols incompatible with LIANA resource (simulated/non-HGNC names). "
            "Use human/mouse symbol datasets; do not emit placeholder LR tables.",
            file=sys.stderr,
        )
        sys.exit(1)

    ln.mt.rank_aggregate(adata, groupby=g, resource_name=resource, use_raw=False, verbose=True)

    res = adata.uns.get("rank_aggregate")
    if res is None:
        res = adata.uns.get("liana_res")
    if res is None:
        for k in adata.uns.keys():
            if "rank" in k.lower() or "liana" in k.lower():
                res = adata.uns[k]
                break
    if res is None:
        print("FAIL: no rank_aggregate in adata.uns", list(adata.uns.keys())[:20])
        sys.exit(1)
    def unwrap_df(x):
        if isinstance(x, pd.DataFrame) and len(x):
            return x
        if isinstance(x, dict):
            for v in x.values():
                got = unwrap_df(v)
                if got is not None:
                    return got
        return None
    got = unwrap_df(res)
    if got is not None:
        res = got
    if not isinstance(res, pd.DataFrame):
        print("FAIL: unexpected result type", type(res))
        sys.exit(1)
    output_temp = Path(f"{output_path}.tmp")
    cells_temp = Path(f"{cells_path}.tmp")
    run_temp = Path(f"{run_path}.tmp")
    pd.DataFrame({"cell_id": adata.obs_names.astype(str)}).to_csv(
        cells_temp,
        index=False,
    )
    run_temp.write_text(
        json.dumps(
            {
                "content_contract": CONFIG["content_contract"],
                "dataset": dataset,
                "canonical_source_path": str(input_path),
                "n_cells": int(adata.n_obs),
                "n_genes": int(adata.n_vars),
                "groupby": g,
                "resource": resource,
            },
            ensure_ascii=True,
            indent=2,
        )
        + "\n",
        encoding="utf-8",
        newline="\n",
    )
    res.to_csv(output_temp, index=False)
    cells_temp.replace(cells_path)
    run_temp.replace(run_path)
    output_temp.replace(output_path)
    print(
        f"OK: {len(res)} rows -> {output_path}; "
        f"{adata.n_obs} cells -> {cells_path}; provenance -> {run_path}"
    )

if __name__ == "__main__":
    main()
