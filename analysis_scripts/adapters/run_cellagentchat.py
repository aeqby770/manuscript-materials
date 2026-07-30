#!/usr/bin/env python
"""CellAgentChat ABM CCI (official Communication.CCI, full permutation background)."""
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

CONFIG = {
    "paths": {
        "vendor_root": Path(os.environ["LLR_VENDOR_ROOT"]).resolve(),
    },
    "threads": os.cpu_count() or 1,
    "seed": 12345,
    "plot": {},
}
CONFIG["paths"]["src"] = CONFIG["paths"]["vendor_root"] / "CellAgentChat" / "src"

os.environ["OMP_NUM_THREADS"] = str(CONFIG["threads"])
np.random.seed(CONFIG["seed"])
SRC = CONFIG["paths"]["src"]
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))
import Communication  # noqa: E402
import permutations  # noqa: E402
from bckground_distribution import get_distribution  # noqa: E402

import _mem_guard
_mem_guard.cap_self(55)


def main():
    if len(sys.argv) < 3:
        print("Usage: run_cellagentchat.py <h5ad> <out_dir>", file=sys.stderr)
        sys.exit(2)
    inpath, outdir = sys.argv[1], sys.argv[2]
    os.makedirs(outdir, exist_ok=True)
    marker = os.path.join(outdir, "cellagentchat_done.txt")
    if os.path.exists(marker) and os.path.getsize(marker) > 0:
        print("CACHED")
        return

    adata = sc.read_h5ad(inpath)

    if "lognorm" in adata.layers:
        x = adata.layers["lognorm"]
        adata.X = x.toarray() if sparse.issparse(x) else np.asarray(x)
    else:
        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)

    gcol = None
    for c in ("cell_type", "cluster", "clusters"):
        if c in adata.obs.columns and adata.obs[c].nunique() >= 2:
            gcol = c
            break
    if gcol is None:
        sc.pp.neighbors(adata, n_neighbors=min(15, adata.n_obs - 1))
        sc.tl.leiden(adata, resolution=0.5, key_added="cac_leiden")
        gcol = "cac_leiden"
    adata.obs["cell_type"] = adata.obs[gcol].astype(str)

    # The vendor ships human_lr_pair.tsv. For mouse data, map symbols by title case
    # when the direct symbols do not overlap the dataset.
    sp = os.environ.get("LLR_DATASET_SPECIES", "").lower()
    lr_path = SRC / "human_lr_pair.tsv"
    lr = pd.read_csv(lr_path, sep="\t")
    genes = set(adata.var_names.astype(str))

    def _match(symbols):
        s = [str(x) for x in symbols]
        direct = [g for g in s if g in genes]
        if sp == "mouse" and len(direct) < 5:
            mapped = [g.title() for g in s if g.title() in genes]
            if len(mapped) > len(direct):
                return mapped
        return direct

    lig_uni = list(dict.fromkeys(_match(lr["ligand"].unique())))
    rec_uni = list(dict.fromkeys(_match(lr["receptor"].unique())))
    if len(lig_uni) < 5 or len(rec_uni) < 5:
        # Return a non-zero status when no supported ligand-receptor genes remain.
        print(f"FAIL: too few matched LR genes (lig={len(lig_uni)}, rec={len(rec_uni)}); "
              f"species={sp or 'unknown'} vs human LR DB", file=sys.stderr)
        sys.exit(1)

    clusters = ["Names"] + sorted(adata.obs["cell_type"].unique().astype(str).tolist())
    rates = np.ones(len(clusters) - 1)
    # threshold is the minimum average number of score values per LR pair before
    # the permutation background stops accumulating. The tutorial recommends
    # at least 10,000; use 10,000 for the full analysis.
    fin, _, _ = permutations.permutation_test(
        threshold=10000,
        N=adata.n_obs,
        adata=adata,
        lig_uni=lig_uni,
        rec_uni=rec_uni,
        rates=rates,
        dist=False,
        tau=2,
        max_steps=1,
        seed=12345,
    )
    dist = get_distribution(fin)

    Communication.CCI(
        N=adata.n_obs,
        adata=adata,
        lig_uni=lig_uni,
        rec_uni=rec_uni,
        rates=rates,
        distribution=dist,
        clusters=clusters,
        dist=1.0,
        delta=1.0,
        max_steps=1,
        tau=2,
        rec_block=False,
        plot_every_step=False,
        path=outdir,
    )
    with open(marker, "w", encoding="utf-8") as fh:
        fh.write("OK\n")
    print(f"OK: CellAgentChat -> {outdir}")


if __name__ == "__main__":
    main()
