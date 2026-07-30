#!/usr/bin/env python
"""Palantir — official scanpy workflow (Palantir_sample_notebook.ipynb).

Usage: run_palantir.py <h5ad> <output_h5ad>
The adapter uses all cells and requires usable counts or expression values.
"""
import os
import sys

import numpy as np
import pandas as pd
import scanpy as sc
import _mem_guard
from _benchmark_contract import expression_root_position
_mem_guard.cap_self(55)

CONFIG = {
    "paths": {},
    "threads": os.cpu_count(),
    "seed": 12345,
    "plot": {},
}

os.environ["OMP_NUM_THREADS"] = str(CONFIG["threads"])
np.random.seed(CONFIG["seed"])


ROOT_CANDIDATES = ("0", "root", "Root", "STEM", "stem", "branch_0", "early", "Early")


def benchmark_track():
    track = os.environ.get("LLR_BENCHMARK_TRACK", "annotation_free")
    if track not in {"annotation_free", "common_root", "legacy_annotated"}:
        raise ValueError(f"Unknown LLR_BENCHMARK_TRACK: {track}")
    return track


def pick_cluster_and_labels(ad):
    """Pick a cluster-like column and its labels from ad.obs.

    Preference order mirrors R-side helpers: clusters, cluster, cell_type, celltype, seurat_clusters.
    A column is accepted only if it has at least two non-empty, non-NA levels.
    """

    candidates = (
        ("clusters", "cluster", "cell_type", "celltype", "seurat_clusters")
        if benchmark_track() == "legacy_annotated"
        else ("benchmark_clusters",)
    )
    for cn in candidates:
        if cn not in ad.obs.columns:
            continue
        v = ad.obs[cn].astype(str)
        valid = v.notna() & (v != "")
        if v[valid].nunique() >= 2:
            return cn, v
    return None, None


def pick_common_root_cell(ad):
    root_cell = os.environ["LLR_BENCHMARK_ROOT_CELL"]
    if root_cell not in ad.obs_names:
        raise ValueError(f"Common-root cell is absent from input: {root_cell}")
    return root_cell


def attach_fallback_clusters(ad):
    """If no annotation clusters are available, derive clusters from PCA space.

    This mirrors the monocle3 fallback of using internal clustering when meta is missing.
    """

    if "X_pca" not in ad.obsm:
        sc.pp.pca(ad, n_comps=min(30, ad.n_vars - 1, ad.n_obs - 1))
    sc.pp.neighbors(ad, use_rep="X_pca", n_neighbors=30)
    sc.tl.leiden(ad, resolution=1.0, key_added="benchmark_clusters")
    v = ad.obs["benchmark_clusters"].astype(str)
    return "benchmark_clusters", v


def determine_root_type(labels):
    """Choose a root cluster consistent with PickMonocle3RootCells semantics.

    1) Prefer well-known labels (0/root/STEM/early...)
    2) Then minimum numeric cluster ID if labels are numeric
    3) Finally, the most frequent label.
    """

    cl = pd.Series(labels.astype(str))
    root = None

    # Named candidates first
    for candid in ROOT_CANDIDATES:
        if (cl == candid).any():
            root = candid
            break

    # Then numeric minimum if possible
    if root is None:
        num = pd.to_numeric(cl, errors="coerce")
        if np.isfinite(num).any():
            root = str(np.nanmin(num.values))

    # Fallback: most frequent label
    if root is None and len(cl) > 0:
        vc = cl.value_counts()
        if len(vc) >= 1:
            root = str(vc.index[0])

    if root is None:
        raise RuntimeError("Palantir: cannot determine root cluster from obs")

    return root


def load_counts_matrix(adata):
    if "counts" in adata.layers:
        x = adata.layers["counts"]
    elif adata.raw is not None:
        x = adata.raw.X
    else:
        x = adata.X
    if hasattr(x, "toarray"):
        x = x.toarray()
    return np.asarray(x, dtype=np.float64)


def main():
    if len(sys.argv) != 3:
        print("Usage: run_palantir.py <h5ad> <output_h5ad>", file=sys.stderr)
        sys.exit(2)
    inpath, outpath = sys.argv[1], sys.argv[2]

    import palantir

    ad = sc.read_h5ad(inpath)
    canonical_cell_ids = ad.obs_names.astype(str).to_numpy(copy=True)
    if canonical_cell_ids.size != ad.n_obs or pd.Index(canonical_cell_ids).has_duplicates:
        raise RuntimeError("Palantir: canonical obs_names must be complete and unique")

    # No subsampling: Palantir runs on the full cell set (official workflow).
    ad.X = load_counts_matrix(ad)
    if not np.any(ad.X > 0):
        raise RuntimeError("Palantir: empty count matrix")

    sc.pp.filter_genes(ad, min_cells=3)
    sc.pp.normalize_total(ad)
    sc.pp.log1p(ad)
    sc.pp.highly_variable_genes(ad, n_top_genes=min(2000, ad.n_vars - 1))
    ad = ad[:, ad.var.highly_variable].copy()
    sc.pp.pca(ad, n_comps=min(30, ad.n_vars - 1, ad.n_obs - 1))

    palantir.utils.run_diffusion_maps(ad, n_components=5)
    palantir.utils.determine_multiscale_space(ad)

    track = benchmark_track()
    if track == "common_root":
        start_cell = pick_common_root_cell(ad)
    else:
        celltype_column, labels = pick_cluster_and_labels(ad)
        if celltype_column is None:
            celltype_column, labels = attach_fallback_clusters(ad)
        if track == "annotation_free":
            start_cell = str(ad.obs_names[expression_root_position(ad, labels.to_numpy())])
        else:
            root_type = determine_root_type(labels)
            start_cell = palantir.utils.early_cell(
                ad,
                root_type,
                celltype_column=celltype_column,
                fallback_seed=CONFIG["seed"],
            )

    run_fn = getattr(palantir.core, "run_palantir", None) or palantir.utils.run_palantir
    kwargs = {"num_waypoints": min(500, max(100, ad.n_obs // 2))}
    pr_res = run_fn(ad, start_cell, **kwargs)

    if not np.array_equal(ad.obs_names.astype(str).to_numpy(), canonical_cell_ids):
        raise RuntimeError("Palantir: preprocessing changed canonical cell order")
    raw_pt = pr_res.pseudotime
    if isinstance(raw_pt, pd.Series):
        pt_index = raw_pt.index.astype(str)
        if pt_index.has_duplicates:
            raise RuntimeError("Palantir: pseudotime contains duplicate cell IDs")
        pal_pt = pd.to_numeric(
            pd.Series(raw_pt.to_numpy(), index=pt_index).reindex(canonical_cell_ids),
            errors="coerce",
        ).to_numpy(dtype=float)
    else:
        pal_pt = np.asarray(raw_pt, dtype=float).reshape(-1)
    if pal_pt.size != canonical_cell_ids.size:
        raise RuntimeError(
            f"Palantir: pseudotime length {pal_pt.size} != canonical cells {canonical_cell_ids.size}"
        )
    if not np.isfinite(pal_pt).all() or np.ptp(pal_pt) <= 0:
        raise RuntimeError("Palantir: pal_pt must be finite for every cell and have non-zero range")

    ad.obs["pal_pt"] = pal_pt
    for c in pr_res.branch_probs.columns:
        ad.obs[f"pal_br_{c}"] = pr_res.branch_probs[c].values
    ad.write_h5ad(outpath)

    written = sc.read_h5ad(outpath, backed="r")
    written_cell_ids = written.obs_names.astype(str).to_numpy(copy=True)
    written_pt = pd.to_numeric(written.obs["pal_pt"], errors="coerce").to_numpy(dtype=float)
    written.file.close()
    if (
        not np.array_equal(written_cell_ids, canonical_cell_ids)
        or pd.Index(written_cell_ids).has_duplicates
        or written_pt.size != canonical_cell_ids.size
        or not np.isfinite(written_pt).all()
        or np.ptp(written_pt) <= 0
    ):
        raise RuntimeError("Palantir: written H5AD failed the full-cell pal_pt contract")
    print(f"OK: Palantir {ad.n_obs} cells [{track}] -> {outpath}")


if __name__ == "__main__":
    main()
