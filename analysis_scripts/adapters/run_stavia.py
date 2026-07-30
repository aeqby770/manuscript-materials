#!/usr/bin/env python
"""StaVia/VIA adapter for canonical h5ad input.

Usage: python run_stavia.py <canonical.h5ad> <output_pseudotime.csv>
"""
import os
import shutil
import sys
import types

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
import _mem_guard
from _benchmark_contract import (
    benchmark_track,
    common_root_position,
    expression_root_position,
    pick_cluster_column,
)
_mem_guard.cap_self(55)

VENDOR_ROOT = os.environ.get(
    "LLR_VENDOR_ROOT",
    r"F:\biotrade\LineaLinkR_benchmark_assets\02_vendor_algorithms",
)
GIT_ROOT = os.path.dirname(os.path.dirname(shutil.which("git")))
CONFIG = {
    "paths": {
        "vendor": os.path.join(VENDOR_ROOT, "StaVia"),
        "mingw_bin": os.path.join(GIT_ROOT, "mingw64", "bin"),
    },
    "threads": os.cpu_count(),
    "seed": 12345,
    "plot": {
        "dpi": 300,
        "font": "Arial",
        "fontsize": 12,
        "palette": "viridis",
    },
    "n_pcs": 30,
    "knn": 30,
}

os.environ["OMP_NUM_THREADS"] = str(CONFIG["threads"])
np.random.seed(CONFIG["seed"])
np.bool8 = np.bool_
np.object0 = np.object_
np.int0 = np.intp
np.uint0 = np.uintp
np.float_ = np.float64
np.longfloat = np.longdouble
np.complex_ = np.complex128
np.cfloat = np.complex128
np.singlecomplex = np.complex64
np.clongfloat = np.clongdouble
np.longcomplex = np.clongdouble
np.string_ = np.bytes_
np.bytes0 = np.bytes_
np.unicode_ = np.str_
np.str0 = np.str_
np.void0 = np.void
if os.path.isdir(CONFIG["paths"]["mingw_bin"]):
    MINGW_DLL_DIR = os.add_dll_directory(CONFIG["paths"]["mingw_bin"])
sys.path.insert(0, CONFIG["paths"]["vendor"])
pyvia_package = types.ModuleType("pyVIA")
pyvia_package.__path__ = [os.path.join(CONFIG["paths"]["vendor"], "VIA")]
sys.modules["pyVIA"] = pyvia_package

import pyVIA.utils_via as via_utils  # noqa: E402

_scipy_csr_matrix = via_utils.csr_matrix
via_utils.csr_matrix = lambda arg1=None, *args, **kwargs: _scipy_csr_matrix((arg1[0], tuple(arg1[1])), *args, **kwargs) if isinstance(arg1, tuple) and len(arg1) == 2 and isinstance(arg1[1], zip) else _scipy_csr_matrix(arg1, *args, **kwargs)
import pyVIA.core as via  # noqa: E402


def main():
    if len(sys.argv) != 3:
        print("Usage: run_stavia.py <canonical.h5ad> <output_pseudotime.csv>", file=sys.stderr)
        sys.exit(2)

    inpath, outpath = sys.argv[1], sys.argv[2]
    adata = sc.read_h5ad(inpath)
    canonical_cell_ids = adata.obs_names.astype(str).to_numpy(copy=True)
    if canonical_cell_ids.size != adata.n_obs or pd.Index(canonical_cell_ids).has_duplicates:
        raise RuntimeError("StaVia: canonical obs_names must be complete and unique")
    _mem_guard.require_under(
        _mem_guard.dense_gb(adata.n_obs, min(adata.n_vars, 3000), copies=1),
        55,
        "StaVia",
    )

    if "lognorm" in adata.layers:
        adata.X = adata.layers["lognorm"].copy()
    elif "counts" in adata.layers:
        adata.X = adata.layers["counts"].copy()
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
    else:
        adata.X = adata.X.copy()

    if sparse.issparse(adata.X):
        adata.X.data = np.nan_to_num(
            adata.X.data, nan=0.0, posinf=0.0, neginf=0.0
        )
        adata.X.eliminate_zeros()
    else:
        adata.X = np.nan_to_num(
            np.asarray(adata.X), nan=0.0, posinf=0.0, neginf=0.0
        )

    n_pcs = min(CONFIG["n_pcs"], adata.n_obs - 1, adata.n_vars - 1)
    sc.pp.pca(adata, n_comps=n_pcs, random_state=CONFIG["seed"])

    label_col = pick_cluster_column(adata.obs)
    if label_col is not None:
        labels = adata.obs[label_col].astype(str).tolist()
    else:
        # Fallback: PCA + Leiden clustering
        sc.pp.neighbors(adata, use_rep="X_pca", n_neighbors=min(30, adata.n_obs - 1))
        sc.tl.leiden(adata, resolution=1.0, key_added="benchmark_clusters")
        labels = adata.obs["benchmark_clusters"].astype(str).tolist()
        label_col = "benchmark_clusters"
        print(f"StaVia: using fallback Leiden clusters")

    # Root cell selection: consistent with DPT/PHLOWER — pick root cluster by
    # label semantics, then select PCA-centroid cell within that cluster.
    # VIA root_user accepts cell indices (int), not group labels.
    labels_arr = np.array(labels)
    root_label = "common_root"
    if benchmark_track() == "common_root":
        root_idx = common_root_position(adata.obs_names)
    elif benchmark_track() == "annotation_free":
        root_label = "expression_only"
        root_idx = expression_root_position(adata, labels_arr)
    else:
        root_candidates = ("0", "root", "Root", "STEM", "stem", "branch_0", "early", "Early")
        matches = [candidate for candidate in root_candidates if (labels_arr == candidate).any()]
        root_label = matches[0] if matches else None
        if root_label is None:
            numeric_labels = pd.to_numeric(pd.Series(labels_arr), errors="coerce")
            root_label = str(np.nanmin(numeric_labels.values)) if np.isfinite(numeric_labels).any() else None
        if root_label is not None:
            root_cells = np.where(labels_arr == root_label)[0]
            pca = adata.obsm["X_pca"][root_cells, :]
            centroid = pca.mean(axis=0)
            root_idx = int(root_cells[np.argmin(np.linalg.norm(pca - centroid, axis=1))])
        else:
            root_idx = 0
            print("StaVia: no root cluster determined, using cell_idx=0")
    print(f"StaVia: root_label='{root_label}', root_user=[{root_idx}]")

    model = via.VIA(
        adata.obsm["X_pca"][:, :n_pcs],
        labels,
        knn=min(CONFIG["knn"], adata.n_obs - 1),
        random_seed=CONFIG["seed"],
        num_threads=CONFIG["threads"],
        root_user=[root_idx],
    )
    model.run_VIA()
    pseudotime = np.asarray(model.single_cell_pt_markov, dtype=float).reshape(-1)
    if not np.array_equal(adata.obs_names.astype(str).to_numpy(), canonical_cell_ids):
        raise RuntimeError("StaVia: preprocessing changed canonical cell order")
    if pseudotime.size != canonical_cell_ids.size:
        raise RuntimeError(
            f"StaVia: pseudotime length {pseudotime.size} != canonical cells {canonical_cell_ids.size}"
        )
    if not np.isfinite(pseudotime).all() or np.ptp(pseudotime) <= 0:
        raise RuntimeError("StaVia: pseudotime must be finite for every cell and have non-zero range")

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
        raise RuntimeError("StaVia: written CSV failed the full-cell pseudotime contract")

    print(f"StaVia OK: {adata.n_obs} cells -> {outpath}", flush=True)


if __name__ == "__main__":
    main()
    sys.exit(0)
