# -*- coding: utf-8 -*-
"""Task B: ligand→receptor peak-bin lag for temporal-only headline pairs."""
from __future__ import annotations


import sys
from pathlib import Path
_LLR_CFG = Path(__file__).resolve().parents[1] / "config"
if str(_LLR_CFG) not in sys.path:
    sys.path.insert(0, str(_LLR_CFG))
from llr_paths import work_dir as _llr_work_dir  # noqa: E402
_LLR_WORK = _llr_work_dir()

import os

import numpy as np
import pandas as pd
import scipy.io
from scipy import sparse

CONFIG = {
    "paths": {
        "workspace": _LLR_WORK / "cache/linealinkr_workspace",
        "temporal_only": _LLR_WORK / "results/paper_tables/supplement/tables/temporal_only_validation",
        "out_dir": _LLR_WORK / "results/paper_tables/supplement/tables/temporal_only_orthogonal_validation",
    },
    "threads": os.cpu_count() or 1,
    "seed": 12345,
    "n_bins": 12,
    "n_perm": 500,
    "datasets": ["pancreas", "dentategyrus"],
}

os.environ["OMP_NUM_THREADS"] = str(CONFIG["threads"])
np.random.seed(CONFIG["seed"])
CONFIG["paths"]["out_dir"].mkdir(parents=True, exist_ok=True)


def load_expr(ds: str):
    extract = CONFIG["paths"]["workspace"] / ds / "extract"
    genes = (extract / f"{ds}_genes.txt").read_text(encoding="utf-8").splitlines()
    cells = (extract / f"{ds}_cells.txt").read_text(encoding="utf-8").splitlines()
    mat = scipy.io.mmread(extract / f"{ds}_matrix.mtx")
    if sparse.issparse(mat):
        mat = mat.tocsr()
    else:
        mat = sparse.csr_matrix(mat)
    if mat.shape[0] != len(genes) and mat.shape[1] == len(genes):
        mat = mat.T.tocsr()
    gene_index = {g: i for i, g in enumerate(genes)}
    return mat, genes, cells, gene_index


def load_pt(ds: str) -> pd.Series:
    p = CONFIG["paths"]["out_dir"] / f"{ds}_llr_pseudotime.csv"
    df = pd.read_csv(p)
    return pd.Series(df["pseudotime"].to_numpy(), index=df["cell"].astype(str))


def bin_means(values: np.ndarray, pt: np.ndarray, n_bins: int) -> np.ndarray:
    ok = np.isfinite(values) & np.isfinite(pt)
    x = values[ok]
    t = pt[ok]
    br = np.unique(np.quantile(t, np.linspace(0, 1, n_bins + 1)))
    out = np.full(n_bins, np.nan)
    if br.size < 3:
        return out
    bins = np.digitize(t, br[1:-1], right=True)
    for b in range(n_bins):
        idx = bins == b
        if np.any(idx):
            out[b] = float(np.mean(x[idx]))
    return out


def peak_bin(v: np.ndarray) -> int:
    if not np.any(np.isfinite(v)):
        return -1
    vv = np.where(np.isfinite(v), v, -np.inf)
    return int(np.argmax(vv)) + 1


def cross_corr_lag(x: np.ndarray, y: np.ndarray, max_lag: int = 5):
    ok = np.isfinite(x) & np.isfinite(y)
    x = x[ok]
    y = y[ok]
    if x.size < 4:
        return np.nan, np.nan
    x = (x - x.mean()) / (x.std() + 1e-12)
    y = (y - y.mean()) / (y.std() + 1e-12)
    lags = np.arange(-max_lag, max_lag + 1)
    cors = []
    # lag L > 0 means ligand leads receptor: corr(ligand[t], receptor[t+L])
    for L in lags:
        if L > 0:
            xs, ys = x[: x.size - L], y[L:]
        elif L < 0:
            xs, ys = x[-L:], y[: y.size + L]
        else:
            xs, ys = x, y
        if xs.size < 3 or xs.std() == 0 or ys.std() == 0:
            cors.append(np.nan)
        else:
            cors.append(float(np.corrcoef(xs, ys)[0, 1]))
    cors = np.asarray(cors, dtype=float)
    i = int(np.nanargmax(cors))
    return int(lags[i]), float(cors[i])


def perm_p_lag_positive(x: np.ndarray, y: np.ndarray, n_perm: int, seed: int) -> float:
    rng = np.random.default_rng(seed)
    obs = peak_bin(y) - peak_bin(x)
    if obs == 0 and (peak_bin(x) < 0 or peak_bin(y) < 0):
        return np.nan
    nulls = []
    for _ in range(n_perm):
        yp = rng.permutation(y)
        nulls.append(peak_bin(yp) - peak_bin(x))
    nulls = np.asarray(nulls)
    return float(np.mean(nulls >= obs))


def gene_profile(mat, gene_index, gene_str: str):
    genes = [g for g in gene_str.replace(";", ",").replace("|", ",").replace("_", ",").split(",") if g]
    genes = [g.strip() for g in genes if g.strip()]
    idxs = [gene_index[g] for g in genes if g in gene_index]
    if not idxs:
        return None
    if len(idxs) == 1:
        return np.asarray(mat[idxs[0], :].todense()).ravel()
    return np.asarray(mat[idxs, :].mean(axis=0)).ravel()


rows = []
profiles = []

for ds in CONFIG["datasets"]:
    to = pd.read_csv(CONFIG["paths"]["temporal_only"] / f"{ds}_temporal_only_vs_external.csv")
    mat, genes, cells, gene_index = load_expr(ds)
    pt = load_pt(ds)
    cell_idx = {c: i for i, c in enumerate(cells)}
    common = [c for c in cells if c in pt.index]
    col_sel = np.array([cell_idx[c] for c in common], dtype=int)
    pt_arr = pt.loc[common].to_numpy(dtype=float)
    mat_c = mat[:, col_sel]

    for i, r in to.iterrows():
        lig = str(r["ligand_genes"])
        rec = str(r["receptor_genes"])
        pair = str(r["pair_name"])
        lig_v = gene_profile(mat_c, gene_index, lig)
        rec_v = gene_profile(mat_c, gene_index, rec)
        if lig_v is None or rec_v is None:
            rows.append(
                {
                    "dataset": ds,
                    "pair_name": pair,
                    "ligand_genes": lig,
                    "receptor_genes": rec,
                    "family": r.get("family"),
                    "temporal_mean_score": r.get("temporal_mean_score"),
                    "ligand_peak_bin": np.nan,
                    "receptor_peak_bin": np.nan,
                    "lag_bins": np.nan,
                    "lag_positive": np.nan,
                    "best_xcorr_lag": np.nan,
                    "best_xcorr": np.nan,
                    "perm_p_lag_positive": np.nan,
                    "status": "gene_missing_in_matrix",
                }
            )
            continue
        lig_bin = bin_means(lig_v, pt_arr, CONFIG["n_bins"])
        rec_bin = bin_means(rec_v, pt_arr, CONFIG["n_bins"])
        lp = peak_bin(lig_bin)
        rp = peak_bin(rec_bin)
        lag = rp - lp if lp > 0 and rp > 0 else np.nan
        xc_lag, xc_cor = cross_corr_lag(lig_bin, rec_bin, max_lag=5)
        pp = perm_p_lag_positive(lig_bin, rec_bin, CONFIG["n_perm"], CONFIG["seed"] + int(i))
        rows.append(
            {
                "dataset": ds,
                "pair_name": pair,
                "ligand_genes": lig,
                "receptor_genes": rec,
                "family": r.get("family"),
                "temporal_mean_score": r.get("temporal_mean_score"),
                "ligand_peak_bin": lp,
                "receptor_peak_bin": rp,
                "lag_bins": lag,
                "lag_positive": bool(np.isfinite(lag) and lag > 0),
                "best_xcorr_lag": xc_lag,
                "best_xcorr": xc_cor,
                "perm_p_lag_positive": pp,
                "status": "ok",
            }
        )
        for b in range(CONFIG["n_bins"]):
            profiles.append(
                {
                    "dataset": ds,
                    "pair_name": pair,
                    "bin": b + 1,
                    "ligand_mean": lig_bin[b],
                    "receptor_mean": rec_bin[b],
                }
            )

lag_dt = pd.DataFrame(rows)
prof_dt = pd.DataFrame(profiles)
lag_dt.to_csv(CONFIG["paths"]["out_dir"] / "temporal_only_ligand_receptor_lag.csv", index=False)
prof_dt.to_csv(CONFIG["paths"]["out_dir"] / "temporal_only_bin_expression_profiles.csv", index=False)

ok = lag_dt[lag_dt["status"] == "ok"]
summary = (
    ok.groupby("dataset")
    .agg(
        n_pairs=("pair_name", "size"),
        n_lag_positive=("lag_positive", "sum"),
        frac_lag_positive=("lag_positive", "mean"),
        median_lag=("lag_bins", "median"),
        n_xcorr_lag_positive=("best_xcorr_lag", lambda s: int(np.sum(np.asarray(s) > 0))),
        median_perm_p=("perm_p_lag_positive", "median"),
    )
    .reset_index()
)
summary.to_csv(CONFIG["paths"]["out_dir"] / "temporal_only_lag_summary_by_dataset.csv", index=False)
print(lag_dt.to_string(index=False))
print(summary.to_string(index=False))
