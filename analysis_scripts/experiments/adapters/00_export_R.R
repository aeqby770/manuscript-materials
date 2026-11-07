# Export h5ad -> R/Python-friendly artifacts under cache/r_export/<dataset>/.
# Used by R-package adapters (TSCAN, monocle3, slingshot, SCORPIUS, …).
# Python adapters typically read h5ad directly via DatasetPath().
# See docs/COMPETITOR_INPUT_CONTRACTS.md for which artifact each algorithm needs.
suppressPackageStartupMessages(library(data.table))

# path bootstrap (relative; set LLR_ASSETS or LLR_PIPELINE to override)
.llr_cmd <- commandArgs(trailingOnly = FALSE)
.llr_file <- sub("^--file=", "", grep("^--file=", .llr_cmd, value = TRUE)[1])
if (!length(.llr_file) || is.na(.llr_file) || !nzchar(.llr_file)) {
  stop("Run this script via Rscript so --file= is set")
}
.llr_here <- dirname(normalizePath(.llr_file, winslash = "/", mustWork = TRUE))
WORK_ROOT <- normalizePath(file.path(.llr_here,"../.."),winslash="/",mustWork=TRUE)
Sys.setenv(LLR_WORK_ROOT=WORK_ROOT,LLR_PIPELINE=file.path(WORK_ROOT,"experiments"),LLR_ASSETS=WORK_ROOT)
source(file.path(.llr_here, "../config/llr_bootstrap.R"), local = FALSE)

llr_source_config()
args <- commandArgs(trailingOnly = TRUE)
ds_name <- args[1]
if (is.na(ds_name)) stop("Usage: Rscript 00_export_R.R <dataset_name>")
outdir <- file.path(CONFIG$dir$r_export, ds_name)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
h5ad <- DatasetPath(ds_name)

sp <- CONFIG$dataset_species[[ds_name]]
if (is.null(sp)) sp <- "human"
input_scale <- unname(CONFIG$dataset_lr_input_scale[ds_name])
input_scale[is.na(input_scale)] <- "counts"

py <-
'import atexit
import json
import os
import shutil
import tempfile
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from scipy.io import mmread, mmwrite

CONFIG = {
    "paths": {
        "h5ad": r"__H5AD__", "output": r"__OUTDIR__",
        "exclusions": r"__EXCLUSIONS__",
    },
    "dataset": "__DATASET__",
    "threads": int(os.cpu_count() or 1),
    "seed": 12345,
    "plot": {},
    "species": "__SPECIES__",
    "input_scale": "__INPUT_SCALE__",
    "contract": "linealinkr064_r_export_sparse",
    "anchor_count": 64,
    "rtol": 1e-6,
    "atol": 1e-10,
    "sum_atol": 1e-3,
}
os.environ["OMP_NUM_THREADS"] = str(CONFIG["threads"])
np.random.seed(CONFIG["seed"])
exclusion_table = pd.read_csv(CONFIG["paths"]["exclusions"], encoding="utf-8")
excluded_genes = exclusion_table.loc[
    exclusion_table["dataset"].eq(CONFIG["dataset"]), "gene"
].tolist()

source_path = str(Path(CONFIG["paths"]["h5ad"]).resolve())
output_dir = Path(CONFIG["paths"]["output"])
output_dir.mkdir(parents=True, exist_ok=True)
required_paths = {
    "counts": output_dir / "expr.mtx",
    "lognorm": output_dir / "expr_lognorm.mtx",
    "meta": output_dir / "meta.csv",
    "genes": output_dir / "genes.txt",
    "cells": output_dir / "cells.txt",
    "embedding": output_dir / "umap.csv",
    "info": output_dir / "info.json",
}
receipt_path = output_dir / "decoded_validation_receipt.json"
cache_present = (
    required_paths["counts"].is_file()
    and required_paths["lognorm"].is_file()
    and required_paths["meta"].is_file()
    and required_paths["genes"].is_file()
    and required_paths["cells"].is_file()
    and required_paths["embedding"].is_file()
    and required_paths["info"].is_file()
)
receipt_valid = False
if receipt_path.is_file() and cache_present:
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    source_stat = Path(source_path).stat()
    artifact_state = {
        "counts": {"size": required_paths["counts"].stat().st_size, "mtime_ns": required_paths["counts"].stat().st_mtime_ns},
        "lognorm": {"size": required_paths["lognorm"].stat().st_size, "mtime_ns": required_paths["lognorm"].stat().st_mtime_ns},
        "meta": {"size": required_paths["meta"].stat().st_size, "mtime_ns": required_paths["meta"].stat().st_mtime_ns},
        "genes": {"size": required_paths["genes"].stat().st_size, "mtime_ns": required_paths["genes"].stat().st_mtime_ns},
        "cells": {"size": required_paths["cells"].stat().st_size, "mtime_ns": required_paths["cells"].stat().st_mtime_ns},
        "embedding": {"size": required_paths["embedding"].stat().st_size, "mtime_ns": required_paths["embedding"].stat().st_mtime_ns},
        "info": {"size": required_paths["info"].stat().st_size, "mtime_ns": required_paths["info"].stat().st_mtime_ns},
    }
    receipt_valid = bool(
        receipt.get("receipt_contract") == "linealinkr064_r_export_decoded_validation_v1"
        and receipt.get("source") == {
            "path": source_path,
            "size": source_stat.st_size,
            "mtime_ns": source_stat.st_mtime_ns,
        }
        and receipt.get("artifacts") == artifact_state
        and receipt.get("excluded_genes", []) == excluded_genes
        and receipt.get("input_scale") == CONFIG["input_scale"]
    )
if receipt_valid:
    print("CACHED_VALIDATED")
    raise SystemExit(0)

adata = ad.read_h5ad(source_path)
adata = adata[:, ~adata.var_names.isin(excluded_genes)].copy()
cell_ids = np.asarray(adata.obs_names.astype(str), dtype=str)
gene_ids = np.asarray(adata.var_names.astype(str), dtype=str)

counts_obs_var = adata.layers["counts"] if "counts" in adata.layers else adata.X
lognorm_obs_var = adata.layers["lognorm"] if "lognorm" in adata.layers else adata.X
counts = sparse.csr_matrix(counts_obs_var).transpose().tocsr()
lognorm = sparse.csr_matrix(lognorm_obs_var).transpose().tocsr()
if CONFIG["input_scale"] == "counts" and (excluded_genes or "lognorm" not in adata.layers):
    library_size = np.asarray(counts.sum(axis=0)).ravel()
    lognorm = counts.multiply(10000 / library_size).tocsr()
    lognorm.data = np.log1p(lognorm.data)
    lognorm_obs_var = lognorm.transpose().tocsr()
counts.sum_duplicates()
counts.eliminate_zeros()
counts.sort_indices()
lognorm.sum_duplicates()
lognorm.eliminate_zeros()
lognorm.sort_indices()
gene_weights = np.ones(counts.shape[0], dtype=np.float64)
cell_weights = np.ones(counts.shape[1], dtype=np.float64)
counts_cell_sums = counts.transpose() @ gene_weights
counts_gene_sums = counts @ cell_weights
lognorm_cell_sums = lognorm.transpose() @ gene_weights
lognorm_gene_sums = lognorm @ cell_weights

obs = adata.obs.copy()
obs = obs.drop(columns=["cell_id"], errors="ignore")
obs.insert(0, "cell_id", cell_ids)
meta_text = obs.to_csv(index=False, lineterminator="\\n")

cache_readable = False
if cache_present:
    try:
        exported_counts = sparse.csr_matrix(mmread(required_paths["counts"])).tocsr()
        exported_lognorm = sparse.csr_matrix(mmread(required_paths["lognorm"])).tocsr()
        exported_counts.sum_duplicates()
        exported_counts.eliminate_zeros()
        exported_counts.sort_indices()
        exported_lognorm.sum_duplicates()
        exported_lognorm.eliminate_zeros()
        exported_lognorm.sort_indices()
        exported_cells = pd.read_csv(
            required_paths["cells"], header=None, dtype=str, keep_default_na=False
        ).iloc[:, 0].to_numpy(dtype=str)
        exported_genes = pd.read_csv(
            required_paths["genes"], header=None, dtype=str, keep_default_na=False
        ).iloc[:, 0].to_numpy(dtype=str)
        exported_meta_text = required_paths["meta"].read_text(encoding="utf-8")
        exported_meta = pd.read_csv(
            required_paths["meta"], dtype=str, keep_default_na=False
        )
        exported_embedding = np.atleast_2d(
            np.loadtxt(required_paths["embedding"], delimiter=",")
        )
        exported_info = json.loads(required_paths["info"].read_text(encoding="utf-8"))
        cache_readable = True
    except (OSError, ValueError, TypeError, KeyError, IndexError, json.JSONDecodeError):
        cache_readable = False

count_positions = np.unique(np.linspace(
    0, max(counts.nnz - 1, 0), min(CONFIG["anchor_count"], max(counts.nnz, 1)), dtype=np.int64
)) if counts.nnz else np.asarray([], dtype=np.int64)
count_rows = np.searchsorted(counts.indptr[1:], count_positions, side="right")
count_columns = counts.indices[count_positions] if counts.nnz else np.asarray([], dtype=np.int64)
count_anchor_readable = cache_readable and exported_counts.shape == counts.shape
count_exported_anchors = np.asarray(
    exported_counts[count_rows, count_columns]
).ravel() if count_anchor_readable and counts.nnz else np.asarray([], dtype=np.float64)
count_source_anchors = counts.data[count_positions] if counts.nnz else np.asarray([], dtype=np.float64)

lognorm_positions = np.unique(np.linspace(
    0, max(lognorm.nnz - 1, 0), min(CONFIG["anchor_count"], max(lognorm.nnz, 1)), dtype=np.int64
)) if lognorm.nnz else np.asarray([], dtype=np.int64)
lognorm_rows = np.searchsorted(lognorm.indptr[1:], lognorm_positions, side="right")
lognorm_columns = lognorm.indices[lognorm_positions] if lognorm.nnz else np.asarray([], dtype=np.int64)
lognorm_anchor_readable = cache_readable and exported_lognorm.shape == lognorm.shape
lognorm_exported_anchors = np.asarray(
    exported_lognorm[lognorm_rows, lognorm_columns]
).ravel() if lognorm_anchor_readable and lognorm.nnz else np.asarray([], dtype=np.float64)
lognorm_source_anchors = lognorm.data[lognorm_positions] if lognorm.nnz else np.asarray([], dtype=np.float64)

cache_valid = bool(
    cache_readable
    and np.array_equal(exported_cells, cell_ids)
    and np.array_equal(exported_genes, gene_ids)
    and exported_meta_text.replace("\\r\\n", "\\n") == meta_text
    and "cell_id" in exported_meta.columns
    and np.array_equal(exported_meta["cell_id"].to_numpy(dtype=str), cell_ids)
    and exported_counts.shape == counts.shape
    and exported_lognorm.shape == lognorm.shape
    and exported_counts.nnz == counts.nnz
    and exported_lognorm.nnz == lognorm.nnz
    and np.allclose(exported_counts.transpose() @ gene_weights, counts_cell_sums, rtol=CONFIG["rtol"], atol=CONFIG["sum_atol"])
    and np.allclose(exported_counts @ cell_weights, counts_gene_sums, rtol=CONFIG["rtol"], atol=CONFIG["sum_atol"])
    and np.allclose(exported_lognorm.transpose() @ gene_weights, lognorm_cell_sums, rtol=CONFIG["rtol"], atol=CONFIG["sum_atol"])
    and np.allclose(exported_lognorm @ cell_weights, lognorm_gene_sums, rtol=CONFIG["rtol"], atol=CONFIG["sum_atol"])
    and np.allclose(count_exported_anchors, count_source_anchors, rtol=CONFIG["rtol"], atol=CONFIG["atol"])
    and np.allclose(lognorm_exported_anchors, lognorm_source_anchors, rtol=CONFIG["rtol"], atol=CONFIG["atol"])
    and exported_embedding.shape == (adata.n_obs, 2)
    and np.isfinite(exported_embedding).all()
    and exported_info.get("content_contract") == CONFIG["contract"]
    and exported_info.get("source_path") == source_path
    and exported_info.get("species") == CONFIG["species"]
    and exported_info.get("n_obs") == int(adata.n_obs)
    and exported_info.get("n_vars") == int(adata.n_vars)
    and exported_info.get("excluded_genes", []) == excluded_genes
    and exported_info.get("input_scale") == CONFIG["input_scale"]
)

if not cache_valid:
    staging_dir = Path(tempfile.mkdtemp(prefix=".r_export_staging_", dir=output_dir))
    staging_cleanup = atexit.register(shutil.rmtree, staging_dir, ignore_errors=True)
    mmwrite(staging_dir / "expr.mtx", counts)
    mmwrite(staging_dir / "expr_lognorm.mtx", lognorm)
    (staging_dir / "meta.csv").write_text(meta_text, encoding="utf-8")
    np.savetxt(staging_dir / "genes.txt", gene_ids, fmt="%s")
    np.savetxt(staging_dir / "cells.txt", cell_ids, fmt="%s")

    embedding = None
    if not excluded_genes and "X_umap" in adata.obsm:
        candidate = np.asarray(adata.obsm["X_umap"][:, :2], dtype=np.float64)
        embedding = candidate if candidate.shape == (adata.n_obs, 2) and np.isfinite(candidate).all() and np.ptp(candidate, axis=0).max() > 0 else None
    if embedding is None and not excluded_genes and "X_pca" in adata.obsm:
        candidate = np.asarray(adata.obsm["X_pca"][:, :2], dtype=np.float64)
        embedding = candidate if candidate.shape == (adata.n_obs, 2) and np.isfinite(candidate).all() and np.ptp(candidate, axis=0).max() > 0 else None
    if embedding is None:
        reduced = ad.AnnData(X=sparse.csr_matrix(lognorm_obs_var).copy())
        reduced.var_names = adata.var_names.copy()
        sc.pp.filter_genes(reduced, min_cells=3)
        n_pcs = int(min(30, max(2, reduced.n_vars - 1), max(2, reduced.n_obs - 1)))
        if reduced.n_vars > 500:
            sc.pp.highly_variable_genes(
                reduced,
                n_top_genes=min(2000, reduced.n_vars - 1),
                flavor="seurat_v3",
            )
            reduced = reduced[:, reduced.var.highly_variable].copy() if int(reduced.var.highly_variable.sum()) > 2 else reduced
        sc.pp.pca(reduced, n_comps=n_pcs, random_state=CONFIG["seed"])
        embedding = np.asarray(reduced.obsm["X_pca"][:, :2], dtype=np.float64)
    np.savetxt(staging_dir / "umap.csv", embedding, delimiter=",")
    (staging_dir / "info.json").write_text(
        json.dumps(
            {
                "content_contract": CONFIG["contract"],
                "source_path": source_path,
                "species": CONFIG["species"],
                "input_scale": CONFIG["input_scale"],
                "n_obs": int(adata.n_obs),
                "n_vars": int(adata.n_vars),
                "excluded_genes": excluded_genes,
            },
            ensure_ascii=False,
            indent=2,
        ) + "\\n",
        encoding="utf-8",
    )

    rebuilt_counts = sparse.csr_matrix(mmread(staging_dir / "expr.mtx")).tocsr()
    rebuilt_lognorm = sparse.csr_matrix(mmread(staging_dir / "expr_lognorm.mtx")).tocsr()
    rebuilt_counts.sum_duplicates()
    rebuilt_counts.eliminate_zeros()
    rebuilt_counts.sort_indices()
    rebuilt_lognorm.sum_duplicates()
    rebuilt_lognorm.eliminate_zeros()
    rebuilt_lognorm.sort_indices()
    rebuilt_cells = pd.read_csv(
        staging_dir / "cells.txt", header=None, dtype=str, keep_default_na=False
    ).iloc[:, 0].to_numpy(dtype=str)
    rebuilt_genes = pd.read_csv(
        staging_dir / "genes.txt", header=None, dtype=str, keep_default_na=False
    ).iloc[:, 0].to_numpy(dtype=str)
    rebuilt_meta = pd.read_csv(staging_dir / "meta.csv", dtype=str, keep_default_na=False)
    rebuilt_embedding = np.atleast_2d(np.loadtxt(staging_dir / "umap.csv", delimiter=","))
    rebuilt_count_anchors = np.asarray(
        rebuilt_counts[count_rows, count_columns]
    ).ravel() if counts.nnz else np.asarray([], dtype=np.float64)
    rebuilt_lognorm_anchors = np.asarray(
        rebuilt_lognorm[lognorm_rows, lognorm_columns]
    ).ravel() if lognorm.nnz else np.asarray([], dtype=np.float64)
    rebuilt_checks = {
        "cell_ids": np.array_equal(rebuilt_cells, cell_ids),
        "gene_ids": np.array_equal(rebuilt_genes, gene_ids),
        "meta_cell_order": np.array_equal(rebuilt_meta["cell_id"].to_numpy(dtype=str), cell_ids),
        "counts_shape": rebuilt_counts.shape == counts.shape,
        "lognorm_shape": rebuilt_lognorm.shape == lognorm.shape,
        "counts_nnz": rebuilt_counts.nnz == counts.nnz,
        "lognorm_nnz": rebuilt_lognorm.nnz == lognorm.nnz,
        "counts_cell_sums": np.allclose(rebuilt_counts.transpose() @ gene_weights, counts_cell_sums, rtol=CONFIG["rtol"], atol=CONFIG["sum_atol"]),
        "counts_gene_sums": np.allclose(rebuilt_counts @ cell_weights, counts_gene_sums, rtol=CONFIG["rtol"], atol=CONFIG["sum_atol"]),
        "lognorm_cell_sums": np.allclose(rebuilt_lognorm.transpose() @ gene_weights, lognorm_cell_sums, rtol=CONFIG["rtol"], atol=CONFIG["sum_atol"]),
        "lognorm_gene_sums": np.allclose(rebuilt_lognorm @ cell_weights, lognorm_gene_sums, rtol=CONFIG["rtol"], atol=CONFIG["sum_atol"]),
        "counts_anchors": np.allclose(rebuilt_count_anchors, count_source_anchors, rtol=CONFIG["rtol"], atol=CONFIG["atol"]),
        "lognorm_anchors": np.allclose(rebuilt_lognorm_anchors, lognorm_source_anchors, rtol=CONFIG["rtol"], atol=CONFIG["atol"]),
        "embedding_shape": rebuilt_embedding.shape == (adata.n_obs, 2),
        "embedding_finite": np.isfinite(rebuilt_embedding).all(),
    }
    rebuilt_valid = all(rebuilt_checks.values())
    if not rebuilt_valid:
        failed_checks = [name for name, passed in rebuilt_checks.items() if not passed]
        lognorm_cell_delta = np.asarray(
            rebuilt_lognorm.transpose() @ gene_weights
            - lognorm_cell_sums
        ).ravel()
        raise RuntimeError(
            "Rebuilt r_export failed decoded content validation: "
            + ", ".join(failed_checks)
            + "; lognorm_cell_max_abs="
            + str(float(np.nanmax(np.abs(lognorm_cell_delta))))
        )

    os.replace(staging_dir / "expr.mtx", required_paths["counts"])
    os.replace(staging_dir / "expr_lognorm.mtx", required_paths["lognorm"])
    os.replace(staging_dir / "meta.csv", required_paths["meta"])
    os.replace(staging_dir / "genes.txt", required_paths["genes"])
    os.replace(staging_dir / "cells.txt", required_paths["cells"])
    os.replace(staging_dir / "umap.csv", required_paths["embedding"])
    os.replace(staging_dir / "info.json", required_paths["info"])
    shutil.rmtree(staging_dir)
    atexit.unregister(staging_cleanup)

source_stat = Path(source_path).stat()
artifact_state = {
    "counts": {"size": required_paths["counts"].stat().st_size, "mtime_ns": required_paths["counts"].stat().st_mtime_ns},
    "lognorm": {"size": required_paths["lognorm"].stat().st_size, "mtime_ns": required_paths["lognorm"].stat().st_mtime_ns},
    "meta": {"size": required_paths["meta"].stat().st_size, "mtime_ns": required_paths["meta"].stat().st_mtime_ns},
    "genes": {"size": required_paths["genes"].stat().st_size, "mtime_ns": required_paths["genes"].stat().st_mtime_ns},
    "cells": {"size": required_paths["cells"].stat().st_size, "mtime_ns": required_paths["cells"].stat().st_mtime_ns},
    "embedding": {"size": required_paths["embedding"].stat().st_size, "mtime_ns": required_paths["embedding"].stat().st_mtime_ns},
    "info": {"size": required_paths["info"].stat().st_size, "mtime_ns": required_paths["info"].stat().st_mtime_ns},
}
receipt_path.write_text(
    json.dumps(
        {
            "receipt_contract": "linealinkr064_r_export_decoded_validation_v1",
            "source": {
                "path": source_path,
                "size": source_stat.st_size,
                "mtime_ns": source_stat.st_mtime_ns,
            },
            "artifacts": artifact_state,
            "excluded_genes": excluded_genes,
            "input_scale": CONFIG["input_scale"],
        },
        ensure_ascii=False,
        indent=2,
    ) + "\\n",
    encoding="utf-8",
)

print("CACHED" if cache_valid else "REBUILT")
'
py <- gsub("__H5AD__", h5ad, py, fixed = TRUE)
py <- gsub("__OUTDIR__", outdir, py, fixed = TRUE)
py <- gsub("__SPECIES__", sp, py, fixed = TRUE)
py <- gsub("__INPUT_SCALE__", input_scale, py, fixed = TRUE)
py <- gsub("__DATASET__", ds_name, py, fixed = TRUE)
py <- gsub("__EXCLUSIONS__", file.path(CONFIG$dir$pipeline, "config/feature_exclusions.csv"), py, fixed = TRUE)

log_out <- file.path(CONFIG$dir$log, sprintf("export_%s.log", ds_name))
python_status <- sys::exec_wait(CONFIG$exe$python, args = c("-X", "utf8", "-c", py),
  std_out = log_out,
  std_err = file.path(CONFIG$dir$log, sprintf("export_%s.err", ds_name)))
if (python_status != 0L) {
  stop("r_export content validation failed for ", ds_name,
       "; see ", file.path(CONFIG$dir$log, sprintf("export_%s.err", ds_name)))
}
export_state <- tail(readLines(log_out, encoding = "UTF-8", warn = FALSE), 1L)
cat(sprintf("%s %s -> %s\n", export_state, ds_name, outdir))
