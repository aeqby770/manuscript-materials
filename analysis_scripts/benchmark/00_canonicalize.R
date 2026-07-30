#!/usr/bin/env Rscript
# =============================================================================
# 00_canonicalize.R
# Convert every registered dataset to the canonical schema (H5AD and RDS).
# Write canonical H5AD files to CONFIG$dir$canonical.
# =============================================================================

# path bootstrap (relative; set LLR_ASSETS to override)
.llr_cmd <- commandArgs(trailingOnly = FALSE)
.llr_file <- sub("^--file=", "", grep("^--file=", .llr_cmd, value = TRUE)[1])
if (!length(.llr_file) || is.na(.llr_file) || !nzchar(.llr_file)) {
  stop("Run via Rscript so --file= is set")
}
.llr_here <- dirname(normalizePath(.llr_file, winslash = "/", mustWork = TRUE))
source(file.path(.llr_here, "../config/llr_bootstrap.R"), local = FALSE)
.llr_assets <- llr_assets_dir()
llr_source_config()

suppressPackageStartupMessages({
  library(data.table)
})

# =============================================================================
# Canonical schema definition
# =============================================================================
# obs fields:
#   cell_id | cell_type | cluster | timepoint | pseudotime |
#   sample_id | condition | n_counts | n_genes | batch
# var fields:
#   gene_symbol
# layers：
#   counts（raw counts, genes x cells）
#   lognorm（log-normalized, genes x cells）
# =============================================================================

CANONICAL_FIELDS <- c("cell_id", "cell_type", "cluster", "timepoint",
                      "pseudotime", "sample_id", "condition",
                      "n_counts", "n_genes")

# Known aliases (obs columns -> canonical fields)
ALIAS_MAP <- list(
  cell_type  = c("cell_type", "celltype", "cao_cell_type", "meta_celltype",
                 "Cluster_Idents", "cell_type_main"),
  cluster    = c("cluster", "clusters", "louvain_cluster", "paul15_clusters",
                 "cao_cluster", "seurat_clusters"),
  timepoint  = c("time_point", "day", "age(days)"),
  pseudotime = c("pseudotime", "palantir_pseudotime"),
  sample_id  = c("sample_id", "sample", "orig.ident", "orig_ident",
                 "Sample_Name", "patient"),
  condition  = c("condition", "disease", "group", "diagnosis"),
  n_counts   = c("n_counts", "nCount_RNA", "total_counts"),
  n_genes    = c("n_genes", "nFeature_RNA", "n_genes_by_counts")
)

GuessCanonicalCol <- function(obs_colnames, canonical_field) {
  # Match known aliases.
  aliases <- ALIAS_MAP[[canonical_field]]
  for (a in aliases) {
    m <- match(a, tolower(obs_colnames))
    if (!is.na(m)) return(obs_colnames[m])
  }
  NA_character_
}

# =============================================================================
# Python helper: read an H5AD and write the canonical H5AD.
# =============================================================================
CanonicalizeH5ad <- function(ds_name, h5ad_path, out_dir) {
  # Use Python scanpy to read, canonicalize, and save
  py_script <- r"(
import scanpy as sc, pandas as pd, numpy as np, os, sys

inpath, outdir, dsname = sys.argv[1:4]

adata = sc.read(inpath)
if adata.n_obs == 0 or adata.n_vars == 0:
    sys.exit(1)

# Ensure counts in raw if available
if adata.raw is not None:
    raw_adata = adata.raw.to_adata()
    # Transfer raw counts
    adata.layers['counts'] = raw_adata.X.copy()
elif hasattr(adata, 'layers') and 'counts' in adata.layers:
    pass
else:
    # Assume X is counts
    adata.layers['counts'] = adata.X.copy()

# Log-normalized
import scanpy.preprocessing as scpp
scpp.normalize_total(adata, target_sum=1e4, key_added=None, inplace=False)
scpp.log1p(adata)
adata.layers['lognorm'] = adata.X.copy()

# Restore raw counts to X if raw existed
if adata.raw is not None:
    adata.X = raw_adata.X.copy()
else:
    adata.X = adata.layers['counts'].copy()

# Ensure obs columns
obs_map = {
    'cell_type': None, 'cluster': None, 'timepoint': None,
    'pseudotime': None, 'sample_id': None, 'condition': None
}

aliases = {
    'cell_type': ['cell_type','celltype','cao_cell_type','meta_celltype','Cluster_Idents','cell_type_main','ct'],
    'cluster': ['cluster','clusters','louvain_cluster','paul15_clusters','cao_cluster','seurat_clusters'],
     'timepoint': ['time_point','day','age(days)','Time_points','days','time_info','timepoint'],
    'pseudotime': ['pseudotime','palantir_pseudotime'],
    'sample_id': ['sample_id','sample','orig.ident','orig_ident','Sample_Name','patient'],
    'condition': ['condition','disease','group','diagnosis','status','Class','HF.etiology','Disease.status']
}

for canon, syns in aliases.items():
    obs_cols_lower = [c.lower() for c in adata.obs.columns]
    for s in syns:
        s_lower = s.lower()
        if s_lower in obs_cols_lower:
            idx = obs_cols_lower.index(s_lower)
            adata.obs[canon] = adata.obs.iloc[:, idx].astype(str)
            break
    if canon not in adata.obs.columns:
        adata.obs[canon] = ''

# Ensure key columns exist
for col in ['n_counts','n_genes']:
    if col not in adata.obs.columns:
        adata.obs[col] = 0

# Write canonical
if 'X_umap' not in adata.obsm:
    try:
        sc.pp.highly_variable_genes(adata, n_top_genes=2000)
        adata_sub = adata[:, adata.var.highly_variable].copy()
        sc.pp.pca(adata_sub, n_comps=30)
        sc.pp.neighbors(adata_sub, n_neighbors=30)
        sc.tl.umap(adata_sub)
        adata.obsm['X_umap'] = adata_sub.obsm['X_umap']
    except:
        adata.obsm['X_umap'] = np.zeros((adata.n_obs, 2))

adata.uns['linealinkr_canonical_contract'] = 'linealinkr064_canonical_v1'
adata.uns['linealinkr_canonical_dataset'] = dsname
adata.uns['linealinkr_canonical_source'] = os.path.normcase(os.path.realpath(inpath))
adata.write(os.path.join(outdir, f'{dsname}_canonical.h5ad'), compression='gzip')
print(f'OK: {dsname}')
)"

  status <- sys::exec_wait(
    CONFIG$exe$python,
    args = c("-X", "utf8", "-c", py_script, h5ad_path, out_dir, ds_name),
    std_out = file.path(CONFIG$dir$log, sprintf("canon_%s.log", ds_name)),
    std_err = file.path(CONFIG$dir$log, sprintf("canon_%s.err", ds_name))
  )

  out_path <- file.path(out_dir, sprintf("%s_canonical.h5ad", ds_name))
  if (identical(status, 0L) && file.exists(out_path)) out_path else NA_character_
}

ValidateCanonicalH5ad <- function(ds_name, source_path, canonical_path) {
  validator <- r"(
import anndata as ad
import numpy as np
import os
import sys
from scipy import sparse

dataset, source_path, canonical_path = sys.argv[1:4]
source = ad.read_h5ad(source_path, backed='r')
canonical = ad.read_h5ad(canonical_path, backed='r')

assert canonical.uns['linealinkr_canonical_contract'] == 'linealinkr064_canonical_v1'
assert canonical.uns['linealinkr_canonical_dataset'] == dataset
assert os.path.normcase(os.path.realpath(canonical.uns['linealinkr_canonical_source'])) == os.path.normcase(os.path.realpath(source_path))
assert source.n_obs == canonical.n_obs and source.n_vars == canonical.n_vars
assert np.array_equal(source.obs_names.to_numpy(), canonical.obs_names.to_numpy())
assert np.array_equal(source.var_names.to_numpy(), canonical.var_names.to_numpy())
assert set(['cell_type', 'cluster', 'timepoint', 'pseudotime', 'sample_id', 'condition', 'n_counts', 'n_genes']).issubset(canonical.obs.columns)
assert set(['counts', 'lognorm']).issubset(canonical.layers.keys())

row_index = np.unique(np.linspace(0, source.n_obs - 1, min(source.n_obs, 257), dtype=np.int64))
col_index = np.unique(np.linspace(0, source.n_vars - 1, min(source.n_vars, 257), dtype=np.int64))
source_matrix = source.raw.X if source.raw is not None else source.layers['counts'] if 'counts' in source.layers else source.X
source_anchor = source_matrix[row_index, :][:, col_index]
canonical_anchor = canonical.X[row_index, :][:, col_index]
counts_anchor = canonical.layers['counts'][row_index, :][:, col_index]
lognorm_anchor = canonical.layers['lognorm'][row_index, :][:, col_index]
source_anchor = source_anchor.toarray() if sparse.issparse(source_anchor) else np.asarray(source_anchor)
canonical_anchor = canonical_anchor.toarray() if sparse.issparse(canonical_anchor) else np.asarray(canonical_anchor)
counts_anchor = counts_anchor.toarray() if sparse.issparse(counts_anchor) else np.asarray(counts_anchor)
lognorm_anchor = lognorm_anchor.toarray() if sparse.issparse(lognorm_anchor) else np.asarray(lognorm_anchor)
assert np.array_equal(source_anchor, canonical_anchor, equal_nan=True)
assert np.array_equal(source_anchor, counts_anchor, equal_nan=True)
assert np.isfinite(lognorm_anchor).all()
print(f'TRUSTED\t{canonical.n_obs}\t{canonical.n_vars}')
)"
  validation <- sys::exec_internal(
    CONFIG$exe$python,
    args = c("-X", "utf8", "-c", validator, ds_name, source_path, canonical_path),
    error = FALSE
  )
  fields <- strsplit(trimws(rawToChar(validation$stdout)), "\t", fixed = TRUE)[[1L]]
  trusted <- identical(validation$status, 0L) && length(fields) == 3L && identical(fields[[1L]], "TRUSTED")
  list(
    trusted = trusted,
    n_obs = ifelse(trusted, as.integer(fields[[2L]]), NA_integer_),
    n_vars = ifelse(trusted, as.integer(fields[[3L]]), NA_integer_)
  )
}

# =============================================================================
# Main loop over the registered datasets.
# =============================================================================
cat("\n===== 00_canonicalize.R =====\n")
cat(sprintf("Datasets to canonicalize: %d\n\n", length(datasets)))

results <- data.table(
  dataset   = character(),
  status    = character(),
  source_path = character(),
  canonical_path = character(),
  n_obs = integer(),
  n_vars = integer()
)

for (nm in names(datasets)) {
  cat(sprintf("[%s] ... ", nm))

  out_path <- file.path(CONFIG$dir$canonical, sprintf("%s_canonical.h5ad", nm))
  source_path <- DatasetPath(nm)

  # Reuse only after decoding source and output and checking their content contract.
  cached <- if (file.exists(out_path)) {
    ValidateCanonicalH5ad(nm, source_path, out_path)
  } else {
    list(trusted = FALSE, n_obs = NA_integer_, n_vars = NA_integer_)
  }
  if (cached$trusted) {
    cat("CACHED (content contract)\n")
    results <- rbind(results, data.table(
      dataset = nm, status = "cached", source_path = source_path,
      canonical_path = normalizePath(out_path, winslash = "/", mustWork = TRUE),
      n_obs = cached$n_obs, n_vars = cached$n_vars
    ))
    next
  }

  ok <- CanonicalizeH5ad(nm, source_path, CONFIG$dir$canonical)
  verified <- if (!is.na(ok)) {
    ValidateCanonicalH5ad(nm, source_path, ok)
  } else {
    list(trusted = FALSE, n_obs = NA_integer_, n_vars = NA_integer_)
  }

  if (!verified$trusted) {
    cat("FAIL\n")
    results <- rbind(results, data.table(
      dataset = nm, status = "fail", source_path = source_path,
      canonical_path = NA_character_, n_obs = NA_integer_, n_vars = NA_integer_
    ))
  } else {
    cat(sprintf("OK → %s\n", basename(ok)))
    results <- rbind(results, data.table(
      dataset = nm, status = "ok", source_path = source_path,
      canonical_path = normalizePath(ok, winslash = "/", mustWork = TRUE),
      n_obs = verified$n_obs, n_vars = verified$n_vars
    ))
  }
}

# =============================================================================
# Summary.
# =============================================================================
cat("\n===== Canonicalization Summary =====\n")
print(results)
fwrite(results, file.path(CONFIG$dir$canonical, "canonical_manifest.csv"))
cat(sprintf("Done. %d/%d succeeded.\n",
    sum(results$status %in% c("ok", "cached")), nrow(results)))
