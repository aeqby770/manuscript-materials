#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))

# path bootstrap (relative; set LLR_ASSETS or LLR_PIPELINE to override)
.llr_cmd <- commandArgs(trailingOnly = FALSE)
.llr_file <- sub("^--file=", "", grep("^--file=", .llr_cmd, value = TRUE)[1])
if (!length(.llr_file) || is.na(.llr_file) || !nzchar(.llr_file)) {
  stop("Run this script via Rscript so --file= is set")
}
.llr_here <- dirname(normalizePath(.llr_file, winslash = "/", mustWork = TRUE))
source(file.path(.llr_here, "../config/llr_bootstrap.R"), local = FALSE)

llr_source_config()
llr_source_rel("config/panel_skip_rules.R")
Sys.setenv(
  LLR_VENDOR_ROOT = CONFIG$dir$algorithms,
  LLR_TICCI_ROOT = CONFIG$ticci$root,
  LLR_BENCHMARK_TRACK = "annotation_free",
  LLR_ROOT_CELL = ""
)

# Optional smoke mode uses an explicitly supplied isolated output tree.
smoke_root <- Sys.getenv("LLR_SMOKE_ROOT", unset = NA_character_)
if (!is.na(smoke_root) && nzchar(smoke_root)) {
  CONFIG$dir$algorithm_out <- file.path(smoke_root, "by_algorithm")
  CONFIG$dir$algorithm_misc <- file.path(smoke_root, "misc")
  CONFIG$dir$log <- file.path(CONFIG$dir$assets, "09_logs", "smoke_competitors")
  dir.create(CONFIG$dir$algorithm_out, recursive = TRUE, showWarnings = FALSE)
  dir.create(CONFIG$dir$algorithm_misc, recursive = TRUE, showWarnings = FALSE)
  dir.create(CONFIG$dir$log, recursive = TRUE, showWarnings = FALSE)
}

args <- commandArgs(trailingOnly = TRUE)
filter_datasets <- NULL; filter_algos <- NULL; skip_datasets <- character(0)
skip_by_algo <- list(); force_rerun <- FALSE
reuse_existing_results <- FALSE
register_validated <- FALSE
split_ds <- function(x) {
  x <- trimws(x)
  x <- x[nzchar(x)]
  unique(x)
}
read_ds_file <- function(f) {
  lines <- readLines(f, warn = FALSE)
  lines <- trimws(lines)
  lines <- lines[!grepl("^#", lines) & nzchar(lines)]
  if (length(lines) == 1L && grepl(",", lines[[1L]], fixed = TRUE)) {
    return(split_ds(unlist(strsplit(lines[[1L]], ",", fixed = TRUE))))
  }
  split_ds(lines)
}
read_skip_by_algo_tsv <- function(f) {
  if (!file.exists(f)) return(list())
  lines <- readLines(f, warn = FALSE)
  lines <- trimws(lines)
  lines <- lines[!grepl("^#", lines) & nzchar(lines)]
  if (!length(lines)) return(list())
  out <- list()
  for (ln in lines[-1L]) {
    parts <- strsplit(ln, "\t", fixed = TRUE)[[1L]]
    if (length(parts) < 2L) next
    algo <- parts[[1L]]
    ds <- split_ds(unlist(strsplit(parts[[2L]], ",", fixed = TRUE)))
    if (length(ds)) out[[algo]] <- unique(c(out[[algo]], ds))
  }
  out
}
PanelSkip <- function(algo, ds) {
  panel_should_skip(algo, ds, skip_datasets, skip_by_algo)$skip
}
PanelSkipReason <- function(algo, ds) {
  panel_should_skip(algo, ds, skip_datasets, skip_by_algo)$reason
}
for (a in args) {
  if (grepl("^--datasets-file=", a)) {
    f <- sub("^--datasets-file=", "", a)
    if (!file.exists(f)) stop("datasets file not found: ", f)
    filter_datasets <- read_ds_file(f)
  }
  if (grepl("^--datasets=", a)) {
    raw <- sub("^--datasets=", "", a)
    sep <- if (grepl(";", raw, fixed = TRUE)) ";" else ","
    filter_datasets <- split_ds(unlist(strsplit(raw, sep, fixed = TRUE)))
  }
  if (grepl("^--skip-datasets-file=", a)) {
    f <- sub("^--skip-datasets-file=", "", a)
    if (!file.exists(f)) stop("skip datasets file not found: ", f)
    skip_datasets <- unique(c(skip_datasets, read_ds_file(f)))
  }
  if (grepl("^--skip-datasets=", a)) {
    raw <- sub("^--skip-datasets=", "", a)
    sep <- if (grepl(";", raw, fixed = TRUE)) ";" else ","
    skip_datasets <- unique(c(skip_datasets, split_ds(unlist(strsplit(raw, sep, fixed = TRUE)))))
  }
  if (grepl("^--skip-by-algo-file=", a)) {
    f <- sub("^--skip-by-algo-file=", "", a)
    if (!file.exists(f)) stop("skip-by-algo file not found: ", f)
    skip_by_algo <- read_skip_by_algo_tsv(f)
  }
  if (grepl("^--algos=", a)) filter_algos <- split_ds(strsplit(sub("^--algos=", "", a), ",", fixed = TRUE)[[1]])
  if (a == "--force") force_rerun <- TRUE
  if (a %in% c("--reuse-existing", "--reuse-existing-results")) {
    reuse_existing_results <- TRUE
  }
  if (a == "--register-validated") register_validated <- TRUE
}

# CLI exclusions are explicit. Programmatic applicability rules are evaluated
# by panel_skip_rule_reason().
if (length(skip_datasets)) {
  cat("Panel skip (global):", paste(skip_datasets, collapse = ", "), "\n")
}
if (length(skip_by_algo)) {
  cat("Panel skip (by algo):", length(skip_by_algo), "algorithms (+ programmatic rules)\n")
}

# Official outputs use the result root consumed by evaluation. Smoke runs use
# an isolated output tree.
formal_competitor_root <- normalizePath(
  CONFIG$dir$competitor_baseline,
  winslash = "/",
  mustWork = TRUE
)
OUT <- if (!is.na(smoke_root) && nzchar(smoke_root)) {
  CONFIG$dir$algorithm_out
} else {
  file.path(formal_competitor_root, "by_algorithm")
}
OUT_MISC <- if (!is.na(smoke_root) && nzchar(smoke_root)) {
  CONFIG$dir$algorithm_misc
} else {
  file.path(formal_competitor_root, "misc")
}
ADAPTERS <- file.path(CONFIG$dir$pipeline, "adapters")
# The Job Object launcher applies the memory guard and preserves child exit status.
MEM_GUARD_PY <- file.path(ADAPTERS, "_mem_guard.py")
MEM_BUDGET_GB <- 55
dir.create(OUT_MISC, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ---- Output checks --------------------------------------------------------
# Reuse requires a matching run-registry record and a decoded primary artifact.
DiskOK <- function(path, min_size = 100) {
  file.exists(path) && file.info(path)$size > min_size
}

RUN_REGISTRY_PATH <- file.path(
  CONFIG$dir$paper_tables,
  "_registry",
  sprintf("benchmark_results_%s.csv", CONFIG$llr$benchmark_protocol)
)
RUN_REGISTRY_REQUIRED_COLUMNS <- c(
  "algorithm", "dataset", "status", "benchmark_protocol",
  "result_revision"
)
RUN_REGISTRY <- if (file.exists(RUN_REGISTRY_PATH)) {
  fread(RUN_REGISTRY_PATH, showProgress = FALSE)
} else {
  data.table()
}

ReuseRegistered <- function(algo, ds_name) {
  if (!reuse_existing_results || force_rerun) return(FALSE)
  if (!all(RUN_REGISTRY_REQUIRED_COLUMNS %in% names(RUN_REGISTRY))) return(FALSE)
  nrow(RUN_REGISTRY[
    algorithm == algo &
      dataset == ds_name &
      status == "OK" &
      benchmark_protocol == CONFIG$llr$benchmark_protocol &
      result_revision == CONFIG$llr$result_revision
  ]) == 1L
}

REFERENCE_IDS_PY <- paste(c(
  "import h5py",
  "import pandas as pd",
  "import sys",
  "h = h5py.File(sys.argv[1], 'r')",
  "key = h['obs'].attrs['_index']",
  "key = key.decode('utf-8') if isinstance(key, bytes) else str(key)",
  "cell_ids = h['obs'][key].asstr()[...]",
  "h.close()",
  "pd.DataFrame({'cell_id': cell_ids}).to_csv(sys.stdout, index=False, lineterminator='\\n')"
), collapse = "\n")
REFERENCE_ID_CACHE <- new.env(parent = emptyenv())

ReadCanonicalCellIds <- function(ds_name) {
  if (exists(ds_name, envir = REFERENCE_ID_CACHE, inherits = FALSE)) {
    return(get(ds_name, envir = REFERENCE_ID_CACHE, inherits = FALSE))
  }
  decoded <- sys::exec_internal(
    CONFIG$exe$python,
    args = c("-X", "utf8", "-c", REFERENCE_IDS_PY, DatasetPath(ds_name)),
    error = FALSE
  )
  stopifnot(identical(decoded$status, 0L))
  reference <- fread(
    text = rawToChar(decoded$stdout),
    sep = ",",
    showProgress = FALSE
  )
  stopifnot(
    identical(names(reference), "cell_id"),
    nrow(reference) > 0L,
    !anyNA(reference$cell_id),
    all(nzchar(reference$cell_id)),
    !anyDuplicated(reference$cell_id)
  )
  cell_ids <- as.character(reference$cell_id)
  assign(ds_name, cell_ids, envir = REFERENCE_ID_CACHE)
  cell_ids
}

PseudotimeCsvOK <- function(path, ds_name) {
  if (!DiskOK(path)) return(FALSE)
  x <- fread(path, colClasses = c(cell_id = "character"))
  if (!identical(names(x), c("cell_id", "pseudotime"))) return(FALSE)
  reference_ids <- ReadCanonicalCellIds(ds_name)
  output_ids <- as.character(x$cell_id)
  pt <- as.numeric(x$pseudotime)
  length(output_ids) == length(reference_ids) &&
    !anyNA(output_ids) && all(nzchar(output_ids)) &&
    !anyDuplicated(output_ids) && setequal(output_ids, reference_ids) &&
    all(is.finite(pt)) && diff(range(pt)) > 0
}

R_PSEUDOTIME_OUTPUT <- c(
  slingshot = "slingshot_pseudotime.csv",
  SCORPIUS = "scorpius_pseudotime.csv",
  monocle3 = "monocle3_pseudotime.csv",
  TSCAN = "tscan_pseudotime.csv"
)

# Reuse is opt-in. A directory is reusable only when the adapter wrote its
# completion marker and the primary artifact still covers the requested cells.
# `--force` always takes precedence over reuse.
DirCached <- function(algo, ds_name) {
  full <- file.path(OUT, algo, ds_name)
  if (!dir.exists(full)) return(FALSE)
  marker_ok <- file.exists(file.path(full, "DONE"))
  primary_name <- unname(R_PSEUDOTIME_OUTPUT[algo])
  primary_name <- primary_name[!is.na(primary_name)]
  marker_ok && (
    !length(primary_name) ||
      PseudotimeCsvOK(file.path(full, primary_name), ds_name)
  )
}

MarkDirDone <- function(algo, ds_name, out_dir) {
  writeLines(as.character(Sys.time()), file.path(out_dir, "DONE"))
}

# Remove a pre-existing flat output before the child run.
MarkFlatFail <- function(path) {
  if (file.exists(path)) unlink(path)
}

RunPython <- function(script, log_name) {
  sys::exec_wait(CONFIG$exe$python, args = c(script),
    std_out = file.path(CONFIG$dir$log, sprintf("%s.log", log_name)),
    std_err = file.path(CONFIG$dir$log, sprintf("%s.err", log_name)))
}

# cellrank2 runs inside a WSL conda env (cr2) that ships slepc4py/petsc4py, so
# GPCCA uses sparse Krylov-Schur instead of densifying the transition matrix
# (Windows has no slepc wheels). Same algorithm/parameters — only the Schur
# backend changes, so results are equivalent but the O(n^3) dense step is gone.
WSL_DISTRO <- "Ubuntu-22.04"
WSL_CR2_PY <- "/root/anaconda3/envs/cr2/bin/python"
ToWSL <- function(p) {
  p <- gsub("\\\\", "/", p)
  if (grepl("^[A-Za-z]:/", p)) {
    p <- paste0("/mnt/", tolower(substr(p, 1, 1)), substring(p, 3))
  }
  p
}
WriteInlinePy <- function(py, tag) {
  d <- file.path(CONFIG$dir$cache, "inline_py")
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  tag <- gsub("[^A-Za-z0-9_.-]", "_", tag)
  tf <- file.path(d, paste0(tag, ".py"))
  writeLines(py, tf, useBytes = TRUE)
  tf
}
RunPythonWSL <- function(script, log_name) {
  script_wsl <- ToWSL(normalizePath(script, winslash = "/", mustWork = FALSE))
  system2("wsl", args = c("-d", WSL_DISTRO, "--", WSL_CR2_PY, script_wsl),
    stdout = file.path(CONFIG$dir$log, sprintf("%s.log", log_name)),
    stderr = file.path(CONFIG$dir$log, sprintf("%s.err", log_name)))
}

RunRScript <- function(algo, ds_name, script_name) {
  cat(sprintf("  [%s] %s ... ", algo, ds_name))
  out_dir <- file.path(OUT, algo, ds_name); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  force_this <- force_rerun
  if (ReuseRegistered(algo, ds_name) && DirCached(algo, ds_name)) {
    out_marker <- file.path(out_dir, "DONE")
    if (!file.exists(out_marker)) writeLines(as.character(Sys.time()), out_marker)
    cat("CACHED\n"); return(out_dir)
  }
  old_entries <- list.files(out_dir, all.files = TRUE, no.. = TRUE, full.names = TRUE)
  if (length(old_entries)) unlink(old_entries, recursive = TRUE, force = TRUE)
  script <- file.path(ADAPTERS, script_name)
  if (!file.exists(script)) { cat(sprintf("FAIL (no adapter: %s)\n", script_name)); return(NA) }
  status <- system2(CONFIG$exe$python,
    args = c("-X", "utf8", MEM_GUARD_PY, "--limit-gb", MEM_BUDGET_GB, "--",
             CONFIG$exe$rscript, script, ds_name, out_dir),
    stdout = file.path(CONFIG$dir$log, sprintf("%s_%s.log", algo, ds_name)),
    stderr = file.path(CONFIG$dir$log, sprintf("%s_%s.err", algo, ds_name)))
  out_marker <- file.path(out_dir, "DONE")
  produced <- setdiff(list.files(out_dir, all.files = TRUE, no.. = TRUE, full.names = TRUE), out_marker)
  produced <- produced[file.exists(produced) & file.info(produced)$size > 100]
  primary_name <- unname(R_PSEUDOTIME_OUTPUT[algo])
  primary_name <- primary_name[!is.na(primary_name)]
  produced_ok <- if (length(primary_name)) {
    PseudotimeCsvOK(file.path(out_dir, primary_name), ds_name)
  } else {
    length(produced) > 0L
  }
  if (status == 0 && produced_ok) {
    MarkDirDone(algo, ds_name, out_dir)
    cat("OK\n"); out_dir
  } else if (status == 0) {
    if (file.exists(out_marker)) unlink(out_marker)
    cat("FAIL (invalid primary output)\n"); NA
  } else {
    # Remove the completion marker after a non-zero exit.
    if (file.exists(out_marker)) unlink(out_marker)
    cat("FAIL\n"); NA
  }
}

# ===== Python inline runners =====

RunPalantir <- function(ds_name) {
  cat(sprintf("  [Palantir] %s ... ", ds_name))
  out_file <- FlatOut(sprintf("palantir_%s_result.h5ad", ds_name))
  if (ReuseRegistered("Palantir", ds_name) && DiskOK(out_file)) { cat("CACHED\n"); return(out_file) }
  MarkFlatFail(out_file)
  adapter <- file.path(ADAPTERS, "run_palantir.py")
  status <- system2(CONFIG$exe$python,
    args = c("-X", "utf8", adapter, DatasetPath(ds_name), out_file),
    stdout = file.path(CONFIG$dir$log, sprintf("palantir_%s.log", ds_name)),
    stderr = file.path(CONFIG$dir$log, sprintf("palantir_%s.err", ds_name)))
  ok <- status == 0 && file.exists(out_file) && file.info(out_file)$size > 100
  if (ok) { cat("OK\n"); out_file } else { MarkFlatFail(out_file); cat("FAIL\n"); NA }
}

RunLIANA <- function(ds_name) {
  cat(sprintf("  [LIANA] %s ... ", ds_name))
  out_file <- FlatOut(sprintf("liana_%s_result.csv", ds_name))
  cells_file <- FlatOut(sprintf("liana_%s_cells.csv", ds_name))
  run_file <- FlatOut(sprintf("liana_%s_run.json", ds_name))
  liana_output_ok <- DiskOK(out_file, min_size = 80) &&
    DiskOK(cells_file, min_size = 20) && DiskOK(run_file, min_size = 80)
  if (ReuseRegistered("LIANA", ds_name) && liana_output_ok) { cat("CACHED\n"); return(out_file) }
  MarkFlatFail(out_file)
  MarkFlatFail(cells_file)
  MarkFlatFail(run_file)
  sp <- CONFIG$dataset_species[[ds_name]]; if (is.null(sp)) sp <- "human"
  resource <- if (sp == "mouse") "mouseconsensus" else "consensus"
  adapter <- file.path(ADAPTERS, "run_liana.py")
  status <- system2(CONFIG$exe$python,
    args = c("-X", "utf8", adapter, DatasetPath(ds_name), out_file, resource),
    stdout = file.path(CONFIG$dir$log, sprintf("liana_%s.log", ds_name)),
    stderr = file.path(CONFIG$dir$log, sprintf("liana_%s.err", ds_name)))
  liana_output_ok <- DiskOK(out_file, min_size = 80) &&
    DiskOK(cells_file, min_size = 20) && DiskOK(run_file, min_size = 80)
  if (status == 0 && liana_output_ok) {
    cat("OK\n"); out_file
  } else {
    MarkFlatFail(out_file)
    MarkFlatFail(cells_file)
    MarkFlatFail(run_file)
    cat("FAIL\n"); NA
  }
}

RunCellRank2 <- function(ds_name) {
  cat(sprintf("  [cellrank2] %s ... ", ds_name))
  out_file <- FlatOut(sprintf("cellrank2_%s_result.h5ad", ds_name))
  if (ReuseRegistered("cellrank2", ds_name) && DiskOK(out_file)) { cat("CACHED\n"); return(out_file) }
  py <- sprintf('import scanpy as sc, cellrank as cr, numpy as np, os, sys, scipy.sparse
os.environ.setdefault("UCX_TLS", "tcp,self")
a = sc.read(r"%s")
# No subsampling: run on the full cell set. The cr2 WSL env ships slepc4py/petsc4py,
# so GPCCA uses the sparse Krylov-Schur backend (no dense transition matrix) and
# fate probabilities are solved with the sparse PETSc GMRES iterative solver below.
sc.pp.filter_genes(a, min_cells=1)
# Robust raw-vs-processed input handling (CellRank expects log-normalized X).
# Prefer explicit layers; otherwise infer from X by integrality: raw counts are
# non-negative integers -> normalize_total + log1p; non-integer X is taken as
# already log-normalized. This replaces the fragile "a.X.max() < 50" heuristic
# that mislabelled shallow raw-count datasets as log-normalized.
if "lognorm" in a.layers:
  a.X = a.layers["lognorm"].copy()
elif "counts" in a.layers:
  a.X = a.layers["counts"].copy(); sc.pp.normalize_total(a); sc.pp.log1p(a)
else:
  _xd = a.X.data if scipy.sparse.issparse(a.X) else np.asarray(a.X).ravel()
  if _xd.size and np.all(_xd >= 0) and np.allclose(_xd, np.round(_xd)):
    sc.pp.normalize_total(a); sc.pp.log1p(a)
# Manual HVG by top variance: avoids seurat expm1 overflow (lognorm max ~3223),
# skmisc dependency (not in WSL cr2 env), and cell_ranger duplicate bin edges.
_Xh = a.X
if scipy.sparse.issparse(_Xh):
  _mh = np.asarray(_Xh.mean(axis=0)).ravel()
  _vh = np.asarray(_Xh.power(2).mean(axis=0)).ravel() - _mh ** 2
else:
  _mh = np.asarray(_Xh).mean(axis=0)
  _vh = np.square(np.asarray(_Xh)).mean(axis=0) - _mh ** 2
_vh = np.nan_to_num(_vh, nan=0.0, posinf=0.0, neginf=0.0)
_nh = min(2000, a.n_vars - 1)
_hvg = np.argsort(_vh)[-_nh:]
a.var["highly_variable"] = np.zeros(a.n_vars, dtype=bool)
a.var["highly_variable"].values[_hvg] = True
a = a[:, a.var["highly_variable"]].copy()
a_hvg = a.copy()
_CR2_SEEDS = [%d]
_last_err = None
g = None
for _seed in _CR2_SEEDS:
  a = a_hvg.copy()
  sc.pp.pca(a, random_state=_seed)
  sc.pp.neighbors(a, random_state=_seed)
  ck = cr.kernels.ConnectivityKernel(a)
  ck.compute_transition_matrix()
  from scipy.sparse import csgraph as _csgraph
  _T = ck.transition_matrix
  _n_comp, _ = _csgraph.connected_components(_T, directed=False)
  g = cr.estimators.GPCCA(ck)
  g.compute_schur(n_components=20)
  _eig = getattr(g, "eigendecomposition", None)
  _eg = _eig.get("eigengap") if isinstance(_eig, dict) else None
  n_states = min(12, max(2, int(_eg))) if _eg is not None else 5
  n_states = max(n_states, min(12, _n_comp))
  g.compute_macrostates(n_states=n_states)
  g.predict_terminal_states()
  try:
    if a.n_obs <= 3000:
      g.compute_fate_probabilities(solver="direct", use_petsc=False, show_progress_bar=False)
    else:
      try:
        g.compute_fate_probabilities(solver="gmres", use_petsc=True, show_progress_bar=False)
      except ValueError as _fp_err:
        if "do not sum to 1" in str(_fp_err) and a.n_obs <= 8000:
          g.compute_fate_probabilities(solver="direct", use_petsc=False, show_progress_bar=False)
        else:
          raise
    break
  except ValueError as _fp_err:
    if "do not sum to 1" not in str(_fp_err):
      raise
    _last_err = _fp_err
    g = None
else:
  raise _last_err
a.obs["cr2_term"] = g.terminal_states.astype(str)
a.write(r"%s")
print(f"OK: {g.fate_probabilities.shape}")
', ToWSL(DatasetPath(ds_name)), CONFIG$seed, ToWSL(out_file))
  tf <- WriteInlinePy(py, sprintf("cr2_%s", ds_name))
  status <- RunPythonWSL(tf, sprintf("cr2_%s", ds_name))
  unlink(tf)
  ok <- status == 0 && DiskOK(out_file)
  if (ok) { cat("OK\n"); out_file } else { MarkFlatFail(out_file); cat("FAIL\n"); NA }
}

FlatOut <- function(name) file.path(OUT_MISC, name)

RunDPT <- function(ds_name) {
  cat(sprintf("  [DPT] %s ... ", ds_name))
  out_file <- FlatOut(sprintf("dpt_%s_pseudotime.csv", ds_name))
  if (ReuseRegistered("DPT", ds_name) && PseudotimeCsvOK(out_file, ds_name)) { cat("CACHED\n"); return(out_file) }
  MarkFlatFail(out_file)
  adapter <- file.path(ADAPTERS, "run_dpt.py")
  status <- system2(CONFIG$exe$python,
    args = c("-X", "utf8", adapter, DatasetPath(ds_name), out_file, "auto"),
    stdout = file.path(CONFIG$dir$log, sprintf("dpt_%s.log", ds_name)),
    stderr = file.path(CONFIG$dir$log, sprintf("dpt_%s.err", ds_name)))
  if (status == 0 && PseudotimeCsvOK(out_file, ds_name)) { cat("OK\n"); out_file } else { MarkFlatFail(out_file); cat("FAIL\n"); NA }
}

RunPAGA <- function(ds_name) {
  cat(sprintf("  [PAGA] %s ... ", ds_name))
  out_file <- FlatOut(sprintf("paga_%s_connectivities.csv", ds_name))
  if (ReuseRegistered("PAGA", ds_name) && DiskOK(out_file)) { cat("CACHED\n"); return(out_file) }
  adapter <- file.path(ADAPTERS, "run_paga.py")
  status <- system2(CONFIG$exe$python,
    args = c(adapter, DatasetPath(ds_name), out_file, "metadata"),
    stdout = file.path(CONFIG$dir$log, sprintf("paga_%s.log", ds_name)),
    stderr = file.path(CONFIG$dir$log, sprintf("paga_%s.err", ds_name)))
  if (status == 0 && DiskOK(out_file)) { cat("OK\n"); out_file } else { MarkFlatFail(out_file); cat("FAIL\n"); NA }
}

# STREAM-like (kNN Dijkstra pseudotime)
RunSTREAM <- function(ds_name) {
  cat(sprintf("  [STREAM] %s ... ", ds_name))
  out_file <- file.path(OUT, sprintf("stream_%s_pseudotime.csv", ds_name))
  if (ReuseRegistered("STREAM", ds_name) && DiskOK(out_file)) { cat("CACHED\n"); return(out_file) }
  py <- sprintf('import scanpy as sc, pandas as pd, numpy as np, os
a = sc.read(r"%s")
if "lognorm" in a.layers: a.X = a.layers["lognorm"].copy()
sc.pp.highly_variable_genes(a, n_top_genes=2000); a = a[:, a.var.highly_variable].copy()
sc.pp.pca(a); sc.pp.neighbors(a); sc.tl.umap(a)
from sklearn.neighbors import NearestNeighbors
umap = a.obsm["X_umap"]; nn = NearestNeighbors(n_neighbors=15).fit(umap)
dist, idx = nn.kneighbors(umap); root = np.argmin(umap[:, 0])
import heapq
n = a.n_obs; d = np.full(n, np.inf); d[root] = 0; visited = np.zeros(n, dtype=bool)
pq = [(0, root)]
while pq:
    du, u = heapq.heappop(pq)
    if visited[u]: continue
    visited[u] = True
    for vi, w in enumerate(idx[u]):
        v = idx[u][vi]; wd = dist[u][vi]
        if wd > 0 and d[u] + wd < d[v]:
            d[v] = d[u] + wd; heapq.heappush(pq, (d[v], v))
pd.DataFrame({"pseudotime": d}, index=a.obs_names).to_csv(r"%s")
print("OK")
', DatasetPath(ds_name), out_file)
  tf <- WriteInlinePy(py, sprintf("stream_%s", ds_name))
  status <- RunPython(tf, sprintf("stream_%s", ds_name))
  unlink(tf)
  ok <- status == 0 && DiskOK(out_file)
  if (ok) { cat("OK\n"); out_file } else { MarkFlatFail(out_file); cat("FAIL\n"); NA }
}

# Wishbone-like (DPT + branch score)
RunWishbone <- function(ds_name) {
  cat(sprintf("  [Wishbone] %s ... ", ds_name))
  out_file <- file.path(OUT, sprintf("wishbone_%s_pseudotime.csv", ds_name))
  if (ReuseRegistered("Wishbone", ds_name) && DiskOK(out_file)) { cat("CACHED\n"); return(out_file) }
  py <- 'import scanpy as sc, pandas as pd, numpy as np, os
inpath = r"PATH_PLACEHOLDER"
outpath = r"OUT_PLACEHOLDER"
a = sc.read(inpath)
if "lognorm" in a.layers: a.X = a.layers["lognorm"].copy()
sc.pp.highly_variable_genes(a, n_top_genes=2000)
a = a[:, a.var.highly_variable].copy()
sc.pp.pca(a); sc.pp.neighbors(a); sc.tl.diffmap(a)
# Cluster-based root selection (consistent with DPT/StaVia, no data leakage)
ROOT_CANDIDATES = ("0", "root", "Root", "STEM", "stem", "branch_0", "early", "Early")
cluster_col = None
for cn in ("clusters", "cluster", "cell_type", "celltype", "seurat_clusters"):
    if cn not in a.obs.columns: continue
    v = a.obs[cn].astype(str); valid = v.notna() & (v != "")
    if v[valid].nunique() >= 2: cluster_col = cn; break
if cluster_col is None:
    sc.tl.leiden(a, resolution=1.0, key_added="benchmark_clusters")
    cluster_col = "benchmark_clusters"
labels = a.obs[cluster_col].astype(str).values
root_label = None
for c in ROOT_CANDIDATES:
    if (labels == c).any(): root_label = c; break
if root_label is None:
    num = pd.to_numeric(pd.Series(labels), errors="coerce")
    if np.isfinite(num).any(): root_label = str(np.nanmin(num.values))
if root_label is not None:
    mask = labels == root_label; root_cells = np.where(mask)[0]
    if "X_pca" in a.obsm and len(root_cells) > 0:
        pca = a.obsm["X_pca"][root_cells, :]; centroid = pca.mean(axis=0)
        dists = np.linalg.norm(pca - centroid, axis=1)
        a.uns["iroot"] = int(root_cells[np.argmin(dists)])
    else: a.uns["iroot"] = int(root_cells[0]) if len(root_cells) > 0 else 0
else: a.uns["iroot"] = 0
sc.tl.dpt(a)
pt = a.obs["dpt_pseudotime"].values
from sklearn.neighbors import NearestNeighbors
x = a.obsm.get("X_umap", a.obsm.get("X_pca", np.zeros((a.n_obs,2))))
nn = NearestNeighbors(n_neighbors=10).fit(x)
_, idx = nn.kneighbors(x)
bs = np.array([np.std(pt[i]) for i in idx])
pd.DataFrame({"pseudotime": pt, "branch_score": bs}, index=a.obs_names).to_csv(outpath)
print("OK")'
  py <- gsub("PATH_PLACEHOLDER", DatasetPath(ds_name), py, fixed=TRUE)
  py <- gsub("OUT_PLACEHOLDER", out_file, py, fixed=TRUE)
  tf <- WriteInlinePy(py, sprintf("wishbone_%s", ds_name))
  status <- RunPython(tf, sprintf("wishbone_%s", ds_name))
  unlink(tf)
  ok <- status == 0 && DiskOK(out_file)
  if (ok) { cat("OK\n"); out_file } else { MarkFlatFail(out_file); cat("FAIL\n"); NA }
}

# MultistageOT (OT-based pseudotime + fate probabilities)
RunMultistageOT <- function(ds_name) {
  cat(sprintf("  [MultistageOT] %s ... ", ds_name))
  out_file <- FlatOut(sprintf("multistageot_%s_pseudotime.csv", ds_name))
  out_dir  <- file.path(OUT, "MultistageOT", ds_name)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  if (ReuseRegistered("MultistageOT", ds_name) && PseudotimeCsvOK(out_file, ds_name)) { cat("CACHED\n"); return(out_file) }
  MarkFlatFail(out_file)
  adapter <- file.path(ADAPTERS, "run_multistageot.py")
  status <- system2(CONFIG$exe$python,
    args = c("-X", "utf8", adapter, DatasetPath(ds_name), out_file),
    stdout = file.path(CONFIG$dir$log, sprintf("multistageot_%s.log", ds_name)),
    stderr = file.path(CONFIG$dir$log, sprintf("multistageot_%s.err", ds_name)))
  if (status == 0 && PseudotimeCsvOK(out_file, ds_name)) { cat("OK\n"); out_file } else { MarkFlatFail(out_file); cat("FAIL\n"); NA }
}

# scFates (principal tree pseudotime via ElPiGraph/PPT)
RunScFates <- function(ds_name) {
  cat(sprintf("  [scFates] %s ... ", ds_name))
  out_file <- FlatOut(sprintf("scfates_%s_pseudotime.csv", ds_name))
  if (ReuseRegistered("scFates", ds_name) && PseudotimeCsvOK(out_file, ds_name)) { cat("CACHED\n"); return(out_file) }
  MarkFlatFail(out_file)
  stale_sidecars <- file.path(OUT_MISC, sprintf(c(
    "scfates_lineage_edges_%s.csv",
    "scfates_lineage_segments_%s.csv",
    "scfates_lineage_cell_segments_%s.csv"
  ), ds_name))
  unlink(stale_sidecars[file.exists(stale_sidecars)])
  adapter <- file.path(ADAPTERS, "run_scfates.py")
  status <- system2(CONFIG$exe$python,
    args = c("-X", "utf8", adapter, DatasetPath(ds_name), out_file),
    stdout = file.path(CONFIG$dir$log, sprintf("scfates_%s.log", ds_name)),
    stderr = file.path(CONFIG$dir$log, sprintf("scfates_%s.err", ds_name)))
  if (status == 0 && PseudotimeCsvOK(out_file, ds_name)) { cat("OK\n"); out_file } else { MarkFlatFail(out_file); cat("FAIL\n"); NA }
}

RunPHLOWER <- function(ds_name) {
  cat(sprintf("  [PHLOWER] %s ... ", ds_name))
  out_file <- file.path(OUT, sprintf("phlower_%s_pseudotime.csv", ds_name))
  if (ReuseRegistered("PHLOWER", ds_name) && DiskOK(out_file)) { cat("CACHED\n"); return(out_file) }
  py <- "E:/conda/envs/phlower-official/python.exe"
  adapter <- file.path(ADAPTERS, "run_phlower.py")
  if (!file.exists(py)) { cat("SKIP (no official PHLOWER python)\n"); return(NA) }
  status <- system2(py,
    args = c(adapter, DatasetPath(ds_name), out_file),
    stdout = file.path(CONFIG$dir$log, sprintf("phlower_%s.log", ds_name)),
    stderr = file.path(CONFIG$dir$log, sprintf("phlower_%s.err", ds_name)))
  if (status == 0 && DiskOK(out_file)) { cat("OK\n"); out_file } else { MarkFlatFail(out_file); cat("FAIL\n"); NA }
}

# TICCI, CellPhoneDB, wot — keeping existing implementations short for brevity
RunTICCI <- function(ds_name) {
  cat(sprintf("  [TICCI] %s ... ", ds_name))
  out_dir <- file.path(OUT, "TICCI", ds_name)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  marker <- file.path(out_dir, "ticci_summary.csv")
  if (ReuseRegistered("TICCI", ds_name) && DiskOK(marker, min_size = 50)) { cat("CACHED\n"); return(out_dir) }
  adapter <- file.path(ADAPTERS, "run_ticci_benchmark.py")
  status <- system2(CONFIG$exe$python,
    args = c(adapter, DatasetPath(ds_name), out_dir),
    stdout = file.path(CONFIG$dir$log, sprintf("ticci_%s.log", ds_name)),
    stderr = file.path(CONFIG$dir$log, sprintf("ticci_err_%s.log", ds_name)))
  if (status == 0 && DiskOK(marker, min_size = 50)) {
    MarkDirDone("TICCI", ds_name, out_dir)
    cat("OK\n"); out_dir
  } else { MarkFlatFail(marker); cat("FAIL\n"); NA }
}
RunCellPhoneDB <- function(ds_name) {
  cat(sprintf("  [CellPhoneDB] %s ... ", ds_name))
  out_dir <- file.path(OUT, "CellPhoneDB", ds_name)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  marker <- file.path(out_dir, "cpdb_done.txt")
  if (ReuseRegistered("CellPhoneDB", ds_name) && DirCached("CellPhoneDB", ds_name)) { cat("CACHED\n"); return(out_dir) }
  adapter <- file.path(ADAPTERS, "run_cellphonedb.py")
  status <- system2(CONFIG$exe$python,
    args = c(adapter, DatasetPath(ds_name), out_dir),
    stdout = file.path(CONFIG$dir$log, sprintf("cellphonedb_%s.log", ds_name)),
    stderr = file.path(CONFIG$dir$log, sprintf("cellphonedb_%s.err", ds_name)))
  if (status == 0 && file.exists(marker)) { MarkDirDone("CellPhoneDB", ds_name, out_dir); cat("OK\n"); out_dir } else { cat("FAIL\n"); NA }
}
RunNichenetr <- function(ds_name) {
  cat(sprintf("  [nichenetr] %s ... ", ds_name))
  out_file <- FlatOut(sprintf("nichenetr_%s_ligand_activity.csv", ds_name))
  if (ReuseRegistered("nichenetr", ds_name) && DiskOK(out_file, min_size = 80)) { cat("CACHED\n"); return(out_file) }
  adapter <- file.path(ADAPTERS, "run_nichenetr.R")
  status <- system2(CONFIG$exe$python, args = c(MEM_GUARD_PY, "--limit-gb", MEM_BUDGET_GB, "--", CONFIG$exe$rscript, adapter, ds_name, out_file),
    stdout = file.path(CONFIG$dir$log, sprintf("nichenetr_%s.log", ds_name)),
    stderr = file.path(CONFIG$dir$log, sprintf("nichenetr_%s.err", ds_name)))
  if (status == 0 && DiskOK(out_file, min_size = 80)) { cat("OK\n"); out_file } else { MarkFlatFail(out_file); cat("FAIL\n"); NA }
}
RunCytoTalk <- function(ds_name) {
  cat(sprintf("  [CytoTalk] %s ... ", ds_name))
  out_dir <- file.path(OUT, "CytoTalk", ds_name)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  marker <- file.path(out_dir, "cytotalk_done.txt")
  if (ReuseRegistered("CytoTalk", ds_name) && DirCached("CytoTalk", ds_name)) { cat("CACHED\n"); return(out_dir) }
  adapter <- file.path(ADAPTERS, "run_cytotalk.R")
  status <- system2(CONFIG$exe$python, args = c(MEM_GUARD_PY, "--limit-gb", MEM_BUDGET_GB, "--", CONFIG$exe$rscript, adapter, ds_name, out_dir),
    stdout = file.path(CONFIG$dir$log, sprintf("cytotalk_%s.log", ds_name)),
    stderr = file.path(CONFIG$dir$log, sprintf("cytotalk_%s.err", ds_name)))
  if (status == 0 && file.exists(marker)) { MarkDirDone("CytoTalk", ds_name, out_dir); cat("OK\n"); out_dir } else { cat("FAIL\n"); NA }
}
RunWOT <- function(ds_name) {
  cat(sprintf("  [wot] %s ... ", ds_name))
  out_dir <- file.path(OUT, "wot", ds_name)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  marker <- file.path(out_dir, "wot_done.txt")
  if (ReuseRegistered("wot", ds_name) && DirCached("wot", ds_name)) { cat("CACHED\n"); return(out_dir) }
  adapter <- file.path(ADAPTERS, "run_wot.py")
  status <- system2(CONFIG$exe$python,
    args = c(adapter, DatasetPath(ds_name), out_dir),
    stdout = file.path(CONFIG$dir$log, sprintf("wot_%s.log", ds_name)),
    stderr = file.path(CONFIG$dir$log, sprintf("wot_%s.err", ds_name)))
  if (status == 0 && file.exists(marker)) { MarkDirDone("wot", ds_name, out_dir); cat("OK\n"); out_dir } else { cat("FAIL\n"); NA }
}
RunNATMI <- function(ds_name) {
  cat(sprintf("  [NATMI] %s ... ", ds_name))
  out_file <- FlatOut(sprintf("natmi_%s_edges.csv", ds_name))
  if (ReuseRegistered("NATMI", ds_name) && DiskOK(out_file, min_size = 80)) { cat("CACHED\n"); return(out_file) }
  sp <- CONFIG$dataset_species[[ds_name]]; if (is.null(sp)) sp <- "human"
  adapter <- file.path(ADAPTERS, "run_natmi.py")
  status <- system2(CONFIG$exe$python,
    args = c(adapter, DatasetPath(ds_name), out_file, sp),
    stdout = file.path(CONFIG$dir$log, sprintf("natmi_%s.log", ds_name)),
    stderr = file.path(CONFIG$dir$log, sprintf("natmi_%s.err", ds_name)))
  if (status == 0 && DiskOK(out_file, min_size = 80)) { cat("OK\n"); out_file } else { MarkFlatFail(out_file); cat("FAIL\n"); NA }
}
RunGeneTrajectory <- function(ds_name) {
  cat(sprintf("  [GeneTrajectory] %s ... ", ds_name))
  out_file <- FlatOut(sprintf("genetrajectory_%s_gene_pt.csv", ds_name))
  if (ReuseRegistered("GeneTrajectory", ds_name) && DiskOK(out_file, min_size = 80)) { cat("CACHED\n"); return(out_file) }
  adapter <- file.path(ADAPTERS, "run_genetrajectory.py")
  status <- system2(CONFIG$exe$python,
    args = c(adapter, DatasetPath(ds_name), out_file),
    stdout = file.path(CONFIG$dir$log, sprintf("genetrajectory_%s.log", ds_name)),
    stderr = file.path(CONFIG$dir$log, sprintf("genetrajectory_%s.err", ds_name)))
  if (status == 0 && DiskOK(out_file, min_size = 80)) { cat("OK\n"); out_file } else { MarkFlatFail(out_file); cat("FAIL\n"); NA }
}
RunScTour <- function(ds_name) {
  cat(sprintf("  [scTour] %s ... ", ds_name))
  out_file <- FlatOut(sprintf("sctour_%s_pseudotime.csv", ds_name))
  if (ReuseRegistered("scTour", ds_name) && PseudotimeCsvOK(out_file, ds_name)) { cat("CACHED\n"); return(out_file) }
  MarkFlatFail(out_file)
  adapter <- file.path(ADAPTERS, "run_sctour.py")
  status <- system2(CONFIG$exe$python,
    args = c("-X", "utf8", adapter, DatasetPath(ds_name), out_file),
    stdout = file.path(CONFIG$dir$log, sprintf("sctour_%s.log", ds_name)),
    stderr = file.path(CONFIG$dir$log, sprintf("sctour_%s.err", ds_name)))
  if (status == 0 && PseudotimeCsvOK(out_file, ds_name)) { cat("OK\n"); out_file } else { MarkFlatFail(out_file); cat("FAIL\n"); NA }
}
RunStaVia <- function(ds_name) {
  cat(sprintf("  [StaVia] %s ... ", ds_name))
  out_file <- FlatOut(sprintf("stavia_%s_pseudotime.csv", ds_name))
  if (ReuseRegistered("StaVia", ds_name) && PseudotimeCsvOK(out_file, ds_name)) { cat("CACHED\n"); return(out_file) }
  MarkFlatFail(out_file)
  adapter <- file.path(ADAPTERS, "run_stavia.py")
  log_file <- file.path(CONFIG$dir$log, sprintf("stavia_%s.log", ds_name))
  err_file <- file.path(CONFIG$dir$log, sprintf("stavia_%s.err", ds_name))
  status <- system2(CONFIG$exe$python,
      args = c("-X", "utf8", adapter, DatasetPath(ds_name), out_file),
    stdout = log_file,
    stderr = err_file)
  ok <- (status == 0) && PseudotimeCsvOK(out_file, ds_name)
  if (ok) { cat("OK\n"); out_file } else { MarkFlatFail(out_file); cat("FAIL\n"); NA }
}
RunTraSig <- function(ds_name) {
  cat(sprintf("  [TraSig] %s ... ", ds_name))
  out_dir <- file.path(OUT, "TraSig", ds_name)
  run_trasig_adapter <- file.path(ADAPTERS, "run_trasig.R")
  if (!file.exists(run_trasig_adapter)) { cat("SKIP\n"); return(NA) }
  status <- system2(CONFIG$exe$python, args = c(MEM_GUARD_PY, "--limit-gb", MEM_BUDGET_GB, "--", CONFIG$exe$rscript, run_trasig_adapter, ds_name, out_dir),
    stdout = file.path(CONFIG$dir$log, sprintf("trasig_%s.log", ds_name)),
    stderr = file.path(CONFIG$dir$log, sprintf("trasig_%s.err", ds_name)))
  prep_ok <- file.exists(file.path(out_dir, "hid_var.pickle"))
  res_ok <- dir.exists(file.path(out_dir, "results"))
  if (status == 0 && res_ok) {
    cat("OK\n"); out_dir
  } else if (status == 0 && prep_ok) {
    # Prep produced hid_var.pickle but the actual TraSig permutation run did not
    # yield results/ -> this is NOT a completed competitor result. Report PARTIAL
    # (NA) so it is not counted as a successful run.
    cat("PARTIAL (prep only, no results)\n"); NA
  } else { cat("FAIL\n"); NA }
}

# ===== Main dispatch =====
RunBenchmark <- function() {
  cat("\n===== 01_run_other_algos.R =====\n")
  ds_names <- names(datasets)
  if (!is.null(filter_datasets)) ds_names <- intersect(filter_datasets, ds_names)
  cat(sprintf("Datasets: %d\n", length(ds_names)))

  # Execution order: CONFIG$algorithm_run_order (light → heavy by time/space complexity).
  formal_pool <- unique(c(
    CONFIG$algorithm_groups$ti_de_other,
    CONFIG$algorithm_groups$ti_trajectory,
    CONFIG$algorithm_groups$ccc_static,
    CONFIG$algorithm_groups$ccc_temporal))
  requested <- if (is.null(filter_algos)) {
    CONFIG$algorithm_groups$ti_trajectory
  } else {
    filter_algos
  }
  stopifnot(all(requested %in% formal_pool))
  order <- CONFIG$algorithm_run_order
  all_algos <- intersect(order, requested)
  tail <- setdiff(requested, all_algos)
  if (length(tail)) all_algos <- c(all_algos, tail)
  cat(sprintf("Algorithm order (%d): %s\n", length(all_algos), paste(all_algos, collapse = ", ")))

  runner_map <- list(
    Palantir = RunPalantir, LIANA = RunLIANA, cellrank2 = RunCellRank2,
    DPT = RunDPT, PAGA = RunPAGA, STREAM = RunSTREAM, Wishbone = RunWishbone,
    MultistageOT = RunMultistageOT,
    scFates = RunScFates,
    PHLOWER = RunPHLOWER,
    TICCI = RunTICCI, CellPhoneDB = RunCellPhoneDB, wot = RunWOT, NATMI = RunNATMI,
    GeneTrajectory = RunGeneTrajectory, scTour = RunScTour, StaVia = RunStaVia,
    TraSig = RunTraSig,
    nichenetr = RunNichenetr, CytoTalk = RunCytoTalk,
    monocle3 = function(d) RunRScript("monocle3", d, "run_monocle3.R"),
    slingshot = function(d) RunRScript("slingshot", d, "run_slingshot.R"),
    CellChat = function(d) RunRScript("CellChat", d, "run_cellchat.R"),
    Tempora = function(d) RunRScript("Tempora", d, "run_tempora.R"),
    SingleCellSignalR = function(d) RunRScript("SingleCellSignalR", d, "run_scsignalr.R"),
    tradeSeq = function(d) RunRScript("tradeSeq", d, "run_tradeseq.R"),
    TimeTalk = function(d) RunRScript("TimeTalk", d, "run_timetalk.R"),
    SCORPIUS = function(d) RunRScript("SCORPIUS", d, "run_scorpius.R"),
    TSCAN = function(d) RunRScript("TSCAN", d, "run_tscan.R"),
    SLICER = function(d) RunRScript("SLICER", d, "run_slicer.R"),
    FateID = function(d) RunRScript("FateID", d, "run_fateid.R"),
    pCreode = function(d) {
      cat(sprintf("  [pCreode] %s ... ", d))
      out_file <- FlatOut(sprintf("pcreode_%s_pseudotime.csv", d))
      if (ReuseRegistered("pCreode", d) && PseudotimeCsvOK(out_file, d)) { cat("CACHED\n"); return(out_file) }
      MarkFlatFail(out_file)
      adapter <- file.path(ADAPTERS, "run_pcreode.py")
      status <- system2(CONFIG$exe$python,
        args = c("-X", "utf8", adapter, DatasetPath(d), out_file),
        stdout = file.path(CONFIG$dir$log, sprintf("pcreode_%s.log", d)),
        stderr = file.path(CONFIG$dir$log, sprintf("pcreode_%s.err", d)))
      if (status == 0 && PseudotimeCsvOK(out_file, d)) { cat("OK\n"); out_file } else { MarkFlatFail(out_file); cat("FAIL\n"); NA }
    }
  )

  results <- data.table(
    algorithm = character(), dataset = character(), status = character(),
    registration_basis = character()
  )
  for (algo in all_algos) {
    runner <- runner_map[[algo]]
    if (is.null(runner)) { cat(sprintf("\n[%s] SKIP (no runner)\n", algo)); next }
    cat(sprintf("\n>>> %s <<<\n", algo))
    for (ds in ds_names) {
      if (PanelSkip(algo, ds)) {
        reason <- PanelSkipReason(algo, ds)
        cat(sprintf("  [%s] %s ... SKIP (%s)\n", algo, ds, reason))
        results <- rbind(results, data.table(algorithm=algo, dataset=ds, status="SKIP_PANEL"))
        next
      }
      if (!file.exists(DatasetPath(ds))) {
        # Record the missing dataset and continue the algorithm-by-dataset sweep.
        cat(sprintf("  [%s] %s ... SKIP (dataset file not found)\n", algo, ds))
        results <- rbind(results, data.table(algorithm=algo, dataset=ds, status="SKIP_NO_DATASET"))
        next
      }
      out <- runner(ds)
      status <- if (length(out) != 1L || is.na(out) || !nzchar(as.character(out))) "FAIL" else "OK"
      results <- rbind(
        results,
        data.table(
          algorithm = algo,
          dataset = ds,
          status = status,
          registration_basis = if (status == "OK") {
            "runner_verified_primary_artifact"
          } else {
            "runner_failed"
          }
        )
      )
    }
  }
  results[, `:=`(
    benchmark_protocol = CONFIG$llr$benchmark_protocol,
    result_revision = CONFIG$llr$result_revision
  )]
  registry_path <- file.path(
    CONFIG$dir$paper_tables,
    "_registry",
    sprintf("benchmark_results_%s.csv", CONFIG$llr$benchmark_protocol)
  )
  registry <- if (file.exists(registry_path)) fread(registry_path) else results[0L]
  registry <- registry[
    benchmark_protocol == CONFIG$llr$benchmark_protocol &
      result_revision == CONFIG$llr$result_revision
  ]
  registry <- rbindlist(
    list(registry[!results, on = .(algorithm, dataset)], results),
    use.names = TRUE,
    fill = TRUE
  )
  setorder(registry, algorithm, dataset)
  fwrite(registry, registry_path)
  cat("\n===== Summary =====\n")
  print(results[, .N, by = .(algorithm, status)])
}

RegisterValidatedResults <- function() {
  register_algos <- c(
    "slingshot", "DPT", "monocle3", "SCORPIUS", "MultistageOT", "scFates",
    "pCreode", "TSCAN", "Palantir", "StaVia", "scTour"
  )
  log_stems <- c(
    slingshot = "slingshot", DPT = "dpt", monocle3 = "monocle3",
    SCORPIUS = "SCORPIUS", MultistageOT = "multistageot", scFates = "scfates",
    pCreode = "pcreode", TSCAN = "TSCAN", Palantir = "palantir",
    StaVia = "stavia", scTour = "sctour"
  )
  success_markers <- c(
    slingshot = "OK:", DPT = "DPT OK", monocle3 = "OK:", SCORPIUS = "OK:",
    MultistageOT = "MultistageOT OK", scFates = "scFates OK", pCreode = "pCreode OK",
    TSCAN = "OK:", Palantir = "OK: Palantir", StaVia = "StaVia OK", scTour = "scTour OK"
  )
  LogOK <- function(algo, ds_name) {
    log_path <- file.path(CONFIG$dir$log, sprintf("%s_%s.log", log_stems[[algo]], ds_name))
    if (!DiskOK(log_path, min_size = 1)) return(FALSE)
    txt <- paste(readLines(log_path, encoding = "UTF-8", warn = FALSE), collapse = "\n")
    grepl(success_markers[[algo]], txt, fixed = TRUE) && !grepl("FAIL", txt, fixed = TRUE)
  }
  ArtifactOK <- function(algo, ds_name) {
    if (algo %in% names(R_PSEUDOTIME_OUTPUT)) {
      return(DirCached(algo, ds_name))
    }
    if (algo == "Palantir") {
      return(DiskOK(FlatOut(sprintf("palantir_%s_result.h5ad", ds_name))))
    }
    PseudotimeCsvOK(FlatOut(sprintf("%s_%s_pseudotime.csv", tolower(algo), ds_name)), ds_name)
  }
  rows <- rbindlist(lapply(register_algos, function(algo) {
    rbindlist(lapply(BENCHMARK_COMPLETE_CORE_DATASETS, function(ds_name) {
      data.table(
        algorithm = algo,
        dataset = ds_name,
        status = if (ArtifactOK(algo, ds_name) && LogOK(algo, ds_name)) "OK" else "UNVERIFIED"
      )
    }))
  }))
  rows[, `:=`(
    benchmark_protocol = CONFIG$llr$benchmark_protocol,
    result_revision = CONFIG$llr$result_revision,
    registration_basis = "validated_primary_artifact_and_success_log"
  )]
  registry <- if (file.exists(RUN_REGISTRY_PATH)) fread(RUN_REGISTRY_PATH) else rows[0L]
  registry <- registry[
    benchmark_protocol == CONFIG$llr$benchmark_protocol &
      result_revision == CONFIG$llr$result_revision
  ]
  registry <- rbindlist(
    list(registry[!rows, on = .(algorithm, dataset)], rows),
    use.names = TRUE,
    fill = TRUE
  )
  setorder(registry, algorithm, dataset)
  fwrite(registry, RUN_REGISTRY_PATH)
  print(rows[, .N, by = .(algorithm, status)][order(algorithm, status)])
  stopifnot(all(rows$status == "OK"))
}

if (register_validated) {
  RegisterValidatedResults()
} else {
  RunBenchmark()
}
