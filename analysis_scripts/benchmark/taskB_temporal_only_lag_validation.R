# Task B: ligand→receptor peak-bin lag for temporal-only headline pairs
# Interpreters: D:\R\R-4.5.3\bin\Rscript.exe --encoding=UTF-8
# path bootstrap (relative; set LLR_ASSETS to override)
.llr_cmd <- commandArgs(trailingOnly = FALSE)
.llr_file <- sub("^--file=", "", grep("^--file=", .llr_cmd, value = TRUE)[1])
if (!length(.llr_file) || is.na(.llr_file) || !nzchar(.llr_file)) {
  stop("Run via Rscript so --file= is set")
}
.llr_here <- dirname(normalizePath(.llr_file, winslash = "/", mustWork = TRUE))
source(file.path(.llr_here, "../config/llr_bootstrap.R"), local = FALSE)
.llr_assets <- llr_assets_dir()

set.seed(12345)
Sys.setenv(OMP_NUM_THREADS = parallel::detectCores())
data.table::setDTthreads(parallel::detectCores())

library(data.table)
library(Matrix)

CONFIG <- list(
  paths = list(
    workspace = file.path(.llr_assets, "cache/linealinkr_workspace"),
    temporal_only = file.path(.llr_assets, "results/paper_tables/supplement/tables/temporal_only_validation"),
    out_dir = file.path(.llr_assets, "results/paper_tables/supplement/tables/temporal_only_orthogonal_validation")
  ),
  threads = parallel::detectCores(),
  seed = 12345,
  n_bins = 12L,
  n_perm = 500L,
  datasets = c("pancreas", "dentategyrus")
)

dir.create(CONFIG$paths$out_dir, recursive = TRUE, showWarnings = FALSE)

LoadExpr <- function(ds) {
  extract_dir <- file.path(CONFIG$paths$workspace, ds, "extract")
  genes <- readLines(file.path(extract_dir, paste0(ds, "_genes.txt")), warn = FALSE)
  cells <- readLines(file.path(extract_dir, paste0(ds, "_cells.txt")), warn = FALSE)
  mat <- Matrix::readMM(file.path(extract_dir, paste0(ds, "_matrix.mtx")))
  if (nrow(mat) != length(genes) && ncol(mat) == length(genes)) {
    mat <- t(mat)
  }
  rownames(mat) <- genes
  colnames(mat) <- cells
  mat
}

LoadPt <- function(ds) {
  rds <- readRDS(file.path(CONFIG$paths$workspace, ds, "linealinkr_complete.rds"))
  pt <- as.numeric(rds$pseudotime)
  names(pt) <- as.character(rds$cells)
  list(pt = pt, n_lineages = rds$n_lineages, lineage_weights = rds$lineage_weights)
}

BinMeans <- function(expr_row, pt, n_bins) {
  ok <- is.finite(expr_row) & is.finite(pt)
  x <- expr_row[ok]
  t <- pt[ok]
  br <- quantile(t, probs = seq(0, 1, length.out = n_bins + 1L), na.rm = TRUE)
  br <- unique(br)
  if (length(br) < 3L) {
    return(rep(NA_real_, n_bins))
  }
  bin <- as.integer(cut(t, breaks = br, include.lowest = TRUE, labels = FALSE))
  vapply(seq_len(n_bins), function(b) {
    idx <- which(bin == b)
    if (length(idx) == 0L) NA_real_ else mean(x[idx])
  }, numeric(1))
}

PeakBin <- function(v) {
  if (all(!is.finite(v))) return(NA_integer_)
  as.integer(which.max(replace(v, !is.finite(v), -Inf)))
}

CrossCorrLag <- function(x, y, max_lag = 5L) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  if (length(x) < 4L) {
    return(list(best_lag = NA_integer_, best_cor = NA_real_))
  }
  x <- as.numeric(scale(x))
  y <- as.numeric(scale(y))
  lags <- (-max_lag):max_lag
  cors <- vapply(lags, function(L) {
    if (L < 0L) {
      xs <- x[seq_len(length(x) + L)]
      ys <- y[seq_len(length(y) + L) - L]
    } else if (L > 0L) {
      xs <- x[seq_len(length(x) - L) + L]
      ys <- y[seq_len(length(y) - L)]
    } else {
      xs <- x
      ys <- y
    }
    if (length(xs) < 3L || sd(xs) == 0 || sd(ys) == 0) return(NA_real_)
    cor(xs, ys)
  }, numeric(1))
  i <- which.max(replace(cors, !is.finite(cors), -Inf))
  list(best_lag = as.integer(lags[i]), best_cor = as.numeric(cors[i]))
}

PermPLagPositive <- function(x, y, n_perm, seed) {
  set.seed(seed)
  obs <- PeakBin(y) - PeakBin(x)
  if (!is.finite(obs)) return(NA_real_)
  nulls <- replicate(n_perm, {
    yp <- sample(y)
    PeakBin(yp) - PeakBin(x)
  })
  mean(nulls >= obs)
}

rows <- list()
profile_rows <- list()

for (ds in CONFIG$datasets) {
  to_path <- file.path(CONFIG$paths$temporal_only, paste0(ds, "_temporal_only_vs_external.csv"))
  to <- fread(to_path, encoding = "UTF-8")
  expr <- LoadExpr(ds)
  pt_obj <- LoadPt(ds)
  pt <- pt_obj$pt
  common <- intersect(colnames(expr), names(pt))
  expr <- expr[, common, drop = FALSE]
  pt <- pt[common]

  for (i in seq_len(nrow(to))) {
    lig <- as.character(to$ligand_genes[i])
    rec <- as.character(to$receptor_genes[i])
    pair <- as.character(to$pair_name[i])
    lig_genes <- unlist(strsplit(lig, "[,_|; ]+"))
    rec_genes <- unlist(strsplit(rec, "[,_|; ]+"))
    lig_genes <- lig_genes[lig_genes %in% rownames(expr)]
    rec_genes <- rec_genes[rec_genes %in% rownames(expr)]
    if (length(lig_genes) == 0L || length(rec_genes) == 0L) {
      rows[[length(rows) + 1L]] <- data.table(
        dataset = ds, pair_name = pair, ligand_genes = lig, receptor_genes = rec,
        family = to$family[i], temporal_mean_score = to$temporal_mean_score[i],
        ligand_peak_bin = NA_integer_, receptor_peak_bin = NA_integer_,
        lag_bins = NA_integer_, lag_positive = NA, best_xcorr_lag = NA_integer_,
        best_xcorr = NA_real_, perm_p_lag_positive = NA_real_,
        status = "gene_missing_in_matrix"
      )
      next
    }
    lig_vec <- if (length(lig_genes) == 1L) {
      as.numeric(expr[lig_genes, ])
    } else {
      as.numeric(Matrix::colMeans(expr[lig_genes, , drop = FALSE]))
    }
    rec_vec <- if (length(rec_genes) == 1L) {
      as.numeric(expr[rec_genes, ])
    } else {
      as.numeric(Matrix::colMeans(expr[rec_genes, , drop = FALSE]))
    }
    lig_bin <- BinMeans(lig_vec, pt, CONFIG$n_bins)
    rec_bin <- BinMeans(rec_vec, pt, CONFIG$n_bins)
    lp <- PeakBin(lig_bin)
    rp <- PeakBin(rec_bin)
    lag <- as.integer(rp - lp)
    xc <- CrossCorrLag(lig_bin, rec_bin, max_lag = 5L)
    pp <- PermPLagPositive(lig_bin, rec_bin, CONFIG$n_perm, CONFIG$seed + i)
    rows[[length(rows) + 1L]] <- data.table(
      dataset = ds, pair_name = pair, ligand_genes = lig, receptor_genes = rec,
      family = to$family[i], temporal_mean_score = to$temporal_mean_score[i],
      ligand_peak_bin = lp, receptor_peak_bin = rp,
      lag_bins = lag, lag_positive = is.finite(lag) & lag > 0L,
      best_xcorr_lag = xc$best_lag, best_xcorr = xc$best_cor,
      perm_p_lag_positive = pp, status = "ok"
    )
    profile_rows[[length(profile_rows) + 1L]] <- data.table(
      dataset = ds, pair_name = pair, bin = seq_len(CONFIG$n_bins),
      ligand_mean = lig_bin, receptor_mean = rec_bin
    )
  }
}

lag_dt <- rbindlist(rows, fill = TRUE)
prof_dt <- rbindlist(profile_rows, fill = TRUE)
fwrite(lag_dt, file.path(CONFIG$paths$out_dir, "temporal_only_ligand_receptor_lag.csv"))
fwrite(prof_dt, file.path(CONFIG$paths$out_dir, "temporal_only_bin_expression_profiles.csv"))

summary_dt <- lag_dt[status == "ok", .(
  n_pairs = .N,
  n_lag_positive = sum(lag_positive, na.rm = TRUE),
  frac_lag_positive = mean(lag_positive, na.rm = TRUE),
  median_lag = median(lag_bins, na.rm = TRUE),
  n_xcorr_lag_positive = sum(best_xcorr_lag > 0, na.rm = TRUE),
  median_perm_p = median(perm_p_lag_positive, na.rm = TRUE)
), by = dataset]
fwrite(summary_dt, file.path(CONFIG$paths$out_dir, "temporal_only_lag_summary_by_dataset.csv"))

cat("Wrote lag table n=", nrow(lag_dt), "\n")
print(lag_dt[, .(dataset, pair_name, ligand_peak_bin, receptor_peak_bin, lag_bins, lag_positive, best_xcorr_lag, perm_p_lag_positive, status)])
print(summary_dt)
