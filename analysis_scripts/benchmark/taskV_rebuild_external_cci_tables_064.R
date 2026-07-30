.llr_cmd <- commandArgs(trailingOnly = FALSE)
.llr_file <- sub("^--file=", "", grep("^--file=", .llr_cmd, value = TRUE)[1])
.llr_here <- dirname(normalizePath(.llr_file, winslash = "/", mustWork = TRUE))
.llr_root <- normalizePath(
  file.path(.llr_here, "..", ".."),
  winslash = "/",
  mustWork = TRUE
)
.llr_config_env <- new.env(parent = globalenv())
sys.source(
  file.path(.llr_root, "pipeline", "config", "00_config.R"),
  envir = .llr_config_env
)

CONFIG <- list(
  paths = list(
    root = .llr_root,
    rlib = file.path(.llr_root, "rlib"),
    workspace = .llr_config_env$CONFIG$dir$linealinkr,
    cellchat = file.path(
      .llr_config_env$CONFIG$dir$algorithm_out, "CellChat"
    ),
    liana = .llr_config_env$CONFIG$dir$algorithm_misc,
    canonical = .llr_config_env$CONFIG$dataset_registry[c(
      "pancreas", "dentategyrus"
    )],
    output = .llr_config_env$CONFIG$dir$supplement_tables
  ),
  protocol = .llr_config_env$CONFIG$llr$benchmark_protocol,
  revision = .llr_config_env$CONFIG$llr$result_revision,
  external_source = "canonical_native_cellchat_liana_outputs",
  threads = parallel::detectCores(),
  seed = 12345L
)
rm(.llr_config_env)

set.seed(CONFIG$seed)
Sys.setenv(OMP_NUM_THREADS = CONFIG$threads)
data.table::setDTthreads(CONFIG$threads)
.libPaths(c(CONFIG$paths$rlib, .libPaths()))
library(data.table)

required_inputs <- c(
  file.path(CONFIG$paths$workspace, "pancreas", "linealinkr_complete.rds"),
  file.path(CONFIG$paths$workspace, "dentategyrus", "linealinkr_complete.rds"),
  file.path(CONFIG$paths$workspace, "GSE210152_athero", "linealinkr_complete.rds"),
  file.path(CONFIG$paths$cellchat, "pancreas", "cellchat.rds"),
  file.path(CONFIG$paths$cellchat, "dentategyrus", "cellchat.rds"),
  file.path(CONFIG$paths$liana, "liana_pancreas_result.csv"),
  file.path(CONFIG$paths$liana, "liana_pancreas_cells.csv"),
  file.path(CONFIG$paths$liana, "liana_pancreas_run.json"),
  file.path(CONFIG$paths$liana, "liana_dentategyrus_result.csv"),
  file.path(CONFIG$paths$liana, "liana_dentategyrus_cells.csv"),
  file.path(CONFIG$paths$liana, "liana_dentategyrus_run.json")
)
stopifnot(all(file.exists(required_inputs)))

pan_res <- readRDS(required_inputs[[1L]])
dent_res <- readRDS(required_inputs[[2L]])
gse_res <- readRDS(required_inputs[[3L]])
pan_cellchat_object <- readRDS(required_inputs[[4L]])
dent_cellchat_object <- readRDS(required_inputs[[5L]])
stopifnot(
  identical(pan_res$benchmark_protocol, CONFIG$protocol),
  identical(pan_res$result_revision, CONFIG$revision),
  identical(dent_res$benchmark_protocol, CONFIG$protocol),
  identical(dent_res$result_revision, CONFIG$revision),
  identical(gse_res$benchmark_protocol, CONFIG$protocol),
  identical(gse_res$result_revision, CONFIG$revision),
  identical(as.character(pan_res$cells), as.character(colnames(pan_cellchat_object@data))),
  identical(as.character(dent_res$cells), as.character(colnames(dent_cellchat_object@data)))
)

out_cmp <- file.path(CONFIG$paths$output, "external_cci_comparison")
out_temporal <- file.path(CONFIG$paths$output, "temporal_only_validation")
out_case <- file.path(CONFIG$paths$output, "case_study_extension")
out_narrative <- file.path(CONFIG$paths$output, "linealinkr_narrative")
dir.create(out_cmp, recursive = TRUE, showWarnings = FALSE)
dir.create(out_temporal, recursive = TRUE, showWarnings = FALSE)
dir.create(out_case, recursive = TRUE, showWarnings = FALSE)
dir.create(out_narrative, recursive = TRUE, showWarnings = FALSE)

cellchat_columns <- c(
  "source", "target", "ligand", "receptor", "prob", "pval"
)
liana_columns <- c(
  "source", "target", "ligand_complex", "receptor_complex",
  "magnitude_rank", "specificity_rank", "lr_means", "expr_prod",
  "scaled_weight", "lr_logfc", "spec_weight", "lrscore"
)
liana_method_columns <- c(
  "lr_means", "expr_prod", "scaled_weight", "lr_logfc", "spec_weight",
  "lrscore"
)

pan_scores_all <- as.data.table(pan_res$temporal_lr$scores)
pan_scores <- pan_scores_all[
  is.finite(score) & expression_supported == TRUE
]
pan_scores <- merge(
  pan_scores,
  as.data.table(pan_res$pair_gene_map),
  by = "pair_idx",
  all.x = TRUE
)
pan_scores[, pair_name := paste(ligand_genes, receptor_genes, sep = "_")]
pan_scores[, pair_key := paste0(
  tolower(gsub("[^A-Za-z0-9]+", "", ligand_genes)),
  "___",
  tolower(gsub("[^A-Za-z0-9]+", "", receptor_genes))
)]
pan_llr <- pan_scores[, .(
  pair_name = pair_name[1],
  ligand_genes = ligand_genes[1],
  receptor_genes = receptor_genes[1],
  temporal_mean_score = mean(score),
  temporal_max_score = max(score),
  n_bins_detected = uniqueN(bin[score > 0]),
  n_lineages_detected = uniqueN(lineage[score > 0]),
  n_source_target = uniqueN(paste(source, target, sep = "->"))
), by = pair_key]
pan_lineage <- pan_scores[, .(
  score = mean(score),
  n_rows = .N,
  n_bins = uniqueN(bin),
  n_source_target = uniqueN(paste(source, target, sep = "->"))
), by = .(
  pair_key, pair_name, ligand_genes, receptor_genes, lineage
)]
pan_lineage[, `:=`(
  dataset = "pancreas",
  benchmark_protocol = CONFIG$protocol,
  result_revision = CONFIG$revision
)]
setcolorder(pan_lineage, c(
  "dataset", "pair_key", "pair_name", "ligand_genes", "receptor_genes",
  "lineage", "score", "n_rows", "n_bins", "n_source_target",
  "benchmark_protocol", "result_revision"
))
setorder(pan_lineage, pair_key, lineage)
fwrite(
  pan_lineage,
  file.path(out_narrative, "pancreas_temporal_lr_lineage_scores.csv")
)
pan_llr[, family := fcase(
  grepl("^(Wnt)|_(Fzd)", pair_name, ignore.case = TRUE), "Wnt-Fzd",
  grepl("^(Dll|Jag)|_(Notch)", pair_name, ignore.case = TRUE), "Notch",
  grepl("^(Fgf)|_(Fgfr)", pair_name, ignore.case = TRUE), "FGF",
  grepl("^(Bmp|Gdf|Tgfb)|_(Bmpr|Acvr|Tgfbr)", pair_name, ignore.case = TRUE), "BMP-TGF",
  grepl("^(Efna|Efnb)|_(Epha|Ephb)", pair_name, ignore.case = TRUE), "Ephrin-Eph",
  grepl("^(Sema)|_(Plxn|Nrp)", pair_name, ignore.case = TRUE), "Semaphorin",
  grepl("^(Mdk)|_(Ncl|Sdc)", pair_name, ignore.case = TRUE), "MDK-axis",
  grepl("^(Thbs)", ligand_genes, ignore.case = TRUE), "Thrombospondin",
  grepl("^(Lama|Lamb|Lamc|Agrn|Hspg2|Col)", ligand_genes, ignore.case = TRUE), "ECM-adhesion",
  grepl("^(Gas6)|_(Mertk|Tyro3|Axl)", pair_name, ignore.case = TRUE), "GAS6-TAM",
  grepl("^(Ccl|Cxcl)|_(Ackr|Ccr|Cxcr)", pair_name, ignore.case = TRUE), "Chemokine",
  grepl("^(Guca2)|_(Gucy)", pair_name, ignore.case = TRUE), "Guanylin",
  grepl("^(App)|_(Cd74)", pair_name, ignore.case = TRUE), "APP-CD74",
  grepl("^(Cdh|Dsc|Dsg|F11r|Ocln|Mpz|Mpzl)", pair_name, ignore.case = TRUE), "Adhesion-junction",
  default = "other"
)]

stopifnot(identical(
  as.character(pan_res$cells),
  as.character(colnames(pan_cellchat_object@data))
))
pan_cellchat <- as.data.table(CellChat::subsetCommunication(
  pan_cellchat_object
))
stopifnot(all(cellchat_columns %in% names(pan_cellchat)))
pan_cellchat[, pair_key := paste0(
  tolower(gsub("[^A-Za-z0-9]+", "", ligand)),
  "___",
  tolower(gsub("[^A-Za-z0-9]+", "", receptor))
)]
pan_cellchat <- pan_cellchat[
  is.finite(prob) & is.finite(pval) & pval <= 0.05,
  .(
    cellchat_mean_prob = mean(prob),
    cellchat_max_prob = max(prob),
    cellchat_n_source_target = uniqueN(paste(source, target, sep = "->"))
  ),
  by = pair_key
]

pan_liana_cells <- fread(
  file.path(CONFIG$paths$liana, "liana_pancreas_cells.csv"),
  colClasses = c(cell_id = "character"),
  encoding = "UTF-8"
)
pan_liana_run <- jsonlite::fromJSON(
  file.path(CONFIG$paths$liana, "liana_pancreas_run.json"),
  simplifyVector = TRUE
)
stopifnot(
  identical(names(pan_liana_cells), "cell_id"),
  !anyDuplicated(pan_liana_cells$cell_id),
  identical(as.character(pan_res$cells), pan_liana_cells$cell_id),
  identical(
    as.character(pan_liana_run$content_contract),
    "liana_rank_aggregate_native_v1"
  ),
  identical(as.character(pan_liana_run$dataset), "pancreas"),
  identical(
    normalizePath(
      as.character(pan_liana_run$canonical_source_path),
      winslash = "/",
      mustWork = TRUE
    ),
    normalizePath(
      CONFIG$paths$canonical[["pancreas"]],
      winslash = "/",
      mustWork = TRUE
    )
  ),
  identical(as.integer(pan_liana_run$n_cells), as.integer(pan_res$n_cells)),
  identical(as.integer(pan_liana_run$n_genes), as.integer(pan_res$n_genes)),
  identical(as.integer(pan_liana_run$n_cells), nrow(pan_liana_cells)),
  nzchar(as.character(pan_liana_run$groupby)),
  nzchar(as.character(pan_liana_run$resource))
)
pan_liana <- fread(
  file.path(CONFIG$paths$liana, "liana_pancreas_result.csv"),
  encoding = "UTF-8"
)
stopifnot(all(liana_columns %in% names(pan_liana)))
pan_liana[, pair_key := paste0(
  tolower(gsub("[^A-Za-z0-9]+", "", ligand_complex)),
  "___",
  tolower(gsub("[^A-Za-z0-9]+", "", receptor_complex))
)]
pan_liana[, methods_supported := rowSums(
  is.finite(as.matrix(.SD))
), .SDcols = liana_method_columns]
pan_liana <- pan_liana[is.finite(magnitude_rank), .(
  liana_best_magnitude_rank = min(magnitude_rank),
  liana_best_specificity_rank = min(c(
    specificity_rank[is.finite(specificity_rank)], Inf
  )),
  liana_n_source_target = uniqueN(paste(source, target, sep = "->")),
  liana_max_methods_supported = max(methods_supported)
), by = pair_key]
pan_liana[!is.finite(liana_best_specificity_rank),
  liana_best_specificity_rank := NA_real_]
setorder(pan_liana, liana_best_magnitude_rank, liana_best_specificity_rank)
pan_liana[, liana_pair_rank := seq_len(.N)]
pan_liana[, liana_top10pct :=
  liana_pair_rank <= max(1L, ceiling(.N * 0.10))]
pan_liana <- pan_liana[, .(
  pair_key, liana_best_magnitude_rank, liana_best_specificity_rank,
  liana_n_source_target, liana_max_methods_supported,
  liana_pair_rank, liana_top10pct
)]

pan_cmp <- merge(pan_llr, pan_cellchat, by = "pair_key", all.x = TRUE)
pan_cmp <- merge(pan_cmp, pan_liana, by = "pair_key", all.x = TRUE)
pan_cmp[, cellchat_supported := is.finite(cellchat_mean_prob)]
pan_cmp[, liana_supported := is.finite(liana_best_magnitude_rank)]
pan_cmp[is.na(liana_top10pct), liana_top10pct := FALSE]
pan_cmp[, external_support_count :=
  as.integer(cellchat_supported) + as.integer(liana_supported)]
pan_cmp[, comparator_membership := fifelse(
  cellchat_supported & liana_supported,
  "cellchat_and_liana",
  fifelse(
    cellchat_supported,
    "cellchat_only",
    fifelse(liana_supported, "liana_only", "none")
  )
)]
setorder(
  pan_cmp,
  -temporal_mean_score,
  -temporal_max_score,
  -n_bins_detected,
  -n_lineages_detected
)
pan_cmp[, temporal_pair_rank := seq_len(.N)]
pan_cmp[, temporal_top10pct := temporal_pair_rank <= max(1L, ceiling(.N * 0.10))]
pan_cmp[, temporal_broad := n_bins_detected >= 2L & n_lineages_detected >= 1L]
pan_cmp[, temporal_candidate := temporal_top10pct & temporal_broad]
pan_cmp[, temporal_only_vs_external :=
  temporal_candidate & !cellchat_supported & !liana_top10pct]
pan_cmp[, benchmark_protocol := CONFIG$protocol]
pan_cmp[, result_revision := CONFIG$revision]
pan_cmp[, external_source := CONFIG$external_source]
fwrite(pan_cmp, file.path(out_cmp, "pancreas_llr_vs_external_pair_table.csv"))

pan_temporal_only <- pan_cmp[temporal_only_vs_external == TRUE]
fwrite(
  pan_temporal_only,
  file.path(out_temporal, "pancreas_temporal_only_vs_external.csv")
)
pan_family <- pan_temporal_only[, .(
  n_pairs = .N,
  mean_temporal_score = mean(temporal_mean_score),
  max_temporal_score = max(temporal_max_score),
  dataset = "pancreas"
), by = family]
fwrite(
  pan_family,
  file.path(out_temporal, "pancreas_temporal_only_family_summary.csv")
)
pan_case <- pan_cmp[temporal_candidate == TRUE, .(
  pair_name, family, temporal_mean_score, temporal_max_score,
  n_bins_detected, n_lineages_detected, temporal_pair_rank,
  comparator_membership, cellchat_supported, liana_supported, liana_top10pct,
  benchmark_protocol, result_revision, external_source
)]
fwrite(pan_case, file.path(out_case, "pancreas_case_support_pairs.csv"))

pan_peak <- pan_scores[order(-score), .SD[1L], by = pair_key][, .(
  pair_key,
  peak_score = score,
  peak_bin = bin,
  peak_lineage = lineage,
  peak_source = source,
  peak_target = target
)]
pan_top <- merge(pan_cmp, pan_peak, by = "pair_key", all.x = TRUE)
setorder(pan_top, temporal_pair_rank)
pan_top <- head(pan_top, 50L)
pan_top[, dataset := "pancreas"]
fwrite(pan_top, file.path(out_narrative, "pancreas_top_temporal_lr_pairs.csv"))
pan_summary <- data.table(
  dataset = "pancreas",
  benchmark_protocol = pan_res$benchmark_protocol,
  result_revision = pan_res$result_revision,
  n_cells = as.integer(pan_res$n_cells),
  n_genes = as.integer(pan_res$n_genes),
  n_types = as.integer(pan_res$n_types),
  n_lineages = as.integer(pan_res$n_lineages),
  n_lr_pairs = as.integer(nrow(pan_res$lr_scores)),
  n_temporal_rows = as.integer(nrow(pan_scores_all)),
  n_temporal_pairs = as.integer(uniqueN(pan_scores$pair_idx)),
  mean_temporal_score = mean(pan_scores$score),
  max_temporal_score = max(pan_scores$score)
)
fwrite(pan_summary, file.path(out_narrative, "pancreas_dataset_summary.csv"))

dent_scores_all <- as.data.table(dent_res$temporal_lr$scores)
dent_scores <- dent_scores_all[
  is.finite(score) & expression_supported == TRUE
]
dent_scores <- merge(
  dent_scores,
  as.data.table(dent_res$pair_gene_map),
  by = "pair_idx",
  all.x = TRUE
)
dent_scores[, pair_name := paste(ligand_genes, receptor_genes, sep = "_")]
dent_scores[, pair_key := paste0(
  tolower(gsub("[^A-Za-z0-9]+", "", ligand_genes)),
  "___",
  tolower(gsub("[^A-Za-z0-9]+", "", receptor_genes))
)]
dent_llr <- dent_scores[, .(
  pair_name = pair_name[1],
  ligand_genes = ligand_genes[1],
  receptor_genes = receptor_genes[1],
  temporal_mean_score = mean(score),
  temporal_max_score = max(score),
  n_bins_detected = uniqueN(bin[score > 0]),
  n_lineages_detected = uniqueN(lineage[score > 0]),
  n_source_target = uniqueN(paste(source, target, sep = "->"))
), by = pair_key]
dent_llr[, family := fcase(
  grepl("^(Wnt)|_(Fzd)", pair_name, ignore.case = TRUE), "Wnt-Fzd",
  grepl("^(Dll|Jag)|_(Notch)", pair_name, ignore.case = TRUE), "Notch",
  grepl("^(Fgf)|_(Fgfr)", pair_name, ignore.case = TRUE), "FGF",
  grepl("^(Bmp|Gdf|Tgfb)|_(Bmpr|Acvr|Tgfbr)", pair_name, ignore.case = TRUE), "BMP-TGF",
  grepl("^(Efna|Efnb)|_(Epha|Ephb)", pair_name, ignore.case = TRUE), "Ephrin-Eph",
  grepl("^(Sema)|_(Plxn|Nrp)", pair_name, ignore.case = TRUE), "Semaphorin",
  grepl("^(Mdk)|_(Ncl|Sdc)", pair_name, ignore.case = TRUE), "MDK-axis",
  grepl("^(Thbs)", ligand_genes, ignore.case = TRUE), "Thrombospondin",
  grepl("^(Lama|Lamb|Lamc|Agrn|Hspg2|Col)", ligand_genes, ignore.case = TRUE), "ECM-adhesion",
  grepl("^(Gas6)|_(Mertk|Tyro3|Axl)", pair_name, ignore.case = TRUE), "GAS6-TAM",
  grepl("^(Ccl|Cxcl)|_(Ackr|Ccr|Cxcr)", pair_name, ignore.case = TRUE), "Chemokine",
  grepl("^(Guca2)|_(Gucy)", pair_name, ignore.case = TRUE), "Guanylin",
  grepl("^(App)|_(Cd74)", pair_name, ignore.case = TRUE), "APP-CD74",
  grepl("^(Cdh|Dsc|Dsg|F11r|Ocln|Mpz|Mpzl)", pair_name, ignore.case = TRUE), "Adhesion-junction",
  default = "other"
)]

stopifnot(identical(
  as.character(dent_res$cells),
  as.character(colnames(dent_cellchat_object@data))
))
dent_cellchat <- as.data.table(CellChat::subsetCommunication(
  dent_cellchat_object
))
stopifnot(all(cellchat_columns %in% names(dent_cellchat)))
dent_cellchat[, pair_key := paste0(
  tolower(gsub("[^A-Za-z0-9]+", "", ligand)),
  "___",
  tolower(gsub("[^A-Za-z0-9]+", "", receptor))
)]
dent_cellchat <- dent_cellchat[
  is.finite(prob) & is.finite(pval) & pval <= 0.05,
  .(
    cellchat_mean_prob = mean(prob),
    cellchat_max_prob = max(prob),
    cellchat_n_source_target = uniqueN(paste(source, target, sep = "->"))
  ),
  by = pair_key
]

dent_liana_cells <- fread(
  file.path(CONFIG$paths$liana, "liana_dentategyrus_cells.csv"),
  colClasses = c(cell_id = "character"),
  encoding = "UTF-8"
)
dent_liana_run <- jsonlite::fromJSON(
  file.path(CONFIG$paths$liana, "liana_dentategyrus_run.json"),
  simplifyVector = TRUE
)
stopifnot(
  identical(names(dent_liana_cells), "cell_id"),
  !anyDuplicated(dent_liana_cells$cell_id),
  identical(as.character(dent_res$cells), dent_liana_cells$cell_id),
  identical(
    as.character(dent_liana_run$content_contract),
    "liana_rank_aggregate_native_v1"
  ),
  identical(as.character(dent_liana_run$dataset), "dentategyrus"),
  identical(
    normalizePath(
      as.character(dent_liana_run$canonical_source_path),
      winslash = "/",
      mustWork = TRUE
    ),
    normalizePath(
      CONFIG$paths$canonical[["dentategyrus"]],
      winslash = "/",
      mustWork = TRUE
    )
  ),
  identical(
    as.integer(dent_liana_run$n_cells),
    as.integer(dent_res$n_cells)
  ),
  identical(
    as.integer(dent_liana_run$n_genes),
    as.integer(dent_res$n_genes)
  ),
  identical(as.integer(dent_liana_run$n_cells), nrow(dent_liana_cells)),
  nzchar(as.character(dent_liana_run$groupby)),
  nzchar(as.character(dent_liana_run$resource))
)
dent_liana <- fread(
  file.path(CONFIG$paths$liana, "liana_dentategyrus_result.csv"),
  encoding = "UTF-8"
)
stopifnot(all(liana_columns %in% names(dent_liana)))
dent_liana[, pair_key := paste0(
  tolower(gsub("[^A-Za-z0-9]+", "", ligand_complex)),
  "___",
  tolower(gsub("[^A-Za-z0-9]+", "", receptor_complex))
)]
dent_liana[, methods_supported := rowSums(
  is.finite(as.matrix(.SD))
), .SDcols = liana_method_columns]
dent_liana <- dent_liana[is.finite(magnitude_rank), .(
  liana_best_magnitude_rank = min(magnitude_rank),
  liana_best_specificity_rank = min(c(
    specificity_rank[is.finite(specificity_rank)], Inf
  )),
  liana_n_source_target = uniqueN(paste(source, target, sep = "->")),
  liana_max_methods_supported = max(methods_supported)
), by = pair_key]
dent_liana[!is.finite(liana_best_specificity_rank),
  liana_best_specificity_rank := NA_real_]
setorder(dent_liana, liana_best_magnitude_rank, liana_best_specificity_rank)
dent_liana[, liana_pair_rank := seq_len(.N)]
dent_liana[, liana_top10pct :=
  liana_pair_rank <= max(1L, ceiling(.N * 0.10))]
dent_liana <- dent_liana[, .(
  pair_key, liana_best_magnitude_rank, liana_best_specificity_rank,
  liana_n_source_target, liana_max_methods_supported,
  liana_pair_rank, liana_top10pct
)]

dent_cmp <- merge(dent_llr, dent_cellchat, by = "pair_key", all.x = TRUE)
dent_cmp <- merge(dent_cmp, dent_liana, by = "pair_key", all.x = TRUE)
dent_cmp[, cellchat_supported := is.finite(cellchat_mean_prob)]
dent_cmp[, liana_supported := is.finite(liana_best_magnitude_rank)]
dent_cmp[is.na(liana_top10pct), liana_top10pct := FALSE]
dent_cmp[, external_support_count :=
  as.integer(cellchat_supported) + as.integer(liana_supported)]
dent_cmp[, comparator_membership := fifelse(
  cellchat_supported & liana_supported,
  "cellchat_and_liana",
  fifelse(
    cellchat_supported,
    "cellchat_only",
    fifelse(liana_supported, "liana_only", "none")
  )
)]
setorder(
  dent_cmp,
  -temporal_mean_score,
  -temporal_max_score,
  -n_bins_detected,
  -n_lineages_detected
)
dent_cmp[, temporal_pair_rank := seq_len(.N)]
dent_cmp[, temporal_top10pct := temporal_pair_rank <= max(1L, ceiling(.N * 0.10))]
dent_cmp[, temporal_broad := n_bins_detected >= 2L & n_lineages_detected >= 1L]
dent_cmp[, temporal_candidate := temporal_top10pct & temporal_broad]
dent_cmp[, temporal_only_vs_external :=
  temporal_candidate & !cellchat_supported & !liana_top10pct]
dent_cmp[, benchmark_protocol := CONFIG$protocol]
dent_cmp[, result_revision := CONFIG$revision]
dent_cmp[, external_source := CONFIG$external_source]
fwrite(dent_cmp, file.path(out_cmp, "dentategyrus_llr_vs_external_pair_table.csv"))

dent_temporal_only <- dent_cmp[temporal_only_vs_external == TRUE]
fwrite(
  dent_temporal_only,
  file.path(out_temporal, "dentategyrus_temporal_only_vs_external.csv")
)
dent_family <- dent_temporal_only[, .(
  n_pairs = .N,
  mean_temporal_score = mean(temporal_mean_score),
  max_temporal_score = max(temporal_max_score),
  dataset = "dentategyrus"
), by = family]
fwrite(
  dent_family,
  file.path(out_temporal, "dentategyrus_temporal_only_family_summary.csv")
)
dent_case <- dent_cmp[temporal_candidate == TRUE, .(
  pair_name, family, temporal_mean_score, temporal_max_score,
  n_bins_detected, n_lineages_detected, temporal_pair_rank,
  comparator_membership, cellchat_supported, liana_supported, liana_top10pct,
  benchmark_protocol, result_revision, external_source
)]
fwrite(dent_case, file.path(out_case, "dentategyrus_case_support_pairs.csv"))

dent_peak <- dent_scores[order(-score), .SD[1L], by = pair_key][, .(
  pair_key,
  peak_score = score,
  peak_bin = bin,
  peak_lineage = lineage,
  peak_source = source,
  peak_target = target
)]
dent_top <- merge(dent_cmp, dent_peak, by = "pair_key", all.x = TRUE)
setorder(dent_top, temporal_pair_rank)
dent_top <- head(dent_top, 50L)
dent_top[, dataset := "dentategyrus"]
fwrite(
  dent_top,
  file.path(out_narrative, "dentategyrus_top_temporal_lr_pairs.csv")
)
dent_summary <- data.table(
  dataset = "dentategyrus",
  benchmark_protocol = dent_res$benchmark_protocol,
  result_revision = dent_res$result_revision,
  n_cells = as.integer(dent_res$n_cells),
  n_genes = as.integer(dent_res$n_genes),
  n_types = as.integer(dent_res$n_types),
  n_lineages = as.integer(dent_res$n_lineages),
  n_lr_pairs = as.integer(nrow(dent_res$lr_scores)),
  n_temporal_rows = as.integer(nrow(dent_scores_all)),
  n_temporal_pairs = as.integer(uniqueN(dent_scores$pair_idx)),
  mean_temporal_score = mean(dent_scores$score),
  max_temporal_score = max(dent_scores$score)
)
fwrite(
  dent_summary,
  file.path(out_narrative, "dentategyrus_dataset_summary.csv")
)

gse_scores_all <- as.data.table(gse_res$temporal_lr$scores)
gse_scores <- gse_scores_all[
  is.finite(score) & expression_supported == TRUE
]
gse_scores <- merge(
  gse_scores,
  as.data.table(gse_res$pair_gene_map),
  by = "pair_idx",
  all.x = TRUE
)
gse_scores[, pair_name := paste(ligand_genes, receptor_genes, sep = "_")]
gse_scores[, pair_key := paste0(
  tolower(gsub("[^A-Za-z0-9]+", "", ligand_genes)),
  "___",
  tolower(gsub("[^A-Za-z0-9]+", "", receptor_genes))
)]
gse_top <- gse_scores[, .(
  ligand_genes = ligand_genes[1],
  receptor_genes = receptor_genes[1],
  temporal_mean_score = mean(score),
  temporal_median_score = median(score),
  temporal_max_score = max(score),
  n_bins_detected = uniqueN(bin[score > 0]),
  n_lineages_detected = uniqueN(lineage[score > 0]),
  n_source_target = uniqueN(paste(source, target, sep = "->")),
  peak_score = max(score),
  peak_bin = bin[which.max(score)],
  peak_lineage = lineage[which.max(score)],
  peak_source = source[which.max(score)],
  peak_target = target[which.max(score)]
), by = .(pair_name, pair_key)]
gse_top[, family := fcase(
  grepl("^(Wnt)|_(Fzd)", pair_name, ignore.case = TRUE), "Wnt-Fzd",
  grepl("^(Dll|Jag)|_(Notch)", pair_name, ignore.case = TRUE), "Notch",
  grepl("^(Fgf)|_(Fgfr)", pair_name, ignore.case = TRUE), "FGF",
  grepl("^(Bmp|Gdf|Tgfb)|_(Bmpr|Acvr|Tgfbr)", pair_name, ignore.case = TRUE), "BMP-TGF",
  grepl("^(Ccl|Cxcl)|_(Ackr|Ccr|Cxcr)", pair_name, ignore.case = TRUE), "Chemokine",
  grepl("^(Lama|Lamb|Lamc|Agrn|Hspg2|Col)", ligand_genes, ignore.case = TRUE), "ECM-adhesion",
  default = "other"
)]
setorder(gse_top, -temporal_mean_score, -temporal_max_score)
gse_top[, temporal_pair_rank := seq_len(.N)]
gse_top <- head(gse_top, 50L)
gse_top[, dataset := "GSE210152_athero"]
gse_top[, benchmark_protocol := CONFIG$protocol]
gse_top[, result_revision := CONFIG$revision]
fwrite(
  gse_top,
  file.path(out_narrative, "GSE210152_athero_top_temporal_lr_pairs.csv")
)
gse_summary <- data.table(
  dataset = "GSE210152_athero",
  benchmark_protocol = gse_res$benchmark_protocol,
  result_revision = gse_res$result_revision,
  n_cells = as.integer(gse_res$n_cells),
  n_genes = as.integer(gse_res$n_genes),
  n_types = as.integer(gse_res$n_types),
  n_lineages = as.integer(gse_res$n_lineages),
  n_lr_pairs = as.integer(nrow(gse_res$lr_scores)),
  n_temporal_rows = as.integer(nrow(gse_scores_all)),
  n_temporal_pairs = as.integer(uniqueN(gse_scores$pair_idx)),
  mean_temporal_score = mean(gse_scores$score),
  max_temporal_score = max(gse_scores$score)
)
fwrite(
  gse_summary,
  file.path(out_narrative, "GSE210152_athero_dataset_summary.csv")
)

summary_all <- rbindlist(list(
  data.table(
    dataset = "pancreas",
    llr_pairs_total = nrow(pan_cmp),
    llr_pairs_candidates = sum(pan_cmp$temporal_candidate),
    cellchat_pairs = sum(pan_cmp$cellchat_supported),
    liana_pairs = sum(pan_cmp$liana_supported),
    llr_supported_by_cellchat = sum(
      pan_cmp$temporal_candidate & pan_cmp$cellchat_supported
    ),
    llr_supported_by_liana = sum(
      pan_cmp$temporal_candidate & pan_cmp$liana_supported
    ),
    llr_supported_by_both = sum(
      pan_cmp$temporal_candidate & pan_cmp$cellchat_supported &
        pan_cmp$liana_supported
    ),
    llr_temporal_only_vs_external = sum(pan_cmp$temporal_only_vs_external),
    benchmark_protocol = CONFIG$protocol,
    result_revision = CONFIG$revision,
    external_source = CONFIG$external_source
  ),
  data.table(
    dataset = "dentategyrus",
    llr_pairs_total = nrow(dent_cmp),
    llr_pairs_candidates = sum(dent_cmp$temporal_candidate),
    cellchat_pairs = sum(dent_cmp$cellchat_supported),
    liana_pairs = sum(dent_cmp$liana_supported),
    llr_supported_by_cellchat = sum(
      dent_cmp$temporal_candidate & dent_cmp$cellchat_supported
    ),
    llr_supported_by_liana = sum(
      dent_cmp$temporal_candidate & dent_cmp$liana_supported
    ),
    llr_supported_by_both = sum(
      dent_cmp$temporal_candidate & dent_cmp$cellchat_supported &
        dent_cmp$liana_supported
    ),
    llr_temporal_only_vs_external = sum(dent_cmp$temporal_only_vs_external),
    benchmark_protocol = CONFIG$protocol,
    result_revision = CONFIG$revision,
    external_source = CONFIG$external_source
  )
))
fwrite(summary_all, file.path(out_cmp, "external_cci_comparator_summary.csv"))

family_all <- rbindlist(list(pan_family, dent_family), fill = TRUE)
family_wide <- dcast(family_all, family ~ dataset, value.var = "n_pairs", fill = 0)
family_wide[, datasets_supported := rowSums(.SD > 0), .SDcols = c(
  "pancreas", "dentategyrus"
)]
setorder(family_wide, -datasets_supported, -pancreas, -dentategyrus)
fwrite(
  family_wide,
  file.path(out_temporal, "cross_dataset_temporal_family_recurrence.csv")
)

pair_all <- rbindlist(list(
  pan_temporal_only[, .(
    dataset = "pancreas", pair_name, pair_key, family, temporal_mean_score
  )],
  dent_temporal_only[, .(
    dataset = "dentategyrus", pair_name, pair_key, family, temporal_mean_score
  )]
), fill = TRUE)
pair_recurrence <- pair_all[, .(
  datasets_supported = uniqueN(dataset),
  datasets = paste(sort(unique(dataset)), collapse = ";"),
  family = family[1],
  mean_temporal_score = mean(temporal_mean_score)
), by = .(pair_name, pair_key)]
setorder(pair_recurrence, -datasets_supported, -mean_temporal_score)
fwrite(
  pair_recurrence,
  file.path(out_temporal, "cross_dataset_temporal_pair_recurrence.csv")
)

fwrite(
  rbindlist(list(pan_summary, dent_summary, gse_summary), fill = TRUE),
  file.path(out_narrative, "linealinkr_narrative_dataset_summary.csv")
)
fwrite(
  rbindlist(list(pan_top, dent_top, gse_top), fill = TRUE),
  file.path(out_narrative, "linealinkr_narrative_top_temporal_lr_pairs.csv")
)

cat("EXTERNAL_CCI_TABLES_064_REBUILT\n")
