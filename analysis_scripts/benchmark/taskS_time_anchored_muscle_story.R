.llr_cmd <- commandArgs(trailingOnly = FALSE)
.llr_file <- normalizePath(sub("^--file=", "", grep("^--file=", .llr_cmd, value = TRUE)[1L]), winslash = "/", mustWork = TRUE)
.llr_pipeline <- normalizePath(file.path(dirname(.llr_file), ".."), winslash = "/", mustWork = TRUE)
source(file.path(.llr_pipeline, "config", "llr_bootstrap.R"), local = TRUE, encoding = "UTF-8")
.llr_assets <- llr_assets_dir()
.llr_config_env <- new.env(parent = globalenv())
sys.source(file.path(.llr_assets, "pipeline", "config", "00_config.R"), envir = .llr_config_env)
PROJECT_CONFIG <- .llr_config_env$CONFIG
rm(.llr_config_env)
.libPaths(c(file.path(.llr_assets, "rlib"), .libPaths()))

library(data.table)
library(Matrix)
library(LineaLinkR)

set.seed(12345)
Sys.setenv(OMP_NUM_THREADS = parallel::detectCores())
data.table::setDTthreads(parallel::detectCores())

CONFIG <- list(
  paths = list(
    result = file.path(.llr_assets, "cache/linealinkr_workspace/GSE143437_muscle_regeneration/linealinkr_complete.rds"),
    matrix = file.path(.llr_assets, "cache/linealinkr_workspace/GSE143437_muscle_regeneration/extract/GSE143437_muscle_regeneration_matrix.mtx"),
    metadata = file.path(.llr_assets, "cache/linealinkr_workspace/GSE143437_muscle_regeneration/extract/GSE143437_muscle_regeneration_meta.csv"),
    genes = file.path(.llr_assets, "cache/linealinkr_workspace/GSE143437_muscle_regeneration/extract/GSE143437_muscle_regeneration_genes.txt"),
    depletion = file.path(.llr_assets, "results/biology/muscle/fap_depletion_summary_064.csv"),
    output = file.path(.llr_assets, "results/biology/muscle"),
    cache = file.path(.llr_assets, "cache/linealinkr_workspace/GSE143437_muscle_regeneration/time_anchored_muscle_story_064.rds")
  ),
  threads = parallel::detectCores(),
  seed = 12345L,
  package_version = PROJECT_CONFIG$llr$package_version,
  benchmark_protocol = PROJECT_CONFIG$llr$benchmark_protocol,
  result_revision = PROJECT_CONFIG$llr$result_revision,
  n_bins = 4L,
  n_perms = 199L,
  phase_days = c(0, 2, 5, 7),
  candidate_ligands = c(
    "C3", "Fgf7", "Il6", "Il33", "Tgfb1", "Csf1", "Vegfa",
    "Col1a1", "Col1a2", "Fn1", "Thbs2"
  ),
  target_labels = c(
    "MuSCs and progenitors",
    "Monocytes/Macrophages/Platelets",
    "Pro-inflammatory macrophages",
    "Anti-inflammatory macrophages",
    "Resident Macrophages/APCs"
  )
)

dir.create(CONFIG$paths$output, recursive = TRUE, showWarnings = FALSE)

result <- readRDS(CONFIG$paths$result)
stopifnot(identical(as.character(packageVersion("LineaLinkR")), CONFIG$package_version))
stopifnot(identical(as.character(result$package_version), CONFIG$package_version))
stopifnot(identical(
  as.character(result$result_revision),
  CONFIG$result_revision
))
stopifnot(identical(
  as.character(result$benchmark_protocol),
  CONFIG$benchmark_protocol
))
metadata <- fread(CONFIG$paths$metadata, encoding = "UTF-8")
genes <- readLines(CONFIG$paths$genes, encoding = "UTF-8")
expr <- as(readMM(CONFIG$paths$matrix), "CsparseMatrix")
rownames(expr) <- genes

type_levels <- sort(unique(as.character(metadata$clusters)))
celltypes <- as.integer(factor(metadata$clusters, levels = type_levels))
attr(celltypes, "levels") <- seq_along(type_levels)
type_map <- data.table(type_idx = seq_along(type_levels), label = type_levels)

pair_map <- as.data.table(result$pair_gene_map)[order(pair_idx)]
lr_db <- list(
  ligand_rows = lapply(
    strsplit(pair_map$ligand_genes, "_", fixed = TRUE),
    match,
    table = genes
  ),
  receptor_rows = lapply(
    strsplit(pair_map$receptor_genes, "_", fixed = TRUE),
    match,
    table = genes
  )
)

expr_lr <- expr %*% Diagonal(x = 10000 / colSums(expr))
expr_lr@x <- log1p(expr_lr@x)
observed_time <- metadata$time_days / max(CONFIG$phase_days)
lineage_weights <- matrix(1, nrow = ncol(expr_lr), ncol = 1L)

source_idx <- type_map[label == "FAPs", type_idx]
target_idx <- type_map[label %in% CONFIG$target_labels, type_idx]
requested_contexts <- CJ(
  pair_idx = pair_map$pair_idx,
  source = source_idx,
  target = target_idx,
  unique = TRUE
)

temporal <- ComputeTemporalLRWithIntermediate(
  expr_lr,
  celltypes,
  lr_db,
  observed_time,
  lineage_weights,
  n_bins = CONFIG$n_bins,
  requested_contexts = requested_contexts
)

scores <- as.data.table(temporal$scores)
scores[, source_label := type_map$label[source]]
scores[, target_label := type_map$label[target]]
scores[, phase_day := CONFIG$phase_days[bin]]
scores <- merge(scores, pair_map, by = "pair_idx", all.x = TRUE)
scores[, valid_score := fifelse(
  structural_valid & expression_supported & is.finite(score),
  score,
  NA_real_
)]

static_scores <- as.data.table(result$lr_scores)
static_scores[, source_label := type_map$label[source]]
static_scores[, target_label := type_map$label[target]]
static_scores <- merge(static_scores, pair_map, by = "pair_idx", all.x = TRUE)
static_scores <- static_scores[
  source_label == "FAPs" & target_label %in% CONFIG$target_labels,
  .(pair_idx, target_label, static_score = score, static_pvalue = pvalue)
]

candidate_scores_all <- scores[ligand_genes %in% CONFIG$candidate_ligands]
candidate_scores <- candidate_scores_all[is.finite(valid_score)]
candidate_peak <- candidate_scores[
  , .SD[which.max(valid_score)],
  by = .(pair_idx, ligand_genes, receptor_genes, target_label)
][, .(
  pair_idx,
  ligand_genes,
  receptor_genes,
  target_label,
  peak_day = phase_day,
  peak_score = valid_score
)]
candidate_n <- candidate_scores[
  , .(n_valid_phases = sum(is.finite(valid_score))),
  by = .(pair_idx, ligand_genes, receptor_genes, target_label)
]
candidate_summary <- merge(
  candidate_peak,
  candidate_n,
  by = c("pair_idx", "ligand_genes", "receptor_genes", "target_label"),
  all = TRUE
)
candidate_summary <- merge(
  candidate_summary,
  static_scores,
  by = c("pair_idx", "target_label"),
  all.x = TRUE
)

candidate_trajectories <- unique(candidate_scores[
  , .(pair_idx, source = source_idx, target, lineage = 1L)
])
temporal_null <- ComputeTemporalPermutation(
  expr_lr,
  celltypes,
  lr_db,
  observed_time,
  lineage_weights,
  n_bins = CONFIG$n_bins,
  n_perms = CONFIG$n_perms,
  seed = CONFIG$seed,
  stat = "amplitude",
  min_perm_observed = 100L,
  requested_trajectories = candidate_trajectories,
  permutation_stratum = "block_celltype"
)

hill_scale <- temporal_null$hill_scale
c3_cell_expression <- data.table(
  sampleID = metadata$sampleID,
  time_days = metadata$time_days,
  cell_type = metadata$clusters,
  C3 = as.numeric(expr_lr["C3", ]),
  C3ar1 = as.numeric(expr_lr["C3ar1", ])
)
c3_source_sample <- c3_cell_expression[
  cell_type == "FAPs",
  .(
    ligand_mean = mean(C3),
    ligand_pct = mean(C3 > 0),
    n_source_cells = .N
  ),
  by = .(sampleID, time_days)
]
c3_target_sample <- c3_cell_expression[
  cell_type == "Resident Macrophages/APCs",
  .(
    receptor_mean = mean(C3ar1),
    receptor_pct = mean(C3ar1 > 0),
    n_target_cells = .N
  ),
  by = .(sampleID, time_days)
]
c3_sample_scores <- merge(
  c3_source_sample,
  c3_target_sample,
  by = c("sampleID", "time_days"),
  all = FALSE
)
c3_sample_scores[, geometric_mean := sqrt(ligand_mean * receptor_mean)]
c3_sample_scores[, sample_lr_score := geometric_mean / (hill_scale + geometric_mean)]
c3_sample_scores[, `:=`(
  hill_scale = hill_scale,
  source = "FAPs",
  target = "Resident Macrophages/APCs",
  pair = "C3-C3ar1",
  package_version = CONFIG$package_version,
  benchmark_protocol = CONFIG$benchmark_protocol,
  result_revision = CONFIG$result_revision
)]
c3_sample_scores <- c3_sample_scores[order(time_days, sampleID)]
c3_sample_day_summary <- c3_sample_scores[
  , .(
    mean_sample_lr_score = mean(sample_lr_score),
    sd_sample_lr_score = sd(sample_lr_score),
    min_sample_lr_score = min(sample_lr_score),
    max_sample_lr_score = max(sample_lr_score),
    n_samples = .N
  ),
  by = time_days
][order(time_days)]
c3_sample_anova <- summary(aov(sample_lr_score ~ factor(time_days), data = c3_sample_scores))[[1]][
  "factor(time_days)", "Pr(>F)"
]
c3_sample_kruskal <- kruskal.test(
  sample_lr_score ~ factor(time_days),
  data = c3_sample_scores
)$p.value
c3_sample_inference <- data.table(
  endpoint = "FAP C3-C3ar1 to resident macrophages",
  n_samples = nrow(c3_sample_scores),
  n_timepoints = uniqueN(c3_sample_scores$time_days),
  anova_p = c3_sample_anova,
  kruskal_p = c3_sample_kruskal,
  evidence_scope = "sample-level expression score across recorded injury days"
)
depletion_summary <- fread(CONFIG$paths$depletion, encoding = "UTF-8")
depletion_c3 <- depletion_summary[
  evidence_type == "prespecified candidate gene" & feature == "C3",
  .(
    dataset = "GSE291091",
    evidence_axis = "whole-population FAP depletion",
    context,
    effect = effect,
    p_value,
    adjusted_p,
    inference_unit = analysis_unit,
    evidence_scope = "whole-population FAP-depletion perturbation contrast"
  )
]
cross_dataset_evidence <- rbindlist(
  list(
    data.table(
      dataset = "GSE143437",
      evidence_axis = "recorded injury-day C3-C3ar1 expression score",
      context = "FAPs to resident macrophages/APCs",
      effect = NA_real_,
      p_value = c3_sample_inference$kruskal_p,
      adjusted_p = NA_real_,
      inference_unit = "biological sample",
      evidence_scope = "recorded injury-day expression score"
    ),
    depletion_c3
  ),
  use.names = TRUE,
  fill = TRUE
)
permutation_summary <- as.data.table(temporal_null$trajectories)
permutation_summary[, target_label := type_map$label[target]]
permutation_summary <- merge(
  permutation_summary,
  pair_map,
  by = "pair_idx",
  all.x = TRUE
)[
  , .(
    pair_idx,
    ligand_genes,
    receptor_genes,
    target_label,
    amplitude = temporal_stat,
    p_temporal,
    q_temporal,
    n_perm_observed,
    perm_coverage
  )
]
candidate_summary <- merge(
  candidate_summary,
  permutation_summary,
  by = c("pair_idx", "ligand_genes", "receptor_genes", "target_label"),
  all.x = TRUE
)

missing_fbn1 <- data.table(
  pair_idx = NA_integer_,
  ligand_genes = "Fbn1",
  receptor_genes = NA_character_,
  target_label = CONFIG$target_labels,
  peak_day = NA_real_,
  peak_score = NA_real_,
  n_valid_phases = 0L,
  static_score = NA_real_,
  static_pvalue = NA_real_,
  amplitude = NA_real_,
  p_temporal = NA_real_,
  q_temporal = NA_real_,
  n_perm_observed = 0L,
  perm_coverage = NA_real_
)
candidate_summary <- rbindlist(
  list(candidate_summary, missing_fbn1),
  use.names = TRUE,
  fill = TRUE
)[order(target_label, peak_day, -peak_score)]

candidate_db_coverage <- CJ(
  ligand_genes = c(CONFIG$candidate_ligands, "Fbn1"),
  target_label = CONFIG$target_labels,
  unique = TRUE
)
candidate_db_counts <- pair_map[
  ligand_genes %in% CONFIG$candidate_ligands,
  .(n_lr_database_pairs = uniqueN(pair_idx)),
  by = ligand_genes
]
candidate_evaluable_counts <- candidate_scores[
  , .(n_evaluable_pairs = uniqueN(pair_idx)),
  by = .(ligand_genes, target_label)
]
candidate_db_coverage <- merge(
  candidate_db_coverage,
  candidate_db_counts,
  by = "ligand_genes",
  all.x = TRUE
)
candidate_db_coverage <- merge(
  candidate_db_coverage,
  candidate_evaluable_counts,
  by = c("ligand_genes", "target_label"),
  all.x = TRUE
)
candidate_db_coverage[is.na(n_lr_database_pairs), n_lr_database_pairs := 0L]
candidate_db_coverage[is.na(n_evaluable_pairs), n_evaluable_pairs := 0L]
candidate_db_coverage[, status := fifelse(
  n_lr_database_pairs == 0L,
  "not represented in dataset-specific scored LR map",
  fifelse(n_evaluable_pairs == 0L, "represented but not evaluable", "evaluable")
)][order(target_label, match(ligand_genes, c(CONFIG$candidate_ligands, "Fbn1")))]

top_fap_programs <- scores[
  is.finite(valid_score),
  .(
    peak_score = max(valid_score),
    peak_day = phase_day[which.max(valid_score)],
    mean_score = mean(valid_score),
    n_valid_phases = .N
  ),
  by = .(pair_idx, ligand_genes, receptor_genes, target_label)
][order(target_label, -peak_score)][, head(.SD, 20L), by = target_label]

metadata[, blind_pseudotime := result$pseudotime]
blind_signed_spearman <- cor(
  metadata$blind_pseudotime,
  metadata$time_days,
  method = "spearman",
  use = "complete.obs"
)
blind_ordering_audit <- data.table(
  dataset = "GSE143437",
  n_cells = nrow(metadata),
  direction_agnostic_abs_spearman = abs(blind_signed_spearman),
  cell_spearman = blind_signed_spearman,
  blind_root_day = metadata$time_days[result$start_cell],
  blind_root_cell_type = metadata$clusters[result$start_cell],
  evaluation_axis = "recorded injury day withheld from blind inference",
  application_axis = "observed injury day"
)
sample_ordering <- metadata[
  , .(
    n_cells = .N,
    median_blind_pseudotime = median(blind_pseudotime),
    mean_blind_pseudotime = mean(blind_pseudotime)
  ),
  by = .(sampleID, time_days)
][order(time_days, sampleID)]
phase_cell_counts <- metadata[
  , .(n_cells = .N),
  by = .(time_days, sampleID, cell_type = clusters)
][order(time_days, sampleID, cell_type)]

fwrite(
  candidate_summary,
  file.path(CONFIG$paths$output, "gse143437_time_anchored_candidate_pairs.csv")
)
fwrite(
  candidate_scores_all,
  file.path(CONFIG$paths$output, "gse143437_time_anchored_candidate_trajectories.csv")
)
fwrite(
  candidate_db_coverage,
  file.path(CONFIG$paths$output, "gse143437_candidate_database_coverage.csv")
)
fwrite(
  top_fap_programs,
  file.path(CONFIG$paths$output, "gse143437_top_fap_programs.csv")
)
fwrite(
  blind_ordering_audit,
  file.path(CONFIG$paths$output, "gse143437_blind_ordering_audit.csv")
)
fwrite(
  sample_ordering,
  file.path(CONFIG$paths$output, "gse143437_sample_ordering.csv")
)
fwrite(
  phase_cell_counts,
  file.path(CONFIG$paths$output, "gse143437_phase_cell_counts.csv")
)
fwrite(
  c3_sample_scores,
  file.path(CONFIG$paths$output, "gse143437_sample_level_c3_c3ar1.csv")
)
fwrite(
  c3_sample_scores,
  file.path(CONFIG$paths$output, "muscle_c3_sample_scores_064.csv")
)
fwrite(
  c3_sample_day_summary,
  file.path(CONFIG$paths$output, "gse143437_sample_level_c3_c3ar1_by_day.csv")
)
fwrite(
  c3_sample_inference,
  file.path(CONFIG$paths$output, "gse143437_sample_level_c3_c3ar1_inference.csv")
)
fwrite(
  cross_dataset_evidence,
  file.path(CONFIG$paths$output, "muscle_cross_dataset_evidence_064.csv")
)
saveRDS(
  list(
    scores = scores,
    temporal = temporal,
    temporal_null = temporal_null,
    candidate_summary = candidate_summary,
    candidate_scores = candidate_scores_all,
    candidate_db_coverage = candidate_db_coverage,
    top_fap_programs = top_fap_programs,
    blind_ordering_audit = blind_ordering_audit,
    phase_cell_counts = phase_cell_counts,
    c3_sample_scores = c3_sample_scores,
    c3_sample_day_summary = c3_sample_day_summary,
    c3_sample_inference = c3_sample_inference,
    phase_days = CONFIG$phase_days,
    protocol = "native LineaLinkR temporal LR with observed sampling time in [0,1], one lineage, four bins"
  ),
  CONFIG$paths$cache
)

print(blind_ordering_audit)
print(candidate_summary)
