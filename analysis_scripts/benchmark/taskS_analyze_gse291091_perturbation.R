library(data.table)
library(edgeR)
library(fgsea)
library(msigdbr)

set.seed(12345)
Sys.setenv(OMP_NUM_THREADS = parallel::detectCores())
data.table::setDTthreads(parallel::detectCores())

CONFIG <- list(
  paths = list(
    counts = paste0(
      "F:/biotrade/0.64LLR/results/biology/muscle/",
      "gse291091_pseudobulk_counts.csv"
    ),
    metadata = paste0(
      "F:/biotrade/0.64LLR/results/biology/muscle/",
      "gse291091_pseudobulk_metadata.csv"
    ),
    output = "F:/biotrade/0.64LLR/results/biology/muscle"
  ),
  threads = parallel::detectCores(),
  seed = 12345L,
  candidate_genes = c(
    "C3", "Fgf7", "Il6", "Il33", "Fbn1", "Tgfb1", "Csf1", "Vegfa",
    "Col1a1", "Col1a2", "Fn1", "Thbs2"
  ),
  focus_pathways = c(
    "HALLMARK_IL6_JAK_STAT3_SIGNALING",
    "HALLMARK_INFLAMMATORY_RESPONSE",
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
    "HALLMARK_INTERFERON_GAMMA_RESPONSE",
    "HALLMARK_TGF_BETA_SIGNALING",
    "HALLMARK_ANGIOGENESIS",
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"
  )
)

counts_dt <- fread(CONFIG$paths$counts, encoding = "UTF-8")
metadata <- fread(CONFIG$paths$metadata, encoding = "UTF-8")
genes <- counts_dt$gene
count_matrix <- as.matrix(counts_dt[, -"gene"])
rownames(count_matrix) <- genes

mac_meta <- metadata[
  context == "Mononuclear phagocytes" & injury_status == "3dpi"
][order(factor(fap_status, levels = c("control", "depleted")), replicate)]
mac_group <- factor(mac_meta$fap_status, levels = c("control", "depleted"))
mac_y <- DGEList(
  counts = count_matrix[, mac_meta$group_id, drop = FALSE],
  group = mac_group
)
mac_keep <- filterByExpr(mac_y, group = mac_group)
mac_y <- calcNormFactors(mac_y[mac_keep, , keep.lib.sizes = FALSE])
mac_design <- model.matrix(~0 + mac_group)
colnames(mac_design) <- levels(mac_group)
mac_y <- estimateDisp(mac_y, mac_design, robust = TRUE)
mac_fit <- glmQLFit(mac_y, mac_design, robust = TRUE)
mac_test <- glmQLFTest(mac_fit, contrast = c(-1, 1))
mac_de <- as.data.table(topTags(mac_test, n = Inf, sort.by = "none")$table, keep.rownames = "gene")
mac_de[, `:=`(
  context = "Mononuclear phagocytes",
  contrast = "FAP-depleted minus control at 3 dpi"
)]

fap_3dpi_meta <- metadata[
  context == "FAP compartment" & injury_status == "3dpi"
][order(factor(fap_status, levels = c("control", "depleted")), replicate)]
fap_3dpi_group <- factor(
  fap_3dpi_meta$fap_status,
  levels = c("control", "depleted")
)
fap_3dpi_y <- DGEList(
  counts = count_matrix[, fap_3dpi_meta$group_id, drop = FALSE],
  group = fap_3dpi_group
)
fap_3dpi_keep <- filterByExpr(fap_3dpi_y, group = fap_3dpi_group)
fap_3dpi_y <- calcNormFactors(
  fap_3dpi_y[fap_3dpi_keep, , keep.lib.sizes = FALSE]
)
fap_3dpi_design <- model.matrix(~0 + fap_3dpi_group)
colnames(fap_3dpi_design) <- levels(fap_3dpi_group)
fap_3dpi_y <- estimateDisp(fap_3dpi_y, fap_3dpi_design, robust = TRUE)
fap_3dpi_fit <- glmQLFit(fap_3dpi_y, fap_3dpi_design, robust = TRUE)
fap_3dpi_test <- glmQLFTest(fap_3dpi_fit, contrast = c(-1, 1))
fap_3dpi_de <- as.data.table(
  topTags(fap_3dpi_test, n = Inf, sort.by = "none")$table,
  keep.rownames = "gene"
)
fap_3dpi_de[, `:=`(
  context = "FAP compartment",
  contrast = "FAP-depleted minus control at 3 dpi"
)]

fap_injury_meta <- metadata[
  context == "FAP compartment" & fap_status == "control"
][order(factor(injury_status, levels = c("uninjured", "3dpi")), replicate)]
fap_injury_group <- factor(
  fap_injury_meta$injury_status,
  levels = c("uninjured", "3dpi")
)
fap_injury_y <- DGEList(
  counts = count_matrix[, fap_injury_meta$group_id, drop = FALSE],
  group = fap_injury_group
)
fap_injury_keep <- filterByExpr(fap_injury_y, group = fap_injury_group)
fap_injury_y <- calcNormFactors(
  fap_injury_y[fap_injury_keep, , keep.lib.sizes = FALSE]
)
fap_injury_design <- model.matrix(~0 + fap_injury_group)
colnames(fap_injury_design) <- levels(fap_injury_group)
fap_injury_y <- estimateDisp(fap_injury_y, fap_injury_design, robust = TRUE)
fap_injury_fit <- glmQLFit(fap_injury_y, fap_injury_design, robust = TRUE)
fap_injury_test <- glmQLFTest(fap_injury_fit, contrast = c(-1, 1))
fap_injury_de <- as.data.table(
  topTags(fap_injury_test, n = Inf, sort.by = "none")$table,
  keep.rownames = "gene"
)
fap_injury_de[, `:=`(
  context = "FAP compartment",
  contrast = "3 dpi minus uninjured in controls"
)]

hallmark <- as.data.table(msigdbr(species = "Mus musculus", collection = "H"))
hallmark_pathways <- split(hallmark$gene_symbol, hallmark$gs_name)
mac_ranks <- setNames(mac_de$logFC, mac_de$gene)
mac_ranks <- sort(mac_ranks[is.finite(mac_ranks)], decreasing = TRUE)
mac_fgsea <- as.data.table(fgseaMultilevel(
  pathways = hallmark_pathways,
  stats = mac_ranks,
  minSize = 10L,
  maxSize = 500L,
  eps = 0
))
mac_fgsea[, `:=`(
  context = "Mononuclear phagocytes",
  contrast = "FAP-depleted minus control at 3 dpi"
)]

fap_injury_ranks <- setNames(fap_injury_de$logFC, fap_injury_de$gene)
fap_injury_ranks <- sort(
  fap_injury_ranks[is.finite(fap_injury_ranks)],
  decreasing = TRUE
)
fap_injury_fgsea <- as.data.table(fgseaMultilevel(
  pathways = hallmark_pathways,
  stats = fap_injury_ranks,
  minSize = 10L,
  maxSize = 500L,
  eps = 0
))
fap_injury_fgsea[, `:=`(
  context = "FAP compartment",
  contrast = "3 dpi minus uninjured in controls"
)]

all_de <- rbindlist(list(mac_de, fap_3dpi_de, fap_injury_de), fill = TRUE)
all_de[, `:=`(
  analysis_unit = "sample pseudobulk",
  n_biological_replicates_per_group = 2L,
  inference_scope = "sample-level pseudobulk perturbation contrast"
)]
candidate_effects <- all_de[gene %in% CONFIG$candidate_genes][
  order(contrast, context, match(gene, CONFIG$candidate_genes))
]
pathway_results <- rbindlist(list(mac_fgsea, fap_injury_fgsea), fill = TRUE)
pathway_results[, `:=`(
  analysis_unit = "sample pseudobulk",
  n_biological_replicates_per_group = 2L,
  inference_scope = "sample-level pseudobulk perturbation contrast"
)]
focus_pathways <- pathway_results[pathway %in% CONFIG$focus_pathways][
  order(contrast, context, pathway)
]

de_summary <- all_de[
  contrast == "FAP-depleted minus control at 3 dpi",
  .(
    evidence_type = "differential-expression summary",
    feature = "all tested genes",
    effect = NA_real_,
    p_value = NA_real_,
    adjusted_p = NA_real_,
    n_tested = .N,
    n_fdr_005 = sum(FDR < 0.05, na.rm = TRUE),
    n_up_fdr_005 = sum(FDR < 0.05 & logFC > 0, na.rm = TRUE),
    n_down_fdr_005 = sum(FDR < 0.05 & logFC < 0, na.rm = TRUE)
  ),
  by = .(context, contrast)
]
candidate_summary <- candidate_effects[
  contrast == "FAP-depleted minus control at 3 dpi",
  .(
    evidence_type = "prespecified candidate gene",
    feature = gene,
    effect = logFC,
    p_value = PValue,
    adjusted_p = FDR,
    n_tested = NA_integer_,
    n_fdr_005 = NA_integer_,
    n_up_fdr_005 = NA_integer_,
    n_down_fdr_005 = NA_integer_
  ),
  by = .(context, contrast)
]
pathway_summary <- focus_pathways[
  contrast == "FAP-depleted minus control at 3 dpi",
  .(
    evidence_type = "prespecified Hallmark pathway",
    feature = pathway,
    effect = NES,
    p_value = pval,
    adjusted_p = padj,
    n_tested = NA_integer_,
    n_fdr_005 = NA_integer_,
    n_up_fdr_005 = NA_integer_,
    n_down_fdr_005 = NA_integer_
  ),
  by = .(context, contrast)
]
fap_depletion_summary <- rbindlist(
  list(de_summary, candidate_summary, pathway_summary),
  use.names = TRUE,
  fill = TRUE
)
fap_depletion_summary[, `:=`(
  analysis_unit = "sample pseudobulk",
  n_biological_replicates_per_group = 2L,
  evidence_scope = "sample-level pseudobulk perturbation contrast"
)]

fap_all_meta <- metadata[context == "FAP compartment"][order(group_id)]
fap_all_y <- calcNormFactors(DGEList(
  counts = count_matrix[, fap_all_meta$group_id, drop = FALSE]
))
fap_logcpm <- as.data.table(
  cpm(fap_all_y, log = TRUE, prior.count = 2),
  keep.rownames = "gene"
)
fap_logcpm <- melt(
  fap_logcpm[gene %in% CONFIG$candidate_genes],
  id.vars = "gene",
  variable.name = "group_id",
  value.name = "logCPM"
)
fap_logcpm <- merge(fap_logcpm, fap_all_meta, by = "group_id", all.x = TRUE)[
  order(match(gene, CONFIG$candidate_genes), injury_status, fap_status, replicate)
]

fwrite(
  all_de,
  file.path(CONFIG$paths$output, "gse291091_pseudobulk_differential_expression.csv")
)
fwrite(
  all_de,
  file.path(CONFIG$paths$output, "fap_depletion_pseudobulk_064.csv")
)
fwrite(
  fap_depletion_summary,
  file.path(CONFIG$paths$output, "fap_depletion_summary_064.csv")
)
fwrite(
  candidate_effects,
  file.path(CONFIG$paths$output, "gse291091_candidate_ligand_effects.csv")
)
fwrite(
  pathway_results,
  file.path(CONFIG$paths$output, "gse291091_hallmark_enrichment.csv")
)
fwrite(
  focus_pathways,
  file.path(CONFIG$paths$output, "gse291091_focus_pathway_enrichment.csv")
)
fwrite(
  fap_logcpm,
  file.path(CONFIG$paths$output, "gse291091_candidate_ligand_logcpm.csv")
)

print(candidate_effects)
print(focus_pathways)
