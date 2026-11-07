# End-to-end lrsim evaluation with independent target x variant denominators.
# No truth target is created from a returned prediction. Wrong forks remain in
# the candidate list as false candidates; no-call targets remain misses.
#
# Evaluation objects are biological sets, not local indices:
#   - LR pairs join on the gene identity in each result's own pair_gene_map,
#     never on pair_idx or file row order.
#   - A prediction is one (sender x receiver x focal side x sibling side)
#     contrast. For L0 runs the focal side is the focal terminal lineage
#     (rows where child != the unit's focal are duplicate emissions and are
#     dropped). For L1 runs the focal side is the child node's descendant
#     terminal set from lineage_paths; the sibling side uses the recorded
#     contributing children when available. For baseline/B runs the focal side is
#     `lineage` and the sibling side is sibling_map at the same fork_cluster.
#   - A truth target is hit when focal set == {truth focal_branch} and the
#     truth other_branch is covered by the sibling set; wrong event, mirror,
#     and non-evaluable rows stay in the candidate list as misses.
set.seed(12345)
Sys.setenv(OMP_NUM_THREADS = as.integer(Sys.getenv("LLR_THREADS", as.character(parallel::detectCores()))))
data.table::setDTthreads(as.integer(Sys.getenv("LLR_THREADS", as.character(parallel::detectCores()))))

REV_ROOT <- normalizePath(Sys.getenv("LLR_WORK_ROOT",getwd()),winslash="/",mustWork=TRUE)
.libPaths(c(file.path(REV_ROOT, "library"), .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(pROC)
  library(PRROC)
})

EXTRACT_ROOT <- Sys.getenv("LLR_SIM_EXTRACT_ROOT",file.path(REV_ROOT,"cache/final065"))
READONLY_REVISION_ROOT <- file.path(REV_ROOT,"inputs/revision")
CFG <- list(
  sim_dir = file.path(READONLY_REVISION_ROOT, "sim"),
  ws = NULL,
  out_dir = file.path(REV_ROOT, "results", "evaluation", "final065"),
  datasets = c("lrsim_shared_trend", "lrsim_focal_branch",
               "lrsim_donor_confounded", "lrsim_dropout", "lrsim_nested_fork",
               "lrsim_sampling_thinning", "lrsim_transient_focal",
               "lrsim_donor_transient_focal", "lrsim_donor_ligand_only",
               "lrsim_donor_receptor_only", "lrsim_donor_nested_fork")
)
SCORE_OVERRIDE_FILE <- Sys.getenv("LLR_EVAL_SCORE_OVERRIDE", "")
SCORE_OVERRIDES <- switch(as.character(nzchar(SCORE_OVERRIDE_FILE)),
  "TRUE" = fread(SCORE_OVERRIDE_FILE), "FALSE" = NULL)
variants <- fread(file.path(REV_ROOT, "experiments/config/variants.csv"))
CFG$ws <- setNames(file.path(REV_ROOT, "cache", variants$workspace), variants$variant)
eval_variant <- Sys.getenv("LLR_EVAL_VARIANT", "Final065")
eval_workspace <- Sys.getenv("LLR_EVAL_WORKSPACE", file.path(REV_ROOT, "cache", "final065"))
switch(as.character(nzchar(eval_variant) && nzchar(eval_workspace)),
  "TRUE" = CFG$ws[eval_variant] <- eval_workspace,
  "FALSE" = NULL)
CFG$datasets <- switch(as.character(nzchar(eval_variant)),
  "TRUE" = trimws(strsplit(Sys.getenv("LLR_EVAL_DATASET", paste(CFG$datasets, collapse = ",")), ",", fixed = TRUE)[[1L]]),
  "FALSE" = CFG$datasets)
CFG$out_dir <- switch(as.character(nzchar(eval_variant)),
  "TRUE" = Sys.getenv("LLR_EVAL_OUT_DIR", CFG$out_dir), "FALSE" = CFG$out_dir)

dir.create(CFG$out_dir, recursive = TRUE, showWarnings = FALSE)

SIMPLE_BASELINES <- c("simple_lr", "simple_lr_smooth",
  "simple_lr_L1D0", "simple_lr_L1D0_smooth")
selected_versions <- switch(as.character(nzchar(eval_variant)),
  "TRUE" = eval_variant, "FALSE" = c(names(CFG$ws), SIMPLE_BASELINES))
jobs <- CJ(dataset = CFG$datasets,
           version = selected_versions)
jobs[, result_status := "not_read"]
baseline_path <- file.path(CFG$out_dir, "baseline_simple", "simple_lr_contrasts.csv")
baseline_path_L1D0 <- file.path(CFG$out_dir, "baseline_simple_L1D0",
  "simple_lr_contrasts.csv")

# Independent truth target base. The key is biological, not a local pair_idx:
# ligand, receptor, sender state, receiver state, focal branch and the sibling
# branch that must sit on the other side of the same event.
truth_all <- rbindlist(lapply(CFG$datasets, function(ds) {
  t <- fread(file.path(CFG$sim_dir, sprintf("%s_truth.csv", ds)))
  t[, .(dataset = ds,
        ligand = as.character(ligand),
        receptor = as.character(receptor),
        role = as.character(role),
        true_positive = true_positive %in% TRUE,
        donor_confounded = donor_confounded %in% TRUE,
        truth_focal_branch = as.character(focal_branch),
        truth_other_branch = as.character(other_branch))]
}), fill = TRUE)
truth_base <- unique(truth_all[true_positive %in% TRUE, .(
  dataset, ligand, receptor,
  source_label = "SRC", target_label = "TGT",
  focal_branch = truth_focal_branch, other_branch = truth_other_branch)])
truth_base[, target_key := paste(dataset, ligand, receptor, source_label,
                                 target_label, focal_branch, other_branch,
                                 sep = "|")]
truth_counts <- truth_base[, .(n_targets = uniqueN(target_key)), by = dataset]
variant_grid <- merge(
  CJ(dataset = CFG$datasets, version = selected_versions),
  truth_counts, by = "dataset", all.x = TRUE
)
variant_grid[is.na(n_targets), n_targets := 0L]

# Per-dataset structural maps from each run's own saved objects. Terminal leaf
# branches come from the cell metadata majority label of each leaf cluster,
# identical across LRs because scores are never consulted.
BuildRunMaps <- function(x, meta) {
  leaf_cl <- vapply(x$lineage_paths, function(p) tail(p, 1L), numeric(1L))
  lc <- x$lineage_clusters
  tb <- as.character(meta$true_branch)
  leaf_branch <- vapply(leaf_cl, function(cl) {
    vals <- tb[lc == cl]
    names(sort(table(vals), decreasing = TRUE))[1L]
  }, character(1L))

  paths <- x$lineage_paths
  child_desc <- function(ev, ch) {
    which(vapply(paths, function(p) {
      pos <- match(ev, p)
      !is.na(pos) && pos < length(p) && p[pos + 1L] == ch
    }, logical(1L)))
  }

  smap <- as.data.table(x$sibling_map)
  foc_order <- sort(unique(smap$focal_lineage))
  sibs_by_foc <- lapply(foc_order, function(f) {
    sort(unique(smap$sibling_lineage[smap$focal_lineage == f]))
  })
  names(sibs_by_foc) <- as.character(foc_order)
  conflict_focs <- smap[, uniqueN(fork_cluster), by = focal_lineage][
    V1 > 1L, focal_lineage]
  smap_siblings <- tapply(
    seq_len(nrow(smap)),
    paste(smap$focal_lineage, smap$fork_cluster),
    function(ii) smap$sibling_lineage[ii]
  )

  list(leaf_branch = leaf_branch, child_desc = child_desc,
       sibs_by_foc = sibs_by_foc, foc_order = foc_order,
       conflict_focs = conflict_focs, smap_siblings = smap_siblings)
}

set_label <- function(idx, leaf_branch) {
  if (length(idx) == 0L) return(NA_character_)
  paste(sort(unique(leaf_branch[idx])), collapse = "+")
}

# Attach focal_set / other_set for every returned contrast row.
# L0 ("nearest"): one unit per focal lineage; the unit emits rows for every
#   child, only the row with child == focal is the real object.
# L1 ("all"): one unit per fork event; every emitted child is a real object
#   whose identity is the descendant terminal set of (event -> child).
AttachSets <- function(fc, maps, ver) {
  fc[, `:=`(focal_set = NA_character_, other_set = NA_character_,
            duplicate_emission = FALSE, event_conflict = FALSE)]
  is_event_run <- all(c("event_node", "child", "unit_id") %in% names(fc))
  if (is_event_run && !isTRUE(all(fc$lineage == fc$child))) {
    # L1: child is an internal tree node id.
    ev_children <- lapply(split(fc$child, fc$unit_id), unique)
    desc_cache <- new.env(parent = emptyenv())
    desc_of <- function(ev, ch) {
      k <- paste(ev, ch, sep = "_")
      v <- desc_cache[[k]]
      if (is.null(v)) {
        v <- maps$child_desc(ev, ch)
        assign(k, v, envir = desc_cache)
      }
      v
    }
    row_sets <- lapply(seq_len(nrow(fc)), function(i) {
      u <- fc$unit_id[i]; ch <- fc$child[i]; ev <- fc$event_node[i]
      foc_idx <- desc_of(ev, ch)
      others <- if ("sibling_children" %in% names(fc)) {
        selected <- as.integer(strsplit(
          as.character(fc$sibling_children[i]), ",", fixed = TRUE)[[1L]])
        selected[is.finite(selected)]
      } else {
        setdiff(ev_children[[as.character(u)]], ch)
      }
      oth_idx <- unique(unlist(lapply(others, function(o) desc_of(ev, o))))
      list(focal = set_label(foc_idx, maps$leaf_branch),
           other = set_label(oth_idx, maps$leaf_branch))
    })
    fc[, focal_set := vapply(row_sets, `[[`, character(1L), "focal")]
    fc[, other_set := vapply(row_sets, `[[`, character(1L), "other")]
    fc[, unit_kind := "event_child"]
  } else if (is_event_run) {
    # L0: child == lineage index; unit_id indexes sorted focal lineages.
    foc_of_unit <- maps$foc_order
    unit_foc <- foc_of_unit[fc$unit_id]
    fc[, event_conflict := unit_foc %in% maps$conflict_focs]
    kids_by_unit <- lapply(split(fc$child, fc$unit_id), unique)
    bad <- fc$child != unit_foc
    fc[bad, duplicate_emission := TRUE]
    kept <- which(!bad)
    row_sets <- lapply(kept, function(i) {
      u <- fc$unit_id[i]
      foc <- foc_of_unit[u]
      oth <- maps$sibs_by_foc[[as.character(foc)]]
      list(focal = set_label(foc, maps$leaf_branch),
           other = set_label(oth, maps$leaf_branch))
    })
    fc[kept, focal_set := vapply(row_sets, `[[`, character(1L), "focal")]
    fc[kept, other_set := vapply(row_sets, `[[`, character(1L), "other")]
    fc[, unit_kind := "focal_terminal"]
  } else {
    # Baseline/B/0.6.x and simple_lr: focal is `lineage`, siblings come from
    # sibling_map at the same fork_cluster. A row whose (lineage,
    # fork_cluster) has no mapped siblings keeps an undefined other set and
    # cannot hit a target.
    sib_idx <- lapply(seq_len(nrow(fc)), function(i) {
      f <- fc$lineage[i]
      s <- if ("fork_cluster" %in% names(fc)) {
        unique(maps$smap_siblings[[paste(f, fc$fork_cluster[i])]])
      } else {
        maps$sibs_by_foc[[as.character(f)]]
      }
      if (is.null(s)) integer(0L) else s
    })
    fc[, focal_set := maps$leaf_branch[lineage]]
    fc[, other_set := vapply(sib_idx, function(s) {
      set_label(s, maps$leaf_branch)
    }, character(1L))]
    fc[, unit_kind := "focal_terminal"]
  }
  invisible(fc)
}

# Build predictions while retaining wrong-fork and mirror candidates.
detail <- rbindlist(lapply(seq_len(nrow(jobs)), function(j) {
  ds <- jobs$dataset[j]
  ver <- jobs$version[j]
  is_baseline <- !ver %in% names(CFG$ws)
  ws_key <- data.table::fcase(
    ver %in% c("simple_lr_L1D0", "simple_lr_L1D0_smooth"), "L1D0",
    is_baseline, "0.6.5", default = ver)
  rds_path <- file.path(CFG$ws[[ws_key]], ds, "linealinkr_complete.rds")
  truth_path <- file.path(CFG$sim_dir, sprintf("%s_truth.csv", ds))
  if (!file.exists(rds_path)) {
    jobs[j, result_status := "missing_input"]
    return(NULL)
  }
  x <- readRDS(rds_path)
  meta_path <- file.path(EXTRACT_ROOT,ds,"extract",sprintf("%s_meta.csv",ds))
  if (!file.exists(meta_path)) {
    jobs[j, result_status := "missing_input"]
    return(NULL)
  }
  jobs[j, result_status := "available"]

  meta <- fread(meta_path)
  truth <- fread(truth_path)
  fc <- if (is_baseline) {
    bl <- fread(data.table::fcase(
      ver %in% c("simple_lr_L1D0", "simple_lr_L1D0_smooth"),
      baseline_path_L1D0, default = baseline_path))
    as.data.table(bl[dataset == ds & method == ver])
  } else {
    as.data.table(x$fork_contrasts)
  }
  if (!"sibling_contrast" %in% names(fc)) {
    jobs[j, result_status := "missing_fork_output"]
    return(NULL)
  }
  if (!nrow(fc)) return(NULL)
  rank_field <- intersect(c("event_ranking_score", "donor_stratified_score", "branch_balance_score", "persistent_positive_score",
                            "sibling_contrast"),
                          names(fc))[1L]
  rank_field <- Sys.getenv("LLR_EVAL_SCORE_FIELD", rank_field)
  fc[, ranking_score := get(rank_field)]
  switch(as.character(nzchar(SCORE_OVERRIDE_FILE)), "TRUE" = {
    override <- SCORE_OVERRIDES[dataset == ds & version == ver][order(source_row)]
    identity_fields <- c("unit_id", "pair_idx", "source", "target", "event_node", "child")
    stopifnot(nrow(override) == nrow(fc),
      identical(as.integer(override$source_row), seq_len(nrow(fc))),
      identical(unname(as.matrix(override[, ..identity_fields])),
        unname(as.matrix(fc[, ..identity_fields]))))
    fc[, ranking_score := override$ranking_score_override]
  })

  pmap <- as.data.table(x$pair_gene_map)[, .(
    pair_idx,
    ligand = as.character(ligand_genes),
    receptor = as.character(receptor_genes)
  )]
  fc <- merge(fc, pmap, by = "pair_idx", all.x = TRUE, sort = FALSE)
  truth_pair <- truth[, .(
    ligand = as.character(ligand), receptor = as.character(receptor),
    role, true_positive, donor_confounded,
    truth_focal_branch = as.character(focal_branch),
    truth_other_branch = as.character(other_branch)
  )]
  fc <- merge(fc, truth_pair, by = c("ligand", "receptor"),
              all.x = TRUE, sort = FALSE)

  ct_levels <- sort(unique(as.character(meta$clusters)))
  source_as_label <- ifelse(
    as.character(fc$source) %in% ct_levels,
    as.character(fc$source),
    ct_levels[as.integer(fc$source)]
  )
  target_as_label <- ifelse(
    as.character(fc$target) %in% ct_levels,
    as.character(fc$target),
    ct_levels[as.integer(fc$target)]
  )
  fc[, `:=`(
    dataset = ds,
    version = ver,
    source_label = source_as_label,
    target_label = target_as_label
  )]
  fc[, sender_direction := source_label == "SRC" & target_label == "TGT"]

  maps <- BuildRunMaps(x, meta)
  AttachSets(fc, maps, ver)

  fc[, hit_target_key := ifelse(
    true_positive %in% TRUE & !is.na(ligand),
    paste(dataset, ligand, receptor, "SRC", "TGT",
          truth_focal_branch, truth_other_branch, sep = "|"),
    NA_character_
  )]
  fc[, truth_row := true_positive %in% TRUE & sender_direction &
       !is.na(focal_set) & focal_set == truth_focal_branch &
       mapply(function(o, need) isTRUE(need %in% strsplit(o, "+", fixed = TRUE)[[1]]),
              other_set, truth_other_branch)]
  fc[, mirror_row := true_positive %in% TRUE & sender_direction &
       !is.na(focal_set) & focal_set == truth_other_branch &
       mapply(function(o, need) isTRUE(need %in% strsplit(o, "+", fixed = TRUE)[[1]]),
              other_set, truth_focal_branch)]
  fc[, resource_extra := is.na(role)]
  fc[, null_row := sender_direction & !(truth_row %in% TRUE) &
       !(resource_extra %in% TRUE)]
  fc[, false_fork_row := null_row & !(mirror_row %in% TRUE)]
  fc[, prediction_key := ifelse(
    !is.na(ligand) & !is.na(receptor) & !is.na(focal_set) &
      !is.na(other_set) & !is.na(source_label) & !is.na(target_label),
    paste(dataset, ligand, receptor, source_label, target_label,
          focal_set, other_set, sep = "|"),
    paste(dataset, "unmapped_row", seq_len(.N), sep = "|")
  )]
  if (!"evaluable" %in% names(fc)) {
    fc[, evaluable := is.finite(sibling_contrast)]
  }
  if (!"reason" %in% names(fc)) fc[, reason := NA_character_]
  fc
}), fill = TRUE)

fwrite(detail, file.path(CFG$out_dir, "lrsim_contrast_detail.csv"))

# One prediction per biological candidate key. Duplicate emissions for the
# same object keep the highest score under one fixed rule; wrong-fork and
# mirror keys stay as their own candidates.
scored_all <- detail[
  duplicate_emission %in% c(FALSE, NA) &
    evaluable %in% TRUE & is.finite(ranking_score) &
    sender_direction %in% TRUE
]
setorder(scored_all, dataset, version, prediction_key,
         -ranking_score, lineage, pair_idx)
scored_all <- scored_all[, .SD[1L], by = .(dataset, version, prediction_key)]

# Fixed denominator: every dataset x variant receives the same independent
# truth target count, including variants with no returned predictions.
# Recovery is counted on the hit truth targets, not on candidate objects:
# several distinct returned objects can satisfy the same target.
hit_targets <- scored_all[truth_row %in% TRUE,
                          .(hit_target_key = unique(hit_target_key)),
                          by = .(dataset, version)]
target_availability <- merge(
  variant_grid[, .(dataset, version, n_targets)],
  hit_targets[, .(n_targets_available = uniqueN(hit_target_key)),
              by = .(dataset, version)],
  by = c("dataset", "version"), all.x = TRUE
)
target_availability[is.na(n_targets_available), n_targets_available := 0L]
target_availability[, target_availability := ifelse(n_targets > 0, n_targets_available / n_targets, NA_real_)]

# Direction per target: sign of the target's top-scored correct candidate.
dir_targets <- scored_all[
  truth_row %in% TRUE,
  .SD[which.max(ranking_score)],
  by = .(dataset, version, hit_target_key)
][, .(n_direction_correct = sum(sibling_contrast > 0),
  n_targets_positive_score = sum(ranking_score > 0)),
  by = .(dataset, version)]

# Metric universe: candidates for the 12 truth-defined LR pairs. Resource
# pairs outside the truth table are reported separately as undefined.
scored_main <- scored_all[resource_extra %in% c(FALSE, NA)]
metric_groups <- split(scored_main,
                       paste(scored_main$dataset, scored_main$version))
metrics <- rbindlist(lapply(metric_groups, function(d) {
  pos <- d[truth_row %in% TRUE]
  nul <- d[null_row %in% TRUE]
  rankset <- d[truth_row %in% TRUE | null_row %in% TRUE]
  auroc <- if (nrow(pos) >= 2 && nrow(nul) >= 2) {
    as.numeric(pROC::auc(pROC::roc(rankset$truth_row %in% TRUE,
                                   rankset$ranking_score, quiet = TRUE,
                                   direction = "<")))
  } else NA_real_
  auprc <- if (nrow(pos) >= 2 && nrow(nul) >= 2) {
    PRROC::pr.curve(
      scores.class0 = rankset$ranking_score[rankset$truth_row %in% TRUE],
      scores.class1 = rankset$ranking_score[!(rankset$truth_row %in% TRUE)],
      curve = FALSE
    )$auc.integral
  } else NA_real_
  data.table(
    dataset = d$dataset[1], version = d$version[1],
    n_scored = nrow(d), n_pos_scored = nrow(pos), n_null_scored = nrow(nul),
    n_mirror_scored = sum(d$mirror_row %in% TRUE),
    auroc = auroc, auprc = auprc,
    pos_median = if (nrow(pos)) median(pos$sibling_contrast) else NA_real_,
    pos_direction_accuracy = if (nrow(pos)) mean(pos$sibling_contrast > 0) else NA_real_,
    null_median_abs = if (nrow(nul)) median(abs(nul$sibling_contrast)) else NA_real_,
    null_max_abs = if (nrow(nul)) max(abs(nul$sibling_contrast)) else NA_real_,
    null_rank_median = if (nrow(nul)) median(nul$ranking_score) else NA_real_,
    null_rank_max_abs = if (nrow(nul)) max(abs(nul$ranking_score)) else NA_real_,
    null_rank_positive = sum(nul$ranking_score > 0)
  )
}), fill = TRUE)

# Oracle K comes from the independent truth grid, never from returned positives.
pk <- rbindlist(lapply(split(scored_main, paste(scored_main$dataset,
                                                scored_main$version)), function(d) {
  k_budget <- variant_grid[dataset == d$dataset[1] & version == d$version[1],
                           n_targets][1]
  d <- d[order(-ranking_score, prediction_key)]
  k_used <- min(k_budget, nrow(d))
  data.table(
    dataset = d$dataset[1], version = d$version[1],
    k_budget = k_budget, k_used = k_used, n_returned = nrow(d),
    precision_at_k = if (k_budget > 0) {
      sum(d$truth_row[seq_len(k_used)]) / k_budget
    } else NA_real_
  )
}), fill = TRUE)

false_forks <- detail[
  sender_direction %in% TRUE & false_fork_row %in% TRUE &
    evaluable %in% TRUE & is.finite(sibling_contrast),
  .(n_false_fork_returns = .N,
    false_fork_median_abs = median(abs(sibling_contrast)),
    false_fork_max_abs = max(abs(sibling_contrast))),
  by = .(dataset, version)
]
coverage <- detail[, .(
  n_total = .N,
  n_evaluable = sum(evaluable %in% TRUE),
  n_duplicate_emission = sum(duplicate_emission %in% TRUE),
  n_resource_extra = sum(resource_extra %in% TRUE & sender_direction %in% TRUE),
  coverage = mean(evaluable %in% TRUE)
), by = .(dataset, version)]
reasons <- detail[
  sender_direction %in% TRUE & !(evaluable %in% TRUE),
  .(n = .N), by = .(dataset, version, reason)
]
metric_grid <- variant_grid[, .(dataset, version)]
out <- Reduce(function(a, b) merge(a, b, by = c("dataset", "version"), all = TRUE),
              list(metric_grid, target_availability, coverage, metrics, dir_targets,
                   false_forks, pk))
out <- merge(out, detail[, .(n_returned_rows = .N), by = .(dataset, version)],
             by = c("dataset", "version"), all = TRUE)

out <- merge(out, jobs, by = c("dataset", "version"), all.x = TRUE)
# A completed method returning no candidates has zero counts and zero P@K.
# Missing runs stay missing; conditional AUROC/AUPRC remain undefined without
# both positive and negative scores.
count_columns <- c("n_total", "n_evaluable", "n_duplicate_emission",
  "n_resource_extra", "n_scored", "n_pos_scored", "n_null_scored",
  "n_mirror_scored", "n_direction_correct", "n_targets_positive_score",
  "null_rank_positive", "n_false_fork_returns",
  "k_used", "n_returned", "n_returned_rows")
out[result_status == "available", (count_columns) := lapply(.SD, nafill, fill = 0),
  .SDcols = count_columns]
out[result_status == "available", k_budget := n_targets]
out[result_status == "available" & n_targets > 0 & n_returned == 0,
  precision_at_k := 0]
out[result_status == "available" & n_total == 0, coverage := 0]
out[result_status != "available", `:=`(n_targets_available = NA_integer_, target_availability = NA_real_,
                                         n_direction_correct = NA_integer_,
                                         n_targets_positive_score = NA_integer_)]
target_availability <- merge(target_availability, jobs, by = c("dataset", "version"), all.x = TRUE)
target_availability[result_status != "available", `:=`(n_targets_available = NA_integer_, target_availability = NA_real_)]
target_availability <- merge(target_availability,out[,.(dataset,version,n_targets_positive_score)],
  by=c("dataset","version"),all.x=TRUE)
# Per-target output uses the same candidate ordering and planted-target budget
# as native fork-event precision. Availability is not called detection.
target_ranked <- copy(scored_main)
setorderv(target_ranked,c("dataset","version","ranking_score","prediction_key"),c(1L,1L,-1L,1L))
target_ranked[,candidate_rank:=seq_len(.N),by=.(dataset,version)]
best_correct <- target_ranked[truth_row %in% TRUE,.SD[1L],
  by=.(dataset,version,hit_target_key)][,.(dataset,version,target_key=hit_target_key,
    best_correct_prediction=prediction_key,best_correct_rank=candidate_rank,
    ranking_score,local_contrast=sibling_contrast)]
target_rows <- merge(truth_base,variant_grid[,.(dataset,version,k_budget=n_targets)],
  by="dataset",allow.cartesian=TRUE)
target_rows <- merge(target_rows,best_correct,by=c("dataset","version","target_key"),all.x=TRUE)
target_rows[,c("correct_candidate_available","positive_ranking_score",
  "correct_local_direction","hit_in_budget"):=.(
    is.finite(ranking_score),is.finite(ranking_score) & ranking_score>0,
    is.finite(local_contrast) & local_contrast>0,
    is.finite(best_correct_rank) & best_correct_rank<=k_budget)]
target_rows <- merge(target_rows,jobs,by=c("dataset","version"),all.x=TRUE)
target_rows[result_status!="available",c("correct_candidate_available",
  "positive_ranking_score","correct_local_direction","hit_in_budget"):=NA]
budget_hits <- target_rows[,.(n_targets_in_budget=sum(hit_in_budget %in% TRUE)),
  by=.(dataset,version)]
out <- merge(out,budget_hits,by=c("dataset","version"),all.x=TRUE)
out[result_status=="available" & n_targets==0,n_targets_in_budget:=0L]
out[,budget_recall:=ifelse(n_targets>0,n_targets_in_budget/n_targets,NA_real_)]
out[,n_targets_correct_direction:=n_direction_correct]
out[,target_direction_accuracy:=ifelse(n_targets>0,n_targets_correct_direction/n_targets,NA_real_)]
out[result_status!="available",c("n_targets_in_budget","budget_recall",
  "n_targets_correct_direction","target_direction_accuracy"):=NA]

setorder(target_rows,dataset,version,target_key)
fwrite(target_rows, file.path(CFG$out_dir, "lrsim_target_recovery.csv"))
fwrite(out, file.path(CFG$out_dir, "lrsim_metrics.csv"))
fwrite(reasons, file.path(CFG$out_dir, "lrsim_nocall_reasons.csv"))
cat("\nSaved:", file.path(CFG$out_dir, "lrsim_metrics.csv"), "\n")
