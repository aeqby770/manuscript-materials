# Programmatic method-data applicability rules (algorithm x dataset).
# Applicability rules depend only on required modality, species, and reference scope.

panel_dataset_species <- function(ds) {
  sp <- CONFIG$dataset_species[[ds]]
  if (is.null(sp)) "human" else sp
}

panel_is_simulated <- function(ds) {
  sp <- panel_dataset_species(ds)
  sp == "simulated" ||
    grepl("^sim_", ds) ||
    grepl("^simulated", ds)
}

panel_is_human_mouse <- function(ds) {
  panel_dataset_species(ds) %in% c("human", "mouse")
}

# Return skip reason string, or NA if this pair should run.
panel_skip_rule_reason <- function(algo, ds) {
  sp <- panel_dataset_species(ds)
  sim <- panel_is_simulated(ds)
  hm <- panel_is_human_mouse(ds)

  if (algo == "CellPhoneDB" && sp != "human") {
    return("non-human (no ortholog DB)")
  }

  ccc_hm_only <- c("NATMI", "LIANA", "nichenetr", "CytoTalk", "CellChat")
  if (algo %in% ccc_hm_only) {
    if (sim) return("simulated (no human/mouse L-R DB)")
    if (!hm) return("non-human/mouse species")
  }

  if (algo == "SingleCellSignalR") {
    if (sim) return("simulated (no applicable L-R)")
    if (!hm) return("non-human/mouse species")
  }

  if (algo %in% c("TimeTalk", "TraSig")) {
    if (sim) return("simulated (no DoRothEA temporal CCC)")
    if (!hm) return("non-human/mouse species (no DoRothEA/CCC DB)")
    tp_ok <- panel_has_timepoints(ds)
    if (isFALSE(tp_ok)) return("no temporal metadata (need >= 2 timepoints)")
  }

  if (algo == "wot") {
    tp_ok <- panel_has_timepoints(ds)
    if (isFALSE(tp_ok)) return("requires discrete timepoint metadata")
  }

  NA_character_
}

panel_has_timepoints <- function(ds, min_n = 2L) {
  mp <- file.path(CONFIG$dir$r_export, ds, "meta.csv")
  if (!file.exists(mp)) return(NA)
  hdr <- names(data.table::fread(mp, nrows = 0))
  candidates <- c("timepoint", "time_point", "Timepoint", "time_info", "stage", "theiler",
                  "age(days)", "age_days", "age")
  for (col in candidates) {
    if (!col %in% hdr) next
    v <- as.character(data.table::fread(mp, select = col, showProgress = FALSE)[[col]])
    v <- unique(v[!is.na(v) & v != ""])
    if (length(v) >= min_n) return(TRUE)
  }
  FALSE
}

merge_skip_by_algo <- function(base, extra) {
  algos <- unique(c(names(base), names(extra)))
  out <- list()
  for (a in algos) {
    out[[a]] <- unique(c(base[[a]], extra[[a]]))
  }
  out
}

panel_should_skip <- function(algo, ds, skip_datasets, skip_by_algo) {
  if (length(skip_datasets) && ds %in% skip_datasets) {
    return(list(skip = TRUE, reason = "panel skip list (global)"))
  }
  algo_skip <- skip_by_algo[[algo]]
  if (length(algo_skip) && ds %in% algo_skip) {
    return(list(skip = TRUE, reason = "panel skip list"))
  }
  rule <- panel_skip_rule_reason(algo, ds)
  if (!is.na(rule)) {
    return(list(skip = TRUE, reason = rule))
  }
  list(skip = FALSE, reason = NA_character_)
}
