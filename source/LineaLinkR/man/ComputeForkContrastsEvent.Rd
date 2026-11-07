\name{ComputeForkContrastsEvent}
\alias{ComputeForkContrastsEvent}
\title{Estimate ligand--receptor changes at native fork events}
\usage{
ComputeForkContrastsEvent(
  expr, celltypes, lr_db, pseudotime, lineage_weights,
  lineage_paths = NULL, sibling_map = NULL, lineage_clusters = NULL,
  fork_bins = NULL, bin_bounds = NULL,
  event_scope = c("all", "nearest_sibling"),
  boundary = c("exact", "fixed_bin", "tree_segment_interval"),
  window_bins = 2L, hill_scale = NULL, min_n_eff = 5,
  min_pct = 0.05, activity_floor = 1e-4, bif_df = NULL,
  units_override = NULL, descendant_window_children = integer(),
  anchor_rule = c("min", "mean", "child_mass_entropy"),
  scoring_mode = c("bin", "interval"),
  expression_mean_mode = c("all_cells", "positive_cells"),
  activity_mode = c("thresholded", "continuous"),
  sibling_selection = c("all", "time_matched", "program_matched"),
  sibling_program_profiles = NULL,
  activity_transform = c("hill", "log_exposure"),
  measurement_error = FALSE, measurement_median_bounds = NULL,
  prefork_weighting = c("shared", "descendant")
)
}
\arguments{
\item{expr}{Log-normalized sparse gene-by-cell expression matrix.}
\item{celltypes}{Integer receiver/sender identity for each cell.}
\item{lr_db}{Validated ligand--receptor database for the expression matrix.}
\item{pseudotime}{One state coordinate per cell.}
\item{lineage_weights}{Cell-by-terminal-lineage soft assignments.}
\item{lineage_paths}{Root-to-terminal cluster paths.}
\item{sibling_map}{Nearest-sibling comparison map.}
\item{lineage_clusters}{Integer cluster membership for each cell.}
\item{fork_bins,bin_bounds,bif_df}{Saved fork crossings and temporal bin boundaries.}
\item{event_scope}{All tree fork events or nearest-sibling units.}
\item{boundary}{Exact, fixed-bin or native tree-segment windows.}
\item{window_bins}{Number of bins around a fork for bin-based windows.}
\item{hill_scale}{Positive dataset-level Hill scale; required.}
\item{min_n_eff,min_pct,activity_floor}{Effective-cell, detection and activity support thresholds.}
\item{units_override}{The pooled result's event-unit attribute, for sample-matched estimation.}
\item{descendant_window_children}{Children whose post segment includes later descendants.}
\item{anchor_rule}{How child crossing coordinates determine the common event anchor.}
\item{scoring_mode}{Binned or interval scoring.}
\item{expression_mean_mode}{Mean over all or positive cells.}
\item{activity_mode}{Thresholded or continuous activity.}
\item{sibling_selection}{Supported sibling panel selection.}
\item{sibling_program_profiles}{Profiles used by program-matched sibling selection.}
\item{activity_transform}{Native Hill activity or log exposure.}
\item{measurement_error}{Estimate conditional delete-one-cell measurement error on the selected Hill or log-exposure scale, with geometry, support and gates held fixed.}
\item{measurement_median_bounds}{Optional fixed bounds for propagating that error through a sample median.}
\item{prefork_weighting}{For native intervals, shared observed parent weights (default) or retrospective descendant weights.}
}
\description{Scores each focal child and its supported siblings on the same native event windows.}
\details{With shared parent weights, every child uses the same observed pre-fork population. The sibling difference estimates divergence after a common parent; it does not identify true ancestry or a causal effect. Unsupported rows retain an explicit reason.}
\value{A table with event, child, molecular-context identity, activity, contrast and support fields. The event-unit attribute records the comparison windows and weights.}
\seealso{\code{AddLateForkPersistence}}
