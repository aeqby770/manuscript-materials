\name{AddLateForkPersistence}
\alias{AddLateForkPersistence}
\title{Attach matched late support to fork contrasts}
\usage{
AddLateForkPersistence(
  fork_contrasts, temporal_lr, lineage_paths,
  expression_mean_mode = c("all_cells", "positive_cells"),
  sibling_scope = c("native", "pooled_native", "all_descendants", "pairwise"),
  persistence_mode = c("signed", "positive"),
  fallback_temporal_lr = temporal_lr$late_fallback
)
}
\arguments{
\item{fork_contrasts}{Event table from \code{ComputeForkContrastsEvent}.}
\item{temporal_lr}{Binned temporal ligand--receptor activity and bin boundaries.}
\item{lineage_paths}{Root-to-terminal lineage paths.}
\item{expression_mean_mode}{Use all-cell or positive-cell temporal activity.}
\item{sibling_scope}{Observed cell-mass pooling within the fixed local sibling panel (default); a strict common horizon with equal child weights; historical unrestricted descendants; or pairwise observed horizons with equal sibling weights.}
\item{fallback_temporal_lr}{Optional adjacent-bin temporal grid scored from the same original cells, pooled Hill scale and unchanged effective-count threshold. Fine observations take priority.}
\item{persistence_mode}{Symmetric signed persistence (default), or the historical positive-only cap.}
}
\description{Attaches signed late support to a shared-parent local contrast. The default pools observed cells within the fixed local sibling panel; a short sibling does not set the horizon for every other sibling.}
\details{For local contrast theta and late change contrast lambda, signed persistence is sign(theta) times min(abs(theta), max(0, sign(theta) times lambda)). When fine late support is absent, a supplied adjacent-bin grid can restore observed support before the missing-support rule is used. persistence_support records late_observed_adjacent_bins for those rows. Without support on either grid, signed mode uses sign(theta) times min(abs(theta), max(0, sign(theta) times focal local change)). This preserves local direction while requiring its magnitude to be supported by the focal change; missing late support is not treated as full persistence. persistence_support distinguishes local_focal_only, late_observed and not_evaluable. The positive-only historical mode retains its local-contrast fallback. The historical field name persistent_positive_score contains signed values in the default mode. persistence_phase_start and persistence_phase_end describe the bins actually used by persistence; late_*_native fields retain the strict common-horizon diagnostic. phase_pattern follows the reference actually used for persistence. The pairwise option is retained as a tested mechanism control. This function does not fit the sample-conditioned pulse model or produce calibrated event probabilities.}
\value{The event table with late activity, baseline-relative change, persistence and phase-support fields.}
\seealso{\code{ComputeForkContrastsEvent}}
