\name{SummarizeForkDynamics}
\alias{SummarizeForkDynamics}
\title{Describe branch-associated temporal activity on observed support}
\usage{
SummarizeForkDynamics(
  fork_contrasts,
  temporal_lr,
  lineage_paths,
  expression_mean_mode = c("all_cells", "positive_cells"),
  background = c("branch_conditioned", "focal_only"),
  sample_labels = rep.int("pooled", nrow(temporal_lr$scores)),
  condition_on_sample = TRUE,
  time_coordinate = c("state_pseudotime", "observed_calendar")
)
}
\arguments{
\item{fork_contrasts}{Native event contrast table.}
\item{temporal_lr}{Temporal LR scores and their original bin bounds.}
\item{lineage_paths}{Root-to-terminal paths for those scores.}
\item{expression_mean_mode}{Activity field used by the native event analysis.}
\item{background}{Child and shared-time effects, or a focal-only comparison.}
\item{sample_labels}{Biological sample corresponding to each score row.}
\item{condition_on_sample}{Retain within-sample backgrounds in state mode.}
\item{time_coordinate}{State pseudotime by default; observed calendar time
requires a finite \code{observed_time} field and one collection time per
biological sample.}
}
\description{
Fits constant, linear, plateau, delayed and two-slope profiles on the same
observed native child panel. Relative and focal shapes are reported separately.
}
\details{
State mode preserves sample-by-child intercepts and within-sample time effects.
Calendar mode uses observed sampling times, sample plus child intercepts and
equal total weight per biological sample. A sample-by-child intercept would
absorb calendar interactions when each sample is collected once.
Weighted partial regression profiles out the baseline. Actual design rank and
remaining sample degrees of freedom determine whether a model is estimable.
Calendar profiles use their own observed LR support, independently of
state-coordinate pre/post event eligibility.
BIC gains are descriptive model support, without calibrated significance or a
causal interpretation. State coordinates and calendar time remain distinct;
missing support is retained. This function does not change event rankings.
}
\value{
The input table with \code{dynamic_*} shape, support, slopes and BIC fields,
\code{pulse_bic_gain} and \code{dynamic_time_coordinate}. Branch-conditioned
output also contains \code{focal_dynamic_*}, \code{focal_pulse_bic_gain} and
\code{focal_peak_time}, defined only for a focal pulse.
}
