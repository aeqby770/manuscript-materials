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
  time_coordinate = c("state_pseudotime", "observed_calendar"),
  condition_labels = NULL,
  reference_conditions = character(),
  sample_contrasts = NULL,
  transfer_pulse = TRUE,
  include_level = TRUE
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
\item{condition_labels}{Observed-condition labels, one per temporal score row.}
\item{reference_conditions}{Reference labels included in each condition-specific fit.}
\item{sample_contrasts}{Optional native sample contrasts for supplementary cohort summaries.}
\item{transfer_pulse}{Allow pulse completion in the optional cohort summary.}
\item{include_level}{Estimate a sustained post-parent branch difference together with temporal shape.}
}
\description{
Fits constant, linear, plateau, delayed, two-slope and flat-top pulse profiles on the same
observed native child panel. Relative and focal shapes are reported separately.
}
\details{
Native shared-parent state fits use the supplied fork coordinate as a common
origin and retain within-sample time backgrounds. A free branch level preserves
sustained differences. Sample-specific branch offsets have zero mean, preserving
the cohort mean while adjusting between-sample heterogeneity. The descendant-parent control
retains sample-by-child and sample-by-time intercepts on its observed horizon.
Continuous breakpoints are obtained from fitted line intersections within
adjacent observed-time intervals.
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
\code{focal_peak_time}, defined only for a fitted focal pulse.
\code{dynamic_effect} and \code{dynamic_effect_time} report the signed largest
fitted branch difference on observed support and its coordinate in native
shared-parent state mode. Focal and calendar controls retain fitted-change
semantics. The model level is recorded in \code{dynamic_level_effect}. \code{event_peak_time} reports the
maximum of the focal activity curve, with sample-intercept adjustment in
state mode. A state-coordinate interior maximum is refined by quadratic
interpolation over its two neighboring supported times. Calendar maxima
remain at observed sampling times. \code{event_peak_lower} and
\code{event_peak_upper} are neighboring support coordinates, rather than
confidence limits. A peak coordinate and support for a transient shape are
reported separately.
}
