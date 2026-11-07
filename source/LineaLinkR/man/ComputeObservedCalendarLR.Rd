\name{ComputeObservedCalendarLR}
\alias{ComputeObservedCalendarLR}
\title{Observed calendar LR activity in native branch subtrees}
\usage{
ComputeObservedCalendarLR(
  expr, celltypes, lr_db, observed_time, sample_labels,
  lineage_clusters, lineage_paths, hill_scale, min_pct = 0.05,
  min_n_eff = 5, requested_contexts = NULL,
  condition_labels = rep.int("all", ncol(expr))
)
}
\arguments{
\item{expr, celltypes, lr_db, min_pct, min_n_eff}{See \code{ComputeTemporalLR}.}
\item{observed_time}{Numeric sampling time for every cell.}
\item{sample_labels}{Observation or biological-sample labels, with one time per label.}
\item{lineage_clusters}{Native tree-node label for every cell.}
\item{lineage_paths}{Root-to-leaf native tree paths.}
\item{hill_scale}{Fixed pooled Hill saturation scale.}
\item{requested_contexts}{Optional pair/source/target requests.}
\item{condition_labels}{Experimental-condition labels for every cell.}
}
\value{
A temporal-LR list with scores, shared calendar boundaries, native subtree
paths, profile-node labels, and the observed-population scope.
}
\description{
Each child profile uses its observed descendant cells. Expression and detection
are pooled before the nonlinear LR transform. Nested subtrees retain their raw
memberships, without distributing ancestor cells to future children.
\code{SummarizeForkDynamics} selects the direct child profiles for each event.
}

