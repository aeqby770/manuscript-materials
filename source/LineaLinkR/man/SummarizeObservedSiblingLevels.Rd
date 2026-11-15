\name{SummarizeObservedSiblingLevels}
\alias{SummarizeObservedSiblingLevels}
\title{Sample-paired levels across observed receiver subtrees}
\usage{
SummarizeObservedSiblingLevels(temporal_lr, event_units,
  activity_transform = c("log_exposure", "hill"))
}
\arguments{
\item{temporal_lr}{Output of ComputeObservedCalendarLR, including observations,
requested contexts and native profile-node identities.}
\item{event_units}{The event-unit attribute of a fork-contrast result.}
\item{activity_transform}{The original Hill scale or log exposure.}
}
\value{A list containing cohort effects, individual observed-sample differences,
coverage and the declared sampling unit.}
\description{Each focal subtree is compared with the equal mean of its supported
structural sisters in the same observed sample. The cohort effect averages these
paired differences. A sample-wide sender population can be declared when
computing the upstream observed profiles. The output retains unobserved requested
contexts with missing effects and zero paired-sample support. Finite zero scores
are observed opportunities on the declared upstream activity scale and remain
in the paired comparison.}
