\name{RankTransientForkEvents}
\alias{RankTransientForkEvents}
\title{Rank transient profiles within fork contexts}
\usage{
RankTransientForkEvents(
  fork_dynamics, component = c("relative", "focal"), n = 1L,
  group_by = c("event_node", "child", "source", "target"),
  aggregation = c("event", "mean_children")
)
}
\arguments{
\item{fork_dynamics}{Output from \code{SummarizeForkDynamics}.}
\item{component}{Relative branch activity or focal-child activity.}
\item{n}{Maximum number of returned molecular queries per context.}
\item{group_by}{Columns defining the context and query budget.}
\item{aggregation}{Individual event ranks, or mean bounded pulse evidence times amplitude over planned children.}
}
\description{Orders fitted pulse profiles by their BIC advantage over nonpulse
profiles. Selection uses RNA-derived profiles and a specified context budget.}
\details{The score summarizes temporal model support. Ranking alone provides
no significance threshold or false-discovery calibration.}
\value{Selected events or molecular queries with a \code{temporal_score} column. Aggregated queries also report planned, observed and pulse-supporting children.}
