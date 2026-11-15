\name{PredictLineageComposition}
\alias{PredictLineageComposition}
\title{Conditional descendant composition from terminal RNA references}
\usage{
PredictLineageComposition(query_memberships, reference_memberships,
  reference_outcomes, neighbors, neighbor_weights, return_components = FALSE)
}
\arguments{
\item{query_memberships}{Query cells by lineage matrix.}
\item{reference_memberships}{Terminal RNA reference cells by lineage matrix.}
\item{reference_outcomes}{Reference cells by outcome indicator matrix.}
\item{neighbors}{Matrix of RNA neighbor indices in the reference rows.}
\item{neighbor_weights}{Nonnegative RNA neighbor similarity weights.}
\item{return_components}{Return the likelihood, prior, effective count and context-prior control.}
}
\description{Combines local terminal-RNA composition with a weak lineage-derived
Dirichlet prior. The local effective count is the inverse sum of squared
normalized neighbor weights; one-half per outcome defines the prior concentration.}
\value{An outcome probability matrix, or a list of components.}
