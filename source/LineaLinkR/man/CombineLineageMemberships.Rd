\name{CombineLineageMemberships}
\alias{CombineLineageMemberships}
\title{Combine geometric support with local branch assignment}
\usage{
CombineLineageMemberships(reference_weights, assignment_weights)
}
\arguments{
\item{reference_weights}{Cell-by-lineage geometric memberships.}
\item{assignment_weights}{Matching local branch memberships.}
}
\value{Cell-by-lineage normalized propensities in the original path order.}
\description{
Normalize the componentwise maximum of the two row-normalized inputs.
The result preserves their union of positive support, and retains at least
half of each input's mass. It is a branch propensity, not a calibrated
individual-cell fate posterior. Geometric memberships remain the reference
population weights for communication contrasts.
}
