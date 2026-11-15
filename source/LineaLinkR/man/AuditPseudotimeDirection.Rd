\name{AuditPseudotimeDirection}
\alias{AuditPseudotimeDirection}
\title{Audit whether pseudotime direction is supported by blind diagnostics}
\usage{
AuditPseudotimeDirection(
  pseudotime,
  potency = NULL,
  embedding = NULL,
  n_pc = 3L,
  endpoint_fraction = 0.05,
  min_potency_delta = 0.1,
  min_blind_delta = 0.05
)
}
\arguments{
\item{pseudotime}{Numeric vector of per-cell pseudotime values.}

\item{potency}{Optional numeric vector; higher values are interpreted as more root-like for this audit.}

\item{embedding}{Optional numeric matrix used for the blind PC score.}

\item{n_pc}{Number of PCs passed to \code{BlindEmbeddingScore()}.}

\item{endpoint_fraction}{Lower/upper pseudotime fraction used as root and terminal endpoint summaries.}

\item{min_potency_delta}{Minimum scaled endpoint potency difference.}

\item{min_blind_delta}{Minimum blind-score margin for a PC-convention call.}
}
\value{
Named list with potency, PC-convention and combined status fields.
}
\description{
This diagnostic does not flip pseudotime. It separates a potency prior supplied by the caller from the sign-anchored PC convention used by \code{BlindEmbeddingScore()}. The PC convention can detect internally conflicting orientations, but it is not a biological ground truth.
}
