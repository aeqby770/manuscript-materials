\name{ComputeSourceRemovalLR}
\alias{ComputeSourceRemovalLR}
\title{Fixed-composition RNA source-removal query}
\usage{
ComputeSourceRemovalLR(expr, celltypes, lr_db, sample_labels, source,
  target, hill_scale, min_n_eff = 5)
}
\arguments{
\item{expr}{Genes by cells transformed RNA matrix used for LR scoring.}
\item{celltypes}{Integer population code per cell.}
\item{lr_db}{Lists of ligand and receptor subunit row indices.}
\item{sample_labels}{Biological sample identity per cell.}
\item{source,target}{Source and receiving population codes.}
\item{hill_scale}{Common positive LR saturation scale.}
\item{min_n_eff}{Minimum source and receiving cells per sample.}
}
\description{Removes the contribution of one source population from a fixed
RNA mixture while holding receptor expression and the original total-cell
denominator fixed. Complexes use their limiting averaged subunit.}
\details{Source loss is the difference between total and source-removed
log-exposure. Source loss fraction divides this difference by total activity.
The query describes RNA availability under the fixed-composition operation.}
\value{Per-sample, per-pair availability, loss, loss fraction and cell support.}
