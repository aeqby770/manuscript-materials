\name{RefineStateGraph}
\alias{RefineStateGraph}
\title{Observed-neighbourhood RNA-state geometry}
\usage{
RefineStateGraph(embedding, knn_conn, start_cell, k = 15L,
  n_neighbors = 50L, seed = 12345L, project_root = TRUE)
}
\arguments{
\item{embedding}{Numeric cells-by-coordinates RNA embedding.}
\item{knn_conn}{Observed sparse cell connectivity.}
\item{start_cell}{Root-prior cell index.}
\item{k}{Number of coarse RNA states.}
\item{n_neighbors}{Maximum within-component neighbour count.}
\item{seed}{Random seed for state clustering.}
\item{project_root}{Project an expression-derived root onto an endpoint.}
}
\value{Refined geometry, observed graph, states, root, paths and weights.}
\description{Applies one graph-conditional neighbourhood average to all RNA coordinates. Original observed edges reconnect newly separated components below the retained refined supports. Originally disconnected components, including singleton states, retain their own paths and memberships.}
