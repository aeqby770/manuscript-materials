#' Lineage priming from within-state standardized count contrasts.
#'
#' Retains the complete state-discriminative subspace of log-normalized RNA.
#' State identities and paths are unchanged; clonal outcomes are not inputs.
#'
#' @param counts Gene-by-cell count matrix.
#' @param cluster_labels State label per cell.
#' @param lineage_paths Rooted paths in output column order.
#' @param expression_scale Total-count normalization scale.
#' @return Cell-by-lineage weights with the fitted count projection attached.
ComputeLineageWeightsCountResidual <- function(counts, cluster_labels,
                                               lineage_paths,
                                               expression_scale = 1e4) {
  expression <- as(counts, "CsparseMatrix")
  expression <- expression %*% Matrix::Diagonal(
    x = expression_scale / Matrix::colSums(expression))
  expression@x <- log1p(expression@x)
  states <- sort(unique(cluster_labels))
  state_index <- match(cluster_labels, states)
  state_size <- tabulate(state_index, length(states))
  membership <- Matrix::sparseMatrix(i = seq_along(state_index),
    j = state_index, x = 1,
    dims = c(length(state_index), length(states)))
  profile <- as.matrix(expression %*% membership) /
    rep(state_size, each = nrow(expression))
  within_variance <- (Matrix::rowSums(expression^2) -
    rowSums(profile^2 * rep(state_size, each = nrow(profile)))) /
    (ncol(expression) - length(states))
  retained <- which(within_variance > 0)
  noise_sd <- sqrt(within_variance[retained])
  contrast <- sweep(profile[retained, , drop = FALSE], 1L,
    rowMeans(profile[retained, , drop = FALSE]), "-") / noise_sd
  basis <- svd(contrast, nu = length(states) - 1L, nv = 0L)$u / noise_sd
  embedding <- as.matrix(Matrix::t(expression[retained, , drop = FALSE]) %*% basis)
  weights <- ComputeLineageWeightsLocalBranchResidual(embedding, cluster_labels,
    lineage_paths, priming = FALSE)
  attr(weights, "count_projection") <- list(genes = rownames(counts)[retained],
    loadings = basis, expression_scale = expression_scale,
    within_variance = within_variance[retained],
    method = "complete_within_state_standardized_mean_contrast_subspace")
  weights
}
