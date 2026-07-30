# Shared 55 GiB process budget for R benchmark adapters.
# Estimate peak RAM from data dimensions before allocating.
# The hard backstop (Windows Job Object) is applied by the orchestrator, which runs
# each Rscript child through adapters/_mem_guard.py --limit-gb <budget>.
MEM_BUDGET_GB <- 55

mem_dense_gb <- function(n_cells, n_genes, bytes_per = 8, copies = 1) {
  as.numeric(n_cells) * as.numeric(n_genes) * bytes_per * copies / (1024^3)
}
mem_square_gb <- function(n, bytes_per = 8, copies = 1) {
  as.numeric(n) * as.numeric(n) * bytes_per * copies / (1024^3)
}
mem_require_under_gb <- function(est_gb, limit_gb = MEM_BUDGET_GB, what = "") {
  if (is.finite(est_gb) && est_gb > limit_gb) {
    stop(sprintf(paste0("FAIL: estimated peak memory %.1f GB exceeds %g GB budget for %s."),
                 est_gb, limit_gb, what), call. = FALSE)
  }
}
