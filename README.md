# LineaLinkR

LineaLinkR is an R package for RNA-state trajectories, lineage memberships and branch-relative ligand–receptor activity. Functions return state assignments, lineage weights, temporal profiles and sample-matched interaction contrasts.

## Installation

Source installation requires R 4.2 or later, C++17, and Rust/Cargo 1.80 or later. On Windows, install the Rtools version matching R and the `x86_64-pc-windows-gnu` Rust toolchain. Make `cargo` and `rustc` available on `PATH` before installing.

For Windows HTTPS downloads, run `options(download.file.method = "curl")` in the R session before installation.

```r
install.packages("remotes")
remotes::install_github(
  "aeqby770/manuscript-materials",
  subdir = "source/LineaLinkR",
  auth_token = NULL,
  upgrade = "never"
)
library(LineaLinkR)
```

## Example

The embedding has one row per cell. The sparse connectivity matrix uses the same cell order.

```r
embedding <- rbind(
  c(0.0, 0.0), c(0.4, 0.0), c(0.9, 0.0), c(1.4, 0.6),
  c(1.8, 1.0), c(1.4, -0.6), c(1.8, -1.0), c(2.2, -1.2)
)
knn_conn <- Matrix::sparseMatrix(
  i = c(1, 2, 3, 4, 3, 6, 7),
  j = c(2, 3, 4, 5, 6, 7, 8),
  x = 1, dims = c(8, 8), symmetric = TRUE
)
states <- RefineStateGraph(
  embedding, knn_conn, start_cell = 1L,
  k = 4L, n_neighbors = 3L, seed = 12345L,
  project_root = FALSE
)
data.frame(cell = seq_len(nrow(embedding)), state = states$labels)
states$lineage_paths
states$lineage_weights
```

For ligand–receptor comparisons, construct the embedding and graph from RNA features excluding the LR subunits being evaluated. Expression matrices use genes as rows and cells as columns; expression, graph and cell metadata must share the same cell order.

## Functions

| Function | Output |
| --- | --- |
| `RefineStateGraph()` | RNA-state graph, paths and cell memberships |
| `ComputeLineageWeightsLocalBranchResidual()` | Local assignment or priming weights along supplied paths |
| `ComputeForkContrastsEvent()` | Focal-child and sibling contrasts on shared event windows |
| `SummarizeForkDynamics()` | Temporal effects, profile evidence and supported peak times |
| `ComputeSourceRemovalLR()` | Sender contributions to modeled LR activity |
| `PredictLineageComposition()` | Descendant-composition estimates from RNA neighbours and lineage information |

Argument definitions and return values are available in R:

```r
help(package = "LineaLinkR")
?RefineStateGraph
?ComputeForkContrastsEvent
?SummarizeForkDynamics
```

## License

MIT. See [LICENSE](LICENSE).
