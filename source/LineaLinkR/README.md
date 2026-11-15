# LineaLinkR

RNA-state trajectories, lineage memberships and branch-relative ligand–receptor activity in R.

Source installation requires R >= 4.2, C++17 and Rust/Cargo >= 1.80. On Windows, use Rtools matching R and the `x86_64-pc-windows-gnu` Rust toolchain.

For Windows HTTPS downloads, run `options(download.file.method = "curl")` in the R session before installation.

Install from GitHub:

```r
install.packages("remotes")
remotes::install_github("aeqby770/manuscript-materials", subdir = "source/LineaLinkR", auth_token = NULL)
library(LineaLinkR)
help(package = "LineaLinkR")
```

To install this source directory with its dependencies already available, run `R CMD INSTALL .` from the package directory.

See the [installation guide and example](https://github.com/aeqby770/manuscript-materials#installation) for input formats and a complete example.
