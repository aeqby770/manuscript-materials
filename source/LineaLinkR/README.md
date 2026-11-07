# LineaLinkR 0.6.5 source package

This archive contains the R/Rust implementation, function documentation and package tests.

Source installation requires R >= 4.2, the R packages listed in DESCRIPTION, C++17, and Rust/Cargo >= 1.80. On Windows, use Rtools matching R and the GNU Rust host toolchain x86_64-pc-windows-gnu. Put that toolchain's cargo.exe and rustc.exe on PATH before installation.

From the directory containing the manuscript attachments, create a private R library and install the dependencies:

~~~r
dir.create("library",showWarnings=FALSE)
Sys.setenv(R_LIBS_USER=normalizePath("library"),R_LIBS_SITE="")
.libPaths(c(Sys.getenv("R_LIBS_USER"),.Library))
install.packages(c("data.table","RANN","igraph","RSpectra","RcppEigen","matrixStats"),
  lib=Sys.getenv("R_LIBS_USER"),repos="https://cloud.r-project.org")
~~~

With R on PATH, install the file under its distributed attachment name:

~~~text
R CMD INSTALL --preclean --clean --no-multiarch --library=library SupplementarySoftware1_LineaLinkR.tar.gz
~~~

Load the installed package from the same private library:

~~~r
.libPaths(c(normalizePath("library"),.Library))
library(LineaLinkR)
packageVersion("LineaLinkR")
~~~

For numerical reproduction of Table 1 and the uLIPSTIC peak readout, use reproduce_saved_results.R in SupplementarySoftware2_Evaluation.zip, with SourceData.zip as the publication-table comparison. That calculation uses the supplied individual candidate and event records. The evaluation-archive README gives commands for that calculation and for fitting the supplied 1,200-cell, 1,200-gene lrsim_focal_branch dataset from its H5AD. On Windows 10 Pro with R 4.5.3, a Ryzen 7 5700X3D (8 cores, 16 threads) and 128 GiB RAM, installing this source package in a fresh library took 84 s with the compiler and dependencies in place. The single-dataset fit took 62 s and the saved-record calculation took 0.6 s, excluding downloads and archive extraction.
