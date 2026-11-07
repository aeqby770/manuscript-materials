# LineaLinkR

LineaLinkR evaluates ligand-receptor expression changes across related receiver states. It compares a focal daughter with observed sisters within each biological sample and records support for local and late observations.

## Installation

The package requires R >= 4.2, a C++17 compiler and Rust/Cargo >= 1.80. On Windows, use the Rtools version matching R and the GNU Rust toolchain (x86_64-pc-windows-gnu).

~~~r
install.packages("remotes")
remotes::install_github(
  "aeqby770/manuscript-materials",
  subdir = "source/LineaLinkR"
)
library(LineaLinkR)
packageVersion("LineaLinkR")
help(package = "LineaLinkR")
~~~

## Repository contents

- source/LineaLinkR/ contains the R package, native Rust source, documentation and package tests.
- analysis_scripts/experiments/ contains the scientific input-preparation, fitting, evaluation and comparison code.
- data_sources/public_input_manifest.csv identifies the fixed inputs and their original study accessions.
- results/ contains the displayed main and supplementary result tables.

## Reproduce the main comparison

The paper's SupplementarySoftware2_Evaluation.zip supplies candidate scores, reference labels, target records and fixed simulation matrices. Extract it to a directory named evaluation. Keep SourceData.zip from the paper's attachments alongside it. With data.table and PRROC installed, run:

~~~text
Rscript --vanilla analysis_scripts/experiments/evaluation/reproduce_saved_results.R evaluation SourceData.zip
~~~

The command calculates the twelve rows of Table 1 from individual records and checks the measured contact-peak profiles. Outputs are written under evaluation/reproduced/.

For refitting, experiments/run/run_linealinkr.R is the main method entry in an assembled analysis workspace. The evaluation archive README provides the expression-input layout and fitting commands. External-method predictions retain their task-specific candidate and coverage definitions.

The package is distributed under the license in source/LineaLinkR/DESCRIPTION and source/LineaLinkR/LICENSE.
