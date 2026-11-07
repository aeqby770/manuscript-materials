#!/usr/bin/env Rscript

# path bootstrap (relative; set LLR_ASSETS or LLR_PIPELINE to override)
.llr_cmd <- commandArgs(trailingOnly = FALSE)
.llr_file <- sub("^--file=", "", grep("^--file=", .llr_cmd, value = TRUE)[1])
if (!length(.llr_file) || is.na(.llr_file) || !nzchar(.llr_file)) {
  stop("Run this script via Rscript so --file= is set")
}
.llr_here <- dirname(normalizePath(.llr_file, winslash = "/", mustWork = TRUE))
WORK_ROOT <- normalizePath(file.path(.llr_here,"../.."),winslash="/",mustWork=TRUE)
Sys.setenv(LLR_WORK_ROOT=WORK_ROOT,LLR_PIPELINE=file.path(WORK_ROOT,"experiments"),LLR_ASSETS=WORK_ROOT)
source(file.path(.llr_here, "../config/llr_bootstrap.R"), local = FALSE)

# TraSig: slingshot prep (Python, no rpy2) + optional TraSig main.py
llr_source_config()
args <- commandArgs(trailingOnly = TRUE)
reuse_native <- "--reuse-native" %in% args
args <- args[args!="--reuse-native"]
ds_name <- args[1]
if (is.na(ds_name)) stop("Usage: Rscript run_trasig.R <dataset_name> [output_dir] [group_column]")
out_dir <- if (length(args) >= 2 && !is.na(args[2])) args[2] else file.path(CONFIG$dir$algorithm_out, "TraSig", ds_name)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

h5ad <- DatasetPath(ds_name)
lr_resource <- switch(CONFIG$dataset_species[[ds_name]],
  mouse = system.file("extdata/Mouse-2020-Shao-LR-pairs.txt", package = "TimeTalk"),
  human = file.path(CONFIG$dir$source_assets, "tools/competitors/original/TraSig/ligand_receptor_lists/ligand_receptor_FANTOM.pickle"),
  stop("TraSig requires an explicit LR resource for this species")
)
stopifnot(nzchar(lr_resource))
binding_path <- file.path(CONFIG$dir$algorithm_out,"TraSig","native_input_binding.csv")
registered <- if (file.exists(binding_path)) data.table::fread(binding_path)[dataset==ds_name] else data.frame()
group_arg <- if (length(args)>=3L) args[3L] else if (nrow(registered)==1L) registered$grouping else ""
if (reuse_native && nrow(registered)==1L) {
  pt <- file.path(CONFIG$dir$algorithm_out,"slingshot",ds_name,"slingshot_lineage_pseudotime.csv")
  files <- c(h5ad,pt,lr_resource)
  expected <- c(registered$h5ad_sha256,registered$slingshot_pt_sha256,registered$fantom_resource_sha256)
  hashes <- vapply(files,digest::digest,character(1),algo="sha256",file=TRUE,USE.NAMES=FALSE)
  stopifnot(identical(group_arg,registered$grouping),identical(unname(hashes),unname(expected)))
}
input_scale <- unname(CONFIG$dataset_lr_input_scale[ds_name])
input_scale[is.na(input_scale)] <- "counts"
status <- sys::exec_wait(
  CONFIG$exe$python,
  args = c(
    file.path(llr_pipeline_dir(), "adapters/run_trasig_prep.py"),
    h5ad, out_dir, ds_name, CONFIG$dir$source_assets, lr_resource, group_arg, input_scale,
    if (reuse_native) "--prepare-only" else character()
  ),
  std_out = file.path(CONFIG$dir$log, sprintf("trasig_prep_%s.log", ds_name)),
  std_err = file.path(CONFIG$dir$log, sprintf("trasig_prep_%s.err", ds_name))
)
if (status!=0L) stop("TraSig preparation or native scoring failed for ",ds_name)
stopifnot(file.exists(file.path(out_dir,"trasig_scores.csv")))
status <- sys::exec_wait(CONFIG$exe$python,args=c("-X","utf8",
  file.path(llr_pipeline_dir(),"adapters/export_trasig_temporal.py"),ds_name,
  "TraSig","complete",out_dir))
if (status!=0L) stop("TraSig curve export failed for ",ds_name)
cat("TraSig scores and temporal curves: ",out_dir,"\n",sep="")
