# Shared path bootstrap for LineaLinkR benchmark pipeline scripts.
# Prefer this over hardcoded G:/ or F:/ absolute paths.
#
# Usage (Rscript entry points):
#   source(file.path(dirname(llr_this_file()), "../config/llr_bootstrap.R"), local = TRUE)
#   # or from adapters/:
#   source(file.path(dirname(llr_this_file()), "../config/llr_bootstrap.R"), local = TRUE)
#   llr_source_config()
#   llr_source_rel("../adapters/_mem_guard.R")  # example

# Capture the bootstrap path while source() still exposes its `ofile` frame.
# Calling sys.frames() later loses that frame and breaks location discovery.
.llr_bootstrap_frames <- sys.frames()
.llr_bootstrap_ofiles <- vapply(.llr_bootstrap_frames, function(fr) {
  if (exists("ofile", envir = fr, inherits = FALSE)) {
    of <- get("ofile", envir = fr, inherits = FALSE)
    if (is.character(of) && length(of) == 1L && nzchar(of)) return(of)
  }
  ""
}, character(1L))
.llr_bootstrap_file <- tail(.llr_bootstrap_ofiles[nzchar(.llr_bootstrap_ofiles)], 1L)
.llr_bootstrap_file <- normalizePath(.llr_bootstrap_file, winslash = "/", mustWork = TRUE)

llr_this_file <- function() {
  cmd <- commandArgs(trailingOnly = FALSE)
  hit <- grep("^--file=", cmd, value = TRUE)
  if (length(hit) < 1L) {
    stop("llr_this_file(): run via Rscript so --file= is available", call. = FALSE)
  }
  normalizePath(sub("^--file=", "", hit[[1L]]), winslash = "/", mustWork = TRUE)
}

llr_pipeline_dir <- function() {
  env <- Sys.getenv("LLR_PIPELINE", unset = "")
  if (nzchar(env)) {
    return(normalizePath(env, winslash = "/", mustWork = TRUE))
  }
  # this file lives in pipeline/config/
  normalizePath(file.path(dirname(llr_this_file_bootstrap()), ".."), winslash = "/", mustWork = TRUE)
}

# When bootstrap itself is sourced, capture its path via sys.frames ofile if present.
llr_this_file_bootstrap <- function() {
  if (exists(".llr_bootstrap_file", inherits = TRUE) &&
      is.character(.llr_bootstrap_file) && length(.llr_bootstrap_file) == 1L &&
      nzchar(.llr_bootstrap_file)) {
    return(.llr_bootstrap_file)
  }
  for (i in rev(seq_len(sys.nframe()))) {
    fr <- sys.frame(i)
    if (exists("ofile", envir = fr, inherits = FALSE)) {
      of <- get("ofile", envir = fr, inherits = FALSE)
      if (is.character(of) && length(of) == 1L && nzchar(of)) {
        return(normalizePath(of, winslash = "/", mustWork = TRUE))
      }
    }
  }
  # Fallback: locate from the isolated work root.
  work <- Sys.getenv("LLR064_WORK_ROOT", unset = "")
  if (nzchar(work)) {
    return(normalizePath(file.path(work, "pipeline", "config", "llr_bootstrap.R"),
                         winslash = "/", mustWork = TRUE))
  }
  env <- Sys.getenv("LLR_ASSETS", unset = "")
  if (nzchar(env)) {
    return(normalizePath(file.path(env, "pipeline", "config", "llr_bootstrap.R"),
                         winslash = "/", mustWork = TRUE))
  }
  stop("Cannot locate llr_bootstrap.R; set LLR_ASSETS or LLR_PIPELINE", call. = FALSE)
}

llr_assets_dir <- function() {
  work <- Sys.getenv("LLR064_WORK_ROOT", unset = "")
  if (nzchar(work)) {
    return(normalizePath(work, winslash = "/", mustWork = TRUE))
  }
  env <- Sys.getenv("LLR_ASSETS", unset = "")
  if (nzchar(env)) {
    return(normalizePath(env, winslash = "/", mustWork = TRUE))
  }
  normalizePath(file.path(llr_pipeline_dir(), ".."), winslash = "/", mustWork = TRUE)
}

llr_source_config <- function() {
  cfg <- file.path(llr_pipeline_dir(), "config", "00_config.R")
  source(cfg, local = FALSE)
  invisible(TRUE)
}

llr_source_rel <- function(rel) {
  path <- normalizePath(file.path(llr_pipeline_dir(), rel), winslash = "/", mustWork = TRUE)
  source(path, local = FALSE)
  invisible(TRUE)
}
