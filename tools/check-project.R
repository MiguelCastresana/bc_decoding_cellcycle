#!/usr/bin/env Rscript

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || is.na(x)) y else x
}

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg) > 0) {
  project_root <- normalizePath(file.path(dirname(sub("^--file=", "", file_arg[1])), ".."), mustWork = TRUE)
} else {
  project_root <- normalizePath(getwd(), mustWork = TRUE)
}

required_files <- c(
  "README.md",
  "install_packages.R",
  "python_environment.yaml",
  "src/data/data_generation_discovery.R",
  "src/data/data_generation_validation.R",
  "src/DEG/DGE_analysis.R",
  "src/GRN/GRN_analysis.ipynb"
)

missing <- required_files[!file.exists(file.path(project_root, required_files))]
if (length(missing) > 0) {
  stop("Missing required files:\n", paste(missing, collapse = "\n"), call. = FALSE)
}

r_files <- list.files(file.path(project_root, "src"), pattern = "\\.R$", recursive = TRUE, full.names = TRUE)
r_files <- c(file.path(project_root, "install_packages.R"), r_files)

parse_one <- function(path) {
  tryCatch(
    {
      parse(path)
      TRUE
    },
    error = function(err) {
      message("Parse failed: ", path)
      message(conditionMessage(err))
      FALSE
    }
  )
}

ok <- vapply(r_files, parse_one, logical(1))
if (!all(ok)) {
  stop("One or more R files failed to parse.", call. = FALSE)
}

message("Project structure and R syntax checks passed.")
