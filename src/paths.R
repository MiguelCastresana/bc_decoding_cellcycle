# Shared paths for breast cancer cell-cycle analysis scripts.

script_args <- commandArgs(trailingOnly = FALSE)
script_file_arg <- grep("^--file=", script_args, value = TRUE)
script_dir <- if (length(script_file_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", script_file_arg[1]), mustWork = TRUE))
} else {
  normalizePath(getwd(), mustWork = TRUE)
}

repo_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = FALSE)
if (!file.exists(file.path(repo_root, "README.md"))) {
  repo_root <- normalizePath(getwd(), mustWork = TRUE)
}

data_root <- Sys.getenv("BC_CELL_CYCLE_DATA_DIR", unset = file.path(repo_root, "data"))
results_root <- Sys.getenv("BC_CELL_CYCLE_RESULTS_DIR", unset = file.path(repo_root, "results"))
figures_root <- Sys.getenv("BC_CELL_CYCLE_FIGURES_DIR", unset = file.path(repo_root, "figures"))

data_file <- function(...) {
  file.path(data_root, ...)
}

results_file <- function(...) {
  path <- file.path(results_root, ...)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  path
}

figures_file <- function(...) {
  path <- file.path(figures_root, ...)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  path
}

repo_file <- function(...) {
  file.path(repo_root, ...)
}

resolve_project_path <- function(path) {
  vapply(path, function(value) {
    if (is.na(value) || value == "") {
      return(value)
    }

    expanded <- path.expand(value)
    is_absolute <- grepl("^([A-Za-z]:)?[\\/]", expanded)
    if (is_absolute) {
      return(normalizePath(expanded, mustWork = FALSE))
    }

    file.path(repo_root, expanded)
  }, character(1), USE.NAMES = FALSE)
}

resolve_config_paths <- function(config, path_names) {
  for (path_name in intersect(path_names, names(config))) {
    config[[path_name]] <- resolve_project_path(config[[path_name]])
  }
  config
}
