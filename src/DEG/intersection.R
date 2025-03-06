



# -------------------------
# Load required libraries
# -------------------------
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(devtools)
library(SCopeLoomR)
library(SCENIC)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(grid)
library(data.table)
library(tidyverse)
library(org.Hs.eg.db)
library(annotables)
library(clusterProfiler)
library(msigdbr)
library(reshape2)
library(rlang)

'%!in%' <- function(x, y)!('%in%'(x, y))





# Read the configuration file (adjust the path if necessary)
config_file_path <- "~/src/repository/repository_final/DEG/config_intersection.yml"
config <- yaml::read_yaml(config_file_path)$default




validate_parameters <- function(config) {
  required_params <- c("analysis","discovery_path", 
                       "validation_path", "output_dir")
  
  # Check if all required parameters exist
  missing_params <- setdiff(required_params, names(config))
  if (length(missing_params) > 0) {
    stop(paste("❌ Missing parameters in YAML config:", paste(missing_params, collapse = ", ")))
  }
  
  # Convert `analysis_var` to lowercase
  if (!is.null(config$analysis) && is.character(config$analysis)) {
    config$analysis <- tolower(config$analysis)
  } else {
    stop("❌ Error: 'analysis' must be a character string.")
  }
  
  # Validate `analysis_var`
  if (!config$analysis %in% c("phase", "no_phase")) {
    stop("❌ Error: 'analysis must be either 'phase' or 'no_phase'.")
  }
  
  
  # Validate file paths
  if (!file.exists(config$discovery_path)) {
    stop(paste("❌ Error: Discovery file does not exist:", config$discovery_path))
  }
  if (!file.exists(config$validation_path)) {
    stop(paste("❌ Error: Validation file does not exist:", config$validation_path))
  }
  # If the output directory does not exist, create it recursively.
  if (!dir.exists(config$output_dir)) {
    message("Output directory does not exist. Creating: ", config$output_dir)
    dir.create(config$output_dir, recursive = TRUE)
  }
  
  message("✅ All configuration checks passed successfully.")
  
  return(config)  # ✅ Return the modified config
}



# Safe validation wrapper
safe_validate_parameters <- function(config) {
  tryCatch(
    {
      return(validate_parameters(config))  # ✅ Return the validated config
    },
    error = function(e) {
      message("❌ Caught an error in validation: ", e$message)
      return(NULL)  # Return NULL in case of an error
    }
  )
}

# Run validation and get the updated `config`
config <- safe_validate_parameters(config)


# Assign validated parameters
analysis <- config$analysis
discovery_path <- config$discovery_path
validation_path <- config$validation_path
output_dir <- config$output_dir



comparaciones_discovery_grn <- load(discovery_path)

comparaciones_validation_grn <- load(validation_path)

# -------------------------
# Determine ordered row names:
# Common names first, then unique names.
# -------------------------
names_discovery <- rownames(comparison_discovery)
names_validation <- rownames(comparison_validation)
common_names <- intersect(names_discovery, names_validation)
unique_names <- setdiff(union(names_discovery, names_validation), common_names)
ordered_names <- c(common_names, unique_names)

# -------------------------
# Tidyverse function: Reorder matrix rows and add missing rows
# -------------------------
reorder_add_missing_rows <- function(mat, ordered_names) {
  mat %>%
    as_tibble(rownames = "row") %>%
    complete(row = ordered_names) %>%      # Ensure all ordered_names appear
    pivot_longer(-row, names_to = "col", values_to = "value") %>%
    pivot_wider(names_from = col, values_from = value) %>%
    arrange(factor(row, levels = ordered_names)) %>%
    { 
      m <- as.matrix(select(., -row))
      rownames(m) <- .$row
      m
    }
}

# Reorder both matrices
comparison_discovery <- reorder_add_missing_rows(comparison_discovery, ordered_names)
comparison_validation  <- reorder_add_missing_rows(comparison_validation, ordered_names)

# -------------------------
# Replace NA values with 0 using tidyverse
# -------------------------
# Convert NA values to 0 using tidyverse (we convert to tibble and back to matrix)
fix_na <- function(mat) {
  mat %>%
    as_tibble(rownames = "row") %>%
    dplyr::mutate(across(-row, ~replace_na(.x, 0))) %>%
    pivot_longer(-row, names_to = "col", values_to = "value") %>%
    pivot_wider(names_from = col, values_from = value) %>%
    arrange(match(row, rownames(mat))) %>%
    { as.matrix(select(., -row)) -> m; 
      rownames(m) <- .$row; m }
}


comparison_discovery <- fix_na(comparison_discovery)
comparison_validation  <- fix_na(comparison_validation)

colnames(comparison_validation) = colnames(comparison_discovery)

# -------------------------
# Create the results matrix based on directionality
# -------------------------
# For each cell: if one matrix shows a hit (abs > 0) and the other is 0 → result 1;
# if both are positive → result 2;
# if both are negative → result -2; else 0.
resultado <- left_join(
  comparison_discovery %>% as_tibble(rownames = "row") %>% pivot_longer(-row, names_to = "col", values_to = "disc"),
  comparison_validation %>% as_tibble(rownames = "row") %>% pivot_longer(-row, names_to = "col", values_to = "val"),
  by = c("row", "col")
) %>%
  mutate(result = case_when(
    (abs(disc) > 0 & val == 0) | (disc == 0 & abs(val) > 0) ~ 1,
    (disc > 0 & val > 0) ~ 2,
    (disc < 0 & val < 0) ~ -2,
    TRUE ~ 0
  )) %>%
  select(row, col, result) %>%
  pivot_wider(names_from = col, values_from = result) %>%
  column_to_rownames("row") %>%
  as.matrix()

# Replace any remaining NA with 0
resultado[is.na(resultado)] <- 0



# -------------------------
# Save the results matrix
# -------------------------
save(resultado, file = paste0(output_dir,"DGE_",analysis,".RData"))




