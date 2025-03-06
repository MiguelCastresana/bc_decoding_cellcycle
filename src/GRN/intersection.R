



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
config_file_path <- "~/src/repository/repository_final/GRN/config_intersection.yml"
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


# Extract common row and column names
common_rownames <- union(rownames(comparaciones_discovery_grn), rownames(comparaciones_validation_grn))
common_colnames <- union(colnames(comparaciones_discovery_grn), colnames(comparaciones_validation_grn))

# A tidyverse function to reformat a matrix to have given rows and columns
reformat_matrix <- function(mat, common_rows, common_cols) {
  mat %>%
    as_tibble(rownames = "row") %>%                     # Convert to tibble with rownames in column "row"
    pivot_longer(-row, names_to = "col", values_to = "value") %>% 
    complete(row = common_rows, col = common_cols, 
             fill = list(value = NA)) %>%               # Ensure all rows and columns are present
    pivot_wider(names_from = col, values_from = value) %>% 
    arrange(match(row, common_rows)) -> df_wide         # Order rows based on common_rows
  
  mat_out <- df_wide %>% select(-row) %>% as.matrix()
  rownames(mat_out) <- df_wide$row
  return(mat_out)
}

# Reformat both matrices
comparaciones_discovery_grn <- reformat_matrix(comparaciones_discovery_grn, common_rownames, common_colnames)
comparaciones_validation_grn <- reformat_matrix(comparaciones_validation_grn, common_rownames, common_colnames)

# Convert NA values to 0 using tidyverse (we convert to tibble and back to matrix)
fix_na <- function(mat) {
  mat %>%
    as_tibble(rownames = "row") %>%
    mutate(across(-row, ~replace_na(.x, 0))) %>%
    pivot_longer(-row, names_to = "col", values_to = "value") %>%
    pivot_wider(names_from = col, values_from = value) %>%
    arrange(match(row, rownames(mat))) %>%
    { as.matrix(select(., -row)) -> m; 
      rownames(m) <- .$row; m }
}

comparaciones_discovery_grn <- fix_na(comparaciones_discovery_grn)
comparaciones_validation_grn <- fix_na(comparaciones_validation_grn)

resultado <- left_join(
  comparaciones_discovery_grn %>% 
    as_tibble(rownames = "row") %>% 
    pivot_longer(-row, names_to = "col", values_to = "disc"),
  comparaciones_validation_grn %>% 
    as_tibble(rownames = "row") %>% 
    pivot_longer(-row, names_to = "col", values_to = "val"),
  by = c("row", "col")
) %>% 
  mutate(result = case_when(
    disc == 1 & val == 1 ~ 2,   # both have a hit
    disc == 1 | val == 1 ~ 1,   # only one has a hit
    TRUE ~ 0
  )) %>% 
  select(row, col, result) %>% 
  pivot_wider(names_from = col, values_from = result) %>% 
  column_to_rownames("row") %>% 
  as.matrix()



# Save the result

save(resultado, file = paste0(output_dir,"GRN_",analysis,".RData"))

