###############################################################################
# SingleR Annotation Pipeline for Single-Cell RNA-Seq Data
# Author: Miguel Castresna
# Date: 2024-12-18
# Description: This script performs annotation of validation single-cell 
#              RNA-seq data using the SingleR package. It loads training 
#              and test datasets, preprocesses them, runs SingleR in chunks,
#              and combines the results.
###############################################################################

library(tidyverse)            # Loads dplyr, tidyr, purrr, readr, tibble, stringr, etc.
library(Seurat)
library(SingleR)
library(SummarizedExperiment)
library(scuttle)

# ---------------------------
# 2. Define File Paths
# ---------------------------
discovery_seurat_path    <- "~/data/discovery/combined_discovery"
discovery_metadata_path  <- "~/data/discovery/Whole_miniatlas_meta.csv"
validation_seurat_path   <- "~/data/validation/validation_data_all_normalized"
results_dir              <- "~/data/validation/singleR_results_new/"
final_genes_path         <- "~/data/final_genes_both.csv"

if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
  message("Created results directory: ", results_dir)
}

# ---------------------------
# 3. Define Utility Functions
# ---------------------------
# Tidyverse function to load a Seurat object from an RData file
load_seurat_object <- function(filepath, object_name = "combined") {
  if (!file.exists(filepath)) stop("File not found: ", filepath)
  env <- new.env()
  load(filepath, envir = env)
  if (!exists(object_name, envir = env)) stop("Object '", object_name, "' not found in file: ", filepath)
  get(object_name, envir = env)
}

# ---------------------------
# 4. Load and Preprocess Training Data (Discovery)
# ---------------------------
seurat_discovery <- load_seurat_object(discovery_seurat_path, object_name = "combined")
message("Successfully loaded discovery Seurat object from: ", discovery_seurat_path)

their_metadata <- read_delim(discovery_metadata_path, delim = ",", col_names = TRUE) %>%
  tibble::column_to_rownames("NAME") %>%
  rownames_to_column(var = "NAME")

metadata_discovery <- seurat_discovery@meta.data %>% 
  rownames_to_column("NAME") %>% 
  select(NAME, orig.ident) %>% 
  left_join(their_metadata, by = "NAME") %>% 
  column_to_rownames("NAME")
seurat_discovery@meta.data <- metadata_discovery

counts_discovery <- GetAssayData(seurat_discovery, slot = "counts")

remove(seurat_discovery)
gc()

se_train <- SummarizedExperiment(assays = list(counts = counts_discovery), colData = metadata_discovery)
se_train <- as(se_train, "SingleCellExperiment")
se_train <- logNormCounts(se_train)
message("Training data preprocessed and normalized.")
# ---------------------------
# 5. Load and Preprocess Test Data (Validation)
# ---------------------------
seurat_validation <- load_seurat_object(validation_seurat_path, object_name = "combined")
message("Successfully loaded validation Seurat object from: ", validation_seurat_path)

metadata_validation <- seurat_validation@meta.data
rownames(metadata_validation) = metadata_validation$NAME
counts_validation <- GetAssayData(seurat_validation, slot = "counts")
colnames(counts_validation) <- rownames(metadata_validation)

se_test <- SummarizedExperiment(assays = list(counts = counts_validation), colData = metadata_validation)
se_test <- as(se_test, "SingleCellExperiment")
se_test <- logNormCounts(se_test)
message("Test data preprocessed and normalized.")

# ---------------------------
# 6. Run SingleR in Chunks on Test Data
# ---------------------------
chunk_size <- 100
total_cells <- ncol(se_test)
chunk_starts <- seq(1, total_cells, by = chunk_size)
chunk_ends <- c(chunk_starts[-1] - 1, total_cells)

remove(seurat_validation,counts_discovery,counts_validation)
gc()

singleR_results <- map2(chunk_starts, chunk_ends, function(start_idx, end_idx) {
  se_test_chunk <- se_test[, start_idx:end_idx]
  pred <- SingleR(
    test = se_test_chunk,
    ref = se_train,
    assay.type.test = 1,
    labels = se_train@colData$celltype_major
  )
  message("Processed cells: ", start_idx, " to ", end_idx)
  # Save each chunk's result as an RData file
  result_filename <- paste0("result_", start_idx, "_", end_idx, ".RData")
  save(pred, file = file.path(results_dir, result_filename))
  pred
})
message("Completed running SingleR on all chunks.")

# ---------------------------
# 7. Combine SingleR Results
# ---------------------------
result_files <- list.files(results_dir, pattern = "^result_.*\\.RData$", full.names = TRUE)
result_files <- result_files[order(file.info(result_files)$ctime)]

singleR_predictions <- map(result_files, function(f) {
  load(f)  # loads object 'pred'
  if (!exists("pred")) {
    warning("Object 'pred' not found in file: ", f)
    return(NULL)
  }
  pred$labels
}) %>% compact()

combined_predictions <- reduce(singleR_predictions, c)
prediction_table <- table(combined_predictions)
message("Distribution of SingleR Predicted Labels:")
print(prediction_table)

# =============================================================================
# End of Script
# =============================================================================
