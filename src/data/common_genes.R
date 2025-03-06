# =============================================================================
# Gene Intersection Between Discovery and Validation Datasets
# Author: Miguel Castresana
# Date: 2024-12-18
# Description: This script loads two Seurat objects (discovery and validation),
#              extracts their gene features, finds the intersection of genes,
#              and saves the intersected genes to a CSV file.
# =============================================================================

# ---------------------------
# 1. Setup Environment
# ---------------------------

# Load necessary libraries
library(Seurat)   # For handling Seurat objects
library(dplyr)    # For data manipulation functions like left_join
library(readr)    # For reading and writing data

# ---------------------------
# 2. Define File Paths
# ---------------------------

# Define the file paths for discovery and validation datasets
# Consider using relative paths or environment variables for better portability
discovery_path <- "~/data/discovery/combined_discovery.RData"
validation_path <- "~/data/validation/validation_data_all_normalized.RData"

# Define the output path for the intersected genes
output_genes_path <- "~/data/final_genes_both.csv"

# ---------------------------
# 3. Define Utility Functions
# ---------------------------

# Function to load a Seurat object from an RData file without overwriting existing objects
load_seurat_object <- function(filepath, object_name = "combined") {
  if (!file.exists(filepath)) {
    stop("File not found: ", filepath)
  }
  
  # Create a new environment to load the RData file
  env <- new.env()
  load(filepath, envir = env)
  
  # Check if the expected object exists in the loaded environment
  if (!exists(object_name, envir = env)) {
    stop("Object '", object_name, "' not found in the file: ", filepath)
  }
  
  # Return the loaded Seurat object
  return(get(object_name, envir = env))
}

# ---------------------------
# 4. Load Seurat Objects
# ---------------------------

# Load the discovery Seurat object
seurat_discovery <- load_seurat_object(discovery_path, object_name = "combined")
message("Successfully loaded discovery Seurat object from: ", discovery_path)

# Load the validation Seurat object
seurat_validation <- load_seurat_object(validation_path, object_name = "combined")
message("Successfully loaded validation Seurat object from: ", validation_path)

# ---------------------------
# 5. Extract Gene Features
# ---------------------------

# Function to extract gene features from a Seurat object
extract_gene_features <- function(seurat_obj, assay = "RNA") {
  if (!assay %in% Assays(seurat_obj)) {
    stop("Assay '", assay, "' not found in the Seurat object.")
  }
  
  # Extract gene names from the specified assay
  # Assuming 'Features()' is a custom function; using 'rownames()' for standard extraction
  genes <- rownames(seurat_obj[[assay]])
  
  return(genes)
}

# Extract genes from discovery dataset
genes_discovery <- extract_gene_features(seurat_discovery, assay = "RNA")
message("Extracted ", length(genes_discovery), " genes from discovery dataset.")

# Extract genes from validation dataset
genes_validation <- extract_gene_features(seurat_validation, assay = "RNA")
message("Extracted ", length(genes_validation), " genes from validation dataset.")

# ---------------------------
# 6. Find Intersection of Genes
# ---------------------------

# Find intersecting genes between discovery and validation datasets
intersect_genes <- intersect(genes_discovery, genes_validation)
message("Found ", length(intersect_genes), " intersecting genes between discovery and validation datasets.")

# Convert the intersecting genes to a data frame
genes_df <- data.frame(genes = intersect_genes, stringsAsFactors = FALSE)

# ---------------------------
# 7. Save Intersected Genes to CSV
# ---------------------------

# Save the intersected genes to a CSV file
# Using 'write_csv' from the 'readr' package for better performance and handling
write_csv(genes_df, output_genes_path)
message("Successfully saved intersected genes to: ", output_genes_path)

# =============================================================================
# End of Script
# =============================================================================
