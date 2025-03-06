# =============================================================================
# Single-Cell RNA-Seq Data Processing and Subtype Classification Pipeline
# Author: Miguel Castresana
# Date: 2024-12-18
# Description: This script processes single-cell RNA-seq data using Seurat,
#              merges metadata, calculates SC-subtype scores, and outputs
#              classification results, all using tidyverse.
# =============================================================================

# ---------------------------
# 1. Setup Environment
# ---------------------------
library(tidyverse)
library(Seurat)
library(Matrix)

# Function to install missing packages using tidyverse style
install_if_missing <- function(pkgs) {
  pkgs %>% 
    setdiff(rownames(installed.packages())) %>%
    walk(~ install.packages(.x, dependencies = TRUE))
}

required_packages <- c("Seurat", "Matrix", "dplyr")
install_if_missing(required_packages)

# Define the "not in" operator for convenience
`%!in%` <- function(x, y) !(x %in% y)

# Define data directory
data_dir <- "~/data/discovery"

# ---------------------------
# 2. Load Gene Expression Data
# ---------------------------
load_gene_matrix <- function(filepath) {
  if (!file.exists(filepath)) stop("Gene matrix file not found: ", filepath)
  gene_matrix <- Matrix::readMM(filepath)
  message("Loaded gene expression matrix from: ", filepath)
  return(gene_matrix)
}

# Define file paths
gene_matrix_path <- file.path(data_dir, "gene_sorted-matrix.mtx")
features_path <- file.path(data_dir, "features.tsv.gz")
barcodes_path <- file.path(data_dir, "barcodes.tsv.gz")

# Load gene expression matrix
gene_matrix <- load_gene_matrix(gene_matrix_path)
message("Gene matrix dimensions: ", paste(dim(gene_matrix), collapse = " x "))

# Load features and barcodes using tidyverse functions
features <- read_tsv(features_path, col_names = FALSE)
barcodes_cells <- read_tsv(barcodes_path, col_names = FALSE)
message("Loaded features and barcodes.")

# Assign row and column names
rownames(gene_matrix) <- features$X1
colnames(gene_matrix) <- barcodes_cells$X1

# Create a Seurat object
seurat_object <- CreateSeuratObject(counts = gene_matrix)
seurat_object[["RNA"]] <- CreateAssayObject(data = as.matrix(gene_matrix))
message("Created Seurat object with RNA assay.")

# ---------------------------
# 3. Load and Merge Metadata
# ---------------------------
load_metadata <- function(filepath) {
  if (!file.exists(filepath)) stop("Metadata file not found: ", filepath)
  metadata <- read_csv(filepath)
  metadata <- metadata %>% column_to_rownames(var = "NAME")
  message("Loaded metadata from: ", filepath)
  return(metadata)
}

metadata_path <- file.path(data_dir, "Whole_miniatlas_meta.csv")
their_metadata <- load_metadata(metadata_path)
message("Metadata dimensions: ", paste(dim(their_metadata), collapse = " x "))

# Merge metadata with Seurat object
seurat_object@meta.data <- their_metadata
message("Merged metadata with Seurat object.")

# Filter to Cancer Epithelial cells using tidyverse
cancer_epithelial_metadata <- their_metadata %>% filter(celltype_major == "Cancer Epithelial")
seurat_object <- subset(seurat_object, cells = cancer_epithelial_metadata %>% rownames())
message("Subset Seurat object to Cancer Epithelial cells. Dimensions: ",
        paste(dim(seurat_object@assays$RNA@data), collapse = " x "))

# ---------------------------
# 4. Check Patient Samples
# ---------------------------
unique_patients <- seurat_object@meta.data %>% pull(Patient) %>% unique()
message("Unique patients identified: ", length(unique_patients))
print(unique_patients)

# ---------------------------
# 5. Load SC-Subtype Signatures
# ---------------------------
load_sc_subtype_signatures <- function(filepath) {
  if (!file.exists(filepath)) stop("SC-subtype signatures file not found: ", filepath)
  sigdat <- read_csv(filepath, col_types = cols())
  
  # Replace periods with dashes in all columns
  sigdat = as.data.frame(apply(sigdat,2,function(x) sub("\\.", "-", x)))
  
  sigdat[sigdat == "AP000769-1"] <- "AP000769.1"
  
  # Extract unique genes from selected subtype columns and remove blanks
  signature_genes <- sigdat %>%
    select(Basal_SC, Her2E_SC, LumA_SC, LumB_SC) %>%
    pivot_longer(everything(), values_to = "gene") %>%
    pull(gene) %>%
    discard(~ is.na(.x) || .x == "") %>%
    unique()
  
  message("Loaded and processed SC-subtype signatures from: ", filepath)
  return(list(sigdat = sigdat, signature_genes = signature_genes))
}

sigdat_path <- file.path(data_dir, "NatGen_Supplementary_table_S4.csv")
sc_subtype_data <- load_sc_subtype_signatures(sigdat_path)
sigdat <- sc_subtype_data$sigdat
signature_genes <- sc_subtype_data$signature_genes

# Check for missing genes in the dataset using tidyverse set operations
missing_genes <- signature_genes[signature_genes %!in% Features(seurat_object, assay = "RNA")]
if (length(missing_genes) == 0) {
  message("No missing genes in the dataset.")
} else {
  warning("Missing genes: ", paste(missing_genes, collapse = ", "))
}

# ---------------------------
# 6. Calculate SC-Subtype Scores
# ---------------------------
# Scale data for the signature genes
seurat_object <- ScaleData(seurat_object, features = signature_genes)
message("Scaled data for signature genes.")

# Extract scaled data (for the signature genes)
scaled_data <- GetAssayData(seurat_object, assay = "RNA", slot = "scale.data")
scaled_data <- scaled_data[signature_genes, ]
message("Extracted scaled data for signature genes.")

# Calculate mean SC-subtype scores using tidyverse and purrr
# Convert selected signature columns to a list (each element is a vector of genes)
sig_list <- sigdat %>% 
  select(Basal_SC, Her2E_SC, LumA_SC, LumB_SC) %>% 
  as.list()

# For each subtype, compute the column means over signature genes present
subtype_scores <- sig_list %>% map(~ {
  genes <- .x[.x != ""]
  genes_present <- intersect(genes, rownames(scaled_data))
  if (length(genes_present) > 0) {
    colMeans(scaled_data[genes_present, , drop = FALSE], na.rm = TRUE)
  } else {
    rep(0, ncol(scaled_data))
  }
})

# Convert the list of subtype scores into a matrix (rows: subtypes, columns: cells)
score_matrix <- do.call(rbind, subtype_scores)

# Remove rows with all zero scores and round numerical values
final_scores <- score_matrix[rowSums(score_matrix, na.rm = TRUE) != 0, ]
final_scores <- as.data.frame(final_scores) %>%
  mutate(across(where(is.numeric), ~ round(., 4)))
final_scores_matrix <- as.matrix(final_scores)
message("Calculated SC-subtype scores.")

# Define a scaling function to center scores
center_sweep <- function(x, row.w = rep(1, nrow(x)) / nrow(x)) {
  get_average <- function(v) sum(v * row.w) / sum(row.w)
  average <- apply(x, 2, get_average)
  sweep(x, 2, average)
}

# Obtain the highest call per cell
finalmt <- as.data.frame(t(final_scores_matrix))
scaled_scores_df <- center_sweep(finalmt)
Finalnames <- colnames(scaled_scores_df)[max.col(scaled_scores_df, ties.method = "first")]

# Convert results to a tibble and add SCSubtypeCall
scaled_scores_df <- scaled_scores_df %>% 
  as_tibble(rownames = "Cell") %>%
  mutate(SCSubtypeCall = Finalnames)
message("Assigned SC-subtype calls.")

# ---------------------------
# 7. Output Results
# ---------------------------
pred_table <- table(scaled_scores_df$SCSubtypeCall)
print("Distribution of SC-Subtype Calls:")
print(pred_table)

message("SC-subtype classification completed.")

# =============================================================================
# End of Script
# =============================================================================
