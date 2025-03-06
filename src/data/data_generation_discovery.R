###############################################################################
# SETUP
###############################################################################
library(tidyverse)    # loads dplyr, tidyr, purrr, readr, tibble, etc.
library(Seurat)
library(Matrix)
library(scCustomize)
# Define “not in” operator (if needed)
`%!in%` <- function(x, y) !('%in%'(x, y))

###############################################################################
# CELL CYCLE GENES
###############################################################################
# Use the built‑in cc.genes list from Seurat
s.genes <- cc.genes$s.genes %>% 
  if_else(. == "MLF1IP", "CENPU", .)
g2m.genes <- cc.genes$g2m.genes %>% 
  if_else(. == "MLF1IP", "CENPU", .)

###############################################################################
# LOAD AND PROCESS INDIVIDUAL SAMPLES
###############################################################################
# List tumor sample directories
sample_dirs <- list.files("~/data/tumor/data/raw_data/", full.names = FALSE)

# Process each sample using purrr::map
seurat_object_list <- sample_dirs %>% 
  set_names() %>%   # names will be the sample directory names
  map(function(sample) {
    base_path <- file.path("~/data/tumor/data/raw_data", sample)
    
    # Load raw count matrix (MTX)
    exp_data <- Matrix::readMM(file.path(base_path, "count_matrix_sparse.mtx"))
    
    # Load gene names and assign them to rows
    genes <- read_delim(file.path(base_path, "count_matrix_genes.tsv"), delim = "\t", col_names = FALSE)
    rownames(exp_data) <- genes$X1
    
    # Load cell barcodes and assign them to columns
    cells <- read_delim(file.path(base_path, "count_matrix_barcodes.tsv"), delim = "\t", col_names = FALSE)
    colnames(exp_data) <- cells$X1
    
    # Create Seurat object
    so <- CreateSeuratObject(
      counts = exp_data,
      min.cells = 1,
      project = sample,
      min.features = 0
    )
    
    # Calculate percent mitochondrial genes
    so[["percent.mito"]] <- PercentageFeatureSet(so, pattern = "^MT-")
    
    # Filtering: keep cells with >200 features and <20% mitochondrial reads
    so <- so %>% subset(subset = nFeature_RNA > 200) %>% subset(subset = percent.mito < 20)
    
    # Normalize data and find variable features
    so <- so %>% NormalizeData() %>% FindVariableFeatures()
    
    # Cell cycle scoring
    so <- CellCycleScoring(so,
                           s.features = s.genes,
                           g2m.features = g2m.genes,
                           set.ident = TRUE)
    
    message("Processed sample: ", sample)
    return(so)
  })

###############################################################################
# SAVE SEPARATE LIST OF SAMPLES
###############################################################################
save(seurat_object_list, file = "~/data/discovery/tumor_data_normalized_list_separate")

###############################################################################
# MERGE AND NORMALIZE COMBINED SEURAT OBJECT
###############################################################################
combined <- Merge_Seurat_List(list_seurat = seurat_object_list) %>% 
  NormalizeData() %>% 
  JoinLayers()  # Join layers if needed by scCustomize

# Clean up cell names by removing underscores
old_cell_names <- Cells(combined)
new_cell_names <- sub("_", "", old_cell_names)

# Check for duplicates after removing underscores
if (anyDuplicated(new_cell_names)) {
  stop("Error: Removing underscores results in duplicate cell names.")
}

# Rename cells in the Seurat object
cell_name_map <- setNames(new_cell_names, old_cell_names)
combined <- RenameCells(object = combined, new.names = cell_name_map)

###############################################################################
# ADD SCSUBTYPE CALLING
###############################################################################
# Load calls from CSV using readr
calls <- read_csv("~/data/discovery/results_mycalls.csv") %>% 
  mutate(my_calls = as_factor(my_calls) %>% 
           fct_recode(Her2_SC = "Her2E_SC")) %>% 
  rename(Cell = 1)   # assuming the first column contains cell names

# Convert combined metadata to a tibble and join calls data by "Cell"
metadata <- combined@meta.data %>% 
  rownames_to_column("Cell") %>% 
  left_join(calls %>% select(Cell, my_calls), by = "Cell") %>% 
  rename(Calls = my_calls)%>%
  drop_na(Calls)

# Save the updated metadata of only cancer epithelial cells to CSV
write_csv(metadata, "~/data/discovery/metadata_discovery.csv")

# Optionally check updated cell names
print(head(Cells(combined)))

###############################################################################
# SAVE COMBINED OBJECT
###############################################################################
save(combined, file = "~/data/discovery/combined_discovery")
