###############################################################################
# SETUP
###############################################################################
library(tidyverse)    # loads dplyr, tidyr, purrr, readr, tibble, stringr, etc.
library(Seurat)
library(Matrix)
library(scCustomize)
library(readxl)

# Define a “not in” operator (if needed)
`%!in%` <- function(x, y) !('%in%'(x, y))

###############################################################################
# LOAD SAMPLE INFORMATION
###############################################################################
info <- read_excel("~/data/validation/sample_info.xlsx")

# List expression data file names from the validation directory
list_files <- list.files("~/data/validation/exp_data/", full.names = FALSE)

# Remove normal samples (first 56)
list_files <- list_files[-(1:56)]

# Remove male samples (those matching patterns "0068-T", "0068-LN", "0178")
list_files <- list_files %>% 
  discard(~ str_detect(.x, "0068-T|0068-LN|0178"))

# Remove involved lymph node samples using info from the Excel file
ln_samples <- info %>% filter(Condition == "Involved LN") %>% pull(`Sample Name`)
ln_patterns <- str_replace(ln_samples, "^[^-]+-(.*)$", "\\1")
pattern <- paste(ln_patterns, collapse = "|")
list_files <- list_files %>% discard(~ str_detect(.x, pattern))

# Function to extract substring before the last dash ("-")
extract_substring <- function(x) {
  last_dash <- str_locate_all(x, "-")[[1]]
  if(nrow(last_dash) > 0) {
    str_sub(x, 1, last_dash[nrow(last_dash), 1] - 1)
  } else {
    x
  }
}
# Clean sample names using map_chr
list_files <- map_chr(list_files, extract_substring)

###############################################################################
# PREPARE CELL CYCLE GENES
###############################################################################
# Use the built-in cc.genes list from Seurat and adjust names with if_else
s.genes <- cc.genes$s.genes %>% if_else(. == "MLF1IP", "CENPU", .)
g2m.genes <- cc.genes$g2m.genes %>% if_else(. == "MLF1IP", "CENPU", .)

###############################################################################
# CREATE SEURAT OBJECTS PER SAMPLE
###############################################################################
# Define a function to process one sample
process_sample <- function(sample) {
  base_path <- file.path("~/data/validation/exp_data", sample)
  mtx_file <- paste0(base_path, "-matrix.mtx.gz")
  barcode_file <- paste0(base_path, "-barcodes.tsv.gz")
  gene_file <- "~/data/validation/gene_features/features.tsv"
  
  # Load raw count matrix using Matrix and read_delim from readr
  exp_data <- Matrix::readMM(mtx_file)
  
  # Load cell barcodes and assign column names
  muestra <- read_delim(barcode_file, delim = "\t", col_names = FALSE)
  colnames(exp_data) <- muestra$X1
  
  # Load gene names and assign row names (first column for IDs, second for gene symbols)
  genes <- read_delim(gene_file, delim = "\t", col_names = FALSE)
  rownames(exp_data) <- genes$X1
  # Convert to standard matrix and switch row names to gene symbols (column 2)
  exp_data <- as.matrix(exp_data)
  rownames(exp_data) <- genes$X2
  
  # Handle duplicated gene symbols: keep the row with the highest sum
  dup_genes <- rownames(exp_data)[duplicated(rownames(exp_data))]
  if (length(dup_genes) > 0) {
    dup_rows <- exp_data[rownames(exp_data) %in% dup_genes, ]
    rest_data <- exp_data[!rownames(exp_data) %in% dup_genes, ]
    
    # For duplicated genes, compute row sums and keep the row with the highest sum per gene
    dup_df <- as.data.frame(dup_rows) %>%
      rownames_to_column("gene") %>%
      mutate(rowsum = rowSums(across(everything(), ~ as.numeric(.x)))) %>%
      group_by(gene) %>%
      arrange(desc(rowsum)) %>%
      slice(1) %>%
      ungroup() %>%
      select(-rowsum)
    dup_mat <- dup_df %>% column_to_rownames("gene") %>% as.matrix()
    
    exp_data <- rbind(rest_data, dup_mat)
  }
  
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
  so <- so %>% subset(subset = nFeature_RNA > 200 & percent.mito < 20)
  
  # Normalize data and identify variable features
  so <- so %>% NormalizeData() %>% FindVariableFeatures()
  
  # Cell cycle scoring
  so <- CellCycleScoring(so, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  
  message("Processed sample: ", sample)
  return(so)
}

# Process samples in pairs (using every other sample as in the original code)
# Here we assume that samples come in pairs, so we use seq(1, length(list_files), by = 2)
seurat_object_list <- map(seq(1, length(list_files), by = 2), function(i) {
  process_sample(list_files[i])
})
# Optionally, you can also extract meta data from each object
phases_epithelial <- map(seurat_object_list, ~ .x@meta.data)

# Name the list using a unique substring from list_files
names(seurat_object_list) <- unique(sub("^[^-]+-(.*)$", "\\1", list_files))

###############################################################################
# SAVE PROCESSED DATA
###############################################################################
save(seurat_object_list, file = "~/data/validation/tumor_data_normalized_list_separate")

# Merge Seurat objects from the list, normalize combined object, and join layers
combined <- Merge_Seurat_List(list_seurat = seurat_object_list, add.cell.ids = names(seurat_object_list)) %>%
  NormalizeData() %>%
  JoinLayers()
metadata <- combined@meta.data %>% rownames_to_column("NAME")
combined@meta.data <- metadata
save(combined, file = "~/data/validation/validation_data_all_normalized")

###############################################################################
# ADD SC SUBTYPE CALLS
###############################################################################
calls <- read_delim("~/data/validation/validation_scsubtype", delim = ",")
metadata <- combined@meta.data %>% 
  left_join(calls %>% select(NAME, SCSubtypeCall), by = "NAME") %>% 
  rename(Calls = SCSubtypeCall)
combined@meta.data <- metadata

###############################################################################
# ADD CELL TYPE FROM SINGLE-R RESULTS
###############################################################################
dir_path <- "~/data/validation/singleR_results/"
files <- list.files(dir_path, full.names = TRUE)
file_info <- file.info(files)
res_files <- files[order(file_info$ctime)]

list_subtypes <- map(res_files, function(f) {
  load(f)
  tibble(NAME = pred.hesc@rownames, celltype_major = pred.hesc$labels)
}) %>% bind_rows() %>%
  mutate(NAME = str_remove(NAME, "_RNA")) %>%
  arrange(NAME) %>%
  distinct(NAME, .keep_all = TRUE)

# Recode Calls and remove calls for cells not in list_subtypes
metadata <- combined@meta.data %>% 
  mutate(Calls = if_else(Calls == "Her2E_SC", "Her2_SC", Calls),
         Calls = if_else(!NAME %in% list_subtypes$NAME, NA_character_, Calls)) %>%
  left_join(list_subtypes, by = "NAME")
combined@meta.data <- metadata

# Filter metadata for Cancer Epithelial cells and rename column for output
metadata <- combined@meta.data %>% 
  filter(celltype_major == "Cancer Epithelial") %>% 
  rename(Cell = NAME)
write_csv(metadata, "~/data/validation/metadata_validation.csv")

###############################################################################
# SAVE COMBINED OBJECT
###############################################################################
save(combined, file = "~/data/validation/validation_data_all_normalized")
