# =============================================================================
# Single-Cell RNA-Seq Data Processing and Annotation Pipeline
# Author: Miguel Castresana
# Date: 2024-12-18
# Description: This script processes single-cell RNA-seq data for tumor samples,
#              performs normalization, splits data into training and test sets,
#              runs SingleR for cell type prediction, and evaluates the results,
#              using tidyverse principles throughout.
# =============================================================================

# ---------------------------
# 1. Setup Environment
# ---------------------------

library(tidyverse)
library(caret)
library(Seurat)
library(SingleR)
library(SummarizedExperiment)
library(scuttle)
library(scCustomize)
library(readxl)

# Function to install missing packages using purrr
install_if_missing <- function(pkgs) {
  pkgs %>%
    setdiff(rownames(installed.packages())) %>%
    walk(~ install.packages(.x, dependencies = TRUE))
}
required_packages <- c("dplyr", "caret", "Seurat", "SingleR",
                       "SummarizedExperiment", "scuttle", "scCustomize", "readxl")
install_if_missing(required_packages)

# ---------------------------
# 2. Load Data
# ---------------------------

# Function to load a Seurat object from an RData file into a new environment
load_seurat_object <- function(filepath, object_name = "combined") {
  if (!file.exists(filepath)) stop("File not found: ", filepath)
  env <- new.env()
  load(filepath, envir = env)
  if (!exists(object_name, envir = env)) stop("Object '", object_name, "' not found in file: ", filepath)
  get(object_name, envir = env)
}

# Load integrated Seurat object for discovery
discovery_seurat_path <- "~/data/discovery/combined_discovery"
seurat_discovery <- load_seurat_object(discovery_seurat_path, object_name = "combined")
message("Successfully loaded discovery Seurat object from: ", discovery_seurat_path)

# Load metadata and set rownames using tibble functions
their_metadata_path <- "~/data/discovery/Whole_miniatlas_meta.csv"
if (!file.exists(their_metadata_path)) stop("Metadata file not found: ", their_metadata_path)
their_metadata <- read_delim(their_metadata_path, delim = ",", col_types = cols()) %>%
  column_to_rownames(var = "NAME")

# Merge metadata with Seurat object meta.data using tidyverse piping
metadata_discovery <- seurat_discovery@meta.data %>%
  rownames_to_column(var = "NAME") %>%
  select(NAME, orig.ident) %>%
  left_join(their_metadata %>% rownames_to_column(var = "NAME"), by = "NAME") %>%
  column_to_rownames(var = "NAME")
seurat_discovery@meta.data <- metadata_discovery
message("Merged discovery metadata with Seurat object.")

# ---------------------------
# 3. Create SummarizedExperiment and Normalize Counts
# ---------------------------

counts_discovery <- GetAssayData(object = seurat_discovery, assay = "RNA", slot = "counts")
se_train <- SummarizedExperiment(assays = list(counts = counts_discovery),
                                 colData = metadata_discovery) %>%
  as("SingleCellExperiment")
assayNames(se_train)[1] <- "counts"
se_train <- logNormCounts(se_train)
message("Training data converted to SingleCellExperiment and normalized.")

# ---------------------------
# 4. Load Additional Data
# ---------------------------

# Load supplementary tables from Excel using readxl and purrr
supplementary_tables_path <- "~/data/discovery/supplementary_tables.xlsx"
if (!file.exists(supplementary_tables_path)) stop("Supplementary tables file not found: ", supplementary_tables_path)
sheet_names <- excel_sheets(supplementary_tables_path)
mylist <- map(sheet_names, ~ read_excel(supplementary_tables_path, sheet = .x))
tumor_labels <- mylist[[3]]
message("Loaded supplementary tables from: ", supplementary_tables_path)

# Load normalized list of separate Seurat objects (assumes object named 'seurat_object_list')
normalized_seurat_list_path <- "~/data/discovery/discovery_normalized_list_separate"
if (!file.exists(normalized_seurat_list_path)) stop("Normalized Seurat object list file not found: ", normalized_seurat_list_path)
load(normalized_seurat_list_path)
if (is.null(seurat_object_list)) stop("Seurat object list not found in file: ", normalized_seurat_list_path)
message("Loaded normalized Seurat object list from: ", normalized_seurat_list_path)

# Filter metadata for Cancer Epithelial cells using tidyverse
their_metadata <- their_metadata %>% filter(celltype_major == "Cancer Epithelial")
message("Filtered metadata for Cancer Epithelial cells.")

# ---------------------------
# 5. Prepare Training Set
# ---------------------------

# Clean tumor_labels data frame
nombres <- tumor_labels %>% dplyr::slice(3) %>% unlist() %>% as.character()
tumor_labels <- tumor_labels %>% dplyr::slice(-1:-3)
colnames(tumor_labels) <- nombres
training_set <- tumor_labels %>% as_tibble()
message("Prepared training set from tumor labels.")

# Identify unique subtypes
subtypes <- training_set %>% pull(3) %>% unique()
message("Identified unique subtypes: ", paste(subtypes, collapse = ", "))

# Filter seurat_object_list to include only samples present in training_set (using first column)
seurat_object_list <- seurat_object_list[names(seurat_object_list) %in% training_set[[1]]]
message("Filtered Seurat object list to include relevant samples.")

# Subset each Seurat object to include only cells present in their_metadata (rownames)
seurat_object_list <- map(seurat_object_list, ~ subset(.x, cells = rownames(their_metadata)))
message("Subset Seurat objects to include only Cancer Epithelial cells.")

# Split seurat_object_list into subtypes based on training_set information
lumA <- seurat_object_list[names(seurat_object_list) %in%
                             (training_set %>% filter(.[[3]] == "LumA") %>% pull(1))]
lumB <- seurat_object_list[names(seurat_object_list) %in%
                             (training_set %>% filter(.[[3]] == "LumB") %>% pull(1))]
her2 <- seurat_object_list[names(seurat_object_list) %in%
                             (training_set %>% filter(.[[3]] == "Her2") %>% pull(1))]
basal <- seurat_object_list[names(seurat_object_list) %in%
                              (training_set %>% filter(.[[3]] == "Basal") %>% pull(1))]
message("Split Seurat objects into subtypes: LumA, LumB, Her2, Basal.")

# ---------------------------
# 6. Define Function to Split Data 70-30 using Tidyverse
# ---------------------------

split_data_7030 <- function(values, groups) {
  df <- tibble(Group = groups, Value = as.numeric(values)) %>%
    arrange(Group, desc(Value)) %>%
    mutate(cum_value = cumsum(Value))
  
  total <- sum(df$Value)
  target_train <- total * 0.7
  
  df <- df %>% mutate(Set = if_else(cum_value <= target_train, "Train", "Test"))
  
  list(TrainSet = df %>% filter(Set == "Train") %>% select(-cum_value, -Set),
       TestSet  = df %>% filter(Set == "Test") %>% select(-cum_value, -Set))
}

message("Defined function to split data into 70% training and 30% testing using tidyverse.")

# ---------------------------
# 7. Calculate Cells per Sample and Split by Subtype
# ---------------------------

luma_g  <- map_int(lumA, ~ ncol(.x))
lumb_g  <- map_int(lumB, ~ ncol(.x))
her2_g  <- map_int(her2, ~ ncol(.x))
basal_g <- map_int(basal, ~ ncol(.x))

lumA_split  <- split_data_7030(luma_g, names(luma_g))
lumB_split  <- split_data_7030(lumb_g, names(lumb_g))
her2_split  <- split_data_7030(her2_g, names(her2_g))
basal_split <- split_data_7030(basal_g, names(basal_g))

message("Split data into training and testing sets for each subtype.")

# ---------------------------
# 8. Manual Sample Selection for Train/Test
# ---------------------------

basal_train <- c("CID4465", "CID4495", "CID4523", "CID3963", "CID4515")
basal_test  <- c("CID4513", "CID44971")
lumA_train  <- c("CID4290A", "CID4463", "CID4530N")
lumA_test   <- c("CID3941", "CID4067")
her2_train  <- c("CID3921", "CID44991")
her2_test   <- c("CID4066", "CID45171")
lumB_train  <- c("CID4535")
lumB_test   <- c("CID3948", "CID4461")

normal_samples <- training_set %>% filter(.[[3]] == "Normal") %>% pull(`Tumour ID`)

train_samples <- c(basal_train, lumA_train, lumB_train, her2_train, normal_samples)
test_samples  <- c(basal_test, lumA_test, lumB_test, her2_test)

discovery_training <- seurat_object_list[train_samples]
message("Selected training and testing samples.")

# ---------------------------
# 9. Clean Up Environment
# ---------------------------

rm(seurat_object_list, tumor_labels, mylist)
gc()
message("Cleaned up environment and freed memory.")

# ---------------------------
# 10. Prepare Train and Test Sets from SingleCellExperiment
# ---------------------------

# Convert colData to tibble for easier filtering (assuming a 'Patient' column exists)
se_train_df <- as_tibble(colData(se_train), rownames = "Cell")
pos_train <- se_train_df %>% filter(Patient %in% train_samples) %>% pull(Cell)
pos_test  <- se_train_df %>% filter(Patient %in% test_samples) %>% pull(Cell)

se_train_subset <- se_train[, pos_train]
se_test_subset  <- se_train[, pos_test]

# (Re)apply log-normalization if needed
se_train_subset <- logNormCounts(se_train_subset)
se_test_subset  <- logNormCounts(se_test_subset)

logcounts_data_train <- assay(se_train_subset, "logcounts")
logcounts_data_test  <- assay(se_test_subset, "logcounts")

cell_labels <- colData(se_train_subset)$celltype_major
message("Prepared training and testing sets for SingleR.")

# ---------------------------
# 11. Run SingleR for Cell Type Prediction
# ---------------------------

pred.hesc <- SingleR(
  test = logcounts_data_test,
  ref = se_train_subset,
  assay.type.test = 1,
  labels = cell_labels
)

singleR_results_path <- "~/data/validation/singleR_celltypes_training_test.RData"
save(pred.hesc, file = singleR_results_path)
message("SingleR analysis completed and results saved to: ", singleR_results_path)

# ---------------------------
# 12. Evaluation
# ---------------------------

predicted_labels <- pred.hesc$labels %>% set_names(rownames(pred.hesc))
test_labels_df <- tibble(Cell = names(predicted_labels), Predicted = predicted_labels)
actual_labels_df <- tibble(Cell = colnames(logcounts_data_test),
                           Actual = colData(se_test_subset)$celltype_major)
fusion <- left_join(test_labels_df, actual_labels_df, by = "Cell")

# Overall accuracy
overall_accuracy_table <- table(fusion$Predicted == fusion$Actual)
print("Overall Accuracy Comparison:")
print(overall_accuracy_table)

correct_overall <- sum(fusion$Predicted == fusion$Actual, na.rm = TRUE)
total_overall <- sum(!is.na(fusion$Predicted) & !is.na(fusion$Actual))
overall_accuracy <- (correct_overall / total_overall) * 100
cat(sprintf("Overall Accuracy: %.2f%%\n", overall_accuracy))

# Accuracy for Cancer Epithelial cells only
fusion_cancer_epi <- fusion %>% filter(Actual == "Cancer Epithelial")
cancer_epi_accuracy_table <- table(fusion_cancer_epi$Predicted == fusion_cancer_epi$Actual)
print("Cancer Epithelial Accuracy Comparison:")
print(cancer_epi_accuracy_table)

correct_cancer_epi <- sum(fusion_cancer_epi$Predicted == fusion_cancer_epi$Actual, na.rm = TRUE)
total_cancer_epi <- sum(!is.na(fusion_cancer_epi$Predicted) & !is.na(fusion_cancer_epi$Actual))
cancer_epi_accuracy <- (correct_cancer_epi / total_cancer_epi) * 100
cat(sprintf("Cancer Epithelial Accuracy: %.2f%%\n", cancer_epi_accuracy))

# =============================================================================
# End of Script
# =============================================================================
