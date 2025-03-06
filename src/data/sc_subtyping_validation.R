# =============================================================================
# Single-Cell RNA-Seq Data Processing and Annotation Pipeline
# Author: Miguel Castresana
# Date: 2024-12-18
# Description: This script processes single-cell RNA-seq data for tumor samples,
#              normalizes the data, splits into training and validation sets,
#              integrates SingleR annotations, computes SC-subtype scores per patient,
#              and outputs the results.
# =============================================================================

library(tidyverse)
library(Seurat)
library(readxl)
library(caret)
library(SummarizedExperiment)
library(scuttle)
library(scCustomize)

# Define the "not in" operator
`%!in%` <- function(x, y) !(x %in% y)

# ---------------------------
# 2. Load Validation Data
# ---------------------------
# Function to load a Seurat object from an RData file
load_seurat_object <- function(filepath, object_name = "combined") {
  if (!file.exists(filepath)) stop("File not found: ", filepath)
  env <- new.env()
  load(filepath, envir = env)
  if (!exists(object_name, envir = env)) {
    stop("Object '", object_name, "' not found in file: ", filepath)
  }
  get(object_name, envir = env)
}

# Load validation data (assumes the file loads an object named "combined")
load("~/data/validation/validation_data_all_normalized")
DefaultAssay(combined) <- 'RNA'

# Clean and update metadata:
metadata_validation <- combined@meta.data %>% 
  as.data.frame() %>% 
  { rownames(.) <- NULL; . } %>% 
  column_to_rownames(var = "NAME") %>% 
  mutate(NAME = rownames(.))
combined@meta.data <- metadata_validation



# ---------------------------
# 4. Filter to Cancer Epithelial Cells
# ---------------------------
metadata_filtered <- metadata_validation %>% filter(celltype_major == "Cancer Epithelial")
validation_data <- subset(combined, cells = metadata_filtered$NAME)

# ---------------------------
# 5. Load Training Data
# ---------------------------
training_excel_path <- "~/data/discovery/41588_2021_911_MOESM4_ESM.xlsx"
if (!file.exists(training_excel_path)) stop("Training Excel file not found: ", training_excel_path)
training_data_excel <- read_excel(training_excel_path, sheet = 3)

training_samples <- training_data_excel %>% 
  filter(`SCTyper dataset` == "Training") %>% 
  pull(`Tumour ID`)

their_metadata_path <- "~/data/discovery/metadata_discovery.csv"
if (!file.exists(their_metadata_path)) stop("Metadata file not found: ", their_metadata_path)
their_metadata <- read_csv(their_metadata_path) %>% 
  column_to_rownames(var = "Cell") %>% 
  mutate(Cell = rownames(.))

their_metadata <- their_metadata %>% 
  filter(Patient %in% training_samples)
train_cells <- their_metadata %>% pull(Cell)

training_seurat_path <- "~/data/discovery/combined_discovery"
seurat_training <- load_seurat_object(training_seurat_path, object_name = "combined")
DefaultAssay(seurat_training) <- 'RNA'
training_data <- subset(seurat_training, cells = train_cells)
message("Training data loaded and subset to Cancer Epithelial cells.")

# ---------------------------
# 6. Load SC-Subtype Signatures
# ---------------------------
sc_subtype_signatures_path <- "~/data/discovery/NatGen_Supplementary_table_S4.csv"
if (!file.exists(sc_subtype_signatures_path)) stop("SC-subtype signatures file not found: ", sc_subtype_signatures_path)
sigdat <- read_csv(sc_subtype_signatures_path, col_types = cols())

# Replace only the first occurrence of a period with a dash in each element,
# then correct specific gene name "AP000769-1" to "AP000769.1"
sigdat = as.data.frame(apply(sigdat,2,function(x) sub("\\.", "-", x)))
sigdat[sigdat == "AP000769-1"] <- "AP000769.1"

# Extract unique signature genes from the four subtype columns, removing empties and NAs
temp_allgenes <- c(sigdat$Basal_SC, sigdat$Her2E_SC, sigdat$LumA_SC, sigdat$LumB_SC) %>% 
  discard(~ is.na(.x) || .x == "") %>% 
  unique()


# ---------------------------
# 7. Prepare Training Data for SingleR
# ---------------------------
counts_train <- GetAssayData(training_data, assay = "RNA", slot = "data")
missing_genes_train <- setdiff(temp_allgenes, rownames(counts_train))
if (length(missing_genes_train) > 0) {
  warning("Missing signature genes in training data: ", paste(missing_genes_train, collapse = ", "))
} else {
  message("All signature genes are present in training data.")
}

# ---------------------------
# 8. Compute SC-Subtypes for Validation Data by Patient
# ---------------------------
patients <- validation_data@meta.data$orig.ident %>% unique()

# Define helper function for centering scores
center_sweep <- function(x, row.w = rep(1, nrow(x)) / nrow(x)) {
  get_average <- function(v) sum(v * row.w) / sum(row.w)
  average <- apply(x, 2, get_average)
  sweep(x, 2, average)
}

sc_subtypes_all <- map(patients, function(patient) {
  
  cells_sub <- metadata_filtered %>% filter(orig.ident == patient) %>% pull(NAME)
  seurat_object2 <- subset(validation_data, cells = cells_sub)
  
  counts1 <- GetAssayData(seurat_object2, assay = "RNA", slot = "data")
  common_genes <- intersect(rownames(counts_train), rownames(counts1))
  
  seurat_object1 <- subset(training_data, features = common_genes)
  seurat_object2 <- subset(seurat_object2, features = common_genes)
  
  merged_seurat_object <- merge(x = seurat_object1, y = seurat_object2)
  merged_seurat_object <- JoinLayers(merged_seurat_object)
  DefaultAssay(merged_seurat_object) <- "RNA"
  
  Mydata <- ScaleData(merged_seurat_object, features = temp_allgenes, slot = "data")
  tocalc <- GetAssayData(Mydata, assay = "RNA", slot = "scale.data") %>% as.data.frame()
  tocalc <- tocalc[rownames(tocalc) %in% temp_allgenes, ]
  
  subtype_names <- c("Basal_SC", "Her2E_SC", "LumA_SC", "LumB_SC")
  # For each subtype, compute the mean score across signature genes
  outdat <- map(subtype_names, function(subtype) {
    genes <- unique(sigdat[[subtype]])
    genes <- genes[genes != ""]
    common <- intersect(rownames(tocalc), genes)
    if (length(common) == 0) {
      rep(0, ncol(tocalc))
    } else {
      colMeans(tocalc[common, , drop = FALSE], na.rm = TRUE)
    }
  })
  outdat <- do.call(rbind, outdat)
  rownames(outdat) <- subtype_names
  
  final <- outdat[rowSums(outdat, na.rm = TRUE) != 0, , drop = FALSE]
  final_df <- as.data.frame(final) %>% mutate(across(where(is.numeric), ~ round(., 4)))
  finalm <- as.matrix(final_df)
  
  finalmt <- as.data.frame(t(finalm))
  finalmt_sweep <- center_sweep(finalmt)
  Finalnames <- colnames(finalmt_sweep)[max.col(finalmt_sweep, ties.method = "first")]
  finalmt_sweep$SCSubtypeCall <- Finalnames
  
  finalmt_sweep %>% rownames_to_column(var = "Cell")
})

final_calls <- bind_rows(sc_subtypes_all) %>%
  mutate(Cell = str_replace(Cell, "^_", "")) %>%
  filter(Cell %in% colnames(validation_data@assays$RNA))

print(table(final_calls$SCSubtypeCall))

# ---------------------------
# 9. Save the Calls
# ---------------------------
final_calls <- final_calls %>% 
  mutate(NAME = Cell,
         patient = str_replace(Cell, "_.*", ""))

# Most common subtype per patient
result <- final_calls %>%
  group_by(patient) %>%
  summarize(MostCommonSubtypeCall = SCSubtypeCall[which.max(tabulate(match(SCSubtypeCall, SCSubtypeCall)))],
            .groups = "drop")

# If desired, join additional metadata (adjust column selections as needed)
final_calls <- left_join(final_calls,
                         metadata_filtered %>% rownames_to_column(var = "NAME") %>% select(NAME, col7, col9),
                         by = "NAME")

write.table(final_calls, "~/data/validation/validation_scsubtype", sep = ",", quote = FALSE, row.names = FALSE)

# =============================================================================
# End of Script
# =============================================================================
