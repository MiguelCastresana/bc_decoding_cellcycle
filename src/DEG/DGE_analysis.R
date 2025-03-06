# =============================================================================
# Single-Cell RNA-Seq Differential Expression Analysis with Command-Line Input
# Author: Miguel Castresana
# Date: 2025-01-23
# Description: This script performs differential expression analysis using limma
#              on single-cell RNA-seq data. It loads a Seurat object and metadata,
#              constructs a design matrix with interaction terms, fits linear models,
#              extracts top differentially expressed genes, and conducts pathway
#              enrichment analysis using FGSEA.
# =============================================================================

# -------------------------
# 1. Load Required Libraries
# -------------------------
library(limma)
library(edgeR)
library(plyr)
library(dplyr)
library(statmod)
library(irlba)
library(doParallel)
library(igraph)
library(magrittr)
library(Seurat)
library(msigdbr)
library(fgsea)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(openxlsx)
library(pheatmap)
library(yaml)

# Define the "not in" operator
'%!in%' <- function(x, y) !(x %in% y)

# -------------------------
# 2. Parse Configuration
# -------------------------
# (This section mirrors the command-line/configuration style of the uploaded script)
config_file_path <- "~/src/repository/repository_final/DEG//config.yaml"
config <- yaml::read_yaml(config_file_path)$default
validate_parameters <- function(config) {
  required_params <- c("data_dir", "results_dir", "seurat_path", 
                       "metadata_path", "genes_file","analysis_var")
  
  # Check if all required parameters exist
  missing_params <- setdiff(required_params, names(config))
  if (length(missing_params) > 0) {
    stop(paste("❌ Missing parameters in YAML config:", paste(missing_params, collapse = ", ")))
  }
  
  # Validate that each parameter is a single character string
  for (param in required_params) {
    if (!is.character(config[[param]]) || length(config[[param]]) != 1) {
      stop(paste("❌ Error:", param, "must be a single character string."))
    }
  }
  
  # Validate `analysis_var`
  if (!tolower(config$analysis_var) %in% c("phase", "subtype","no_phase")) {
    stop("❌ Error: 'analysis_var' must be either 'phase' ,subtype'. or ''no_phase")
  }
  
  
  # Construct full paths for the Seurat object and metadata files (assumed relative to data_dir)
  seurat_full_path <- file.path(config$data_dir, config$seurat_path)
  metadata_full_path <- file.path(config$data_dir, config$metadata_path)
  
  # Validate file existence
  if (!file.exists(seurat_full_path)) {
    stop(paste("❌ Error: Seurat object file does not exist:", seurat_full_path))
  }
  if (!file.exists(metadata_full_path)) {
    stop(paste("❌ Error: Metadata file does not exist:", metadata_full_path))
  }
  if (!file.exists(config$genes_file)) {
    stop(paste("❌ Error: Genes file does not exist:", config$genes_file))
  }
  
  # If the output directory does not exist, create it recursively
  if (!dir.exists(config$results_dir)) {
    message("Output directory does not exist. Creating: ", config$results_dir)
    dir.create(config$results_dir, recursive = TRUE)
  }
  
  message("✅ All configuration checks passed successfully.")
  
  # Add computed full paths to config for later use
  config$seurat_full_path <- seurat_full_path
  config$metadata_full_path <- metadata_full_path
  
  return(config)
}

# Safe validation wrapper
safe_validate_parameters <- function(config) {
  tryCatch(
    {
      return(validate_parameters(config))  # Return the validated config
    },
    error = function(e) {
      message("❌ Caught an error in validation: ", e$message)
      return(NULL)  # Return NULL in case of an error
    }
  )
}


# Run validation and get the updated config
config <- safe_validate_parameters(config)
if (is.null(config)) {
  stop("Configuration validation failed.")
}

# Assign validated parameters
data_dir          <- config$data_dir
results_dir       <- config$results_dir
seurat_full_path  <- config$seurat_full_path
metadata_full_path<- config$metadata_full_path
genes_file        <- config$genes_file
analysis_var      <- tolower(config$analysis_var)

# -------------------------
# 3. Load and Prepare Seurat Object
# -------------------------
load_seurat_object <- function(filepath) {
  if (!file.exists(filepath)) {
    stop("Seurat object file not found: ", filepath)
  }
  env <- new.env()
  load(filepath, envir = env)
  objects_loaded <- ls(env)
  if (length(objects_loaded) == 0) {
    stop("No objects found in file: ", filepath)
  }
  if (length(objects_loaded) > 1) {
    message("Multiple objects found in file. Using the first object: ", objects_loaded[1])
  }
  seurat_obj <- get(objects_loaded[1], envir = env)
  message("Loaded Seurat object '", objects_loaded[1], "' from: ", filepath)
  return(seurat_obj)
}

# Usage:
seurat_obj <- load_seurat_object(seurat_full_path)


# -------------------------
# 4. Load and Merge Metadata
# -------------------------
load_metadata <- function(filepath) {
  if (!file.exists(filepath)) {
    stop("Metadata file not found: ", filepath)
  }
  
  metadata <- read_csv(filepath, col_types = cols()) %>%
    rename_with(tolower)
  
  # If a column named "name" exists, set it as row names
  if ("name" %in% names(metadata)) {
    metadata <- column_to_rownames(metadata, var = "name")
  }
  
  message("Loaded metadata from: ", filepath)
  return(metadata)
}
metadata <- load_metadata(metadata_full_path)


# ---------------------------
# 5. Clean data
# ---------------------------

# Subset Seurat object to Cancer Epithelial cells using metadata Cell IDs
seurat_obj <- subset(seurat_obj, cells = metadata$Cell)
message("Subset Seurat object to Cancer Epithelial cells: ", 
        ncol(GetAssayData(seurat_obj, slot = "data")), " cells.")

cellinfo <- metadata %>%
  mutate(
    Batch   = orig.ident,
    Patient = orig.ident,
    Group   = orig.ident
  ) %>%
  # Select and rename columns; assume original names are in lowercase
  select(Cell = cell, Phase = phase, SC_subtype = subtype) %>%
  # Recode SC_subtype, simplify Phase, and set factor levels explicitly
  mutate(
    SC_subtype = recode(SC_subtype, "Her2E_SC" = "Her2_SC"),
    Phase = factor(str_split(Phase, "_", simplify = TRUE)[, 1],
                   levels = c("G1", "S", "G2M")),
    SC_subtype = factor(SC_subtype, levels = c("LumA_SC", "LumB_SC", "Her2_SC", "Basal_SC"))
  ) %>%
  select(Cell, SC_subtype, Phase)


cellinfo <- cellinfo %>%
      mutate(SC_subtype = factor(SC_subtype, levels = c("LumA_SC", "LumB_SC", "Her2_SC", "Basal_SC")))

# Check the levels
print("SC_subtype levels:")
print(levels(cellinfo$SC_subtype))




# Get the assay data (as a matrix) from the Seurat object in the "data" slot
datmat <- GetAssayData(object = seurat_obj, slot = "data") %>% 
  # Keep only the columns (cells) that are in cellinfo$Cell
  { .[, intersect(colnames(.), cellinfo$Cell)] } %>% 
  # Filter rows: keep only genes that are in your genes file
  { 
    genes_f <- read_csv(genes_file, col_types = cols())
    .[rownames(.) %in% genes_f$genes, ]
  } %>% 
  # Reorder columns so they are in the same order as cellinfo$Cell
  { .[, match(cellinfo$Cell, colnames(.))] }

if (!identical(colnames(datmat), cellinfo$Cell)) {
  stop("Error: Column names of expression data and cellinfo do not match.")
}
message("Expression data and cell information are properly aligned.")

# -------------------------
# 7. Define Design matrix
# -------------------------



if(analysis_var=="phase" | analysis_var == "subtype"){
  
  # Build the design matrix (using interaction terms)
  design_matrix <- model.matrix(~ SC_subtype * Phase, data = cellinfo)

}else if(analysis_var == "no_phase"){
  
  # Create design matrix with interaction terms
  design_matrix <- model.matrix(
    ~ SC_subtype,
    data = cellinfo
  )
}else{
  stop("analysis_var must be either 'phase', 'subtype' or 'no_phase'")
}

colnames(design_matrix) <- gsub(":", "_", colnames(design_matrix))
message("Design Matrix Columns:")
print(colnames(design_matrix))

# -------------------------
# 7. Define Contrasts for Differential Expression
# -------------------------
# Rename intercept for clarity
colnames(design_matrix)[1] <- "Intercept"

if (analysis_var == "phase") {
  contrast_matrix <- makeContrasts(
    # G1 Phase comparisons
    LumAvsOther_G1 = -(SC_subtypeLumB_SC + SC_subtypeHer2_SC + SC_subtypeBasal_SC) / 3,
    LumBvsOther_G1 = SC_subtypeLumB_SC - (0 + SC_subtypeHer2_SC + SC_subtypeBasal_SC) / 3,
    Her2vsOther_G1 = SC_subtypeHer2_SC - (0 + SC_subtypeLumB_SC + SC_subtypeBasal_SC) / 3,
    BasalvsOther_G1 = SC_subtypeBasal_SC - (0 + SC_subtypeLumB_SC + SC_subtypeHer2_SC) / 3,
    
    # S Phase comparisons
    LumAvsOther_S = -((SC_subtypeLumB_SC + SC_subtypeHer2_SC + SC_subtypeBasal_SC) / 3 +
                        (SC_subtypeLumB_SC_PhaseS + SC_subtypeHer2_SC_PhaseS + SC_subtypeBasal_SC_PhaseS) / 3),
    LumBvsOther_S = (SC_subtypeLumB_SC + SC_subtypeLumB_SC_PhaseS) -
      ((0 + SC_subtypeHer2_SC + SC_subtypeBasal_SC + 0 +
          SC_subtypeHer2_SC_PhaseS + SC_subtypeBasal_SC_PhaseS) / 3),
    Her2vsOther_S = (SC_subtypeHer2_SC + SC_subtypeHer2_SC_PhaseS) -
      ((0 + SC_subtypeLumB_SC + SC_subtypeBasal_SC + 0 +
          SC_subtypeLumB_SC_PhaseS + SC_subtypeBasal_SC_PhaseS) / 3),
    BasalvsOther_S = (SC_subtypeBasal_SC + SC_subtypeBasal_SC_PhaseS) -
      ((0 + SC_subtypeLumB_SC + SC_subtypeHer2_SC + 0 +
          SC_subtypeLumB_SC_PhaseS + SC_subtypeHer2_SC_PhaseS) / 3),
    
    # G2M Phase comparisons
    LumAvsOther_G2M = -((SC_subtypeLumB_SC + SC_subtypeHer2_SC + SC_subtypeBasal_SC) / 3 +
                          (SC_subtypeLumB_SC_PhaseG2M + SC_subtypeHer2_SC_PhaseG2M + SC_subtypeBasal_SC_PhaseG2M) / 3),
    LumBvsOther_G2M = (SC_subtypeLumB_SC + SC_subtypeLumB_SC_PhaseG2M) -
      ((0 + SC_subtypeHer2_SC + SC_subtypeBasal_SC + 0 +
          SC_subtypeHer2_SC_PhaseG2M + SC_subtypeBasal_SC_PhaseG2M) / 3),
    Her2vsOther_G2M = (SC_subtypeHer2_SC + SC_subtypeHer2_SC_PhaseG2M) -
      ((0 + SC_subtypeLumB_SC + SC_subtypeBasal_SC + 0 +
          SC_subtypeLumB_SC_PhaseG2M + SC_subtypeBasal_SC_PhaseG2M) / 3),
    BasalvsOther_G2M = (SC_subtypeBasal_SC + SC_subtypeBasal_SC_PhaseG2M) -
      ((0 + SC_subtypeLumB_SC + SC_subtypeHer2_SC + 0 +
          SC_subtypeLumB_SC_PhaseG2M + SC_subtypeHer2_SC_PhaseG2M) / 3),
    
    levels = colnames(design_matrix)
  )
  
} else if (analysis_var == "subtype") {
  contrast_matrix <- makeContrasts(
    # LumA Contrasts
    G1LumA_vs_OtherLumA = - (PhaseS + PhaseG2M) / 2,
    SLumA_vs_OtherLumA = PhaseS - (PhaseG2M) / 2,
    G2MLumA_vs_OtherLumA = PhaseG2M - PhaseS / 2,
    
    # LumB Contrasts
    G1LumB_vs_OtherLumB = - (PhaseS + PhaseG2M + SC_subtypeLumB_SC_PhaseS + SC_subtypeLumB_SC_PhaseG2M) / 2,
    SLumB_vs_OtherLumB = (PhaseS + SC_subtypeLumB_SC_PhaseS) - (PhaseG2M + SC_subtypeLumB_SC_PhaseG2M) / 2,
    G2MLumB_vs_OtherLumB = (PhaseG2M + SC_subtypeLumB_SC_PhaseG2M) - (PhaseS + SC_subtypeLumB_SC_PhaseS) / 2,
    
    # Her2 Contrasts
    G1Her2_vs_OtherHer2 = - (PhaseS + PhaseG2M + SC_subtypeHer2_SC_PhaseS + SC_subtypeHer2_SC_PhaseG2M) / 2,
    SHer2_vs_OtherHer2 = (PhaseS + SC_subtypeHer2_SC_PhaseS) - (PhaseG2M + SC_subtypeHer2_SC_PhaseG2M) / 2,
    G2MHer2_vs_OtherHer2 = (PhaseG2M + SC_subtypeHer2_SC_PhaseG2M) - (PhaseS + SC_subtypeHer2_SC_PhaseS) / 2,
    
    # Basal Contrasts
    G1Basal_vs_OtherBasal = - (PhaseS + PhaseG2M + SC_subtypeBasal_SC_PhaseS + SC_subtypeBasal_SC_PhaseG2M) / 2,
    SBasal_vs_OtherBasal = (PhaseS + SC_subtypeBasal_SC_PhaseS) - (PhaseG2M + SC_subtypeBasal_SC_PhaseG2M) / 2,
    G2MBasal_vs_OtherBasal = (PhaseG2M + SC_subtypeBasal_SC_PhaseG2M) - (PhaseS + SC_subtypeBasal_SC_PhaseS) / 2,
    
    levels = colnames(design_matrix)
  )
  
} else if (analysis_var == "no_phase") {
  
  contrast_matrix <- makeContrasts(
    
    # LumA vs Others
    LumAvsOther = -(SC_subtypeLumB_SC + SC_subtypeHer2_SC + SC_subtypeBasal_SC) / 3,
    
    # LumB vs Others
    LumBvsOther = SC_subtypeLumB_SC - (SC_subtypeHer2_SC + SC_subtypeBasal_SC + 0) / 3,
    
    # Her2 vs Others
    Her2vsOther = SC_subtypeHer2_SC - (SC_subtypeLumB_SC + SC_subtypeBasal_SC + 0) / 3,
    
    # Basal vs Others
    BasalvsOther = SC_subtypeBasal_SC - (SC_subtypeLumB_SC + SC_subtypeHer2_SC + 0) / 3,
    
    levels = colnames(design_matrix)
  )
  
  
} else {
  stop("analysis_var must be either 'phase', 'subtype' or 'no_phase'")
}

# -------------------------
# 8. Fit Linear Models and Apply Contrasts
# -------------------------
lmfit <- lmFit(datmat, design_matrix) %>%
  { message("Fitted linear model using limma."); . } %>%
  eBayes(trend = TRUE, robust = TRUE) %>%
  { message("Applied empirical Bayes moderation."); . }

contrasts_fit <- contrasts.fit(lmfit, contrast_matrix) %>%
  eBayes(trend = TRUE, robust = TRUE) %>%
  { message("Fitted contrasts and applied empirical Bayes moderation."); . }

# -------------------------
# 9. Extract Top Differentially Expressed Genes
# -------------------------
# For each contrast (column of contrast_matrix), extract the top table,
# then combine them by gene.
results_list <- map(seq_len(ncol(contrast_matrix)), function(i) {
  contrast_name <- colnames(contrast_matrix)[i]
  tbl <- topTable(contrasts_fit, coef = i, n = Inf, adjust = "BH")
  message("Extracted topTable for contrast: ", contrast_name)
  tibble(gene = rownames(tbl),
         logFC = tbl$logFC,
         tstats = tbl$t) %>%
    rename_with(~ paste0(., "_", contrast_name), .cols = c("logFC", "tstats"))
})

# Merge all contrast results by the gene column
merged_df <- purrr::reduce(results_list, full_join, by = "gene") %>% 
  filter(rowSums(across(-gene, ~ is.na(.))) < (ncol(.) - 1)) %>%
  column_to_rownames("gene")

# Retain only the t-statistics columns from each contrast
res_table <- merged_df %>% 
  select(contains("tstats")) %>%
  rename_with(~ gsub("_tstats", "", .x))

# Save differential expression results
de_results_path <- file.path(results_dir, "limmatrend_DE_results.RData")
save(res_table, file = de_results_path)
message("Saved differential expression results to: ", de_results_path)
# -------------------------
# 10. Pathway Enrichment Analysis using FGSEA
# -------------------------

# -------------------------
# Load MSigDB Hallmark gene sets and prepare pathway list
# -------------------------
msigdb_h <- msigdbr(species = "Homo sapiens", category = "H") %>%
  distinct(gs_name, gene_symbol, .keep_all = TRUE)

# Create a pathway list: each pathway (gs_name) maps to a character vector of unique gene symbols.
path_df <- msigdb_h %>%
  dplyr::group_by(gs_name) %>%
  dplyr::summarise(genes = list(unique(gene_symbol)), .groups = "drop") %>%
  arrange(gs_name)

# Convert the summarised tibble into a named list
path_list <- setNames(path_df$genes, path_df$gs_name)

message("Prepared MSigDB Hallmark pathways for enrichment analysis.")
print(paste("Number of pathways:", length(path_list)))

# -------------------------
# Perform FGSEA for each contrast using tidyverse/purrr
# -------------------------
names_comparisons <- colnames(res_table)

fgsea_results <- map(names_comparisons, function(contrast_name) {
  # Extract gene list for current contrast, assign names, and sort in decreasing order
  gene_list <- res_table[, contrast_name] %>%
    set_names(rownames(res_table)) %>%
    sort(decreasing = TRUE)
  
  message("Running FGSEA for contrast: ", contrast_name, " with ", length(gene_list), " genes.")
  
  # Run FGSEA
  enrichment <- fgseaMultilevel(
    pathways    = path_list,
    stats       = gene_list,
    nproc       = 1,
    eps         = 1e-50,
    sampleSize  = 2000,
    nPermSimple = 100000
  )
  
  # Check if the expected 'pathway' column exists
  if (!"pathway" %in% colnames(enrichment)) {
    stop("FGSEA output does not contain a 'pathway' column for contrast: ", contrast_name)
  }
  
  enrichment %>% arrange(pval)
})
names(fgsea_results) <- names_comparisons

# -------------------------
# Save FGSEA results to Excel and RData file
# -------------------------
fgsea_output_path <- file.path(results_dir, "hallmarks_limma_DE_results.xlsx")
openxlsx::write.xlsx(fgsea_results, file = fgsea_output_path)
message("Saved FGSEA results to Excel: ", fgsea_output_path)

fgsea_rdata_path <- file.path(results_dir, "hallmarks_limma_DE_results.RData")
save(fgsea_results, file = fgsea_rdata_path)
message("Saved FGSEA results to RData: ", fgsea_rdata_path)
# =============================================================================
# End of Script
# =============================================================================
