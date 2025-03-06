
# =============================================================================
# Combined SCENIC Downstream Analysis Script
# Author: Miguel Castresana
# Date: 2025-01-23
# Description:Combined SCENIC Downstream Analysis Script (AUCell, RSS, pathway enrichment)


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




# -------------------------
# Parse command-line arguments
# -------------------------

# Read the configuration file (adjust the path if necessary)
config_file_path <- "~/src/repository/repository_final/GRN/config.yaml"
config <- yaml::read_yaml(config_file_path)$default

validate_parameters <- function(config) {
  required_params <- c("analysis_var", "z_threshold", "relative_activity_cutoff", 
                       "min_qvalue", "loom_base_path", "metadata_path", 
                       "final_genes_path", "output_dir")
  
  # Check if all required parameters exist
  missing_params <- setdiff(required_params, names(config))
  if (length(missing_params) > 0) {
    stop(paste("❌ Missing parameters in YAML config:", paste(missing_params, collapse = ", ")))
  }
  
  # Convert `analysis_var` to lowercase
  if (!is.null(config$analysis_var) && is.character(config$analysis_var)) {
    config$analysis_var <- tolower(config$analysis_var)
  } else {
    stop("❌ Error: 'analysis_var' must be a character string.")
  }
  
  # Validate `analysis_var`
  if (!tolower(config$analysis_var) %in% c("phase", "subtype")) {
    stop("❌ Error: 'analysis_var' must be either 'phase' or 'subtype'.")
  }
  
  # Validate numeric parameters
  numeric_params <- c("z_threshold", "relative_activity_cutoff", "min_qvalue")
  numeric_limits <- list(
    z_threshold = c(0, Inf),
    relative_activity_cutoff = c(0, 1),
    min_qvalue = c(0, 1)
  )
  
  for (param in numeric_params) {
    value <- config[[param]]
    if (!is.numeric(value) || length(value) != 1 || value < numeric_limits[[param]][1] || value > numeric_limits[[param]][2]) {
      stop(paste("❌ Error:", param, "must be a numeric value in range", numeric_limits[[param]][1], "to", numeric_limits[[param]][2]))
    }
  }
  
  # Validate file paths
  if (!dir.exists(config$loom_base_path)) {
    stop(paste("❌ Error: Loom base directory does not exist:", config$loom_base_path))
  }
  if (!file.exists(config$metadata_path)) {
    stop(paste("❌ Error: Metadata file does not exist:", config$metadata_path))
  }
  if (!file.exists(config$final_genes_path)) {
    stop(paste("❌ Error: Final genes file does not exist:", config$final_genes_path))
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
analysis_var <- config$analysis_var
z_threshold <- config$z_threshold
relative_activity_cutoff <- config$relative_activity_cutoff
min_qvalue <- config$min_qvalue
loom_base_path <- config$loom_base_path
metadata_path <- config$metadata_path
final_genes_path <- config$final_genes_path
output_dir <- config$output_dir


# Define analysis groups; for discovery mode we use cell phases, for validation mode you might switch to cell types.
if(analysis_var == "phase"){
  network_groups <- c("G1", "S", "G2M")           
  comparison_group <- c("LumA", "LumB", "Her2", "Basal")  # not used in discovery mode
  analysis_comparison = "subtype"
} else if(analysis_var == "subtype"){
  network_groups <- c("LumA", "LumB", "Her2", "Basal")
  comparison_group <- c("G1", "S", "G2M")
  analysis_comparison = "phase"
}else{
  print("Try again")
}



process_loom <- function(loom_path,current_group) {
  
  
  # loom_path = current_loom_path
  loom <- open_loom(loom_path)
  
  # Extract regulon information
  regulons_incidMat <- get_regulons(loom, column.attr.name = "Regulons")
  regulons <- regulonsToGeneLists(regulons_incidMat) %>%
    enframe(name = "TranscriptionFactor", value = "Gene") %>%
    unnest(Gene)
  regulons$comparison_group <- rep(current_group, nrow(regulons))
  
  # Get AUcell results
  regulonAUC <- get_regulons_AUC(loom, column.attr.name = "RegulonsAUC")
  total_cells <- dim(regulonAUC@colData)[1]
  regulonAucThresholds <- get_regulon_thresholds(loom)
  embeddings <- get_embeddings(loom)
  
  # Load metadata and filter by current group using analysis_var
  meta <- read_csv(metadata_path, col_types = cols())
  names(meta) <- tolower(names(meta))
  meta <- meta %>% 
    column_to_rownames("cell") %>%
    mutate(subtype = str_replace_all(subtype, "_SC", "")) %>%
    filter(.data[[analysis_var]] %in% current_group)
  
  # Build cellClusters from AUcell data
  # Build cellClusters from AUcell data
  cellClusters <- tibble(cell = colnames(regulonAUC))
  # Here, decide which column to merge on; if analysis_var is "subtype", use that,
  # otherwise use the provided column (here we assume "phase")
  merge_column <- analysis_comparison
  cellClusters <- cellClusters %>% 
    left_join(meta %>% rownames_to_column("cell") %>% select(cell, all_of(merge_column)),
              by = "cell") %>%
    mutate(clusters = factor(.data[[merge_column]], levels =  comparison_group)) %>%
    column_to_rownames("cell")
  
  close_loom(loom)
  
  # Remove duplicate cells from regulonAUC
  regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)), ]
  
  # Split cells by cluster and calculate average activity per group
  cellsPerCluster <- split(rownames(cellClusters), cellClusters$clusters)
  regulonActivity_bycellType <- sapply(cellsPerCluster, function(cells) rowMeans(getAUC(regulonAUC)[, cells]))
  
  
  topRegulators <- melt(regulonActivity_bycellType) %>%
    mutate(geneset = paste(Var2, "__", Var1, sep = ""),
           !!sym(analysis_var) := current_group) %>%
    rename(Regulon = Var1,
           !!sym(analysis_comparison) := Var2,
           RelativeActivity = value)
  
  # Optionally, select the top N regulons per cell type (if desired)
  topRegulators <- topRegulators %>%
    group_by(!!sym(analysis_var)) %>%
    arrange(desc(RelativeActivity)) %>%
    ungroup()
  
  # Compute RSS in a tidyverse manner
  selectedResolution <- "clusters"
  rss <- calcRSS(AUC = getAUC(regulonAUC),
                 cellAnnotation = cellClusters[colnames(regulonAUC), selectedResolution])
  rssNorm <- scale(rss)
  rssNorm[rssNorm < 0] <- 0
  rssSubset <- rssNorm
  rssSubset[rssSubset < z_threshold] <- 0
  
  # Convert rss into a long-format tibble with a dynamic column name
  rss.df <- as.data.frame(rss) %>%
    rownames_to_column("Topic") %>%
    pivot_longer(
      cols = -Topic,
      names_to = analysis_comparison,   # dynamically defined column name
      values_to = "RSS"
    ) %>%
    mutate(geneset = str_c(.data[[analysis_comparison]], "__", Topic))
  
  # Convert rssNorm into a long-format tibble with the same dynamic column name
  rssNorm.df <- as.data.frame(rssNorm) %>%
    rownames_to_column("Topic") %>%
    pivot_longer(
      cols = -Topic,
      names_to = analysis_comparison,   # dynamically defined column name
      values_to = "Z"
    ) %>%
    mutate(geneset = str_c(.data[[analysis_comparison]], "__", Topic))
  
  # Left join the two data frames by the dynamic key 'geneset' and add a column named by analysis_var
  sub_list <- left_join(rss.df, rssNorm.df %>% select(geneset, Z), by = "geneset") %>%
    mutate(!!sym(analysis_var) := current_group)
  
  return(list(regulons = regulons, 
              regulonAUC = regulonAUC, 
              cellClusters = cellClusters, 
              topRegulators = topRegulators,
              total_cells = total_cells,
              rss = sub_list))
}


dfToNamedMatrix <- function(df) {
  # For now, simply return a matrix version of the data frame.
  return(as.matrix(df))
}

# -------------------------
# Main Pipeline: Process all groups
# -------------------------
all_topRegulators <- list()
all_regulons_list <- list()
all_total_cells <- numeric()
rss_list <- list()

for(current_group in network_groups) {
  current_loom_path <- file.path(loom_base_path, current_group, "pbmc10k_scenic_integrated-output.loom")
  res <- process_loom(current_loom_path, current_group)
  
  all_topRegulators[[current_group]] <- res$topRegulators
  all_regulons_list[[current_group]] <- res$regulons
  all_total_cells[current_group] <- res$total_cells
  rss_list[[current_group]] <- res$rss
}

regulons_allasone <- bind_rows(all_topRegulators)
regulons_allasone$geneset <- paste(regulons_allasone$comparison_group, "__", regulons_allasone$Regulon, sep = "")

regulons_list_df <- bind_rows(all_regulons_list)
regulons_list_df$geneset <- paste(regulons_list_df$comparison_group, "__", regulons_list_df$TranscriptionFactor, sep = "")

topRegulators <- bind_rows(all_topRegulators)

# RSS
rss_list <- bind_rows(rss_list)

write.table(rss_list, file = paste0(output_dir, "RSS_results.csv"),
            sep = ";", row.names = FALSE, quote = FALSE)

# -------------------------
# Gene Mapping
# -------------------------
ensembl <- mapIds(org.Hs.eg.db, keys = regulons_list_df$Gene,
                  column = "ENSEMBL", keytype = "SYMBOL", multiVals = "first")
regulons_list_df$ensembl <- ensembl

regulons_list_df <- regulons_list_df[regulons_list_df$TranscriptionFactor %in% topRegulators$Regulon, ]
regulons_list_df <- regulons_list_df[!is.na(regulons_list_df$ensembl), c("ensembl", "geneset")]

cutoff <- quantile(topRegulators$RelativeActivity, probs = relative_activity_cutoff)

# -------------------------
# HALLMARKS GENE MAPPING
# -------------------------
ens2symbol <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(regulons_list_df$ensembl),
  columns = c("SYMBOL", "ENTREZID"),
  keytype = "ENSEMBL"
) %>%
  filter(complete.cases(.)) %>%
  distinct(ENSEMBL, .keep_all = TRUE) %>%
  rename(ensembl = ENSEMBL, symbol = SYMBOL, entrez = ENTREZID)

grch38_mapping <- grch38 %>%
  filter(ensgene %in% regulons_list_df$ensembl) %>%
  select(ensembl = ensgene, symbol, entrez) %>%
  filter(complete.cases(.), !symbol %in% c("", NA)) %>%
  mutate(entrez = as.character(entrez))

combined_mappings <- bind_rows(ens2symbol, grch38_mapping) %>%
  group_by(ensembl) %>%
  summarise(symbol = coalesce(symbol[1], symbol[2]),
            entrez = coalesce(entrez[1], entrez[2])) %>%
  ungroup() %>%
  filter(complete.cases(.))

regulons_list_new <- regulons_list_df %>%
  inner_join(combined_mappings, by = "ensembl") %>%
  select(entrez, geneset) %>%
  filter(complete.cases(.)) %>%
  distinct(entrez, geneset)

# -------------------------
# Pathway Enrichment Analysis
# -------------------------
msig_h <- msigdbr(species = "Homo sapiens", category = "H") %>%
  select(gs_name, entrez_gene, gene_symbol) %>%
  rename(ont = gs_name, symbol = gene_symbol, entrezid = entrez_gene)
msig_h <- msig_h[!duplicated(msig_h[c("ont", "entrezid")]), ]

total_genes <- read.delim(final_genes_path, sep = ";")
total_genes <- mapIds(org.Hs.eg.db,
                      keys = total_genes$genes,
                      column = "ENTREZID",
                      keytype = "SYMBOL",
                      multiVals = "first")
total_genes <- total_genes[!is.na(total_genes)]

# For each unique geneset, run the enrichment and combine results
gea_output <- regulons_list_new %>%
  distinct(geneset, entrez) %>% 
  group_by(geneset) %>%
  nest() %>%
  mutate(enrich = map(data, ~ {
    clusterProfiler::enricher(gene = .x$entrez,
                              universe = as.vector(total_genes),
                              TERM2GENE = msig_h)
  })) %>%
  ungroup() %>%
  # Extract enrichment results (if available) and add the geneset column
  mutate(results = map(enrich, function(x) {
    if (!is.null(x)) {
      res <- x@result[, c("ID", "GeneRatio", "pvalue", "p.adjust")]
      res$geneset <- unique(x@result$geneset)  # should be identical to the grouping geneset
      return(as_tibble(res))
    } else {
      return(NULL)
    }
  })) %>%
  filter(!map_lgl(results, is.null)) %>%
  select(geneset, results) %>%
  unnest(results)

# Rename columns to match your original naming convention
gea_output <- gea_output %>%
  rename(Geneset = geneset,
         Pathway = ID,
         Overlap = GeneRatio,
         `p-value` = pvalue,
         `q-value` = p.adjust) %>%
  mutate(FWER = p.adjust(`p-value`, method = "bonferroni"))


# First, build a combined identifier if not already done:
# If topRegulators or rss_list contain separate phase and subtype columns, combine them.
# For example, if topRegulators has "comparison" (which might be cell phase) and topRegulators has "cellType" (e.g., LumA):
topRegulators <- topRegulators %>% 
  mutate(both = paste(subtype, phase, sep = " "))

# Similarly, ensure rss_list has a combined "both" column:
rss_list <- rss_list %>%
  mutate(both = paste(subtype, phase, sep = " "))

# Now, perform the fusion step:
# For each group in unique combined identifiers (unique_cols), we:
# - Subset topRegulators and rss_list for that group.
# - Filter rss_list so that only entries with Z > z_threshold remain.
# - Left join the two by the common 'geneset' column.
fus_tibble <- tibble(both = unique(topRegulators$both)) %>%
  mutate(fusion = map(both, function(b) {
    topregulon_sub <- topRegulators %>% filter(both == b)
    rss_list_sub <- rss_list %>% filter(both == b, Z > z_threshold)
    
    # Left join topregulon_sub and rss_list_sub on "geneset"
    fus <- left_join(topregulon_sub %>% select(Regulon, geneset, RelativeActivity),
                     rss_list_sub %>% select(geneset, subtype, RSS, Z, phase),
                     by = "geneset")
    # Ensure RSS is numeric and filter complete cases:
    fus <- fus %>% mutate(RSS = as.numeric(RSS)) %>% drop_na()
    fus
  }))

# Combine all fusion results into one data frame:
fus_df <- fus_tibble %>% unnest(fusion)

# Adjust column names and factor levels:
fus_df <- fus_df %>%
  mutate(reg_comparison = paste(.data[[analysis_var]], "__", Regulon, sep = ""))

new_cols <- c(analysis_var,"Regulon")
gea_output <- gea_output %>%
  separate(Geneset, into = new_cols, sep = "__", remove = FALSE)

gea_output <- gea_output %>%
  rename(!!analysis_var := !!sym(analysis_var)) %>%
  mutate(!!analysis_var := factor(.data[[analysis_var]], levels = network_groups),
         reg_comparison = paste(.data[[analysis_var]], "__", Regulon, sep = ""))

# -------------------------
# HALLMARKS Enrichment Section
# -------------------------

# First, filter gea_output1 by q-value threshold:
gea_output1 <- gea_output %>% 
  filter(`q-value` < min_qvalue)

# Get the unique reg_comparison values from fus_df
unique_reg <- unique(fus_df$reg_comparison)

# For each unique reg_comparison, filter gea_output1_filtered and combine results
gea_cases <- unique_reg %>% 
  map_df(~ gea_output1 %>% filter(reg_comparison == .x)) %>% 
  select(-!!analysis_var, -Regulon)



# Perform a left join using 'reg_phase' as the join key:
final <- inner_join(gea_cases, fus_df, by = "reg_comparison",relationship = "many-to-many")



# Further processing: modify pathway names, filter by RelativeActivity cutoff, and create a Gene column.
analysis <- final %>%
  mutate(both = paste(subtype, phase, sep = " "),
         Pathway = gsub("_", " ", Pathway),
         Pathway = sub("HALLMARK ", "", Pathway)) %>%
  filter(RelativeActivity >= cutoff) %>%
  mutate(Gene = sub("_\\(.*\\)", "", Regulon))


write_csv(analysis, file = paste0(output_dir, "GRN_pathways_enrichment.csv"))


# -------------------------
# Optional: Create Summary Matrix for Pathway Overlaps
# -------------------------
# Generate all combinations and order them by phase
unique_cols <- expand_grid(
  !!sym(analysis_comparison) := network_groups,
  !!sym(analysis_var) := comparison_group
) %>%
  mutate(Combination = paste(.data[[analysis_comparison]], .data[[analysis_var]])) %>%
  pull(Combination)

unique_paths = unique(analysis$Pathway)

# Process the analysis data:
# 1. Remove duplicate rows based on 'both' and 'Pathway'
# 2. Create a new column 'sig' that is 1 when q-value < 0.05 and 0 otherwise.
# 3. Pivot the data so that each row is a pathway and each column is a combination.
comparaciones1 <- analysis %>%
  distinct(both, Pathway, .keep_all = TRUE) %>%
  mutate(sig = as.integer(`q-value` < 0.05)) %>%
  select(both, Pathway, sig) %>%
  pivot_wider(names_from = both, values_from = sig, values_fill = list(sig = 0)) %>%
  # Ensure that the columns appear in the desired order (unique_cols):
  select(any_of(unique_cols)) %>%
  # Add the pathway names as a column so we can reassemble later:
  mutate(Pathway = unique_paths[match(row_number(), row_number())]) %>%  # This line is optional if row order is already correct
  relocate(Pathway)

# If some pathways in unique_paths did not appear in the data, we need to ensure they are included.
# Here we ensure that the final table has a row for every pathway in unique_paths:
comparaciones1 <- tibble(Pathway = unique_paths) %>%
  left_join(comparaciones1, by = "Pathway") %>%
  mutate(across(-Pathway, ~ replace_na(.x, 0)))  # replace any missing values with 0

# Now, compute the sum of 1s for each pathway (row) and sort in descending order
comparaciones1 <- comparaciones1 %>%
  rowwise() %>%
  mutate(sum_of_1s = sum(c_across(-Pathway))) %>%
  ungroup() %>%
  arrange(desc(sum_of_1s)) %>%
  select(-sum_of_1s)

# Optionally, set Pathway as rownames if you prefer a data frame that way:
comparaciones1 <- as.data.frame(comparaciones1)
rownames(comparaciones1) <- comparaciones1$Pathway
comparaciones1$Pathway <- NULL


# Save the summary matrix to an RData file
save(comparaciones1, file = paste0(output_dir, "comparaciones.RData"))

# End of Script

