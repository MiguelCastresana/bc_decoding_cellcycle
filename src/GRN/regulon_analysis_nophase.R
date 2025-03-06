
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
library(yaml)

'%!in%' <- function(x, y) !('%in%'(x, y))

# -------------------------
# Parse command-line arguments via YAML
# -------------------------
config_file_path <- "~/src/repository/repository_final/GRN/config_nophase.yaml"
config <- yaml::read_yaml(config_file_path)$default

validate_parameters <- function(config) {
  required_params <- c( "z_threshold", "relative_activity_cutoff", 
                       "min_qvalue", "loom_base_path", "metadata_path", 
                       "final_genes_path", "output_dir")
  missing_params <- setdiff(required_params, names(config))
  if (length(missing_params) > 0) {
    stop(paste("❌ Missing parameters in YAML config:", paste(missing_params, collapse = ", ")))
  }
  numeric_params <- c("z_threshold", "relative_activity_cutoff", "min_qvalue")
  numeric_limits <- list(z_threshold = c(0, Inf),
                         relative_activity_cutoff = c(0, 1),
                         min_qvalue = c(0, 1))
  for (param in numeric_params) {
    value <- config[[param]]
    if (!is.numeric(value) || length(value) != 1 ||
        value < numeric_limits[[param]][1] || value > numeric_limits[[param]][2]) {
      stop(paste("❌ Error:", param, "must be a numeric value in range", numeric_limits[[param]][1], "to", numeric_limits[[param]][2]))
    }
  }
  if (!dir.exists(config$loom_base_path)) stop(paste("❌ Loom base directory does not exist:", config$loom_base_path))
  if (!file.exists(config$metadata_path)) stop(paste("❌ Metadata file does not exist:", config$metadata_path))
  if (!file.exists(config$final_genes_path)) stop(paste("❌ Final genes file does not exist:", config$final_genes_path))
  if (!dir.exists(config$output_dir)) {
    message("Output directory does not exist. Creating: ", config$output_dir)
    dir.create(config$output_dir, recursive = TRUE)
  }
  message("✅ All configuration checks passed successfully.")
  return(config)
}

safe_validate_parameters <- function(config) {
  tryCatch(
    {
      return(validate_parameters(config))
    },
    error = function(e) {
      message("❌ Caught an error in validation: ", e$message)
      return(NULL)
    }
  )
}

config <- safe_validate_parameters(config)

# Assign validated parameters
z_threshold <- config$z_threshold
relative_activity_cutoff <- config$relative_activity_cutoff
min_qvalue <- config$min_qvalue
loom_base_path <- config$loom_base_path
metadata_path <- config$metadata_path
final_genes_path <- config$final_genes_path
output_dir <- config$output_dir

# -------------------------
# Process Loom Files Function
# -------------------------
process_loom <- function(loom_path) {
  
  # loom_path =current_loom_path
  loom <- open_loom(loom_path)
  
  # Extract regulon information
  regulons_incidMat <- get_regulons(loom, column.attr.name = "Regulons")
  regulons <- regulonsToGeneLists(regulons_incidMat) %>%
    enframe(name = "TranscriptionFactor", value = "Gene") %>%
    unnest(Gene)
  
  # Get AUCell results and related info
  regulonAUC <- get_regulons_AUC(loom, column.attr.name = "RegulonsAUC")
  total_cells <- dim(regulonAUC@colData)[1]
  regulonAucThresholds <- get_regulon_thresholds(loom)
  embeddings <- get_embeddings(loom)
  
  # Load metadata and filter by current group using analysis_var
  meta <- read_csv(metadata_path, col_types = cols()) %>%
    rename_all(tolower) %>%
    mutate(subtype = str_replace_all(subtype, "_SC", ""))
  
  cellClusters <- tibble(cell = colnames(regulonAUC)) %>%
    left_join(meta %>% select(cell, subtype), by = "cell") %>%
    mutate(clusters = factor(subtype, levels = c("LumA", "LumB", "Her2", "Basal"))) %>%
    select(cell, clusters) %>%
    column_to_rownames("cell")
  
  close_loom(loom)
  
  # Remove duplicate cells from regulonAUC
  regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)), ]
  
  # Split cells by cluster and calculate average regulon activity
  cellsPerCluster <- split(rownames(cellClusters), cellClusters$clusters)
  regulonActivity_bycellType <- sapply(cellsPerCluster,
                                       function(cells) rowMeans(getAUC(regulonAUC)[, cells]))
  
  # Dynamic topRegulators using melt and renaming with dynamic columns
  topRegulators <- reshape2::melt(regulonActivity_bycellType) %>%
    mutate(geneset = paste(Var2, "__", Var1, sep = "")) %>%
    rename(Regulon = Var1,
           subtype := Var2,
           RelativeActivity = value) %>%
    arrange(desc(RelativeActivity)) %>%
    ungroup()
  
  # Calculate RSS and reshape dynamically using pivot_longer
  selectedResolution <- "clusters"
  rss <- calcRSS(AUC = getAUC(regulonAUC),
                 cellAnnotation = cellClusters[colnames(regulonAUC), selectedResolution])
  rssNorm <- scale(rss)
  rssNorm[rssNorm < 0] <- 0
  rssSubset <- rssNorm
  rssSubset[rssSubset < z_threshold] <- 0
  
  rss.df <- as.data.frame(rss) %>%
    rownames_to_column("Topic") %>%
    pivot_longer(cols = -Topic,
                 names_to = "subtype",
                 values_to = "RSS") %>%
    mutate(geneset = str_c(.data[["subtype"]], "__", Topic))
  
  rssNorm.df <- as.data.frame(rssNorm) %>%
    rownames_to_column("Topic") %>%
    pivot_longer(cols = -Topic,
                 names_to = analysis_comparison,
                 values_to = "Z") %>%
    mutate(geneset = str_c(.data[["subtype"]], "__", Topic))
  
  sub_list <- left_join(rss.df, rssNorm.df %>% select(geneset, Z), by = "geneset")
  
  return(list(regulons = regulons, 
              regulonAUC = regulonAUC, 
              cellClusters = cellClusters, 
              topRegulators = topRegulators,
              total_cells = total_cells,
              rss = sub_list))
}

# -------------------------
# Helper: Convert DataFrame to Named Matrix
# -------------------------
dfToNamedMatrix <- function(df) {
  return(as.matrix(df))
}

# -------------------------
# Main Pipeline: Process all groups
# -------------------------


current_loom_path <- file.path(loom_base_path, "pbmc10k_scenic_integrated-output.loom")
res <- process_loom(current_loom_path)

all_topRegulators <- res$topRegulators
all_regulons_list <- res$regulons
all_total_cells <- res$total_cells
rss_list <- res$rss



regulons_list_df <- bind_rows(all_regulons_list)

topRegulators <- bind_rows(all_topRegulators)
rss_list <- bind_rows(rss_list)

write.table(rss_list, file = paste0(output_dir, "RSS_results.csv"),
            sep = ";", row.names = FALSE, quote = FALSE)

# -------------------------
# Gene Mapping
# -------------------------
ensembl <- mapIds(org.Hs.eg.db, keys = regulons_list_df$Gene,
                  column = "ENSEMBL", keytype = "SYMBOL", multiVals = "first")
regulons_list_df$ensembl <- ensembl



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


regulons_list_new <- regulons_list %>%
  inner_join(combined_mappings, by = "ensembl") %>%  # merge by "ensembl"
  select(entrez, geneset) %>%                        # select only the desired columns
  filter(complete.cases(.)) %>%                      # remove rows with NA values
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

result_list <- list()
groups <- unique(regulons_list_new$geneset)
c <- 1
for(g in groups) {
  sub <- regulons_list_new %>% filter(geneset == g)
  msig_ora <- clusterProfiler::enricher(gene = sub$entrez,
                                        universe = as.vector(total_genes),
                                        TERM2GENE = msig_h)
  if (!is.null(msig_ora)) {
    results <- msig_ora@result[, c("ID", "GeneRatio", "pvalue", "p.adjust")]
    results$geneset <- rep(g, nrow(results))
    result_list[[c]] <- results
    c <- c + 1
  }
  message("Processing geneset: ", g)
}

result_list <- do.call(rbind.data.frame, result_list)
gea_output <- result_list
rownames(gea_output) <- NULL
gea_output <- gea_output[, c("geneset", "ID", "GeneRatio", "pvalue", "p.adjust")]
names(gea_output) <- c("geneset", "Pathway", "Overlap", "p-value", "q-value")
gea_output <- gea_output %>% mutate(FWER = p.adjust(`p-value`, method = "bonferroni"))

# -------------------------
# Fusion of Top Regulators with RSS Results
# -------------------------

# Fuse top regulators with RSS list
fus_list <- topRegulators %>% 
  select(Regulon, geneset, RelativeActivity) %>% 
  left_join(
    rss_list %>% 
      filter(Z > z_threshold) %>% 
      select(geneset, subtype, RSS, Z), 
    by = "geneset"
  ) %>% 
  mutate(RSS = as.numeric(RSS)) %>% 
  drop_na() %>% 
  mutate(subtype = factor(subtype, levels = c("LumA", "LumB", "Her2", "Basal")))


# Get unique Regulon names from the fusion result
unique_reg <- fus_list %>% 
  pull(Regulon) %>% 
  unique()

# Filter gea_output by those regulators (assuming gea_output has a 'geneset' column)
gea_output <- gea_output %>% 
  filter(geneset %in% topRegulators$Regulon)

# Further filter gea_output by q-value threshold and by unique regulators,
# then rename the first column (geneset) to "Regulon"
gea_output1 <- gea_output %>% 
  filter(`q-value` < 0.05, geneset %in% unique_reg) %>% 
  rename(Regulon = geneset)

# Join the filtered enrichment output with the fusion list on "Regulon"
final <- left_join(gea_output1, fus_list, by = "Regulon",relationship = "many-to-many")

# Build the analysis data frame:
# - Rename column 2 to "Pathway" (if not already)
# - Create a 'both' column from Subtype and phase
# - Clean up the 'Pathway' names
# - Filter rows based on RelativeActivity cutoff
analysis <- final %>% 
  mutate(Pathway = as.character(Pathway),  # ensure Pathway is character
         Pathway = gsub("_", " ", Pathway),
         Pathway = sub("HALLMARK ", "", Pathway)) %>% 
  filter(RelativeActivity >= cutoff)

# Define the unique set of subtypes (columns) and the unique set of pathways (rows)
unique_cols <- c("LumA", "LumB", "Her2", "Basal")
unique_paths <- analysis %>% 
  pull(Pathway) %>% 
  unique()

# Create a summary matrix using tidyverse
comparaciones1 <- analysis %>%
  # Remove duplicate entries based on Pathway and Subtype
  distinct(Pathway, subtype, .keep_all = TRUE) %>%
  # Create a binary indicator for significance (q-value < 0.05)
  mutate(sig = as.integer(`q-value` < 0.05)) %>%
  select(Pathway, subtype, sig) %>%
  # Pivot so that each Subtype becomes a column and Pathway stays as rows
  pivot_wider(names_from = subtype, values_from = sig, values_fill = list(sig = 0)) %>%
  # Ensure all pathways are present (fill missing rows with zeros)
  complete(Pathway = unique_paths, fill = list(sig = 0)) %>%
  arrange(factor(Pathway, levels = unique_paths)) %>%
  # Set Pathway as row names and order columns as defined in unique_cols
  column_to_rownames("Pathway") %>%
  select(any_of(unique_cols)) %>%
  as.data.frame()

save(comparaciones1, file = paste0(output_dir, "comparaciones.RData"))

# End of Script
