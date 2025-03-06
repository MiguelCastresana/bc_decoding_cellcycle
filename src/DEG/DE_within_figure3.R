

# -------------------------
# Load required libraries
# -------------------------

library(tidyverse)
require(pheatmap)

'%!in%' <- function(x, y)!('%in%'(x, y))




# ─── Load Data ─────────────────────────────────────────────────────────────
# Load discovery data
# Create a temporary environment and load the file into it
tmp <- new.env()
load("~/results/results_DEG/discovery/within/binary_limma_discovery_within.RData", envir = tmp)


# Get the name of the first object in that environment and assign it to a variable
comparison_discovery <- get(ls(tmp)[1], envir = tmp)


# validation

# Create a temporary environment and load the file into it
tmp <- new.env()
load("~/results/results_DEG/validation//within/binary_limma_validation_within.RData", envir = tmp)

comparison_validation <- get(ls(tmp)[1], envir = tmp)




# -------------------------
# Determine ordered row names:
# Common names first, then unique names.
# -------------------------
names_discovery <- rownames(comparison_discovery)
names_validation <- rownames(comparison_validation)
common_names <- intersect(names_discovery, names_validation)
unique_names <- setdiff(union(names_discovery, names_validation), common_names)
ordered_names <- c(common_names, unique_names)

# -------------------------
# Tidyverse function: Reorder matrix rows and add missing rows
# -------------------------
reorder_add_missing_rows <- function(mat, ordered_names) {
  mat %>%
    as_tibble(rownames = "row") %>%
    complete(row = ordered_names) %>%      # Ensure all ordered_names appear
    pivot_longer(-row, names_to = "col", values_to = "value") %>%
    pivot_wider(names_from = col, values_from = value) %>%
    arrange(factor(row, levels = ordered_names)) %>%
    { 
      m <- as.matrix(select(., -row))
      rownames(m) <- .$row
      m
    }
}

# Reorder both matrices
comparison_discovery <- reorder_add_missing_rows(comparison_discovery, ordered_names)
comparison_validation  <- reorder_add_missing_rows(comparison_validation, ordered_names)

# -------------------------
# Replace NA values with 0 using tidyverse
# -------------------------
# Convert NA values to 0 using tidyverse (we convert to tibble and back to matrix)
fix_na <- function(mat) {
  mat %>%
    as_tibble(rownames = "row") %>%
    dplyr::mutate(across(-row, ~replace_na(.x, 0))) %>%
    pivot_longer(-row, names_to = "col", values_to = "value") %>%
    pivot_wider(names_from = col, values_from = value) %>%
    arrange(match(row, rownames(mat))) %>%
    { as.matrix(select(., -row)) -> m; 
      rownames(m) <- .$row; m }
}


comparison_discovery <- fix_na(comparison_discovery)
comparison_validation  <- fix_na(comparison_validation)

colnames(comparison_validation) = colnames(comparison_discovery)

# -------------------------
# Create the results matrix based on directionality
# -------------------------
# For each cell: if one matrix shows a hit (abs > 0) and the other is 0 → result 1;
# if both are positive → result 2;
# if both are negative → result -2; else 0.
results <- left_join(
  comparison_discovery %>% as_tibble(rownames = "row") %>% pivot_longer(-row, names_to = "col", values_to = "disc"),
  comparison_validation %>% as_tibble(rownames = "row") %>% pivot_longer(-row, names_to = "col", values_to = "val"),
  by = c("row", "col")
) %>%
  mutate(result = case_when(
    (abs(disc) > 0 & val == 0) | (disc == 0 & abs(val) > 0) ~ 1,
    (disc > 0 & val > 0) ~ 2,
    (disc < 0 & val < 0) ~ -2,
    TRUE ~ 0
  )) %>%
  select(row, col, result) %>%
  pivot_wider(names_from = col, values_from = result) %>%
  column_to_rownames("row") %>%
  as.matrix()

# Replace any remaining NA with 0
results[is.na(results)] <- 0



# Adjust the color mapping to reflect directionality
color_mapping <- c("red", "white","grey", "green2")

# Adjust the breaks to ensure proper color mapping
breaks <- c(-2, -1, 0, 1, 2)




# Extracting the phase and subtype from each column name
subtypes_phases <- strsplit(colnames(results), " ")
subtypes <- sapply(subtypes_phases, function(x) gsub("Other", "", x[3]))
phases <- rep(c("G0/G1", "S", "G2/M"),4)

# Define the desired order of subtypes and phases
subtype_order <- c("LumA", "LumB", "Her2", "Basal")
phase_order <- c("G0/G1", "S", "G2/M")

# Create a data frame for ordering
ordering_df <- data.frame(column = colnames(comparison_discovery), subtype = subtypes, phase = phases)

# Order by subtype first, then by phase
ordering_df <- ordering_df[order(match(ordering_df$subtype, subtype_order), match(ordering_df$phase, phase_order)),]

# Reorder the columns of your matrix
comparaciones_ordered <- comparison_discovery[, ordering_df$column]



# Define sample colors
my_sample_col <- data.frame(Subgroup = rep(c("LumA", "LumB", "Her2", "Basal"), each = 3))

colnames(results) <- gsub("G1", "G0/G1", colnames(results))
colnames(results) <- gsub("G2M", "G2/M", colnames(results))

row.names(my_sample_col) <- colnames(results)

my_colour = list(
  Subgroup = c(LumA = "#076fffff", LumB = "#0cd9faff", `Her2` = "#f343d3ff", Basal = "#e1041eff")
)

# Count occurrences of 2 and -2 in each row
count_2 <- rowSums(results == 2)
count_neg2 <- rowSums(results == -2)

# Create a new column that adds the counts
total_counts <- count_2 + count_neg2

# Order rows based on total counts in decreasing order
results <- results[order(total_counts, decreasing = TRUE), ]


# Plot the heatmap
p <- pheatmap(results, cluster_cols = FALSE, cluster_rows = FALSE,
              annotation_colors = my_colour, annotation_col = my_sample_col,
              color = color_mapping,  fontsize = 6, 
              fontsize_row = 9, fontsize_col = 9,
              angle_col = "45",annotation_legend = FALSE,
              legend = FALSE,gaps_col = cumsum(c(3,3,3,3)),cellwidth = 13)




ggsave("~/figures//Figure3_DE_within.pdf",p, width=6, height=8,dpi = 500)





