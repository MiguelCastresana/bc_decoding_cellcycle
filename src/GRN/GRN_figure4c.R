# Load required libraries
library(pheatmap)
library(tidyverse)

# ─── Load Data ─────────────────────────────────────────────────────────────
# Load discovery data
# Create a temporary environment and load the file into it
tmp <- new.env()
load("~/results/results_network/scenic_discovery_within/comparaciones", envir = tmp)

# Get the name of the first object in that environment and assign it to a variable
comparaciones_discovery_grn <- get(ls(tmp)[1], envir = tmp)


# Load validation data
tmp <- new.env()
load("~/results/results_network/scenic_validation_within/comparaciones", envir = tmp)
comparaciones_validation_grn <- get(ls(tmp)[1], envir = tmp)

# ─── Merge Row Names and Reorder Rows ─────────────────────────────────────────
# Get row names from both datasets and form union
names_discovery <- rownames(comparaciones_discovery_grn)
names_validation <- rownames(comparaciones_validation_grn)
all_names <- union(names_discovery, names_validation)
common_names <- intersect(names_discovery, names_validation)
unique_names <- setdiff(all_names, common_names)
ordered_names <- c(common_names, unique_names)

# Function: Reorder rows and add missing ones (filled with NA)
reorder_add_missing_rows <- function(df, ordered_names) {
  missing_rows <- setdiff(ordered_names, rownames(df))
  df_existing <- df[ordered_names[ordered_names %in% rownames(df)], , drop = FALSE]
  if (length(missing_rows) > 0) {
    df_missing <- matrix(NA, nrow = length(missing_rows), ncol = ncol(df),
                         dimnames = list(missing_rows, colnames(df)))
    df_existing <- rbind(df_existing, df_missing)
  }
  df_existing[match(ordered_names, rownames(df_existing)), , drop = FALSE]
}

comparaciones_discovery_grn <- reorder_add_missing_rows(comparaciones_discovery_grn, ordered_names)
comparaciones_validation_grn <- reorder_add_missing_rows(comparaciones_validation_grn, ordered_names)

# Check row names
print(rownames(comparaciones_discovery_grn))
print(rownames(comparaciones_validation_grn))

# Replace NAs with 0 in both datasets
comparaciones_discovery_grn[is.na(comparaciones_discovery_grn)] <- 0
comparaciones_validation_grn[is.na(comparaciones_validation_grn)] <- 0

# ─── Create the Union/Intersection Matrix ─────────────────────────────────────
# Initialize result matrix with zeros (same dimnames as discovery dataset)
resultado <- matrix(0, nrow = nrow(comparaciones_discovery_grn), 
                    ncol = ncol(comparaciones_discovery_grn),
                    dimnames = dimnames(comparaciones_discovery_grn))

# Fill in matrix: 1 if either dataset has a 1, then overwrite cells to 2 if both have 1
resultado[comparaciones_discovery_grn == 1 | comparaciones_validation_grn == 1] <- 1
resultado[comparaciones_discovery_grn == 1 & comparaciones_validation_grn == 1] <- 2

# ─── Adjust Colors and Column Names ───────────────────────────────────────────
color_mapping <- c("white", "gray", "green2")

# Replace "G1" with "G0/G1" and "G2M" with "G2/M"
colnames(resultado) <- colnames(resultado) %>%
  str_replace_all("G1", "G0/G1") %>%
  str_replace_all("G2M", "G2/M")

# ─── Define Sample Annotation ─────────────────────────────────────────────────
# Create a data frame for sample colors; here assume 4 subgroups each with 3 columns
my_sample_col <- data.frame(Subgroup = rep(c("LumA", "LumB", "Her2", "Basal"), each = 3),
                            row.names = colnames(resultado))
my_colour <- list(
  Subgroup = c(LumA = "#076fffff", LumB = "#0cd9faff", Her2 = "#f343d3ff", Basal = "#e1041eff")
)

# ─── Order Rows by Total Count of Intersection (value = 2) ─────────────────────
total_counts <- rowSums(resultado == 2)
resultado <- resultado[order(total_counts, decreasing = TRUE), ]

# ─── Plot the Heatmap ──────────────────────────────────────────────────────────
p <- pheatmap(resultado, 
              cluster_cols = FALSE, 
              cluster_rows = FALSE,
              annotation_col = my_sample_col,
              annotation_colors = my_colour,
              color = color_mapping, 
              fontsize = 9, 
              angle_col = 45,
              main = "Within subtypes", 
              legend = FALSE,
              gaps_col = cumsum(rep(3, 4)),
              annotation_legend = FALSE,
              cellwidth = 13)

ggsave("~/figures/Figure4_GRNs_within.pdf",p, width = 13, height = 8, dpi = 500)

