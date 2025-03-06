# Load required libraries
library(pheatmap)
library(tidyverse)

# ─── Load Data ─────────────────────────────────────────────────────────────
# Load discovery data
# Create a temporary environment and load the file into it
tmp <- new.env()
load("~/results/results_network/scenic_discovery_across/comparaciones", envir = tmp)

# Get the name of the first object in that environment and assign it to a variable
comparaciones_discovery_grn <- get(ls(tmp)[1], envir = tmp)


# Load validation data
tmp <- new.env()
load("~/results/results_network/scenic_validation_across/comparaciones", envir = tmp)
# Get the name of the first object in that environment and assign it to a variable
comparaciones_validation_grn <- get(ls(tmp)[1], envir = tmp)

# ─── Merge and Reorder Row Names ─────────────────────────────────────────────
# Get row names and create a union of names
names_discovery <- rownames(comparaciones_discovery_grn)
names_validation <- rownames(comparaciones_validation_grn)
ordered_names <- c(intersect(names_discovery, names_validation), setdiff(union(names_discovery, names_validation), intersect(names_discovery, names_validation)))

# Function to reorder and add missing rows (filled with NA)
reorder_add_missing_rows <- function(df, ordered_names) {
  df_existing <- df[rownames(df) %in% ordered_names, , drop = FALSE]
  missing_rows <- setdiff(ordered_names, rownames(df))
  if (length(missing_rows) > 0) {
    df_missing <- matrix(NA, nrow = length(missing_rows), ncol = ncol(df),
                         dimnames = list(missing_rows, colnames(df)))
    df_existing <- rbind(df_existing, df_missing)
  }
  df_existing[match(ordered_names, rownames(df_existing)), , drop = FALSE]
}

comparaciones_discovery_grn <- reorder_add_missing_rows(comparaciones_discovery_grn, ordered_names)
comparaciones_validation_grn <- reorder_add_missing_rows(comparaciones_validation_grn, ordered_names)

# Replace NA values with 0
comparaciones_discovery_grn[is.na(comparaciones_discovery_grn)] <- 0
comparaciones_validation_grn[is.na(comparaciones_validation_grn)] <- 0

# ─── Create Union/Intersection Matrix ─────────────────────────────────────────
# Initialize result matrix with zeros (using discovery dimnames)
resultado <- matrix(0, nrow = nrow(comparaciones_discovery_grn), ncol = ncol(comparaciones_discovery_grn),
                    dimnames = dimnames(comparaciones_discovery_grn))

# Set cell value to 1 if either dataset has 1; overwrite with 2 if both have 1.
resultado[(comparaciones_discovery_grn == 1 | comparaciones_validation_grn == 1)] <- 1
resultado[(comparaciones_discovery_grn == 1 & comparaciones_validation_grn == 1)] <- 2

# ─── Adjust Column Names and Define Breaks ─────────────────────────────────────
resultado <- resultado %>%
  { colnames(.) <- gsub("G1", "G0/G1", colnames(.)); . } %>%
  { colnames(.) <- gsub("G2M", "G2/M", colnames(.)); . }
color_mapping <- c("white", "gray", "green2")
breaks <- c(-0.5, 0.5, 1.5, 2.5)

# ─── Sample Annotation ───────────────────────────────────────────────────────
# For the across dataset, assume 3 columns per phase; adjust if needed.
sample_names <- colnames(resultado)
my_sample_col <- data.frame(phase = rep(c("G0/G1", "S", "G2/M"), each = 4))
rownames(my_sample_col) <- sample_names
my_colour <- list(
  phase = c("G0/G1" = "#d93b01ff", "S" = "#4fa601ff", "G2/M" = "#d4aa00ff")
)

# ─── Order Rows by Intersection Count ─────────────────────────────────────────
total_counts <- rowSums(resultado == 2)
resultado <- resultado[order(total_counts, decreasing = TRUE), ]

# ─── Plot and Save the Heatmap ────────────────────────────────────────────────
p <- pheatmap(resultado, 
              cluster_cols = FALSE, 
              cluster_rows = FALSE,
              annotation_col = my_sample_col,
              annotation_colors = my_colour,
              color = color_mapping, 
              fontsize = 9, 
              angle_col = 45,
              main = "Across subtypes", 
              legend = FALSE,
              gaps_col = cumsum(rep(4, 3)),  # Assuming 3 groups with 4 columns each
              cellwidth = 13,
              annotation_legend = FALSE)

ggsave("~/figures/Figure4_GRNs_across.pdf",p, width = 13, height = 8, dpi = 500)
