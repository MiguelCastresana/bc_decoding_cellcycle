

# Load across phase intersection between discovery and validation
load("~/results/results_DEG/intersection/phases")

colnames(results) <- gsub("G1", "G0/G1", colnames(results))  # Replace "G1" with "G0/G1"
colnames(results) <- gsub("G2M", "G2/M", colnames(results))  # Replace "G2M" with "G2/M"


# Define the color mapping (colors should match the number of breaks - 1)
color_mapping <- c("red", "white", "grey", "green2")

# Define the breaks to properly map to the color scale
breaks <- c(-2, -0.5, 0.5, 1.5, 2)



my_sample_col <- data.frame(phase = rep(c("G0/G1", "S","G2/M"), each = 4))
row.names(my_sample_col) <- colnames(results)

my_colour = list(
  phase = c('G0/G1' = "#CD6600", S = "#4fa601ff",'G2/M' = "#d4aa00ff")

)
# colnames(results) = trimws(colnames(results),"left")
 valores = sub(" .*", "", colnames(results))
# row.names(my_sample_col) = colnames(results)
 
# Count occurrences of 2 and -2 in each row
count_2 <- rowSums(results == 2)
count_neg2 <- rowSums(results == -2)

# Create a new column that adds the counts
total_counts <- count_2 + count_neg2

# Order rows based on total counts in decreasing order
results <- results[order(total_counts, decreasing = TRUE), ]


# Plot the heatmap
p <- pheatmap(results, 
              cluster_cols = FALSE, 
              cluster_rows = FALSE,
              annotation_colors = my_colour, 
              annotation_col = my_sample_col,
              color = color_mapping, 
              breaks = breaks, 
              fontsize = 8, 
              fontsize_row = 9, fontsize_col = 9,
              angle_col = "45", 
              legend = FALSE,
              annotation_legend = FALSE,
              gaps_col = cumsum(c(4,4,4)),
              cellwidth = 13)  # Adjust to make boxes shorter



ggsave("~/figures//Figure3_DE_across.pdf",p, width=6, height=8,dpi = 500)



