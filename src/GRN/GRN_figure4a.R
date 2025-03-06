library(dplyr)
library(tidyr)
library(purrr)
library(ggvenn)
library(eulerr)
library(ggplot2)

# ─── Load the phase and no-phase datasets ────────────────────────────────
tmp <- new.env()
load("~/results/results_network/intersection/GRN_phase_rss_1_5_auc_005", envir = tmp)

# Get the name of the first object in that environment and assign it to a variable
comparaciones_phase <- get(ls(tmp)[1], envir = tmp)

tmp <- new.env()
load("~/results/results_network/intersection/GRN_nophase_rss_1_5_auc_005", envir = tmp)

# Get the name of the first object in that environment and assign it to a variable
comparaciones_nophase <- get(ls(tmp)[1], envir = tmp)

# Add row names as a column and merge datasets
comparaciones_phase <- as.data.frame(comparaciones_phase) %>% 
  rownames_to_column("row_name")
comparaciones_nophase <- as.data.frame(comparaciones_nophase) %>% 
  rownames_to_column("row_name")

comparaciones_combined <- full_join(comparaciones_nophase, comparaciones_phase, by = "row_name") %>% 
  column_to_rownames("row_name")

# ─── Reorder columns based on expected groups ─────────────────────────────
target_order <- c("LumA", "LumA G1", "LumA S", "LumA G2M",
                  "LumB", "LumB G1", "LumB S", "LumB G2M",
                  "Her2", "Her2 G1", "Her2 S", "Her2 G2M",
                  "Basal", "Basal G1", "Basal S", "Basal G2M")
target_order <- intersect(target_order, colnames(comparaciones_combined))
comparaciones_final <- comparaciones_combined %>% 
  select(all_of(target_order)) %>% 
  mutate(across(everything(), ~ replace_na(.x, 0))) %>% 
  # Convert all 1's to 0 (if that is your intention)
  mutate(across(everything(), ~ if_else(.x == 1, 0, .x)))

# ─── Define functions for applying phase conditions ───────────────────────

# Function: For a set of 4 columns, if the first value is nonzero and the other three are zero, return 1;
# if the first value is zero and any of the other three are nonzero, return -1; else NA.
apply_conditions <- function(df) {
  df %>% 
    pmap_dbl(~ {
      row <- c(...)
      if(row[1] != 0 && all(row[2:4] == 0)) {
        1
      } else if(row[1] == 0 && any(row[2:4] != 0)) {
        -1
      } else {
        NA_real_
      }
    })
}

# For "NO PHASE": return indices where the (single-column) value is nonzero.
apply_conditions_nophase <- function(vec) {
  which(vec != 0)
}

# For "PHASE": return indices where any value in the provided columns is nonzero.
apply_conditions_phase <- function(df) {
  df %>% 
    mutate(any_nonzero = if_else(rowSums(across(everything(), ~ if_else(. != 0, 1, 0)) ) > 0, 1, 0)) %>% 
    pull(any_nonzero) %>% 
    { which(. == 1) }
}

# ─── Apply conditions for each subgroup ─────────────────────────────────────
# (Assuming comparaciones_final has 16 columns arranged as in target_order)

# Number of rows in comparaciones_final
n_rows <- nrow(comparaciones_final)

result1 <- comparaciones_final %>% select(1:4) %>% as_tibble() %>% mutate(cond = apply_conditions(cur_data_all())) %>% pull(cond)
result2 <- comparaciones_final %>% select(5:8) %>% as_tibble() %>% mutate(cond = apply_conditions(cur_data_all())) %>% pull(cond)
result3 <- comparaciones_final %>% select(9:12) %>% as_tibble() %>% mutate(cond = apply_conditions(cur_data_all())) %>% pull(cond)
result4 <- comparaciones_final %>% select(13:16) %>% as_tibble() %>% mutate(cond = apply_conditions(cur_data_all())) %>% pull(cond)

results <- tibble(
  Comparison1 = result1,
  Comparison2 = result2,
  Comparison3 = result3,
  Comparison4 = result4
)
print(results)

# ─── NO PHASE: Get indices of nonzero elements ─────────────────────────────
result1_np <- comparaciones_final %>% select(1) %>% pull() %>% { which(. != 0) }
result2_np <- comparaciones_final %>% select(5) %>% pull() %>% { which(. != 0) } + n_rows
result3_np <- comparaciones_final %>% select(9) %>% pull() %>% { which(. != 0) } + 2 * n_rows
result4_np <- comparaciones_final %>% select(13) %>% pull() %>% { which(. != 0) } + 3 * n_rows

nophase <- c(result1_np, result2_np, result3_np, result4_np)

# ─── PHASE: Get indices where any value is nonzero ───────────────────────────
result1_phase <- comparaciones_final %>% select(2:4) %>% apply_conditions_phase()
result2_phase <- comparaciones_final %>% select(6:8) %>% apply_conditions_phase() + n_rows
result3_phase <- comparaciones_final %>% select(10:12) %>% apply_conditions_phase() + 2 * n_rows
result4_phase <- comparaciones_final %>% select(14:16) %>% apply_conditions_phase() + 3 * n_rows

phase <- c(result1_phase, result2_phase, result3_phase, result4_phase)

# ─── Generate Venn Diagram ───────────────────────────────────────────────────
x <- list(A = nophase, B = phase)

ggvenn(x, 
       fill_color = c("#0073C2FF", "#EFC000FF"),
       stroke_size = 0.5, set_name_size = 2)

# ─── Generate Euler Diagram ───────────────────────────────────────────────────
fit <- euler(c("No phases" = 3, "With Phases" = 48, "No phases&With Phases" = 24))
plot(fit,
     fills = list(fill = c("yellow1", "dodgerblue", "#66CD00"), alpha = 0.7),
     labels = list(col = "black", font = 4, cex = 1.6),
     quantities = list(cex = 2))
