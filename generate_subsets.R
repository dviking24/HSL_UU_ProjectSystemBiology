#install.packages("dplyr")
library(dplyr)

# =========================================================================
# DATA PREPARATION & SAVING
# Run this once to generate the 3 TSV formats (for both vivo & vitro).
# =========================================================================

# ----- 1A. Load Raw Data -----
vivo_raw <- read.csv("matrices/vivo_matrix.tsv", header = TRUE, sep = "\t", row.names = 1)
vitro_raw <- read.csv("matrices/vitro_matrix.tsv", header = TRUE, sep = "\t", row.names = 1)

# Replace . with |
colnames(vivo_raw)  <- gsub("\\.", "|", colnames(vivo_raw))
colnames(vitro_raw) <- gsub("\\.", "|", colnames(vitro_raw))

# ----- 1B. Filtering Helpers -----
is_embryo <- function(x) grepl("^(Zo|SW)_", x)
is_endo   <- function(x) grepl("^(NP|PR)_", x)

# Function to filter matrix based on direction
filter_combinations <- function(mat, direction = c("both", "emb2endo", "endo2emb")) {
  direction <- match.arg(direction)
  if (direction == "both") return(mat)
  
  split_cols <- strsplit(colnames(mat), "\\|")
  left_part  <- sapply(split_cols, `[`, 1)
  right_part <- sapply(split_cols, `[`, 2)
  
  if (direction == "emb2endo") {
    keep_idx <- is_embryo(left_part) & is_endo(right_part)
  } else if (direction == "endo2emb") {
    keep_idx <- is_endo(left_part) & is_embryo(right_part)
  }
  
  return(mat[, keep_idx, drop = FALSE])
}

# ----- 1C. Processing and Saving Function -----
process_and_save_matrices <- function(raw_mat, prefix) {
  directions <- c("both", "emb2endo", "endo2emb")
  
  for (dir in directions) {
    # Filter columns
    filtered_mat <- filter_combinations(raw_mat, direction = dir)
    
    # Transpose data.frame
    df <- as.data.frame(t(filtered_mat))
    
    # Add PregnancyStatus
    df <- cbind(
      PregnancyStatus = factor(
        ifelse(grepl("NP", rownames(df)), "NOT_PREGNANT", "PREGNANT")
      ),
      df
    )
    
    # Save to TSV
    filename <- paste0("matrices/", prefix, "_matrix_", dir, ".tsv")
    write.table(df, file = filename, sep = "\t", col.names = NA, quote = FALSE)
    message("Saved: ", filename)
  }
}

# Generate and save all 6 datasets (3 vivo, 3 vitro)
process_and_save_matrices(vivo_raw, "vivo")
process_and_save_matrices(vitro_raw, "vitro")
