# ===== IMPORTS =====
# install.packages("dplyr")
library(dplyr)

# =========================================================================
# CONFIGURATION
# =========================================================================

# 1. Define the list of summary files you want to analyze.
# You can mix base feature importance and SHAP importance files.
summary_files <- c(
  "model_results/10setseeds_emb2endo_rf_feature_importance_summary_2026-06-14_11h54m09s.csv",
  "model_results/10setseeds_emb2endo_rf_shap_importance_summary_2026-06-14_11h54m09s.csv"
  # Add as many paths as you want to compare here
)

# 2. Dynamic Variables for Feature Selection
top_n_features <- 20          # How many top features to extract per dataset, per file?
min_file_appearances <- 2     # Out of the files listed above, how many times must it appear to be kept?

# =========================================================================
# PROCESS & AGGREGATE
# =========================================================================

extract_top_features <- function(file_path, top_n) {
  # Read the CSV
  df <- read.csv(file_path, stringsAsFactors = FALSE)
  
  # Standardize column names dynamically based on whether it's SHAP or Base Importance
  if ("mean_shap" %in% colnames(df)) {
    df <- df %>% rename(score = mean_shap)
  } else if ("mean_importance" %in% colnames(df)) {
    df <- df %>% rename(score = mean_importance)
  } else {
    stop(paste("Could not identify the importance score column in:", file_path))
  }
  
  # Group by vivo/vitro, rank them by score, and pull the top N
  top_df <- df %>%
    group_by(dataset) %>%
    slice_max(order_by = score, n = top_n, with_ties = FALSE) %>%
    mutate(
      source_file = basename(file_path),
      rank_in_file = row_number()
    ) %>%
    ungroup() %>%
    select(dataset, feature, score, rank_in_file, source_file)
  
  return(top_df)
}

# 1. Read and extract top features from all files in the vector
all_top_features_list <- lapply(summary_files, function(file) {
  if (file.exists(file)) {
    extract_top_features(file, top_n_features)
  } else {
    warning(paste("File not found and skipped:", file))
    NULL
  }
})

# Combine into one master dataframe
all_top_features <- bind_rows(all_top_features_list)

# 2. Find the CONSISTENT features
consensus_features <- all_top_features %>%
  group_by(dataset, feature) %>%
  summarise(
    appearance_count = n(),
    avg_rank = mean(rank_in_file),
    avg_score = mean(score),
    .groups = "drop"
  ) %>%
  # Filter by our consistency threshold
  filter(appearance_count >= min_file_appearances) %>%
  arrange(dataset, desc(appearance_count), avg_rank)

# =========================================================================
# SAVE RESULTS
# =========================================================================

# Ensure output directory exists
if(!dir.exists("feature_lists")) dir.create("feature_lists") 

timestamp <- format(Sys.time(), "%Y-%m-%d_%Hh%Mm%Ss")
save_path <- paste0("feature_lists/consensus_top_", top_n_features, "_features_", timestamp, ".csv")

write.csv(consensus_features, save_path, row.names = FALSE)

# Print summary to console
print(paste("Analysis complete. Saved consensus features to:", save_path))
print("Overview of kept features per dataset:")
print(table(consensus_features$dataset))

# View the top ones in the console
print(head(consensus_features, 15))
