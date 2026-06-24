library(dplyr)

# 1. Load your all-features SHAP summary file
shap_summary_path <- "model_results/endo2emb/10setseeds_endo2emb_rf_shap_importance_summary_2026-06-14_12h58m05s.csv"
shap_data <- read.csv(shap_summary_path)

# 2. Extract the top N (e.g., 20) features per dataset based on mean_shap
N_features <- 10
top_shap_list <- shap_data %>%
  group_by(dataset) %>%
  slice_max(order_by = mean_shap, n = N_features, with_ties = FALSE) %>%
  select(feature, dataset) # Keep only the columns the pipeline expects

# 3. Save it to your feature_lists folder
output_path <- file.path("feature_lists", paste0("top_", N_features, "_shap_features_endo2emb.csv"))
write.csv(top_shap_list, output_path, row.names = FALSE)

message("Saved top ", N_features, " SHAP features to: ", output_path)


# ===== ALTERNATIVE BY FEATIMPORTANCE =====
# 1. Load your all-features feature importance summary file
feat_summary_path <- "model_results/endo2emb/10setseeds_endo2emb_rf_feature_importance_summary_2026-06-14_12h58m05s.csv"
feat_data <- read.csv(feat_summary_path)

# 2. Extract the top N (e.g., 20) features per dataset based on mean_shap
N_features <- 10
top_feat_list <- feat_data %>%
  group_by(dataset) %>%
  slice_max(order_by = mean_importance, n = N_features, with_ties = FALSE) %>%
  select(feature, dataset) # Keep only the columns the pipeline expects

# 3. Save it to your feature_lists folder
output_path <- file.path("feature_lists", paste0("top_", N_features, "_featImp_features_endo2emb.csv"))
write.csv(top_feat_list, output_path, row.names = FALSE)

message("Saved top ", N_features, " featImportance features to: ", output_path)

