library(caret)
library(ranger)
library(treeshap)
library(shapviz)

# =========================================================================
# CONFIGURATION
# =========================================================================
# Select the specific variation, dataset, and seed model you want to dissect
TARGET_VARIATION <- "emb2endo"
TARGET_DATASET   <- "vivo"       # "vivo" or "vitro"
TARGET_SEED      <- 64

# Point to your auto-saved model directory (Update timestamps to match your folder!)
MODEL_PATH <- "saved_models/emb2endo_Top20SHAP_Only_emb2endo_subset_2026-06-14_18h41m15s/emb2endo_Top20SHAP_Only_emb2endo_subset_vivo_model_seed_64.rds"
MATRIX_PATH <- paste0("matrices/", TARGET_DATASET, "_matrix_", TARGET_VARIATION, ".tsv")

# =========================================================================
# LOAD AND PREPARE DATA
# =========================================================================
# Load Model and Data Matrix
rf_caret <- readRDS(MODEL_PATH)
matrix_df <- read.table(MATRIX_PATH, header = TRUE, sep = "\t", row.names = 1)

# Extract features used in training (excluding label)
x_data <- matrix_df[, colnames(matrix_df) != "PregnancyStatus", drop = FALSE]

# =========================================================================
# GENERATE DIRECTED SHAP VISUALIZATIONS
# =========================================================================
# Unify the model and compute raw directional SHAP values via treeshap
unified <- treeshap::unify(rf_caret$finalModel, x_data)
treeshap_res <- treeshap::treeshap(unified, x_data, verbose = FALSE)

# Convert to a shapviz object for elite plotting capabilities
shp <- shapviz(treeshap_res, X = x_data)

# 1. Global Beeswarm Summary Plot
message("--> Rendering Beeswarm Plot...")
sv_importance(shp, kind = "beeswarm", max_display = 15)

# 2. Dynamic Feature Dependence Plot
# Programmatically find the #1 feature with the highest mean absolute SHAP value
mean_shaps <- colMeans(abs(treeshap_res$shaps))
top_feature_name <- names(sort(mean_shaps, decreasing = TRUE))[1]

message("--> Code dynamically selected the top SHAP feature: ", top_feature_name)
message("--> Rendering Dependence Plot...")

# Render the dependence plot safely
sv_dependence(shp, v = top_feature_name, alpha = 0.8)
