# =========================================================================
# ===== 1. IMPORTS =====
# =========================================================================
library(caret)
library(ranger)
library(dplyr)
library(treeshap)
library(pROC)

# =========================================================================
# ===== 2. USER CONFIGURATION BLOCK =====
# =========================================================================
# Matrix directional type to load. Options: "both", "emb2endo", "endo2emb"
SELECTED_VARIATION <- "endo2emb" 

# Path to consensus feature list (.csv) to apply feature selection/trimming.
# Set this to NULL if you want to run the model on ALL features.
# FEATURE_LIST_PATH  <- NULL 
FEATURE_LIST_PATH <- "feature_lists/top_10_shap_features_endo2emb.csv"

# Output naming management
BASE_JOB_NAME      <- "endo2emb_Top10SHAP_Only"

# OPTIONAL: Save trained caret/ranger models to .rds files for future testing?
SAVE_MODELS        <- TRUE  # Set to FALSE to skip saving model objects

# Reproducibility settings
# Set to either the vector with the given seeds, or generate vector with 10 random seeds from range 1-1000
SEEDS_TO_RUN       <- c(64, 28, 21, 94, 41, 12, 53, 22, 17, 62)
# SEEDS_TO_RUN       <- sample.int(1000, 10)

# Experimental train/test sample configurations
VIVO_TRAIN_N       <- 7
VIVO_TEST_N        <- 4

VITRO_TRAIN_N      <- 5
VITRO_TEST_N       <- 4


# =========================================================================
# ===== 3. DATA LOADING & PREPROCESSING FUNCTIONS =====
# =========================================================================

#' Load Matrices and Enforce Factor Formatting
load_target_matrices <- function(variation) {
  vivo_path  <- file.path("matrices", paste0("vivo_matrix_", variation, ".tsv"))
  vitro_path <- file.path("matrices", paste0("vitro_matrix_", variation, ".tsv"))
  
  if (!file.exists(vivo_path) || !file.exists(vitro_path)) {
    stop("Target data matrices not found for variation: ", variation)
  }
  
  vivo  <- read.table(vivo_path, header = TRUE, sep = "\t", row.names = 1)
  vitro <- read.table(vitro_path, header = TRUE, sep = "\t", row.names = 1)
  
  vivo$PregnancyStatus  <- factor(vivo$PregnancyStatus)
  vitro$PregnancyStatus <- factor(vitro$PregnancyStatus)
  
  list(vivo = vivo, vitro = vitro)
}

#' Apply Dynamic Feature Filtering From Consensus Lists
apply_feature_selection <- function(matrices, feature_list_path) {
  if (is.null(feature_list_path) || !file.exists(feature_list_path)) {
    message("--> No valid feature subset provided. Running on ALL features.")
    return(list(matrices = matrices, mode_tag = "_all_"))
  }
  
  message("--> Loading consensus feature selection file from: ", feature_list_path)
  consensus_features <- read.csv(feature_list_path, stringsAsFactors = FALSE)
  
  # Trim Vivo Matrix
  vivo_feats <- consensus_features$feature[consensus_features$dataset == "vivo"]
  if (length(vivo_feats) > 0) {
    vivo_keep <- intersect(c("PregnancyStatus", vivo_feats), colnames(matrices$vivo))
    matrices$vivo <- matrices$vivo[, vivo_keep, drop = FALSE]
    message("    Trimmed vivo_matrix down to ", length(vivo_keep) - 1, " features.")
  }
  
  # Trim Vitro Matrix
  vitro_feats <- consensus_features$feature[consensus_features$dataset == "vitro"]
  if (length(vitro_feats) > 0) {
    vitro_keep <- intersect(c("PregnancyStatus", vitro_feats), colnames(matrices$vitro))
    matrices$vitro <- matrices$vitro[, vitro_keep, drop = FALSE]
    message("    Trimmed vitro_matrix down to ", length(vitro_keep) - 1, " features.")
  }
  
  return(list(matrices = matrices, mode_tag = "_subset_"))
}


# =========================================================================
# ===== 4. SAMPLE-ORDER AGNOSTIC SPLITTING ENGINE =====
# =========================================================================

#' Parse Entity Types Uniformly from Rownames Regardless of Pipe Sequence Order
extract_entities <- function(df) {
  rn <- rownames(df)
  parts <- strsplit(rn, "\\|")
  all_samples <- unique(unlist(parts))
  
  embryos <- all_samples[grepl("^(Zo|SW)_", all_samples)]
  endos   <- all_samples[grepl("^(NP|PR)_", all_samples)]
  
  list(embryos = embryos, endos = endos)
}

split_embryos <- function(embryos) {
  list(
    NP = embryos[grepl("_NP_", embryos)],
    PR = embryos[grepl("_PR_", embryos)]
  )
}

split_endos <- function(endos) {
  list(
    NP = endos[grepl("^NP_", endos)],
    PR = endos[grepl("^PR_", endos)]
  )
}

#' Draw Uniform Resamples Without Category Clashes
sample_entities <- function(embryos, endos, train_n, test_n, seed) {
  set.seed(seed)
  emb <- split_embryos(embryos)
  end <- split_endos(endos)
  
  train_emb_PR <- sample(emb$PR, train_n)
  train_emb_NP <- sample(emb$NP, train_n)
  
  test_emb_PR <- sample(setdiff(emb$PR, train_emb_PR), test_n)
  test_emb_NP <- sample(setdiff(emb$NP, train_emb_NP), test_n)
  
  train_end_PR <- sample(end$PR, train_n)
  train_end_NP <- sample(end$NP, train_n)
  
  test_end_PR <- sample(setdiff(end$PR, train_end_PR), test_n)
  test_end_NP <- sample(setdiff(end$NP, train_end_NP), test_n)
  
  list(
    train_embryos = c(train_emb_PR, train_emb_NP),
    test_embryos  = c(test_emb_PR, test_emb_NP),
    train_endos   = c(train_end_PR, train_end_NP),
    test_endos    = c(test_end_PR, test_end_NP)
  )
}

#' Assign Sample Pairs to Partitions Based on Complete Biological Matches
build_split <- function(df, split_info) {
  rn <- rownames(df)
  parts <- strsplit(rn, "\\|")
  
  train_idx <- sapply(parts, function(p) {
    any(p %in% split_info$train_embryos) & any(p %in% split_info$train_endos)
  })
  
  test_idx <- sapply(parts, function(p) {
    any(p %in% split_info$test_embryos) & any(p %in% split_info$test_endos)
  })
  
  list(
    train = df[train_idx, , drop = FALSE],
    test  = df[test_idx, , drop = FALSE]
  )
}


# =========================================================================
# ===== 5. CORE MODELING AND EVALUATING EXECUTIVE =====
# =========================================================================

#' Run Training, Predictions, ROC Analysis, and SHAP Unification For an Experiment
run_rf_experiment <- function(df, seed, dataset_name, train_n, test_n) { 
  set.seed(seed)
  
  # ---- Data Partitioning ----
  entities   <- extract_entities(df)
  split_info <- sample_entities(entities$embryos, entities$endos, train_n, test_n, seed)
  splits     <- build_split(df, split_info)
  
  train_df   <- splits$train
  test_df    <- splits$test
  
  # ---- Feature Matrix Formatting (The 0/1 SHAP Conversion Trick) ----
  x_train     <- train_df[, names(train_df) != "PregnancyStatus", drop = FALSE]
  y_train_num <- ifelse(train_df$PregnancyStatus == "PREGNANT", 1, 0)
  
  x_test        <- test_df[, names(test_df) != "PregnancyStatus", drop = FALSE]
  y_test_factor <- test_df$PregnancyStatus 
  
  # ---- Model Training (Muting the expected 0/1 caret regression warning) ----
  model <- suppressWarnings(
    train(
      x = x_train, y = y_train_num,
      method = "ranger",
      importance = "impurity"
    )
  )
  
  # ---- Predict Class Probabilities & Evaluate Metrics ----
  preds_prob   <- predict(model, x_test)
  preds_class  <- ifelse(preds_prob > 0.5, "PREGNANT", "NOT_PREGNANT")
  preds_factor <- factor(preds_class, levels = levels(y_test_factor))
  
  cm        <- confusionMatrix(preds_factor, y_test_factor)
  roc_obj   <- roc(response = y_test_factor, predictor = preds_prob, quiet = TRUE)
  auc_value <- as.numeric(auc(roc_obj))
  
  # ---- Gini / Impurity Feature Importance ----
  vi    <- varImp(model)$importance
  vi_df <- data.frame(
    feature    = rownames(vi),
    importance = vi[, 1],
    dataset    = dataset_name,
    seed       = seed
  )
  
  # ---- Unified SHAP Tree Calculation ----
  ranger_model  <- model$finalModel
  unified       <- treeshap::unify(ranger_model, x_train)
  shap_res      <- treeshap::treeshap(unified, x_test, verbose = FALSE)
  mean_abs_shap <- colMeans(abs(shap_res$shaps))
  
  shap_df <- data.frame(
    feature         = names(mean_abs_shap),
    shap_importance = unname(mean_abs_shap),
    dataset         = dataset_name,
    seed            = seed
  )
  
  # ---- Return Consolidated Metrics Packaging ----
  metrics <- data.frame(
    dataset  = dataset_name,
    seed     = seed,
    accuracy = cm$overall["Accuracy"],
    kappa    = cm$overall["Kappa"],
    auc      = auc_value
  )
  
  list(metrics = metrics, importance = vi_df, shap_importance = shap_df, model = model)
}


# =========================================================================
# ===== 6. PIPELINE ORCHESTRATION ENGINE =====
# =========================================================================

execute_pipeline <- function() {
  message("=== STARTING RANDOM FOREST MODELLING PIPELINE ===")
  
  # 1. Setup Time and Directory Metadata Contexts upfront
  timestamp     <- format(Sys.time(), "%Y-%m-%d_%Hh%Mm%Ss")
  
  # 2. Load Data
  matrices <- load_target_matrices(SELECTED_VARIATION)
  
  # 3. Conditionally Filter Features
  filtering_results <- apply_feature_selection(matrices, FEATURE_LIST_PATH)
  matrices <- filtering_results$matrices
  
  # 4. Formulate Output String Metadata Tagging
  job_name      <- paste0(BASE_JOB_NAME, "_", SELECTED_VARIATION, filtering_results$mode_tag)
  clean_job_dir <- sub("_+$", "", job_name)
  
  # 5. Multi-Seed Iterations Processing Lookups
  message("--> Iterating over experimental seeds...")
  results <- lapply(SEEDS_TO_RUN, function(s) {
    vivo_res  <- run_rf_experiment(matrices$vivo, s, "vivo", train_n = VIVO_TRAIN_N, test_n = VIVO_TEST_N)
    vitro_res <- run_rf_experiment(matrices$vitro, s, "vitro", train_n = VITRO_TRAIN_N, test_n = VITRO_TEST_N)
    
    # Optional .rds Model Archiving
    if (SAVE_MODELS) {
      model_dir <- file.path("saved_models", paste0(clean_job_dir, "_", timestamp))
      if(!dir.exists(model_dir)) dir.create(model_dir, recursive = TRUE)
      
      saveRDS(vivo_res$model,  file.path(model_dir, paste0(job_name, "vivo_model_seed_",  s, ".rds")))
      saveRDS(vitro_res$model, file.path(model_dir, paste0(job_name, "vitro_model_seed_", s, ".rds")))
    }
    
    list(
      metrics         = rbind(vivo_res$metrics, vitro_res$metrics),
      importance      = rbind(vivo_res$importance, vitro_res$importance),
      shap_importance = rbind(vivo_res$shap_importance, vitro_res$shap_importance)
    )
  })
  
  # 6. Collapse Sublists to Clean Matrices
  results_df    <- do.call(rbind, lapply(results, `[[`, "metrics"))
  importance_df <- do.call(rbind, lapply(results, `[[`, "importance"))
  shap_df       <- do.call(rbind, lapply(results, `[[`, "shap_importance"))
  
  # 7. Generate Analytical Calculation Summaries
  summary_stats <- results_df %>%
    group_by(dataset) %>%
    summarise(
      mean_accuracy = mean(accuracy), sd_accuracy = sd(accuracy),
      mean_kappa    = mean(kappa),    sd_kappa    = sd(kappa),
      mean_auc      = mean(auc),      sd_auc      = sd(auc)
    )
  
  importance_summary <- importance_df %>%
    group_by(dataset, feature) %>%
    summarise(mean_importance = mean(importance), sd_importance = sd(importance), .groups = "drop") %>%
    arrange(dataset, desc(mean_importance))
  
  shap_summary <- shap_df %>%
    group_by(dataset, feature) %>%
    summarise(mean_shap = mean(shap_importance), sd_shap = sd(shap_importance), .groups = "drop") %>%
    arrange(dataset, desc(mean_shap))
  
  # 8. Print Console Diagnostics
  print(summary_stats)
  print(head(shap_summary))
  
  # 9. Export Result Sheets Safely into Managed Folders
  job_dir <- file.path("model_results", paste0(clean_job_dir, "_", timestamp))
  if(!dir.exists(job_dir)) dir.create(job_dir, recursive = TRUE)
  
  write.csv(results_df,         file.path(job_dir, paste0(job_name, "rf_seed_results_", timestamp, ".csv")), row.names = FALSE)
  write.csv(summary_stats,       file.path(job_dir, paste0(job_name, "rf_summary_stats_", timestamp, ".csv")), row.names = FALSE)
  write.csv(importance_df,      file.path(job_dir, paste0(job_name, "rf_feature_importance_", timestamp, ".csv")), row.names = FALSE)
  write.csv(importance_summary, file.path(job_dir, paste0(job_name, "rf_feature_importance_summary_", timestamp, ".csv")), row.names = FALSE)
  write.csv(shap_df,            file.path(job_dir, paste0(job_name, "rf_shap_importance_", timestamp, ".csv")), row.names = FALSE)
  write.csv(shap_summary,       file.path(job_dir, paste0(job_name, "rf_shap_importance_summary_", timestamp, ".csv")), row.names = FALSE)
  
  message("--> Pipeline Execution Successful.")
  message("    Data metrics saved inside: ", job_dir)
  if (SAVE_MODELS) message("    Trained models saved inside: saved_models/", paste0(clean_job_dir, "_", timestamp))
}

# =========================================================================
# ===== 7. EXECUTION TRIGGER =====
# =========================================================================
execute_pipeline()
