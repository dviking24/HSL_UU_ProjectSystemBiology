# ===== IMPORTS =====
#install.packages("caret")
#install.packages("ranger")
#install.packages("dplyr")
#install.packages("treeshap")
library(caret)
library(ranger)
library(dplyr)
library(treeshap)
library(pROC)



# ===== LOAD DATA =====
vivo_matrix <- read.csv("matrices/vivo_matrix.tsv", header = TRUE, sep = "\t", row.names = 1)
vitro_matrix <- read.csv("matrices/vitro_matrix.tsv", header = TRUE, sep = "\t", row.names = 1)

# Replace . with |
colnames(vivo_matrix)  <- gsub("\\.", "|", colnames(vivo_matrix))
colnames(vitro_matrix) <- gsub("\\.", "|", colnames(vitro_matrix))

# Function to keep only EMBRYO|ENDOMETRIUM columns
keep_embryo_to_endo <- function(mat) {
  cols <- colnames(mat)
  split_cols <- strsplit(cols, "\\|")
  
  left_part  <- sapply(split_cols, `[`, 1)
  right_part <- sapply(split_cols, `[`, 2)
  
  # embryo identifiers
  is_embryo_left <- grepl("^(Zo|SW)_", left_part)
  
  # endometrium identifiers
  is_endo_right <- grepl("^(NP|PR)_", right_part)
  
  mat[, is_embryo_left & is_endo_right, drop = FALSE]
}

vivo_matrix_emb2endo  <- keep_embryo_to_endo(vivo_matrix)
vitro_matrix_emb2endo <- keep_embryo_to_endo(vitro_matrix)

# transpose data.frames
vivo_matrix <- as.data.frame(t(vivo_matrix_emb2endo))
vitro_matrix <- as.data.frame(t(vitro_matrix_emb2endo))



# ===== ADD 'PregnancyStatus' =====
# ----- vivo -----
rn <- rownames(vivo_matrix)
vivo_matrix <- cbind(
  PregnancyStatus = ifelse(
    grepl("NP", rn),
    "NOT_PREGNANT",
    "PREGNANT"
  ),
  vivo_matrix
)

vivo_matrix$PregnancyStatus <- factor(vivo_matrix$PregnancyStatus)

# ----- vitro -----
rn <- rownames(vitro_matrix)
vitro_matrix <- cbind(
  PregnancyStatus = ifelse(
    grepl("NP", rn),
    "NOT_PREGNANT",
    "PREGNANT"
  ),
  vitro_matrix
)

vitro_matrix$PregnancyStatus <- factor(vitro_matrix$PregnancyStatus)



# ===== TRAIN-TEST SPLIT =====
extract_entities <- function(df) {
  rn <- rownames(df)
  parts <- strsplit(rn, "\\|")
  embryos <- unique(sapply(parts, `[`, 1))
  endos   <- unique(sapply(parts, `[`, 2))
  
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

sample_entities <- function(embryos, endos, train_n = 6, test_n = 3, seed = 64) {
  set.seed(seed)
  emb <- split_embryos(embryos)
  end <- split_endos(endos)
  
  train_emb_PR <- sample(emb$PR, train_n)
  train_emb_NP <- sample(emb$NP, train_n)
  
  remaining_emb_PR <- setdiff(emb$PR, train_emb_PR)
  remaining_emb_NP <- setdiff(emb$NP, train_emb_NP)
  
  test_emb_PR <- sample(remaining_emb_PR, test_n)
  test_emb_NP <- sample(remaining_emb_NP, test_n)
  
  train_end_PR <- sample(end$PR, train_n)
  train_end_NP <- sample(end$NP, train_n)
  
  remaining_end_PR <- setdiff(end$PR, train_end_PR)
  remaining_end_NP <- setdiff(end$NP, train_end_NP)
  
  test_end_PR <- sample(remaining_end_PR, test_n)
  test_end_NP <- sample(remaining_end_NP, test_n)
  
  list(
    train_embryos = c(train_emb_PR, train_emb_NP),
    test_embryos  = c(test_emb_PR, test_emb_NP),
    train_endos = c(train_end_PR, train_end_NP),
    test_endos  = c(test_end_PR, test_end_NP)
  )
}

build_split <- function(df, split_info) {
  rn <- rownames(df)
  parts <- strsplit(rn, "\\|")
  embryo <- sapply(parts, `[`, 1)
  endo   <- sapply(parts, `[`, 2)
  
  train_idx <- embryo %in% split_info$train_embryos & endo %in% split_info$train_endos
  test_idx <- embryo %in% split_info$test_embryos & endo %in% split_info$test_endos
  
  list(
    train = df[train_idx, , drop = FALSE],
    test  = df[test_idx, , drop = FALSE]
  )
}



# ===== TRAINING, TESTING AND EVALUATING RANDOM FOREST MODELS =====
run_rf_experiment <- function(df, seed, dataset_name, train_n, test_n) { 
  
  set.seed(seed)
  
  # ---- split ----
  entities <- extract_entities(df)
  
  split_info <- sample_entities(
    entities$embryos,
    entities$endos,
    train_n = train_n, # Pass it to the sampler
    test_n = test_n,   # Pass it to the sampler
    seed = seed
  )
  
  splits <- build_split(df, split_info)
  
  train_df <- splits$train
  test_df  <- splits$test
  
  # ---- PREPARE DATA FOR TREESHAP (The 0/1 Trick) ----
  x_train <- train_df[, names(train_df) != "PregnancyStatus", drop = FALSE]
  
  # Convert Target to Numeric: NOT_PREGNANT = 0, PREGNANT = 1
  y_train_num <- ifelse(train_df$PregnancyStatus == "PREGNANT", 1, 0)
  
  x_test  <- test_df[, names(test_df) != "PregnancyStatus", drop = FALSE]
  
  # Keep original factor for the Confusion Matrix
  y_test_factor <- test_df$PregnancyStatus 
  
  # ---- model ----
  model <- train(
    x = x_train,
    y = y_train_num,      # Train on numeric 0 and 1
    method = "ranger",
    importance = "impurity"
  )
  
  # ---- predictions ----
  # Because we trained on 0/1, the model predicts a numeric probability of PREGNANT 
  # (e.g., 0.85 means 85% confidence it is PREGNANT)
  preds_prob <- predict(model, x_test)
  
  # Threshold at 0.5 to convert the probability back to class labels
  preds_class <- ifelse(preds_prob > 0.5, "PREGNANT", "NOT_PREGNANT")
  preds_factor <- factor(preds_class, levels = levels(y_test_factor))
  
  cm <- confusionMatrix(preds_factor, y_test_factor)
  
  # ---- ROC and AUC ANALYSIS ----
  # Calculate the ROC curve by comparing the true labels to the predicted probabilities
  roc_obj <- roc(response = y_test_factor, predictor = preds_prob, quiet = TRUE)
  
  # Extract the Area Under the Curve (AUC)
  auc_value <- as.numeric(auc(roc_obj))
  
  # ---- base feature importance ----
  vi <- varImp(model)$importance
  vi_df <- data.frame(
    feature = rownames(vi),
    importance = vi[,1],
    dataset = dataset_name,
    seed = seed
  )
  
  # ---- SHAP ANALYSIS ----
  ranger_model <- model$finalModel
  
  # Unify and compute SHAP (This will now run smoothly as a numeric tree)
  unified <- treeshap::unify(ranger_model, x_train)
  shap_res <- treeshap::treeshap(unified, x_test, verbose = FALSE)
  
  # Aggregate Global SHAP Importance
  mean_abs_shap <- colMeans(abs(shap_res$shaps))
  
  shap_df <- data.frame(
    feature = names(mean_abs_shap),
    shap_importance = unname(mean_abs_shap),
    dataset = dataset_name,
    seed = seed
  )
  
  # ---- results ----
  metrics <- data.frame(
    dataset = dataset_name,
    seed = seed,
    accuracy = cm$overall["Accuracy"],
    kappa = cm$overall["Kappa"],
    auc = auc_value
  )
  
  list(
    metrics = metrics,
    importance = vi_df,
    shap_importance = shap_df
  )
}



# ===== RUN PIPELINE =====
seeds <- c(64,28,21,94,41,12,53,22,17,62)
random_seeds <- sample.int(1000, 10) # generate vector with 10 random seeds from range 1-1000

# train, run and evaluate models on given vector of seeds
job_name <- "jobname_"
results <- lapply(seeds, function(s) { # 'seeds' or 'random_seeds'
  
  # Run vivo with a 7/4 split
  vivo_res <- run_rf_experiment(vivo_matrix, s, "vivo", train_n = 7, test_n = 4)
  
  # Run vitro with a 5/4 split
  vitro_res <- run_rf_experiment(vitro_matrix, s, "vitro", train_n = 5, test_n = 4)
  
  list(
    metrics = rbind(vivo_res$metrics, vitro_res$metrics),
    importance = rbind(vivo_res$importance, vitro_res$importance),
    shap_importance = rbind(vivo_res$shap_importance, vitro_res$shap_importance)
  )
})

# combine results
results_df <- do.call(rbind, lapply(results, `[[`, "metrics"))
importance_df <- do.call(rbind, lapply(results, `[[`, "importance"))
shap_df <- do.call(rbind, lapply(results, `[[`, "shap_importance"))

# make summary & save results to files
summary_stats <- results_df %>%
  group_by(dataset) %>%
  summarise(
    mean_accuracy = mean(accuracy),
    sd_accuracy   = sd(accuracy),
    mean_kappa    = mean(kappa),
    sd_kappa      = sd(kappa),
    mean_auc    = mean(auc),
    sd_auc      = sd(auc)
  )

importance_summary <- importance_df %>%
  group_by(dataset, feature) %>%
  summarise(
    mean_importance = mean(importance),
    sd_importance   = sd(importance),
    .groups = "drop"
  ) %>%
  arrange(dataset, desc(mean_importance))

shap_summary <- shap_df %>%
  group_by(dataset, feature) %>%
  summarise(
    mean_shap = mean(shap_importance),
    sd_shap   = sd(shap_importance),
    .groups = "drop"
  ) %>%
  arrange(dataset, desc(mean_shap))

print(summary_stats)
print(head(shap_summary)) # Just printing the top features to keep console clean

timestamp <- format(Sys.time(), "%Y-%m-%d_%Hh%Mm%Ss")
write.csv(results_df, paste0("model_results/", job_name, "rf_seed_results_", timestamp, ".csv"), row.names = FALSE)
write.csv(summary_stats, paste0("model_results/", job_name, "rf_summary_stats_", timestamp, ".csv"), row.names = FALSE)
write.csv(importance_df, paste0("model_results/", job_name, "rf_feature_importance_", timestamp, ".csv"), row.names = FALSE)
write.csv(importance_summary, paste0("model_results/", job_name, "rf_feature_importance_summary_", timestamp, ".csv"), row.names = FALSE)
write.csv(shap_df, paste0("model_results/", job_name, "rf_shap_importance_", timestamp, ".csv"), row.names = FALSE)
write.csv(shap_summary, paste0("model_results/", job_name, "rf_shap_importance_summary_", timestamp, ".csv"), row.names = FALSE)
