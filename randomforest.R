# ===== IMPORTS =====
#install.packages("caret")
#install.packages("ranger")
library(caret)
library(ranger)
library(dplyr)



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
# add a 'PregnancyStatus'-column to the new data.frames, and fill this column for every row
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

# ----- vitro -----
vivo_matrix$PregnancyStatus <- factor(vivo_matrix$PregnancyStatus) # converts 'PregnancyStatus' target variable to a factor for classification

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
# ----- functions -----
extract_entities <- function(df) {
  
  rn <- rownames(df)
  
  parts <- strsplit(rn, "\\|")
  
  embryos <- unique(sapply(parts, `[`, 1))
  endos   <- unique(sapply(parts, `[`, 2))
  
  list(
    embryos = embryos,
    endos = endos
  )
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


sample_entities <- function(embryos,
                            endos,
                            train_n = 5,
                            test_n = 4,
                            seed = 64) {
  
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
  
  train_idx <- embryo %in% split_info$train_embryos &
    endo   %in% split_info$train_endos
  
  test_idx <- embryo %in% split_info$test_embryos &
    endo   %in% split_info$test_endos
  
  list(
    train = df[train_idx, , drop = FALSE],
    test  = df[test_idx, , drop = FALSE]
  )
}



# ===== TRAINING, TESTING AND EVALUATING RANDOM FOREST MODELS =====
run_rf_experiment <- function(df, seed, dataset_name) {
  
  set.seed(seed)
  
  # ---- split ----
  entities <- extract_entities(df)
  
  split_info <- sample_entities(
    entities$embryos,
    entities$endos,
    seed = seed
  )
  
  splits <- build_split(df, split_info)
  
  train_df <- splits$train
  test_df  <- splits$test
  
  # ---- model ----
  model <- train(
    PregnancyStatus ~ .,
    data = train_df,
    method = "ranger",
    importance = "impurity"
  )
  
  # ---- predictions ----
  preds <- predict(model, test_df)
  cm <- confusionMatrix(preds, test_df$PregnancyStatus)
  
  # ---- feature importance ----
  vi <- varImp(model)$importance
  vi_df <- data.frame(
    feature = rownames(vi),
    importance = vi[,1],
    dataset = dataset_name,
    seed = seed
  )
  
  # ---- results ----
  metrics <- data.frame(
    dataset = dataset_name,
    seed = seed,
    accuracy = cm$overall["Accuracy"],
    kappa = cm$overall["Kappa"]
  )
  
  list(
    metrics = metrics,
    importance = vi_df
  )
}


# ===== RUN PIPELINE =====
seeds <- c(64,28,21,94,41,12,53,22,17,62)
random_seeds <- sample.int(1000, 10) # generate vector with 10 random seeds from range 1-1000

# train, run and evaluate models on given vector of seeds
results <- lapply(seeds, function(s) {
  
  vivo_res <- run_rf_experiment(vivo_matrix, s, "vivo")
  vitro_res <- run_rf_experiment(vitro_matrix, s, "vitro")
  
  list(
    metrics = rbind(vivo_res$metrics, vitro_res$metrics),
    importance = rbind(vivo_res$importance, vitro_res$importance)
  )
})

# combine results
results_df <- do.call(rbind, lapply(results, `[[`, "metrics"))
importance_df <- do.call(rbind, lapply(results, `[[`, "importance"))

# make summary & save results to file
aggregate(accuracy ~ dataset, data = results_df, mean)
aggregate(kappa ~ dataset, data = results_df, mean)

summary_stats <- results_df %>%
  group_by(dataset) %>%
  summarise(
    mean_accuracy = mean(accuracy),
    sd_accuracy   = sd(accuracy),
    mean_kappa    = mean(kappa),
    sd_kappa      = sd(kappa)
  )

importance_summary <- importance_df %>%
  group_by(dataset, feature) %>%
  summarise(
    mean_importance = mean(importance),
    sd_importance   = sd(importance),
    .groups = "drop"
  ) %>%
  arrange(dataset, desc(mean_importance))

print(results_df)
print(summary_stats)
print(importance_df)
print(importance_summary)

timestamp <- format(Sys.time(), "%Y-%m-%d_%Hh%Mm%Ss")
write.csv(results_df,
          paste0("model_results/rf_seed_results_", timestamp, ".csv"),
          row.names = FALSE)
write.csv(summary_stats,
          paste0("model_results/rf_summary_stats_", timestamp, ".csv"),
          row.names = FALSE)
write.csv(importance_df,
          paste0("model_results/rf_feature_importance_", timestamp, ".csv"),
          row.names = FALSE)
write.csv(importance_summary,
          paste0("model_results/rf_feature_importance_summary_", timestamp, ".csv"),
          row.names = FALSE)
