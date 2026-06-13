# ===== IMPORTS =====
#install.packages("caret")
#install.packages("ranger")
library(caret)
library(ranger)



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


# ----- vivo train/test split -----
vivo_entities  <- extract_entities(vivo_matrix)

vivo_split_info <- sample_entities(
  vivo_entities$embryos,
  vivo_entities$endos,
  seed = 64
)

vivo_sets <- build_split(
  vivo_matrix,
  vivo_split_info
)

vivo_train <- vivo_sets$train
vivo_test  <- vivo_sets$test

# ----- vitro train/test split -----
vitro_entities <- extract_entities(vitro_matrix)

vitro_split_info <- sample_entities(
  vitro_entities$embryos,
  vitro_entities$endos,
  seed = 64
)

vitro_sets <- build_split(
  vitro_matrix,
  vitro_split_info
)

vitro_train <- vitro_sets$train
vitro_test  <- vitro_sets$test



# ===== TRAINING, TESTING AND EVALUATING RANDOM FOREST MODELS =====
# ----- vivo random forest -----
#model_rf_vivo <- train(PregnancyStatus ~ ., data = vivo_training_set, method = "ranger", importance = "impurity")
#saveRDS(model_rf_vivo, "models/rf_model_vivo.rds") # save the trained model as .rds
model_rf_vivo <- readRDS(file = "models/rf_model_vivo.rds") # load trained vivo rf-model from RDS-file instead of training manually

predict_rf_vivo <- predict(model_rf_vivo, vivo_test_set) # predicting the vivo test-set using Random Forest model

cm_rf_vivo <- confusionMatrix(predict_rf_vivo, vivo_test_set$PregnancyStatus) # evaluate performance using a confusion matrix
print(cm_rf_vivo)

print(varImp(model_rf_vivo)) # display feature importance

# ----- vitro random forest -----
#model_rf_vitro <- train(PregnancyStatus ~ ., data = vitro_training_set, method = "ranger", importance = "impurity")
#saveRDS(model_rf_vitro, "models/rf_model_vitro.rds")
model_rf_vitro <- readRDS(file = "models/rf_model_vitro.rds")

predict_rf_vitro <- predict(model_rf_vitro, vitro_test_set)

cm_rf_vitro <- confusionMatrix(predict_rf_vitro, vitro_test_set$PregnancyStatus)
print(cm_rf_vitro)

print(varImp(model_rf_vitro))


