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

# add a 'PregnancyStatus'-column to the new data.frames, and fill this column for every row
# ----- vivo -----
rn <- rownames(vivo_matrix)
vivo_matrix2 <- cbind(
  PregnancyStatus = ifelse(
    grepl("NP", rn),
    "NOT_PREGNANT",
    "PREGNANT"
  ),
  vivo_matrix
)

# ----- vitro -----
rn <- rownames(vitro_matrix)
vitro_matrix2 <- cbind(
  PregnancyStatus = ifelse(
    grepl("NP", rn),
    "NOT_PREGNANT",
    "PREGNANT"
  ),
  vitro_matrix
)



# ===== TRAIN-TEST SPLIT =====
set.seed(64) # setting seed

# ----- vivo -----
vivo_matrix2$PregnancyStatus <- factor(vivo_matrix2$PregnancyStatus) # converts 'PregnancyStatus' target variable to a factor for classification
vivo_train_index <- createDataPartition(vivo_matrix2$PregnancyStatus, p = 0.7, list = FALSE) # 70/30 split
vivo_training_set <- vivo_matrix2[vivo_train_index, ]
vivo_test_set <- vivo_matrix2[-vivo_train_index, ]

# ----- vitro -----
vitro_matrix2$PregnancyStatus <- factor(vitro_matrix2$PregnancyStatus)
vitro_train_index <- createDataPartition(vitro_matrix2$PregnancyStatus, p = 0.7, list = FALSE)
vitro_training_set <- vitro_matrix2[vitro_train_index, ]
vitro_test_set <- vitro_matrix2[-vitro_train_index, ]



# ===== TRAINING, TESTING AND EVALUATING RANDOM FOREST MODELS =====
# ----- vivo -----
#model_rf_vivo <- train(PregnancyStatus ~ ., data = vivo_training_set, method = "ranger", importance = "impurity")
#saveRDS(model_rf_vivo, "models/rf_model_vivo.rds") # save the trained model as .rds
model_rf_vivo <- readRDS(file = "models/rf_model_vivo.rds") # load trained vivo rf-model from RDS-file instead of training manually

predict_rf_vivo <- predict(model_rf_vivo, vivo_test_set) # predicting the vivo test-set using Random Forest model

cm_rf_vivo <- confusionMatrix(predict_rf_vivo, vivo_test_set$PregnancyStatus) # evaluate performance using a confusion matrix
print(cm_rf_vivo)

print(varImp(model_rf_vivo)) # display feature importance

# ----- vitro -----
#model_rf_vitro <- train(PregnancyStatus ~ ., data = vitro_training_set, method = "ranger", importance = "impurity")
#saveRDS(model_rf_vitro, "models/rf_model_vitro.rds")
model_rf_vitro <- readRDS(file = "models/rf_model_vitro.rds")

predict_rf_vitro <- predict(model_rf_vitro, vitro_test_set)

cm_rf_vitro <- confusionMatrix(predict_rf_vitro, vitro_test_set$PregnancyStatus)
print(cm_rf_vitro)

print(varImp(model_rf_vitro))


