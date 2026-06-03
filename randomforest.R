#install.packages("caret")
#install.packages("ranger")
library(caret)
library(ranger)

# load data
vivo_matrix <- read.csv("matrices/vivo_matrix.tsv", header = TRUE, sep = "\t", row.names = 1)
vitro_matrix <- read.csv("matrices/vitro_matrix.tsv", header = TRUE, sep = "\t", row.names = 1)


# transpose data.frames
vivo_matrix <- as.data.frame(t(vivo_matrix))
vitro_matrix <- as.data.frame(t(vitro_matrix))


# add 'PregnancyStatus'-column to data.frames, and automatically fill column for every row
# vivo
rn <- rownames(vivo_matrix)

vivo_matrix2 <- cbind(
  PregnancyStatus = ifelse(
    grepl("NP", rn),
    "NOT_PREGNANT",
    "PREGNANT"
  ),
  vivo_matrix
)

# vitro
rn <- rownames(vitro_matrix)

vitro_matrix2 <- cbind(
  PregnancyStatus = ifelse(
    grepl("NP", rn),
    "NOT_PREGNANT",
    "PREGNANT"
  ),
  vitro_matrix
)


# train-test split
set.seed(64) # setting seed

# vivo
vivo_matrix2$PregnancyStatus <- factor(vivo_matrix2$PregnancyStatus) # converts 'PregnancyStatus' target variable to a factor for classification

vivo_train_index <- createDataPartition(vivo_matrix2$PregnancyStatus, p = 0.8, list = FALSE) # 80/20 split
vivo_training_set <- vivo_matrix2[vivo_train_index, ]
vivo_test_set <- vivo_matrix2[-vivo_train_index, ]

# vitro
vitro_matrix2$PregnancyStatus <- factor(vitro_matrix2$PregnancyStatus) # converts 'PregnancyStatus' target variable to a factor for classification

vitro_train_index <- createDataPartition(vitro_matrix2$PregnancyStatus, p = 0.8, list = FALSE) # 80/20 split
vitro_training_set <- vitro_matrix2[vitro_train_index, ]
vitro_test_set <- vitro_matrix2[-vitro_train_index, ]


# training Random Forest model
# vivo, ~40-45 min.
# model_rf_vivo <- train(PregnancyStatus ~ ., data = vivo_training_set, method = "ranger", importance = "impurity")

# save the trained model as .rds
saveRDS(model_rf_vivo, "models/rf_model_vivo.rds")

# load trained model from RDS-file instead of training manually
model_rf_vivo <- readRDS(file = "models/rf_model_vivo.rds")



