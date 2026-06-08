library(class)
library(caret)
library(ggplot2)

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

features <- vivo_matrix2[, -1] # select all but first column ('PregnancyStatus') as features
labels <- factor(vivo_matrix2$PregnancyStatus)

# train-test split
set.seed(64) # setting seed

# vivo
vivo_train_index <- createDataPartition(vivo_matrix2$PregnancyStatus, p = 0.8, list = FALSE) # 80/20 split
vivo_train_features <- features[vivo_train_index, ]
vivo_test_features <- features[-vivo_train_index, ]
vivo_train_labels <- labels[vivo_train_index]
vivo_test_labels <- labels[-vivo_train_index]

# vitro
vitro_train_index <- createDataPartition(vitro_matrix2$PregnancyStatus, p = 0.8, list = FALSE) # 80/20 split
vitro_train_features <- features[vitro_train_index, ]
vitro_test_features <- features[-vitro_train_index, ]
vitro_train_labels <- labels[vitro_train_index]
vitro_test_labels <- labels[-vitro_train_index]


# training, predicting and evaluating KNN-Classifier models:
# ----- vivo -----
knn_pred_vivo <- knn(vivo_train_features, vivo_test_features, vivo_train_labels, k = 5)

# Confusion Matrix
print(confusionMatrix(knn_pred_vivo, vivo_test_labels))

# Plot
test_data_vivo <- cbind(vivo_test_features, Predicted = knn_pred_vivo)
ggplot(
  test_data_vivo,
  aes(
    x = `ENSBTAG00000017664-ENSBTAG00000019712`,
    y = `ENSBTAG00000001141-ENSBTAG00000013745`,
    color = Predicted
  )
) +
  geom_point(size = 3) +
  labs(title = "K-NN Predictions on Vivo Test Data") +
  theme_minimal()

# ----- vitro -----
knn_pred_vitro <- knn(vitro_train_features, vitro_test_features, vitro_train_labels, k = 5)

# Confusion Matrix
print(confusionMatrix(knn_pred_vitro, vitro_test_labels))

# plot
test_data_vitro <- cbind(vitro_test_features, Predicted = knn_pred_vitro)
ggplot(
  test_data_vitro,
  aes(
    x = `ENSBTAG00000017664-ENSBTAG00000019712`,
    y = `ENSBTAG00000001141-ENSBTAG00000013745`,
    color = Predicted
  )
) +
  geom_point(size = 3) +
  labs(title = "K-NN Predictions on Vitro Test Data") +
  theme_minimal()
