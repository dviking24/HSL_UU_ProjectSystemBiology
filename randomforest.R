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
# setting seed
set.seed(64)

# vivo
# converts 'PregnancyStatus' target variable to a factor for classification
vivo_matrix2$PregnancyStatus <- factor(vivo_matrix2$PregnancyStatus)

vivo_train_index <- createDataPartition(vivo_matrix2$PregnancyStatus, p = 0.8, list = FALSE)
vivo_training_set <- vivo_matrix2[vivo_train_index, ]
vivo_test_set <- vivo_matrix2[-vivo_train_index, ]






