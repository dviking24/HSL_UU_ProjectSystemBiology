#install.packages("caret")
library(caret)
library(class)
library(ggplot2)
#install.packages("pROC")
library(pROC)
#install.packages("e1071")
library(e1071)



naam_tsv <- "matrices/vivo_matrix.tsv"

matrix <- read.table(
  naam_tsv,
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE
)

matrix <- as.data.frame(t(matrix))

rn <- rownames(matrix)

matrix$PregnancyStatus <- factor(
  ifelse(
    grepl("NP", rn),
    "NOT_PREGNANT",
    "PREGNANT"
  )
)

X <- matrix[, !colnames(matrix) %in% "PregnancyStatus"]
y <- matrix$PregnancyStatus

set.seed(64)

train_index <- createDataPartition(
  y,
  p = 0.8,
  list = FALSE
)

X_train <- X[train_index, ]
X_test  <- X[-train_index, ]

y_train <- y[train_index]
y_test  <- y[-train_index]


# crossvalidatie
ctrl <- trainControl(
  method = "repeatedcv",
  number = 5,
  repeats = 3
)

# SVM
set.seed(64)

svm_model <- train(
  x = X_train,
  y = y_train,
  method = "svmLinear",
  trControl = ctrl,
  tuneLength = 10
)

print(svm_model)

#voorspellen
pred <- predict(
  svm_model,
  X_test
)

#resultaten
conf <- confusionMatrix(
  pred,
  y_test
)

train_pred <- predict(
  svm_model,
  X_train
)

train_conf <- confusionMatrix(
  train_pred,
  y_train
)

test_pred <- predict(
  svm_model,
  X_test
)

test_conf <- confusionMatrix(
  test_pred,
  y_test
)

probs <- predict(
  svm_model,
  X_test,
  type = "prob"
)

cat(
  "\nTrain accuracy:",
  round(train_conf$overall["Accuracy"], 4),
  "\n"
)

cat(
  "\nTest accuracy:",
  round(test_conf$overall["Accuracy"], 4),
  "\n"
)

cm <- confusionMatrix(pred, y_test)

cm_df <- as.data.frame(cm$table)

ggplot(cm_df, aes(Prediction, Reference, fill = Freq)) +
  geom_tile() +
  geom_text(aes(label = Freq)) +
  scale_fill_gradient(low = "white", high = "red") +
  theme_minimal() +
  ggtitle("Confusion Matrix - SVM")

pca <- prcomp(X_test, scale. = TRUE)

pca_df <- data.frame(
  PC1 = pca$x[,1],
  PC2 = pca$x[,2],
  label = y_test,
  pred = pred
)

ggplot(pca_df, aes(PC1, PC2, color = label)) +
  geom_point(size = 3) +
  theme_minimal() +
  ggtitle("PCA of LR interaction space")



top_features <- names(sort(apply(X_train, 2, var), decreasing = TRUE))[1:2]

X_small <- X_train[, top_features]

model_small <- svm(X_small, y_train, kernel = "linear")
View(model_small)



imp <- varImp(svm_model)
plot(imp, top = 20)




