#install.packages("caret")
library(caret)
library(class)
library(ggplot2)
#install.packages("pROC")
library(pROC)

voorbereiding_machinelearning <- function(naam_tsv) {
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
  
  return(list(
    X = X, 
    y = y
  ))
}


train_test_split <- function(X, y) {
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
  return(list(
    X_train = X_train,
    X_test = X_test,
    y_train = y_train,
    y_test = y_test
  ))
}


svm_machine_learning <- function(train_and_test) {
  X_train <- train_and_test$X_train
  X_test  <- train_and_test$X_test
  y_train <- train_and_test$y_train
  y_test  <- train_and_test$y_test
  
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
    method = "svmLinear2",
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
  
  return(list(
    train_conf = train_conf,
    test_conf = test_conf,
    probs = probs,
    y_test = y_test,
    pred = pred,
    X_test = X_test
  ))
}

accuracy <-  function(train_conf, test_conf) {
  
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
}

AUC <- function(probs, y_test) {
  class(probs)
  str(probs)
  head(probs)
  
  method = "svmLinear2"
  head(probs)
  
  roc_obj <- roc(
    response = y_test,
    predictor = probs$PREGNANT
  )
  
  auc(roc_obj)
  
  plot(
    roc_obj,
    print.auc = TRUE
  )
}

heatmap <- function(pred, y_test) {
  cm <- confusionMatrix(pred, y_test)
  
  cm_df <- as.data.frame(cm$table)
  
  ggplot(cm_df, aes(Prediction, Reference, fill = Freq)) +
    geom_tile() +
    geom_text(aes(label = Freq)) +
    scale_fill_gradient(low = "white", high = "red") +
    theme_minimal() +
    ggtitle("Confusion Matrix - SVM")
}


pca <- function(X_test, y_test, pred) {
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
}

function_training_testing <- function(naam_tsv) {
  X_and_y <- voorbereiding_machinelearning(naam_tsv)
  train_and_test <- train_test_split(
    X_and_y$X,
    X_and_y$y
  )
  lijst <- svm_machine_learning(train_and_test)
  train_conf <- lijst$train_conf
  test_conf <- lijst$test_conf
  probs <- lijst$probs
  y_test <- lijst$y_test
  pred <-  lijst$pred
  X_test <- lijst$X_test
  
  accuracy(train_conf, test_conf)
  
  #AUC(probs, y_test)
  heatmap(pred, y_test)
  pca(X_test, y_test, pred)
}

function_training_testing("matrices/vivo_matrix.tsv")
function_training_testing("matrices/vitro_matrix.tsv")




