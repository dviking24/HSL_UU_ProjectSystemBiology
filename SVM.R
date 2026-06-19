##                      PACKAGES                            ##

#install.packages("caret")
library(caret)
library(class)
library(ggplot2)
#install.packages("pROC")
library(pROC)
#install.packages("e1071")
library(e1071)
#install.packages("kernelshap")
library(kernelshap)


SELECTED_VARIATION <- "both" 


FEATURE_LIST_PATH <- ""

BASE_JOB_NAME      <- ""

SAVE_MODELS        <- TRUE

SEEDS_TO_RUN       <- c(64, 28, 21, 94, 41, 12, 53, 22, 17, 62)

VIVO_TRAIN_N       <- 7
VIVO_TEST_N        <- 4

VITRO_TRAIN_N      <- 5
VITRO_TEST_N       <- 4


##                     LOADING DATA                         ##

pregnancy_status <- function(matrix) {
  matrix <- as.data.frame(t(matrix))
  rn <- rownames(matrix)
  matrix$PregnancyStatus <- factor(
    ifelse(
      grepl("PR", rn),
      "PREGNANT",
      "NOT_PREGNANT"
    ),
    levels = c("NOT_PREGNANT", "PREGNANT")
  )
  return(matrix)
}

# Inladen matrix
load_target_matrices <- function(vivo_file_name, vitro_file_name) {
  vivo_path  <- file.path("matrices", paste0(vivo_file_name))
  vitro_path <- file.path("matrices", paste0(vitro_file_name))
  
  vivo  <- read.table(vivo_path, header = TRUE, sep = "\t", row.names = 1)
  vitro <- read.table(vitro_path, header = TRUE, sep = "\t", row.names = 1)
  
  vivo <- pregnancy_status(vivo)
  vitro <- pregnancy_status(vitro)
 
  return(
    list(vivo = vivo, 
         vitro = vitro)
  )
}



##               TRAIN TEST SET                    ##


extract_entities <- function(martix) {
  rn <- rownames(martix)
  parts <- strsplit(rn, "\\.")
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
build_split <- function(matrix, split_info) {
  print(head(rownames(vivo_matrix)))
  rn <- rownames(matrix)
  parts <- strsplit(rn, "\\.")
  print(head(parts))
  
  train_idx <- sapply(parts, function(p) {
    any(p %in% split_info$train_embryos) & any(p %in% split_info$train_endos)
  })
  
  test_idx <- sapply(parts, function(p) {
    any(p %in% split_info$test_embryos) & any(p %in% split_info$test_endos)
  })
  print(head(test_idx))
  
  list(
    train = matrix[train_idx, , drop = FALSE],
    test  = matrix[test_idx, , drop = FALSE]
  )
}


##              SPLIT UITVOEREN                 ##

# Features en labels
features_and_labels <- function(split_data) {
  train_matrix <- split_data$train
  test_matrix <- split_data$test
  
  
  X_train <- train_matrix[, names(train_matrix) != "PregnancyStatus", drop = FALSE]
  y_train <- train_matrix$PregnancyStatus
  
  X_test <- test_matrix[, names(test_matrix) != "PregnancyStatus", drop = FALSE]
  y_test <- test_matrix$PregnancyStatus
  
  # Controle
  cat("Train rows:", nrow(train_matrix), "\n")
  cat("Test rows :", nrow(test_matrix), "\n")
  
  table(y_train)
  table(y_test)
  
  return(list(
    X_train = X_train,
    X_test = X_test,
    y_train = y_train,
    y_test = y_test
  ))
}


preprocessing <- function(matrix, TRAIN_N, TEST_N, seed) {
  set.seed(64)
  entities <- extract_entities(matrix)
  
  split_info <- sample_entities(
    embryos = entities$embryos,
    endos   = entities$endos,
    train_n = TRAIN_N,
    test_n  = TEST_N,
    seed    = seed
  )
  
  print(head(split_info))
  
  split_data <- build_split(matrix, split_info)
  test_train_sets <- features_and_labels(split_data)
  
  return(test_train_sets)
}


##                          MODEL TRAINING                         ##

SVM_model_running <-  function(method, X_train, y_train) {
  # crossvalidatie
  ctrl <- trainControl(
    method = "repeatedcv",
    number = 5,
    repeats = 3,
    classProbs = TRUE,
    summaryFunction = twoClassSummary
  )
  
  grid <- expand.grid(
    C = c(
      0.001,
      0.01,
      0.1,
      1,
      10,
      100
    )
  )
  
  # SVM
  set.seed(64)
  
  svm_model <- train(
    x = X_train,
    y = y_train,
    method = method,
    metric = "ROC",
    preProcess = c("center", "scale"),
    trControl = ctrl,
    tuneLength = grid
  )
  
  return(svm_model)
}

##                              PREDICTION VARIABLES                        ##

# Train and Test predictions of SVM model

evaluate_predictions <- function(model, X_train, y_train, X_test, y_test) {
  
  train_pred <- predict(model, X_train)
  test_pred  <- predict(model, X_test)
  
  train_conf <- confusionMatrix(train_pred, y_train)
  test_conf  <- confusionMatrix(test_pred, y_test)
  
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
  
  return(test_pred)
}


# ROC AUC curve
plot_roc_curve <- function(model, X_test, y_test) {
  
  probs <- predict(
    model,
    X_test,
    type = "prob"
  )
  
  roc_obj <- roc(
    response = y_test,
    predictor = probs$PREGNANT,
    levels = c("NOT_PREGNANT", "PREGNANT")
  )
  
  plot(
    roc_obj,
    print.auc = TRUE,
    main = "ROC curve"
  )
  
  return(roc_obj)
}


# Confusion Matrix
plot_confusion_matrix <- function(model, X_test, y_test) {
  
  pred <- predict(model, X_test)
  
  cm <- confusionMatrix(pred, y_test)
  
  cm_df <- as.data.frame(cm$table)
  
  ggplot(cm_df,
         aes(Prediction, Reference, fill = Freq)) +
    geom_tile() +
    geom_text(aes(label = Freq)) +
    scale_fill_gradient(low = "white", high = "red") +
    theme_minimal() +
    ggtitle("Confusion Matrix")
}


# PCA of the testset data
plot_pca <- function(X_test, y_test, pred) {
  
  pca <- prcomp(X_test, scale. = TRUE)
  
  pca_df <- data.frame(
    PC1 = pca$x[,1],
    PC2 = pca$x[,2],
    label = y_test,
    pred = pred
  )
  
  ggplot(
    pca_df,
    aes(PC1, PC2, color = label)
  ) +
    geom_point(size = 3) +
    theme_minimal() +
    ggtitle("PCA")
}


# Feature Importance 
plot_feature_importance <- function(model) {
  
  imp <- varImp(model)
  
  plot(
    imp,
    top = 20
  )
}

# Permutation Importance
# Bereken permutation importance met log-loss
permutation_importance <- function(model, X, y) {
  
  y_num <- ifelse(
    y == "PREGNANT",
    1,
    0
  )
  
  baseline_prob <- predict(
    model,
    X,
    type = "prob"
  )[,"PREGNANT"]
  
  baseline_loss <- -mean(
    y_num * log(baseline_prob + 1e-15) +
      (1 - y_num) * log(1 - baseline_prob + 1e-15)
  )
  
  imp <- numeric(ncol(X))
  
  for (i in seq_len(ncol(X))) {
    
    X_perm <- X
    
    X_perm[, i] <- sample(
      X_perm[, i]
    )
    
    perm_prob <- predict(
      model,
      X_perm,
      type = "prob"
    )[,"PREGNANT"]
    
    perm_loss <- -mean(
      y_num * log(perm_prob + 1e-15) +
        (1 - y_num) * log(1 - perm_prob + 1e-15)
    )
    
    imp[i] <- perm_loss - baseline_loss
  }
  
  data.frame(
    feature = colnames(X),
    importance = imp
  )[
    order(imp, decreasing = TRUE),
  ]
}

# Visualiseer permutation importance
plot_permutation_importance <- function(perm_df) {
  
  ggplot(
    perm_df,
    aes(
      x = reorder(feature, importance),
      y = importance
    )
  ) +
    geom_col(fill = "steelblue") +
    coord_flip() +
    theme_minimal() +
    labs(
      title = "Permutation Importance",
      x = "Feature",
      y = "Importance"
    )
}

# Visualiseer decision boundary met 2 features
plot_decision_boundary <- function(X_train, y_train) {
  
  top_features <- names(
    sort(
      apply(X_train, 2, var),
      decreasing = TRUE
    )
  )[1:2]
  
  X_small <- X_train[, top_features]
  
  model_small <- svm(
    X_small,
    y_train,
    kernel = "linear"
  )
  
  plot_df <- data.frame(
    X1 = X_small[,1],
    X2 = X_small[,2],
    label = y_train
  )
  
  grid <- expand.grid(
    X1 = seq(min(X_small[,1]),
             max(X_small[,1]),
             length.out = 200),
    X2 = seq(min(X_small[,2]),
             max(X_small[,2]),
             length.out = 200)
  )
  
  grid$pred <- predict(
    model_small,
    grid
  )
  
  ggplot() +
    geom_tile(
      data = grid,
      aes(X1, X2, fill = pred),
      alpha = 0.3
    ) +
    geom_point(
      data = plot_df,
      aes(X1, X2, color = label),
      size = 2
    ) +
    theme_minimal() +
    ggtitle("SVM Decision Boundary")
}


prediction_figures <- function(model, X_train, X_test, y_train, y_test) {
  test_pred <- evaluate_predictions(model, X_train, y_train, X_test, y_test)
  
  # roc_obj <- plot_roc_curve(model, X_test, y_test)
  # 
  # plot_confusion_matrix(model, X_test, y_test)
  # 
  # plot_pca(X_test, y_test, test_pred)
  # 
  # plot_feature_importance(model)
  # 
  # perm_df <- permutation_importance(model, X_test, y_test)
  # plot_permutation_importance(perm_df)
  # 
  # plot_decision_boundary(X_train, y_train)
}



##                        PERMUTATIONS IMPORTANCE                                ##


permutation_importance <- function(model, X, y) {
  
  baseline_pred <- predict(model, X)
  
  baseline_acc <- mean(baseline_pred == y)
  
  imp <- numeric(ncol(X))
  
  for (i in seq_len(ncol(X))) {
    
    X_perm <- X
    X_perm[, i] <- sample(X_perm[, i])
    
    perm_pred <- predict(model, X_perm)
    
    perm_acc <- mean(perm_pred == y)
    
    imp[i] <- baseline_acc - perm_acc
  }
  
  data.frame(
    feature = colnames(X),
    importance = imp
  )[order(imp, decreasing = TRUE), ]
}


permutation_importance_auc <- function(model, X, y) {
  
  base_prob <- predict(model, X, type = "prob")[,"PREGNANT"]
  
  base_auc <- auc(roc(y, base_prob))
  
  imp <- numeric(ncol(X))
  
  for (i in seq_len(ncol(X))) {
    
    Xp <- X
    Xp[, i] <- sample(Xp[, i])
    
    p <- predict(model, Xp, type = "prob")[,"PREGNANT"]
    
    auc_perm <- auc(roc(y, p))
    
    imp[i] <- base_auc - auc_perm
  }
  
  data.frame(
    feature = colnames(X),
    importance = imp
  )[order(imp, decreasing = TRUE), ]
}


##                                 SHAP                                          ##

# Kernal shap

shap_kernel_function <- function(model, X_train) {
  
  shap_values <- kernelshap::kernelshap(
    object = model,
    X = X_train[1:50, ],
    bg_X = X_train,
    pred_fun = function(object, newdata) {
      predict(object, newdata, type = "prob")[, "PREGNANT"]
    }
  )
  
  return(shap_values)
}

# SHAP results
shap_results <- function(shaps, feature_names) {
  
  shap_importance <- colMeans(abs(shaps))
  
  shap_df <- data.frame(
    feature = feature_names,
    importance = shap_importance
  )
  
  shap_df <- shap_df[order(shap_df$importance, decreasing = TRUE), ]
  
  head(shap_df, 10)
  ggplot(head(shap_df, 10), aes(x = reorder(feature, importance), y = importance)) +
    geom_col(fill = "steelblue") +
    coord_flip() +
    theme_minimal() +
    labs(title = "Top 10 SHAP features", x = "Gene", y = "Mean |SHAP|")
}

##                                PIPELINE                                       ##

vivo_file_name <- "vivo_lig-Endo_rec-Enbryo_matrix.tsv"
vitro_file_name <- "vitro_lig-Endo_rec-Enbryo_matrix.tsv"

matrices <- load_target_matrices(vivo_file_name, vitro_file_name)

vivo_matrix  <- matrices$vivo
#vitro_matrix <- matrices$vitro

test_train_sets <- preprocessing(vivo_matrix, VIVO_TRAIN_N, VIVO_TEST_N, 64)
#test_train_sets <- preprocessing(vitro_matrix)

X_train <- test_train_sets$X_train
X_test <- test_train_sets$X_test
y_train <- test_train_sets$y_train
y_test <- test_train_sets$y_test
  
svm_model_lin <- SVM_model_running("svmLinear", X_train, y_train)
# Radial has worse predictions then linear
# svm_model_rad <- SVM_model_running("svmRadial", X_train, y_train)
prediction_figures(svm_model_lin, X_train, X_test, y_train, y_test)

# permutations
perm_df <- permutation_importance(svm_model_lin, X_train, y_train)
top_features <- perm_df$feature[1:10]
X_train_small <- X_train[, top_features, drop = FALSE]
X_test_small  <- X_test[, top_features, drop = FALSE]
svm_permutation <- SVM_model_running("svmLinear", X_train_small, y_train)
prediction_figures(svm_permutation, X_train_small, X_test_small, y_train, y_test)

# AUC (niet goed)
# perm_df <- permutation_importance_auc(svm_model_lin, X_train, y_train)
# top_features <- perm_df$feature[1:10]
# X_train_small <- X_train[, top_features, drop = FALSE]
# X_test_small  <- X_test[, top_features, drop = FALSE]
# svm_permutation <- SVM_model_running("svmLinear", X_train_small, y_train)
# prediction_figures(svm_permutation, X_train_small, X_test_small, y_train, y_test)


#VarImp
imp <- varImp(svm_model_lin, scale = TRUE)
imp_df <- imp$importance
imp_df$Overall <- rowMeans(abs(imp_df))
imp_df$feature <- rownames(imp_df)
imp_df <- imp_df[order(imp_df$Overall, decreasing = TRUE), ]
top_features <- head(imp_df$feature, 30)
X_train_small <- X_train[, top_features, drop = FALSE]
X_test_small  <- X_test[, top_features, drop = FALSE]
svm_Imp <- SVM_model_running("svmLinear",X_train_small,y_train)
prediction_figures(svm_Imp, X_train_small, X_test_small, y_train, y_test)

shap_values <- shap_kernel_function(svm_Imp, X_train_small)
shap_results(shap_values$S, colnames(X_train_small))


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

head(probs)

levels(y_test)

roc_obj <- roc(
  response = y_test,
  predictor = probs$PREGNANT,
  levels = c("NOT_PREGNANT", "PREGNANT")
)

plot(
  roc_obj,
  print.auc = TRUE,
  main = "ROC curve - Pregnancy prediction"
)

roc_df <- data.frame(
  FPR = 1 - roc_obj$specificities,
  TPR = roc_obj$sensitivities
)

library(ggplot2)

ggplot(roc_df, aes(FPR, TPR)) +
  geom_line(linewidth = 1) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed"
  ) +
  labs(
    title = paste(
      "ROC curve (AUC =",
      round(as.numeric(auc(roc_obj)), 3),
      ")"
    ),
    x = "False Positive Rate",
    y = "True Positive Rate"
  ) +
  theme_minimal()

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

plot_df <- data.frame(
  X1 = X_small[,1],
  X2 = X_small[,2],
  label = y_train
)

x1_range <- seq(min(plot_df$X1), max(plot_df$X1), length.out = 200)
x2_range <- seq(min(plot_df$X2), max(plot_df$X2), length.out = 200)

grid <- expand.grid(
  X1 = x1_range,
  X2 = x2_range
)

grid_pred <- predict(model_small, grid)
grid$pred <- grid_pred

sv <- model_small$index

plot_df$sv <- FALSE
plot_df$sv[sv] <- TRUE

ggplot() +
  
  # decision regions
  geom_tile(
    data = grid,
    aes(X1, X2, fill = pred),
    alpha = 0.3
  ) +
  
  # data points
  geom_point(
    data = plot_df,
    aes(X1, X2, color = label),
    size = 2
  ) +
  
  scale_fill_manual(values = c("NOT_PREGNANT" = "lightblue",
                               "PREGNANT" = "pink")) +
  
  scale_color_manual(values = c("NOT_PREGNANT" = "blue",
                                "PREGNANT" = "red")) +
  
  theme_minimal() +
  ggtitle("SVM decision boundary (top 2 LR features)")

