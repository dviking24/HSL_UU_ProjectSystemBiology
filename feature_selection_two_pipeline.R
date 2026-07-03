#------------------------------------------------------------------------------#
#                                PACKAGES                                      #
#------------------------------------------------------------------------------#

# Data handling & manipulation
#install.packages(c("data.table", "tibble", "glue"))
library(data.table)
library(tibble)
library(glue)
library(dplyr)
library(purrr)

# Machine learning & modelling
#install.packages(c("caret", "e1071", "class"))
library(caret)
library(e1071)
library(class)

# Evaluations & metrics
#install.packages("pROC")
library(pROC)

# Explainability (SHAP)
#install.packages("kernelshap")
library(kernelshap) 

# Overig
#install.packages("TDM")
library(TDM)
library(ggplot2)
library(png)
library(readr)
library(pheatmap)


#------------------------------------------------------------------------------#
#                              GLOBAL VARIABLES                                #
#------------------------------------------------------------------------------#

DATA_ENDOM        <- "Datasets/DataEndom.txt"
DATA_VIVO         <- "Datasets/Data_BlastoIVV.txt"
DATA_VITRO        <- "Datasets/Data_BlastoIVT.txt"
DATA_LR           <- "Datasets/LRdb_bovine_ENSEMBL.txt"

META_ENDOM        <- "Datasets/SampleInfo_Endom.txt"
META_VIVO         <- "Datasets/SampleInfo_BlastoIVV.txt"
META_VITRO        <- "Datasets/SampleInfo_BlastoIVT.txt"

VIVO_LIG_ENDO_REC_EMBRYO <- "matrices/vivo_lig-Endo_rec-Embryo_matrix.tsv"
VITRO_LIG_ENDO_REC_EMBRYO <- "matrices/vitro_lig-Endo_rec-Embryo_matrix.tsv"

VIVO_LIG_EMBRYO_REC_ENDO <- "matrices/vivo_lig-Embryo_rec-Endo_matrix.tsv"
VITRO_LIG_EMBRYO_REC_ENDO <- "matrices/vitro_lig-Embryo_rec-Endo_matrix.tsv"

WHICH             <- "VIVO_VITRO_lig-endo"

SEEDS_TO_RUN      <- c(64, 28, 21, 94, 41, 12, 53, 22, 17, 62)
# SEEDS_TO_RUN <- c(53)

COMPARE_VAR_IMP   <- c(10, 15, 20, 25, 30, 40, 60)
SHAP_AMT_FEATURE  <- 10

VIVO_TRAIN_A      <- 7
VIVO_TEST_A       <- 4

VITRO_TRAIN_A     <- 5
VITRO_TEST_A      <- 4







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

##-------------------------------------------------------------------------------##
##                   PREVENTING DATA LEAKAGE (PREPROCESSING)                     ##
##-------------------------------------------------------------------------------##


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

# Draw Uniform Resamples Without Category Clashes
sample_entities <- function(embryos, endos, train_n, test_n, seed) {
  cat("\nSeperating entitities\n")
  set.seed(seed)
  emb <- split_embryos(embryos)
  end <- split_endos(endos)
  
  cat("Emb PR:", length(emb$PR), "\n")
  cat("Emb NP:", length(emb$NP), "\n")
  cat("End PR:", length(end$PR), "\n")
  cat("End NP:", length(end$NP), "\n")
  
  cat("train_n =", train_n, "\n")
  cat("test_n =", test_n, "\n")
  
  
  train_emb_PR <- sample(emb$PR, train_n)
  train_emb_NP <- sample(emb$NP, train_n)
  
  cat("Remaining emb PR:",
      length(setdiff(emb$PR, train_emb_PR)), "\n")
  
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

# Assign Sample Pairs to Partitions Based on Complete Biological Matches
build_split <- function(matrix, split_info) {
  cat("\nBuilding new entities\n")
  rn <- rownames(matrix)
  parts <- strsplit(rn, "\\.")
  
  train_idx <- sapply(parts, function(p) {
    any(p %in% split_info$train_embryos) & any(p %in% split_info$train_endos)
  })
  
  test_idx <- sapply(parts, function(p) {
    any(p %in% split_info$test_embryos) & any(p %in% split_info$test_endos)
  })
  
  list(
    train = matrix[train_idx, , drop = FALSE],
    test  = matrix[test_idx, , drop = FALSE]
  )
}


##-------------------------------------------------------------------------------##
##                          SPLITTING TRAIN AND TEST SET                         ##
##-------------------------------------------------------------------------------##

# Features en labels
features_and_labels <- function(split_data) {
  cat("\nSplitting train and test data\n")
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


preprocessing <- function(matrix, train_n, test_n, seed) {
  cat("Starting preprocessing\n")
  set.seed(seed)
  
  entities <- extract_entities(matrix)
  
  split_info <- sample_entities(
    embryos = entities$embryos,
    endos   = entities$endos,
    train_n = train_n,
    test_n  = test_n,
    seed    = seed
  )
  
  split_data <- build_split(matrix, split_info)
  test_train_sets <- features_and_labels(split_data)
  
  return(test_train_sets)
}


##-------------------------------------------------------------------------------##
##                              MODEL TRAINING                                   ##
##-------------------------------------------------------------------------------##

SVM_model_running <-  function(method, X_train, y_train, seed) {
  cat("\nRunning SVM model\n")
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
  set.seed(seed)
  
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

##-------------------------------------------------------------------------------##
##                      PREDICTION VARIABLES AND VISUALISATION                   ##
##-------------------------------------------------------------------------------##

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
  
  return(list(
    train_pred = train_pred,
    test_pred = test_pred,
    train_conf = train_conf,
    test_conf = test_conf
  ))
}


# ROC AUC curve
plot_roc_curve <- function(model, X_test, y_test, seed_dir) {
  
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
  
  png(
    filename = file.path(seed_dir, "roc_curve.png"),
    width = 1200,
    height = 1000,
    res = 150
  )
  
  p <- plot(
    roc_obj,
    print.auc = TRUE,
    main = "ROC curve"
  )
  
  print(p)
  
  dev.off()
  
  return(list(
    roc = roc_obj,
    probabilities = probs
  ))
}


# Confusion Matrix
plot_confusion_matrix <- function(model, X_test, y_test, test_pred, seed_dir) {
  cm <- confusionMatrix(test_pred, y_test)
  cm_df <- as.data.frame(cm$table)
  
  png(
    filename = file.path(seed_dir, "confusion_matrix.png"),
    width = 1200,
    height = 1000,
    res = 150
  )
  
  p <- ggplot(cm_df,
              aes(Prediction, Reference, fill = Freq)) +
    geom_tile() +
    geom_text(aes(label = Freq)) +
    scale_fill_gradient(low = "white", high = "red") +
    theme_minimal() +
    ggtitle("Confusion Matrix")
  
  print(p)
  
  dev.off()
  
  return(list(
    cm = cm,
    cm_df = cm_df
  ))
}


# PCA of the testset data
plot_pca <- function(X_test, y_test, test_pred, seed_dir) {
  pca <- prcomp(X_test, scale. = TRUE)
  
  pca_df <- data.frame(
    PC1 = pca$x[,1],
    PC2 = pca$x[,2],
    label = y_test,
    pred = test_pred
  )
  
  png(
    filename = file.path(seed_dir, "pca.png"),
    width = 1200,
    height = 1000,
    res = 150
  )
  
  p <- ggplot(
    pca_df,
    aes(PC1, PC2, color = label)
  ) +
    geom_point(size = 3) +
    theme_minimal() +
    ggtitle("PCA")
  
  print(p)
  
  dev.off()
  
  return(list(
    pca = pca,
    pca_df = pca_df
  ))
}


# Feature Importance 
plot_feature_importance <- function(model, seed_dir) {
  imp <- varImp(model)
  
  png(
    filename = file.path(seed_dir, "varImportance.png"),
    width = 1200,
    height = 1000,
    res = 150
  )
  
  p <- plot(
    imp,
    top = 25
  )
  
  print(p)
  
  dev.off()
  
  return(imp)
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
  
  perm_df <- data.frame(
    feature = colnames(X),
    importance = imp
  )[
    order(imp, decreasing = TRUE),
  ]
  
  return(perm_df)
}

# Visualiseer permutation importance
plot_permutation_importance <- function(perm_df, seed_dir) {
  png(
    filename = file.path(seed_dir, "permutation_importance.png"),
    width = 1200,
    height = 1000,
    res = 150
  )
  
  p <- ggplot(
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
  
  print(p)
  
  dev.off()
}

# Visualiseer decision boundary met 2 features
plot_decision_boundary <- function(X_train, y_train, seed_dir) {
  
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
  
  png(
    filename = file.path(seed_dir, "decision_boundry.png"),
    width = 1200,
    height = 1000,
    res = 150
  )
  
  p <- ggplot() +
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
  
  print(p)
  
  dev.off()
  
  return(list(
    model = model_small,
    grid = grid,
    plot_df = plot_df,
    features = top_features
  ))
}


prediction_figures <- function(model, X_train, X_test, y_train, y_test, seed_dir) {
  pred_list <- evaluate_predictions(model, X_train, y_train, X_test, y_test)
  train_pred <- pred_list$train_pred
  test_pred <- pred_list$test_pred
  train_conf <- pred_list$train_conf
  test_conf <- pred_list$test_conf
  
  roc_prob <- plot_roc_curve(model, X_test, y_test, seed_dir)
  
  conf <- plot_confusion_matrix(model, X_test, y_test, test_pred, seed_dir)
  
  pca_data <- plot_pca(X_test, y_test, test_pred, seed_dir)
  
  imp <- plot_feature_importance(model, seed_dir)
  
  perm_df <- permutation_importance(model, X_test, y_test)
  plot_permutation_importance(perm_df, seed_dir)
  
  db_list <- plot_decision_boundary(X_train, y_train, seed_dir)
  # model_db <- db_list$model
  # grid_db <- db_list$grid
  # plot_df_db <- db_list$plot_df
  # features_db <- db_list$features
  
  
  return(list(
    model = model,
    X_train = X_train,
    
    train_conf = train_conf,
    test_conf = test_conf,
    train_pred = train_pred,
    test_pred = test_pred,
    
    train_accuracy = train_conf$overall["Accuracy"],
    test_accuracy  = test_conf$overall["Accuracy"],
    
    train_kappa = train_conf$overall["Kappa"],
    test_kappa  = test_conf$overall["Kappa"],
    
    sensitivity = test_conf$byClass["Sensitivity"],
    specificity = test_conf$byClass["Specificity"],
    balanced_accuracy = test_conf$byClass["Balanced Accuracy"],
    
    roc = roc_prob$roc,
    probabilities = roc_prob$probabilities,
    
    cm = conf$cm,
    confusion_df = conf$cm_df,
    
    pca = pca_data$pca,
    pca_df = pca_data$pca_df,
    
    feature_importance = imp,
    
    permutation_importance = perm_df,
    
    descision_boundary = db_list
  ))
}

##-------------------------------------------------------------------------------##
##                            SAVING IMPORTANT DATA                              ##
##-------------------------------------------------------------------------------##


save_important_data <- function(pred_results, y_test, seed, name, seed_dir) {
  cat("Saving important data")
  # Important metrics
  metrics <- data.frame(
    seed = seed,
    dataset = name,
    
    train_accuracy = pred_results$train_accuracy,
    test_accuracy  = pred_results$test_accuracy,
    
    train_kappa = pred_results$train_kappa,
    test_kappa  = pred_results$test_kappa,
    
    sensitivity = pred_results$sensitivity,
    specificity = pred_results$specificity,
    balanced_accuracy = pred_results$balanced_accuracy,
    
    auc = as.numeric(pred_results$roc$auc)
  )
  
  # Feature importance
  write.csv(
    pred_results$feature_importance$importance,
    file.path(seed_dir, "feature_importance.csv"),
    row.names = FALSE
  )
  
  write.csv(
    pred_results$permutation_importance,
    file.path(seed_dir, "permutation_importance.csv"),
    row.names = FALSE
  )
  
  # Predictions
  write.csv(
    data.frame(
      truth = y_test,
      prediction = pred_results$test_pred
    ),
    file.path(seed_dir, "test_predictions.csv"),
    row.names = FALSE
  )
  return(metrics)
}


##-------------------------------------------------------------------------------##
##                                 VAR IMP                                       ##
##-------------------------------------------------------------------------------##

varImp_function <- function(model, nb_of_feat, X_train, X_test, y_train, y_test, seed) {
  imp <- varImp(model, scale = TRUE)
  imp_df <- imp$importance
  imp_df$Overall <- rowMeans(abs(imp_df))
  imp_df$feature <- rownames(imp_df)
  imp_df <- imp_df[order(imp_df$Overall, decreasing = TRUE), ]
  top_features <- head(imp_df$feature, nb_of_feat)
  X_train_small <- X_train[, top_features, drop = FALSE]
  X_test_small  <- X_test[, top_features, drop = FALSE]
  svm_Imp <- SVM_model_running("svmLinear", X_train_small, y_train, seed)
  pred_list <- evaluate_predictions(svm_Imp, X_train_small, y_train, X_test_small, y_test)
  
  train_pred <- pred_list$train_pred
  test_pred <- pred_list$test_pred
  train_conf <- pred_list$train_conf
  test_conf <- pred_list$test_conf
  
  return(list(
    model = svm_Imp,
    X_train = X_train_small,
    train_conf = train_conf,
    test_conf = test_conf,
    X_test = X_test_small
  ))
}

accuracy_function <- function(variable) {
  accuracy <- as.numeric(round(variable$overall["Accuracy"], 4))
  return(accuracy)
}

compare_varImp <- function(model, X_train, X_test, y_train, y_test,
                           list_amt_feat, seed, shap_amt_feature) {
  
  results <- vector("list", length(list_amt_feat))
  
  for (i in seq_along(list_amt_feat)) {
    
    number <- list_amt_feat[i]
    
    varImp_variables <- varImp_function(
      model, number,
      X_train, X_test,
      y_train, y_test,
      seed
    )
    
    results[[i]] <- list(
      amount_of_features = number,
      test_accuracy = accuracy_function(varImp_variables$test_conf),
      train_accuracy = accuracy_function(varImp_variables$train_conf),
      model = varImp_variables$model,
      X_train = varImp_variables$X_train,
      X_test = varImp_variables$X_test
    )
  }
  
  results_df <- tibble::tibble(
    amount_of_features = sapply(results, `[[`, "amount_of_features"),
    test_accuracy = sapply(results, `[[`, "test_accuracy"),
    train_accuracy = sapply(results, `[[`, "train_accuracy"),
    model = lapply(results, `[[`, "model"),
    X_train = lapply(results, `[[`, "X_train"),
    X_test = lapply(results, `[[`, "X_test")
  )
  
  max_acc <- max(results_df$test_accuracy, na.rm = TRUE)
  if (all(is.na(results_df$test_accuracy))) {
    stop("All accuracies are NA — model failed")
  }
  
  candidate_idx <- which(results_df$test_accuracy == max_acc)
  
  best_idx <- candidate_idx[
    which.min(results_df$amount_of_features[candidate_idx])
  ]
  
  print(results_df[best_idx,
                   c("amount_of_features",
                     "test_accuracy",
                     "train_accuracy")])
  
  return(list(
    best_model = results_df$model[[best_idx]],
    best_X_train = results_df$X_train[[best_idx]],
    best_X_test = results_df$X_test[[best_idx]]
  ))
}

##-------------------------------------------------------------------------------##
##                                  SHAP                                         ##
##-------------------------------------------------------------------------------##

# Kernal shap

shap_kernel_function <- function(model, X) {
  
  bg_X <- X[sample(nrow(X), min(70, nrow(X))), ]
  
  kernelshap(
    object = model,
    X = X,
    bg_X = bg_X,
    pred_fun = function(object, newdata) {
      p <- predict(object, newdata, type = "prob")
      p[["PREGNANT"]]
    }
  )
}

# SHAP results
shap_results <- function(shaps, seed_dir) {
  
  shap_importance <- colMeans(abs(shaps))
  
  shap_df <- data.frame(
    feature = colnames(shaps),
    importance = shap_importance
  )
  
  shap_df <- shap_df[order(shap_df$importance, decreasing = TRUE), ]
  
  print(head(shap_df, 10))
  
  png(
    filename = file.path(seed_dir, "shap_results.png"),
    width = 1200,
    height = 1000,
    res = 150
  )
  
  p <- ggplot(head(shap_df, 10), aes(x = reorder(feature, importance), y = importance)) +
    geom_col(fill = "steelblue") +
    coord_flip() +
    theme_minimal() +
    labs(title = "Top 10 SHAP features", x = "Gene", y = "Mean |SHAP|")
  
  print(p)
  
  dev.off()
}



top_10_shap_per_seed <-  function(name, seeds_to_run) {
  results <- lapply(seeds_to_run, function(seed) {
    
    base_dir <- file.path("final_results", name)
    seed_dir <- file.path(base_dir, sprintf("seed_%03d", seed))
    shap_file <- file.path(seed_dir, "best_shap_values.rds")
    
    best_shap_values <- readRDS(shap_file)
    S <- best_shap_values$S
    
    importance <- colMeans(abs(S))
    top10 <- sort(importance, decreasing = TRUE)[1:10]
    
    data.frame(
      seed = seed,
      rank = 1:10,
      feature = names(top10),
      mean_abs_shap = as.numeric(top10)
    )
  })
  results_df <- bind_rows(results)
  
  write.csv(results_df,
            file = file.path("final_results", name, "top10_shap.csv"),
            row.names = FALSE)
  
} 


amt_feature <- function(n) {
  vec <- seq(5, n, by = 2)
  
  if (tail(vec, 1) != n) {
    vec <- c(vec, n)
  }
  
  return(vec)
}


calculate_metric_means <- function(df) {
  
  # seed en dataset uitsluiten
  metric_cols <- setdiff(names(df), c("seed", "dataset"))
  
  # naar numeriek forceren
  df[metric_cols] <- lapply(df[metric_cols], function(x)
    as.numeric(as.character(x)))
  
  means <- sapply(df[metric_cols], mean, na.rm = TRUE)
  
  data.frame(
    metric = names(means),
    mean = means,
    row.names = NULL
  )
}


small_best_pipeline <- function(matrix, train_sample_amount, test_sample_amount, seeds_to_run, name, len_features) {
  
  base_dir <- file.path("feature_selection_two_pipeline_results", name)
  # base_dir <- file.path("rf_results", name)
  dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)
  cat("\nIteration over multiple seeds\n\n\n")
  results <- lapply(seeds_to_run, function(seed) { 
    seed_dir <- file.path(base_dir, sprintf("seed_%03d", seed))
    cat(glue("\ncurrently seed {seed}\n\n"))
    dir.create(seed_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Data split
    test_train_sets <- preprocessing(matrix, train_sample_amount, test_sample_amount, seed)
    
    X_train <- test_train_sets$X_train
    X_test  <- test_train_sets$X_test
    y_train <- test_train_sets$y_train
    y_test  <- test_train_sets$y_test
    
    # Model training
    svm_model_lin <- SVM_model_running("svmLinear", X_train, y_train, seed)
    # evaluate_predictions(svm_model_lin, X_train, y_train, X_test, y_test)
    
    
    
    # varImp + best model
    best_list <- compare_varImp(
      svm_model_lin,
      X_train, X_test,
      y_train, y_test,
      amt_feature(len_features),
      seed,
      SHAP_AMT_FEATURE
    )
    
    best_model <- best_list$best_model
    best_X_train <- best_list$best_X_train
    best_X_test <- best_list$best_X_test
    
    pred_results <- prediction_figures(best_model, best_X_train, best_X_test, y_train, y_test, seed_dir)
    metrics <- save_important_data(pred_results, y_test, seed, name, seed_dir)
    
    # SHAP
    shap_file <- file.path(seed_dir, "best_shap_values.rds")
    
    if (file.exists(shap_file)) {
      cat("\nSHAP file already exitsts\n")
      shap_values <- readRDS(shap_file)
      
    } else {
      shap_values <- shap_kernel_function(best_model, best_X_train)
      saveRDS(shap_values, shap_file)
    }
    shap_results(shap_values$S, seed_dir)
    
    # return object
    list(
      metrics = metrics,
      # model = svm_model_lin,
      predictions = data.frame(
        truth = y_test,
        pred = pred_results$test_pred
      ))
    })
  
    # Combine all seeds
    metrics_df <- dplyr::bind_rows(lapply(results, `[[`, "metrics"))
    
    write.csv(
      metrics_df,
      file.path(base_dir, "metrics_summary.csv"),
      row.names = FALSE
    )
    
    
    file <- file.path(base_dir, "metrics_summary.csv")
    M <- read.csv(file, header = TRUE)
    results <- calculate_metric_means(M)
    
    write.csv(
      results,
      file.path(base_dir, "summary.csv"),
      row.names = FALSE
    )
    
    create_top10_shap(base_dir)
    
}

load_top10_matrix <- function(result_name, matrix_file) {
  
  # Lees top 10 features
  top10 <- read.csv(file.path("final_results", result_name, "top10_shap.csv"))
  features <- unique(top10$feature)
  len_features <-  length(features)
  cat("Length of top features", len_features, "\n\n")
  
  # Lees matrix
  mat <- read.table(
    matrix_file,
    header = TRUE,
    sep = "\t",
    row.names = 1
  )
  
  # Filter op features
  mat <- mat[
    rownames(mat) %in% features,
    ,
    drop = FALSE
  ]
  
  return(mat)
}

create_top10_shap <- function(matrix_dir){
  
  shap_files <- list.files(
    matrix_dir,
    pattern = "best_shap_values\\.rds$",
    recursive = TRUE,
    full.names = TRUE
  )
  
  if(length(shap_files) == 0){
    stop("Geen best_shap_values.rds bestanden gevonden.")
  }
  
  # Top 10 per seed
  all_seed_top10 <- map_dfr(shap_files, function(file){
    
    shap <- readRDS(file)
    
    if(!inherits(shap, "kernelshap")){
      stop(paste(file, "is geen kernelshap object"))
    }
    
    seed <- basename(dirname(file))
    
    tibble(
      gene = colnames(shap$S),
      shap = colMeans(abs(shap$S))
    ) %>%
      arrange(desc(shap)) %>%
      slice_head(n = 10) %>%
      mutate(seed = seed)
  })
  
  # Alle unieke genen uit de top10's
  unique_genes <- all_seed_top10 %>%
    distinct(gene)
  
  # Gemiddelde SHAP over alle seeds
  top10 <- all_seed_top10 %>%
    group_by(gene) %>%
    summarise(
      mean_shap = mean(shap),
      sd_shap = sd(shap),
      n_seeds = n(),
      .groups = "drop"
    ) %>%
    arrange(desc(mean_shap)) %>%
    slice_head(n = 10)
  
  write_csv(all_seed_top10,
            file.path(matrix_dir, "all_seed_top10.csv"))
  
  write_csv(unique_genes,
            file.path(matrix_dir, "best_unique_genes.csv"))
  
  write_csv(top10,
            file.path(matrix_dir, "top10_shap.csv"))
  
  invisible(list(
    all_seed_top10 = all_seed_top10,
    unique_genes = unique_genes,
    top10 = top10
  ))
}


load_svm_matrix <- function(result_name, name, matrix_file) {
  file_path <- file.path("rf_results", result_name, 
                         paste0("svm_gem10features.txt"))
  
  # features inlezen
  features_df <- read.delim(file_path, header = TRUE, stringsAsFactors = FALSE)
  features <- unique(as.character(features_df$feature))
  len_features <-  length(features)
  cat("Length of top features", len_features, "\n\n")
  
  # Lees matrix
  mat <- read.table(
    matrix_file,
    header = TRUE,
    sep = "\t",
    row.names = 1
  )
  
  # Filter op features
  mat <- mat[
    rownames(mat) %in% features,
    ,
    drop = FALSE
  ]
  
  return(mat)
}

pheatmap_function <-  function(name) {
  lig_Endo_mat <- load_top10_matrix(
    result_name = name,
    matrix_file = VITRO_LIG_ENDO_REC_EMBRYO
  )
  
  mat <- as.matrix(lig_Endo_mat)
  mat_scaled <- t(scale(t(mat)))
  
  status <- pregnancy_status(lig_Endo_mat)$PregnancyStatus
  
  
  status <- trimws(status)
  status <- toupper(status)
  
  ord <- order(status)
  
  mat_scaled <- mat_scaled[, ord]
  status <- status[ord]
  
  annotation_col <- data.frame(
    PregnancyStatus = factor(status, levels = c("NOT_PREGNANT", "PREGNANT"))
  )
  
  rownames(annotation_col) <- colnames(mat_scaled)
  
  ann_colors <- list(
    PregnancyStatus = c(
      NOT_PREGNANT = "skyblue",
      PREGNANT = "firebrick"
    )
  )
  
  
  pheatmap(
    mat_scaled,
    annotation_col = annotation_col,
    annotation_colors = ann_colors,
    
    cluster_rows = TRUE,
    cluster_cols = FALSE,   # IMPORTANT: anders breekt je sorting
    
    show_colnames = FALSE,
    show_rownames = TRUE,
    
    fontsize_row = 6,
    main = "Vivo - Top features (sorted by Pregnancy Status)"
  )
  
}


full_pipeline <- function() {
  start_time <- Sys.time()
  cat("Filtering matrices on top10 best features\n")
  
  
  if (choice == "ALL" | choice == "VIVO") {
    # Vivo matrix
    vivo_lig_Endo_mat <- load_top10_matrix(
      result_name = "vivo_lig-endo",
      matrix_file = VIVO_LIG_ENDO_REC_EMBRYO
    )
    
    cat(glue("\nRUNNING SVM PIPELINE FOR {WHICH} matrices\n\n"))
    
    vivo_n <-  "vivo_lig-endo"
    vivo <- pregnancy_status(vivo_lig_Endo_mat)
    len_features_vivo <- ncol(vivo)
    
    small_best_pipeline(vivo, VIVO_TRAIN_A, VIVO_TEST_A, SEEDS_TO_RUN, vivo_n, len_features_vivo)
    
    timestamp     <- format(Sys.time(), "%Y-%m-%d_%Hh%Mm%Ss")
    cat(glue("\n\nENDING FULL PIPELINE RUN FOR matrices at {timestamp}"))
    
    end_time <- Sys.time()
    
    total_time <- end_time - start_time
    cat(glue("Total time run {total_time}"))
    
    pheatmap_function(vivo_n)
  
  } else if (choice == "ALL" | choice == "VITRO"){
    # Vitro matrix
    vitro_lig_Endo_mat <- load_top10_matrix(
      result_name = "vitro_lig-endo",
      matrix_file = VITRO_LIG_ENDO_REC_EMBRYO
    )
    
    vitro_n <-  "vitro_lig-endo"
    vitro <- pregnancy_status(vitro_lig_Endo_mat)
    len_features_vitro <- ncol(vitro)
    
    
    small_best_pipeline(vitro, VITRO_TRAIN_A, VITRO_TEST_A, SEEDS_TO_RUN, vitro_n, len_features_vitro)
    
    timestamp     <- format(Sys.time(), "%Y-%m-%d_%Hh%Mm%Ss")
    cat(glue("\n\nENDING FULL PIPELINE RUN FOR matrices at {timestamp}"))
    
    end_time <- Sys.time()
    
    total_time <- end_time - start_time
    cat(glue("Total time run {total_time}"))
    
    pheatmap_function(vitro_n)
  }
  
}

full_pipeline()