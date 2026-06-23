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
library(glue)

# Which dataset to run
# Options are:
# VIVO_VITRO_lig-embryo
# VIVO_VITRO_lig-endo
# ALL
WHICH <- "VIVO_VITRO_lig-endo"

SEEDS_TO_RUN <- c(64, 28, 21, 94, 41, 12, 53, 22, 17, 62)
SEED <- 64

COMPARE_VAR_IMP <- c(10, 15, 20, 25, 30, 40)

VIVO_TRAIN_A <- 7
VIVO_TEST_A <- 4

VITRO_TRAIN_A <- 5
VITRO_TEST_A <- 4

VIVO_LIG_ENDO_REC_EMBRYO <- "matrices/vivo_lig-Endo_rec-Embryo_matrix.tsv"
VITRO_LIG_ENDO_REC_EMBRYO <- "matrices/vitro_lig-Endo_rec-Embryo_matrix.tsv"

VIVO_LIG_EMBRYO_REC_ENDO <- "matrices/vivo_lig-Embryo_rec-Endo_matrix.tsv"
VITRO_LIG_EMBRYO_REC_ENDO <- "matrices/vitro_lig-Embryo_rec-Endo_matrix.tsv"



##-------------------------------------------------------------------------------##
##                             LOADING MATRIX DATA                               ##
##-------------------------------------------------------------------------------##

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
load_matrices <- function(vivo_file_name, vitro_file_name) {
  cat("\nLoading vivo and vitro matrices\n")
  vivo_path  <- file.path(vivo_file_name)
  vitro_path <- file.path(vitro_file_name)
  
  vivo  <- read.table(vivo_path, header = TRUE, sep = "\t", row.names = 1)
  vitro <- read.table(vitro_path, header = TRUE, sep = "\t", row.names = 1)
  
  vivo <- pregnancy_status(vivo)
  vitro <- pregnancy_status(vitro)
 
  return(
    list(vivo = vivo, 
         vitro = vitro)
  )
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


preprocessing <- function(matrix, TRAIN_N, TEST_N, seed) {
  set.seed(seed)
  entities <- extract_entities(matrix)
  
  split_info <- sample_entities(
    embryos = entities$embryos,
    endos   = entities$endos,
    train_n = TRAIN_N,
    test_n  = TEST_N,
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
  
  return(list(
    roc = roc_obj,
    probabilities = probs
  ))
}


# Confusion Matrix
plot_confusion_matrix <- function(model, X_test, y_test, test_pred) {
  cm <- confusionMatrix(test_pred, y_test)
  cm_df <- as.data.frame(cm$table)
  
  ggplot(cm_df,
         aes(Prediction, Reference, fill = Freq)) +
    geom_tile() +
    geom_text(aes(label = Freq)) +
    scale_fill_gradient(low = "white", high = "red") +
    theme_minimal() +
    ggtitle("Confusion Matrix")
  
  return(list(
    cm = cm,
    cm_df = cm_df
  ))
}


# PCA of the testset data
plot_pca <- function(X_test, y_test, test_pred) {
  pca <- prcomp(X_test, scale. = TRUE)
  
  pca_df <- data.frame(
    PC1 = pca$x[,1],
    PC2 = pca$x[,2],
    label = y_test,
    pred = test_pred
  )
  
  ggplot(
    pca_df,
    aes(PC1, PC2, color = label)
  ) +
    geom_point(size = 3) +
    theme_minimal() +
    ggtitle("PCA")
  
  return(list(
    pca = pca,
    pca_df = pca_df
  ))
}


# Feature Importance 
plot_feature_importance <- function(model) {
  imp <- varImp(model)
  plot(
    imp,
    top = 25
  )
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
  
  return(list(
    model = model_small,
    grid = grid,
    plot_df = plot_df,
    features = top_features
  ))
}


prediction_figures <- function(model, X_train, X_test, y_train, y_test) {
  pred_list <- evaluate_predictions(model, X_train, y_train, X_test, y_test)
  train_pred <- pred_list$train_pred
  test_pred <- pred_list$test_pred
  train_conf <- pred_list$train_conf
  test_conf <- pred_list$test_conf
  
  roc_prob <- plot_roc_curve(model, X_test, y_test)

  conf <- plot_confusion_matrix(model, X_test, y_test, test_pred)

  pca_data <- plot_pca(X_test, y_test, test_pred)

  imp <- plot_feature_importance(model)

  perm_df <- permutation_importance(model, X_test, y_test)
  plot_permutation_importance(perm_df)

  db_list <- plot_decision_boundary(X_train, y_train)
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
##                            PERMUTATION IMPORTANCE                             ##
##-------------------------------------------------------------------------------##


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

compare_varImp <- function(model, X_train, X_test, y_train, y_test, list_amt_feat, seed) {
  test_accuracies <- c()
  train_accuracies <- c()
  models <- list()
  X_trains <- list()
  X_tests <- list()
  for (number in list_amt_feat) {
    varImp_variables <- varImp_function(model, number, X_train, X_test, y_train, y_test, seed)
    
    train_accuracies <-  c(train_accuracies, accuracy_function(varImp_variables$train_conf))
    test_accuracies <- c(test_accuracies, accuracy_function(varImp_variables$test_conf))
    models[[length(models) + 1]] <- varImp_variables$model
    X_trains[[length(X_trains) + 1]] <- varImp_variables$X_train
    X_tests[[length(X_tests) + 1]] <- varImp_variables$X_test
    
  }
  
  accuracy_df <- data.frame(
    amount_of_features = list_amt_feat,
    test_accuracy = test_accuracies,
    train_accuracy = train_accuracies
  )
  
  max_acc <- max(accuracy_df$test_accuracy)
  candidate_idx <- which(accuracy_df$test_accuracy == max_acc)
  best_idx <- candidate_idx[
    which.min(accuracy_df$amount_of_features[candidate_idx])
  ]
  
  best_model <- models[[best_idx]]
  best_X_train <- X_trains[[best_idx]]
  best_X_test <- X_tests[[best_idx]]
  print(accuracy_df[best_idx, ])
  
  return(list(
    best_model = best_model,
    best_X_train = best_X_train,
    best_X_test = best_X_test
  ))
}

##-------------------------------------------------------------------------------##
##                                  SHAP                                         ##
##-------------------------------------------------------------------------------##

# Kernal shap

shap_kernel_function <- function(model, X) {
  kernelshap(
    object = model,
    X = X_small,
    bg_X = X_small,
    pred_fun = function(object, newdata) {
      predict(object, newdata, type = "prob")[, "PREGNANT"]
    }
  )
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


##-------------------------------------------------------------------------------##
##                                PIPELINE                                       ##
##-------------------------------------------------------------------------------##

pipeline <- function(matrix, train_sample_amount, test_sample_amount, seeds_to_run, name) {
  
  base_dir <- file.path("results_SVM", name)
  dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)
  cat("\nIteration over multiple seeds\n")
  
  results <- lapply(seeds_to_run, function(seed) {
    cat(glue("\ncurrently seed {seed}\n"))
    
    seed_dir <- file.path(base_dir, sprintf("seed_%03d", seed))
    dir.create(seed_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Data split
    test_train_sets <- preprocessing(matrix, train_sample_amount, test_sample_amount, seed)
    
    X_train <- test_train_sets$X_train
    X_test  <- test_train_sets$X_test
    y_train <- test_train_sets$y_train
    y_test  <- test_train_sets$y_test
    
    # Model
    ## Radial has a worse predictions then linear
    # svm_model_rad <- SVM_model_running("svmRadial", X_train, y_train)
    svm_model_lin <- SVM_model_running("svmLinear", X_train, y_train, seed)
    evaluate_predictions(svm_model_lin, X_train, y_train, X_test, y_test)
    
    ## permutations, slower and performs worse then varImp
    # perm_df <- permutation_importance(svm_model_lin, X_train, y_train)
    # top_features <- perm_df$feature[1:25]
    # X_train_small <- X_train[, top_features, drop = FALSE]
    # X_test_small  <- X_test[, top_features, drop = FALSE]
    # svm_permutation <- SVM_model_running("svmLinear", X_train_small, y_train)
    # prediction_figures(svm_permutation, X_train_small, X_test_small, y_train, y_test)
    
    ## AUC did not perform well
    # perm_df <- permutation_importance_auc(svm_model_lin, X_train, y_train)
    # top_features <- perm_df$feature[1:10]
    # X_train_small <- X_train[, top_features, drop = FALSE]
    # X_test_small  <- X_test[, top_features, drop = FALSE]
    # svm_permutation <- SVM_model_running("svmLinear", X_train_small, y_train)
    # prediction_figures(svm_permutation, X_train_small, X_test_small, y_train, y_test)
    
    
    # varImp + best model
    best_list <- compare_varImp(
      svm_model_lin,
      X_train, X_test,
      y_train, y_test,
      COMPARE_VAR_IMP,
      seed
    )
    
    best_model <- best_list$best_model
    best_X_train <- best_list$best_X_train
    best_X_test <- best_list$best_X_test
    
    
    ## Saving all the important information
    saveRDS(
      best_model,
      file.path(seed_dir, "model.rds")
    )
    
    # Evaluation + plots + features
    pred_results <- prediction_figures(best_model, best_X_train, best_X_test, y_train, y_test)
    
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
    
    write.csv(
      metrics,
      file.path(seed_dir, "metrics.csv"),
      row.names = FALSE
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
    
    # SHAP
    shap_values <- shap_kernel_function(best_model, best_X_train)
    
    saveRDS(shap_values,file.path(seed_dir, "shap_values.rds"))
    
    shap_results(shap_values$S, colnames(best_X_train))
    
    # return object
    list(
      metrics = metrics,
      model = svm_model_lin,
      shap = shap_values,
      predictions = data.frame(
        truth = y_test,
        pred = pred_results$test_pred
      )
    )
  })
  
  # Combine all seeds
  metrics_df <- dplyr::bind_rows(lapply(results, `[[`, "metrics"))
  
  write.csv(
    metrics_df,
    file.path(base_dir, "metrics_summary.csv"),
    row.names = FALSE
  )
  
  return(results)
}


run_pipeline_data <- function() {
  timestamp     <- format(Sys.time(), "%Y-%m-%d_%Hh%Mm%Ss")
  
  glue("RUNNING SVM PIPELINE FOR {WHICH} matrices at {timestamp}")
  
  if (WHICH == "ALL" | WHICH == "VIVO_VITRO_lig-embryo") {
    
    vivo_file_name <- VITRO_LIG_EMBRYO_REC_ENDO 
    vitro_file_name <- VITRO_LIG_EMBRYO_REC_ENDO
    
    matrices <- load_matrices(vivo_file_name, vitro_file_name)
    
    pipeline(matrices$vivo, VIVO_TRAIN_A, VIVO_TEST_A, SEEDS_TO_RUN, "vivo_lig-embryo")
    pipeline(matrices$vitro, VITRO_TRAIN_A, VITRO_TEST_A, SEEDS_TO_RUN, "vitro_lig-embryo")
    
  } else if (WHICH == "ALL" | WHICH == "VIVO_VITRO_lig-endo") {
    
    vivo_file_name <- VIVO_LIG_ENDO_REC_EMBRYO 
    vitro_file_name <- VITRO_LIG_ENDO_REC_EMBRYO
    
    matrices <- load_matrices(vivo_file_name, vitro_file_name)
    
    pipeline(matrices$vivo, VIVO_TRAIN_A, VIVO_TEST_A, SEED, "vivo_lig-endo")
    # pipeline(matrices$vitro, VITRO_TRAIN_A, VITRO_TEST_A, SEED, "vitro_lig-endo")
  }
}

run_pipeline_data()


