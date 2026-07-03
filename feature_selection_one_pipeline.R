#------------------------------------------------------------------------------#
#                                PACKAGES                                      #
#------------------------------------------------------------------------------#

# Data handling & manipulation
#install.packages(c("data.table", "tibble", "glue"))
library(data.table)
library(tibble)
library(glue)
library(dplyr)

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

WHICH             <- "VIVO_VITRO_lig-embryo"

SEEDS_TO_RUN      <- c(64, 28, 21, 94, 41, 12, 53, 22, 17, 62)
# SEEDS_TO_RUN      <- c(17)

COMPARE_VAR_IMP   <- c(10, 15, 20, 25, 30, 40, 60, 80, 100, 120, 200, 300, 400)
SHAP_AMT_FEATURE  <- 10

VIVO_TRAIN_A      <- 7
VIVO_TEST_A       <- 4

VITRO_TRAIN_A     <- 5
VITRO_TEST_A      <- 4


##============================================================================##
##                          NORMALISATION PIPELINE                            ##
##============================================================================##

#------------------------------------------------------------------------------#
#                                LOADING DATA                                  #
#------------------------------------------------------------------------------#

load_data <- function(path, row_names = TRUE) {
  read.csv(
    path,
    header = TRUE,
    sep = "\t",
    row.names = if (row_names) 1 else NULL
  )
}


#------------------------------------------------------------------------------#
#                              TDM NORMALISATION                               #
#------------------------------------------------------------------------------#


tdm_aligment <- function(data_embryo, data_endom) {
  gemeenschappelijk <- intersect(rownames(data_embryo), rownames(data_endom))
  embryo_sub   <- data_embryo[gemeenschappelijk, ]
  endom_sub_v <- data_endom[gemeenschappelijk, ]
  
  endom_tdm_v <- data.table(gene = rownames(endom_sub_v), endom_sub_v, check.names = FALSE)
  embryo_tdm   <- data.table(gene = rownames(embryo_sub), embryo_sub,   check.names = FALSE)
  
  embryo_norm_tdm <- tdm_transform(ref_data = endom_tdm_v, target_data = embryo_tdm)
  
  embryo_norm <- as.matrix(embryo_norm_tdm[, -1])
  rownames(embryo_norm) <- embryo_norm_tdm$gene
  endom_norm_v <- as.matrix(endom_sub_v)
  
  return(list(
    embryo_norm = embryo_norm,
    endom_norm = endom_norm_v
  ))
}


#------------------------------------------------------------------------------#
#                         PCA AFTER TDM NORMALISATION                          #
#------------------------------------------------------------------------------#


run_pca_plot <- function(mat, groups, title) {
  
  mat <- mat[apply(mat, 1, sd) != 0, ]
  
  pca <- prcomp(t(mat), center = TRUE, scale. = TRUE)
  
  var_exp <- round(100 * summary(pca)$importance[2, 1:2], 1)
  
  scores <- data.frame(
    PC1 = pca$x[,1],
    PC2 = pca$x[,2],
    group = groups
  )
  
  ggplot(scores, aes(PC1, PC2, color = group)) +
    geom_point(size = 3) +
    theme_minimal() +
    labs(
      title = title,
      x = paste0("PC1 (", var_exp[1], "%)"),
      y = paste0("PC2 (", var_exp[2], "%)")
    )
}

save_pca <- function(matrix_norm, endom, name) {
  png(
    filename = file.path("feature_selection_one_pipeline_results", paste0(name,".png")),
    width = 1200,
    height = 1000,
    res = 150
  )
  
  p <- run_pca_plot(
    cbind(matrix_norm, endom),
    c(rep(name, ncol(matrix_norm)), rep("endom", ncol(endom))),
    glue("{name} - after TDM")
  )
  
  print(p)
  
  dev.off()
}

#------------------------------------------------------------------------------#
#                         FULL NORMALISATION PIPELINE                          #
#------------------------------------------------------------------------------#


normalisation_pipeline <- function() {
  cat("Beginning normalisation\n")
  dir.create("feature_selection_one_pipeline_results", showWarnings = FALSE, recursive = TRUE)
  # Loading data
  data_endom <- load_data(DATA_ENDOM)
  data_vivo  <- load_data(DATA_VIVO)
  data_vitro <- load_data(DATA_VITRO)
  data_lr <- load_data(DATA_LR, row_names = FALSE)
  
  # Vivo + Endom
  vivo_list <- tdm_aligment(data_vivo, data_endom)
  
  # Vitro + Endom
  vitro_list <- tdm_aligment(data_vitro, data_endom)
  
  # PCA before and after normalisation
  save_pca(vivo_list$embryo_norm, vivo_list$endom_norm, "vivo")
  save_pca(vitro_list$embryo_norm, vitro_list$endom_norm, "vitro")
  
  cat("Succesfull normalisation\n")
  
  return(list(
    vivo_norm = vivo_list$embryo_norm,
    vitro_norm = vitro_list$embryo_norm,
    endom_vivo = vivo_list$endom_norm,
    endom_vitro = vitro_list$endom_norm,
    data_lr = data_lr
  ))
  
}


##============================================================================##
##                             MATRICES PIPELINE                              ##
##============================================================================##

#------------------------------------------------------------------------------#
#                              MAKING THE MATRIX                               #
#------------------------------------------------------------------------------#


standardize_names <- function(x) {
  x <- gsub("nonPR", "NP", x)
  x
}

build_lr_matrix <- function(ligand_m,
                            receptor_m,
                            lr_table) {
  # Takes the unique ligand-receptor pairs
  data <- unique(
    lr_table[, c("ligand_ensembl", "receptor_ensembl")]
  )
  
  # Gives the ligand receptor pairs a name
  # ENSG1, ENSG123 --> ENSG1-ENSG123
  LR <- paste(
    data$ligand_ensembl,
    data$receptor_ensembl,
    sep = "-"
  )
  
  # Gives all possible sample pairs a name
  # SampleEmbryo01, SampleEndo01 --> SampleEmbryo01-SampleEndo01
  
  sample_name_pairs <- expand.grid(
    lig = colnames(ligand_m),
    rec = colnames(receptor_m),
    stringsAsFactors = FALSE
  )
  
  sample_name_pairs$pair <- paste(
    sample_name_pairs$lig,
    sample_name_pairs$rec,
    sep = "-"
  )
  
  sample_name_pairs$type <- ifelse(
    (grepl("NP", sample_name_pairs$lig) & grepl("PR", sample_name_pairs$rec)) |
      (grepl("PR", sample_name_pairs$lig) & grepl("NP", sample_name_pairs$rec)),
    "INVALID",
    "VALID"
  )
  
  valid_pairs <- sample_name_pairs[sample_name_pairs$type == "VALID", ]
  
  sample_names <- valid_pairs$pair
  
  col_index <- setNames(seq_along(sample_names), sample_names)
  
  valid_lookup <- setNames(rep(TRUE, length(sample_names)), sample_names)
  
  
  # Makes an empty matrix
  M <- matrix(
    NA_real_,
    nrow = nrow(data),
    ncol = length(sample_names),
    dimnames = list(LR, sample_names)
  )
  
  
  
  # Calculates the interaction score for echt ligand-receptor pair
  
  for (i in seq_len(nrow(data))) {
    
    ligand <- data$ligand_ensembl[i]
    receptor <- data$receptor_ensembl[i]
    
    # skip als genen niet bestaan in matrices
    if (!(ligand %in% rownames(ligand_m))) next
    if (!(receptor %in% rownames(receptor_m))) next
    
    # expression ophalen (1x per ligand-receptor pair per sample-combo)
    lig_expr_all <- ligand_m[ligand, , drop = TRUE]
    rec_expr_all <- receptor_m[receptor, , drop = TRUE]
    
    for (lig_col in names(lig_expr_all)) {
      
      lig_expr <- lig_expr_all[lig_col]
      
      for (rec_col in names(rec_expr_all)) {
        
        rec_expr <- rec_expr_all[rec_col]
        
        # sample pair naam
        pair_name <- paste(lig_col, rec_col, sep = "-")
        
        # alleen valid pairs gebruiken
        if (is.na(valid_lookup[pair_name])) next
        
        # kolom index opzoeken
        idx <- col_index[[pair_name]]
        
        # safety check
        if (is.null(idx)) next
        
        # score berekenen
        M[i, idx] <- lig_expr * rec_expr
      }
    }
  }
  M <- checks(M)
  return(M)
}

matrix_sanity_check <- function(matrix) {
  cat(
    "Totaal aantal cellen:", length(matrix), "\n",
    "NA:", sum(is.na(matrix)), "\n",
    "Niet-NA:", sum(!is.na(matrix)), "\n",
    "Percentage NA:", round(100 * mean(is.na(matrix)), 2), "%\n",
    "Aantal volledige NA-kolommen:", sum(colSums(!is.na(matrix)) == 0), "\n",
    "Aantal volledige NA-rijen:", sum(rowSums(!is.na(matrix)) == 0),
    "Aantal rijnamen:", nrow(matrix),
    "Aantal unieke rijnamen:", length(unique(rownames(matrix))),
    "Aantal kolomnamen:", ncol(matrix),
    "Aantal unieke kolomnamen", length(unique(colnames(matrix))), "\n\n"
  )
}

# Remove columns or rows which are all NA's
filteren_matrix <- function(matrix) {
  matrix2 <- matrix[
    rowSums(!is.na(matrix)) > 0,
    colSums(!is.na(matrix)) > 0,
    drop = FALSE
  ]
  return(matrix2)
}

# Control if all pairs are PR-PR or NP-NP
control_PR_NP <- function(matrix) {
  PR_pairs <- colnames(matrix)[
    (grepl("PR", colnames(matrix)))
  ]
  
  NP_pairs <- colnames(matrix)[
    (grepl("NP", colnames(matrix)))
  ]
  
  
  cat("length PR pairs:", length(PR_pairs), "\n")
  head(PR_pairs)
  
  cat("length NP pairs:", length(NP_pairs), "\n")
  head(NP_pairs)
  
  cat("length total pairs:", length(colnames(matrix)), "\n")
  head(colnames(matrix))
  cat("sum of PR and NP pairs", length(PR_pairs)+length(NP_pairs))
}

# Building the matrix
checks <- function(M) {
  matrix_sanity_check(M)
  M <- filteren_matrix(M)
  matrix_sanity_check(M)
  control_PR_NP(M)
  return(M)
}

make_matrix_pipeline <- function(datasets) {
  cat("Beginning designing the matrix\n")
  vivo_norm <- as.matrix(datasets$vivo_norm)
  vitro_norm <- as.matrix(datasets$vitro_norm)
  endom_vivo <- as.matrix(datasets$endom_vivo)
  endom_vitro <- as.matrix(datasets$endom_vitro)
  colnames(endom_vivo) <- standardize_names(colnames(endom_vivo))
  colnames(endom_vitro) <- standardize_names(colnames(endom_vitro))
  # Vivo lig-Embryo rec-Endometrium
  cat("\n\nMaking Vivo lig-Embryo rec-Endometrium\n")
  vivo_lig_Embryo_matrix <- build_lr_matrix(
    ligand_m = vivo_norm, 
    receptor_m = endom_vivo,
    lr_table = datasets$data_lr
  )
  write.table(
    vivo_lig_Embryo_matrix,
    file = VIVO_LIG_EMBRYO_REC_ENDO,
    sep = "\t",
    quote = FALSE,
    col.names = NA
  )
  
  # Vivo lig-Endometrium rec-Embryo
  cat("\n\nMaking Vivo lig-Endometrium rec-Embryo\n")
  vivo_lig_Endo_matrix <- build_lr_matrix(
    ligand_m = endom_vivo, 
    receptor_m = vivo_norm,
    lr_table = datasets$data_lr
  )
  write.table(
    vivo_lig_Endo_matrix,
    file = VIVO_LIG_ENDO_REC_EMBRYO,
    sep = "\t",
    quote = FALSE,
    col.names = NA
  )
  
  # Vitro lig-Embryo rec-Endometrium
  cat("\n\nMaking Vitro lig-Embryo rec-Endometrium\n")
  vitro_lig_Embryo_matrix <- build_lr_matrix(
    ligand_m = vitro_norm, 
    receptor_m = endom_vitro,
    lr_table = datasets$data_lr
  )
  write.table(
    vitro_lig_Embryo_matrix,
    file = VITRO_LIG_EMBRYO_REC_ENDO,
    sep = "\t",
    quote = FALSE,
    col.names = NA
  )
  
  # Vitro lig-Endometrium rec-Embryo
  cat("\n\nMaking Vitro lig-Endometrium rec-Embryo\n")
  vitro_lig_Endo_matrix <- build_lr_matrix(
    ligand_m = endom_vitro, 
    receptor_m = vitro_norm,
    lr_table = datasets$data_lr
  )
  write.table(
    vitro_lig_Endo_matrix,
    file = VITRO_LIG_ENDO_REC_EMBRYO,
    sep = "\t",
    quote = FALSE,
    col.names = NA
  )
  
  cat("Designing the matrix was succesfull")
  
  return(list(
    vivo_lig_Embryo_matrix = vivo_lig_Embryo_matrix,
    vivo_lig_Endo_matrix = vivo_lig_Endo_matrix,
    vitro_lig_Embryo_matrix = vitro_lig_Embryo_matrix,
    vitro_lig_Endo_matrix = vitro_lig_Endo_matrix
  ))
}

##============================================================================##
##                               SVM MODEL                                    ##
##============================================================================##

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
  cat("Starting preprocessing\n")
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


save_important_data <- function(basic_pred_list, pred_results, y_test, seed, name, seed_dir) {
  # Important metrics
  metrics <- data.frame(
    seed = seed,
    dataset = name,
    
    basic_train_pred <- basic_pred_list$train_conf$overall["Accuracy"],
    basic_test_conf <- basic_pred_list$test_conf$overall["Accuracy"],
    
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
  X_test <- read.csv(file.path(seed_dir, "X_test.csv"))
  write.csv(
    data.frame(
      samples = X_test$Sample,
      truth = y_test,
      prediction = pred_results$test_pred
    ),
    file.path(seed_dir, "test_predictions.csv"),
    row.names = FALSE
  )
  return(metrics)
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
  
  max_acc <- max(results_df$test_accuracy)
  
  candidate_idx <- which(results_df$test_accuracy == max_acc)
  
  best_idx <- candidate_idx[
    which.min(results_df$amount_of_features[candidate_idx])
  ]
  
  shap_idx <- which(results_df$amount_of_features == shap_amt_feature)
  
  print(results_df[best_idx,
                   c("amount_of_features",
                     "test_accuracy",
                     "train_accuracy")])
  
  return(list(
    best_model = results_df$model[[best_idx]],
    best_X_train = results_df$X_train[[best_idx]],
    best_X_test = results_df$X_test[[best_idx]],
    shap_model = results_df$model[[shap_idx]],
    shap_X_train = results_df$X_train[[shap_idx]]
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
    
    base_dir <- file.path("feature_selection_one_pipeline_results", name)
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
            file = file.path("feature_selection_one_pipeline_results", name, "top10_shap.csv"),
            row.names = FALSE)
}



##-------------------------------------------------------------------------------##
##                                PIPELINE                                       ##
##-------------------------------------------------------------------------------##

pipeline <- function(matrix, train_sample_amount, test_sample_amount, seeds_to_run, name) {
  
  base_dir <- file.path("feature_selection_one_pipeline_results", name)
  dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)
  cat("\nIteration over multiple seeds\n")
  gene_file <- file.path(base_dir, "best_genes.csv")
  write.csv(data.frame(Seed = integer(), Gene = character(), Accuracy= numeric()), gene_file, row.names = FALSE)
  
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
    
    # Model
    ## Radial has a worse predictions then linear
    # svm_model_rad <- SVM_model_running("svmRadial", X_train, y_train)
    svm_model_lin <- SVM_model_running("svmLinear", X_train, y_train, seed)
    
    basic_pred_list <- evaluate_predictions(svm_model_lin, X_train, y_train, X_test, y_test)

    basic_train_pred <- basic_pred_list$train_pred
    basic_test_pred <- basic_pred_list$test_pred
    basic_train_conf <- basic_pred_list$train_conf
    basic_test_conf <- basic_pred_list$test_conf
    
    ## permutations, slower and performs worse then varImp
    # perm_df <- permutation_importance(svm_model_lin, X_train, y_train)
    # top_features <- perm_df$feature[1:25]
    # X_train_small <- X_train[, top_features, drop = FALSE]
    # X_test_small  <- X_test[, top_features, drop = FALSE]
    # svm_permutation <- SVM_model_running("svmLinear", X_train_small, y_train)
    # evaluate_predictions(svm_permutation, X_train_small, y_train, X_test_small, y_test)
    
    ## AUC did not perform well
    # perm_df <- permutation_importance_auc(svm_model_lin, X_train, y_train)
    # top_features <- perm_df$feature[1:10]
    # X_train_small <- X_train[, top_features, drop = FALSE]
    # X_test_small  <- X_test[, top_features, drop = FALSE]
    # svm_permutation <- SVM_model_running("svmLinear", X_train_small, y_train)
    # evaluate_predictions(svm_permutation, X_train_small, y_train, X_test_small, y_test)
    
    # varImp + best model
    best_list <- compare_varImp(
      svm_model_lin,
      X_train, X_test,
      y_train, y_test,
      COMPARE_VAR_IMP,
      seed,
      SHAP_AMT_FEATURE
    )
    
    best_model <- best_list$best_model
    best_X_train <- best_list$best_X_train
    best_X_test <- best_list$best_X_test
    
    
    ## Saving all the important information
    saveRDS(
      best_model,
      file.path(seed_dir, "best_model.rds")
    )
    
    best_X_train_out <- cbind(Sample = rownames(best_X_train), best_X_train)
    write.csv(
      best_X_train_out,
      file.path(seed_dir, "X_train.csv"),
      row.names = FALSE
    )
    
    best_X_test_out <- cbind(Sample = rownames(best_X_test), best_X_test)
    write.csv(
      best_X_test_out,
      file.path(seed_dir, "X_test.csv"),
      row.names = FALSE
    )
    
    # Evaluation + plots + features
    pred_results <- prediction_figures(best_model, best_X_train, best_X_test, y_train, y_test, seed_dir)
    
    metrics <- save_important_data(basic_pred_list, pred_results, y_test, seed, name, seed_dir)
    
    gene_list <- data.frame(
      Seed = seed,
      Gene = colnames(best_X_train),
      Accuracy = accuracy_function(pred_results$test_conf)
    )
    
    write.table(gene_list,
                gene_file,
                sep = ",",
                row.names = FALSE,
                col.names = FALSE,  # geen header opnieuw schrijven
                append = TRUE)


    # SHAP
    shap_file <- file.path(seed_dir, "best_shap_values.rds")

    if (file.exists(shap_file)) {
      shap_values <- readRDS(shap_file)

    } else {
      shap_values <- shap_kernel_function(best_model, best_X_train)
      saveRDS(shap_values, shap_file)
    }
    shap_results(shap_values$S, seed_dir)
    return(list(
      metrics = metrics,
      model = svm_model_lin,
      shap = shap_values,
      predictions = data.frame(
        truth = y_test,
        pred = pred_results$test_pred
      )
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
  
  return(results)
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



full_pipeline <- function() {
  matrix_file <- file.path(VITRO_LIG_EMBRYO_REC_ENDO)
  # matrix_file <- "pyhton.py"
  timestamp     <- format(Sys.time(), "%Y-%m-%d_%Hh%Mm%Ss")
  print(glue("RUNNING FULL PIPELINE FOR matrices at {timestamp}"))
  
  if (file.exists(matrix_file)) {
    cat("Files already exist\n")
    vivo_lig_Embryo_matrix <- read.table(VIVO_LIG_EMBRYO_REC_ENDO, header = TRUE, sep = "\t", row.names = 1)
    vivo_lig_Endo_matrix <- read.table(VIVO_LIG_ENDO_REC_EMBRYO, header = TRUE, sep = "\t", row.names = 1)
    vitro_lig_Embryo_matrix <- read.table(VITRO_LIG_EMBRYO_REC_ENDO, header = TRUE, sep = "\t", row.names = 1)
    vitro_lig_Endo_matrix <- read.table(VITRO_LIG_ENDO_REC_EMBRYO, header = TRUE, sep = "\t", row.names = 1)
  } else {
    cat("Making the files\n")
    # Normalisation
    datasets <- normalisation_pipeline()
    
    # Making the matrices
    matrices <- make_matrix_pipeline(datasets)
    vivo_lig_Embryo_matrix <- matrices$vivo_lig_Embryo_matrix
    vivo_lig_Endo_matrix <- matrices$vivo_lig_Endo_matrix
    vitro_lig_Embryo_matrix <- matrices$vitro_lig_Embryo_matrix
    vitro_lig_Endo_matrix <- matrices$vitro_lig_Endo_matrix
  }
  
  cat(glue("\nRUNNING SVM PIPELINE FOR {WHICH} matrices\n\n"))

  if (WHICH == "ALL" | WHICH == "VIVO_VITRO_lig-embryo") {
    vivo_n <-  "vivo_lig-embryo"
    vitro_n <-  "vitro_lig-embryo"

    vivo <- pregnancy_status(vivo_lig_Embryo_matrix)
    vitro <- pregnancy_status(vitro_lig_Embryo_matrix)

    pipeline(vivo, VIVO_TRAIN_A, VIVO_TEST_A, SEEDS_TO_RUN, vivo_n)
    pipeline(vitro, VITRO_TRAIN_A, VITRO_TEST_A, SEEDS_TO_RUN, vitro_n)

    top_10_shap_per_seed(vivo_n, SEEDS_TO_RUN)
    top_10_shap_per_seed(vitro_n, SEEDS_TO_RUN)


  } else if (WHICH == "ALL" | WHICH == "VIVO_VITRO_lig-endo") {
    vivo_n <-  "vivo_lig-endo"
    vitro_n <-  "vitro_lig-endo"

    vivo <- pregnancy_status(vivo_lig_Endo_matrix)
    vitro <- pregnancy_status(vitro_lig_Endo_matrix)

    pipeline(vivo, VIVO_TRAIN_A, VIVO_TEST_A, SEEDS_TO_RUN, vivo_n)
    pipeline(vitro, VITRO_TRAIN_A, VITRO_TEST_A, SEEDS_TO_RUN, vitro_n)

    top_10_shap_per_seed(vivo_n, SEEDS_TO_RUN)
    top_10_shap_per_seed(vitro_n, SEEDS_TO_RUN)
  }

  timestamp     <- format(Sys.time(), "%Y-%m-%d_%Hh%Mm%Ss")
  cat(glue("\n\nENDING FULL PIPELINE RUN FOR matrices at {timestamp}"))
  
  
}


full_pipeline()