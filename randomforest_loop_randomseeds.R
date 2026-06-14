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



# function for training, testing and evaluating different datasets on various seeds
run_rf_experiment <- function(data, dataset_name, seeds){
  
  results <- data.frame()
  
  for(seed in seeds){
    
    cat("Running seed:", seed, "\n")
    
    set.seed(seed)
    
    
    # train/test split
    train_index <- createDataPartition(
      data$PregnancyStatus,
      p = 0.7,
      list = FALSE
    )
    
    train_set <- data[train_index, ]
    test_set  <- data[-train_index, ]
    
    
    # train model
    model <- train(
      PregnancyStatus ~ .,
      data = train_set,
      method = "ranger",
      importance = "impurity"
    )
    
    
    # save model RDS to folder
    saveRDS(
      model,
      paste0(
        "models/rf_",
        dataset_name,
        "_seed_",
        seed,
        ".rds"
      )
    )
    
    
    # make predictions on test_set
    pred_test <- predict(model, test_set)
    
    
    # confusion matrix
    cm <- confusionMatrix(
      pred_test,
      test_set$PregnancyStatus
    )
    
    
    # add run-statistics to results
    results <- rbind(
      results,
      data.frame(
        Dataset = dataset_name,
        Seed = seed,
        Accuracy = cm$overall["Accuracy"],
        Kappa = cm$overall["Kappa"]
      )
    )
  }
  
  
  return(results)
}


seeds <- c(64,28,21,94,41,12,53,22,17,62)

vivo_results <- run_rf_experiment(
  vivo_matrix2,
  "vivo",
  seeds
)

vitro_results <- run_rf_experiment(
  vitro_matrix2,
  "vitro",
  seeds
)


