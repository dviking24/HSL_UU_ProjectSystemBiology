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



seeds <- c(64,28,21,94,41,12,53,22,17,62)


run_rf_experiment <- function(data, dataset_name, seeds){
  
  results <- data.frame()
  
  for(seed in seeds){
    
    cat("Running seed:", seed, "\n")
    
    set.seed(seed)
    
    # do something
    
    results <- rbind(
      results,
      data.frame(
        Dataset = dataset_name,
        Seed = seed,
        Accuracy = "TEST"
      )
    )
  }
  
  return(results)
}



vivo_results <- run_rf_experiment(
  vivo_matrix2,
  "vivo",
  seeds
)

