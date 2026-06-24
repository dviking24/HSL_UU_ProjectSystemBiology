#------------------------------------------------------------------------------#
#                                PACKAGES                                      #
#------------------------------------------------------------------------------#

# Data handling & manipulation
#install.packages(c("data.table", "tibble", "glue"))
library(data.table)
library(tibble)
library(glue)

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
# SEED <- 64

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
#                            LIGAND-RECEPTOR DATA                              #
#------------------------------------------------------------------------------#


get_ligand_receptor_sets <- function(lr_data) {
  
  liganden  <- unique(lr_data$ligand_ensembl)
  receptors <- unique(lr_data$receptor_ensembl)
  
  list(
    liganden = liganden,
    receptors = receptors,
    ligreceptors = unique(c(liganden, receptors))
  )
}


#------------------------------------------------------------------------------#
#                              TDM NORMALISATION                               #
#------------------------------------------------------------------------------#


tdm_normalize_pair <- function(target, reference, common_genes) {
  
  target_sub <- target[common_genes, , drop = FALSE]
  ref_sub    <- reference[common_genes, , drop = FALSE]
  
  target_tdm <- data.table::data.table(
    gene = rownames(target_sub),
    target_sub,
    check.names = FALSE
  )
  
  ref_tdm <- data.table::data.table(
    gene = rownames(ref_sub),
    ref_sub,
    check.names = FALSE
  )
  
  norm_tdm <- TDM::tdm_transform(
    ref_data = ref_tdm,
    target_data = target_tdm
  )
  
  mat <- as.matrix(norm_tdm[, -1])
  rownames(mat) <- norm_tdm$gene
  
  mat
}


#------------------------------------------------------------------------------#
#                   PCA BEFORE AND AFTER TDM NORMALISATION                     #
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
    filename = file.path("final_results", paste0(name,".png")),
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
  dir.create("final_results", showWarnings = FALSE, recursive = TRUE)
  # Loading data
  data_endom <- load_data(DATA_ENDOM)
  data_vivo  <- load_data(DATA_VIVO)
  data_vitro <- load_data(DATA_VITRO)
  data_lr    <- load_data(DATA_LR, row_names = FALSE)
  
  meta_endom <- load_data(META_ENDOM)
  meta_vivo  <- load_data(META_VIVO)
  meta_vitro <- load_data(META_VITRO)
  
  # Ligand-receptor data
  lr_sets <- get_ligand_receptor_sets(data_lr)
  
  # Vivo alignment
  common_vivo <- intersect(rownames(data_vivo), rownames(data_endom))
  vivo_norm <- tdm_normalize_pair(
    target = data_vivo,
    reference = data_endom,
    common_genes = common_vivo
  )
  endom_vivo <- data_endom[common_vivo, ]
  
  
  # Vitro aligment
  common_vitro <- intersect(rownames(data_vitro), rownames(data_endom))
  vitro_norm <- tdm_normalize_pair(
    target = data_vitro,
    reference = data_endom,
    common_genes = common_vitro
  )
  endom_vitro <- data_endom[common_vitro, ]
  
  
  # PCA before and after normalisation
  save_pca(vivo_norm, endom_vivo, "vivo")
  save_pca(vitro_norm, endom_vitro, "vitro")
  
  cat("Succesfull normalisation\n")
  
  return(list(
    vivo_norm = vivo_norm,
    vitro_norm = vivo_norm,
    endom_vivo = endom_vivo,
    endom_vitro = endom_vitro,
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

build_lr_matrix <- function(ligand_matr,
                            receptor_matr,
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
    lig = colnames(ligand_matr),
    rec = colnames(receptor_matr),
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
    if (!(ligand %in% rownames(ligand_matr))) next
    if (!(receptor %in% rownames(receptor_matr))) next
    
    # expression ophalen (1x per ligand-receptor pair per sample-combo)
    lig_expr_all <- ligand_matr[ligand, , drop = TRUE]
    rec_expr_all <- receptor_matr[receptor, , drop = TRUE]
    
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
  checks(M)
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
}

make_matrix_pipeline <- function(datasets) {
  cat("Beginning designing the matrix\n")
  vivo_norm <- as.matrix(datasets$vivo_norm)
  vitro_norm <- as.matrix(datasets$vivo_norm)
  endom_vivo <- as.matrix(datasets$vivo_norm)
  endom_vitro <- as.matrix(datasets$vivo_norm)
  # Vivo lig-Embryo rec-Endometrium
  cat("\nMaking Vivo lig-Embryo rec-Endometrium\n")
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
  cat("\nMaking Vivo lig-Endometrium rec-Embryo\n")
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
  cat("\nMaking Vitro lig-Embryo rec-Endometrium\n")
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
  cat("\nMaking Vitro lig-Endometrium rec-Embryo\n")
  vitro_lig_Endo_matrix <- build_lr_matrix(
    ligand_m = endom_vitro, 
    receptor_m = vitro_norm,
    lr_table = datasets$data_lr
  )
  write.table(
    vitro_lig_Endo_matrix,
    file = VITRO_LIG_EMBRYO_REC_ENDO,
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
##                                   SVM                                      ##
##============================================================================##

run_pipeline_data <- function() {
  timestamp     <- format(Sys.time(), "%Y-%m-%d_%Hh%Mm%Ss")
  
  print(glue("RUNNING SVM PIPELINE FOR {WHICH} matrices at {timestamp}"))
  
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
    
    # pipeline(matrices$vivo, VIVO_TRAIN_A, VIVO_TEST_A, SEEDS_TO_RUN, "vivo_lig-endo")
    pipeline(matrices$vitro, VITRO_TRAIN_A, VITRO_TEST_A, SEEDS_TO_RUN, "vitro_lig-endo")
  }
}




full_pipeline <- function() {
  matrix_file <- file.path(VITRO_LIG_EMBRYO_REC_ENDO)
  timestamp     <- format(Sys.time(), "%Y-%m-%d_%Hh%Mm%Ss")
  print(glue("RUNNING FULL PIPELINE FOR matrices at {timestamp}"))
  
  if (file.exists(matrix_file)) {
    cat("Files already exist\n")
    vivo_lig_Embryo_matrix <- read.delim2(VIVO_LIG_EMBRYO_REC_ENDO)
    vivo_lig_Endo_matrix <- read.delim2(VIVO_LIG_ENDO_REC_EMBRYO)
    vitro_lig_Embryo_matrix <- read.delim2(VITRO_LIG_EMBRYO_REC_ENDO)
    vitro_lig_Endo_matrix <- read.delim2(VITRO_LIG_ENDO_REC_EMBRYO)
  } else {
    cat("Making the files\n")
    # Normalisation
    datasets <- normalisation_pipeline()
    
    # Making the matrices
    matrices <- make_matrix_pipeline(datasets)
  }
  
  run_pipeline_data(matrices)
  
  timestamp     <- format(Sys.time(), "%Y-%m-%d_%Hh%Mm%Ss")
  print(glue("ENDING FULL PIPELINE RUN FOR matrices at {timestamp}"))
}

full_pipeline()














