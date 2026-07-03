#devtools::install_github("greenelab/TDM", build_vignettes = TRUE)
#install.packages("data.table")
#install.packages("glue")
library(glue)
library(TDM)
library(ggplot2)
library(data.table)
library(pheatmap)
library(png)

# core data
data_endom <- read.csv("Datasets/DataEndom.txt", header = TRUE, sep = "\t", row.names = 1)
data_blastovivo <- read.csv("Datasets/Data_BlastoIVV.txt", header = TRUE, sep = "\t", row.names = 1)
data_blastovitro <- read.csv("Datasets/Data_BlastoIVT.txt", header = TRUE, sep = "\t", row.names = 1)
data_ligrecep <- read.csv("Datasets/LRdb_bovine_ENSEMBL.txt", header = TRUE, sep = "\t")

# metadata
meta_endom <- read.csv("Datasets/SampleInfo_Endom.txt", header = TRUE, sep = "\t")
meta_blastovivo <- read.csv("Datasets/SampleInfo_BlastoIVV.txt", header = TRUE, sep = "\t")
meta_blastovitro <- read.csv("Datasets/SampleInfo_BlastoIVT.txt", header = TRUE, sep = "\t")

# --- Data inladen ---
data_endom       <- read.csv("Datasets/DataEndom.txt", header = TRUE, sep = "\t", row.names = 1)
data_blastovivo  <- read.csv("Datasets/Data_BlastoIVV.txt", header = TRUE, sep = "\t", row.names = 1)
data_blastovitro <- read.csv("Datasets/Data_BlastoIVT.txt", header = TRUE, sep = "\t", row.names = 1)
data_ligrecep    <- read.csv("Datasets/LRdb_bovine_ENSEMBL.txt", header = TRUE, sep = "\t")

# metadata
meta_endom       <- read.csv("Datasets/SampleInfo_Endom.txt", header = TRUE, sep = "\t")
meta_blastovivo  <- read.csv("Datasets/SampleInfo_BlastoIVV.txt", header = TRUE, sep = "\t")
meta_blastovitro <- read.csv("Datasets/SampleInfo_BlastoIVT.txt", header = TRUE, sep = "\t")





# Gen-lijsten samenstellen
liganden     <- unique(data_ligrecep$ligand_ensembl)
receptors    <- unique(data_ligrecep$receptor_ensembl)
ligreceptors <- intersect(liganden, receptors)

# Cross-data normalisatie met TDM 
# endom (microarray) is referentie, vivo/vitro (RNA-seq) worden naar die schaal getrokken

# Vivo + Endom
gemeenschappelijk_vivo <- intersect(rownames(data_blastovivo), rownames(data_endom))
vivo_sub  <- data_blastovivo[gemeenschappelijk_vivo, ]
endom_sub <- data_endom[gemeenschappelijk_vivo, ]

# Omzetten naar TDM-formaat (gene als kolom)
endom_tdm <- data.table(gene = rownames(endom_sub), endom_sub, check.names = FALSE)
vivo_tdm  <- data.table(gene = rownames(vivo_sub),  vivo_sub,  check.names = FALSE)

# TDM toepassen
vivo_norm_tdm <- tdm_transform(ref_data = endom_tdm, target_data = vivo_tdm)

# Terug naar matrix-formaat
vivo_norm <- as.matrix(vivo_norm_tdm[, -1])
rownames(vivo_norm) <- vivo_norm_tdm$gene
endom_norm <- as.matrix(endom_sub)   # referentie blijft ongewijzigd


# Vitro + Endom
gemeenschappelijk_vitro <- intersect(rownames(data_blastovitro), rownames(data_endom))
vitro_sub   <- data_blastovitro[gemeenschappelijk_vitro, ]
endom_sub_v <- data_endom[gemeenschappelijk_vitro, ]

endom_tdm_v <- data.table(gene = rownames(endom_sub_v), endom_sub_v, check.names = FALSE)
vitro_tdm   <- data.table(gene = rownames(vitro_sub),   vitro_sub,   check.names = FALSE)

vitro_norm_tdm <- tdm_transform(ref_data = endom_tdm_v, target_data = vitro_tdm)

vitro_norm <- as.matrix(vitro_norm_tdm[, -1])
rownames(vitro_norm) <- vitro_norm_tdm$gene
endom_norm_v <- as.matrix(endom_sub_v)


# PCA check voor en na 
pca_check <- function(mat, groep, titel) {
  mat <- mat[apply(mat, 1, function(x) sd(x) != 0), ]
  pca <- prcomp(t(mat), center = TRUE, scale. = TRUE)
  var <- round(100 * summary(pca)$importance[2, 1:2], 1)
  scores <- as.data.frame(pca$x[, 1:2])
  scores$groep <- groep
  ggplot(scores, aes(PC1, PC2, color = groep)) +
    geom_point(size = 3) +
    labs(title = titel, x = paste0("PC1 (", var[1], "%)"), y = paste0("PC2 (", var[2], "%)")) +
    theme_minimal()
}

n_vivo  <- ncol(vivo_sub)
n_vitro <- ncol(vitro_sub)
groep_vivo  <- c(rep("vivo",  n_vivo),  rep("endom", ncol(endom_sub)))
groep_vitro <- c(rep("vitro", n_vitro), rep("endom", ncol(endom_sub_v)))

# voor TDM
p1 <- pca_check(as.matrix(cbind(vivo_sub, endom_sub)), groep_vivo, "Vivo - before TDM")
p2 <- pca_check(as.matrix(cbind(vitro_sub, endom_sub_v)), groep_vitro, "Vitro - before TDM")

# na TDM
p3 <- pca_check(cbind(vivo_norm, endom_norm), groep_vivo, "Vivo - after TDM")
p4 <- pca_check(cbind(vitro_norm, endom_norm_v), groep_vitro, "Vitro - after TDM")

print(p1)
print(p2)
print(p3)
print(p4)

ggsave("PCA_vivo_voor_TDM.png", plot = p1, width = 6, height = 5, dpi = 300)
ggsave("PCA_vitro_voor_TDM.png", plot = p2, width = 6, height = 5, dpi = 300) 
ggsave("PCA_vivo_na_TDM.png", plot = p3, width = 6, height = 5, dpi = 300) 
ggsave("PCA_vitro_na_TDM.png", plot = p4, width = 6, height = 5, dpi = 300)

#ligand, receptor en ligrec

# unieke liganden en receptoren uit LR database

liganden <- unique(data_ligrecep$ligand_ensembl)

receptors <- unique(data_ligrecep$receptor_ensembl)

ligreceptors <- unique(c(liganden, receptors))

length(liganden)
length(receptors)
length(ligreceptors)

head(liganden)
head(receptors)


# Filteren naar liganden/receptors <-- dit wordt niet gebruikt
df_blastovivo_liganden     <- vivo_norm[rownames(vivo_norm) %in% liganden, ]
df_blastovivo_receptors    <- vivo_norm[rownames(vivo_norm) %in% receptors, ]
df_blastovivo_ligreceptors <- vivo_norm[rownames(vivo_norm) %in% ligreceptors, ]

df_endom_liganden_voor_vivo     <- endom_norm[rownames(endom_norm) %in% liganden, ]
df_endom_receptors_voor_vivo    <- endom_norm[rownames(endom_norm) %in% receptors, ]
df_endom_ligreceptors_voor_vivo <- endom_norm[rownames(endom_norm) %in% ligreceptors, ]

df_blastovitro_liganden     <- vitro_norm[rownames(vitro_norm) %in% liganden, ]
df_blastovitro_receptors    <- vitro_norm[rownames(vitro_norm) %in% receptors, ]
df_blastovitro_ligreceptors <- vitro_norm[rownames(vitro_norm) %in% ligreceptors, ]

df_endom_liganden_voor_vitro     <- endom_norm_v[rownames(endom_norm_v) %in% liganden, ]
df_endom_receptors_voor_vitro    <- endom_norm_v[rownames(endom_norm_v) %in% receptors, ]
df_endom_ligreceptors_voor_vitro <- endom_norm_v[rownames(endom_norm_v) %in% ligreceptors, ]

length(colnames(endom_norm))
length(colnames(vivo_norm))
#indentifying rows and columns
head(rownames(df_endom_liganden_voor_vivo))
head(rownames(df_blastovivo_receptors))



#//------------------------------------MAKING THE MATRICES--------------------------------------------------------//
#making the matrice

standardize_names <- function(x) {
  x <- gsub("nonPR", "NP", x)
  x
}

head(colnames(endom_norm))
colnames(endom_norm) <- standardize_names(colnames(endom_norm))
head(colnames(endom_norm))


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
  return(M)
}

# Sanity check voor de matrix, checken op NA's en structuur
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
    "Aantal unieke kolomnamen", length(unique(colnames(matrix)))
  )
}

# Kolommen en rijen filteren als ze volledig NA zijn
filteren_matrix <- function(matrix) {
  matrix2 <- matrix[
    rowSums(!is.na(matrix)) > 0,
    colSums(!is.na(matrix)) > 0,
    drop = FALSE
  ]
  return(matrix2)
}

# Vivo matrix

make_matrix <- function(ligand_m, receptor_m, naam_heatmap){
  #matrix maken
  matrix <- build_lr_matrix(
    ligand_matr = ligand_m,
    receptor_matr   = receptor_m,
    lr_table    = data_ligrecep
  )
  
  #sanity checks
  matrix_sanity_check(matrix)
  
  #filteren op NA kolommen en rijen
  matrix2 <- filteren_matrix(matrix)
  
  #sanity check
  matrix_sanity_check(matrix2)
  
  #Volledige heatmap maken
  pheatmap(
    matrix2, 
    na_col = "grey90",
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    filename = glue("figures/{naam_heatmap}/heatmap.pdf")
  )
  
  #Meest en minst variabele rijen laten zien
  row_var <- apply(matrix2, 1, var, na.rm = TRUE)
  
  #Meest en minst variabele kolommen laten zien
  col_var <- apply(matrix2, 2, var, na.rm = TRUE)
  
  top_rows <- names(sort(row_var, decreasing = TRUE))[1:100]
  bottom_rows <- names(sort(row_var, decreasing = FALSE))[1:100]
  
  top_cols <- names(sort(col_var, decreasing = TRUE))[1:100]
  bottom_cols <- names(sort(col_var, decreasing = FALSE))[1:100]
  
  pheatmap(
    matrix2[top_rows, ],
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    fontsize_row = 6,
    filename = glue("figures/{naam_heatmap}/top_rows_var_heatmap.pdf")
  )
  
  pheatmap(
    matrix2[bottom_rows, ],
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    fontsize_row = 6,
    filename = glue("figures/{naam_heatmap}/bottom_rows_var_heatmap.pdf")
  )
  
  pheatmap(
    matrix2[, top_cols],
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    fontsize_row = 6,
    filename = glue("figures/{naam_heatmap}/top_cols_var_heatmap.pdf")
  )
  
  pheatmap(
    matrix2[, bottom_cols],
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    fontsize_row = 6,
    filename = glue("figures/{naam_heatmap}/bottom_cols_var_heatmap.pdf")
  )
  
  return(matrix2)
}

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

# Grote vivo matrix maken

vivo_lig_matrix <- make_matrix(
  ligand_m = vivo_norm, 
  receptor_m = endom_norm,
  naam_heatmap = "vivo_lig"
)

vivo_rec_matrix <- make_matrix(
  ligand_m = endom_norm, 
  receptor_m = vivo_norm,
  naam_heatmap = "vivo_rec"
)

vivo_matrix <- cbind(vivo_lig_matrix, vivo_rec_matrix)

# pheatmap(
#   vivo_matrix,
#   cluster_rows = TRUE,
#   cluster_cols = TRUE,
#   fontsize_row = 6,
#   filename = glue("figures/vivo/heatmap.pdf")
# )

control_PR_NP(vivo_matrix)

# opslaan als matrices/vivo_matrix.tsv

write.table(
  vivo_matrix,
  file = "matrices/vivo_matrix.tsv",
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

# opslaan vivo ligand-embryo receptor-endometrium
write.table(
  vivo_lig_matrix,
  file = "matrices/vivo_lig-Embryo_rec-Endo_matrix.tsv",
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

# opslaan vivo ligand-endometrium receptor-embryo
write.table(
  vivo_rec_matrix,
  file = "matrices/vivo_lig-Endo_rec-Embryo_matrix.tsv",
  sep = "\t",
  quote = FALSE,
  col.names = NA
)













# Grote vitro matrix maken

vitro_lig_matrix <- make_matrix(
  ligand_m = vitro_norm, 
  receptor_m = endom_norm,
  naam_heatmap = "vitro_lig"
)

vitro_rec_matrix <- make_matrix(
  ligand_m = endom_norm, 
  receptor_m = vitro_norm,
  naam_heatmap = "vitro_rec"
)

vitro_matrix <- cbind(vitro_lig_matrix, vitro_rec_matrix)

# pheatmap(
#   vitro_matrix,
#   cluster_rows = TRUE,
#   cluster_cols = TRUE,
#   fontsize_row = 6,
#   filename = glue("figures/vitro/heatmap.pdf")
# )

control_PR_NP(vitro_matrix)

# Zijn er nog NP-PR paren aanwezig?
sum(
  grepl("NP", colnames(vitro_matrix)) &
    grepl("-.*PR", colnames(vitro_matrix))
)

# opslaan als matrices/vitro_matrix.tsv
write.table(
  vitro_matrix,
  file = "matrices/vitro_matrix.tsv",
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

# opslaan vitro ligand-embryo receptor-endometrium
write.table(
  vitro_lig_matrix,
  file = "matrices/vitro_lig-Embryo_rec-Endo_matrix.tsv",
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

# opslaan vitro ligand-endometrium receptor-embryo
write.table(
  vitro_rec_matrix,
  file = "matrices/vitro_lig-Endo_rec-Embryo_matrix.tsv",
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

# Meest en minst variabele rijen laten zien
# Meest en minst variabele rijen laten zien
# row_var <- apply(vivo_matrix, 1, var, na.rm = TRUE)
# 
# Meest en minst variabele kolommen laten zien
# col_var <- apply(vivo_matrix, 2, var, na.rm = TRUE)
# 
# top_rows <- names(sort(row_var, decreasing = TRUE))[1:40]
# bottom_rows <- names(sort(row_var, decreasing = FALSE))[1:40]
# 
# top_cols <- names(sort(col_var, decreasing = TRUE))[1:40]
# bottom_cols <- names(sort(col_var, decreasing = FALSE))[1:40]
# 
# 
# pheatmap(
#   vitro_matrix[top_rows, top_cols],
#   cluster_rows = TRUE,
#   cluster_cols = FALSE,
#   fontsize_row = 6,
#   fontsize_col = 6,
#   filename = glue("figures/vitro/top_rows_top_cols_var_heatmap.pdf")
# )
# 
# top_40_vitro <- vitro_matrix[top_rows, top_cols]
# top_40_vivo <- vivo_matrix[top_rows, top_cols]
# View(top_40_vitro)
# View(top_40_vivo)
# 
# View(data_ligrecep)


