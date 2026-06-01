#devtools::install_github("greenelab/TDM", build_vignettes = TRUE)
#install.packages("data.table")
library(TDM)
library(ggplot2)
library(data.table)

# core data
data_endom <- read.csv("Datasets/DataEndom.txt", header = TRUE, sep = "\t", row.names = 1)
data_blastovivo <- read.csv("Datasets/Data_BlastoIVV.txt", header = TRUE, sep = "\t", row.names = 1)
data_blastovitro <- read.csv("Datasets/Data_BlastoIVT.txt", header = TRUE, sep = "\t", row.names = 1)
data_ligrecep <- read.csv("Datasets/LRdb_bovine_ENSEMBL.txt", header = TRUE, sep = "\t")

# metadata
meta_endom <- read.csv("Datasets/SampleInfo_Endom.txt", header = TRUE, sep = "\t")
meta_blastovivo <- read.csv("Datasets/SampleInfo_BlastoIVV.txt", header = TRUE, sep = "\t")
meta_blastovitro <- read.csv("Datasets/SampleInfo_BlastoIVT.txt", header = TRUE, sep = "\t")

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
pca_check(as.matrix(cbind(vivo_sub,  endom_sub)),   groep_vivo,  "Vivo - voor TDM")
pca_check(as.matrix(cbind(vitro_sub, endom_sub_v)), groep_vitro, "Vitro - voor TDM")

# na TDM
pca_check(cbind(vivo_norm,  endom_norm),   groep_vivo,  "Vivo - na TDM")
pca_check(cbind(vitro_norm, endom_norm_v), groep_vitro, "Vitro - na TDM")

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


# Filteren naar liganden/receptors 
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


#indentifying rows and columns
head(rownames(df_endom_liganden_voor_vivo))
head(rownames(df_blastovivo_receptors))

#making the matrices 4 different

#1. Endometrium receptor ← Embryo ligand (vivo)

comm_vivo_embLig_endoRec <- matrix(
  0,
  nrow = ncol(vivo_norm),
  ncol = ncol(endom_norm),
  dimnames = list(
    colnames(vivo_norm),
    colnames(endom_norm)
  )
)

for(i in seq_len(nrow(data_ligrecep))) {
  
  ligand <- data_ligrecep$ligand_ensembl[i]
  receptor <- data_ligrecep$receptor_ensembl[i]
  
  if(!(ligand %in% rownames(vivo_norm))) next
  if(!(receptor %in% rownames(endom_norm))) next
  
  comm_vivo_embLig_endoRec <- comm_vivo_embLig_endoRec +
    outer(
      as.numeric(vivo_norm[ligand, ]),
      as.numeric(endom_norm[receptor, ]),
      "*"
    )
}

write.table(
  comm_vivo_embLig_endoRec,
  file = "matrices/comm_vivo_embLig_endoRec.txt",
  sep = "\t",
  quote = FALSE,
  row.names = TRUE,
  col.names = NA
)

#2. Endometrium ligand → Embryo receptor (vivo)

comm_vivo_endoLig_embRec <- matrix(
  0,
  nrow = ncol(endom_norm),
  ncol = ncol(vivo_norm),
  dimnames = list(
    colnames(endom_norm),
    colnames(vivo_norm)
  )
)

for(i in seq_len(nrow(data_ligrecep))) {
  
  ligand <- data_ligrecep$ligand_ensembl[i]
  receptor <- data_ligrecep$receptor_ensembl[i]
  
  if(!(ligand %in% rownames(endom_norm))) next
  if(!(receptor %in% rownames(vivo_norm))) next
  
  comm_vivo_endoLig_embRec <- comm_vivo_endoLig_embRec +
    outer(
      as.numeric(endom_norm[ligand, ]),
      as.numeric(vivo_norm[receptor, ]),
      "*"
    )
}

write.table(
  comm_vivo_endoLig_embRec,
  file = "matrices/comm_vivo_endoLig_embRec.txt",
  sep = "\t",
  quote = FALSE,
  row.names = TRUE,
  col.names = NA
)

# 3. Endometrium receptor ← Embryo ligand (vitro)

comm_vitro_embLig_endoRec <- matrix(
  0,
  nrow = ncol(vitro_norm),
  ncol = ncol(endom_norm_v),
  dimnames = list(
    colnames(vitro_norm),
    colnames(endom_norm_v)
  )
)

for(i in seq_len(nrow(data_ligrecep))) {
  
  ligand <- data_ligrecep$ligand_ensembl[i]
  receptor <- data_ligrecep$receptor_ensembl[i]
  
  if(!(ligand %in% rownames(vitro_norm))) next
  if(!(receptor %in% rownames(endom_norm_v))) next
  
  comm_vitro_embLig_endoRec <- comm_vitro_embLig_endoRec +
    outer(
      as.numeric(vitro_norm[ligand, ]),
      as.numeric(endom_norm_v[receptor, ]),
      "*"
    )
}

write.table(
  comm_vitro_embLig_endoRec,
  file = "matrices/comm_vitro_embLig_endoRec.txt",
  sep = "\t",
  quote = FALSE,
  row.names = TRUE,
  col.names = NA
)

# 4. Endometrium ligand → Embryo receptor (vitro)

comm_vitro_endoLig_embRec <- matrix(
  0,
  nrow = ncol(endom_norm_v),
  ncol = ncol(vitro_norm),
  dimnames = list(
    colnames(endom_norm_v),
    colnames(vitro_norm)
  )
)

for(i in seq_len(nrow(data_ligrecep))) {
  
  ligand <- data_ligrecep$ligand_ensembl[i]
  receptor <- data_ligrecep$receptor_ensembl[i]
  
  if(!(ligand %in% rownames(endom_norm_v))) next
  if(!(receptor %in% rownames(vitro_norm))) next
  
  comm_vitro_endoLig_embRec <- comm_vitro_endoLig_embRec +
    outer(
      as.numeric(endom_norm_v[ligand, ]),
      as.numeric(vitro_norm[receptor, ]),
      "*"
    )
}

write.table(
  comm_vitro_endoLig_embRec,
  file = "matrices/comm_vitro_endoLig_embRec.txt",
  sep = "\t",
  quote = FALSE,
  row.names = TRUE,
  col.names = NA
)


# Matrix                      | Ligand bron  | Receptor bron |
# --------------------------- | ------------ | ------------- |
# comm_vivo_embLig_endoRec    | embryo vivo  | endometrium   |
# comm_vivo_endoLig_embRec    | endometrium  | embryo vivo   |
# comm_vitro_embLig_endoRec   | embryo vitro | endometrium   |
# comm_vitro_endoLig_embRec   | endometrium  | embryo vitro  |



#making the matrices like belen wants

# 1. vivo-endometrium

expr_vivo <- cbind(endom_norm, vivo_norm)

comm_vivo_all <- matrix(
  0,
  nrow = ncol(expr_vivo),
  ncol = ncol(expr_vivo),
  dimnames = list(
    colnames(expr_vivo),
    colnames(expr_vivo)
  )
)

for(i in seq_len(nrow(data_ligrecep))) {
  
  ligand <- data_ligrecep$ligand_ensembl[i]
  receptor <- data_ligrecep$receptor_ensembl[i]
  
  if(!(ligand %in% rownames(expr_vivo))) next
  if(!(receptor %in% rownames(expr_vivo))) next
  
  comm_vivo_all <- comm_vivo_all +
    outer(
      as.numeric(expr_vivo[ligand, ]),
      as.numeric(expr_vivo[receptor, ]),
      "*"
    )
}

write.table(
  comm_vivo_all,
  file = "matrices/comm_vivo_all.txt",
  sep = "\t",
  quote = FALSE,
  row.names = TRUE,
  col.names = NA
)

# 2. vitro-endometrium
expr_vitro <- cbind(endom_norm_v, vitro_norm)

comm_vitro_all <- matrix(
  0,
  nrow = ncol(expr_vitro),
  ncol = ncol(expr_vitro),
  dimnames = list(
    colnames(expr_vitro),
    colnames(expr_vitro)
  )
)

for(i in seq_len(nrow(data_ligrecep))) {
  
  ligand <- data_ligrecep$ligand_ensembl[i]
  receptor <- data_ligrecep$receptor_ensembl[i]
  
  if(!(ligand %in% rownames(expr_vitro))) next
  if(!(receptor %in% rownames(expr_vitro))) next
  
  comm_vitro_all <- comm_vitro_all +
    outer(
      as.numeric(expr_vitro[ligand, ]),
      as.numeric(expr_vitro[receptor, ]),
      "*"
    )
}

write.table(
  comm_vitro_all,
  file = "matrices/comm_vitro_all.txt",
  sep = "\t",
  quote = FALSE,
  row.names = TRUE,
  col.names = NA
)



