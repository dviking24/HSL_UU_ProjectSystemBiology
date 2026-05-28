# --- load data --- 
# core data
data_endom <- read.csv("Datasets/DataEndom.txt", header = TRUE, sep = "\t", row.names = 1)
data_blastovivo <- read.csv("Datasets/Data_BlastoIVV.txt", header = TRUE, sep = "\t", row.names = 1)
data_blastovitro <- read.csv("Datasets/Data_BlastoIVT.txt", header = TRUE, sep = "\t", row.names = 1)
data_ligrecep <- read.csv("Datasets/LRdb_bovine_ENSEMBL.txt", header = TRUE, sep = "\t")

# metadata
meta_endom <- read.csv("Datasets/SampleInfo_Endom.txt", header = TRUE, sep = "\t")
meta_blastovivo <- read.csv("Datasets/SampleInfo_BlastoIVV.txt", header = TRUE, sep = "\t")
meta_blastovitro <- read.csv("Datasets/SampleInfo_BlastoIVT.txt", header = TRUE, sep = "\t")


# Genen inladen
liganden <- unique(data_ligrecep$ligand_ensembl)
receptors <- unique(data_ligrecep$receptor_ensembl)
ligreceptors <- intersect(liganden, receptors)
genen_waar_we_naar_zoeken <- c(liganden, receptors)
alle_genen_in_vivo <- row.names(data_blastovivo)
alle_genen_in_vitro <- row.names(data_blastovitro)
alle_genen_in_endom <- row.names(data_endom)


# Filteren van genen in blastovivo-dataset
onnodige_genen_in_vivo_data <- setdiff(alle_genen_in_vivo, genen_waar_we_naar_zoeken)
goede_genen_in_vivo_data <- setdiff(alle_genen_in_vivo, onnodige_genen_in_vivo_data)

# Blastovivo-dataframe filteren zodat het alleen nog maar rijen met goede genen overhoudt
data_blastovivo2 <- data_blastovivo[goede_genen_in_vivo_data, ]


# Filteren van genen in blastovitro-dataset
onnodige_genen_in_vitro_data <- setdiff(alle_genen_in_vitro, genen_waar_we_naar_zoeken)
goede_genen_in_vitro_data <- setdiff(alle_genen_in_vitro, onnodige_genen_in_vitro_data)

# Blastovitro-dataframe filteren zodat het alleen nog maar rijen met goede genen overhoudt
data_blastovitro2 <- data_blastovitro[goede_genen_in_vitro_data, ]


# Filteren van genen in endom-dataset
onnodige_genen_in_endom_data <- setdiff(alle_genen_in_endom, genen_waar_we_naar_zoeken)
goede_genen_in_endom_data <- setdiff(alle_genen_in_endom, onnodige_genen_in_endom_data)

# Endom-dataframe filteren zodat het alleen nog maar rijen met goede genen overhoudt
data_endom2 <- data_endom[goede_genen_in_endom_data, ]




library(preprocessCore)
library(ggplot2)

# genen die in beide zitten (vivo + endom)
gemeenschappelijk <- intersect(rownames(data_blastovivo2), rownames(data_endom2))

vivo_sub  <- data_blastovivo2[gemeenschappelijk, ]
endom_sub <- data_endom2[gemeenschappelijk, ]

gecombineerd <- as.matrix(cbind(vivo_sub, endom_sub))

# quantile normalisatie
genorm <- normalize.quantiles(gecombineerd)
dimnames(genorm) <- dimnames(gecombineerd)

# terugsplitsen
n_vivo <- ncol(vivo_sub)
vivo_norm  <- genorm[, 1:n_vivo]
endom_norm <- genorm[, (n_vivo + 1):ncol(genorm)]

# groep-label voor de PCA
groep <- c(rep("vivo", n_vivo), rep("endom", ncol(endom_sub)))

# --- PCA functie (zelfde als bij vitro) ---
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

# before en after
pca_check(gecombineerd, groep, "Vivo - voor normalisatie")
pca_check(genorm, groep, "Vivo - na normalisatie")




# genen die in beide zitten (vitro + endom)
gemeenschappelijk_v <- intersect(rownames(data_blastovitro2), rownames(data_endom2))

vitro_sub  <- data_blastovitro2[gemeenschappelijk_v, ]
endom_sub_v <- data_endom2[gemeenschappelijk_v, ]

gecombineerd_v <- as.matrix(cbind(vitro_sub, endom_sub_v))

# quantile normalisatie
genorm_v <- normalize.quantiles(gecombineerd_v)
dimnames(genorm_v) <- dimnames(gecombineerd_v)

# terugsplitsen
n_vitro <- ncol(vitro_sub)
vitro_norm <- genorm_v[, 1:n_vitro]
endom_norm_v <- genorm_v[, (n_vitro + 1):ncol(genorm_v)]

# groep-label voor de PCA
groep_v <- c(rep("vitro", n_vitro), rep("endom", ncol(endom_sub_v)))

# --- PCA functie zodat je before/after makkelijk draait ---
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

# before en after
pca_check(gecombineerd_v, groep_v, "Vitro - voor normalisatie")
pca_check(genorm_v, groep_v, "Vitro - na normalisatie")







# Dataframes van Liganden, Receptors en LigandReceptor genen in Vivo-dataset
liganden_vivo <- intersect(goede_genen_in_vivo_data, liganden)
receptors_vivo <- intersect(goede_genen_in_vivo_data, receptors)
ligreceptors_vivo <- intersect(goede_genen_in_vivo_data, ligreceptors)

df_blastovivo_liganden <- data_blastovivo2[liganden_vivo, ]
df_blastovivo_receptors <- data_blastovivo2[receptors_vivo, ]
df_blastovivo_ligreceptors <- data_blastovivo2[ligreceptors_vivo, ]


# Dataframes van Liganden, Receptors en LigandReceptor genen in Vitro-dataset
liganden_vitro <- intersect(goede_genen_in_vitro_data, liganden)
receptors_vitro <- intersect(goede_genen_in_vitro_data, receptors)
ligreceptors_vitro <- intersect(goede_genen_in_vitro_data, ligreceptors)

df_blastovitro_liganden <- data_blastovitro2[liganden_vitro, ]
df_blastovitro_receptors <- data_blastovitro2[receptors_vitro, ]
df_blastovitro_ligreceptors <- data_blastovitro2[ligreceptors_vitro, ]


# Dataframes van Liganden, Receptors en LigandReceptor genen in Endometrium-dataset
liganden_endom <- intersect(goede_genen_in_endom_data, liganden)
receptors_endom <- intersect(goede_genen_in_endom_data, receptors)
ligreceptors_endom <- intersect(goede_genen_in_endom_data, ligreceptors)

df_endom_liganden <- data_endom2[liganden_endom, ]
df_endom_receptors <- data_endom2[receptors_endom, ]
df_endom_ligreceptors <- data_endom2[ligreceptors_endom, ]




