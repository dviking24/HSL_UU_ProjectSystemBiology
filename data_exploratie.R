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



# --- data exploration --- 
# samplecolumns in endom:
ncol(data_endom)

# samplecolumns in blastovivo:
ncol(data_blastovivo)

# samplecolumns in blastovitro:
ncol(data_blastovitro)


# NA's in endom:
sum(is.na(data_endom))

# NA's in blastovivo:
sum(is.na(data_blastovivo))

# NA's in blastovitro:
sum(is.na(data_blastovitro))



# --- normality check --- 
# endom normality
shapiro.test(data_endom$nonPR_Hols_1[0:5000])
shapiro.test(data_endom$nonPR_Hols_1[5001:9840])
#hist(data_endom$nonPR_Hols_1)
#qqnorm(data_endom$nonPR_Hols_1)

endom_shapiro <- apply(data_endom[0:5000,], 2, function(x) {
  
  if(sd(x, na.rm = TRUE) == 0){
    return(NA)
  }
  
  shapiro.test(x)$p.value
})

# totaal
length(endom_shapiro)

# niet_normaal
sum(endom_shapiro < 0.05, na.rm = TRUE)

# percentage
mean(endom_shapiro < 0.05, na.rm = TRUE) * 100



# blastovivo normality
shapiro.test(data_blastovivo$Zo_NP_1[0:5000])
shapiro.test(data_blastovivo$Zo_NP_1[5001:10000])
#hist(data_blastovivo$Zo_NP_1)
#qqnorm(data_blastovivo$Zo_NP_1)

blastovivo_shapiro <- apply(data_blastovivo[0:5000,], 2, function(x) {
  
  if(sd(x, na.rm = TRUE) == 0){
    return(NA)
  }
  
  shapiro.test(x)$p.value
})

# totaal
length(blastovivo_shapiro)

# niet_normaal
sum(blastovivo_shapiro < 0.05, na.rm = TRUE)

# percentage
mean(blastovivo_shapiro < 0.05, na.rm = TRUE) * 100



# blastovitro normality
shapiro.test(data_blastovitro$Zo_IVT_NP_3[0:5000])
shapiro.test(data_blastovitro$Zo_IVT_NP_3[5001:9840])
#hist(data_blastovitro$Zo_IVT_NP_3)
#qqnorm(data_blastovitro$Zo_IVT_NP_3)

vitro_shapiro <- apply(data_blastovitro[0:5000,], 2, function(x) {
  
  if(sd(x, na.rm = TRUE) == 0){
    return(NA)
  }
  
  shapiro.test(x)$p.value
})

# totaal
length(vitro_shapiro)

# niet_normaal
sum(vitro_shapiro < 0.05, na.rm = TRUE)

# percentage
mean(vitro_shapiro < 0.05, na.rm = TRUE) * 100



# --- combining blasto datasets --- 
data_blasto <- cbind(data_blastovivo, data_blastovitro)

identical(rownames(data_blastovivo), rownames(data_blastovitro))


# pca 
pca_plot <- function(data, conditie, titel) {
  mat <- t(data)
  mat <- mat[, apply(mat, 2, function(x) sd(x) != 0)]
  pca <- prcomp(mat, center = TRUE, scale. = TRUE)
  var <- round(100 * summary(pca)$importance[2, 1:2], 1)
  scores <- as.data.frame(pca$x[, 1:2])
  scores$conditie <- conditie
  ggplot(scores, aes(PC1, PC2, color = conditie)) +
    geom_point(size = 3) +
    labs(title = titel, x = paste0("PC1 (", var[1], "%)"), y = paste0("PC2 (", var[2], "%)")) +
    theme_minimal()
}

pca_plot(data_endom, ifelse(grepl("^PR", colnames(data_endom)), "PR", "nonPR"), "PCA endom")
pca_plot(data_blastovivo, ifelse(grepl("_PR_", colnames(data_blastovivo)), "PR", "NP"), "PCA blasto in vivo")
pca_plot(data_blastovitro, ifelse(grepl("_PR_", colnames(data_blastovitro)), "PR", "NP"), "PCA blasto in vitro")

#boxplot
boxplot(data_endom, las = 2, cex.axis = 0.5, main = "Endom")
boxplot(data_blastovivo, las = 2, cex.axis = 0.5, main = "Blasto in vivo")
boxplot(data_blastovitro, las = 2, cex.axis = 0.5, main = "Blasto in vitro")