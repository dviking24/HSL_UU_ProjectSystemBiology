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






