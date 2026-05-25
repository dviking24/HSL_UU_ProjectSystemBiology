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
liganden <- data_ligrecep$ligand_ensembl
receptors <- data_ligrecep$receptor_ensembl
liganden_plus_receptors <- c(liganden, receptors)
genen_waar_we_naar_zoeken <- unique(liganden_plus_receptors)
alle_genen_in_vivo <- row.names(data_blastovivo)
alle_genen_in_endom <- row.names(data_endom)


# Filteren van genen in blastovivo-dataset
onnodige_genen_in_vivo_data <- setdiff(alle_genen_in_vivo, genen_waar_we_naar_zoeken)
goede_genen_in_vivo_data <- setdiff(alle_genen_in_vivo, onnodige_genen_in_vivo_data)

# Blastovivo-dataframe filteren zodat het alleen nog maar rijen met goede genen overhoudt
data_blastovivo2 <- data_blastovivo[goede_genen_in_vivo_data, ]


# Filteren van genen in endom-dataset
onnodige_genen_in_endom_data <- setdiff(alle_genen_in_endom, genen_waar_we_naar_zoeken)
goede_genen_in_endom_data <- setdiff(alle_genen_in_endom, onnodige_genen_in_endom_data)

# Endom-dataframe filteren zodat het alleen nog maar rijen met goede genen overhoudt
data_endom2 <- data_blastovivo[goede_genen_in_endom_data, ]






