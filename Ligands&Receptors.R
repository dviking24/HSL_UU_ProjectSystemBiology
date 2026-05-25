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
data_endom2 <- data_endom[goede_genen_in_endom_data, ]


# Dataframes van Liganden, Receptors en LigandReceptor genen in Vivo-dataset
liganden_vivo <- intersect(goede_genen_in_vivo_data, liganden)
receptors_vivo <- intersect(goede_genen_in_vivo_data, receptors)
ligreceptors_vivo <- intersect(goede_genen_in_vivo_data, ligreceptors)

df_blastovivo_liganden <- data_blastovivo2[liganden_vivo, ]
df_blastovivo_receptors <- data_blastovivo2[receptors_vivo, ]
df_blastovivo_ligreceptors <- data_blastovivo2[ligreceptors_vivo, ]


# Dataframes van Liganden, Receptors en LigandReceptor genen in Endometrium-dataset
liganden_endom <- intersect(goede_genen_in_endom_data, liganden)
receptors_endom <- intersect(goede_genen_in_endom_data, receptors)
ligreceptors_endom <- intersect(goede_genen_in_endom_data, ligreceptors)

df_endom_liganden <- data_endom2[liganden_endom, ]
df_endom_receptors <- data_endom2[receptors_endom, ]
df_endom_ligreceptors <- data_endom2[ligreceptors_endom, ]


