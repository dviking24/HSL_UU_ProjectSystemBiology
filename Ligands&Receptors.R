library(preprocessCore)
library(ggplot2)


# Cross-data normalisatie op alle genen eerst

# Vivo + Endom
gemeenschappelijk_vivo <- intersect(rownames(data_blastovivo), rownames(data_endom))
vivo_sub  <- data_blastovivo[gemeenschappelijk_vivo, ]
endom_sub <- data_endom[gemeenschappelijk_vivo, ]
gecombineerd_vivo <- as.matrix(cbind(vivo_sub, endom_sub))

genorm <- normalize.quantiles(gecombineerd_vivo)
dimnames(genorm) <- dimnames(gecombineerd_vivo)

n_vivo <- ncol(vivo_sub)
vivo_norm  <- genorm[, 1:n_vivo]
endom_norm <- genorm[, (n_vivo + 1):ncol(genorm)]

# Vitro + Endom
gemeenschappelijk_vitro <- intersect(rownames(data_blastovitro), rownames(data_endom))
vitro_sub  <- data_blastovitro[gemeenschappelijk_vitro, ]
endom_sub_v <- data_endom[gemeenschappelijk_vitro, ]
gecombineerd_vitro <- as.matrix(cbind(vitro_sub, endom_sub_v))

genorm_v <- normalize.quantiles(gecombineerd_vitro)
dimnames(genorm_v) <- dimnames(gecombineerd_vitro)

n_vitro <- ncol(vitro_sub)
vitro_norm  <- genorm_v[, 1:n_vitro]
endom_norm_v <- genorm_v[, (n_vitro + 1):ncol(genorm_v)]


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

groep_vivo  <- c(rep("vivo", n_vivo), rep("endom", ncol(endom_sub)))
groep_vitro <- c(rep("vitro", n_vitro), rep("endom", ncol(endom_sub_v)))

pca_check(gecombineerd_vivo,  groep_vivo,  "Vivo - voor normalisatie")
pca_check(genorm,             groep_vivo,  "Vivo - na normalisatie")
pca_check(gecombineerd_vitro, groep_vitro, "Vitro - voor normalisatie")
pca_check(genorm_v,           groep_vitro, "Vitro - na normalisatie")



# Filter de genormaliseerde matrices op alleen relevante genen
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