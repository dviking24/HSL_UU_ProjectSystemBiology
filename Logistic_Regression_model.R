library(glmnet)
library(pROC)
library(ggplot2)



#Inladen van de matrixen
vivo_matrix  <- read.table("matrices/vivo_matrix_emb2endo.tsv",  header = TRUE,
                           row.names = 1, sep = "\t", check.names = FALSE)
vitro_matrix <- read.table("matrices/vitro_matrix_emb2endo.tsv", header = TRUE,
                           row.names = 1, sep = "\t", check.names = FALSE)

#Elke iteraties moet een andere seed gebruiken
seeds <- c(64,28,21,94,41,12,53,22,17,62)

#Model functie die groepeert op embryo en endometrium. Test bevat alleen samples
#waarvan zowel embryo als endometrium in de test-groep zitten (gebalanceerd per
#klasse) om leakage te voorkomen.
group_both <- function(mat, train_n, test_n){
  #Samples staan in de rijen; PregnancyStatus-kolom eruit halen, rest = features
  feat_cols <- colnames(mat) != "PregnancyStatus"
  X <- as.matrix(mat[, feat_cols, drop=FALSE])
  
  #Labels maken, dus NP en PR (uit de samplenaam = rijnaam)
  y <- factor(ifelse(grepl("NP", rownames(X)), "NP", "PR"), levels = c("NP","PR"))
  
  #Embryo en endometrium per sample uit de naam halen (scheidingsteken is | )
  parts  <- strsplit(rownames(X), "\\|")
  embryo <- sapply(parts, function(p) p[grepl("Zo|SW|IVT", p)][1])
  endom  <- sapply(parts, function(p) p[grepl("Hols|Cont|Jap|Sim", p)][1])
  
  #Embryo's per klasse
  emb_label <- tapply(as.character(y), embryo, function(x) unique(x)[1])
  pr_embs <- names(emb_label)[emb_label == "PR"]
  np_embs <- names(emb_label)[emb_label == "NP"]
  
  #Endometria per klasse (label uit de endometrium-naam: NP_ of PR_)
  end_label <- ifelse(grepl("^NP", unique(endom)), "NP", "PR")
  names(end_label) <- unique(endom)
  pr_ends <- names(end_label)[end_label == "PR"]
  np_ends <- names(end_label)[end_label == "NP"]
  
  #Vectoren om resulaat opteslaan
  auc_vec <- c()
  acc_vec <- c()
  coef_mat <- NULL
  shap_per_seed <- NULL
  
  #Loop voor de 10 verschillende iteraties
  for(i in seq_along(seeds)){
    set.seed(seeds[i])
    
    #Embryo's gebalanceerd per klasse in train en test
    train_emb <- c(sample(pr_embs, size = train_n), sample(np_embs, size = train_n))
    test_emb  <- c(sample(setdiff(pr_embs, train_emb), size = test_n),
                   sample(setdiff(np_embs, train_emb), size = test_n))
    
    #Endometria gebalanceerd per klasse in train en test (zelfde logica als embryo)
    train_end <- c(sample(pr_ends, size = train_n), sample(np_ends, size = train_n))
    test_end  <- c(sample(setdiff(pr_ends, train_end), size = test_n),
                   sample(setdiff(np_ends, train_end), size = test_n))
    
    #Train = embryo EN endometrium in train-groep ; test = embryo EN endometrium in test-groep
    train_idx <- which(embryo %in% train_emb & endom %in% train_end)
    test_idx  <- which(embryo %in% test_emb  & endom %in% test_end)
    
    #Train en test data ready
    X_train <- X[train_idx, , drop=FALSE]; y_train <- y[train_idx]
    X_test  <- X[test_idx,  , drop=FALSE]; y_test  <- y[test_idx]
    
    #Aantal samples in train en test laten zien
    cat("iter",i,"seed",seeds[i],"| train n:",length(train_idx),"| test n:",length(test_idx),"\n")
    
    #Skip als train of test maar 1 klasse heeft
    if(length(unique(y_train)) < 2 || length(unique(y_test)) < 2){next}
    
    #Logistic regression en cross validatie (LASSO)
    model <- cv.glmnet(X_train, y_train, family="binomial", alpha=1, nfolds=5)
    
    # Coëfficiënten opslaan
    coef_full <- coef(model, s="lambda.min")
    coef_i <- as.numeric(coef_full)[-1] # haalt intercept weg
    feat_names <- rownames(coef_full)[-1] # namen features
    if(is.null(coef_mat)) coef_mat <- matrix(NA, nrow=length(seeds), ncol=length(coef_i),
                                             dimnames=list(NULL, feat_names))
    coef_mat[i, ] <- coef_i
    
    #Probability dat sample PR is, value tussen 0-1
    prob  <- predict(model, X_test, s="lambda.min", type="response")
    
    #Classificatie, als de kans groter is dan 0.5 is het PR
    pred_class <- factor(ifelse(prob > 0.5, "PR", "NP"), levels=levels(y_test))
    acc_vec[i] <- mean(pred_class == y_test)
    
    #Vergelijkt voorspelling met echte labels en maakt ROC curve
    roc_obj    <- roc(y_test, as.numeric(prob), levels=c("NP","PR"), direction="<", quiet=TRUE)
    auc_vec[i] <- as.numeric(auc(roc_obj))
    
    #SHAP per iteratie: analytisch want LASSO lineair is (beta * (x - train gemiddelde))
    betas <- coef_i; names(betas) <- feat_names
    X_test_centered <- sweep(X_test, 2, colMeans(X_train), "-")  #centreren op train gemiddelde
    shap_test <- sweep(X_test_centered, 2, betas, "*")           #SHAP per sample per feature
    shap_imp  <- colMeans(abs(shap_test))                        #mean SHAP per feature deze seed
    if(is.null(shap_per_seed)) shap_per_seed <- matrix(NA, nrow=length(seeds), ncol=length(shap_imp),
                                                       dimnames=list(NULL, names(shap_imp)))
    shap_per_seed[i, ] <- shap_imp
    
    cat("iter",i,"seed",seeds[i],"acc:",round(acc_vec[i],3),"AUC:",round(auc_vec[i],3),"\n")
  }
  cat("Mean ACC:", mean(acc_vec, na.rm=TRUE), " Mean AUC:", mean(auc_vec, na.rm=TRUE), "\n")
  
  freq <- colSums(coef_mat != 0, na.rm=TRUE) #telt hoe vaak non zero
  top10 <- names(sort(freq, decreasing=TRUE))[1:10]
  cat("Top 10 features:\n"); print(top10)
  
  #Gemiddelde mean|SHAP| over de 10 iteraties
  shap_importance <- colMeans(shap_per_seed, na.rm=TRUE)
  
  return(list(acc=acc_vec, auc=auc_vec, freq=freq, top10=top10,
              shap_importance=shap_importance))
}
vivo_both  <- group_both(vivo_matrix,  train_n = 7, test_n = 4)
vitro_both <- group_both(vitro_matrix, train_n = 5, test_n = 4)



#Frequentie plot over hoe vaak de feautres zijn geselecteerd over de 10 iteraties

#vivo
vivo_freq_df <- data.frame(
  feature = names(sort(vivo_both$freq, decreasing=TRUE)[1:10]),
  freq = sort(vivo_both$freq, decreasing=TRUE)[1:10]
)

ggplot(vivo_freq_df, aes(x=reorder(feature, freq), y=freq)) +
  geom_bar(stat="identity", fill="steelblue") +
  coord_flip() +
  labs(title="Vivo top 10 features (emb2endo)", x="Ligand-receptor paar", y="Frequentie (van 10 iteraties)") +
  theme_minimal()

#vitro
vitro_freq_df <- data.frame(
  feature = names(sort(vitro_both$freq, decreasing=TRUE)[1:10]),
  freq = sort(vitro_both$freq, decreasing=TRUE)[1:10]
)

ggplot(vitro_freq_df, aes(x=reorder(feature, freq), y=freq)) +
  geom_bar(stat="identity", fill="darkorange") +
  coord_flip() +
  labs(title="Vitro top 10 features (emb2endo)", x="Ligand-receptor paar", y="Frequentie (van 10 iteraties)") +
  theme_minimal()



#SHAP importance plot van de top 10 features (gemiddelde SHAP over 10 iteraties)

#vivo
vivo_shap_df <- data.frame(
  feature = names(sort(vivo_both$shap_importance, decreasing=TRUE)[1:10]),
  importance = sort(vivo_both$shap_importance, decreasing=TRUE)[1:10]
)

ggplot(vivo_shap_df, aes(x=reorder(feature, importance), y=importance)) +
  geom_bar(stat="identity", fill="steelblue") +
  coord_flip() +
  labs(title="Vivo SHAP importance top 10 (emb2endo)", x="Ligand-receptor paar", y="Gemiddelde SHAP") +
  theme_minimal()

#vitro
vitro_shap_df <- data.frame(
  feature = names(sort(vitro_both$shap_importance, decreasing=TRUE)[1:10]),
  importance = sort(vitro_both$shap_importance, decreasing=TRUE)[1:10]
)

ggplot(vitro_shap_df, aes(x=reorder(feature, importance), y=importance)) +
  geom_bar(stat="identity", fill="darkorange") +
  coord_flip() +
  labs(title="Vitro SHAP importance top 10 (emb2endo)", x="Ligand-receptor paar", y="Gemiddelde SHAP") +
  theme_minimal()



#Top features exporteren naar CSV 
N_TOP <- 10  #aantal top features per dataset

vivo_top  <- names(sort(vivo_both$shap_importance,  decreasing=TRUE))[1:N_TOP]
vitro_top <- names(sort(vitro_both$shap_importance, decreasing=TRUE))[1:N_TOP]

shap_features <- rbind(
  data.frame(feature = vivo_top,  dataset = "vivo"),
  data.frame(feature = vitro_top, dataset = "vitro")
)

#Map aanmaken als die nog niet bestaat en wegschrijven
if(!dir.exists("feature_list_LR")) dir.create("feature_list_LR")
write.csv(shap_features, "feature_list_LR/top_10_shap_features_emb2endo.csv", row.names = FALSE)



#Matrix filteren op alleen de top SHAP-features (per dataset)
shap_feats <- read.csv("feature_list_LR/top_10_shap_features_emb2endo.csv", stringsAsFactors = FALSE)

vivo_feats  <- shap_feats$feature[shap_feats$dataset == "vivo"]
vitro_feats <- shap_feats$feature[shap_feats$dataset == "vitro"]

#De matrices hebben features als rijen
vivo_matrix_shap  <- vivo_matrix[,  c("PregnancyStatus", intersect(vivo_feats,  colnames(vivo_matrix))),  drop=FALSE]
vitro_matrix_shap <- vitro_matrix[, c("PregnancyStatus", intersect(vitro_feats, colnames(vitro_matrix))), drop=FALSE]

#Model opnieuw draaien op alleen de SHAP-features
vivo_shap_only  <- group_both(vivo_matrix_shap,  train_n = 7, test_n = 4)
vitro_shap_only <- group_both(vitro_matrix_shap, train_n = 5, test_n = 4)