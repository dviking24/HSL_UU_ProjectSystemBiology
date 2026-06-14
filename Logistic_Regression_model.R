library(glmnet)
library(pROC)
library(ggplot2)

#Inladen van de matrixen
vivo_matrix <- as.matrix(read.table("matrices/vivo_matrix.tsv", header = TRUE,
                                    row.names = 1, sep = "\t", check.names = FALSE))
vitro_matrix <- as.matrix(read.table("matrices/vitro_matrix.tsv", header = TRUE,
                                     row.names = 1, sep = "\t", check.names = FALSE))

#Elke iteraties moet een andere seed gebruiken
seeds <- c(64,28,21,94,41,12,53,22,17,62)


#Model functie
nogroup <- function(mat){

  #De matrix transponeren, dus rijen worden kolommen en andersom.
  X <- t(mat) 
  
  #Labels maken, dus NP en PR
  y <- factor(ifelse(grepl("NP", rownames(X)), "NP", "PR"), levels = c("NP","PR"))
  
  #Vectoren om resulaat opteslaan
  auc_vec <- c(); acc_vec <- c()
  coef_mat <- NULL
  
  
  #Loop voor de 10 verschillende iteraties
  for(i in seq_along(seeds)){
    set.seed(seeds[i])
    #Train test split 70%train en 30% test random selectie
    n <- nrow(X)
    train_idx <- sample(1:n, size = floor(0.7 * n))
    
    #Train en test data ready
    X_train <- X[train_idx, , drop=FALSE]; y_train <- y[train_idx]
    X_test  <- X[-train_idx, , drop=FALSE]; y_test  <- y[-train_idx]
    
    #Skip als er maar 1 klasse in de train of test data zit
    if(length(unique(y_train)) < 2 || length(unique(y_test)) < 2){next}
    
    #Logistic regression en cross validatie (LASSO)
    model <- cv.glmnet(X_train, y_train, family="binomial", alpha=1, nfolds=5)
    #Probability dat sample PR is, value tussen 0-1
    prob  <- predict(model, X_test, s="lambda.min", type="response") 
    
    #Classificatie, als de kans groter is dan 0.5 is het PR
    pred_class <- factor(ifelse(prob > 0.5, "PR", "NP"), levels=levels(y_test))
    acc_vec[i] <- mean(pred_class == y_test)
    
    #Vergelijkt voorspelling met echte labels en maakt ROC curve
    roc_obj    <- roc(y_test, as.numeric(prob), levels=c("NP","PR"), direction="<", quiet=TRUE)
    auc_vec[i] <- as.numeric(auc(roc_obj))
    
    cat("iter",i,"seed",seeds[i],"acc:",round(acc_vec[i],3),"AUC:",round(auc_vec[i],3),"\n")
  }
  cat("Mean ACC:", mean(acc_vec, na.rm=TRUE), " Mean AUC:", mean(auc_vec, na.rm=TRUE), "\n")
  return(list(acc=acc_vec, auc=auc_vec))
}

vivo_nogroup  <- nogroup(vivo_matrix)
vitro_nogroup <- nogroup(vitro_matrix)



#Model functie die groepeert op embryo en endometrium. Ttest bevat alleen samples
#waarvan zowel embryo als endometrium nieuw zijn (niet in train) om leakage te voorkomen
group_both <- function(mat, n_pr_test, n_np_test){
  #De matrix transponeren, dus rijen worden kolommen en andersom.
  X <- t(mat)
  
  #Labels maken, dus NP en PR
  y <- factor(ifelse(grepl("NP", rownames(X)), "NP", "PR"), levels = c("NP","PR"))
  
  #Embryo en endometrium per sample uit de naam halen
  parts  <- strsplit(rownames(X), "-")
  embryo <- sapply(parts, function(p) p[grepl("Zo|SW|IVT", p)][1])
  endom  <- sapply(parts, function(p) p[grepl("Hols|Cont|Jap|Sim", p)][1])
  ends <- unique(endom)
  emb_label <- tapply(as.character(y), embryo, function(x) unique(x)[1])
  pr_embs <- names(emb_label)[emb_label == "PR"]
  np_embs <- names(emb_label)[emb_label == "NP"]
  
  #Vectoren om resulaat opteslaan
  auc_vec <- c()
  acc_vec <- c() 
  coef_mat <- NULL
  
  
  #Loop voor de 10 verschillende iteraties
  for(i in seq_along(seeds)){
    set.seed(seeds[i])
    #70% van de embryo's en endometria voor train selecteren
    test_emb  <- c(sample(pr_embs, size = n_pr_test), sample(np_embs, size = n_np_test))
    train_emb <- setdiff(c(pr_embs, np_embs), test_emb)
    train_end <- sample(ends, size = floor(0.7 * length(ends)))
    
    #Train = embryo en endometrium bekend ; test = allebei nieuw
    train_idx <- which( embryo %in% train_emb  &  endom %in% train_end)
    test_idx  <- which(embryo %in% test_emb & !endom %in% train_end)
    
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
    
    cat("iter",i,"seed",seeds[i],"acc:",round(acc_vec[i],3),"AUC:",round(auc_vec[i],3),"\n")
  }
  cat("Mean ACC:", mean(acc_vec, na.rm=TRUE), " Mean AUC:", mean(auc_vec, na.rm=TRUE), "\n")
  
  freq <- colSums(coef_mat != 0, na.rm=TRUE) #telt hoe vaak non zero
  top10 <- names(sort(freq, decreasing=TRUE))[1:10]
  cat("Top 10 features:\n"); print(top10)

  return(list(acc=acc_vec, auc=auc_vec, freq=freq, top10=top10))
}
vivo_both  <- group_both(vivo_matrix,  n_pr_test = 4, n_np_test = 4)
vitro_both <- group_both(vitro_matrix, n_pr_test = 3, n_np_test = 4)





#Frequentie plot over hoe vaak de feautres zijn geselecteerd over de 10 iteraties

#vivo
vivo_freq_df <- data.frame(
  feature = names(sort(vivo_both$freq, decreasing=TRUE)[1:10]),
  freq = sort(vivo_both$freq, decreasing=TRUE)[1:10]
)

ggplot(vivo_freq_df, aes(x=reorder(feature, freq), y=freq)) +
  geom_bar(stat="identity", fill="steelblue") +
  coord_flip() +
  labs(title="Vivo top 10 features", x="Ligand-receptor paar", y="Frequentie (van 10 iteraties)") +
  theme_minimal()

#vitro
vitro_freq_df <- data.frame(
  feature = names(sort(vitro_both$freq, decreasing=TRUE)[1:10]),
  freq = sort(vitro_both$freq, decreasing=TRUE)[1:10]
)

ggplot(vitro_freq_df, aes(x=reorder(feature, freq), y=freq)) +
  geom_bar(stat="identity", fill="darkorange") +
  coord_flip() +
  labs(title="Vitro top 10 features", x="Ligand-receptor paar", y="Frequentie (van 10 iteraties)") +
  theme_minimal()