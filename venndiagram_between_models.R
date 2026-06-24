# =========================================================================
# IMPORTS
# =========================================================================
# install.packages("VennDiagram")
# install.packages("dplyr")
# install.packages("scales")
library(VennDiagram)
library(dplyr)
library(scales)

# =========================================================================
# 1. SETTINGS & LOAD DATA
# =========================================================================
file_model1 <- "feature_lists/top_10_shap_features_endo2emb.csv"
file_model2 <- "feature_lists/top_10_featImp_features_endo2emb.csv"

name_model1 <- "RANDOM FOREST"
name_model2 <- "RF - featIMP"

direction <- "endo2emb" 

data_model1 <- read.csv(file_model1, stringsAsFactors = FALSE)
data_model2 <- read.csv(file_model2, stringsAsFactors = FALSE)

# =========================================================================
# 2. ANALYSIS & PLOTTING FUNCTION
# =========================================================================
analyze_between_models <- function(target_dataset) {
  
  # Filter de data om alleen "vivo" of "vitro" op te halen per model
  features_m1 <- data_model1 %>% filter(dataset == target_dataset) %>% pull(feature)
  features_m2 <- data_model2 %>% filter(dataset == target_dataset) %>% pull(feature)
  
  # --- Zoek naar Overeenkomsten en Verschillen ---
  overlapping_features <- intersect(features_m1, features_m2)
  unique_to_m1         <- setdiff(features_m1, features_m2)
  unique_to_m2         <- setdiff(features_m2, features_m1)
  
  # Print de resultaten overzichtelijk naar de console
  message("\n=======================================================")
  message("OVERLAP ANALYSIS FOR: ", toupper(target_dataset), " (", direction, ")")
  message("COMPARING: ", name_model1, " vs ", name_model2)
  message("=======================================================")
  
  message("Shared Features (", length(overlapping_features), "):")
  if(length(overlapping_features) > 0) print(overlapping_features) else message("  None")
  
  message("\nUnique to ", name_model1, " (", length(unique_to_m1), "):")
  if(length(unique_to_m1) > 0) print(unique_to_m1) else message("  None")
  
  message("\nUnique to ", name_model2, " (", length(unique_to_m2), "):")
  if(length(unique_to_m2) > 0) print(unique_to_m2) else message("  None")
  
  img_filename <- paste0("venn_overlap_models_", direction, "_", target_dataset, ".png")
  
  venn_list <- list()
  venn_list[[name_model1]] <- features_m1
  venn_list[[name_model2]] <- features_m2
  
  venn.diagram(
    x = venn_list,
    filename = img_filename,
    output = TRUE,
    disable.logging = TRUE, 
    
    # --- Styling Options ---
    main = paste("Top 10 SHAP Overlap:", name_model1, "vs", name_model2),
    sub = paste("Dataset:", toupper(target_dataset), "| Direction:", direction),
    main.cex = 1.3,
    main.fontface = "bold",
    sub.cex = 1.1,
    col = c("#440154ff", "#fdae61ff"),
    fill = c(alpha("#440154ff", 0.4), alpha("#fdae61ff", 0.4)),
    cex = 2,           
    cat.cex = 1.2,     
    cat.pos = c(-20, 20) 
  )
  
  message("\n--> Venn diagram saved as: ", img_filename)
}

# =========================================================================
# 3. RUN THE ANALYSIS
# =========================================================================
analyze_between_models("vivo")
analyze_between_models("vitro")