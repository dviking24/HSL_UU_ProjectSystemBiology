# =========================================================================
# IMPORTS
# =========================================================================
# install.packages("VennDiagram")
# install.packages("dplyr")
library(VennDiagram)
library(dplyr)

# =========================================================================
# 1. LOAD DATA
# =========================================================================
# Read in the top 10 SHAP feature files
emb2endo_data <- read.csv("feature_lists/top_10_shap_features_emb2endo.csv", stringsAsFactors = FALSE)
endo2emb_data <- read.csv("feature_lists/top_10_shap_features_endo2emb.csv", stringsAsFactors = FALSE)

# =========================================================================
# 2. ANALYSIS & PLOTTING FUNCTION
# =========================================================================
analyze_and_plot_overlap <- function(target_dataset) {
  
  # Filter the data to isolate just "vivo" or just "vitro"
  set_emb2endo <- emb2endo_data %>% filter(dataset == target_dataset) %>% pull(feature)
  set_endo2emb <- endo2emb_data %>% filter(dataset == target_dataset) %>% pull(feature)
  
  # --- Find Overlaps and Differences ---
  overlapping_features <- intersect(set_emb2endo, set_endo2emb)
  unique_to_emb2endo   <- setdiff(set_emb2endo, set_endo2emb)
  unique_to_endo2emb   <- setdiff(set_endo2emb, set_emb2endo)
  
  # Print the actual gene names to the console
  message("\n=======================================================")
  message("OVERLAP ANALYSIS FOR: ", toupper(target_dataset))
  message("=======================================================")
  
  message("Shared Features (", length(overlapping_features), "):")
  if(length(overlapping_features) > 0) print(overlapping_features) else message("  None")
  
  message("\nUnique to emb2endo (", length(unique_to_emb2endo), "):")
  print(unique_to_emb2endo)
  
  message("\nUnique to endo2emb (", length(unique_to_endo2emb), "):")
  print(unique_to_endo2emb)
  
  # --- Create the Venn Diagram ---
  # Define the filename
  img_filename <- paste0("venn_overlap_", target_dataset, ".png")
  
  # Generate the plot (VennDiagram writes this directly to your working directory)
  venn.diagram(
    x = list(
      "emb2endo" = set_emb2endo,
      "endo2emb" = set_endo2emb
    ),
    filename = img_filename,
    output = TRUE,
    disable.logging = TRUE, # Keeps the console clean from log files
    
    # --- Styling Options ---
    main = paste("Top 10 SHAP Overlap -", toupper(target_dataset)),
    main.cex = 1.5,
    main.fontface = "bold",
    col = c("#440154ff", "#21908dff"), # Circle border colors
    fill = c(alpha("#440154ff", 0.4), alpha("#21908dff", 0.4)), # Transparent fill
    cex = 2,           # Size of the numbers inside the diagram
    cat.cex = 1.2,     # Size of the category labels
    cat.pos = c(-20, 20) # Slightly angles the labels away from the top center
  )
  
  message("\n--> Venn diagram saved as: ", img_filename)
}

# =========================================================================
# 3. RUN THE ANALYSIS
# =========================================================================
# Run the pipeline for both vivo and vitro
analyze_and_plot_overlap("vivo")
analyze_and_plot_overlap("vitro")
