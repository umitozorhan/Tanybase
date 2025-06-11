##### DE Genes Calculation for Individual Cells (No pseudobulk) #####
tany <-tany_integrated
tany$label.Age.Diet <- paste(tany$label,tany$Age,tany$Diet, sep="_")
unique(tany$label.Age.Diet)

# Set identity
Idents(tany) <- "label.Age.Diet"

# List of unique identity groups
groups <- unique(Idents(tany))

# Initialize list to collect all DE results
all_results <- list()

# Function to check if two strings share any component when split by "_"
has_shared_component <- function(name1, name2) {
  parts1 <- unlist(strsplit(as.character(name1), "_"))
  parts2 <- unlist(strsplit(as.character(name2), "_"))
  return(length(intersect(parts1, parts2)) >= 2)
}

# Loop through all pairwise combinations
for (i in 1:(length(groups) - 1)) {
  for (j in (i + 1):length(groups)) {
    
    group1 <- groups[i]
    group2 <- groups[j]
    
    # Only compare if there's at least one shared component
    if (!has_shared_component(group1, group2)) next
    
    # Perform DE analysis
    de_res <- tryCatch({
      FindMarkers(tany, ident.1 = group1, ident.2 = group2,
                  logfc.threshold = 0, min.pct = 0.1)
    }, error = function(e) NULL)
    
    # Skip if empty or failed
    if (is.null(de_res) || nrow(de_res) == 0) next
    
    # Filter for meaningful DE genes
    de_res_filtered <- de_res[
      de_res$p_val_adj <= 0.05 & abs(de_res$avg_log2FC) >= 0.3, 
    ]
    
    if (nrow(de_res_filtered) > 0) {
      # Add gene and group info
      de_res_filtered$gene <- rownames(de_res_filtered)
      de_res_filtered$group1 <- group1
      de_res_filtered$group2 <- group2
      
      # Add to results list
      all_results[[paste(group1, group2, sep = "_vs_")]] <- de_res_filtered
    }
  }
}

# Combine all DE results into a single data frame
final_df <- do.call(rbind, all_results)
final_df <- subset(final_df, pct.1 >= 0.1 & pct.2 >= 0.2)

has_2_shared_components <- function(name1, name2) {
  parts1 <- unlist(strsplit(as.character(name1), "_"))
  parts2 <- unlist(strsplit(as.character(name2), "_"))
  return(length(intersect(parts1, parts2)) >= 2)
}

final_df <- final_df[mapply(has_2_shared_components, final_df$group1, final_df$group2), ]

# Save to CSV
write.csv(final_df, file = "all_DE_filtered_shared_comparisons.csv", row.names = FALSE)
