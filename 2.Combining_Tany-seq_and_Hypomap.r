##### Combining Tany-seq and Hypomap Data #####

#To observe the sources of the original data in both dataset
table(tanyseq$Origin)
table(tany_hypomap$Dataset)

#In hypomap we have already have Campbell and Chen databases. These are removed from Tany-seq.
tanyseq <- subset(x=tanyseq, subset = Origin %in% c("Blackshaw","Deng"))

#Now I will add the following columns for tanyseq data: Age, Diet and Inferred Sex columns will be kept for later.

tanyseq$Diet <- tanyseq$State #Diet is chow for all cells
tanyseq$Age <- tanyseq$orig.ident #Age
tanyseq$inferred_sex <- tanyseq$Origin #Inferred Sex

#Age
tanyseq$Age[tanyseq$orig.ident == "D8C1"] <- "0-3 weeks"
tanyseq$Age[tanyseq$orig.ident == "D17C1"] <- "0-3 weeks"
tanyseq$Age[tanyseq$orig.ident == "D45C1"] <- "6+ weeks"
tanyseq$Age[tanyseq$orig.ident == "TanycyteAnalysis"] <- "6+ weeks" 

#Inferred sex
tanyseq$inferred_sex[tanyseq$Origin == "Blackshaw"] <- "U" #This was not clear, other experiments were done on male mice but this has no clear answer.
tanyseq$inferred_sex[tanyseq$Origin == "Deng"] <- "M" #All mice used in this experiments were male

#Same thing for hypomap,normal chow in Diet column was renamed as chow.
unique(tany_hypomap$Diet)
tany_hypomap$Diet[tany_hypomap$Diet == "Normal chow"] <- "Chow"

#There are some NA cells in hypomap and look for the NA age groups to find the ages?
na_cells <- rownames(tany_hypomap@meta.data[!(tany_hypomap@meta.data$Age %in% c("3-6 weeks", "6+ weeks", "0-3 weeks")), ])
tany_NA <- subset(tany_hypomap, cells = na_cells)
unique(tany_NA$Dataset) #These cells were coming from Affinati10x, Anderson10x.

#Affinati10x and Anderson10x mice all were over 6 weeks old.
tany_hypomap$Age[tany_hypomap$Dataset == "Affinati10x"] <- "6+ weeks"
tany_hypomap$Age[tany_hypomap$Dataset == "Anderson10x"] <- "6+ weeks"

#Now I will create Cell ID column for tanyseq, because Cell_IDs were embedded in rownames.
# Create a Cell_ID column in tanyseq using the rownames. 
tanyseq@meta.data$Cell_ID <- rownames(tanyseq@meta.data)

#To check the label mismatch I had to look common barcodes. I transform barcodes to bare minimum.
#Lets extract the barcodes.
extract_opt <- function(cell_id) {
  # Split by underscore into parts
  parts <- str_split(cell_id, "_", simplify = TRUE)
  
  # If no underscore was found, return the original string
  if(ncol(parts) == 1) {
    return(cell_id)
  } else {
    # Check if any part ends with "-1"
    for(part in parts) {
      if(str_ends(part, "-1")) {
        return(part)
      }
    }
    # If no part ends with "-1", choose the part with the maximum length
    lens <- nchar(parts)
    return(parts[which.max(lens)])
  }
}

#Removes _ and rakes the longest bits
tany_hypomap@meta.data$barcode <- sapply(tany_hypomap@meta.data$Cell_ID, extract_opt)
tanyseq@meta.data$barcode      <- sapply(tanyseq@meta.data$Cell_ID, extract_opt)

#Lets see the common barcodes
common_barcodes <- intersect(tany_hypomap@meta.data$barcode, tanyseq@meta.data$barcode)

#[1] "TATTCCACACACCGAC-1" there is one it is weird but it is okay to keep

#Barcode number is not matching to cell ids. These extra cell ids were not present in common_barcodes. So it is fine.
length(unique(tanyseq@meta.data$barcode))
length(unique(tanyseq@meta.data$Cell_ID))

length(unique(tany_hypomap@meta.data$barcode))
length(unique(tany_hypomap@meta.data$Cell_ID))

#There are some duplicate barcodes but it is not important for now.

#Lets merge them now.
#Before that add a column for the cells origin
tany_hypomap@meta.data$DataOrigin <- "Hypomap"
tanyseq@meta.data$DataOrigin <- "Tanyseq"

# Merge the two Seurat objects
combined <- merge(tany_hypomap, y = tanyseq)

# If there are duplicates in the Cell_ID column on the combined object:
dup_ids <- combined@meta.data %>%
  group_by(Cell_ID) %>%
  filter(n() > 1)

print(dup_ids)

#There is none so no need to do this below.
  
if(nrow(dup_ids) > 0){
  message("There are duplicates. Removing them...")
  # Remove duplicates based on the Cell_ID. Here we keep the first occurrence.
  combined@meta.data <- combined@meta.data %>% distinct(Cell_ID, .keep_all = TRUE)
  # Note: After this, you may need to subset the combined Seurat object accordingly.
  combined <- subset(combined, cells = rownames(combined@meta.data))
} else {
  message("No duplicates found.")
}
