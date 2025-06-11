
##### Tany-seq Data Collection and Labelling #####

#Sullivan, A.I.; Potthoff, M.J.; Flippo, K.H. Tany-Seq: Integrated Analysis of the Mouse Tanycyte Transcriptome. Cells 2022, 11, 1565. https://doi.org/10.3390/cells11091565
#Data obtained from Flippo, Kyle (2022), “Tany-Seq: Integrated Analysis of the Mouse Tanycyte Transcriptome”, Mendeley Data, V1, doi: 10.17632/w8yw2c92jg.1
#File name: IntegratedTanycyteDataset.rds

library(plyr)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(docopt)
library(tibble)
library(tidyr)
library(stringr)
library(Matrix)
library(Seurat)
library(SeuratDisk)
library(future)
library(anndata)
library(cowplot)
library(patchwork)

dir.create("/data/pharma_macrophages/Umit/20250602/", showWarnings = TRUE, recursive = FALSE, mode = "0777")
setwd("/data/pharma_macrophages/Umit/20250602/")

l <-16*1024^3
options(future.globals.maxSize = l)

tanyseq <- readRDS("/data/pharma_macrophages/IntegratedTanycyteDataset.rds")

#Labels were done according to the original paper suggested.

tanyseq$label <- as.character(tanyseq$seurat_clusters)
tanyseq$label[tanyseq$seurat_clusters == "5"] <- "alpha 1"
tanyseq$label[tanyseq$seurat_clusters %in% c("1","4")] <- "alpha 2"
tanyseq$label[tanyseq$seurat_clusters %in% c("0","3")] <- "beta 1"
tanyseq$label[tanyseq$seurat_clusters == "2"] <- "beta 2"

head(tanyseq$label)

#To save labeled tanycytes from Tany-seq
saveRDS(tanyseq, "/data/pharma_macrophages/IntegratedTanycyteDataset.rds")



##### Hypomap Data Collection and Labelling #####

#Steuernagel, L., Lam, B.Y.H., Klemm, P. et al. HypoMap—a unified single-cell gene expression atlas of the murine hypothalamus. Nat Metab 4, 1402–1419 (2022). https://doi.org/10.1038/s42255-022-00657-y
#Data obtained from Lam, Y. H., Yeo, G., Steuernagel, L., Klemm, P., & Brüning, J. (2022). Research data supporting “HypoMap – a unified single cell gene expression atlas of the murine hypothalamus”. Apollo - University of Cambridge Repository. https://doi.org/10.17863/CAM.87955
#File name: hypoMap.rds

hypomap <- readRDS("/data/pharma_macrophages/hypoMap.rds")

#Subsetting the cells that were labelled as Tanycyte
# Define the function
is_tanycyte <- function(x) {
  str_detect(tolower(x), "tany") #tany as keyword
}

#Check which candidate columns exist in the metadata, Original research, C66 label and hypomap labels were checked
metadata_cols <- colnames(hypomap@meta.data)
cols_to_check <- c("Author_Class_Curated", "C66_named", "Author_CellType")
existing_cols <- cols_to_check[cols_to_check %in% metadata_cols]

#Now create a logical vector that is TRUE if any of the existing columns contain "tanycyte"
tany_logical <- Reduce("|", lapply(existing_cols, function(col) {
  is_tanycyte(hypomap@meta.data[[col]])
}))

#Subset the object using the logical vector
tany_hypomap <- subset(hypomap, cells = rownames(hypomap@meta.data)[tany_logical])

#To save tanycyte data from hypomap
saveRDS(tany_hypomap, "/data/pharma_macrophages/tany_hypomap.rds")