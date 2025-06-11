##### Relabelling the unlabelled Tanycytes that carried over from Hypomap #####

#Creating module scores from the markers that were obtained from Fig 1E. of Tany-seq paper.

alpha1 <- list(c("Mafb","Necab2","Gadd45g","Stmn1","Mef2c","Ramp1","Pgnmc1","Cryab","Tgfb2"))
alpha2 <- list(c("Pdzph1","Fabp5","Cpe","S100a6","Slc7a11","Ptprz1","Itm2b","Rcn1"))
beta1 <- list(c("Frzb","Crym","Aldoc","Ndn","Ptn","Sparcl1","Sat1","Ly6h"))
beta2 <- list(c("Cntfr","Gpc3","Col27a1","Grlk3","Deptor","Ppp1r1b","Scn7a","Adm","Mest","Cldn10","Zfand6","Tmem176a","Sub1","Igfbp4"))

#Module scores to check tanycytic subtypes of individual cells.
tany_integrated <- AddModuleScore(
  tany_integrated,
  alpha1,
  nbin = 24,
  ctrl = 100,
  assay = "RNA",
  name = "alpha1",
  seed = 1,
  slot = "data",
)

tany_integrated <- AddModuleScore(
  tany_integrated,
  alpha2,
  nbin = 24,
  ctrl = 100,
  assay = "RNA",
  name = "alpha2",
  seed = 1,
  slot = "data",
)

tany_integrated <- AddModuleScore(
  tany_integrated,
  beta1,
  nbin = 24,
  ctrl = 100,
  assay = "RNA",
  name = "beta1",
  seed = 1,
  slot = "data",
)

tany_integrated <- AddModuleScore(
  tany_integrated,
  beta2,
  nbin = 24,
  ctrl = 100,
  assay = "RNA",
  name = "beta2",
  seed = 1,
  slot = "data",
)


head(tany_integrated) #To check if it worked or not

#Print the modules

p1 <- DimPlot(tany_integrated, group.by = "label", reduction = "umap", shuffle = TRUE, label = TRUE, raster=FALSE) + ggtitle("Subclusters") + theme(axis.title = element_text(size = 14)) + NoLegend()
p2 <- FeaturePlot(tany_integrated, "alpha11", reduction = "umap",  label = FALSE, raster=FALSE) + ggtitle("Alpha1 markers") + theme(axis.title = element_text(size = 14)) 
p3 <- FeaturePlot(tany_integrated, "alpha21", reduction = "umap",  label = FALSE, raster=FALSE) + ggtitle("Alpha2 markers") + theme(axis.title = element_text(size = 14)) 
p4 <- FeaturePlot(tany_integrated, "beta11", reduction = "umap",  label = FALSE, raster=FALSE) + ggtitle("Beta1 markers") + theme(axis.title = element_text(size = 14)) 
p5 <- FeaturePlot(tany_integrated, "beta21", reduction = "umap",  label = FALSE, raster=FALSE) + ggtitle("Beta2 markers") + theme(axis.title = element_text(size = 14)) 
ggsave(filename = "Tany_Scores.png", 
       height = 15, width = 10, dpi=300, 
       plot = plot_grid(p1, p2, p3, p4,p5, nrow = 3))

#Fig 1C. of Tany-seq paper.    
p1 <- DimPlot(tany_integrated, group.by = "label", reduction = "umap", shuffle = TRUE, label = TRUE, raster=FALSE) + ggtitle("Subclusters") + theme(axis.title = element_text(size = 14)) + NoLegend()
p2 <- DimPlot(tany_integrated, group.by = "seurat_clusters", reduction = "umap", shuffle = TRUE, label = TRUE, raster=FALSE) + ggtitle("Subclusters") + theme(axis.title = element_text(size = 14)) + NoLegend()       
p3 <- FeaturePlot(tany_integrated, "Mafb", reduction = "umap",  label = FALSE, raster=FALSE) + theme(axis.title = element_text(size = 14)) 
p4 <- FeaturePlot(tany_integrated, "Pdzph1", reduction = "umap",  label = FALSE, raster=FALSE) +  theme(axis.title = element_text(size = 14)) 
p5 <- FeaturePlot(tany_integrated, "Frzb", reduction = "umap",  label = FALSE, raster=FALSE) +  theme(axis.title = element_text(size = 14)) 
p6 <- FeaturePlot(tany_integrated, "Scn7a", reduction = "umap",  label = FALSE, raster=FALSE) + theme(axis.title = element_text(size = 14))  
ggsave(filename = "Tany_markers.png", 
       height = 15, width = 10, dpi=300, 
       plot = plot_grid(p1, p2, p3, p4,p5,p6, nrow = 3))      
       

#Fix the column names
tany_integrated$alpha1 <- tany_integrated$alpha11 
tany_integrated$alpha2 <- tany_integrated$alpha21 
tany_integrated$beta1 <- tany_integrated$beta11 
tany_integrated$beta2 <- tany_integrated$beta21 

tany_integrated$alpha11 <-NULL
tany_integrated$alpha21  <-NULL
tany_integrated$beta11  <-NULL
tany_integrated$beta21  <-NULL

#Save the integrated dataset
saveRDS(tany_integrated, "/data/pharma_macrophages/Tanycyte_integrated_080525.rds")



#Subset nonlabelled or incorrectly labelled cells from tany_integrated.
tany_integrated$keep <- ifelse(tany_integrated$label %in% c("alpha 1", "alpha 2", "beta 1", "beta 2"), FALSE, TRUE)
tanycyte1 <- subset(tany_integrated, subset = keep) #tanycyte1 only contains incorrectly labelled cells
tany_integrated$keep <- NULL

#Here is the cells that labelled incorrectly
table(tanycyte1$label)

#Astrocytes   Dividing   Endothelial   Ependymal   Fibroblast   Mural   Neurons   Oligodendrocytes   Tanycytes
#     9           1           1           6             3         2        7              1              40

#Modulescores for alpha1, alpha2, beta1 and beta2 were plotted for incorrectly labelled cells and clusters for each subtype of tanycytes were determined
p1 <- DimPlot(tanycyte1, group.by = "label", reduction = "umap", shuffle = TRUE, label = TRUE, raster=FALSE) + ggtitle("Subclusters") + theme(axis.title = element_text(size = 14)) + NoLegend()
p2 <- DimPlot(tanycyte1, group.by = "seurat_clusters", reduction = "umap", shuffle = TRUE, label = TRUE, raster=FALSE) + ggtitle("Subclusters") + theme(axis.title = element_text(size = 14)) + NoLegend()
p3 <- FeaturePlot(tanycyte1, "alpha1", reduction = "umap",  label = FALSE, raster=FALSE) + ggtitle("Alpha1 markers") + theme(axis.title = element_text(size = 14)) 
p4 <- FeaturePlot(tanycyte1, "alpha2", reduction = "umap",  label = FALSE, raster=FALSE) + ggtitle("Alpha2 markers") + theme(axis.title = element_text(size = 14)) 
p5 <- FeaturePlot(tanycyte1, "beta1", reduction = "umap",  label = FALSE, raster=FALSE) + ggtitle("Beta1 markers") + theme(axis.title = element_text(size = 14)) 
p6 <- FeaturePlot(tanycyte1, "beta2", reduction = "umap",  label = FALSE, raster=FALSE) + ggtitle("Beta2 markers") + theme(axis.title = element_text(size = 14)) 
ggsave(filename = "Tany_Scores_notlabelled.png", 
       height = 15, width = 10, dpi=300, 
       plot = plot_grid(p1, p2, p3, p4,p5, p6, nrow = 3))

#These clusters in tanycyte1 (incorrectly labelled tanycytes) renamed as such:
tanycyte1$label[tanycyte1$seurat_clusters %in% c("8","10","14")] <- "alpha 1"
tanycyte1$label[tanycyte1$seurat_clusters == "0"] <- "alpha 2"
tanycyte1$label[tanycyte1$seurat_clusters %in% c("4","18")] <- "beta 1"
tanycyte1$label[tanycyte1$seurat_clusters %in% c("2","19")] <- "beta 2"

#Label transfer were performed by this newly renamed cells
unique(tanycyte1$label)
unique(tany_integrated$label)

#making sure the names are the cell barcodes:
labels_to_transfer <- tanycyte1$label

#Identify which of those cells are also in tany_integrated:
common_cells <- intersect(names(labels_to_transfer), colnames(tany_integrated))

#Overwrite the labels in tany_integrated for those cells:
tany_integrated$label[common_cells] <- labels_to_transfer[common_cells]

# Set cluster order
desired_order <- c("alpha 1","alpha 2","beta 1","beta 2") 
tany_integrated$label <- factor(tany_integrated$label, levels = desired_order)

#Save the integrated dataset
saveRDS(tany_integrated, "/data/pharma_macrophages/Tanybase_integrated_080525.rds")