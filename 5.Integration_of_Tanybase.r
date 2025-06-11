##### Integration of the Combined Dataset #####

# Remove all assays except RNA
DefaultAssay(tanycyte) <- "RNA"
assays_to_keep <- "RNA"
tany_integrated@assays <- tanycyte@assays[names(tanycyte@assays) %in% assays_to_keep]

#QC filtering and measuring the percentage of mt- and Rps, Rpl genes
tanycyte[["percent.mt"]] <- PercentageFeatureSet(tanycyte, pattern = "^mt-")
tanycyte[["percent.rb"]] <- PercentageFeatureSet(tanycyte, pattern = "^Rps|^Rpl")

tanycyte <- subset(
  tanycyte,
  subset = percent.mt <= 15 &
           nCount_RNA >= 600 &
           nFeature_RNA >= 200
)

#Combining low cell number (<100 cells) datasets 
dataset_counts <- table(tanycyte$Dataset) # To check cell numbers from each dataset
low_count_datasets <- names(dataset_counts[dataset_counts < 100])
tanycyte$Dataset_grouped <- ifelse(tanycyte$Dataset %in% low_count_datasets, "Others", tanycyte$Dataset)
table(tanycyte$Dataset_grouped)

#Data were converted to seurat v5
tanycyte <- UpdateSeuratObject(object = tanycyte) 
tanycyte[["RNA"]] <- as(object = tanycyte[["RNA"]], Class = "Assay5")
class(tanycyte[["RNA"]]) #Says Assay5
DefaultAssay(tanycyte) <- "RNA"
Assays(tanycyte) #Says "RNA"

#Integration Preperation
obj1<-tanycyte
obj1[["RNA"]] <- split(obj1[["RNA"]], f = obj1$Dataset_grouped)
obj1 <- ScaleData(obj1)
obj1 <- NormalizeData(obj1)
obj1 <- FindVariableFeatures(obj1)
obj1 <- RunPCA(obj1, npcs = 30, verbose = T)

#CCA Integration
obj1 <- IntegrateLayers(
  object = obj1, method = CCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.cca",
  verbose = TRUE
)

#Rejoining the layers after integration
obj1[["RNA"]] <- JoinLayers(obj1[["RNA"]]) 
DefaultAssay(obj1) <- "RNA"

#Scaling, and Dimensionality Reduction
obj1 <- ScaleData(obj1)
obj1 <- RunPCA(obj1, npcs = 50)
obj1 <- FindNeighbors(obj1, reduction = "integrated.cca", dims = 1:30)
obj1 <- FindClusters(obj1, resolution = 1)
obj1 <- RunUMAP(obj1, reduction = "integrated.cca", dims = 1:30, reduction.name = "umap")

tany_integrated<- obj1#Rename obj1 as tany_integrated 

#QC Plots to check integration quality.
p1<-ElbowPlot(obj1)
ggsave(filename = "elbowplot.pdf", 
       height = 10, width = 10, p1)

p2<-DimPlot(obj1, group.by = "Dataset_grouped")
ggsave(filename = "umap.pdf", 
       height = 5, width = 7, p2)

p3<-DimPlot(obj1, group.by = "label")
ggsave(filename = "umap1.pdf", 
       height = 5, width = 7, p3)

p4<-DimPlot(obj1, group.by = "Age")
ggsave(filename = "umap2.pdf", 
       height = 5, width = 7, p4)

p5<-DimPlot(obj1, group.by = "Diet")
ggsave(filename = "umap3.pdf", 
       height = 5, width = 7, p5)

#Additional plots for QC

p1 <- DimPlot(tany_integrated, group.by = "label", reduction = "umap", shuffle = TRUE, label = TRUE, raster=FALSE) + ggtitle("1") + theme(axis.title = element_text(size = 14)) + NoLegend()
p2 <- DimPlot(tany_integrated, group.by = "cluster_name", reduction = "umap", shuffle = TRUE, label = TRUE, raster=FALSE) + ggtitle("2") + theme(axis.title = element_text(size = 14)) + NoLegend()
ggsave(filename = "Tany_UMAP_mapped.png", 
       height = 5, width = 10, dpi=300, 
       plot = plot_grid(p1, p2,nrow = 1))

p1 <- DimPlot(tany_integrated, group.by = "label", reduction = "umap", shuffle = TRUE, label = TRUE, raster=FALSE) + ggtitle("1") + theme(axis.title = element_text(size = 14)) + NoLegend()
p2 <- DimPlot(tany_integrated, group.by = "seurat_clusters", reduction = "umap", shuffle = TRUE, label = TRUE, raster=FALSE) + ggtitle("2") + theme(axis.title = element_text(size = 14)) + NoLegend()
ggsave(filename = "Tany_UMAP_1res.png", 
       height = 5, width = 10, dpi=300, 
       plot = plot_grid(p1, p2,nrow = 1))

p1 <- DimPlot(tany_integrated, group.by = "label", reduction = "umap", shuffle = TRUE, label = TRUE, raster=FALSE) + ggtitle("1") + theme(axis.title = element_text(size = 14)) + NoLegend()
p2 <- DimPlot(tany_integrated, group.by = "cluster_name", reduction = "umap", shuffle = TRUE, label = TRUE, raster=FALSE) + ggtitle("2") + theme(axis.title = element_text(size = 14)) + NoLegend()
p3 <- DimPlot(tany_integrated, group.by = "Age", reduction = "umap", shuffle = TRUE, label = TRUE, raster=FALSE) + ggtitle("2") + theme(axis.title = element_text(size = 14)) + NoLegend()
p4 <- DimPlot(tany_integrated, group.by = "Diet", reduction = "umap", shuffle = TRUE, label = TRUE, raster=FALSE) + ggtitle("2") + theme(axis.title = element_text(size = 14)) + NoLegend()
p5 <- DimPlot(tany_integrated, group.by = "DataOrigin", reduction = "umap", shuffle = TRUE, label = TRUE, raster=FALSE) + ggtitle("2") + theme(axis.title = element_text(size = 14)) + NoLegend()
p6 <- DimPlot(tany_integrated, group.by = "Dataset", reduction = "umap", shuffle = TRUE, label = TRUE, raster=FALSE) + ggtitle("2") + theme(axis.title = element_text(size = 14)) + NoLegend()
ggsave(filename = "Tany_Int_QualityCheck.png", 
       height = 15, width = 10, dpi=300, 
       plot = plot_grid(p1, p2, p3, p4,p5, p6, nrow = 3))

#Save the integrated dataset
saveRDS(obj1, "/data/pharma_macrophages/Tanycyte_integrated_080525.rds")