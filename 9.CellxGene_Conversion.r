##### Converting for CellxGene Compatible Format #####
library(Seurat)
library(reticulate)
library(anndata)
library(Matrix)
library(dplyr)

anndata <- import("anndata")

#Load Seurat object
tanybase <- readRDS("/data/pharma_macrophages/Tanybase_integrated_080525.rds")

#Get normalized data (data slot)
norm_mat <- GetAssayData(tanybase, assay = "RNA", slot = "data")

#Get raw counts
raw_mat <- GetAssayData(tanybase, assay = "RNA", slot = "counts")

#Metadata (obs)
obs <- tanybase@meta.data
rownames(obs) <- colnames(norm_mat)

#Genes (var)
var <- data.frame(symbol = rownames(norm_mat), row.names = rownames(norm_mat))

#Embeddings (obsm)
embeddings <- list()
for (reduction in c("umap", "pca")) {
  if (reduction %in% names(tanybase@reductions)) {
    embeddings[[paste0("X_", reduction)]] <- Embeddings(tanybase, reduction = reduction)
  }
}

#Create AnnData
adata <- anndata$AnnData(
  X = Matrix::t(norm_mat),   # cells × genes
  obs = obs,
  var = var,
  obsm = reticulate::r_to_py(embeddings)
)

#Add raw counts in .raw (cells × genes)
adata$raw <- anndata$AnnData(
  X = Matrix::t(raw_mat),
  var = var
)

#Save to h5ad
adata$write_h5ad("tanybase.h5ad", compression = "gzip")