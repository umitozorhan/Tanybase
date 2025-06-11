#####Tanycyte Mapping to Remove Nontanycytic Cells Through MapMyCells (RRID:SCR_024672) #####

#Software were given at https://portal.brain-map.org/atlases-and-data/bkp/mapmycells

#Convert the tanycyte .rds file to .h5ad file
# Extract the raw count matrix
data_matrix <- tanycyte@assays$RNA@counts

# Transpose the data
data_matrixt <- Matrix::t(data_matrix)

# Convert to AnnData format
ad <- AnnData(
  X = data_matrixt,
  obs = data.frame(group = rownames(data_matrixt), row.names = rownames(data_matrixt)),
  var = data.frame(type = colnames(data_matrixt), row.names = colnames(data_matrixt))
)

# Write to compressed .h5ad file
write_h5ad(ad, 'tanycyte.h5ad', compression = 'gzip')

# Check file size (must be <500MB)
print(paste("Size in MB:", round(file.size("tanycyte.h5ad") / 2^20)))

#H5AD file uploaded to https://knowledge.brain-map.org/mapmycells/process/ and software run through following parameters: 
#Reference Taxonomy: 10xWholeMouseBrain(CCN20230722), Mapping Algorithm: HierarchicalMapping

#Results were read from csv file.
mapping <- read.csv("tanycyte_10xWholeMouseBrain(CCN20230722)_HierarchicalMapping_UTC_1746707160353.csv",comment.char="#")
head(data.frame(mapping))

#Since the query data corresponds to dataQC above, we will call it dataQC again
dataQC_h5ad <- read_h5ad('tanycyte.h5ad')
dataQC <- t(as.matrix(dataQC_h5ad$X))
rownames(dataQC) <- rownames(dataQC_h5ad$var)
colnames(dataQC) <- rownames(dataQC_h5ad$obs)

#Put row.names as data colnames and the order to match the data
rownames(mapping) <- mapping$cell_id
mapping <- mapping[colnames(dataQC),]

#Create the Seurat object using mapping info and h5ad file
dataSeurat <- CreateSeuratObject(counts = dataQC, meta.data = mapping)

unique(dataSeurat$subclass_name) #Tanycytes were labelled as "322 Tanycyte NN"

#Move these labels in tanycyte object and only keep the cells labelled with subclass_name "322 Tanycyte NN"

tanycyte$class_name <- dataSeurat$class_name
tanycyte$subclass_name <- dataSeurat$subclass_name
tanycyte$supertype_name <- dataSeurat$supertype_name
tanycyte$cluster_name <- dataSeurat$cluster_name

tanycyte <- subset(x=tanycyte, subset = subclass_name == "322 Tanycyte NN")

# Check the result and save
head(tanycyte)

saveRDS(tanycyte, "/data/pharma_macrophages/combined_Tanycyte_190225.rds") #This dataset includes only the cells that named as tanycyte from MapMyCells, Hypomap and Tany-seq.