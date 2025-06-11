##### Unifying Tany-seq and Hypomap Column Names #####

tanycyte <- combined

#Before everything else we should fix column names.

tanycyte$orig.ident <- tanycyte$Sample_ID
tanycyte$orig.ident[tanycyte$Sample == "Blackshaw_Chow_1"] <- "Blackshaw_Chow_1"
tanycyte$orig.ident[tanycyte$Sample == "Blackshaw_Chow_2"] <- "Blackshaw_Chow_2"
tanycyte$orig.ident[tanycyte$Sample == "Blackshaw_Chow_3"] <- "Blackshaw_Chow_3"

tanycyte$orig.ident[tanycyte$Sample == "Deng_Chow_1"] <- "Deng_Chow_1"
tanycyte$orig.ident[tanycyte$Sample == "Deng_Chow_2"] <- "Deng_Chow_2"
tanycyte$orig.ident[tanycyte$Sample == "Deng_Chow_3"] <- "Deng_Chow_3"

tanycyte$Dataset[tanycyte$Origin == "Blackshaw"] <- "Blackshaw"
tanycyte$Dataset[tanycyte$Origin == "Deng"] <- "Deng"

meta_cols <- c("orig.ident", "Cell_ID", "DataOrigin","Dataset","Age","Diet","inferred_sex")  # Adjust as needed
tanycyte@meta.data <- tanycyte@meta.data[, meta_cols, drop = FALSE]

# Remove all assays except RNA
DefaultAssay(tanycyte) <- "RNA"
assays_to_keep <- "RNA"
tanycyte@assays <- tanycyte@assays[names(tanycyte@assays) %in% assays_to_keep]

#Label transfer from tany-seq and hypomap databases.

#In Hypomap, tanycytes were labelled as following: "C185-132: Tgfb2.Tanycytes","C185-133: Prr16.Tanycytes","C185-134: Frzb.Tanycytes", "C185-131: Scn7a.Tanycytes". These are analogous to alpha 1, alpha 2, beta 1 and beta 2 respectively.

hypomap$label <- hypomap$Author_Class_Curated

hypomap$label[hypomap$C185_named == "C185-132: Tgfb2.Tanycytes"] <- "alpha 1"
hypomap$label[hypomap$C185_named == "C185-133: Prr16.Tanycytes"] <- "alpha 2"
hypomap$label[hypomap$C185_named == "C185-134: Frzb.Tanycytes"] <- "beta 1"
hypomap$label[hypomap$C185_named == "C185-131: Scn7a.Tanycytes"] <- "beta 2"


#Labels were transferred from tany-seq as it is.
tanycyte$label <- hypomap$label
tanycyte$label <- tanyseq$label

# Check the result and save
head(tanycyte)

saveRDS(tanycyte, "/data/pharma_macrophages/combined_Tanycyte_190225.rds")