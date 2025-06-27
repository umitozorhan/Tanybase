# Tanybase - Combined and Curated Tanycyte Database

# Introduction
Tanycytes are ependymal-like glial cells that play a crucial role in regulating the neuroendocrine system, metabolism, and glucose homeostasis. Tanycytes lie in the wall of the third ventricle and coordinate various functions. Traditionally, tanycytes were classified into four subtypes based on their position in the ventricular wall: alpha 1, alpha 2, beta 1 and beta 2 tanycytes. Recent scRNA-seq studies have investigated the role of tanycytes at the transcriptome level, but they fall short due to the limited number of tanycytic cells in the datasets. To address this shortcoming, two prominent databases, Tany-seq and Hypomap, were compiled by combining leading hypothalamus and tanycyte scRNA-seq studies on this topic. However, the number of tanycytes for each subtype and condition was still limited, even after these noble works. To eliminate this pitfall and help future researchers in their quest to identify the role of tanycytes to a greater extent, we combined these two databases into a single mouse tanycyte-specific database called tanybase. Tanybase have a total of 20610 cells from 4 tanycytic subtypes, 5 different diets and 3 different ages. 

In short, tanybase is formed by combining two scRNAseq databases of mouse tanycytes, Hypomap and Tany-seq. Hypomap and Tany-seq databases were combined using shared parameters, including diet, age, tanycyte subclusters, original identity and dataset. Unique cells in each database were combined, and cells with more than 15% mitochondrial transcripts, a total RNA count below 600, and fewer than 200 unique features were removed. We included only cells that were labelled as tanycyte according to the Allen Institute MapMyCells software (RRID:SCR_024672). For integration, tanybase was split into different layers for each dataset containing more than 100 cells and datasets that have less than 100 cells were combined (‘others’). Tanybase was scaled, normalized and variable features were determined using Seurat functions. Principle-component analysis was performed as well as integration for batch correction using canonical correlation analysis. After the integration, layers were re-joined. For dimensionality reduction, we performed UMAP clustering with FindClusters and FindNeighbors functions. Labels were transferred from corresponding databases; unlabelled cells were labelled using module scores and the markers provided in Tany-seq. For Tanybase, we used Seurat (v.5.2.0), AnnData, Matrix packages. 

![Tany_UMAP_mapped](https://github.com/user-attachments/assets/ae4f8358-aebf-4bc0-b83e-431592ad7182)

Figure 1. UMAP graph of Tanybase database. 1 shows the canonical tanycyte subtypes, 2 shows the MapMyCell subclusters of tanycytes.

# System and Software Requirements
  # System
Standard modern computer 16+ GB RAM and 8 Cores with 20 gb space. Runtime depends on RAM and Core speeds. Whole analysis could be done under 1 hour.
  # Software
Conda environment or R studio. Conda environment packages and requirements were given at Enviromentlist.txt.
Alternatively R studio with R base >= 4.4 should work. This work run previously on debian based Linux machine.
  # Necessary packages
seurat (5.1.0), stringr (1.5.1), dplyr (1.1.4), matrix(1.7_2), reticulate (1.40.0), anndata (0.11.4), ggplot2 (3.5.1), cowplot (1.1.3).

1.install.packages("Seurat", version='5.1.0')
2.install.packages("stringr", version='1.5.1')
3.install.packages("dplyr", version='1.1.4')
4.if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_version("Matrix", version = "1.7-2")
5.install.packages("reticulate", version='1.40.0')
6.install.packages("anndata", version='0.11.4')
7.install.packages("ggplot2", version='3.5.1')
8.install.packages("cowplot", version='1.1.3')

This will take 5-10 minutes depending on internet connection stability and speed. For conda environment that could take up to 30 min because of extra packages attached to the  Enviromentlist.txt.

# Data Availability
# Source Data
HypoMap Mouse https://www.repository.cam.ac.uk/items/8f9c3683-29fd-44f3-aad5-7acf5e963a75
Tany-Seq https://data.mendeley.com/datasets/w8yw2c92jg/1
# Final Data
özorhan, ümit (2025), “Tanybase - Combined and Curated Tanycyte Database”, Mendeley Data, V2, doi: 10.17632/p6jkzkpdd6.2

# Instructions
Follow the instruction given at each section. These should reproduce the same or very similar results.
1.Data_Collection,
2.Combining_Tany-seq_and_Hypomap,
3.Unifying_Tany-seq_and_Hypomap,
4.Tanycyte_Mapping _and_CleanUp,
5.Integration_of_Tanybase,
6.Labelling_unlabelled_Tanycytes,
7.DE_Genes_Calculation_Tanybase,
8.Plots_for_the_Publication,
9.CellxGene_Conversion.

