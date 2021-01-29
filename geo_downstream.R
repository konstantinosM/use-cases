if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GOSim")

ines_pathways <- readRDS("use_case_data/geo_data/comparison/selected_ines_pathways.rds")

example1 <- ines_pathways$`K-3-L1-10`
#Obtain the GO terms and their description for a list of genes.
getGOInfo(geneIDs = example1@nodes,)

