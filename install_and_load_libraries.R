# Step 0: Install and load required packages ####
# Install packages ####
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("DESeq2", quietly = TRUE)) BiocManager::install("DESeq2")
if (!requireNamespace("edgeR", quietly = TRUE)) BiocManager::install("edgeR")
if (!requireNamespace("GEOquery", quietly = TRUE)) BiocManager::install("GEOquery")
if (!requireNamespace("EnhancedVolcano", quietly = TRUE)) BiocManager::install("EnhancedVolcano")
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
if (!requireNamespace("KeyPathwayMineR", quietly = TRUE)) devtools::install_github("baumbachlab/keypathwayminer-R", build_vignettes = TRUE)
if (!requireNamespace("simpIntLists", quietly = TRUE)) BiocManager::install("simpIntLists")
if (!requireNamespace("STRINGdb", quietly = TRUE)) BiocManager::install("STRINGdb")
if (!requireNamespace("igraph", quietly = TRUE)) install.packages("igraph")
if (!requireNamespace("tidyverse", quietly = TRUE)) install.packages("tidyverse")
if (!requireNamespace("ggalt", quietly = TRUE)) install.packages("ggalt")
if (!requireNamespace("TCGAutils", quietly = TRUE)) BiocManager::install("TCGAutils")
if (!requireNamespace("TCGAbiolinks", quietly = TRUE)) BiocManager::install("TCGAbiolinks")
if (!requireNamespace("curatedTCGAData", quietly = TRUE)) BiocManager::install("curatedTCGAData")
if (!requireNamespace("biomaRt", quietly = TRUE)) BiocManager::install("biomaRt")
if (!requireNamespace("UpSetR", quietly = TRUE)) install.packages("UpSetR")
if (!requireNamespace("ggpubr", quietly = TRUE)) install.packages("ggpubr")
if (!requireNamespace("topGO", quietly = TRUE)) BiocManager::install("topGO")
if (!requireNamespace("AnnotationFuncs", quietly = TRUE)) BiocManager::install("AnnotationFuncs")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("Rgraphviz")

# Load libraries ####
# Gene expression datasets
library("GEOquery")
# BioGRID interactions
library("simpIntLists")
# STRING interactions and enrichment analysis
library("STRINGdb")
# Gene set enrichment analysis
library("topGO")
# To save the biological interactions
library("igraph")
# Differential expression analysis
library("DESeq2")
library("edgeR")
# Downstream Analysis
library("KeyPathwayMineR")
# Visualization and data management
library("EnhancedVolcano")
library("tidyverse")
library("UpSetR")
library("ggpubr")
library("ggalt")
# Working with data from TCGA
library("TCGAutils")
library("TCGAbiolinks")
library("curatedTCGAData")
# Convert Ensembl ids to hgnc symbols
library("biomaRt")
# Annotation functions
library("AnnotationFuncs")
library("org.Hs.eg.db")
