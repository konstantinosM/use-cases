# Step 0: Install and load required packages ####
# Install packages
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
if (!requireNamespace("TCGAutils", quietly = TRUE)) BiocManager::install("TCGAutils")
if (!requireNamespace("TCGAbiolinks", quietly = TRUE)) BiocManager::install("TCGAbiolinks")
if (!requireNamespace("curatedTCGAData", quietly = TRUE)) BiocManager::install("curatedTCGAData")


# Load libraries
# Gene expression datasets
library("GEOquery")
# BioGRID interactions
library("simpIntLists")
# STRING interactions
library("STRINGdb")
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
# Working with data from TCGA
library("TCGAutils")
library("TCGAbiolinks")
library("curatedTCGAData")

