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
if (!requireNamespace("igraph", quietly = TRUE)) install.packages("igraph")
if (!requireNamespace("igraph", quietly = TRUE)) install.packages("igraph")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")



# Load libraries
# Gene expression datasets
library("GEOquery")
# BioGRID interactions
library("simpIntLists")
# To save the biological interactions
library("igraph")
# Differential expression analysis
library("DESeq2")
library("edgeR")
# Downstream Analysis
library("KeyPathwayMineR")
# Visualization
library("EnhancedVolcano")
library("ggplot2")

