# SARS-CoV-2 data from GEO
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

# Step 1: Fetch records from GEO ####
# In this case a GEO series record is fetched
gse_147507 <- getGEO("GSE147507", GSEMatrix = TRUE)
show(gse_147507)

# Sometimes the data owners do not load the data into GEO (when assayData has 0 features)
# but only provides supplementary files as in this case.
# In that case the supplementrary files can be download as follows.
getGEOSuppFiles("GSE147507")

# The function saved the supplementary data of the GEO series in the current working directory as a folder with the name of the series
# No we just have to read the correct file as a data frame to obtain the raw counts
gse_147507_raw_counts_human <- as.data.frame.matrix(read.delim("GSE147507/GSE147507_RawReadCounts_Human.tsv.gz"), )

# Step 2: Select datasets ####
# Samples with SARS-CoV-2 infected NHBE cells and mock treated NHBE cells
NHBE_raw_counts <- gse_147507_raw_counts_human[, c(
  "X",
  "Series1_NHBE_Mock_1",
  "Series1_NHBE_Mock_2",
  "Series1_NHBE_Mock_3",
  "Series1_NHBE_SARS.CoV.2_1",
  "Series1_NHBE_SARS.CoV.2_2",
  "Series1_NHBE_SARS.CoV.2_3"
)]
# Samples with HealthyLungBiopsy vs Covid Lung
lung_biopsy_raw_counts <- gse_147507_raw_counts_human[, c(
  "X",
  "Series15_HealthyLungBiopsy_1",
  "Series15_HealthyLungBiopsy_2",
  "Series15_COVID19Lung_1",
  "Series15_COVID19Lung_2"
)]
# Step 3.1: DE analysis: Run DESeq2 on NHBE_raw_counts counts to get DEGs ####
# Set gene name as rownames and remove gene name col
rownames(NHBE_raw_counts) <- NHBE_raw_counts$X
NHBE_raw_counts <- NHBE_raw_counts[-1]

# Create a DESeq dataset object from the count matrix and the colData
coldata <- data.frame(condition = factor(c(rep("Mock", 3), rep("SARS.CoV.2", 3))), type = factor(c(rep("single-end", 3 + 3))))
rownames(coldata) <- colnames(NHBE_raw_counts)
dds <- DESeqDataSetFromMatrix(countData = NHBE_raw_counts, colData = coldata, design = ~condition)

# Filter out genes that have less than two reads
dds <- dds[rowSums(DESeq2::counts(dds)) > 1, ]

# Run DESeq and get results
dds <- DESeq(dds)

# Contrast the samples infected with SARS.CoV.2 to the Mock infected samples
# Important: The order in which the conditions are specified is important
results <- results(dds, contrast = c("condition", "SARS.CoV.2", "Mock"))
# Step 3.2: DE analysis: TMM normalization of raw counts using edgeR and creation of z-score matrix ####
# Save count matrix as dge_list
dge_list <- DGEList(counts = lung_biopsy_raw_counts[, -1], group = c("HealthyLung", "HealthyLung", "COVID19Lung", "COVID19Lung"))
# Compute normalization factors with TMM
tmm_normalization_factors <- calcNormFactors(dge_list, method = "TMM")
# Normalize counts
norm_counts <- cpm(tmm_normalization_factors)
rownames(norm_counts) <- lung_biopsy_raw_counts$X
# Filter out genes that have very low counts across the samples
keep <- filterByExpr(y = norm_counts, min.count = 2, group = c("HealthyLung", "HealthyLung", "COVID19Lung", "COVID19Lung"))
norm_counts <- norm_counts[keep, ]
# Compute z-score of case samples
z_score_matrix <- compute_z_scores(norm_counts, controls = c(1, 2), cases = c(3, 4))
# Step 4.1: Find cutoffs for genes in NHBE samples ####
# Shrink LFC values
results <- lfcShrink(dds, contrast = c("condition", "SARS.CoV.2", "Mock"), res = results, type = "normal")
# Visualize with volcano plots to determine good cutoffss
fc_cutoff <- 1
pCutoff <- 0.0001
# TODO automatization
for (variable in vector) {
  EnhancedVolcano(results,
    lab = rownames(results),
    title = "SARS-CoV-2 versus Mock infected NHBE cells",
    subtitle = paste("P_ADJ ≤,", pCutoff, " and ", "|Log2(FoldChange)| ≥", fc_cutoff, sep = ""),
    x = "log2FoldChange",
    y = "padj",
    FCcutoff = fc_cutoff,
    pCutoff = pCutoff,
    ylab = bquote(~ -Log[10] ~ "(" ~ italic(P_ADJ) ~ ")"),
    xlab = bquote(~ -Log[2] ~ "(" ~ italic(FoldChange) ~ ")")
  )
}

# Exrtact results for specific log2FoldChange and p_adj cutoff
# Get genes that have a p_adj <= 0.0001
p_adjusted_vals <- results$padj <= pCutoff
p_adjusted_vals[is.na(p_adjusted_vals)] <- FALSE
# Get differentially expressed genes with
deg_deseq <- rownames(results[(results$log2FoldChange <= -fc_cutoff | results$log2FoldChange >= fc_cutoff) & p_adjusted_vals, ])

# Step 5: Prepare biological interaction network ####
# Prepare biological interaction network
human_biogrid_network <- findInteractionList(organism = "human", idType = "Official")
# KeyPathwayMineR support two types of input either a sif file or an iGraph object
# For this use case we will utilie iGraph objects
from <- c()
to <- c()
for (entry in human_biogrid_network) {
  from <- c(from, rep(entry$name, length(entry$interactors)))
  to <- c(to, entry$interactors)
} 
edges <-data.frame(from = from, to = to)
human_biogrid_network <- graph_from_data_frame(d = edges, directed=FALSE)
# Remove duplicate edges and loops
human_biogrid_network <- simplify(human_biogrid_network)
# Step 5: downstream analysis with KPM #####
# Prepare indicator matrix
nhbe_indicator_matrix <- data.frame(names = row.names(NHBE_raw_counts), de = c(0))
# Set DEGs to 1
nhbe_indicator_matrix[nhbe_indicator_matrix$names %in% deg_deseq, ][, 2] <- 1
# Set kpm options
kpm_options(
  execution = "Local",
  strategy = "INES",
  algorithm = "Greedy",
  l_min = 20,
  k_min = 5
)
# Execute remote run by using a custom graph_file
local_example_1 <- kpm(graph_file = sample_network, indicator_matrices = huntington_disease_up)

# Visualize the results with shiny
visualize_result(local_example_1)
