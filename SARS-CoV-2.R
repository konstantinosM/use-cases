# SARS-CoV-2 data from GEO
#### Step 0: Install and load required packages ####
# Install packages 
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("DESeq2", quietly = TRUE)) BiocManager::install("DESeq2")
if (!requireNamespace("edgeR", quietly = TRUE)) BiocManager::install("edgeR")
if (!requireNamespace("GEOquery", quietly = TRUE)) BiocManager::install("GEOquery")
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
if (!requireNamespace("KeyPathwayMineR", quietly = TRUE)) devtools::install_github("baumbachlab/keypathwayminer-R", build_vignettes = TRUE)
# Load libraries
library("KeyPathwayMineR")
library("GEOquery")
library("DESeq2")
library("edgeR")
#### Step 1: Fetch records from GEO ####
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

#### Step 2: Select datasets ####
# Samples with SARS-CoV-2 infected NHBE cells and mock treated NHBE cells 
NHBE_raw_counts <- gse_147507_raw_counts_human[,c("X", 
                                                  "Series1_NHBE_Mock_1", 
                                                  "Series1_NHBE_Mock_2", 
                                                    "Series1_NHBE_Mock_3",
                                                    "Series1_NHBE_SARS.CoV.2_1", 
                                                    "Series1_NHBE_SARS.CoV.2_2",
                                                  "Series1_NHBE_SARS.CoV.2_3")]
# Samples with HealthyLungBiopsy vs Covid Lung
biopsy_raw_counts <- gse_147507_raw_counts_human[,c("X", 
                                                    "Series15_HealthyLungBiopsy_1", 
                                                    "Series15_HealthyLungBiopsy_2", 
                                                    "Series15_COVID19Lung_1",
                                                    "Series15_COVID19Lung_2")]
# Step 3: Differential expression analysis ####
# Step 3.1: Run DESeq2 on NHBE_raw_counts counts to get DEGs ####
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
results <- results(dds, contrast = c("condition", 'SARS.CoV.2', 'Mock'))

# Step 3.2 TMM normalization of raw counts using edgeR and creation of z-score matrix ####


library("edgeR")
# Save count matrix as dge_list
dge_list <- DGEList(counts = NHBE_raw_counts[,-1], group = c("control", "control", "control", "sars-cov-2", "sars-cov-2", "sars-cov-2"))
# Get counts and sample information
dge_list
# Compute normalization factors with TMM
tmm_normalization_factors <- calcNormFactors(dge_list, method = "TMM")
# Normalize counts
norm_counts <- data.frame(EnsemblId = NHBE_raw_counts$EnsemblId, cpm(tmm_normalization_factors))
controls <- c(2, 3, 4)
cases <- c(5, 6, 7)















padj <- result$padj < 0.05
padj[is.na(padj)]<- FALSE
rownames(NHBE_raw_counts)[padj]
sum(result$padj < 0.05, na.rm=TRUE)


#### Step 4: Compute z score matrix between control and sars-cov-2 infected  cells ####
p_vals <- computePValues(countMatrix = biopsy_raw_counts[-1] , n = 3, m = 3, method = "DeSeq2" )
# plot on same grid, each series colored differently -- 
# good if the series have same scale
ggplot(df, aes(time, value)) + geom

boxplot(logcounts, 
        xlab="", 
        ylab="Log2(Counts)",
        las=2,
        col=statusCol)
# Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(as.matrix(logcounts)), col="blue")

