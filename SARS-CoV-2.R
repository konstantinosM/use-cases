# SARS-CoV-2 data from GEO
#### Step 1: Fetch records from GEO ####
library("GEOquery")

# In this case a GEO series record is fetched
gse_147507 <- getGEO("GSE147507", GSEMatrix = TRUE)
show(gse_147507)

# Sometimes the data owners do not load the data into GEO (when assayData has 0 features)
# but only provides supplementary files as in this case.
# In that case the supplementrary files can be download as follows.
getGEOSuppFiles("GSE147507")

# The function saved the supplementary data of the GEO series in the current working directory as a folder with the name of the series
# No we just have to read the correct file as a data frame to obtain the raw counts
gse_147507_raw_counts_human <- as.data.frame.matrix(read.delim("GSE147507/GSE147507_RawReadCounts_Human.tsv.gz"),)

#### Step 2: Select samples with SARS-Cov-2 infected NHBE cells and mock treated NHBE cells ####
library("dplyr")
NHBE_raw_counts <- gse_147507_raw_counts_human %>% select(EnsemblId = "X", "Series1_NHBE_Mock_1","Series1_NHBE_Mock_2","Series1_NHBE_Mock_3",
                                                          "Series1_NHBE_SARS.CoV.2_1","Series1_NHBE_SARS.CoV.2_2","Series1_NHBE_SARS.CoV.2_3")
#### Step 3: TMM normalization of raw counts ####
library("edgeR")
# Save count matrix as dge_list
dge_list <- DGEList(counts = NHBE_raw_counts,group = c("control", "control", "control", "sars-cov-2","sars-cov-2","sars-cov-2"))
# Get counts and sample information
dge_list
# Compute normalization factors with TMM
tmm_normalization_factors <- calcNormFactors(dge_list, method = "TMM")
# Normalize counts
norm_counts <- cpm(tmm_normalization_factors)

#### Step 4: Compute z score matrix between control and sars-cov-2 infected  cells ####




