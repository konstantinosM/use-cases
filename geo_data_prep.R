# SARS-CoV-2 data from GEO
# Install and load R packages
source("install_and_load_libraries.R")
# Step 1: Fetch study GSE147507 from GEO with 7 comparisons ####
getGEOSuppFiles("GSE147507")
# The function saved the supplementary data of the GEO series in the current working directory as
# a folder with the name of the series. No we just have to read the correct file as a data frame to obtain the raw counts
gse_147507_raw_counts_human <-read.delim("GSE147507/GSE147507_RawReadCounts_Human.tsv.gz")
colnames(gse_147507_raw_counts_human)[1] <- "hgnc_symbol"

# Step 1.1: Samples with SARS-CoV-2 infected NHBE cells and mock treated NHBE cells ####
NHBE_series1_raw_counts <- gse_147507_raw_counts_human[, c(
  "hgnc_symbol",
  "Series1_NHBE_Mock_1",
  "Series1_NHBE_Mock_2",
  "Series1_NHBE_Mock_3",
  "Series1_NHBE_SARS.CoV.2_1",
  "Series1_NHBE_SARS.CoV.2_2",
  "Series1_NHBE_SARS.CoV.2_3"
)]

# Step 1.2.: Samples with SARS-CoV-2 infected and mock treated A549 cells. Series 2####
A549_series2_raw_counts <- gse_147507_raw_counts_human[, c(
  "hgnc_symbol",
  "Series2_A549_Mock_1",
  "Series2_A549_Mock_2",
  "Series2_A549_Mock_3",
  "Series2_A549_SARS.CoV.2_1",
  "Series2_A549_SARS.CoV.2_2",
  "Series2_A549_SARS.CoV.2_3"
)]
# Step 1.3: Samples with SARS-CoV-2 infected and mock treated A549 cells. Series 5####
A549_series5_raw_counts <- gse_147507_raw_counts_human[, c(
  "hgnc_symbol",
  "Series5_A549_Mock_1",
  "Series5_A549_Mock_2",
  "Series5_A549_Mock_3",
  "Series5_A549_SARS.CoV.2_1",
  "Series5_A549_SARS.CoV.2_2",
  "Series5_A549_SARS.CoV.2_3"
)]

# Step 1.4: Samples with SARS-CoV-2 infected and mock treated A549-ACE2 cells. Series 6####
A549_ACE2_series6_raw_counts <- gse_147507_raw_counts_human[, c(
  "hgnc_symbol",
  "Series6_A549.ACE2_Mock_1",
  "Series6_A549.ACE2_Mock_2",
  "Series6_A549.ACE2_Mock_3",
  "Series6_A549.ACE2_SARS.CoV.2_1",
  "Series6_A549.ACE2_SARS.CoV.2_2",
  "Series6_A549.ACE2_SARS.CoV.2_3"
)]
# Step 1.5: Samples with SARS-CoV-2 infected and mock treated A549-ACE2 cells. Series 16####
A549_ACE2_series16_raw_counts <- gse_147507_raw_counts_human[, c(
  "hgnc_symbol",
  "Series16_A549.ACE2_Mock_1",
  "Series16_A549.ACE2_Mock_2",
  "Series16_A549.ACE2_Mock_3",
  "Series16_A549.ACE2_SARS.CoV.2_1",
  "Series16_A549.ACE2_SARS.CoV.2_2",
  "Series16_A549.ACE2_SARS.CoV.2_3"
)]
# Step 1.6: Samples with SARS-CoV-2 infected and mock treated Calu3 cells. Series 7####
Calu3_raw_count_series7 <- gse_147507_raw_counts_human[, c(
  "hgnc_symbol",
  "Series7_Calu3_Mock_1",
  "Series7_Calu3_Mock_2",
  "Series7_Calu3_Mock_3",
  "Series7_Calu3_SARS.CoV.2_1",
  "Series7_Calu3_SARS.CoV.2_2",
  "Series7_Calu3_SARS.CoV.2_3"
)]

# Step 1.7: Samples with HealthyLungBiopsy vs Covid Lung ####
lung_biopsy_raw_counts_series15 <- gse_147507_raw_counts_human[, c(
  "hgnc_symbol",
  "Series15_HealthyLungBiopsy_1",
  "Series15_HealthyLungBiopsy_2",
  "Series15_COVID19Lung_1",
  "Series15_COVID19Lung_2"
)]

# Step 2: Fetch study GSE148729 from GEO with 2 comparisons####
getGEOSuppFiles("GSE148729", filter_regex = "*Calu3_totalRNA_readcounts*")
gse_148729_raw_counts_human <- as.data.frame.matrix(read.delim("GSE148729/GSE148729_Calu3_totalRNA_readcounts.tsv.gz"))
# Step 2.1: Samples with SARS-CoV-2 infected and mock treated Calu3 cells. 4 hours after infection ####
Calu3_raw_count_4h <- gse_148729_raw_counts_human[, c(
  "gene_id",
  "Calu3_totalRNA.S2.4h.A",
  "Calu3_totalRNA.S2.4h.B",
  "Calu3_totalRNA.mock.4h.A",
  "Calu3_totalRNA.mock.4h.B"
)]
# Remove version number from Ensembl IDS
Calu3_raw_count_4h$gene_id <- unlist(lapply(X = Calu3_raw_count_4h$gene_id, FUN = sub, pattern = "\\.\\d+$", replacement = ""))

# Step 2.2: Samples with SARS-CoV-2 infected and mock treated Calu3 cells. 24 hours after infection ####
Calu3_raw_count_24h <- gse_148729_raw_counts_human[, c(
  "gene_id",
  "Calu3_totalRNA.S2.24h.A",
  "Calu3_totalRNA.S2.24h.B",
  "Calu3_totalRNA.mock.24h.A",
  "Calu3_totalRNA.mock.24h.B"
)]
# Remove version number from Ensembl IDS
Calu3_raw_count_24h$gene_id <- unlist(lapply(X = Calu3_raw_count_24h$gene_id, FUN = sub, pattern = "\\.\\d+$", replacement = ""))

# Step 2.3: Convert Ensembl ids to hgnc symbol
# Select a BioMart databses
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
ensembl_ids <- Calu3_raw_count_4h$gene_id
identifier_conversion_table <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = ensembl_ids,
  mart = ensembl
)
# Since Calu3_raw_count_24h$gene_id == Calu3_raw_count_4h$gene_id we can use the same identifier_conversion_table for both datasets
Calu3_raw_count_4h <- merge(x = Calu3_raw_count_4h, y = identifier_conversion_table, by.x = "gene_id", by.y = "ensembl_gene_id")[-1]
Calu3_raw_count_24h <- merge(x = Calu3_raw_count_24h, y = identifier_conversion_table, by.x = "gene_id", by.y = "ensembl_gene_id")[-1]
# Remove all entries that were not mapped to a hgnc symbol
Calu3_raw_count_4h <- Calu3_raw_count_4h[Calu3_raw_count_4h$hgnc_symbol != "", ] %>% select(hgnc_symbol, everything())
Calu3_raw_count_24h <- Calu3_raw_count_24h[Calu3_raw_count_24h$hgnc_symbol != "", ] %>% select(hgnc_symbol, everything())

# Step 3: Fetch study GSE148697 and GSE148696 from GEO with 2 comparisons ####
getGEOSuppFiles("GSE148697")
getGEOSuppFiles("GSE148696")
# Step 3.1: Samples with SARS-CoV-2 infected and mock treated HPSC-derived Lung organoids. GSE148697 ####
gse_148697_raw_counts_human <- as.data.frame.matrix(read.delim("GSE148697/GSE148697_counts.txt.gz")) 
colnames(gse_148697_raw_counts_human)[1] <- "hgnc_symbol"
# Step 3.2: Samples with SARS-CoV-2 infected and mock treated HPSC-derived Colonic organoids. GSE148696 ####
gse_148696_raw_counts_human <- as.data.frame.matrix(read.delim("GSE148696/GSE148696_counts.txt.gz"))
colnames(gse_148696_raw_counts_human)[1] <- "hgnc_symbol"
# Step 4:Fetch study GSE164073 from GEO with 3 comparisons ####
getGEOSuppFiles("GSE164073")
gse_164073_raw_counts_human <- as.data.frame.matrix(read.delim("GSE164073/GSE164073_Eye_count_matrix.csv.gz", sep = ","))
colnames(gse_164073_raw_counts_human)[1] <- "hgnc_symbol"
# Step 4.1: Cornea Samples SARS-CoV-2 infected and mock treated ####
cornea_raw_counts <- gse_164073_raw_counts_human[, c(
  "hgnc_symbol",
  "MW1_cornea_mock_1",
  "MW2_cornea_mock_2",
  "MW3_cornea_mock_3",
  "MW4_cornea_CoV2_1",
  "MW5_cornea_CoV2_2",
  "MW6_cornea_CoV2_3"
)]

# Step 4.2: Limbus Samples SARS-CoV-2 infected and mock treated ####
limbus_raw_counts <- gse_164073_raw_counts_human[, c(
  "hgnc_symbol",
  "MW7_limbus_mock_1",
  "MW8_limbus_mock_2",
  "MW9_limbus_mock_3",
  "MW10_limbus_CoV2_1",
  "MW11_limbus_CoV2_2",
  "MW12_limbus_CoV2_3"
)]

# Step 4.3: Sdera Samples SARS-CoV-2 infected and mock treated ####
sclera_raw_counts <- gse_164073_raw_counts_human[, c(
  "hgnc_symbol",
  "MW13_sclera_mock_1",
  "MW14_sclera_mock_2",
  "MW15_sclera_mock_3",
  "MW16_sclera_CoV2_1",
  "MW17_sclera_CoV2_2",
  "MW18_sclera_CoV2_3"
)]

# Step 5: Unify hgnc_symbol ####
# List with all comparisons
comparisons <- list(
  NHBE_series1_raw_counts, A549_series2_raw_counts, A549_series5_raw_counts, A549_ACE2_series6_raw_counts, A549_ACE2_series16_raw_counts, Calu3_raw_count_series7, lung_biopsy_raw_counts_series15,
#  Calu3_raw_count_4h, Calu3_raw_count_24h,
  gse_148697_raw_counts_human, gse_148696_raw_counts_human,
  cornea_raw_counts, limbus_raw_counts, sclera_raw_counts
)
# Get ids that are available in every dataset
ids <- comparisons[[1]]$hgnc_symbol

j <- 2
while (j <= length(comparisons)) {
  ids <- intersect(x = ids,y = comparisons[[j]]$hgnc_symbol)
  j <- j + 1
}

# For every matrix only keep genes from the list ids
# And move hgnc_symbols to rownames
j <- 1
while (j <= length(comparisons)) {
  comparisons[[j]] <- comparisons[[j]][comparisons[[j]]$hgnc_symbol%in%ids,]
  rownames(comparisons[[j]]) <- comparisons[[j]]$hgnc_symbol
  comparisons[[j]] <- comparisons[[j]][-1]
  j <- j + 1
}
# Step 6: DE-Analysis with DESeq2 ####
# Run DESeq2 on NHBE_raw_counts counts to get DEGs 
# Create a DESeq dataset object from the count matrix and the colData
comparisons_results <- list()
j <- 1
while (j <= length(comparisons)) {
  coldata <- data.frame(condition = factor(c(rep("Mock", ncol(comparisons[[j]])/2 ), rep("SARS.CoV.2", ncol(comparisons[[j]])/2))))
  rownames(coldata) <- colnames(comparisons[[j]])
  dds <- DESeqDataSetFromMatrix(countData = comparisons[[j]], colData = coldata, design = ~condition)
  # Run DESeq and get results
  dds <- DESeq(dds)
  # Filter out genes that have less than two reads
  #dds <- dds[rowSums(DESeq2::counts(dds)) > 1, ]
  # Contrast the samples infected with SARS.CoV.2 to the Mock infected samples
  # Important: The order in which the conditions are specified is important
  comparisons_results[j] <- results(dds, contrast = c("condition", "SARS.CoV.2", "Mock"))

  j <- j + 1
}

