# SARS-CoV-2 data from GEO
# Install and load R packages
source("install_and_load_libraries.R")
# Step 1: Fetch study GSE147507 from GEO with 7 comparisons ####
getGEOSuppFiles("GSE147507")
# The function saved the supplementary data of the GEO series in the current working directory as
# a folder with the name of the series. No we just have to read the correct file as a data frame to obtain the raw counts
gse_147507_raw_counts_human <- read.delim("GSE147507/GSE147507_RawReadCounts_Human.tsv.gz")
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

# Step 4.3: Sclera Samples SARS-CoV-2 infected and mock treated ####
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
names(comparisons) <- c(
  "SARS-CoV-2 versus Mock infected NHBE cells",
  "SARS-CoV-2 versus Mock infected A549 cells (Series 2)",
  "SARS-CoV-2 versus Mock infected A549 cells (Series 5)",
  "SARS-CoV-2 versus Mock infected A549-ACE2 cells (Series 6)",
  "SARS-CoV-2 versus Mock infected A549-ACE2 cells (Series 16)",
  "SARS-CoV-2 versus Mock infected Calu3 cells",
  "Covid 19 lung versus Healthy lung biopsy",
  "SARS-CoV-2 versus Mock infected HPSC-derived lung organoids",
  "SARS-CoV-2 versus Mock infected HPSC-derived colonic organoids",
  "SARS-CoV-2 versus Mock infected Cornea sample",
  "SARS-CoV-2 versus Mock infected Limbus sample",
  "SARS-CoV-2 versus Mock infected Sclera sample"
)

# Get ids that are available in every dataset
ids <- comparisons[[1]]$hgnc_symbol

j <- 2
while (j <= length(comparisons)) {
  ids <- intersect(x = ids, y = comparisons[[j]]$hgnc_symbol)
  j <- j + 1
}

# For every matrix only keep genes from the list ids
# And move hgnc_symbols to rownames
j <- 1
while (j <= length(comparisons)) {
  comparisons[[j]] <- comparisons[[j]][comparisons[[j]]$hgnc_symbol %in% ids, ]
  rownames(comparisons[[j]]) <- comparisons[[j]]$hgnc_symbol
  comparisons[[j]] <- comparisons[[j]][-1]
  j <- j + 1
}
# Step 6: DE-Analysis with DESeq2 ####
dds_list <- list()
j <- 1
while (j <= length(comparisons)) {
  coldata <- data.frame(condition = factor(c(rep("Mock", ncol(comparisons[[j]]) / 2), rep("SARS.CoV.2", ncol(comparisons[[j]]) / 2))))
  rownames(coldata) <- colnames(comparisons[[j]])
  dds <- DESeqDataSetFromMatrix(countData = comparisons[[j]], colData = coldata, design = ~condition)
  # Filter out genes that have less than two reads
  #dds <- dds[rowSums(DESeq2::counts(dds)) > 1, ]
  
  # Run DESeq and get results
  dds_list[[j]] <- DESeq(dds)

  j <- j + 1
}
names(dds_list) <- names(comparisons)
# Step 7: Cutoffs, Volcano and DE genes ####
volcanos <- list()
degs_deseq_list <- list()
fc_cutoff <- 1
pCutoff <- 0.0001
j <- 1
while (j <= length(dds_list)) {
  # Contrast the samples infected with SARS.CoV.2 to the Mock infected samples
  # Important: The order in which the conditions are specified is important
  results <- results(dds_list[[j]], contrast = c("condition", "SARS.CoV.2", "Mock"))
  #results <- lfcShrink(dds_list[[j]], contrast = c("condition", "SARS.CoV.2", "Mock"), res = results, type = "normal")

  p_adjusted_vals <- results$padj <= pCutoff
  p_adjusted_vals[is.na(p_adjusted_vals)] <- FALSE
  deg_deseq <- rownames(results[(results$log2FoldChange <= -fc_cutoff | results$log2FoldChange >= fc_cutoff) & p_adjusted_vals, ])
  degs_deseq_list[[j]] <- deg_deseq
  # Plot volcano plot to asses good cutoffs
  volcanos[[j]] <- EnhancedVolcano(results,
    lab = rownames(results),
    title = names(comparisons)[j],
    subtitle = paste("P_ADJ ≤", pCutoff, " and ", "|Log2(FoldChange)| ≥", fc_cutoff, sep = ""),
    caption = paste0("[Genes] Total = ", nrow(results), " and DEGs = ", length(deg_deseq)),
    x = "log2FoldChange",
    y = "padj",
    FCcutoff = fc_cutoff,
    pCutoff = pCutoff,
    ylab = bquote(~ -Log[10] ~ "(" ~ italic(P_ADJ) ~ ")"),
    xlab = bquote(~ Log[2] ~ "(" ~ italic(FoldChange) ~ ")")
  )
  j <- j + 1
}
names(volcanos) <- names(dds_list)
names(degs_deseq_list) <- names(dds_list)
for (i in c(1:length(volcanos))) {
  ggsave(plot = volcanos[[i]], filename = paste0("~/Desktop/plots/geo_volcanos/", names(volcanos)[i], ".png"))
}
grid <- ggarrange(plotlist = volcanos, ncol = 3, nrow = 4)
ggsave(plot = grid, filename = "~/Desktop/plots/geo_volcanos/grid.png", width = 30, height = 30)

# Step 8: Create indicator matrix  ####
indicator_matrix <- data.frame(hgnc_symbol = ids)

for (i in c(1:length(degs_deseq_list))) {
  indicator_matrix[names(degs_deseq_list)[i]] <- ifelse(indicator_matrix$hgnc_symbol %in% unlist(degs_deseq_list[i]), 1, 0)
}
colnames(indicator_matrix)[-1][1:7] <- paste0(colnames(indicator_matrix)[-1][1:7], " - GSE147507")
colnames(indicator_matrix)[-1][8] <- paste0(colnames(indicator_matrix)[-1][8], " - GSE148697")
colnames(indicator_matrix)[-1][9] <- paste0(colnames(indicator_matrix)[-1][9], " - GSE148696")
colnames(indicator_matrix)[-1][10:12] <- paste0(colnames(indicator_matrix)[-1][9], " - GSE164073")
upset(indicator_matrix,
  sets = colnames(indicator_matrix)[-1],
  sets.bar.color = "#56B4E9",
  order.by = "freq",
  empty.intersections = "on",
  keep.order = TRUE
)
