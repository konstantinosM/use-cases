# SARS-CoV-2 data from GEO
# Install and load R packages
source("install_and_load_libraries.R")
  # Series 1: Fetch study GSE147507 from GEO with 7 comparisons ####
  # getGeoSuppFiles("GSE147507")
  gse_147507_raw_counts_human <- read.delim("GSE147507/GSE147507_RawReadCounts_Human.tsv.gz")
  colnames(gse_147507_raw_counts_human)[1] <- "hgnc_symbol"
  
  # Comp. 1.1: Samples with SARS-CoV-2 infected NHBE cells and mock treated NHBE cells ####
  NHBE_series1_raw_counts <- gse_147507_raw_counts_human[, c(
    "hgnc_symbol",
    "Series1_NHBE_Mock_1",
    "Series1_NHBE_Mock_2",
    "Series1_NHBE_Mock_3",
    "Series1_NHBE_SARS.CoV.2_1",
    "Series1_NHBE_SARS.CoV.2_2",
    "Series1_NHBE_SARS.CoV.2_3"
  )]
  
  # Comp. 1.2.: Samples with SARS-CoV-2 infected and mock treated A549 cells. Series 2####
  A549_series2_raw_counts <- gse_147507_raw_counts_human[, c(
    "hgnc_symbol",
    "Series2_A549_Mock_1",
    "Series2_A549_Mock_2",
    "Series2_A549_Mock_3",
    "Series2_A549_SARS.CoV.2_1",
    "Series2_A549_SARS.CoV.2_2",
    "Series2_A549_SARS.CoV.2_3"
  )]
  # Comp. 1.3: Samples with SARS-CoV-2 infected and mock treated A549 cells. Series 5####
  A549_series5_raw_counts <- gse_147507_raw_counts_human[, c(
    "hgnc_symbol",
    "Series5_A549_Mock_1",
    "Series5_A549_Mock_2",
    "Series5_A549_Mock_3",
    "Series5_A549_SARS.CoV.2_1",
    "Series5_A549_SARS.CoV.2_2",
    "Series5_A549_SARS.CoV.2_3"
  )]
  
  # Comp. 1.4: Samples with SARS-CoV-2 infected and mock treated A549-ACE2 cells. Series 6####
  A549_ACE2_series6_raw_counts <- gse_147507_raw_counts_human[, c(
    "hgnc_symbol",
    "Series6_A549.ACE2_Mock_1",
    "Series6_A549.ACE2_Mock_2",
    "Series6_A549.ACE2_Mock_3",
    "Series6_A549.ACE2_SARS.CoV.2_1",
    "Series6_A549.ACE2_SARS.CoV.2_2",
    "Series6_A549.ACE2_SARS.CoV.2_3"
  )]
  # Comp. 1.5: Samples with SARS-CoV-2 infected and mock treated A549-ACE2 cells. Series 16####
  A549_ACE2_series16_raw_counts <- gse_147507_raw_counts_human[, c(
    "hgnc_symbol",
    "Series16_A549.ACE2_Mock_1",
    "Series16_A549.ACE2_Mock_2",
    "Series16_A549.ACE2_Mock_3",
    "Series16_A549.ACE2_SARS.CoV.2_1",
    "Series16_A549.ACE2_SARS.CoV.2_2",
    "Series16_A549.ACE2_SARS.CoV.2_3"
  )]
  # Comp. 1.6: Samples with SARS-CoV-2 infected and mock treated Calu3 cells. Series 7####
  Calu3_raw_count_series7 <- gse_147507_raw_counts_human[, c(
    "hgnc_symbol",
    "Series7_Calu3_Mock_1",
    "Series7_Calu3_Mock_2",
    "Series7_Calu3_Mock_3",
    "Series7_Calu3_SARS.CoV.2_1",
    "Series7_Calu3_SARS.CoV.2_2",
    "Series7_Calu3_SARS.CoV.2_3"
  )]
  
  # Comp. 1.7: Samples with HealthyLungBiopsy vs Covid Lung ####
  lung_biopsy_raw_counts_series15 <- gse_147507_raw_counts_human[, c(
    "hgnc_symbol",
    "Series15_HealthyLungBiopsy_1",
    "Series15_HealthyLungBiopsy_2",
    "Series15_COVID19Lung_1",
    "Series15_COVID19Lung_2"
  )]
  
  # Series 2: Fetch study GSE151879 from GEO with 3 comparisons####
  # getGeoSuppFiles("GSE151879")
  # Comp. 2.8: Samples with SARS-CoV-2 infected and mock treated HESC-derived cardiomyocytes cells. ####
  hesc_derived_cardiomyocytes_raw_counts <- read.delim(file = "GSE151879/GSE151879_raw_counts_genes.hESC-derived_cardiomyocytes.txt.gz")
  # Comp. 2.9: Samples with SARS-CoV-2 infected and mock treated Human cardiomyocytes cells. ####
  human_cardiomyocytes_raw_counts <- read.delim(file = "GSE151879/GSE151879_raw_counts_genes.Adult_human_cardiomyocytes.txt.gz")
  # Comp. 2.10: Samples with SARS-CoV-2 infected and mock treated macrophages ####
  macrophages_raw_counts <- read.delim(file = "GSE151879/GSE151879_raw_counts_genes.Macrophages.txt.gz")
  # Processing 2.a: Convert Ensembl ids to hgnc symbol ####
  # Select a BioMart databses
  ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
  ensembl_ids <- macrophages_raw_counts$gene_id
  identifier_conversion_table <- getBM(
    attributes = c("ensembl_gene_id", "hgnc_symbol"),
    filters = "ensembl_gene_id",
    values = ensembl_ids,
    mart = ensembl
  )
  # Remove non mapped entries
  identifier_conversion_table <- identifier_conversion_table[identifier_conversion_table$hgnc_symbol!="",]
  # Remove duplicate entries
  identifier_conversion_table <-  identifier_conversion_table[!identifier_conversion_table$hgnc_symbol%in%identifier_conversion_table$hgnc_symbol[duplicated(identifier_conversion_table$hgnc_symbol)],]
  
  # Since the gene_ids of all datasets are the same we can use the same identifier_conversion_table for all datasets
  hesc_derived_cardiomyocytes_raw_counts <- merge(x = hesc_derived_cardiomyocytes_raw_counts, y = identifier_conversion_table, by.x = "gene_id", by.y = "ensembl_gene_id")[-1]
  human_cardiomyocytes_raw_counts <- merge(x = human_cardiomyocytes_raw_counts, y = identifier_conversion_table, by.x = "gene_id", by.y = "ensembl_gene_id")[-1]
  macrophages_raw_counts <- merge(x = macrophages_raw_counts, y = identifier_conversion_table, by.x = "gene_id", by.y = "ensembl_gene_id")[-1]
  # Remove all entries that were not mapped to a hgnc symbol and remove duplicate entries
  hesc_derived_cardiomyocytes_raw_counts <- hesc_derived_cardiomyocytes_raw_counts[hesc_derived_cardiomyocytes_raw_counts$hgnc_symbol != "", ] %>% dplyr::select(hgnc_symbol, everything())
  human_cardiomyocytes_raw_counts <- human_cardiomyocytes_raw_counts[human_cardiomyocytes_raw_counts$hgnc_symbol != "", ] %>% dplyr::select(hgnc_symbol, everything())
  macrophages_raw_counts <- macrophages_raw_counts[macrophages_raw_counts$hgnc_symbol != "", ] %>% dplyr::select(hgnc_symbol, everything())
  
  # Series 3,4: Fetch study GSE148697 and GSE148696 from GEO with 2 comparisons ####
  # getGeoSuppFiles("GSE148697")
  # getGeoSuppFiles("GSE148696")
  # Comp. 3,4.N: Samples with SARS-CoV-2 infected and mock treated HPSC-derived Lung organoids. GSE148697 ####
  gse_148697_raw_counts_human <- as.data.frame.matrix(read.delim("GSE148697/GSE148697_counts.txt.gz"))
  colnames(gse_148697_raw_counts_human)[1] <- "hgnc_symbol"
  # Comp. 3,4.11: Samples with SARS-CoV-2 infected and mock treated HPSC-derived Colonic organoids. GSE148696 ####
  gse_148696_raw_counts_human <- as.data.frame.matrix(read.delim("GSE148696/GSE148696_counts.txt.gz"))
  colnames(gse_148696_raw_counts_human)[1] <- "hgnc_symbol"
  # Series 5:Fetch study GSE164073 from GEO with 3 comparisons ####
  # getGeoSuppFiles("GSE164073")
  gse_164073_raw_counts_human <- as.data.frame.matrix(read.delim("GSE164073/GSE164073_Eye_count_matrix.csv.gz", sep = ","))
  colnames(gse_164073_raw_counts_human)[1] <- "hgnc_symbol"
  # Comp. 5.12: Cornea Samples SARS-CoV-2 infected and mock treated ####
  cornea_raw_counts <- gse_164073_raw_counts_human[, c(
    "hgnc_symbol",
    "MW1_cornea_mock_1",
    "MW2_cornea_mock_2",
    "MW3_cornea_mock_3",
    "MW4_cornea_CoV2_1",
    "MW5_cornea_CoV2_2",
    "MW6_cornea_CoV2_3"
  )]
  
  # Comp. 5.13: Limbus Samples SARS-CoV-2 infected and mock treated ####
  limbus_raw_counts <- gse_164073_raw_counts_human[, c(
    "hgnc_symbol",
    "MW7_limbus_mock_1",
    "MW8_limbus_mock_2",
    "MW9_limbus_mock_3",
    "MW10_limbus_CoV2_1",
    "MW11_limbus_CoV2_2",
    "MW12_limbus_CoV2_3"
  )]
  
  # Comp. 5.14: Sclera Samples SARS-CoV-2 infected and mock treated ####
  sclera_raw_counts <- gse_164073_raw_counts_human[, c(
    "hgnc_symbol",
    "MW13_sclera_mock_1",
    "MW14_sclera_mock_2",
    "MW15_sclera_mock_3",
    "MW16_sclera_CoV2_1",
    "MW17_sclera_CoV2_2",
    "MW18_sclera_CoV2_3"
  )]
  
  # Series 6: Fetch study GSE160435 from GEO with 1 comparison ####
  # getGeoSuppFiles("GSE160435")
  gse_160435_raw_counts_human <- as.data.frame.matrix(read.delim("GSE160435/GSE160435_count.csv.gz", sep = ","))
  colnames(gse_160435_raw_counts_human)[1] <- "hgnc_symbol"
  # Comp. 6.15:  Samples with SARS-CoV-2 infected and mock treated AT2 cells####
  AT2_raw_counts <- gse_160435_raw_counts_human[, c(
    "hgnc_symbol",
    "cc01.20mock1",
    "cc01.20mock2",
    "cc03.19mock",
    "cc03.20mock1",
    "cc03.20mock2",
    "cc01.20covid1",
    "cc01.20covid2",
    "cc03.19covid",
    "cc03.20covid1",
    "cc03.20covid2"
  )]
  
  # Remove duplicate entries 
  AT2_raw_counts <-  AT2_raw_counts[!AT2_raw_counts$hgnc_symbol%in%AT2_raw_counts$hgnc_symbol[duplicated(AT2_raw_counts$hgnc_symbol)],]
  
  # Series 7: Fetch study GSE157852 from GEO with 2 comparisons ####
  # getGeoSuppFiles("GSE157852")
  gse_157852_raw_counts_human <- as.data.frame.matrix(read.delim("GSE157852/GSE157852_CPO_RawCounts.txt.gz", sep = " ")) %>% rownames_to_column(var = "hgnc_symbol")
  # Comp. 7.16:  Samples with SARS-CoV-2 infected and mock treated choroid plexus organoids cells 72 hours post infection (hpi)####
  choroid_plexus_72_raw_counts <- gse_157852_raw_counts_human[, c(
    "hgnc_symbol",
    "CPO_Mock_72hpi_S1",
    "CPO_Mock_72hpi_S2",
    "CPO_Mock_72hpi_S3",
    "CPO_SARS.CoV.2_72hpi_S7",
    "CPO_SARS.CoV.2_72hpi_S8",
    "CPO_SARS.CoV.2_72hpi_S9"
  )]
  # Comp. 7.17:  Samples with SARS-CoV-2 infected organoids cells 24hpi and mock treated choroid plexus organoids cells 72 hpi####
  choroid_plexus_24_raw_counts <- gse_157852_raw_counts_human[, c(
    "hgnc_symbol",
    "CPO_Mock_72hpi_S1",
    "CPO_Mock_72hpi_S2",
    "CPO_Mock_72hpi_S3",
    "CPO_SARS.CoV.2_24hpi_S4",
    "CPO_SARS.CoV.2_24hpi_S5",
    "CPO_SARS.CoV.2_24hpi_S6"
  )]
  
  # Series 8: Fetch study GSE152075 from GEO with 1 comparison ####
  # getGeoSuppFiles("GSE152075")
  gse_152075_raw_counts_human <- as.data.frame.matrix(read.delim("GSE152075/GSE152075_raw_counts_GEO.txt.gz", sep = " ")) %>% rownames_to_column(var = "hgnc_symbol")
  # Comp. 8.18 Nasopharyngeal swabs from 430 individuals with SARS-CoV-2 and 54 negative controls ####
  control <- colnames(gse_152075_raw_counts_human)[startsWith(x = colnames(gse_152075_raw_counts_human), prefix = "NEG_")]
  sars_cov_2 <- colnames(gse_152075_raw_counts_human)[startsWith(x = colnames(gse_152075_raw_counts_human), prefix = "POS_")]
  nasopharyngeal_swabs_raw_counts <- gse_152075_raw_counts_human[, c("hgnc_symbol", control, sars_cov_2)]
  
  # Series 19: Fetch study GSE150392 from GEO with 1 comparison ####
  # getGeoSuppFiles("GSE150392")
  gse_150392_raw_counts_human <- as.data.frame.matrix(read.delim("GSE150392/GSE150392_Cov_Mock_Raw_COUNTS.csv.gz", sep = ","))
  # Remove entries without HGNC annotation
  gse_150392_raw_counts_human <- gse_150392_raw_counts_human[gse_150392_raw_counts_human$X %in% gse_150392_raw_counts_human$X[grepl("_", gse_150392_raw_counts_human$X, fixed = TRUE)], ]
  #  Remove ensembl identifiers
  gse_150392_raw_counts_human$X <- gsub(".*_", "", gse_150392_raw_counts_human$X)
  colnames(gse_150392_raw_counts_human)[1] <- "hgnc_symbol"
  # Comp. 9.19:  Samples with SARS-CoV-2 infected and mock treated hiPSC-CMs ####
  hiPSC_CMs_raw_counts <- gse_150392_raw_counts_human[, c("hgnc_symbol", "Mock1", "Mock2", "Mock3", "Cov1", "Cov2", "Cov3")]
  # Remove duplicate entries
  hiPSC_CMs_raw_counts <-  hiPSC_CMs_raw_counts[!hiPSC_CMs_raw_counts$hgnc_symbol%in%hiPSC_CMs_raw_counts$hgnc_symbol[duplicated(hiPSC_CMs_raw_counts$hgnc_symbol)],]
  
# Step 1: Unify hgnc_symbol ####
# List with all comparisons
comparisons <- list(
  NHBE_series1_raw_counts, A549_series2_raw_counts, A549_series5_raw_counts, A549_ACE2_series6_raw_counts, A549_ACE2_series16_raw_counts, Calu3_raw_count_series7, lung_biopsy_raw_counts_series15,
  hesc_derived_cardiomyocytes_raw_counts, human_cardiomyocytes_raw_counts, macrophages_raw_counts,
  gse_148696_raw_counts_human,
  cornea_raw_counts, limbus_raw_counts, sclera_raw_counts,
  AT2_raw_counts,
  choroid_plexus_72_raw_counts, choroid_plexus_24_raw_counts,
  nasopharyngeal_swabs_raw_counts,
  hiPSC_CMs_raw_counts
)
names(comparisons) <- c(
  "SARS-CoV-2 versus Mock infected NHBE cells",
  "SARS-CoV-2 versus Mock infected A549 cells (Series 2)",
  "SARS-CoV-2 versus Mock infected A549 cells (Series 5)",
  "SARS-CoV-2 versus Mock infected A549-ACE2 cells (Series 6)",
  "SARS-CoV-2 versus Mock infected A549-ACE2 cells (Series 16)",
  "SARS-CoV-2 versus Mock infected Calu3 cells",
  "SARS-CoV-2 versus Healthy lung biopsies",

  "SARS-CoV-2 versus Mock infected HESC-derived cardiomyocytes cells",
  "SARS-CoV-2 versus Mock infected Human cardiomyocytes cells",
  "SARS-CoV-2 versus Mock infected Macrophages",

  #"SARS-CoV-2 versus Mock infected HPSC-derived lung organoids",

  "SARS-CoV-2 versus Mock infected HPSC-derived colonic organoids",

  "SARS-CoV-2 versus Mock infected Cornea sample",
  "SARS-CoV-2 versus Mock infected Limbus sample",
  "SARS-CoV-2 versus Mock infected Sclera sample",

  "SARS-CoV-2 versus Mock infected AT2 cells",

  "SARS-CoV-2 versus Mock infected Organoid cells 72 hpi",
  "SARS-CoV-2 versus Mock infected Organoid cells 24 hpi",

  "SARS-CoV-2 versus Negative control nasopharyngeal swabs",
  
  "SARS-CoV-2 versus Mock infected hiPSC-CMs"
)

# Get ids that are available in every dataset
ids <- comparisons[[1]]$hgnc_symbol

j <- 2
while (j <= length(comparisons)) {
  ids <- intersect(x = ids, y = comparisons[[j]]$hgnc_symbol)
  j <- j + 1
}

# For every matrix only keep genes that are contained in the list ids
# And move hgnc_symbols to rownames
j <- 1
while (j <= length(comparisons)) {
  comparisons[[j]] <- comparisons[[j]][comparisons[[j]]$hgnc_symbol %in% ids, ]
  rownames(comparisons[[j]]) <- c() 
  comparisons[[j]] <- tibble::column_to_rownames(.data = comparisons[[j]], var = "hgnc_symbol")
  j <- j + 1
}
# Step 2: DE-Analysis with DESeq2 ####
dds_list <- list()
j <- 1
while (j <= length(comparisons)) {
 #message(paste0("Processing: ",names(comparisons)[j]))
  coldata <- data.frame(condition = factor(c(rep("Mock", ncol(comparisons[[j]]) / 2), rep("SARS.CoV.2", ncol(comparisons[[j]]) / 2))))
  rownames(coldata) <- colnames(comparisons[[j]])
  dds <- DESeqDataSetFromMatrix(countData = comparisons[[j]], colData = coldata, design = ~condition)
  # Filter out genes that have less than two reads
  # dds <- dds[rowSums(DESeq2::counts(dds)) > 1, ]

  # Run DESeq and get results
  dds_list[[j]] <- DESeq(dds)

  j <- j + 1
}
names(dds_list) <- names(comparisons)
saveRDS(dds_list, "use_case_data/geo_data/data/DESeq2_result_list.rds")
# Step 3: Cutoffs, Volcano and DE genes ####
dds_list <- readRDS("use_case_data/geo_data/data/DESeq2_result_list.rds")
volcanos <- list()
degs_deseq_list <- list()
# TODO
fc_cutoff <- 2
pCutoff <- 0.001
j <- 1
while (j <= length(dds_list)) {
  # Contrast the samples infected with SARS.CoV.2 to the Mock infected samples
  # Important: The order in which the conditions are specified is important
  results <- results(dds_list[[j]], contrast = c("condition", "SARS.CoV.2", "Mock"))
  # results <- lfcShrink(dds_list[[j]], contrast = c("condition", "SARS.CoV.2", "Mock"), res = results, type = "normal")

  p_adjusted_vals <- results$padj <= pCutoff
  p_adjusted_vals[is.na(p_adjusted_vals)] <- FALSE
  deg_deseq <- rownames(results[(results$log2FoldChange <= -fc_cutoff | results$log2FoldChange >= fc_cutoff) & p_adjusted_vals, ])
  degs_deseq_list[[j]] <- deg_deseq
  # Plot volcano plot to asses good cutoffs
  volcanos[[j]] <- EnhancedVolcano(results,
    lab = rownames(results),
    title = names(dds_list)[j],
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
  # TODO
  ggsave(plot = volcanos[[i]], filename = paste0("~/Desktop/plots/geo_volcanos_upset/", "LFC2_P0,001/", names(volcanos)[i], ".png"), width = 9.5, height = 7)
}

grid <- ggarrange(plotlist = volcanos[1:9], ncol = 3, nrow = 3)
# TODO
ggsave(plot = grid, filename = "~/Desktop/plots/geo_volcanos_upset/LFC2_P0,001/grid_1_9.png", width = 30, height = 30)

grid <- ggarrange(plotlist = volcanos[10:19], ncol = 3, nrow = 4)
# TODO
ggsave(plot = grid, filename = "~/Desktop/plots/geo_volcanos_upset/LFC2_P0,001/grid_10_19.png", width = 30, height = 30)

# Step 4: Create indicator matrix  ####
indicator_matrix <- data.frame(hgnc_symbol = ids)

for (i in c(1:length(degs_deseq_list))) {
  indicator_matrix[names(degs_deseq_list)[i]] <- ifelse(indicator_matrix$hgnc_symbol %in% unlist(degs_deseq_list[i]), 1, 0)
}
colnames(indicator_matrix)[-1][1:7] <- paste0(colnames(indicator_matrix)[-1][1:7])
colnames(indicator_matrix)[-1][8:10] <- paste0(colnames(indicator_matrix)[-1][8:10])
colnames(indicator_matrix)[-1][11] <- paste0(colnames(indicator_matrix)[-1][11])
colnames(indicator_matrix)[-1][12:14] <- paste0(colnames(indicator_matrix)[-1][12:14])
colnames(indicator_matrix)[-1][15] <- paste0(colnames(indicator_matrix)[-1][15])
colnames(indicator_matrix)[-1][16:17] <- paste0(colnames(indicator_matrix)[-1][16:17])
colnames(indicator_matrix)[-1][18] <- paste0(colnames(indicator_matrix)[-1][18])
colnames(indicator_matrix)[-1][19] <- paste0(colnames(indicator_matrix)[-1][19])
colnames(indicator_matrix) <- str_replace(colnames(indicator_matrix), pattern = "SARS-CoV-2", replacement = "S2")
colnames(indicator_matrix) <- str_replace(colnames(indicator_matrix),pattern = "Mock infected", replacement = "MI")
# TODO
saveRDS(indicator_matrix, "use_case_data/geo_data/data/indicator_matrix_lfc2_p0,001.rds")

upset_plot <- upset(indicator_matrix[,c(2:20)],
  sets = colnames(indicator_matrix)[c(2:20)],
  sets.x.label = "Number of differentially expressed genes \nin the dataset",
  main.bar.color = "black",
  sets.bar.color = "blue",
  order.by = "freq")
# TODO
pdf(file="~/Desktop/plots/geo_volcanos_upset/LFC2_P0,001/upset.pdf", width = 10, height = 7) # or other device
upset_plot
dev.off()

# Step 5: Prepare biological interaction network ###
human_biogrid_network <- findInteractionList(organism = "human", idType = "Official")
# KeyPathwayMineR support two types of input either a sif file or an iGraph object
# For this use case we will utilize iGraph objects
from <- c()
to <- c()
for (entry in human_biogrid_network) {
  from <- c(from, rep(entry$name, length(entry$interactors)))
  to <- c(to, entry$interactors)
}
edges <- data.frame(from = from, to = to)
human_biogrid_network <- graph_from_data_frame(d = edges, directed = FALSE)
# Remove duplicate edges and loops
human_biogrid_network <- simplify(human_biogrid_network)
saveRDS(human_biogrid_network, "use_case_data/geo_data/graphs/human_biogrid_network.rds")