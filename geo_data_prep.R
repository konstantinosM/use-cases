# SARS-CoV-2 data from GEO
# Install and load R packages
source("install_and_load_libraries.R")
# Step 1: Fetch study GSE147507 from GEO with 7 comparisons ####
getGEOSuppFiles("GSE147507") 
# The function saved the supplementary data of the GEO series in the current working directory as
# a folder with the name of the series. No we just have to read the correct file as a data frame to obtain the raw counts
gse_147507_raw_counts_human <- as.data.frame.matrix(read.delim("GSE147507/GSE147507_RawReadCounts_Human.tsv.gz"), )

# Step 1.1: Samples with SARS-CoV-2 infected NHBE cells and mock treated NHBE cells ####
NHBE_series1_raw_counts <- gse_147507_raw_counts_human[, c(
  "X",
  "Series1_NHBE_Mock_1",
  "Series1_NHBE_Mock_2",
  "Series1_NHBE_Mock_3",
  "Series1_NHBE_SARS.CoV.2_1",
  "Series1_NHBE_SARS.CoV.2_2",
  "Series1_NHBE_SARS.CoV.2_3"
)]

# Step 1.2.: Samples with SARS-CoV-2 infected and mock treated A549 cells. Series 2#### 
A549_series2_raw_counts <-  gse_147507_raw_counts_human[, c(
  "X",
  "Series2_A549_Mock_1",
  "Series2_A549_Mock_2",
  "Series2_A549_Mock_3",
  "Series2_A549_SARS.CoV.2_1",
  "Series2_A549_SARS.CoV.2_2",
  "Series2_A549_SARS.CoV.2_3"
)]
# Step 1.3: Samples with SARS-CoV-2 infected and mock treated A549 cells. Series 5#### 
A549_series5_raw_counts <-  gse_147507_raw_counts_human[, c(
  "X",
  "Series5_A549_Mock_1",
  "Series5_A549_Mock_2",
  "Series5_A549_Mock_3",
  "Series5_A549_SARS.CoV.2_1",
  "Series5_A549_SARS.CoV.2_2",
  "Series5_A549_SARS.CoV.2_3"
)]

# Step 1.4: Samples with SARS-CoV-2 infected and mock treated A549-ACE2 cells. Series 6#### 
A549_ACE2_series6_raw_counts <-  gse_147507_raw_counts_human[, c(
  "X",
  "Series6_A549.ACE2_Mock_1",
  "Series6_A549.ACE2_Mock_2",
  "Series6_A549.ACE2_Mock_3",
  "Series6_A549.ACE2_SARS.CoV.2_1",
  "Series6_A549.ACE2_SARS.CoV.2_2",
  "Series6_A549.ACE2_SARS.CoV.2_3"
)]
# Step 1.5: Samples with SARS-CoV-2 infected and mock treated A549-ACE2 cells. Series 16#### 
A549_ACE2_series16_raw_counts <-  gse_147507_raw_counts_human[, c(
  "X",
  "Series16_A549.ACE2_Mock_1",
  "Series16_A549.ACE2_Mock_2",
  "Series16_A549.ACE2_Mock_3",
  "Series16_A549.ACE2_SARS.CoV.2_1",
  "Series16_A549.ACE2_SARS.CoV.2_2",
  "Series16_A549.ACE2_SARS.CoV.2_3"
)]
# Step 1.6: Samples with SARS-CoV-2 infected and mock treated Calu3 cells. Series 7#### 
Calu3_raw_count_series7 <- gse_147507_raw_counts_human[, c(
  "X",
  "Series7_Calu3_Mock_1",
  "Series7_Calu3_Mock_2",
  "Series7_Calu3_Mock_3",
  "Series7_Calu3_SARS.CoV.2_1",
  "Series7_Calu3_SARS.CoV.2_2",
  "Series7_Calu3_SARS.CoV.2_3"
)]

# Step 1.7: Samples with HealthyLungBiopsy vs Covid Lung ####
lung_biopsy_raw_counts_series15 <- gse_147507_raw_counts_human[, c(
  "X",
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
# Step 2.2: Samples with SARS-CoV-2 infected and mock treated Calu3 cells. 24 hours after infection ####
Calu3_raw_count_24h <- gse_148729_raw_counts_human[, c(
  "gene_id",
  "Calu3_totalRNA.S2.24h.A",
  "Calu3_totalRNA.S2.24h.B",
  "Calu3_totalRNA.mock.24h.A",
  "Calu3_totalRNA.mock.24h.B"
)]
# Step 3: Fetch study GSE153940 from GEO with 1 comparison ####
getGEOSuppFiles("GSE153940")
# Step 3.1: Samples with SARS-CoV-2 infected and mock treated VeroE6 cells. 24 hours after infection ####
# Unpack tar
untar("GSE153940/GSE153940_RAW.tar", exdir = "GSE153940/")
samples <- paste0("GSE153940/", list.files(path = "GSE153940/", pattern = "*.gz"))s
i <- 1
while(i <= length(samples)){
  if(i==1){
    VeroE6_raw_count_24h <- read.delim(file = samples[[i]])%>% select(gene_id, expected_count)%>% column_to_rownames(var = "gene_id")
  }else{
    VeroE6_raw_count_24h <- data.frame(VeroE6_raw_count_24h, read.delim(file = samples[[i]])%>% select(expected_count))
  }
  i = i + 1
}

colnames(VeroE6_raw_count_24h) <- c("Control1", "Control2", "Control3","SARS-CoV-1", "SARS-CoV-2", "SARS-CoV-3")












# Step 4: Fetch study GSE148697 and GSE148696 from GEO with 2 comparisons ####
getGEOSuppFiles("GSE148697")
getGEOSuppFiles("GSE148696")
# Step 4.1: Samples with SARS-CoV-2 infected and mock treated HPSC-derived Lung organoids. GSE148697 ####
gse_GSE148697_raw_counts_human <- as.data.frame.matrix(read.delim("GSE148697/GSE148697_counts.txt.gz"))

# Step 4.2: Samples with SARS-CoV-2 infected and mock treated HPSC-derived Colonic organoids. GSE148696 ####
gse_GSE148696_raw_counts_human <- as.data.frame.matrix(read.delim("GSE148696/GSE148696_counts.txt.gz"))











# Step 3.1: DE analysis: NHBE raw counts ####
# Run DESeq2 on NHBE_raw_counts counts to get DEGs 
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

# Step 3.2: DE analysis: Lung biopsy raw counts ####
# TMM normalization of raw counts using edgeR and creation of z-score matrix
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

# Step 4.1: Find cutoffs and create indicator matrix (NHBE) ####
# Shrink LFC values
results <- lfcShrink(dds, contrast = c("condition", "SARS.CoV.2", "Mock"), res = results, type = "normal")
# Set cutoff for p_adj and fc
fc_cutoff <- 0.5
pCutoff <- 0.05
# Exrtact results for specific log2FoldChange and p_adj cutoff
p_adjusted_vals <- results$padj <= pCutoff
p_adjusted_vals[is.na(p_adjusted_vals)] <- FALSE
deg_deseq <- rownames(results[(results$log2FoldChange <= -fc_cutoff | results$log2FoldChange >= fc_cutoff) & p_adjusted_vals, ])

# Plot volcano plot to asses good cutoffs
EnhancedVolcano(results,
  lab = rownames(results),
  title = "SARS-CoV-2 versus Mock infected NHBE cells",
  subtitle = paste("P_ADJ ≤", pCutoff, " and ", "|Log2(FoldChange)| ≥", fc_cutoff, sep = ""),
  caption = paste0("[Genes] Total = ", nrow(results), " and DEGs = ", length(deg_deseq)),
  x = "log2FoldChange",
  y = "padj",
  FCcutoff = fc_cutoff,
  pCutoff = pCutoff,
  ylab = bquote(~ -Log[10] ~ "(" ~ italic(P_ADJ) ~ ")"),
  xlab = bquote(~ -Log[2] ~ "(" ~ italic(FoldChange) ~ ")"),
)

# Prepare indicator matrix
nhbe_indicator_matrix <- data.frame(names = row.names(NHBE_raw_counts), de = c(0))
# Set DEGs to 1
nhbe_indicator_matrix[nhbe_indicator_matrix$names %in% deg_deseq, ][, 2] <- 1
# Step 4.2: Find cutoffs and create indicator matrix (lung biopsy) ####
z_score_1.5 <- z_score_matrix
z_score_2 <- z_score_matrix
z_score_3 <- z_score_matrix
z_score_10 <- z_score_matrix
z_score_50 <- z_score_matrix
# Cutoff |temp|> 1.5
z_score_1.5[z_score_1.5 >= 1.5 | z_score_1.5 <= -1.5] <- 1
z_score_1.5[z_score_1.5 != 1] <- 0
# Cutoff |temp|> 2
z_score_2[z_score_2 >= 2 | z_score_2 <= -2] <- 1
z_score_2[z_score_2 != 1] <- 0
# Cutoff |temp|> 3
z_score_3[z_score_3 >= 3 | z_score_3 <= -3] <- 1
z_score_3[z_score_3 != 1] <- 0
# Cutoff |temp|> 10
z_score_10[z_score_10 >= 10 | z_score_10 <= -10] <- 1
z_score_10[z_score_10 != 1] <- 0
# Cutoff |temp|> 50
z_score_50[z_score_50 >= 50 | z_score_50 <= -50] <- 1
z_score_50[z_score_50 != 1] <- 0

cutoff_list <- c("±1.5", "±1.5", "±2", "±2", "±3", "±3","±10.0", "±10.0","±50.0", "±50.0")
z_score_cutoff_comparison <- data.frame(
  Type = c("DE", "NDE", "DE", "NDE", "DE", "NDE", "DE", "NDE", "DE", "NDE"),
  Z_score_cutoff = factor(cutoff_list,levels = unique(cutoff_list)),
  Genes = c(
    table(z_score_1.5)["1"], table(z_score_1.5)["0"],
    table(z_score_2)["1"], table(z_score_2)["0"],
    table(z_score_3)["1"], table(z_score_3)["0"],
    table(z_score_10)["1"], table(z_score_10)["0"],
    table(z_score_50)["1"], table(z_score_50)["0"]
  )
)

ggplot(data = z_score_cutoff_comparison, aes(x = Z_score_cutoff, y = Genes, fill = Type)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_brewer(palette = "Paired") +
  geom_text(aes(label = Genes), vjust = 1.35, color = "black", position = position_dodge(0.9), size = 5) +
  labs(title = "COVID19 versus Healthy lung biopsy", subtitle = "|Z_score| > X", caption = " Considering genes overall 490 primary tumor samlpes.") +
  theme(legend.position = "top", plot.title = element_text(size = 20), text = element_text(size = 15))

lung_biopsy_indicator_matrix <- data_frame(rownames(z_score_50), z_score_50) 

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
edges <- data.frame(from = from, to = to)
human_biogrid_network <- graph_from_data_frame(d = edges, directed = FALSE)
# Remove duplicate edges and loops
human_biogrid_network <- simplify(human_biogrid_network)
# Step 6.1: downstream analysis with KPM for NHBE sample #####
# Greedy GLONE run
kpm_options(
  execution = "Local",
  strategy = "GLONE",
  algorithm = "Greedy",
  use_range_l = TRUE,
  l_min = 10,
  l_step = 5,
  l_max = 50
)

# Execute remote run by using a custom graph_file
glone_results_nhbe <- kpm(graph = human_biogrid_network, indicator_matrices = nhbe_indicator_matrix)

# Visualize the results with shiny
visualize_result(glone_results_nhbe)

reset_options()
# Greedy INES run
kpm_options(
  execution = "Local",
  strategy = "INES",
  algorithm = "Greedy",
  use_range_k = TRUE,
  l_min = 0,
  k_step = 10,
  k_max = 50
)

innes_results_nhbe <- kpm(graph = human_biogrid_network, indicator_matrices = nhbe_indicator_matrix)

# Visualize the results with shiny
visualize_result(innes_results_nhbe)

# Step 6.2: downstream analysis with KPM for lung biopsy samples ####
##GLONE##
reset_options()
kpm_options(
  execution = "Local",
  strategy = "GLONE",
  algorithm = "Greedy",
  use_range_l = TRUE,
  l_min = 10,
  l_step = 5,
  l_max = 50
)

# Execute remote run by using a custom graph_file
glone_results_lung_biopsy <- kpm(graph = human_biogrid_network, indicator_matrices = nhbe_indicator_matrix)

# Visualize the results with shiny
visualize_result(glone_results_lung_biopsy)

##INES##
reset_options()
kpm_options(
  execution = "Local",
  strategy = "INES",
  algorithm = "Greedy",
  use_range_k = TRUE, 
  use_range_l = TRUE,
  l_min = 1,
  l_step = 1,
  l_max = 2,
  k_min = 2,
  k_step = 2,
  k_max = 10
)

# Execute remote run by using a custom graph_file
ines_results_lung_biopsy <- kpm(graph = human_biogrid_network, indicator_matrices = lung_biopsy_indicator_matrix)

# Visualize the results with shiny
visualize_result(ines_results_lung_biopsy)

reset_options()












