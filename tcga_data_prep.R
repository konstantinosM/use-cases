source("install_and_load_libraries.R")
# For this use case we will use the TCGA-PRAD query which contains data from Prostate Adenocarcinoma patients
# Step 1: Create query for gene expression and methylation data quantified with HTSeq from Prostate Adenocarcinoma patients ####
query_TCGA_counts <- GDCquery(
  project = "TCGA-PRAD",
  data.category = "Transcriptome Profiling",
  experimental.strategy = "RNA-Seq",
  workflow.type = "HTSeq - Counts"
)
# Get results of the query
prad_results_counts <- getResults(query_TCGA_counts)

# Step 2: Select data for analysis ####
# For our analysis we will use 498 Primary solid Tumor vs. 52 Solid Tissue Normal
# Primary tumor counts
primary_tumor_counts <- filter(prad_results_counts, sample_type == "Primary Tumor")
# Solid tissue normal counts
solid_tissue_normal_counts <- filter(prad_results_counts, sample_type == "Solid Tissue Normal")
# Step 3: Download data ####
query_primary_tumor_counts <- GDCquery(
  project = "TCGA-PRAD",
  data.category = "Transcriptome Profiling",
  experimental.strategy = "RNA-Seq",
  workflow.type = "HTSeq - Counts",
  barcode = UUIDtoBarcode(id_vector = primary_tumor_counts$id, from_type = "file_id")$associated_entities.entity_submitter_id
)
GDCdownload(query = query_primary_tumor_counts)

query_solid_tissue_normal_counts <- GDCquery(
  project = "TCGA-PRAD",
  data.category = "Transcriptome Profiling",
  experimental.strategy = "RNA-Seq",
  workflow.type = "HTSeq - Counts",
  barcode = UUIDtoBarcode(id_vector = solid_tissue_normal_counts$id, from_type = "file_id")$associated_entities.entity_submitter_id
)
GDCdownload(query = query_solid_tissue_normal_counts)

# Step 4: Load data into R and prepare count matrices ####
# Count data
# Case
data_primary_tumor_counts <- GDCprepare(query_primary_tumor_counts)
data_primary_tumor_counts <- assay(data_primary_tumor_counts)
# Control
data_solid_tissue_normal_counts <- GDCprepare(query_solid_tissue_normal_counts)
data_solid_tissue_normal_counts <- assay(data_solid_tissue_normal_counts)
# Merge in one matrix controls vs. disease
# TODO counts <- merge(x = data_solid_tissue_normal_counts, y = data_primary_tumor_counts, by = "row.names") %>% column_to_rownames(var = "Row.names")
counts <- readRDS("use_case_data/gdc_data/data/raw_counts.rds")
# Step 5: DE analysis: TMM normalization and Z-score computation ####
# TMM normalization of raw counts using edgeR and creation of z-score matrix
# Define contrast groups
contrast <- c(rep("solid_tissue_normal", 52), rep("primary_tumor_counts", 498))
# Save count matrix as dge_list
dge_list <- DGEList(counts = counts, group = contrast)

# Compute normalization factors with TMM
tmm_normalization_factors <- calcNormFactors(dge_list, method = "TMM")

# Normalize counts
norm_counts <- cpm(tmm_normalization_factors)

# Filter out genes that have very low counts across the samples
keep <- filterByExpr(y = norm_counts, min.count = 2, group = contrast)
norm_counts <- norm_counts[keep, ]

# Compute z-score of case samples
# TODO z_score_matrix <- compute_z_scores(norm_counts, controls = c(1:52), cases = c(53:550))
z_score_matrix <- readRDS("use_case_data/gdc_data/data/z_score_counts.rds")
# Since we are going to use the STRING interaction network in this use cases we need to
# substitute the EnsemblIDS with STRING ids
z_score_matrix <- tibble::rownames_to_column(data.frame(z_score_matrix))
# Define which graph to use. In our case we will use 9606(Homo sapiens) version 11
string_db <- STRINGdb$new(version = "11",
                          species = 9606, 
                          input_directory = "",
                          score_threshold = 800)

z_score_matrix <- string_db$map(my_data_frame = z_score_matrix,
                                my_data_frame_id_col_names = "rowname",
                                removeUnmappedRows = TRUE, takeFirst = TRUE)
z_score_matrix["rowname"] <- NULL
# Move STRING_ids to first column
z_score_matrix <- z_score_matrix[, c(499, (1:ncol(z_score_matrix))[-499])]

# Step 6: Determine optimal cutoff ####
z_score_1.5 <- z_score_matrix[-1]
z_score_2 <- z_score_matrix[-1]
z_score_3 <- z_score_matrix[-1]
z_score_4 <- z_score_matrix[-1]
z_score_5 <- z_score_matrix[-1]
# Cutoff |temp|> 1.5
z_score_1.5[z_score_1.5 >= 1.5 | z_score_1.5 <= -1.5] <- 1
z_score_1.5[z_score_1.5 != 1] <- 0
# Cutoff |temp|> 2
z_score_2[z_score_2 >= 2 | z_score_2 <= -2] <- 1
z_score_2[z_score_2 != 1] <- 0
# Cutoff |temp|> 3
z_score_3[z_score_3 >= 3 | z_score_3 <= -3] <- 1
z_score_3[z_score_3 != 1] <- 0
# Cutoff |temp|> 4
z_score_4[z_score_4 >= 4 | z_score_4 <= -4] <- 1
z_score_4[z_score_4 != 1] <- 0
# Cutoff |temp|> 5
z_score_5[z_score_5 >= 5 | z_score_5 <= -5] <- 1
z_score_5[z_score_5 != 1] <- 0

cutoff_list <- c("±1.5", "±2", "±3", "±4", "±5")
z_score_cutoff_comparison <- tibble(
  Z_score_cutoff = factor(cutoff_list, levels = unique(cutoff_list)),
  Genes = c(
    floor(mean(colSums(z_score_1.5))),
    floor(mean(colSums(z_score_2))),
    floor(mean(colSums(z_score_3))),
    floor(mean(colSums(z_score_4))),
    floor(mean(colSums(z_score_5)))
  )
)

z_score_comparison <- ggplot(data = z_score_cutoff_comparison, aes(x = Z_score_cutoff, y = Genes)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_text(aes(label = Genes), vjust = 1.5, color = "black", position = position_dodge(0.9), size = 5) +
  labs(
    title = "Primaty tumor versus Solid tissue normal samples from prostate cancer patients",
    x = "Z-score cutoff",
    y = "Average DEGs per sample",
    subtitle = "|Z_score| > X",
    caption = "Considering genes in the 490 primary tumor samples."
  ) +
  theme(legend.position = "top", plot.title = element_text(size = 18), text = element_text(size = 15))

#ggsave(filename = "~/Desktop/z_score_comparison.png", z_score_comparison, width = 14, height = 8)

# Save count matrices
# Z = 2
counts_matrix_z_2 <- data.frame(id = z_score_matrix$STRING_id, z_score_2)
saveRDS(counts_matrix_z_2, "use_case_data/gdc_data/data/counts_matrix_z_2.rds")
# Z = 3
counts_matrix_z_3 <- data.frame(id = z_score_matrix$STRING_id, z_score_3)
saveRDS(counts_matrix_z_3, "use_case_data/gdc_data/data/counts_matrix_z_3.rds")
# Z = 4
counts_matrix_z_4 <- data.frame(id = z_score_matrix$STRING_id, z_score_4)
saveRDS(counts_matrix_z_4, "use_case_data/gdc_data/data/counts_matrix_z_4.rds")

# Save biological network
human_string_network <- string_db$get_graph()
saveRDS(human_string_network, "use_case_data/gdc_data/graphs/human_string_network_800.rds")