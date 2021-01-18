if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("TCGAutils")
BiocManager::install("TCGAutils")
BiocManager::install("curatedTCGAData")

library("TCGAutils")
library("TCGAbiolinks")
library("curatedTCGAData")
library("tidyverse")


# For this use case we will use the TCGA-PRAD query which cointaint data from Prostate Adenocarcinoma patients
# Step1: Create query for gene expression and methylation data quantified with HTSeq from Prostate Adenocarcinoma patients ####
query_TCGA_counts <- GDCquery(
  project = "TCGA-PRAD",
  data.category = "Transcriptome Profiling",
  experimental.strategy = "RNA-Seq",
  workflow.type = "HTSeq - Counts"
)

# Get results of the query
prad_results_counts <- getResults(query_TCGA_counts)

# Query methylation data quantified with RNA-Seq from Prostate Adenocarcinoma patients
query_TCGA_methylation <- GDCquery(
  project = "TCGA-PRAD",
  data.category = "DNA Methylation"
)

# Get results of the query
# prad_results_methylation <- getResults(query_TCGA_methylation)
# # Step 2: Filtering ####
# # Remove cases which are not in both queries
# keep_counts <- prad_results_counts$cases.submitter_id %in% prad_results_methylation$cases.submitter_id
# keep_methylation <- prad_results_methylation$cases.submitter_id %in% prad_results_counts$cases.submitter_id
# 
# prad_results_counts <- prad_results_counts[keep_counts, ]
# prad_results_methylation <- prad_results_methylation[keep_methylation, ]
# 
# # Now we check if all the cases have the same amount of entries
# equal <- table(prad_results_counts$cases.submitter_id) == table(prad_results_methylation$cases.submitter_id)
# 
# # And we keep all the entries that have the same amount of cases
# keep <- row.names(equal)[equal == TRUE]
# prad_results_counts <- prad_results_counts[prad_results_counts$cases.submitter_id %in% keep, ]
# prad_results_methylation <- prad_results_methylation[prad_results_methylation$cases.submitter_id %in% keep, ]

# Step3 : Select data for analysis ####
# For our analysis we will use 25 Primary solid Tumor vs. 25 Solid Tissue Normal
# Primary tumor counts
primary_tumor_counts <- filter(prad_results_counts, sample_type == "Primary Tumor")[1:25, ]
# Solid tissue normal counts
solid_tissue_normal_counts <- filter(prad_results_counts, sample_type == "Solid Tissue Normal")[1:25, ]
# For the same cases fetch also methylation data
# Primary tumor methylation
primary_tumor_methylation <- filter(prad_results_methylation, sample_type == "Primary Tumor")
primary_tumor_methylation <- primary_tumor_methylation[primary_tumor_methylation$cases.submitter_id %in% primary_tumor_counts$cases.submitter_id, ]

# Solid tissues methylation
solid_tissue_normal_methylation <- filter(prad_results_methylation, sample_type == "Solid Tissue Normal")
solid_tissue_normal_methylation <- solid_tissue_normal_methylation[solid_tissue_normal_methylation$cases.submitter_id %in% solid_tissue_normal_counts$cases.submitter_id, ]

# The primary_tumor_methylation query has one additional count with FFPE == TRUE which has to be removed
primary_tumor_methylation <- filter(primary_tumor_methylation, is_ffpe != TRUE)

# Step 4: Download data ####
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

query_primary_tumor_methylation <- GDCquery(
  project = "TCGA-PRAD",
  data.category = "DNA Methylation",
  barcode = UUIDtoBarcode(id_vector = primary_tumor_methylation$id, from_type = "file_id")$associated_entities.entity_submitter_id
)
GDCdownload(query = query_primary_tumor_methylation)

query_solid_tissue_normal_methylation <- GDCquery(
  project = "TCGA-PRAD",
  data.category = "DNA Methylation",
  barcode = UUIDtoBarcode(id_vector = solid_tissue_normal_methylation$id, from_type = "file_id")$associated_entities.entity_submitter_id
)
GDCdownload(query = query_solid_tissue_normal_methylation)
# Step 5: Load data into R and prepare count matrices####
# Count data
# Case
data_primary_tumor_counts <- GDCprepare(query_primary_tumor_counts, )
data_primary_tumor_counts <- assay(data_primary_tumor_counts)
# Control
data_solid_tissue_normal_counts <- GDCprepare(query_solid_tissue_normal_counts)
data_solid_tissue_normal_counts <- assay(data_solid_tissue_normal_counts)
# Merge in one matrix controls vs. disease
counts <- merge(x = data_solid_tissue_normal_counts, y = data_primary_tumor_counts, by = "row.names") %>%
  column_to_rownames(var = "Row.names")
# Methylation data
# Case
data_primary_tumor_methylation <- GDCprepare(query_primary_tumor_methylation)
data_primary_tumor_methylation <- assay(data_primary_tumor_methylation)
# Control
data_solid_tissue_normal_methylation <- GDCprepare(query_solid_tissue_normal_methylation)
data_solid_tissue_normal_methylation <- assay(data_solid_tissue_normal_methylation)
# Merge controls vs. disease
methylation <- merge(x = data_solid_tissue_normal_methylation, y = data_primary_tumor_methylation, by = "row.names") %>%
  column_to_rownames(var = "Row.names")
