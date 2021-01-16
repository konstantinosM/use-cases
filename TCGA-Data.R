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
prad_results_methylation <- getResults(query_TCGA_methylation)
# Step 2: Filtering ####
# Remove cases which are not in both queries
keep_counts <- prad_results_counts$cases.submitter_id %in% prad_results_methylation$cases.submitter_id
keep_methylation <- prad_results_methylation$cases.submitter_id %in% prad_results_counts$cases.submitter_id

prad_results_counts <- prad_results_counts[keep_counts, ]
prad_results_methylation <- prad_results_methylation[keep_methylation, ]

# Now we check if all the cases have the same amount of entries
equal <- table(prad_results_counts$cases.submitter_id) == table(prad_results_methylation$cases.submitter_id)

# And we keep all the entries that have the same amount of cases
keep <- row.names(equal)[equal == TRUE]
prad_results_counts <- prad_results_counts[prad_results_counts$cases.submitter_id %in% keep, ]
prad_results_methylation <- prad_results_methylation[prad_results_methylation$cases.submitter_id %in% keep, ]

k <- data.frame(names = c("b", "a", "c", "z"), numbers = c(1, 2, 3, 4))
keeper <- c(a = TRUE, b = FALSE, c = FALSE, z = TRUE)
