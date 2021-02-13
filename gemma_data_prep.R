devtools::install_github("PavlidisLab/gemmaAPI.R")
library("gemmaAPI")
# Step 1: Get data from GEMMA ####
# Retrieve available differential expression tests for the given dataset
differential <- datasetInfo("GSE1297", request = "differential")
# Retrieves the differential expression results given the differential_id.
diff_exp_data <- datasetInfo("GSE1297", request = "degs", differential = differential$`59350`$id)
diff_exp_data <- diff_exp_data[["resultset_ID468701.data"]]

# Step 2: Filtering ####

# Remove entries which do not have a NCBI_ID
diff_exp_data <- diff_exp_data[!is.na(diff_exp_data$NCBI_ID), ]

# Remove entries that were mapped to multiple genes
diff_exp_data <- diff_exp_data[str_count(string = diff_exp_data$NCBI_ID, pattern = "\\|") == 0, ]

# Log2 transform fold changes and remove nas
diff_exp_data$FoldChange_incipient_Alzheimer.s.disease <- log2(diff_exp_data$FoldChange_incipient_Alzheimer.s.disease)
diff_exp_data$FoldChange_moderate_Alzheimer.s.disease <- log2(diff_exp_data$FoldChange_moderate_Alzheimer.s.disease)
diff_exp_data$FoldChange_Alzheimer.s.disease_severe <- log2(diff_exp_data$FoldChange_Alzheimer.s.disease_severe)
keep <- !(is.na(diff_exp_data$FoldChange_incipient_Alzheimer.s.disease) | is.na(diff_exp_data$FoldChange_moderate_Alzheimer.s.disease) | is.na(diff_exp_data$FoldChange_Alzheimer.s.disease_severe))
diff_exp_data <- diff_exp_data[keep,]

# Step 3: Determine cutoffs ####
# Set cutoff for p_adj and fc
fc_cutoff <- 0.0
pCutoff <- 0.1
# Get foldchange and pvals for all conditions and determine  degs
# Incipient
incipient_samples <- data.frame(
  ncbi_id = diff_exp_data$NCBI_ID, 
  fold_change = diff_exp_data$FoldChange_incipient_Alzheimer.s.disease,
  pval = diff_exp_data$PValue_incipient_Alzheimer.s.disease
)
pvals <- incipient_samples$pval <= pCutoff
deg_deseq <- incipient_samples[(incipient_samples$log2fold_change <= -fc_cutoff | incipient_samples$log2fold_change  >= fc_cutoff) & pvals, ]$ncbi_id
# Plot volcano plot to asses good cutoffs
EnhancedVolcano(incipient_samples,
                lab = incipient_samples$ncbi_id,
                x = "fold_change",
                y = "pval",
                FCcutoff = fc_cutoff,
                pCutoff = pCutoff
)

# Moderate
fold_change_moderate <- diff_exp_data$FoldChange_moderate_Alzheimer.s.disease
pval_moderate <- diff_exp_data$PValue_moderate_Alzheimer.s.disease

# Severe
fold_change_severe <- diff_exp_data$FoldChange_Alzheimer.s.disease_severe
pval_severe <- diff_exp_data$PValue_Alzheimer.s.disease_severe


# Prepare indicator matrix
nhbe_indicator_matrix <- data.frame(names = row.names(NHBE_raw_counts), de = c(0))
# Set DEGs to 1
nhbe_indicator_matrix[nhbe_indicator_matrix$names %in% deg_deseq, ][, 2] <- 1
