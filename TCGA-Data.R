if (!require("BiocManager")) {
  install.packages("BiocManager")
}
BiocManager::install("GenomicDataCommons")

library(GenomicDataCommons)

# Create manifest file with information on raw data to be used
# Filter gene expression data quantified with HTSeq from Prostate Adenocarcinoma patients
counts = files() %>%
  filter( cases.project.project_id == 'TCGA-PRAD') %>% 
  filter( type == 'gene_expression' ) %>%  
  filter( analysis.workflow_type == 'HTSeq - Counts')  %>%
  manifest()

# Get the first 25 cased ids
cases <- ge_manifest$id[1:25]

# Get counts for the first 25 cases
gdcdata(uuids = cases)

head(ge_manifest)


testin <- as.data.frame.matrix(read.delim(header = FALSE, "/Users/konstantinos/Library/Caches/GenomicDataCommons/10a24f34-ce3b-43c3-8646-6042c2b0e178/ea1ce00b-00f4-4194-97df-6a64227e7bdc.htseq.counts.gz"), )
