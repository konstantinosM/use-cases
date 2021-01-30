# Step 1: Load data ####
pathways_of_interest <- readRDS(file = "use_case_data/geo_data/comparison/selected_ines_pathways.rds")
indicator_matrix <- readRDS("use_case_data/geo_data/data/indicator_matrix_lfc1_p0,001.rds")
# Step 2: Create gene-to-GO mapping annotation ####
gene_universe <- indicator_matrix$hgnc_symbol
GOs <- translate(gene_universe , from = org.Hs.egALIAS2EG, to = org.Hs.egGO)
example <- pathways_of_interest$`K-1-L1-8-union`
myIntersetingGenes <- example@nodes$node
geneList <- factor(as.integer(gene_universe%in%myIntersetingGenes))
names(geneList) <- gene_universe

# Biological process ####
GO_data <- new("topGOdata", ontology = "BP", allGenes = geneList, annot = annFUN.gene2GO, gene2GO = GOs)

classic_fs <- runTest(GO_data, algorithm = "classic",statistic = "fisher")
classic_ks <- runTest(GO_data, algorithm = "classic",statistic = "ks")
#resultKS_elim <- runTest(GO_data, algorithm = "elim", statistic = "ks")
#weight01_ks<- runTest(GO_data, algorithm = "weight01", statistic = "ks")
table_bp <-  GenTable(GO_data, classicFisher = classic_fs, classicKS = classic_ks, orderBy = "classicFisher")
  
# Molecular function ####  
GO_data <- new("topGOdata", ontology = "MF", allGenes = geneList, annot = annFUN.gene2GO, gene2GO =  GOs)
classic_fs <- runTest(GO_data, algorithm = "classic",statistic = "fisher")
classic_ks <- runTest(GO_data, algorithm = "classic",statistic = "ks")
table_bp <-  GenTable(GO_data, classicFisher = classic_fs, classicKS = classic_ks, orderBy = "classicFisher")

  # allRes <- GenTable(GO_data,
  #                  classicFisher = resultFisher,  elimKS = resultKS_elim, weight01KS = weight01_ks,
  #                  orderBy = "weight01KS", ranksOf = "classicKS", topNodes = 10)