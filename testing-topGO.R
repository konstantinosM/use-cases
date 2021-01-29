# Step 1: Load data ####
pathways_of_interest <- readRDS(file = "use_case_data/geo_data/comparison/selected_ines_pathways.rds")
indicator_matrix <- readRDS("use_case_data/geo_data/data/indicator_matrix_lfc1_p0,001.rds")
# Step 2: Create gene-to-GO mapping annotation ####
gene_universe <- indicator_matrix$hgnc_symbol
GOs <- translate(gene_universe , from = org.Hs.egALIAS2EG, to = org.Hs.egGO)
# Step 3: 
pickGO(GOs, category='MF')
myIntersetingGenes <- example@nodes$node
geneList <- factor(as.integer(gene_universe%in%myIntersetingGenes))
names(geneList) <- gene_universe

Godata <- new("topGOdata", ontology = "BP", allGenes = geneList,
              annot = annFUN.gene2GO, gene2GO =  GOs)