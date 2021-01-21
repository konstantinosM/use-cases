files <- list.files(path = "use_case_data/gdc_data/kpm_results/")
results <- list()
for (file in (files)) {
  results <- c(results, readRDS(file = paste0("use_case_data/gdc_data/kpm_results/", file))@configurations)
}
result <- new("Result", configurations = results)
# Visualize the results with shiny
visualize_result(resutls)

# For the division we use length(nodes) and not nrow(pathway_matrix)
# since not all the genes in the pathway are also in the indicator matrix.
pathway_statistics <- function(indicator_matrix, result) {
  configurations <- get_configurations(result)
  # For all configurations in the result object
  for (configuration in configurations) {
    pathways <- get_pathways(result_object = result, configuration = configuration)
    # Pathways
    for (pathway in names(pathways@pathways)) {
      pathway <- get_pathway(configuration = pathways, pathway_name = pathway)
      nodes <- pathway@nodes$node
      pathway_matrix <- indicator_matrix[indicator_matrix$id %in% nodes, ]
      pathway <- set_avg_exp(pathway = pathway, new_avg_exp = sum(rowSums(pathway_matrix[-1])) / length(nodes))
    }
    # Union network
    nodes <- pathways@union_network@nodes$node
    pathway_matrix <- indicator_matrix[indicator_matrix$id %in% nodes, ]
    pathways@union_network <- set_avg_exp(pathway = pathways@union_network, new_avg_exp = sum(rowSums(pathway_matrix[-1])) / length(nodes))
  }
}
