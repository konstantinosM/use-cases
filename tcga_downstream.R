files <- list.files(path = "use_case_data/gdc_data/kpm_results/")
results <- list()
for(file in (files)){
  results <- c(results, readRDS(file = paste0("use_case_data/gdc_data/kpm_results/", file))@configurations)
}
result <- new('Result',configurations = results)
# Visualize the results with shiny
visualize_result(resutls)


pathway_statistics <- function(indicator_matrix, result){
  configurations <- get_configurations(result)
  for(configuration in configurations){
    pathways <- get_pathways(result_object = result, configuration = configuration)
    # Pathways
    for(i in pathways@pathways){
      pathways["i"]
    }
    
    # Union network
    nodes <- pathways@union_network@nodes$node 
    pathway_matrix <- indicator_matrix[indicator_matrix$id%in%nodes,]
    pathways@union_network<- set_avg_exp(pathway =  pathways@union_network, new_avg_exp = sum(rowSums(pathway_matrix[-1]))/nrow(pathway_matrix))
    
  }
  bert <- get_pathways(result_object = result,configuration = ha)
  bert@union_network
  # bento =  configurations[configuration][1]
    # # Union network
    # configuration@union_network@nodes$node
    # configuration@union_network@avg_exp
    # indicator_matrix[indicator_matrix$id=configuration@union_network@nodes$node]
    
  }
