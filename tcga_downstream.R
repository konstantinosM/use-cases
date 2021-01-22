devtools::load_all("../keypathwayminer-R/")
files <- list.files(path = "use_case_data/gdc_data/kpm_results/")
results <- list()
for (file in (files)) {
  results <- c(results, readRDS(file = paste0("use_case_data/gdc_data/kpm_results/", file))@configurations)
}
# Save pathways from several runs in a Result object
result <- new("Result", configurations = results, parameters = list(computed_pathways = 20))
# Compute pathway statistics
counts_matrix_z_4 <- readRDS("use_case_data/gdc_data/data/counts_matrix_z_4.rds")
result <- pathway_statistics(indicator_matrix = counts_matrix_z_4, result = result)
# Visualize the results with shiny
visualize_result(result)


plot_union_network_comparison <- function(result) {
  configurations <- get_configurations(result)
  config <- c()
  numNodes <- c()
  avgDiffExp <- c()
  for (configuration in configurations) {
    pathways <- get_pathways(result_object = result, configuration = configuration)
    union_network <- pathways@union_network
    config <- c(config, configuration)
    numNodes <- c(numNodes, union_network@num_nodes)
    avgDiffExp <- c(avgDiffExp, union_network@avg_exp)
  }
  return(tibble(config=config,numNodes=numNodes,avgDiffExp=avgDiffExp))
}






# install.packages("ggplot2")
# load package and data
options(scipen = 999) # turn-off scientific notation like 1e+48
library(ggplot2)
theme_set(theme_bw()) # pre-set the bw theme.
data("midwest", package = "ggplot2")
# midwest <- read.csv("http://goo.gl/G1K41K")  # bkup data source

# Scatterplot
gg <- ggplot(midwest, aes(x = area, y = poptotal)) +
  geom_point(aes(col = state, size = popdensity)) +
  geom_smooth(method = "loess", se = F) +
  xlim(c(0, 0.1)) +
  ylim(c(0, 500000)) +
  labs(
    subtitle = "Area Vs Population",
    y = "Population",
    x = "Area",
    title = "Scatterplot",
    caption = "Source: midwest"
  )

plot(gg)
