devtools::load_all("../keypathwayminer-R/")
.jaddClassPath("../keypathwayminer-R/inst/java/keypathwayminer-standalone-5.0.jar")

# Step 1: Get KPM results ####
z_score_3_results <- readRDS(file = "use_case_data/geo_data/kpm_results/lfc1_p0,001/INES/innes_greedy_results_sars_cov_2_lfc1_p0001.rds")
z_score_3_indicator_matrix <- 
  
glone_result <- readRDS(file = "use_case_data/geo_data/kpm_results/lfc1_p0,001/GLONE/glone_greedy_results_sars_cov_2_lfc1_p0001.rds")
indicator_matrix <- readRDS("use_case_data/geo_data/data/indicator_matrix_lfc1_p0,001.rds")

# Step 1: Get results and computepathway statistics ####
files <- list.files(path = "use_case_data/gdc_data/kpm_results/z_score_3/")
results <- list()
for (file in (files)) {
  results <- c(results, readRDS(file = paste0("use_case_data/gdc_data/kpm_results/z_score_3/", file))@configurations)
}
# Save pathways from several runs in a Result object
result <- new("Result", configurations = results, parameters = list(computed_pathways = 20))
# Compute pathway statistics
counts_matrix_z_4 <- readRDS("use_case_data/gdc_data/data/counts_matrix_z_4.rds")
result <- pathway_statistics(indicator_matrix = counts_matrix_z_4, result = result)
# Step 2: Visualize and compare results ####
# Shiny
visualize_result(result)

df <- plot_union_network_comparison(result = result)
df_select <- df[df$avgDiffExp > 30 & df$numNodes > 15, ]
# Scatterplot
union_network_comparison # <- ggplot(df, aes(x = numNodes, y = avgDiffExp)) +
geom_point(aes(col = config, size = avgDiffExp)) +
  geom_encircle(data = df_select, aes(x = numNodes, y = avgDiffExp), color = "red", spread = 0.001) +
  geom_text(data = df_select, aes(label = config), hjust = -0.2, vjust = 0) +
  labs(
    subtitle = "Number of nodes Vs Avg. DE. cases per gene",
    y = "Average differentialy expressed cases per gene",
    x = "Nodes in the pathway",
    title = "Union network comparison",
    caption = "Encircled configurations with avg_diff_exp > 30 and num_nodes > 15",
    col = "Configurations",
    size = "Avg. DE cases per configuration"
  ) + theme(legend.position = "right", plot.title = element_text(size = 20), text = element_text(size = 15))

ggsave(filename = "~/Desktop/union_network_comparison.png", union_network_comparison, width = 14, height = 8)

# Step 3: Visualize top networks ####
top1 <- get_pathway(configuration = result@configurations[df_select$config[1]]$`K-4-L1-250`, union = TRUE)
top2 <- get_pathway(configuration = result@configurations[df_select$config[2]]$`K-6-L1-250`, union = TRUE)

top1_nodes <- top1@nodes$node
string_db$plot_network(top1_nodes)
top1_ENRICHMENT <- string_db$get_enrichment(top1_nodes)

top2_nodes <- top2@nodes$node
string_db$plot_network(top2_nodes)
