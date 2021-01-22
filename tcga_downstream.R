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
# visualize_result(result)

df <- plot_union_network_comparison(result = result)
df_select <- df[df$avgDiffExp > 30 & df$numNodes > 15, ]
# Scatterplot
ggplot(df, aes(x = numNodes, y = avgDiffExp)) +
  geom_point(aes(col = config, size = avgDiffExp)) +
  geom_encircle(data = df_select, aes(x = numNodes, y = avgDiffExp), color = "red", spread = 0.001) +
  geom_text(data = df_select, aes(label = config), hjust = -0.2, vjust = 0) +
  scale_x_continuous(breaks = round(seq(0, max(df$numNodes), by = 10), 1)) +
  scale_y_continuous(breaks = round(seq(0, max(df$avgDiffExp), by = 20), 1)) +
  labs(
    subtitle = "Number of nodes Vs Avg. DE. cases per gene",
    y = "Average differentialy expressed cases per gene",
    x = "Nodes in the pathway",
    title = "Union network comparison",
    caption = "Encircled configurations with avg_diff_exp > 30 and num_nodes > 15"
  )
