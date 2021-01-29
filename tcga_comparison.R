devtools::load_all("../keypathwayminer-R/")
.jaddClassPath("../keypathwayminer-R/inst/java/keypathwayminer-standalone-5.0.jar")
# Step 1: Get KPM results ####
# Step 1.1: Z-score = 3 ####
z_score_3_results <- readRDS(file = "use_case_data/gdc_data/kpm_results/ines_greedy_results_prostate_tumor_z_score_3.rds")
z_score_3_indicator_matrix <- readRDS(file = "use_case_data/gdc_data/data/counts_matrix_z_3.rds") 
# Step 1.2: Z-score = 4 ####
z_score_4_results <- readRDS(file = "use_case_data/gdc_data/kpm_results/ines_greedy_results_prostate_tumor_z_score_4.rds")
z_score_4_indicator_matrix <- readRDS(file = "use_case_data/gdc_data/data/counts_matrix_z_4.rds") 

# Step 2: Compute pathway statistics ####
# Step 2.1: Z-score = 3 ####
z_score_3_results <- pathway_statistics(indicator_matrix = z_score_3_indicator_matrix, result = z_score_3_results)
# Step 2.2: Z-score = 4 ####
z_score_4_results <- pathway_statistics(indicator_matrix = z_score_4_indicator_matrix, result = z_score_4_results)

# Step 3: Visualize and browse results with shiny####
# Step 3.1: Z-score = 3 ####
visualize_result(z_score_3_results)
# Step 3.2: Z-score = 4 ####
visualize_result(z_score_4_results)

# Step 4: Compare union networks, and top networks of every configuration for best results ####
# Step 4.1: Z-score = 3 ####
comparison_z_score_3 <- pathway_comparison_plots(z_score_3_results)
# Union network comparison
union_network_comparison_data <- comparison_ines$union_network_comparison$data
# Select best pathways
df_select <- union_network_comparison_data[union_network_comparison_data$avgDiffExp > 6 & union_network_comparison_data$numNodes > 6, ]
df_select_2 <- union_network_comparison_data[union_network_comparison_data$avgDiffExp > 5 & union_network_comparison_data$numNodes > 40, ]
selected_ines_union <-rbind(df_select, df_select_2)
# Encircle pathways in plot
union_network_comparison_plot_ines <- comparison_ines$union_network_comparison$plot +
  geom_encircle(data = df_select, aes(x = numNodes, y = avgDiffExp), color = "red", spread = 0.001) +
  geom_encircle(data = df_select_2, aes(x = numNodes, y = avgDiffExp), color = "red", spread = 0.001) +
  geom_text(data = df_select, aes(label = config), hjust = -0.2, vjust = 0) +
  geom_text(data = df_select_2, aes(label = config), hjust = -0.2, vjust = 0)+
  labs(caption = "Encircled configurations with avgDiffExp > 6 and numNodes > 6, and avgDiffExp > 6 and numNodes > 40") 

# Save plot
#ggsave(filename = "~/Desktop/plots/geo/LFC1_P0,001/network_comparison/union_network_comparison_ines.png", union_network_comparison_plot_ines, width = 14, height = 8)

# Top pathways comparison
top_pathway_comparison_data <- comparison_ines$top_pathway_comparison$data
# Select best pathways
df_select <- top_pathway_comparison_data[top_pathway_comparison_data$avgDiffExp > 8 & top_pathway_comparison_data$numNodes > 5, ]
selected_ines_top <- df_select
# Encircle pathways in plot
top_pathway_comparison_plot_ines <- comparison_ines$top_pathway_comparison$plot +
  geom_encircle(data = df_select, aes(x = numNodes, y = avgDiffExp), color = "red", spread = 0.001) +
  geom_text(data = df_select, aes(label = config), hjust = -0.2, vjust = 0)+ 
  labs(caption = "Encircled configurations with avg_diff_exp > 8 and num_nodes > 5")

# Save plot
#ggsave(filename = "~/Desktop/plots/geo/LFC1_P0,001/network_comparison/top_pathway_comparison_ines.png", top_pathway_comparison_plot_ines, width = 14, height = 8)

# Step 4.2: Z-score = 4 ####

























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
