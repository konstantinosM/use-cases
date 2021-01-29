devtools::load_all("../keypathwayminer-R/")
.jaddClassPath("../keypathwayminer-R/inst/java/keypathwayminer-standalone-5.0.jar")
# Step 1: Get KPM results ####
ines_result <- readRDS(file = "use_case_data/geo_data/kpm_results/lfc1_p0,001/INES/innes_greedy_results_sars_cov_2_lfc1_p0001.rds")
glone_result <- readRDS(file = "use_case_data/geo_data/kpm_results/lfc1_p0,001/GLONE/glone_greedy_results_sars_cov_2_lfc1_p0001.rds")
indicator_matrix <- readRDS("use_case_data/geo_data/data/indicator_matrix_lfc1_p0,001.rds")
# Step 2: Compute pathway statistics ####
ines_result <- pathway_statistics(indicator_matrix = indicator_matrix, result = ines_result)
glone_result <- pathway_statistics(indicator_matrix = indicator_matrix, result = glone_result)
# Step 3: Visualize and browse results with shiny####
visualize_result(ines_result)
visualize_result(glone_result)
# Step 4.1: Compare netwroks INES ####
comparison_ines <- pathway_comparison_plots(ines_result)
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

# Step 4.2: Compare networks GLONE ## ##
comparison_glone <- pathway_comparison_plots(glone_result)
# Union network comparison data
union_network_comparison_data <- comparison_glone$union_network_comparison$data
# Select best pathways
df_select <- union_network_comparison_data[union_network_comparison_data$avgDiffExp > 6 & union_network_comparison_data$numNodes > 30, ]
selected_glone_union <- df_select
# Encircle pathways in plot
union_network_comparison_plot_glone <- comparison_glone$union_network_comparison$plot +
  geom_encircle(data = df_select, aes(x = numNodes, y = avgDiffExp), color = "red", spread = 0.001) +
  geom_text(data = df_select, aes(label = config), hjust = -0.2, vjust = 0) +
  labs(caption = "Encircled configurations with avgDiffExp > 6 and numNodes > 30") 

# Save plot
#ggsave(filename = "~/Desktop/plots/geo/LFC1_P0,001/network_comparison/union_network_comparison_glone.png", union_network_comparison_plot_glone, width = 14, height = 8)

# Top pathways comparison
top_pathway_comparison_data <- comparison_glone$top_pathway_comparison$data
# Select best pathways
df_select <- top_pathway_comparison_data[top_pathway_comparison_data$avgDiffExp > 7.5 & top_pathway_comparison_data$numNodes > 10,]
selected_glone_top <- df_select
# Encircle pathways in plot
top_pathway_comparison_plot_golne <- comparison_glone$top_pathway_comparison$plot +
  geom_encircle(data = df_select, aes(x = numNodes, y = avgDiffExp), color = "red", spread = 0.001) +
  geom_text(data = df_select, aes(label = config), hjust = -0.2, vjust = 0)+ 
  labs(caption = "Encircled configurations with avg_diff_exp > 7.5 and num_nodes > 10")

# Save plot
#ggsave(filename = "~/Desktop/plots/geo/LFC1_P0,001/network_comparison/top_pathway_comparison_glone.png", top_pathway_comparison_plot_golne, width = 14, height = 8)

# Step 4.3 Create grid with all comparison plots ####
plots <- list(union_network_comparison_plot_ines, union_network_comparison_plot_glone, top_pathway_comparison_plot_ines, top_pathway_comparison_plot_golne)
grid <- ggarrange(plotlist = plots,nrow = 2, ncol = 2)
grid <- annotate_figure(grid, left = "INES", right = "GLONE")
#ggsave(filename = "~/Desktop/plots/geo/LFC1_P0,001/network_comparison/grid.png", grid, width = 25, height = 12)

# Step 5 Save pathways that should be analysed further ####
# INES results
ines_pathways <- list()
for(configuration_name in selected_ines_top$config){
  ines_pathways[configuration_name] <- get_pathway(configuration = get_configuration(ines_result, configuration_name), pathway_name = "Pathway-1")
}
for(configuration_name in selected_ines_union$config){
  ines_pathways[paste0(configuration_name,"-union")] <- get_pathway(configuration = get_configuration(ines_result, configuration_name), union = TRUE)
}
saveRDS(ines_pathways, "use_case_data/geo_data/comparison/selected_ines_pathways.rds")
# GLONE results 
glone_pathways <- list()
for(configuration_name in selected_glone_top$config){
  glone_pathways[configuration_name] <- get_pathway(configuration = get_configuration(glone_result, configuration_name), pathway_name = "Pathway-1")
}
for(configuration_name in selected_glone_union$config){
  glone_pathways[paste0(configuration_name,"-union")] <- get_pathway(configuration = get_configuration(glone_result, configuration_name), union = TRUE)
}
saveRDS(glone_pathways, "use_case_data/geo_data/comparison/selected_glone_pathways.rds")



