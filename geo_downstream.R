devtools::load_all("../keypathwayminer-R/")
.jaddClassPath("../keypathwayminer-R/inst/java/keypathwayminer-standalone-5.0.jar")
# Step 1: Get results ####
ines_result <- readRDS(file ="use_case_data/geo_data/kpm_results/lfc1_p0,001/INES/innes_greedy_results_sars_cov_2_lfc1_p0001.rds")
glone_result <- readRDS(file ="use_case_data/geo_data/kpm_results/lfc1_p0,001/GLONE/glone_greedy_results_sars_cov_2_lfc1_p0001.rds")
indicator_matrix <-  readRDS("use_case_data/geo_data/data/indicator_matrix_lfc1_p0,001.rds")
# Step 2: Compute pathway statistics ####
ines_result <- pathway_statistics(indicator_matrix = indicator_matrix, result = ines_result)
glone_result <- pathway_statistics(indicator_matrix = indicator_matrix, result = glone_result)

# Step 3: Visualize and browse results with shiny####
visualize_result(ines_result)
visualize_result(glone_result)

# Step 4.1: Compare netwroks INES ####
comparison_ines <- pathway_comparison_plots(ines_result)

# Union network comparison
union_network_comparison_data <- comparison_plots_ines$union_network_comparison$data
# Select best pathways
df_select <- union_network_comparison_data[union_network_comparison_data$avgDiffExp > 6 & union_network_comparison_data$numNodes > 6, ]
df_select_2 <- union_network_comparison_data[union_network_comparison_data$avgDiffExp > 5 & union_network_comparison_data$numNodes > 40, ]

union_network_comparison_plot <- comparison_plots_ines$union_network_comparison$plot +
  geom_encircle(data = df_select, aes(x = numNodes, y = avgDiffExp), color = "red", spread = 0.001) +
  geom_encircle(data = df_select_2, aes(x = numNodes, y = avgDiffExp), color = "red", spread = 0.001) +
  labs(Caption = "Encircled configurations with avg_diff_exp > 30 and num_nodes > 15") +
  geom_text(data = df_select, aes(label = config), hjust = -0.2, vjust = 0) +
  geom_text(data = df_select_2, aes(label = config), hjust = -0.2, vjust = 0) 

ggsave(filename = "~/Desktop/plots/geo/network_comparison/lfc1_p0,001/union_network_comparison_ines.png", union_network_comparison_plot, width = 14, height = 8)

# Top pathways comparison
top_pathway_comparison_data <- comparison_plots_ines$top_pathway_comparison$data


top_pathway_comparison_plot




# Step 4.1: Compare networks GLONE ####



ggsave(filename = "~/Desktop/union_network_comparison3.png", union_network_comparison,width = 14,height = 8)

union_network_comparison_data <- comparison_plots_ines$union_network_comparison
top_pathway_comparison_plot <- comparison_plots_ines$top_pathway_comparison


# Scatterplot

