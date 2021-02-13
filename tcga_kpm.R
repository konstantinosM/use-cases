library("KeyPathwayMineR")
# Step 1: Read saved files ####
human_string_network <- readRDS("use_case_data/tcga_data/graphs/human_string_network_800.rds")
counts_matrix_z_2 <- readRDS("use_case_data/tcga_data/data/counts_matrix_z_2.rds")
counts_matrix_z_3 <- readRDS("use_case_data/tcga_data/data/counts_matrix_z_3.rds")
counts_matrix_z_4 <- readRDS("use_case_data/tcga_data/data/counts_matrix_z_4.rds")

# Step 2: Run KPM INES and save results for z-score = 2 ####
# Since the datasets and biological network are very big we need to allocate some additional memory
options(java.parameters = "-Xmx64000m")
result_list <- list()

# Greedy INES run
for (l in c(150, 200, 250)) {
  for (k in c(2:20)) {
    kpm_options(
      execution = "Local",
      strategy = "INES",
      algorithm = "Greedy",
      l_min = l,
      k_min = k
    )
    result_list <- c(result_list,  kpm(graph = human_string_network, indicator_matrices = counts_matrix_z_2)@configurations)
  }
}
# Save pathways from several runs in a Result object
result <- new("Result", configurations = result_list, parameters = list(computed_pathways = kpm_options()$computed_pathways, strategy = kpm_options()$strategy))

# Save result
saveRDS(result, "use_case_data/tcga_data/kpm_results/ines_greedy_results_prostate_tumor_z_score_2.rds")

# Step 3: Run KPM INES and save results for z-score = 3 ####
# Since the datasets and biological network are very big we need to allocate some additional memory
options(java.parameters = "-Xmx64000m")
result_list <- list()

# Greedy INES run
for (l in c(150, 200, 250)) {
  for (k in c(2:20)) {
    kpm_options(
      execution = "Local",
      strategy = "INES",
      algorithm = "Greedy",
      l_min = l,
      k_min = k
    )
    result_list <- c(result_list,  kpm(graph = human_string_network, indicator_matrices = counts_matrix_z_3)@configurations)
  }
}
# Save pathways from several runs in a Result object
result <- new("Result", configurations = result_list, parameters = list(computed_pathways = kpm_options()$computed_pathways, strategy = kpm_options()$strategy))

# Save result
saveRDS(result, "use_case_data/tcga_data/kpm_results/ines_greedy_results_prostate_tumor_z_score_3.rds")

# Step 4: Run KPM INES and save results for z-score = 4 ####
# Since the datasets and biological network are very big we need to allocate some additional memory
options(java.parameters = "-Xmx64000m")
result_list <- list()

# Greedy INES run
for (l in c(150, 200, 250)) {
  for (k in c(2:20)) {
    kpm_options(
      execution = "Local",
      strategy = "INES",
      algorithm = "Greedy",
      l_min = l,
      k_min = k
    )
    result_list <- c(result_list,  kpm(graph = human_string_network, indicator_matrices = counts_matrix_z_4)@configurations)
  }
}

# Save pathways from several runs in a Result object
result <- new("Result", configurations = result_list, parameters = list(computed_pathways = kpm_options()$computed_pathways, strategy = kpm_options()$strategy))

# Save result
saveRDS(result, "use_case_data/tcga_data/kpm_results/ines_greedy_results_prostate_tumor_z_score_4.rds")
