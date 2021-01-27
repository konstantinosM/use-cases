devtools::load_all("../keypathwayminer-R/")
.jaddClassPath("../keypathwayminer-R/inst/java/keypathwayminer-standalone-5.0.jar")
# Step 1: Read saved files ####
human_biogrid_network <- readRDS("use_case_data/geo_data/graphs/human_biogrid_network.rds")
counts_matrix_lfc1_p0001 <- readRDS("use_case_data/geo_data/data/indicator_matrix_lfc1_p0,001.rds")

# Step 2: Run KPM INES and save results ####
options(java.parameters = "-Xmx64000m")
result_list <- list()
# Greedy INES run
for (l in c(1:10)) {
  for (k in c(1:10)) {
    kpm_options(
      execution = "Local",
      strategy = "INES",
      algorithm = "Greedy",
      l_min = as.numeric(l),
      k_min = as.numeric(k)
    )
    result_list <- c(result_list, kpm(graph = human_biogrid_network, indicator_matrices = counts_matrix_lfc1_p0001)@configurations)
  }
}

# Save pathways from several runs in a Result object
result <- new("Result", configurations = result_list, parameters = list(computed_pathways = kpm_options()$computed_pathways))

# Save result
saveRDS(result, "use_case_data/geo_data/kpm_results/lfc1_p0,001/INES/innes_greedy_results_sars_cov_2_lfc1_p0001.rds")

# Step 3: Run KPM GLONE and save results ####
options(java.parameters = "-Xmx64000m")
result_list <- list()
# Greedy GLONE run
for (l in seq(20, 400, by= 20)) {
    kpm_options(
      execution = "Local",
      strategy = "GLONE",
      algorithm = "Greedy",
      l_min = as.numeric(l)
    )
    result_list <- c(result_list, kpm(graph = human_biogrid_network, indicator_matrices = counts_matrix_lfc1_p0001)@configurations)
}
# Save pathways from several runs in a Result object
result <- new("Result", configurations = result_list, parameters = list(computed_pathways = kpm_options()$computed_pathways))

# Save result
saveRDS(result, "use_case_data/geo_data/kpm_results/lfc1_p0,001/GLONE/glone_greedy_results_sars_cov_2_lfc1_p0001.rds")
