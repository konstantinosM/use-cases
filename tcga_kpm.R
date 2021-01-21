devtools::load_all("../keypathwayminer-R/")
.jaddClassPath("../keypathwayminer-R/inst/java/keypathwayminer-standalone-5.0.jar")
# Read saved files
human_string_network <- readRDS("use_case_data/gdc_data/graphs/human_string_network_800.rds")
counts_matrix_z_4 <- readRDS("use_case_data/gdc_data/data/counts_matrix_z_4.rds")
# Since the datasets and biological network are very big we need to allocate some additional memory
options(java.parameters = "-Xmx64000m")
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
    innes_results_prostate_tumor_z_4 <- kpm(graph = human_string_network, indicator_matrices = counts_matrix_z_4)
    saveRDS(innes_results_prostate_tumor_z_4, paste("use_case_data/gdc_data/kpm_results/prostate_tumor__z_score_4_innes_K=", k, "_L=",l,".rds", sep = ""))
  }
}
