source("install_and_load_libraries.R")
string_db <- STRINGdb$new(
  version = "11",
  species = 9606,
  input_directory = "",
  score_threshold = 600
)

human_string_network <- string_db$get_graph()
human_string_network <- readRDS("use_case_data/gdc_data/data/counts_matrix_z_4.rds")

# Greedy INES run
kpm_options(
  execution = "Local",
  strategy = "INES",
  algorithm = "Greedy",
  use_range_k = TRUE,
  l_min = 250,
  k_min = 4,
  k_step = 2,
  k_max = 20
)

innes_results_prostate_tumor_z_4 <- kpm(graph = human_string_network, indicator_matrices = counts_matrix_z_4)
saveRDS(innes_results_prostate_tumor_z_4, "use_case_data/gdc_data/kpm_results/prostate_tumor_z_score_4_innes.rds")
