devtools::load_all("../keypathwayminer-R/")
.jaddClassPath("../keypathwayminer-R/inst/java/keypathwayminer-standalone-5.0.jar")
# Read saved files
human_biogrid_network <- readRDS("use_case_data/geo_data/graphs/human_biogrid_network.rds")
counts_matrix_lfc1_p0001 <- readRDS("use_case_data/geo_data/data/indicator_matrix_lfc1_p0,001.rds")

options(java.parameters = "-Xmx64000m")
# Greedy INES run
kpm_options(
  execution = "Local",
  strategy = "INES",
  algorithm = "Greedy",
  use_range_l = TRUE,
  use_range_k = TRUE,
  l_min = 1,
  l_step = 1,
  l_max = 10,
  k_min = 1,
  k_step = 1,
  k_max = 10
)

innes_results_sars_cov_2_lfc1_p0001 <- kpm(graph = human_biogrid_network, indicator_matrices = counts_matrix_lfc1_p0001)
saveRDS(innes_results_sars_cov_2_lfc1_p0001, "use_case_data/geo_data/kpm_results/lfc1_p0,001/INES/innes_results_sars_cov_2_lfc1_p0001.rds")
