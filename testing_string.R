library(STRINGdb)
string_db <- STRINGdb$new(version = "11", species = 9606, input_directory = "")
human_string_network <- string_db$get_graph()
data(diff_exp_example1)

string_db$map(diff_exp_example1, "gene", removeUnmappedRows = TRUE)

example1_mapped <- string_db$map(diff_exp_example1, "gene", removeUnmappedRows = TRUE )