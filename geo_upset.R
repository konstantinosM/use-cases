library("UpSetR")
indicator_matrix <- readRDS("use_case_data/geo_data/data/indicator_matrix_p0,001_lfc0,5.rds")
message("Running upset")
upset_plot <- upset(indicator_matrix,
      sets = colnames(indicator_matrix)[-1],
      sets.bar.color = "#56B4E9",
      order.by = "freq",
      nintersects = 20
)

pdf(file="filename.pdf", width = 20)
upset_plot
dev.off()