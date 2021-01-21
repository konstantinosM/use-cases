devtools::install_github("PavlidisLab/gemmaAPI.R")
library("gemmaAPI")

differential <- datasetInfo("GSE1297", request = "differential")

degs <- datasetInfo("GSE1297", request = "data", filter = FALSE)
