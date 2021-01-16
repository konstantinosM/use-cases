devtools::install_github('PavlidisLab/gemmaAPI.R')
library("gemmaAPI")

data <- datasetInfo('GSE1297', request='data', filter = FALSE)
design <- datasetInfo('GSE1297', request='design')

