library(soil.spec)
library(reshape)


source("1-Output_soil.R")

aggregate_soil()

source("4_predict_soil.R")

predict_properties()
