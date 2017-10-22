library(soil.spec)
library(reshape)


source("1-Output_soil.R")

aggregate_soil()

source("4_predict_soil.R")

#default number of runs is set to 3. No need to change it. Increasing this
#value can potentially increase predictions accuracy but at the cost of
#computational time. DO NOT ENS to TRUE as it's still experimental
predict_properties(numberOfRuns = 3, IND=T, ENS=F) 
