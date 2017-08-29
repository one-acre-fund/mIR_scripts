# customRF <- list(type = "Regression", library = "randomForest", loop = NULL)
# customRF$parameters <- data.frame(parameter = c("mtry", "ntree"), class = rep("numeric", 2), label = c("mtry", "ntree"))
# customRF$grid <- function(x, y, len = NULL, search = "grid") {}
# customRF$fit <- function(x, y, wts, param, lev, last, weights, classProbs, ...) {
#   randomForest(x, y, mtry = param$mtry, ntree=param$ntree, ...)
# }
# customRF$predict <- function(modelFit, newdata, preProc = NULL, submodels = NULL){
#   predict(modelFit, newdata)
# }
# customRF$prob <- function(modelFit, newdata, preProc = NULL, submodels = NULL)
#   predict(modelFit, newdata, type = NULL)
# customRF$sort <- function(x) x[order(x[, 1]),]
# customRF$levels <- function(x) x$classes

HaarTransform=function(DF1,nTimes=1)#re-scale MIR data onto orthonormal basis. This is a sometimes usefule pre-processing step for MIR data
{
  w =function(k)
  {
    s1=dwt(k, filter="haar")
    return (s1@V[[1]])
  }
  Smt=DF1
  for (i in 1:nTimes)
  {
    Smt=t(apply(Smt,1,w))
  }
  return (data.frame(Smt))
}

Derivative=function(DF1, D=1) #smooth transform MIR data to reduce noises.
{
  df1=t(diff(t(DF1), differences = D))	
  return(df1)
}

# my custom error metric
evalerror <- function(preds, dtrain) {
  labels <- getinfo(dtrain, "label")
  err <- 1 - (sum((labels-preds)^2)/sum((labels-mean(labels))^2))
  return(list(metric = "Rsq", value = err))
}

model_in <- function(numberOfRuns, model, M_in,dirin){#this function takes in predictions from a model and computer their average
  for(i in seq(1:numberOfRuns)){#pull all predictions from this model
    fi <- paste(dirin, model, i, sep="/")
    setwd( fi )
    r <- read.csv("All_predictions.csv")
    nam <- r[1]
    ro <- r[-1]
    cols <- names(ro)
    print(dim(ro))
    M_in <- M_in + ro
  }
  #compute average
  M_in <- M_in/numberOfRuns
  names(M_in) <- cols
  M_in$SSN <- nam
  return(M_in)
}

random.sample <- function(ref){# this splits the ref data for testing and training sets
  newm <- 0.35
  testing <- round(runif(round(length(ref$SSN)*newm), min=0, max=length(ref$SSN)))
  testing <- unique(testing)
}
