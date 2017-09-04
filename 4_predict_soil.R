predict_properties <- function(numberOfRuns, ENS=FALSE, IND=TRUE){
  
  start.time <- proc.time()
  
  library(soil.spec)
  library(tcltk)
  #library(mvoutlier)
  #start by setting WD to here
  wd <- getwd()
  #outwd <- choose.dir(default="D:/OneAcre/Google Drive/One Acre Fund/OAF Soil Lab Folder/Projects",caption="please choose project directory")
  outwd <- tk_choose.dir(default="/home/max/oaf/soil_lab/projects",caption="please choose project directory")
  
  
  #ref_path <- paste(outwd,"3_wet_chem/wet_chem.csv",sep="/")
  
  ref_path <- tk_choose.files(default=paste(outwd,"/wet_chem/wet_chem.csv",sep="")
                              , caption="Please select the wet chem CSV",
                              multi = FALSE)
  raw_path <- tk_choose.files(default=paste(outwd,"/1_raw_opus/", sep=""), caption="Please select the mIR data CSV")
  
  #set method, either PLSM or RF
  method=c("PLSM","SVM","XGBM")
  
  #source helper file
  source("misc.R")
  
  raw <- NULL
  ref <- NULL
  #lets try it on OAF data
  
  #mir
  raw <- read.csv (raw_path)
  
  #wetchem
  ref <- read.csv (ref_path)
  
  #remove NA only rows
  ref <- ref[rowSums(is.na(ref)) < (dim(ref)[2]-6), ]
  
  #clean up column data
  ref <- ref[vapply(ref, function(x) length(unique(x)) > 1, logical(1L))]
  
  for(i in colnames(ref)[2:length(colnames(ref))]){
    ref[i] <- gsub("<", "", factor(unlist(ref[i])))
    ref[i] <- gsub(">", "", factor(unlist(ref[i])))
    ref[i] <- as.numeric(unlist(ref[i]))
  }
  
  ref <- Filter(function(x)!all(is.na(x)), ref)
  
  #check we have data left
  if(dim(ref)[1] < 10){
    print("ERROR, wetchem DATA INCOMPLETE")
    stop() }
  #check we have data left
  if(dim(raw)[1] < 10){
    print("ERROR, MIR DATA INCOMPLETE")
    stop() }
  
  # Rename first column in both tables
  colnames(raw) <- c("SSN", colnames(raw[,-1]))
  colnames(ref) <- c("SSN",colnames(ref[,-1]))
  
  ref$SSN <- sapply(ref$SSN, tolower)
  raw$SSN <- sapply(raw$SSN, tolower)
  
  #clean up field names
  ref$SSN <- gsub(" ", "", ref$SSN)
  ref$SSN <- gsub("_", "", ref$SSN)
  ref$SSN <- gsub("\\.0", "", ref$SSN)
  ref$SSN <- gsub("oafoaf", "oaf", ref$SSN)
  
  raw$SSN <- gsub(" ", "", raw$SSN)
  raw$SSN <- gsub("_", "", raw$SSN)
  raw$SSN <- gsub("\\.0", "", raw$SSN)
  raw$SSN <- gsub("oafoaf", "oaf", raw$SSN)
  
  ref$SSN <- as.factor(ref$SSN)
  raw$SSN <- as.factor(raw$SSN)
  
  x <- match(ref$SSN, raw$SSN)
  y <- match(raw$SSN, ref$SSN)
  xn <- !is.na(x)
  yn <- !is.na(y)
  
  for(i in colnames(ref)[2:length(colnames(ref))]){
    ref[[i]] <- as.numeric(ref[[i]])
  }
  
  #test we have matches
  if( sum(!is.na(match(raw$SSN, ref$SSN)))  < dim(ref)[1]*0.9   ){
    print("Found too few matches, something might be wrong, stopping run v2")
    stop()
  }
  
  print(paste("Found",sum(!is.na(match(raw$SSN, ref$SSN))),"matching entries"))
  
  #critical pre-processing step. Using sophisticated approach to outlier
  #detection and removal. Takes about 3 solid minutes to compute outliers.
  # outl <- sign2(raw[, -1]) #Based on the robustly sphered and normed data, robust principal components are computed
  # raw <- raw[outl$wfinal01==1, ]
  
   # ref <- subset(ref, select = -c(X.EC..Salts.,Phosphorus,Potassium,Manganese,Copper,
   #                                Aluminium,X.Sodium,Iron,X.C.E.C,X.Exchangeable.Acidity,
   #                                X.Phosphorus.Sorption.Index..PSI.,Sulphur,Boron,Zinc,
   #                                X.Exchangeable.Aluminium))
  print(names(ref))
  print(" ")
  if(IND){
    print("modeling 3 individual ML algo")
    # Begin by sourcing the scripts
    source('PLS.R')
    source('SVM.R')
    source("XGBM.R")
    
    for( xm in method){
      testing_old <- NULL
      for(i in 1:numberOfRuns){
        testing <- random.sample(ref = ref)
        while(length(testing) == length(testing_old)){# let's do our best to avoid same testing sets 
          testing <- random.sample(ref = ref)
        }
        print(paste("split vector length is ", length(testing)))
        print(paste("iteration", i, "of method", xm))
        newwd <- paste(outwd,"4_predicted", sep="/")
        dir.create(newwd, showWarnings = FALSE) 
        newwd<- paste(newwd, xm ,sep="/")
        dir.create(newwd, showWarnings = FALSE) 
        newwd<- paste(newwd, i, sep="/")
        dir.create(newwd, showWarnings = FALSE) 
        
        #make new folder
        setwd(newwd) 
        print(newwd)
        try(PLS(newwd, raw, ref, testing, method=xm))
        try(SVM(newwd, raw, ref, testing, method=xm))
        try(XGBM(newwd, raw, ref, testing, method=xm))
        testing_old <- testing
      }
    }  
    setwd(wd)
    print("summarising runs")
    source("3-Summarise_soil.R")
    summarise_runs(numberOfRuns=numberOfRuns)
  }
  if(ENS){
    source("ensemble.R")
    print("modeling an ensemble of 3 ML algo")
    newwd <- paste(outwd,"4_predicted", sep="/")
    try(ENS(newwd, raw, ref))
  }
  time.taken <- proc.time() - start.time
  print(paste("time taken is ", round(time.taken[3]/60, 2), " minutes"))
}
