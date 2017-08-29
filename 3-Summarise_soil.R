#summarise results from long runs

#will run on different holdout samps to get a true estimate of accuracy with increasing samp size
#start by setting WD to here
summarise_runs <- function(numberOfRuns){
  setwd("~/oaf/soil_lab/projects/sample_scripts")
  wd <- getwd()
  
  #set methods
  method = c("PLSM","SVM","XGBM")
  
  listwd <- c()
  
  #dirin <- choose.dir(default="D:/OneAcre/Google Drive/One Acre Fund/OAF Soil Lab Folder/Projects",caption = "Choose prediction sub-directory")
  dirin <- tk_choose.dir(default="~/oaf/soil_lab/projects/", caption = "Choose prediction sub-directory")
  #dirin <- "D:/OneAcre/Google Drive/One Acre Fund/OAF Soil Lab Folder/Projects/my_feb_17/4_predicted"
  ##dirin <- "D:\OneAcre\Google Drive\One Acre Fund\OAF Soil Lab Folder\Projects\my_feb_17\4_predicted"
  #get soil props
  newwd <- paste(dirin, "PLSM",1 , sep="/")
  setwd(newwd)
  tin <- read.csv("Model_Summary.csv")
  models <- data.frame("soil property" = tin[1])
  
  for(xm in method ){
    for(nm in seq(1:numberOfRuns)){
      newwd <- paste(dirin, xm ,sep="/")
      newwd <- paste(newwd, nm, sep="/")
      setwd(newwd)
      te <- read.csv("Model_Summary.csv")
      p <- paste(xm, nm, sep = "_")
      df <- data.frame("run" = te$Holdout.Rsquared)
      colnames(df) <- p
      models <- cbind(models, df)
    }
  }
  
  #remember avs is changing the holdout pc
  #fc is changing the actual number of samples run (constant 30% holdout)
  models$meanPLSM <- round(rowMeans(models[grep("PLSM", names(models))]), 2)
  models$sdPLSM   <- round(apply(models[grep("PLSM", names(models))], 1, sd), 2)
  
  models$meanSVM  <- round(rowMeans(models[grep("SVM", names(models))]), 2)
  models$sdSVM    <- round(apply(models[grep("SVM", names(models))], 1, sd), 2)
  
  models$meanXGBM <- round(rowMeans(models[grep("XGBM", names(models))]), 2)
  models$sdXGBM   <- round(apply(models[grep("XGBM", names(models))], 1, sd), 2)
  
  #check if passes R2 and SD test
  models_sub1 <- subset(models, meanSVM  > 0.5 & sdSVM  < 0.1)
  models_sub2 <- subset(models, meanPLSM > 0.5 & sdPLSM < 0.1)
  models_sub3 <- subset(models, meanXGBM > 0.5 & sdXGBM < 0.1)
  
  models_sub <- rbind(models_sub1, models_sub2, models_sub3)
  models_sub <- unique(models_sub)
  
  #find max model - Should probably add sd below as well[Max]
  SVM   <- intersect(which(models_sub$meanSVM >= models_sub$meanPLSM), which(models_sub$meanSVM >= models_sub$meanXGBM))
  PLSM  <- intersect(which(models_sub$meanPLSM > models_sub$meanSVM),  which(models_sub$meanPLSM >= models_sub$meanXGBM))
  XGBM  <- intersect(which(models_sub$meanXGBM > models_sub$meanSVM),  which(models_sub$meanXGBM > models_sub$meanPLSM))
  
  #new df
  models_best <- rbind(models_sub[PLSM, ], models_sub[SVM, ], models_sub[XGBM, ])
  
  #return soil prop names of passing feats
  PLSM_vars <- models_sub[PLSM, ][1]
  SVM_vars  <- models_sub[SVM, ][1]
  XGBM_vars <- models_sub[XGBM, ][1]
  
  PLSM_vars <- as.character(PLSM_vars[,1])
  SVM_vars  <- as.character(SVM_vars[,1])
  XGBM_vars <- as.character(XGBM_vars[,1])
  
  #read in prediction files to merge
  #Initialize result data frame
  fi=paste(dirin, "PLSM", 1, sep="/")
  setwd( fi )
  r <- read.csv("All_predictions.csv")
  r <- r[-1]
  nr <- dim(r)[1]
  nc <- dim(r)[2]
  M_in <- data.frame(matrix(0, ncol = nc, nrow = nr))
  
  PLSM_in <- model_in(model = "PLSM", numberOfRuns = numberOfRuns, M_in = M_in, dirin = dirin)
  SVM_in  <- model_in(model = "SVM",  numberOfRuns = numberOfRuns, M_in = M_in, dirin = dirin)
  XGBM_in <- model_in(model = "XGBM",  numberOfRuns = numberOfRuns, M_in = M_in, dirin = dirin)
  
  ##Okay - now take the best cols from each
  PLSMind <- match(PLSM_vars, names(PLSM_in))
  SVMind  <- match(SVM_vars,  names(SVM_in))
  XGBMind <- match(XGBM_vars, names(XGBM_in))
  
  ###
  newdf <- cbind(PLSM_in[PLSMind], SVM_in[SVMind], XGBM_in[XGBMind])
  
  rs <- NULL
  for(i in seq(1:dim(models_best)[1])) {
    su <- models_best[i, ]
    n <- max(su$meanPLSM, su$meanSVM, su$meanXGBM)
    rs <- c(rs, n)
  }
  
  rsq <- data.frame("new" = c( rs , "mean_rsqaured") )
  newdf <- cbind(newdf, SVM_in$SSN)
  
  fo <- paste(dirin, "best-predictions.csv",sep="/")
  write.csv(newdf, fo)
  
  #slim models best
  models_slim <- models_best
  new_slim <- data.frame("Soil Property" = models_slim[1], "Best R2" =1, "Best SD" = 1)
  for(row in 1:dim(models_slim)[1] ){
    su <- models_slim[row,]
    a <- max(su$meanXGBM, su$meanPLSM, su$meanSVM)
    new_slim[row, 2] <- a
    
    if((su$meanSVM >= su$meanPLSM) | (su$meanSVM >= su$meanXGBM)){
      b <- su$sdSVM
    }
    else if((su$meanSVM < su$meanPLSM) | (su$meanPLSM >= su$meanXGBM)){
      b <- su$sdPLSM
    }
    else{
      b <- su$sdXGBM
    }
    new_slim[row, 3] <- b 
  }
  
  models_ <- models #don't quite understand the logic behind the for-loop above & below. ICRAF is taking max Rsq from one model and min sd from a different model. Will change this later[Max]
  new_ <- data.frame("Soil Property" = models_[1], "Best R2" = 1, "Best SD" = 1) 
  for(row in 1:dim(models_)[1] ){
    s <- models_[row, ]
    a <- max(s$meanSVM, s$meanPLS, s$meanXGBM)
    new_[row, 2] <- a
    if((su$meanSVM >= su$meanPLSM) | (su$meanSVM >= su$meanXGBM)){
      b <- su$sdSVM
    }
    else if((su$meanSVM < su$meanPLSM) | (su$meanPLSM >= su$meanXGBM)){
      b <- su$sdPLS
    }
    else{
      b <- su$sdXGBM
    }
    new_[row, 3] <- b 
  }
  
  foo <- paste(dirin, "Summary_of_best_models.csv", sep="/")
  write.csv(new_slim, foo, row.names = FALSE)
  
  SSNs <- SVM_in$SSN
  #SSNs
  print(paste("SAVED IN:", fo))
  
  #D:\OneAcre\Google Drive\One Acre Fund\OAF Soil Lab Folder\Projects\my_feb_17\4_predicted
  dirin2 <- paste(dirin,"other_summaries", sep="/")
  dir.create(dirin2, showWarnings = FALSE) 
  fo <- paste(dirin2, "PLSM-predictions.csv", sep="/")
  PLSM_in$SSN <- NULL
  PLSM_in <- cbind(PLSM_in, SSNs)
  write.csv(PLSM_in, fo)
  
  fo <- paste(dirin2, "SVM-predictions.csv", sep="/")
  SVM_in$SSN <- NULL
  SVM_in <- cbind(SVM_in, SSNs)
  write.csv(SVM_in, fo)
  
  fo <- paste(dirin2, "XGBM-predictions.csv", sep="/")
  XGBM_in$SSN <- NULL
  XGBM_in <- cbind(XGBM_in, SSNs)
  write.csv(XGBM_in, fo)
  
  fo <- paste(dirin2, "Full-model-summary.csv", sep="/")
  write.csv(new_, fo, row.names = FALSE)
  
  
  ###### spit out all "best" predictions, even if fail cutoff ---------
  SVM_all   <- intersect(which(models$meanSVM >= models$meanPLSM), which(models$meanSVM >= models$meanXGBM))
  PLSM_all  <- intersect(which(models$meanPLSM > models$meanSVM),  which(models$meanPLSM >= models$meanXGBM))
  XGBM_all  <- intersect(which(models$meanXGBM > models$meanSVM),  which(models$meanXGBM > models$meanPLSM))
  
  PLSM_all  <- models[PLSM_all,][1]
  SVM_all   <- models[SVM_all,][1]
  XGBM_all  <- models[XGBM_all,][1]
  
  SVM_all  <- SVM_all[complete.cases(SVM_all),]
  PLSM_all <- PLSM_all[complete.cases(PLSM_all),]
  XGBM_all <- XGBM_all[complete.cases(XGBM_all),]
  
  PLSM_all  <- as.character(PLSM_all)
  SVM_all   <- as.character(SVM_all)
  XGBM_all  <- as.character(XGBM_all)
  
  #match(SVM_all, PLSM_all, XGBM_all)
  
  #pl_in #rf_in are average predictions for each method
  all_out <- cbind(PLSM_in[PLSM_all], SVM_in[SVM_all], XGBM_in[XGBM_all])
  all_out <- cbind(all_out, SSNs)
  fooo <- paste(dirin2, "combined-predictions-including-bad-ones.csv",sep="/")
  write.csv(all_out, fooo)
  
  print("Done!")
}