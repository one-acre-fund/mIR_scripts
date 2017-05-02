#summarise results from long runs

library(ggplot2)
#will run on different holdout samps to get a true estimate of accuracy with increasing samp size

#start by setting WD to here
setwd("D:/OneAcre/Google Drive/One Acre Fund/OAF Soil Lab Folder/Projects/Sample_scripts")
wd <- getwd()


#set method, either PLSM or RF
#method=c("PLSM","RF")
method=c("RF","PLSM")



listwd <- c()

dirin <- choose.dir(default="D:/OneAcre/Google Drive/One Acre Fund/OAF Soil Lab Folder/Projects",caption = "Choose prediction sub-directory")
#dirin <- "D:/OneAcre/Google Drive/One Acre Fund/OAF Soil Lab Folder/Projects/my_feb_17/4_predicted"
##dirin <- "D:\OneAcre\Google Drive\One Acre Fund\OAF Soil Lab Folder\Projects\my_feb_17\4_predicted"
#get soil props
newwd <- paste(dirin, "RF",1 ,sep="/")
setwd(newwd)
tin <- read.csv("Model_Summary.csv")
models <- data.frame("soil property" = tin[1] )



for(xm in method ){
  for(nm in seq(1:3)){
    newwd <- paste(dirin, xm ,sep="/")
    newwd <- paste(newwd, nm, sep="/")
    setwd(newwd)
    te <- read.csv("Model_Summary.csv")
    p <- paste(xm,nm, sep = "_")
    df <- data.frame("run" <- te$Holdout.Rsquared)
    colnames(df) <- p
    models <- cbind(models, df)
    
    
  }}




#remember avs is changing the holdout pc
#fc is changing the actual number of samples run (constant 30% holdout)
models$meanRF <- rowMeans(models[c(4,2,3)])
models$sdRF <- apply(models[, c(4,2,3)],1,sd)

models$meanPLS <- rowMeans(models[c(5,6,7)])
models$sdPLS <- apply(models[, c(5,6,7)],1,sd)


#check if passes R2 and SD test

models_sub1 <- subset(models, meanRF > 0.5 & sdRF <0.1)
models_sub2 <- subset(models, meanPLS > 0.5 & sdPLS < 0.1)

models_sub <- rbind(models_sub1,models_sub2)
models_sub <- unique(models_sub)

#find max model
RF <- which(models_sub$meanRF >= models_sub$meanPLS)
PLS <- which(models_sub$meanRF < models_sub$meanPLS)

#new df
models_best <- rbind(models_sub[PLS,], models_sub[RF,] )


#return soil prop names of passing feats
PLS_vars <- models_sub[PLS,][1]
RF_vars <- models_sub[RF,][1]

PLS_vars <- as.character(PLS_vars[,1])
RF_vars <- as.character(RF_vars[,1])

#read in prediction files to merge
#currently manual enter dim sizes for starter matric - needs to change


#get dims:

fi <- paste(dirin,"RF",1,sep="/")
setwd( fi )
r <- read.csv("All predictions.csv")
r <- r[-1]
nr <- dim(r)[1]
nc <- dim(r)[2]
rf_in <- data.frame(matrix(0, ncol = nc, nrow =nr ) )

for(i in seq(1:3)){
  fi <- paste(dirin,"RF",i,sep="/")
  setwd( fi )
  r <- read.csv("All predictions.csv")
  nam <- r[1]
  ro <- r[-1]
  cols <- names(ro)
  print(dim(ro))
  
  rf_in <- rf_in + ro
}

#average
rf_in <- rf_in/3
names(rf_in) <- cols
rf_in$SSN <- nam


#same for pls
#read in prediction files to merge
#currently manual enter dim sizes for starter matric - needs to change
pl_in <- data.frame(matrix(0, ncol = nc, nrow = nr))
for(i in seq(1:3)){
  fi <- paste(dirin,"PLSM",i,sep="/")
  setwd( fi )
  q <- read.csv("All predictions.csv")
  nam <- q[1]
  qo <- q[-1]
  cols <- names(qo)
  
  pl_in <- pl_in + qo
}

#average
pl_in <- pl_in/3
names(pl_in) <- cols
pl_in$SSN <- nam


##Okay - now take the best cols from each
PLSind <- match(PLS_vars, names(pl_in))
RFind <-match(RF_vars, names(rf_in))

###
newdf <- cbind(pl_in[PLSind], rf_in[RFind])

rs <- NULL
for(i in seq(1:dim(models_best)[1]) ) {
  su <- models_best[i,]
  n <- max(su$meanRF, su$meanPLS)
  rs <- c(rs,n)
  
}

rsq <- data.frame("new" = c( rs , "mean_rsqaured") )


newdf <- cbind(newdf, rf_in$SSN)



fo <- paste(dirin, "best-predictions.csv",sep="/")
write.csv(newdf, fo)



#slim models best

models_slim <- models_best[c(1,8,9,10,11)]
models_slim[c(2:5)] <-round(models_slim[c(2:5)],2)
new_slim <- data.frame("Soil Property" = models_slim[1], "Best R2" =1, "Best SD" = 1)

for(row in 1:dim(models_slim)[1] ){
  su <- models_slim[row,]
  a <- max(su$meanRF, su$meanPLS )
  new_slim[row,2] <- a

  if(su$meanRF >= su$meanPLS){
    b <- su$sdRF
  }
  if(su$meanRF < su$meanPLS){
    b <- su$sdPLS
  }
  
  new_slim[row,3] <- b 
  
}



models_ <- models[c(1,8,9,10,11)]
models_[c(2:5)] <-round(models_[c(2:5)],2)
new_ <- data.frame("Soil Property" = models_[1], "Best R2" =1, "Best SD" = 1)

for(row in 1:dim(models_)[1] ){
  s <- models_[row,]
  a <- max(s$meanRF, s$meanPLS )
  new_[row,2] <- a
  
  if(s$meanRF >= s$meanPLS){
    b <- su$sdRF
  }
  if(s$meanRF < s$meanPLS){
    b <- su$sdPLS
  }
  
  new_[row,3] <- b 
  
}



foo <- paste(dirin, "Summary_of_best_models.csv",sep="/")

write.csv(new_slim, foo, row.names = FALSE)


SSNs <- rf_in$SSN
SSNs
print(paste("SAVED IN:",fo))

#D:\OneAcre\Google Drive\One Acre Fund\OAF Soil Lab Folder\Projects\my_feb_17\4_predicted

dirin2 <- paste(dirin,"other_summaries",sep="/")
dir.create(dirin2, showWarnings = FALSE) 
fo <- paste(dirin2, "PLS-predictions.csv",sep="/")
pl_in$SSN <- NULL
pl_in <- cbind(pl_in, SSNs)
write.csv(pl_in, fo)

fo <- paste(dirin2, "RF-predictions.csv",sep="/")
rf_in$SSN <- NULL
rf_in <- cbind(rf_in, SSNs)
write.csv(rf_in, fo)

fo <- paste(dirin2, "Full-model-summary.csv",sep="/")
write.csv(new_, fo, row.names = FALSE)


###### spit out all "best" predictions, even if fail cutoff ---------

RF_all <- which(models$meanRF >= models$meanPLS)
PLS_all <- which(models$meanRF < models$meanPLS)
PLS_all <- models[PLS_all,][1]
RF_all <- models[RF_all,][1]
RF_all <- RF_all[complete.cases(RF_all),]
PLS_all <- PLS_all[complete.cases(PLS_all),]
PLS_all <- as.character(PLS_all)
RF_all <- as.character(RF_all)
RF_all
PLS_all
match(RF_all, PLS_all)
#pl_in #rf_in are average predictions for each method

all_out <- cbind(pl_in[PLS_all],rf_in[RF_all])

all_out <- cbind(all_out, SSNs)

fooo <- paste(dirin2, "combined-predictions-including-bad-ones.csv",sep="/")
write.csv(all_out, fooo)
stop()

