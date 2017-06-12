predict_properties <- function(){


library(soil.spec)

#start by setting WD to here
wd <- getwd()
outwd <- choose.dir(default="D:/OneAcre/Google Drive/One Acre Fund/OAF Soil Lab Folder/Projects",caption="please choose project directory")


#ref_path <- paste(outwd,"3_wet_chem/wet_chem.csv",sep="/")

ref_path <- choose.files(default=paste(outwd,"\\3_wet_chem\\wet_chem.csv",sep=""), caption="Please select the wet chem CSV")
raw_path <- choose.files(default=paste(outwd,"\\1_raw_opus",sep=""), caption="Please select the mIR data CSV")



#set method, either PLSM or RF
method=c("PLSM","RF")

# Begin by sourcing the scripts
source('./MB_PLS.R')
source('./RF_PLS_optimal.R')

raw<- NULL
ref<-NULL
#lets try it on OAF data

#mir
raw <- read.csv (raw_path)

#wetchem
ref <- read.csv (ref_path)

#remove NA only rows
ref <- ref[rowSums(is.na(ref)) < (dim(ref)[2]-6),]

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

# 
# ref <- ref[xn,]
# raw <- raw[yn,]
# dim(ref)
# dim(raw)

for(i in colnames(ref)[2:length(colnames(ref))]){
  ref[[i]] <- as.numeric(ref[[i]])
  
}


#test we have matches

if( sum(!is.na(match(raw$SSN, ref$SSN)))  < dim(ref)[1]*0.9   ){
  print("Found too few matches, something might be wrong, stopping run v2")
  stop()
}

print(paste("Found",sum(!is.na(match(raw$SSN, ref$SSN))),"matching entries"))

#triple run
###
#best of 3

for( xm in method){
  i
for(i in 1:3){
  print(paste("iteration",i, "of method", xm))
  m<-round(0.3*nrow(ref)) 
  test<-sample(1:nrow(ref),m)
  hout<-ref[test,]

  newwd <- paste(outwd,"4_predicted",sep="/")
  dir.create(newwd, showWarnings = FALSE) 
  newwd<- paste(newwd, xm ,sep="/")
  dir.create(newwd, showWarnings = FALSE) 
  newwd<- paste(newwd, i, sep="/")
  dir.create(newwd, showWarnings = FALSE) 

  #make new folder
  setwd(newwd) 
  print(newwd)
  try(   calibrate(newwd,raw,ref,hout, method = xm))
  #use new method - MB edited ICRAF PLS
  try( MB_PLS(newwd,raw,ref,hout,method=xm))
  } }  




setwd(wd)
print("summarising runs")
source("3-Summarise_soil.R")
summarise_runs(outwd)
}

  
  #repeat with diff holdout
  
  
  ##  #this exists for repeats
  # for( ix in 2:10){
  #   #reset to home
  #   setwd("D:/OneAcre/Google Drive/One Acre Fund/OAF Soil Lab Folder/Projects/Sample_scripts")
  #
  #   #name output folder - will auto create said folder
  #   wd <- "./test/ICRAF_RFtest/"
  #   wd <- paste(wd,ix,sep="")
  #   wd <- paste(wd,"/",sep="")
  #   print(wd)
  #   dir.create(wd, showWarnings = FALSE)
  #
  #
  #
  # # What proportion of calibration set to be held out for validation?
  # newm <- 1/ix
  # m<-round(newm*nrow(ref))
  #
  #
  #
  # test<-sample(1:nrow(ref),m)
  #
  # hout<-ref[test,]
  #
  # #  Finaly do calibrations
  # calibrate(wd,raw,ref,hout, method = "RF")
  # }