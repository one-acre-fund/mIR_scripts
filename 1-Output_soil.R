
aggregate_soil <- function(){
###This will read in OPUS files and output a single CSV with spectra and file name per row
#last checked 7mar17 by MB
library(soil.spec)
library(reshape)
library(tcltk)
#memory.limit(size=12000)

gwd <- getwd()
#get the working directory for path

#pth <- choose.dir(default=getwd(), caption="Select folder with OPUS raw data")
pth <- tk_choose.dir(default=getwd(), caption="Select folder with OPUS raw data")
print(paste("Directory chosen", pth))

#grab files names from path
lst <- as.list(list.files(path=pth, pattern="*.0$", full.names=TRUE))
shortlst <- as.list(list.files(path=pth, pattern="*.0$", full.names=FALSE))
length(lst)

#lifted from R. On code ###############################
spec <- read.opus(lst,speclib = "ICRAF")@data@ab
print("1")

spectraDf<-spec

#lets try taking the 1st deriv of spectra
tmpSpec <- spectraDf[,-1]
names(tmpSpec) <- gsub("^a", "", names(tmpSpec))
row.names(spec)
print("2")
#C:\Users\Michael\Google Drive\One Acre Fund\OAF Soil Lab Folder\Projects\ke_shs_maize_paired\1_raw_opus

# Fix issue with non-float wave numbers
names(tmpSpec) <- gsub("^(2349\\.6)\\.", "\\1",
                       gsub("^(2349.6)\\.([0-9])$", "\\10\\2",
                            colnames(tmpSpec)))
#fix again for non-deriv spectra
print("3")
graphics.off(); par("mar"); par(mar=c(1,1,1,1));
tmpDer <- soil.spec::trans(tmpSpec, tr="derivative", plot.spectrogram=T)
spectraDfDer <- as.data.frame(tmpDer$trans)
rawspec <- as.data.frame(tmpDer$raw)
new.nom <- as.numeric(names(spectraDfDer))
fix.nom <- as.numeric(names(rawspec))
#class(new.nom)
#new.nom <- round(new.nom, 1)

names(spectraDfDer) <- paste0("d", new.nom)
names(rawspec) <- paste0("a", fix.nom)
names(spectraDfDer)
spectraDfDer <- cbind(SAMPLEID=spectraDf$SAMPLEID, spectraDfDer)
spec.fix <- cbind(SAMPLEID=spectraDf$SAMPLEID, rawspec)


#######################################
#date the file
#the following just makes the path happy for windows and dates it in the title
pth <- file.path(pth, .Platform$file.sep)

filnam <- paste(Sys.Date(),"Stored_spectra_d.csv",sep = "_")
filnam <- paste(pth, filnam)
#convert OS types for File out
filnam <- gsub("\\\\","/",filnam)
filnam

filnam2 <- paste(Sys.Date(),"Stored_spectra.csv", sep = "_")
filnam2 <- paste(pth, filnam2)
#convert OS types for File out
filnam2 <- gsub("\\\\","/",filnam2)

randnam <- paste(Sys.Date(),"Random_sample_selection.csv", sep = "_")




print("outputs:")
filnam
filnam2
write.csv(spectraDfDer, file=filnam, row.names = FALSE)
write.csv(spec.fix, file=filnam2, row.names = FALSE)


#output random 20% of samples
setwd(pth)
setwd("..")
dir.create("2_sampled", showWarnings = FALSE) 
setwd("./2_sampled")

set.seed(10101)
ran.df <- do.call("rbind", lapply(sample(shortlst, size=length(shortlst)/10, replace=FALSE ), as.data.frame)) 
colnames(ran.df) <- "Sample ID"
write.csv(ran.df,file=randnam, row.names = TRUE)

setwd(gwd)
}

