# Script for fitting calibration models for infrared spectroscopy data. ----------------------------------------------------------------
# Author: Andrew Sila , April 2016. --------------------------------

# Begin by checking the required packages are installed. --------------------------------
is.installed <- function(anypkg){
	
  is.element(anypkg, installed.packages()[,1])
  
}

required.packages <- c("hexView","caret","readr","mlbench","pROC","rpart","caretEnsemble","soil.spec","FitAR",

"doParallel","ggplot2","dplyr","downloader","prospectr","randomForest","gridExtra")

installp <- which(!is.installed(required.packages) == TRUE)

#Install missing packages

if (length (installp) > 0){
	
install.packages(required.packages[installp])
}
	
# Load packages.

suppressMessages(library(caret))
suppressMessages(library(readr))
suppressMessages(library(mlbench))
suppressMessages(library(pROC))
suppressMessages(library(rpart))
suppressMessages(library(caretEnsemble))
suppressMessages(library(soil.spec))
suppressMessages(library(FitAR))
suppressMessages(library(doParallel))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(downloader))
suppressMessages(library(prospectr))
suppressMessages(library(randomForest))
suppressMessages(library(gridExtra))

# Run in parallel to speed up the validation of the model fitting process.

registerDoParallel()
getDoParWorkers()

calibrate <- function(wd,infrared.data,reference.data,hout,method = c("RF","PLS")){
	
	if(method == "RF"){
		
  # Start by setting working directory. 
  
  setwd(wd)
  
  mir <- infrared.data
  
  ref <- reference.data
  
  hout <- hout
  
  # Exclude metadata variables.
  
  mir1 <- as.matrix(mir[,-1])
  
  wave <- as.numeric(substr(colnames(mir1),2,19))
  
  colnames(mir1) <- wave
  
  # First derivative.
  
  de1 <- trans(mir1,tr = "derivative",order = 1,gap = 23)
  
  der1 <- rev(as.data.frame(de1$trans))
  
  colnames(der1) <- paste0("m",wave)
  
  # Save derivative spectra.
  der1.ssn <- as.data.frame(cbind(as.vector(mir[,1]),der1))
  
  colnames(der1.ssn) <- c("SSN",colnames(der1))
  
  write.table(der1.ssn,file = "First derivative.csv",sep = ",",row.names = FALSE)
  
  der1.ssn<-as.data.frame(read_csv("First derivative.csv"))
  
  # Merge with first derivative preprocessed spectra.
  
  ref.mir <- merge(ref,der1.ssn,by.x = "SSN",by.y = "SSN")
  
  rc <- colnames(ref)
  
  # Which columns contains reference data?
  
  ref <- ref.mir[,rc]
  
  # Extract spectral predictors
  
  mirp <- colnames(der1.ssn)[-1]
  
  spectra <- ref.mir[,mirp]
  
  #Create two new subfolders within the current working using:
  
  b <- getwd()
  
  if(!file.exists("Models")){dir.create("Models")}
  
  if(!file.exists("calibration_plots")){dir.create("calibration_plots")}
  
  # Fit calibration models for the training set and use the testing set to validate
  
  # the models.
  
  set.seed(67523)
  
  testing <- which(ref.mir$SSN%in%hout$SSN)
  
  #with hout
  
  #Loop for calibration of all soil properties in the reference set starts here
  
  msummary <- NULL
  
  hd <- colnames(ref)[-1]
  
  for (q in 1:length(hd)){
  	
  	#Get variable names from reference table
  	
  	refq <- which(colnames(ref)%in%hd[q])
  	
  	ref.q <- ref[,refq]
  	
  	pms.a <- NULL
  	
  	pred.all <- NULL
  	
  	cal <- na.omit(cbind(as.vector(ref.q),spectra)[-testing,])
  	
  	val <- na.omit(cbind(as.vector(ref.q),spectra)[testing,])
  	
  	colnames(cal) <- c(colnames(ref)[refq],colnames(spectra))
  	
  	colnames(val) <- colnames(cal)
  	
  	trainX <- cal[,-1]
  	
  	trainY <- cal[,1]
  	
  	rf.m <- randomForest(data = cal, x = cal[,-1], y = cal[,1],
  	 
  	 ntree = 200, importance = TRUE, na.action = na.omit)
  	 
  	predicted<-predict(rf.m )
  	
  	measured<-trainY
  	
  	# compute RMSE and R-squared values for the calibration set

  	training.parameters <- round(postResample(predicted,measured),2)
  	
  	RSQ <- training.parameters[2]
  	
  	RMSE <- training.parameters[1]
  	
  	predicted.test <- predict(rf.m,val[,-1])
  	
  	measured.test <- val[,1]
  	
  	# compute RMSE and R-squared values for the validation set
  	
  	testing.parameters <- round(postResample(predicted.test,measured.test),2)
  	
  	RSP <- testing.parameters[2]
  	
  	RMSEP <- testing.parameters[1]
  	
  	model.summary <- c(hd[q],training.parameters,testing.parameters)
  	
  	msummary <- rbind(msummary,model.summary)
  	
  	# Save obtained model
  	
  	saveRDS(rf.m,file = paste0(b,"/","models/",hd[q],".rds"))
  	
  	pm <- as.data.frame(cbind(measured,predicted))
  	
  	# Create scatter plot for the predicted versus the measured
  	
	p <- ggplot(pm, aes(x = measured,y = predicted)) + 
	
	geom_point(col = "black",size = 3,alpha = 0.5) + 
	
	ggtitle(paste0("Calibration for ",hd[q])) + 
	
	xlab("Measured") + 
	
	ylab("Predicted")
	
	p <- p + stat_smooth(method = lm, se = FALSE, color = 'black',alpha = 0.1)
	
	p <- p + theme(plot.title = element_text(lineheight = 3, face = "bold", color = "black", size = 20))
	
	p <- p + theme(text = element_text(size = 20)) # this will change all text size
	 
	p <- p + annotate('text', label = paste('R^2 == ',RSQ), parse = TRUE,Inf, 
	
	-Inf,hjust = 2.5, vjust = -7.8) + annotate('text', label = paste('RMSE == ',RMSE),
	
	 parse = TRUE,Inf, -Inf,hjust = 1.8, vjust = -6.4)
	 
	# Centre title
	p <- p + theme(plot.title = element_text(hjust  = 0.5))
	
	p <- p + xlim(range(pm)) + ylim(range(pm))
	
	# Validation data
	
	pmp <- as.data.frame(cbind(measured.test,predicted.test))
	
	p2 <- ggplot(pmp, aes(x = measured.test,y = predicted.test))  + 
	
	geom_point(col = "brown",size = 3,alpha = 0.5)  + 
	
	ggtitle(paste0("Validation for ",hd[q]))  + 
	
	xlab("Measured")  + 
	
	ylab("Predicted")
	
    p2 <- p2 + stat_smooth(method = lm, se = FALSE, color = 'brown',alpha = 0.1)
    
    p2 <- p2 + theme(plot.title = element_text(lineheight = 3, face = "bold",
    
     color = "black", size = 20))
     
    # this will change all text size 

    p2 <- p2 + theme(text = element_text(size = 20)) 
    
    p2 <- p2 + annotate('text', label = paste('R^2 == ',RSP),
    
    parse = TRUE,Inf, -Inf,hjust = 2.5, vjust = -7.8) + 
    
    annotate('text', label = paste('RMSE == ',RMSEP), 
    
    parse = TRUE,Inf, -Inf,hjust = 1.8, vjust = -6.4)
	
	# Centre title
  	p2 <- p2 + theme(plot.title = element_text(hjust  = 0.5))
    
    p2 <- p2 + xlim(range(pmp)) + ylim(range(pmp))

    # Save produced scatter plots inro png files
    
	png(file = paste0(b,"/Calibration_plots/",hd[q],".png"),height = 400,width = 800)
	
	grid.arrange(p,p2,nrow = 1)
	
	dev.off()
	
	}
  
  	# Save the model summary of the calibration and validation sets.
  
  	colnames(msummary)<-c("Soil properties","OOB RMSEC","OOB Rsquared", "Holdout
  	
  	RMSEP","Holdout Rsquared")
  	
  	write.table(msummary,file = "Model_Summary.csv",sep = ",",row.names = FALSE)
  	
  	# Combine validation and calibration sets to get full model.

  	b<-getwd()
  	
  	if(!file.exists("Full_Models")){
  		dir.create("Full_Models")}
  	
  	if(!file.exists("Full_calibration_plots")){
  	
  	dir.create("Full_calibration_plots")}
  	
  	msummary <- NULL
  	
  	hd <- colnames(ref[,-1])#Exclude SSN 
  	
  	all.predicted <- NULL
  	
  	for (q in 1:length(hd)){
  	
  	refq <- which(colnames(ref)%in%hd[q])
  	
  	ref.q <- ref[,refq]
  	
  	pms.a <- NULL
  	
  	pred.all <- NULL
  	
  	cal <- cbind(as.vector(ref.q),spectra)
  	
  	colnames(cal) <- c(colnames(ref)[refq],colnames(spectra))
  	
  	p <- which(is.na(der1.ssn[,1]) ==  TRUE)
  	
  	ifelse(length(p)>0,ssn<-der1.ssn[-p,1],ssn<-der1.ssn[,1])
  	
  	ifelse(length(p)>0,der1.ssn<-der1.ssn[-p,],der1.ssn<-der1.ssn)
  	
  	#Select training and testing sets
  	
  	cal <- na.omit(cal)
  	
  	trainX <- cal[, -1]
  	
  	trainY <- cal[,1]
  	
  	# Begin model
  	
  	rf.m <- randomForest(data = cal, x = cal[,-1], y = cal[,1], ntree = 200, 
  	
  	importance = TRUE, na.action = na.omit)
  	
  	predi <- predict(rf.m)
  	
  	y <- trainY
  	
  	training.parameters <- c(hd[q],round(postResample(predi,y),3))
  	
  	msummary <- rbind(msummary,training.parameters)
  	
  	# Save the model.
  	
  	saveRDS(rf.m,file = paste0(b,"/","Full_Models/",hd[q],".rds"))
  	
  	# Make prediction using the fitted model and all the spectra in the 
  	
  	# training data
  	
  	predicted <- predict(rf.m )
  	
  	measured <- trainY
  	
  	# Computes RMSE and R-squared values for the calibration set
  	
  	training.parameters <- round(postResample(predicted,measured),3)
  	
  	RSQ <- training.parameters[2]
  	
  	RMSE <- training.parameters[1]
  	
  	predicted.test <- predict(rf.m,val[,-1])
  	
  	measured.test <- val[,1]
  	
  	testing.parameters <- round(postResample(predicted.test,measured.test),2)
  	
  	#computes RMSE and R-squared values for the validation set
  	
  	RSP <- testing.parameters[2]
  	
  	RMSEP <- testing.parameters[1]
  	
  	model.summary <- c(hd[q],training.parameters,testing.parameters)
  	
  	msummary <- rbind(msummary,model.summary)
  	
  	pm <- as.data.frame(cbind(y,predicted))
  	
  	colnames(pm) <- c("measured","predicted")
  	
  	# Create scatter plot for the predicted versus the measured
  	
  	p1 <- ggplot(pm, aes(x = measured,y = predicted)) + 
      
   geom_point(col = "brown",size = 3,alpha = 0.2) + 
   
   ggtitle(paste0("Calibration for ",hd[q])) + 
   
   xlab("Measured") + 
   
   ylab("Predicted")
   
   p1 <- p1 + stat_smooth(method = lm, se = FALSE, color = 'brown',alpha = 0.1)
   
   p1 <- p1 + theme(plot.title = element_text(lineheight = 3,
   
   face = "bold", color = "black", size = 20))
   
   # this will change all text size 

   p1 <- p1 + theme(text = element_text(size = 20))
   
   p1 <- p1 + annotate('text', label = paste('R^2 == ',RSQ),
   
   parse = TRUE,Inf, -Inf,hjust = 2.5, vjust = -7.8) + 
   
   annotate('text', label = paste('RMSE == ',RMSE), 
   
   parse = TRUE,Inf, -Inf,hjust = 1.8, vjust = -6.4)

  # Centre title
  p1 <- p1 + theme(plot.title = element_text(hjust  = 0.5))
  
  p1 <- p1 + xlim(range(pm)) + ylim(range(pm))
 
  # Save the created models 
  
  par(mfrow = c(1,1))
  
  ggsave(file = paste0(b,"/","Full_calibration_plots/",hd[q],".png"),
  
  height = 6, width = 6, units = "in",p1)
  
  # Get prediction from the full model
  
  predicted.pq <- predict(rf.m,der1.ssn[,-1])
	
  # Combine all predicted values from the full models
  
  all.predicted <- cbind(all.predicted,predicted.pq)
  
  }

  #Add sample ids (SSN) to the combined predicted values.

  all.predicted.SSN <- cbind(as.vector(ssn),all.predicted)

  colnames(all.predicted.SSN) <- c("SSN",hd)

  colnames(msummary) <- c("Soil properties","RMSEC","Rsquared")

  #Save full models' summary and predicted values
   
   write.table(msummary, file = "Full models summary.csv",
   
   sep = ",",row.names = FALSE)
   
   write.table(all.predicted.SSN, file = "All predictions.csv",
   
   sep = ",", row.names = FALSE)
   
   }

  # PLS Regression method
  # ----------------------------------------------------------------
  if(method ==  "PLS"){
  	
  	setwd(wd)
  	
  	mir <- infrared.data
  	
  	ref <- reference.data
  	
  	#Exclude metadata variables

  	mir1 <- as.matrix(mir[,-1])
  	
  	wave <- as.numeric(substr(colnames(mir1),2,19))
  	
  	colnames(mir1) <- wave
  	
  	#First derivative
  	
  	de1 <- trans(mir1,tr = "derivative",order = 1,gap = 23)
  	
  	der1 <- rev(as.data.frame(de1$trans))
  	
  	colnames(der1) <- paste0("m",wave)
  	
  	#Save derivative spectra
  	
  	der1.ssn <- as.data.frame(cbind(as.vector(mir[,1]),der1))
  	
  	colnames(der1.ssn) <- c("SSN",colnames(der1))
  	
  	write.table(der1.ssn,file = "First derivative.csv",
  	
  	sep = ",",row.names = FALSE)
  	
  	#merge with first derivative preprocessed spectra
  	
  	ref.mir <- merge(ref,der1.ssn,by.x = "SSN",by.y = "SSN")
  	
  	rc <- colnames(ref)
  	
  	#which columns contains reference data?
  	
  	ref<-ref.mir[,rc]
  	
  	#Extract spectral predictors
  	
  	mirp<-colnames(der1.ssn)[-1]
  	
  	spectra<-ref.mir[,mirp]
  	
  	#Create two new subfolders within the current working using:
  	
  	b<-getwd()
  	
  	if(!file.exists("Models")){dir.create("Models")}
  	
  	if(!file.exists("calibration_plots")){dir.create("calibration_plots")}
  	
  	# Fit calibration models for the training set and
  	
  	# use the testing set to validate the models
  	  	
  	set.seed(67523)
  	
  	testing <- which(ref.mir$SSN%in%hout$SSN) #with hout
  	
  	#Use Kennard_Stone.
  	
  	# This is an optional step just to show distribution of spectra in a PCA space.
  	
  	sel <- kenStone(spectra,k = round(0.33*nrow(spectra)),pc = .99)
  	
  	# To view selected samples, remove "#" below two lines to plot
  	
  	# plot(sel$pc[,1:2],xlab = 'PC1',ylab = 'PC2')
  	
  	# points(sel$pc[sel$model,1:2],pch = 19,col = 2)
  	
  	# points selected for calibration
  	
  	#Loop for calibration of all soil properties in the reference set starts here
  	
	msummary <- NULL
	
	hd <- colnames(ref)[-1]
	
	for (q in 1:length(hd)){
		
		refq <- which(colnames(ref)%in%hd[q])
		
		ref.q <- ref[,refq]
		
		pms.a <- NULL
		
		pred.all <- NULL
		
		cal <- cbind(as.vector(ref.q),spectra)[-testing,]
		
		val <- cbind(as.vector(ref.q),spectra)[testing,]
		
		colnames(cal) <- c(colnames(ref)[refq],colnames(spectra))
		
		colnames(val) <- colnames(cal)
		
		cal <- na.omit(cal)
		
		val <- na.omit(val)
		
		trainX <- cal[, -1]
		
		set.seed(100)
		
		colnames(cal) <- c("trainY", colnames(trainX))
		
		cal[,1] <- log(cal[,1])
	
		indx <- createFolds(trainX, returnTrain = TRUE)
		
		ctrl <- trainControl(method = "cv", index = indx)
		
		rf.m <- train(trainY~., method = "pls", data = cal,trControl =
		
		ctrl,tuneGrid = expand.grid(ncomp = 1:10),metric = "RMSE",preProc = 
		
		c("center", "scale"))
		
		# Get final model to compute coefficient for variation explained
		
		predi <- exp(predict(rf.m,rf.m$trainingData))
		
		y <- exp(cal[,1])
		
		#computes RMSE and R-squared values for the calibration set

		training.parameters <- round(postResample(predi,y),2)
		
		RSQ <- training.parameters[2]
		
		RMSE <- training.parameters[1]
		
		# Predict qth soil property of the holdoutset using
		
		# the MIR data and compare with the actual measurement
		
		predi.test <- exp(predict(rf.m,val[,-1]))
		
		y.test <- val[,1]
		
		#Get PCs used
		
		PCs <- rf.m$finalModel$ncomp
		
		#computes RMSE and R-squared values for the validation set

		testing.parameters <- round(postResample(predi.test,y.test),2)
		
		RSP <- testing.parameters[2]
		
		RMSEP <- testing.parameters[1]
		
		model.summary <- c(hd[q],PCs,training.parameters,testing.parameters)
		
		msummary <- rbind(msummary,model.summary)
		
		saveRDS(rf.m,file = paste0(b,"/","models/",hd[q],".rds"))
		
		pm <- as.data.frame(cbind(y,predi))
		
		colnames(pm) <- c("measured","predicted")
		
		# Create scatter plot for the predicted versus the measured - training set
		
		p <- ggplot(pm, aes(x = measured,y = predicted)) + 
		
		geom_point(col = "black",size = 3,alpha = 0.5) + 
		
		ggtitle(paste0("Calibration for ",hd[q])) + 
		
		xlab("Measured") + 
		
		ylab("Predicted")
		
		p <- p + stat_smooth(method = lm, se = FALSE, color = 'black',alpha = 0.1)
		
		p <- p + theme(plot.title = element_text(lineheight = 3, face = "bold",
		
		color = "black", size = 20))
		
		 # this will change all text size 
		
		p <- p + theme(text = element_text(size = 20))
		
		p <- p + annotate('text', label = paste('R^2 == ',RSQ),
		
		parse = TRUE,Inf, -Inf,hjust = 2.5, vjust = -7.8)  + 
		
		annotate('text', label = paste('RMSE == ',RMSE), 
		
		parse = TRUE,Inf, -Inf,hjust = 1.8, vjust = -6.4)
		
		# Centre title
		
      	p <- p + theme(plot.title = element_text(hjust  = 0.5))
      	
      	p <- p + xlim(range(pm)) + ylim(range(pm))
      
		#Validation data
		
		pmp <- as.data.frame(cbind(y.test,predi.test))
		
		colnames(pmp)<-c("measured.test","predicted.test")
		
		# Create scatter plot for the predicted versus the measured
		
		# the validation set
		
		p2 <- ggplot(pmp, aes(x = measured.test,y = predicted.test)) + 
		
		geom_point(col = "brown",size = 3,alpha = 0.5) + 
		
		ggtitle(paste0("Validation for ",hd[q])) + 
		
		xlab("Measured") + 
		
		ylab("Predicted")
		
		p2 <- p2 + stat_smooth(method = lm, se = FALSE, color = 'brown',
		
		alpha = 0.1)
		
		p2 <- p2 + theme(plot.title = element_text(lineheight = 3,
		
		face = "bold", color = "black", size = 20))
		
		# this will change all text size 
		
		p2 <- p2 + theme(text = element_text(size = 20))
	
		p2 <- p2 + annotate('text', label = paste('R^2 == ',RSP),
		
		parse = TRUE,Inf, -Inf,hjust = 2.5, vjust = -7.8)  +
		
		annotate('text', label = paste('RMSE == ',RMSEP),
		
		parse = TRUE,Inf, -Inf,hjust = 1.8, vjust = -6.4)
		
		# Centre title
      	p2 <- p2 + theme(plot.title = element_text(hjust  = 0.5))
      	
      	p2 <- p2 + xlim(range(pmp)) + ylim(range(pmp))

		# Save calibration and validation plots
		png(file = paste0(b,"/Calibration_plots/",hd[q],".png"),
		
		height = 400,width = 800)
		
		grid.arrange(p,p2,nrow = 1)
		
		dev.off()
		
		}
      colnames(msummary) <- c("Soil properties","PCs","LOOCV RMSEC",
      
      "LOOCV Rsquared", "Holdout RMSEP","Holdout Rsquared")
      
      write.table(msummary,file = "Model_Summary.csv",sep = ",",row.names = FALSE)
      
      # All Samples
      
      b<-getwd()
      
      if(!file.exists("Models")){dir.create("Models")}
      
      if(!file.exists("calibration_plots")){dir.create("calibration_plots")}
      
      # Begin calibration 
      
      msummary<-NULL
      
      hd<-colnames(ref[,-1])#Exclude SSN 
      
      all.predicted<-NULL
      
      for (q in 1:length(hd)) {
      	
      	refq<-which(colnames(ref)%in%hd[q])
      	
      	ref.q<-ref[,refq]
      	
      	cal<-cbind(as.vector(ref.q),spectra)
      	
      	cal<-na.omit(cal)
      	
      	trainX <-cal[, -1]
      	
      	colnames (cal) <- c("trainY",colnames(trainX))
      	
      	cal[,1] <-log(cal[,1])
      	
      	#colnames(cal)<-c(colnames(ref)[refq],colnames(spectra))
      	
      	p<-which(is.na(der1.ssn[,1]) == TRUE)
      	
      	ifelse(length(p)>0,ssn<-der1.ssn[-p,1],ssn <- der1.ssn[,1])
      	
      	ifelse(length(p)>0,der1.ssn<-der1.ssn[-p,],der1.ssn<-der1.ssn)
      	
      	#Select training and testing sets
      	
      	set.seed(100)
      	
      	indx <- createFolds(trainY, returnTrain = TRUE)
      	
      	ctrl <- trainControl(method = "cv", index = indx)
      	
      	rf.m <- train(trainY~., method = "pls", data = cal,
      	
      	trControl = ctrl,tuneGrid = expand.grid(ncomp = 1:10),
      	
      	metric = "RMSE",preProc = c("center", "scale"))
      	
      	#Save the model
      	
      	saveRDS(rf.m,file = paste0(b,"/","Full_Models/",hd[q],".rds"))
      	
      	#Get final model to compute coefficient for variation explained
      	
      	predi <- exp(predict(rf.m,rf.m$trainingData))
      	
      	y <- exp(cal[,1])
      	
      	#Get PCs used
      	
      	PCs <- rf.m$finalModel$ncomp
      	
      	training.parameters <- c(hd[q],PCs,round(postResample(predi,y),3))
      	
      	RSQ <- round(as.numeric(training.parameters[4]),2)
      	
      	RMSE <- round(as.numeric(training.parameters[3]),2)
      	
      	msummary <- rbind(msummary,training.parameters)
      	
      	#Training
      	
      	pm <- as.data.frame(cbind(y,predi))
      	
      	colnames(pm) <-c ("measured","predicted")
      	
      	png(file = paste0(b,"/","Full_calibration_plots/",hd[q],".png"),
      	
      	height =  600,width = 600)
      	
      	p1 <- ggplot(pm, aes(x = measured,y = predicted)) + 
      	
      	geom_point(col = "brown",size = 5,alpha = 0.5) + 
      	
      	ggtitle(paste0("Calibration for ",hd[q])) + 
      	
      	xlab("Measured") + 
      	
      	ylab("Predicted")
      	
      	p1 <- p1 + stat_smooth(method = lm, se = FALSE, color = 'brown',
      	
      	alpha = 0.1,size = 1.5) + 
      	
      	theme(plot.title = element_text(lineheight = 3, 
      	
      	face = "bold", color = "black", size = 20)) + 
      	
      	# this will change all text size 

      	theme(text = element_text(size = 20)) +
      	
      	annotate('text', label = paste('R^2 == ',RSQ),
      	
      	parse = TRUE,Inf, -Inf,hjust = 2.5, vjust = -7.8,size = 5)
      	
      	+ annotate('text',label = paste('RMSE == ',RMSE), 
      	
      	parse = TRUE,Inf, -Inf,hjust = 2.0, vjust = -6.8,size = 5) + 
      	
      	annotate('text', label = paste('PCs == ',PCs), 
      	
      	parse = TRUE,Inf, -Inf,hjust = 2.5, vjust = -3.9,size = 5)
      	
      	# Centre title
      	
      	p1 <- p1 + theme(plot.title = element_text(hjust  = 0.5))
      	
      	# Create scatter plot for the predicted versus the measured 
      	
      	# the combined dataset
      	
      	p1 <- p1 + xlim(range(pm)) + ylim(range(pm))

      	ggsave(file = paste0(b,"/","Full_calibration_plots/",hd[q],".png"),
      	
      	height = 6, width = 6, units = "in", p1)
      	
      	prediction.f <- round(exp(predict(rf.m,der1.ssn[,-1])),2)
      	
      	all.predicted <- cbind(all.predicted,prediction.f)
      	
      	}
      	
      	#Combine the predicted values together
      	
      	all.predicted.SSN <- cbind(as.vector(ssn),all.predicted)
      	
      	colnames(all.predicted.SSN) <- c("SSN",hd)
      	
      	colnames(msummary)<-c("Soil properties","PCs","RMSEC","Rsquared")
      	
      	#Save full model summaries
      	
      	write.table(msummary, file = "Full models summary.csv",sep = ",",
      	
      	row.names = FALSE)
      	
      	#Save the linked file
      	
      	write.table(all.predicted.SSN, file = "All predictions.csv",
      	
      	sep = ",", row.names = FALSE)
     }
   }