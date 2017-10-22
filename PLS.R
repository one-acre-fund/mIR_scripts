# Script for fitting calibration models for infrared spectroscopy data. ----------------------------------------------------------------
# Author: Andrew Sila , April 2016. --------------------------------


###
# Begin by checking the required packages are installed. --------------------------------
is.installed <- function(anypkg){
  is.element(anypkg, installed.packages()[,1])
}

required.packages <- c("caret","readr","soil.spec", "doMC",
                       "doParallel","ggplot2","dplyr","prospectr","gridExtra")

installp <- which(!is.installed(required.packages) == TRUE)

#Install missing packages
if (length (installp) > 0){
  install.packages(required.packages[installp])
}

# Load packages.
suppressMessages(library(caret))
suppressMessages(library(readr))
library(soil.spec)
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(prospectr))
suppressMessages(library(gridExtra))

# Run in parallel to speed up the validation of the model fitting process.
suppressMessages(library(doMC)) #Max Sop setup - comment out if you don't a multiple core computer
registerDoMC(detectCores())         #Max Sop setup - this should automatically detetect the #cores in your machine and take advantage of it for faster computation

## RUN FROM HERE
PLS <- function(wd, infrared.data, reference.data, testing, 
                raw_test, ref_test, method = c("PLSM","SVM","XGBM")){
  
  if(method=="PLSM"){
    
    # PLS Regression method
    # ----------------------------------------------------------------
    setwd(wd)
    mir <- infrared.data
    ref <- reference.data
    
    #Exclude metadata variables
    mir1 <- as.matrix(mir[, -1])
    raw_test1 <- as.matrix(raw_test[, -1])
    mir1 <- savitzkyGolay(mir1, p = 3, w = 21, m = 0)
    raw_test1 <- savitzkyGolay(raw_test1, p = 3, w = 21, m = 0)
    wave <- as.numeric(substr(colnames(mir1),2,19))
    wave_test <- as.numeric(substr(colnames(raw_test1),2,19))
    colnames(mir1) <- wave
    colnames(raw_test1) <- wave_test
    
    #First derivative
    de1 <- trans(mir1, tr = "derivative", order = 1, gap = 23)
    der1 <- rev(as.data.frame(de1$trans))
    colnames(der1) <- paste0("m", wave)
    
    de2 <- trans(raw_test1, tr = "derivative", order = 1, gap = 23)
    der2 <- rev(as.data.frame(de2$trans))
    colnames(der2) <- paste0("m", wave_test)
    
    #Save derivative spectra
    der1.ssn <- as.data.frame(cbind(as.vector(mir[, 1]), der1))
    colnames(der1.ssn) <- c("SSN",colnames(der1))
    der2.ssn <- as.data.frame(cbind(as.vector(raw_test[, 1]), der2))
    colnames(der2.ssn) <- c("SSN",colnames(der2))
    
    colnames(ref)[1] <- "SSN"
    
    #colnames(der1.ssn)
    ref.mir <- merge(ref, der1.ssn, by.x = "SSN", by.y = "SSN")
    rc <- colnames(ref)
    
    #which columns contains reference data?
    ref<-ref.mir[, rc]
    
    #Extract spectral predictors
    mirp <- colnames(der1.ssn)[-1]
    spectra <- ref.mir[, mirp]
    
    b <- getwd()
    
    if(!file.exists("Models")){dir.create("Models")}
    
    if(!file.exists("calibration_plots")){dir.create("calibration_plots")}
    
    # Fit calibration models for the training set and
    # use the testing set to validate the models
    msummary <- NULL
    hd <- colnames(ref)[-1]
    print("training calibration models")
    pb <- txtProgressBar(min = 0, max = 2*length(hd), style = 3)	
    best_params <- c()
    for (q in 1:length(hd)){
      setTxtProgressBar(pb, q+length(hd))
      refq <- which(colnames(ref)%in%hd[q])
      ref.q <- ref[, refq]
      pms.a <- NULL
      pred.all <- NULL
      cal <- cbind(as.vector(ref.q),spectra)[-testing,]
      val <- cbind(as.vector(ref.q),spectra)[testing,]
      #check val and cal
      colnames(cal) <- c(colnames(ref)[refq], colnames(spectra))
      colnames(val) <- colnames(cal)
      cal <- na.omit(cal)
      val <- na.omit(val)
      #check val and cal
      trainX <- cal[, -1] #drop target var
      colnames(cal) <- c("trainY", colnames(trainX))
      cal[, 1] <- log(cal[, 1]) 
      
      #set.seed(101)
      indx <- createFolds(cal$trainY, returnTrain = TRUE, k = 10) 
      ctrl <- trainControl(method = "cv", index = indx)
      #set.seed(99)
      pls.fit <- train(trainY ~., method = "kernelpls", data = cal,
                       trControl =ctrl,tuneGrid = expand.grid(ncomp = 2:5),
                       metric = "RMSE", preProc = c("center", "scale"))
      # Get final model to compute coefficient for variation explained
      predi <- exp(predict(pls.fit, pls.fit$trainingData))
      y <- exp(cal[, 1])
      
      #computes RMSE and R-squared values for the calibration set
      training.parameters <- round(postResample(predi, y), 2)
      RSQ <- training.parameters[2]
      RMSE <- training.parameters[1]
      
      # Predict qth soil property of the holdoutset using
      # the MIR data and compare with the actual measurement
      predi.test <- exp(predict(pls.fit, val[,-1]))
      y.test <- val[, 1]
      
      #Get PCs used
      PCs <- pls.fit$finalModel$ncomp
      best_params[q] <- PCs
      #computes RMSE and R-squared values for the validation set
      testing.parameters <- round(postResample(predi.test, y.test), 2)
      RSP <- testing.parameters[2]
      RMSEP <- testing.parameters[1]
      model.summary <- c(hd[q], PCs, training.parameters, testing.parameters)
      msummary <- rbind(msummary, model.summary)
      
      #saveRDS(pls.fit, file = paste0(b,"/","Models/", hd[q], ".rds"))
      
      pm <- as.data.frame(cbind(y, predi))
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
      p <- p + annotate('text', label = paste('R^2 == ', RSQ),
                        parse = TRUE,Inf, -Inf,hjust = 2.5, vjust = -7.8)  + 
        annotate('text', label = paste('RMSE == ',RMSE), 
                 parse = TRUE,Inf, -Inf,hjust = 1.8, vjust = -6.4)
      
      # Centre title
      p <- p + theme(plot.title = element_text(hjust  = 0.5))
      p <- p + xlim(range(pm)) + ylim(range(pm))
      
      #Validation data
      pmp <- as.data.frame(cbind(y.test, predi.test))
      colnames(pmp) <- c("measured.test","predicted.test")
      
      # Create scatter plot for the predicted versus the measured
      # the validation set
      p2 <- ggplot(pmp, aes(x = measured.test,y = predicted.test)) + 
        geom_point(col = "brown",size = 3,alpha = 0.5) + 
        ggtitle(paste0("Validation for ",hd[q])) + 
        xlab("Measured") + 
        ylab("Predicted")
      p2 <- p2 + stat_smooth(method = lm, se = FALSE, color = 'brown', alpha = 0.1)
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
      png(file = paste0(b,"/calibration_plots/", hd[q],".png"),
          height = 400, width = 800)
      grid.arrange(p,p2,nrow = 1)
      
      dev.off()
      gc()
    }
    colnames(msummary) <- c("Soil properties","PCs","CV RMSEC","CV Rsquared",
                            "Holdout RMSEP","Holdout Rsquared")
    write.table(msummary, file = "Model_Summary.csv", 
                sep = ",",row.names = FALSE)
    close(pb)
    # All Samples
    print("training full models")
    pb <- txtProgressBar(min = 0, max = 2*length(hd), style = 3)
    
    b<-getwd()
    
    # if(!file.exists("Full_Models")){
    #   dir.create("Full_Models")}
    # 
    # Begin calibration 
    msummary <- NULL
    msummary2 <- NULL
    hd <- colnames(ref[-1])#Exclude SSN 
    all.predicted <- NULL
    all.predicted.test <- NULL
    for(q in 1:length(hd)) {
      
      setTxtProgressBar(pb, q+length(hd))
      
      refq <- which(colnames(ref) %in% hd[q])
      ref.q <- ref[, refq]
      cal <- cbind(as.vector(ref.q), spectra)
      cal <- na.omit(cal)
      trainX <- cal[, -1]
      #looks like at this point trainX and cal are exactly same but cal includes the soil property to be predicted
      
      #log the property...?
      cal[,1] <- log(cal[,1]) 
      cal <- na.omit(cal)
      #this err, renames trainY back to the original name.... e.g. phosphorus
      colnames(cal) <- c(colnames(ref)[refq], colnames(spectra))
      
      #find empty SSN names
      p <- which(is.na(der1.ssn[, 1]) == TRUE)
      ifelse(length(p)>0,ssn<-der1.ssn[-p,1],ssn <- der1.ssn[,1])
      ifelse(length(p)>0,der1.ssn<-der1.ssn[-p,],der1.ssn<-der1.ssn)
      
      p <- which(is.na(der2.ssn[, 1]) == TRUE)
      ifelse(length(p) > 0, ssn2 <- der2.ssn[-p, 1], ssn2 <- der2.ssn[,1])
      ifelse(length(p) > 0, der2.ssn <- der2.ssn[-p,], der2.ssn <- der2.ssn)
      
      #Select training and testing sets
      colnames (cal) <- c("trainY", colnames(trainX))
      
      #make folds of data trainX - which does NOT have soil prop (phosphorus, pH etc)
      #set.seed(101)
      indx <- createFolds(cal$trainY, returnTrain = TRUE, k = 10) 
      ctrl <- trainControl(method = "cv", index = indx)
      ncomp <- best_params[q]
      #set.seed(99)
      pls.fit <- train(trainY ~ ., method = "kernelpls", data = cal,
                       trControl = ctrl, tuneGrid = expand.grid(ncomp = ncomp),
                       metric = "RMSE", preProc = c("center", "scale"))
      
      #saveRDS(pls.fit, file = paste0(b,"/","Full_Models/", hd[q], ".rds"))
      
      predi <- exp(predict(pls.fit, pls.fit$trainingData))
      y <- exp(cal[, 1])
      
      #Get PCs used
      PCs <- pls.fit$finalModel$ncomp
      training.parameters <- c(hd[q], PCs, round(postResample(predi,y), 3))
      RSQ <- round(as.numeric(training.parameters[4]),2)
      RMSE <- round(as.numeric(training.parameters[3]),2)
      msummary <- rbind(msummary, training.parameters)
      
      pm <- as.data.frame(cbind(y, predi))
      colnames(pm) <-c ("measured","predicted")
      fo <- paste(b,"/","Full_calibration_plots/", sep = "")
      dir.create(fo, showWarnings = FALSE)
      fo <- paste(fo, hd[q],".png",sep = "")
      
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
        annotate('text', label = paste('R^2 == ', RSQ),
                 parse = TRUE,Inf, -Inf,hjust = 2.5, vjust = -7.8,size = 5) + 
        annotate('text', label = paste('PCs == ', PCs), 
                 parse = TRUE,Inf, -Inf,hjust = 2.5, vjust = -3.9, size = 5)
      
      # Centre title
      p1 <- p1 + theme(plot.title = element_text(hjust  = 0.5))
      
      # Create scatter plot for the predicted versus the measured 
      # the combined dataset
      p1 <- p1 + xlim(range(pm)) + ylim(range(pm))
      ggsave(file = paste0(b,"/","Full_calibration_plots/", hd[q],".png"),
             height = 6, width = 6, units = "in", p1)
      
      prediction.f <- round(exp(predict(pls.fit, der1.ssn[, -1])), 2)
      all.predicted <- cbind(all.predicted, prediction.f)
      
      prediction.test <- round(exp(predict(pls.fit, der2.ssn[, -1])), 2)
      all.predicted.test <- cbind(all.predicted.test, prediction.test)
      rsq.test <- c(hd[q] , round(postResample(prediction.test, ref_test[, hd[q]])[2], 3))
      msummary2 <- rbind(msummary2, rsq.test)
    } ##end loop
    
    #Combine the predicted values together
    all.predicted.SSN <- cbind(as.vector(ssn), all.predicted)
    colnames(all.predicted.SSN) <- c("SSN", hd)
    colnames(msummary) <- c("Soil properties","PCs","RMSEC","Rsquared")
    
    all.predicted.test.SSN <- cbind(as.vector(ssn2), all.predicted.test)
    colnames(all.predicted.test.SSN) <- c("SSN", hd)
    colnames(msummary2) <- c("Soil properties","Holdout Rsquared")
    
    #Save full model summaries
    write.table(msummary, file = "Full_models_summary.csv", sep = ",", 
                row.names = FALSE)
    #Save the linked file
    write.table(all.predicted.SSN, file = "All_predictions.csv",
                sep = ",", row.names = FALSE)
    
    external_mir_dir <- "../../../External_MIR"
    sub_dir <- substring(b, 106, 111)
    #if(!file.exists(paste0(external_mir_dir, "/", sub_dir))){
      dir.create(file.path(external_mir_dir, sub_dir), recursive = TRUE)
     # } 
    #print(b)
    write.table(msummary2, file = paste0(external_mir_dir,"/", sub_dir, "/CV_models_summary.csv"), 
                row.names = FALSE, sep = ",")
    write.table(all.predicted.test.SSN, file = paste0(external_mir_dir,"/", sub_dir, "/All_predictions.csv"),
                row.names = FALSE, sep = ",") 
    close(pb)
    gc()
  }
}



