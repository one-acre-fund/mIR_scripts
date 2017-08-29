# Script for fitting calibration models for infrared spectroscopy data. ----------------------------------------------------------------
# Author: Andrew Sila , April 2016. --------------------------------


###
# Begin by checking the required packages are installed. --------------------------------
is.installed <- function(anypkg){
  is.element(anypkg, installed.packages()[,1])
}

required.packages <- c("caret","readr","soil.spec","doMC","xgboost",
                       "ggplot2","dplyr","prospectr","gridExtra", "wavelets")

installp <- which(!is.installed(required.packages) == TRUE)

#Install missing packages

if (length (installp) > 0){
  install.packages(required.packages[installp])
}

# Load packages.

suppressMessages(library(caret))
suppressMessages(library(readr))
suppressMessages(library(wavelets))
library(soil.spec)
suppressMessages(library(xgboost))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(prospectr))
suppressMessages(library(gridExtra))

# Run in parallel to speed up the validation of the model fitting process.
suppressMessages(library(doMC)) #Max Sop setup - comment out if you don't a multiple core computer
registerDoMC(detectCores())         #Max Sop setup - this should automatically detetect the #cores in your machine and take advantage of it for faster computation

## RUN FROM HERE
XGBM <- function(wd, infrared.data, reference.data, testing, method = c("PLSM","SVM","XGBM")){
  
  if(method=="XGBM"){
    
    # PLS Regression method
    # ----------------------------------------------------------------
    setwd(wd)
    mir <- infrared.data
    ref <- reference.data
    
    #Exclude metadata variables
    mir1 <- as.matrix(mir[, -1])
    mir1 <- savitzkyGolay(mir1, p = 3, w = 21, m = 0)
    wave <- as.numeric(substr(colnames(mir1),2,19))
    colnames(mir1) <- wave
    
    #mir pre-processing
    de1 <- trans(mir1, tr = "derivative", order = 1, gap = 23)
    der1 <- rev(as.data.frame(de1$trans))
    #der1 <- Derivative(mir1) #this is a bit faster
    colnames(der1) <- paste0("m", wave)
    der1 <- HaarTransform(der1, 3)
    #der1 <- continuumRemoval(der1, type='R', method='substraction')
    nzv_cols <- nearZeroVar(der1)
    if(length(nzv_cols) > 0) der1 <- der1[, -nzv_cols]
    der1 <- der1[, colMeans(is.na(der1)) <= 0.5]
    colnames(der1) <- paste0("m", colnames(der1))
    
    #Save transformed spectra
    der1.ssn <- as.data.frame(cbind(as.vector(mir[, 1]), der1))
    colnames(der1.ssn) <- c("SSN", colnames(der1))
    #write_csv(der1.ssn, "transformed_mir_data.csv")
    
    #merge with first derivative preprocessed spectra
    colnames(ref)[1] <- "SSN"
    ref.mir <- merge(ref, der1.ssn, by.x = "SSN", by.y = "SSN")
    write_csv(ref.mir, "ref.mir.csv")
    rc <- colnames(ref)
    
    #which columns contains reference data?
    ref <- ref.mir[, rc]
    
    #Extract spectral predictors
    mirp <- colnames(der1.ssn)[-1]
    spectra <- ref.mir[, mirp]
    
    #Create two new subfolders within the current working using:
    b<-getwd()
    
    if(!file.exists("Models")){dir.create("Models")}
    if(!file.exists("calibration_plots")){dir.create("calibration_plots")}
    
    # Fit calibration models for the training set and
    # use the testing set to validate the models
    
    # newm <- 0.35
    # #set.seed(sample(1:1000,1))
    # testing <- round(runif(round(length(ref.mir$SSN)*newm), min=0, max=length(ref.mir$SSN)))
    # testing <- unique(testing)
    # print(paste("data points in testing are ", length(testing)))
    
    msummary <- NULL
    hd <- colnames(ref)[-1]
    print("training calibration models")
    pb <- txtProgressBar(min = 0, max = 2*length(hd), style = 3)	
    best_iterations <- c()
    for(q in 1:length(hd)){
      setTxtProgressBar(pb, q+length(hd))
      refq <- which(colnames(ref)%in%hd[q])
      ref.q <- ref[, refq]
      pms.a <- NULL
      pred.all <- NULL
      cal <- cbind(as.vector(ref.q),spectra)[-testing,]
      val <- cbind(as.vector(ref.q),spectra)[testing,]
      colnames(cal) <- c(colnames(ref)[refq], colnames(spectra))
      colnames(val) <- colnames(cal)
      cal <- na.omit(cal)
      val <- na.omit(val)
      trainX <- as.matrix(cal[, -1]) #drop target var and convert to matrix
      trainY <- cal[, 1]
      testX <- as.matrix(val[, -1])
      testY <- val[, 1]
      
      dtrain <- xgb.DMatrix(data = trainX, label = trainY)
      dtest <- xgb.DMatrix(data = testX, label = testY)
      
      watchlist <- list(train=dtrain, eval=dtest)
      
      params <- list(objective = "reg:linear"
                     , eta = 0.07
                     , subsample = 0.7
                     , colsample_bytree = 0.9
                     , min_child_weight = 10
                     , max_depth = 9)
      nrounds <- 500
      set.seed(71109)
      xgb.fit <- xgb.train(params = params, data=dtrain, watchlist=watchlist,
                     nround = nrounds, eval_metric = evalerror,
                     print_every_n = 10,  early_stopping_rounds = 10, 
                     maximize = TRUE, verbose = 1)
      
      best_iter <- xgb.fit$best_iteration
      train_log <- xgb.fit$evaluation_log
      best_iterations[q] <- best_iter
      #R-squared values for the calibration set
      RSQ <- round(train_log[best_iter, 2], 2)
      RSQ.eval <- round(train_log[best_iter, 3], 2)
      model.summary <- c(hd[q], RSQ, RSQ.eval)
      msummary <- rbind(msummary, model.summary)
      
      saveRDS(xgb.fit, file = paste0(b,"/","Models/", hd[q], ".rds"))
      predi <- round(predict(xgb.fit, trainX), 2)
      pm <- as.data.frame(cbind(trainY, predi))
      colnames(pm) <- c("measured","predicted")
      
      # Create scatter plot for the predicted versus the measured - training set
      p <- ggplot(pm, aes(x = measured,y = predicted)) + 
        geom_point(col = "black",size = 3,alpha = 0.5) + 
        ggtitle(paste0("Calibration for ", hd[q])) + 
        xlab("Measured") + 
        ylab("Predicted")
      p <- p + stat_smooth(method = lm, se = FALSE, color = 'black',alpha = 0.1)
      p <- p + theme(plot.title = element_text(lineheight = 3, face = "bold",
                                               color = "black", size = 20))
      
      # this will change all text size 
      p <- p + theme(text = element_text(size = 20))
      p <- p + annotate('text', label = paste('R^2 == ', RSQ), 
        parse = TRUE,Inf, -Inf,hjust = 2.5, vjust = -7.8) 
      
      # Centre title
      p <- p + theme(plot.title = element_text(hjust  = 0.5))
      p <- p + xlim(range(pm)) + ylim(range(pm))
      
      #Validation data
      predi.test <- round(predict(xgb.fit, testX), 2)
      pmp <- as.data.frame(cbind(testY, predi.test))
      colnames(pmp) <- c("measured.test","predicted.test")
      
      # Create scatter plot for the predicted versus the measured
      # the validation set
      p2 <- ggplot(pmp, aes(x = measured.test,y = predicted.test)) + 
        geom_point(col = "brown",size = 3,alpha = 0.5) + 
        ggtitle(paste0("Validation for ",hd[q])) + 
        xlab("Measured") + 
        ylab("Predicted")
      
      p2 <- p2 + stat_smooth(method = lm, se = FALSE, color = 'brown',alpha = 0.1)
      p2 <- p2 + theme(plot.title = element_text(lineheight = 3, face = "bold", color = "black", size = 20))
      # this will change all text size 
      p2 <- p2 + theme(text = element_text(size = 20))
      p2 <- p2 + annotate('text', label = paste('R^2 == ', RSQ.eval), 
        parse = TRUE,Inf, -Inf,hjust = 2.5, vjust = -7.8) 
      
      # Centre title
      p2 <- p2 + theme(plot.title = element_text(hjust  = 0.5))
      p2 <- p2 + xlim(range(pmp)) + ylim(range(pmp))
      
      # Save calibration and validation plots
      png(file = paste0(b,"/calibration_plots/",hd[q],".png"), height = 400, width = 800)
      grid.arrange(p,p2,nrow = 1)
      
      dev.off()
      gc()
    }
    colnames(msummary) <- c("Soil properties","CV Rsquared", "Holdout Rsquared")
    write.table(msummary, file = "Model_Summary.csv", sep = ",",row.names = FALSE)
    close(pb)
    # All Samples
    print("training full models")
    pb <- txtProgressBar(min = 0, max = 2*length(hd), style = 3)
    
    b<-getwd()
    
    if(!file.exists("Full_Models")){
      dir.create("Full_Models")}
    
    msummary <- NULL
    hd <- colnames(ref[-1])#Exclude SSN 
    all.predicted <- NULL
    print("training full models")
    for(q in 1:length(hd)) {
      
      setTxtProgressBar(pb, q+length(hd))
      refq <- which(colnames(ref) %in% hd[q])
      ref.q <- ref[, refq]
      cal <- cbind(as.vector(ref.q), spectra)
      cal <- na.omit(cal)
      trainX <- cal[, -1]
      trainY <- cal[, 1]
      #looks like at this point trainX and cal are exactly same but cal includes the soil property to be predicted
      #rename property as trainY
      cal <- na.omit(cal)
      colnames(cal) <- c(colnames(ref)[refq], colnames(spectra))
      
      #find empty SSN names
      p <- which(is.na(der1.ssn[, 1]) == TRUE)
      ifelse(length(p)>0,ssn<-der1.ssn[-p,1],ssn <- der1.ssn[,1])
      ifelse(length(p)>0,der1.ssn<-der1.ssn[-p,],der1.ssn<-der1.ssn)
      
      trainX <- as.matrix(trainX) #drop target var and convert to matrix
      dtrain <- xgb.DMatrix(data = trainX, label = trainY)
      watchlist <- list(train=dtrain)
      
      nrounds <- best_iterations[q]
      set.seed(71109)
      xgb.fit <- xgb.train(params = params, data=dtrain, watchlist=watchlist,
                     nround = nrounds, eval_metric = evalerror,
                     print_every_n = 10,  early_stopping_rounds = 10, 
                     maximize = TRUE, verbose = 1)
      
      saveRDS(xgb.fit, file = paste0(b,"/","Full_Models/", hd[q], ".rds"))
      
      best_iter <- xgb.fit$best_iteration
      train_log <- xgb.fit$evaluation_log
      
      #R-squared values for the calibration set
      RSQ <-  round(train_log[best_iter, 2], 2)
      model.summary <- c(hd[q], RSQ)
      msummary <- rbind(msummary, model.summary)
      
      predi <- round(predict(xgb.fit, trainX), 2)
      pm <- as.data.frame(cbind(trainY, predi))
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
        theme(plot.title = element_text(lineheight = 3, face = "bold", color = "black", size = 20)) + 
        
        # this will change all text size 
        theme(text = element_text(size = 20)) + annotate('text', label = paste('R^2 == ', RSQ),
                 parse = TRUE,Inf, -Inf,hjust = 2.5, vjust = -7.8,size = 5)
      
      # Centre title
      p1 <- p1 + theme(plot.title = element_text(hjust  = 0.5))
      
      # Create scatter plot for the predicted versus the measured 
      # the combined dataset
      p1 <- p1 + xlim(range(pm)) + ylim(range(pm))
      ggsave(file = paste0(b,"/","Full_calibration_plots/", hd[q],".png"),
             height = 6, width = 6, units = "in", p1)
      prediction.f <- round(predict(xgb.fit, as.matrix(der1.ssn[, -1])), 2)
      all.predicted <- cbind(all.predicted, prediction.f)
    } ##end loop
    
    #Combine the predicted values together
    all.predicted.SSN <- cbind(as.vector(ssn), all.predicted)
    colnames(all.predicted.SSN) <- c("SSN", hd)
    colnames(msummary) <- c("Soil properties","Rsquared")
    
    #Save full model summaries
    write.table(msummary, file = "Full_models_summary.csv", sep = ",", 
                row.names = FALSE)
    
    #Save the linked file
    write.table(all.predicted.SSN, file = "All_predictions.csv",
                sep = ",", row.names = FALSE)
    close(pb)
    gc()
  }
}



