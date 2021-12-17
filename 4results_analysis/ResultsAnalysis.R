############################################################
################ Simulation: Equal Cases ###################
############################################################
setwd("insert the path to the related results data. ")
library(ggplot2)
library(dplyr)
library(forcats)
library(viridis)
## Percentage of MSE reduction 
modelnames <- c("OW", "TOP3", "RP", "CWM", "SSIN", "SSDE", "REGCOV", "LAS", "SA")
samplesize <- c(16, 32, 64, 128, 256)
#dat <- read.csv("InfStates.csv", header = F, stringsAsFactors = F)
dat <- read.csv("EqualStates.csv", header = T, stringsAsFactors = F)
dat <- dat[ , -3]
colnames(dat) <- c("case", "samplesize", "p1", "p2", "p3", "p4", "p5", "p6", "p7", "p8",
                   "cv1", "cv2", "cv3", "cv4", "cv5", "cv6", "cv7", "cv8",
                   "MSEout1", "MSEout2", "MSEout3", "MSEout4", "MSEout5", "MSEout6", "MSEout7", "MSEout8", "MSEout9",
                   "MSEttest1", "MSEttest2", "MSEttest3", "MSEttest4", "MSEttest5", "MSEttest6", "MSEttest7", "MSEttest8")
head(dat)
#### MSE Reduction: aggregating for each weighting method ####
pvalue <- dat[ , 3:10]
cverror <- dat[ , 11:18]
MSEmean <- dat[ , 19:27]
PMSEtoSA <- (MSEmean[ , 9] - MSEmean)/MSEmean[ , 9]
PMSEdata.test <- matrix(NA, nrow = nrow(PMSEtoSA), ncol = ncol(PMSEtoSA) - 1)
PMSEdata.cv   <- matrix(NA, nrow = nrow(PMSEtoSA), ncol = ncol(PMSEtoSA) - 1)
crit <- 0.05
for(k in 1:(length(modelnames) - 1)){
  test_select_wa_index <- which(pvalue[ , k] <= crit)
  test_select_sa_index <- which(pvalue[ , k] > crit)
  cv_select_wa_index   <- which(cverror[ , k] <= 0)
  cv_select_sa_index   <- which(cverror[ , k] > 0)
  PMSEdata.test[test_select_wa_index, k] <- PMSEtoSA[test_select_wa_index, k]
  PMSEdata.test[test_select_sa_index, k] <- PMSEtoSA[test_select_sa_index, 9]
  PMSEdata.cv[cv_select_wa_index, k]     <- PMSEtoSA[cv_select_wa_index, k]
  PMSEdata.cv[cv_select_sa_index, k]     <- PMSEtoSA[cv_select_sa_index, 9]
}
# remove the outliers for RegCov weighting method 
PMSEtoSA <- PMSEtoSA[ , -9]
mean.indw <- apply(PMSEtoSA[ , -7], 2, function(x) mean(na.omit(x)))
mean.test <- apply(PMSEdata.test[ , -7], 2, function(x) mean(na.omit(x)))
mean.cv <- apply(PMSEdata.cv[ , -7], 2, function(x) mean(na.omit(x)))
se.indw <- apply(PMSEtoSA[ , -7], 2, function(x) sd(na.omit(x))/sqrt(length(na.omit(x))))
se.test <- apply(PMSEdata.test[ , -7], 2, function(x) sd(na.omit(x))/sqrt(length(na.omit(x))))
se.cv <- apply(PMSEdata.cv[ , -7], 2, function(x) sd(na.omit(x))/sqrt(length(na.omit(x))))
# Visualization: including individual weighting methods 
plotdata1 <- data.frame(rbind(mean.indw, mean.test, mean.cv))
plotdata2 <- data.frame(rbind(se.indw, se.test, se.cv))
colnames(plotdata1) <- modelnames[c(1:6, 8)]
print(plotdata1)
print(plotdata2)
print(min(plotdata1, na.rm = T))
print(max(plotdata1, na.rm = T))

#range <- c(-5, 5) # for inflaton rate data 
range <- c(-14, 0) # for unemployment rate data
stepsize <- 2
# plot
setEPS()
postscript("PMSE_model_equal.eps", width = 7, height = 3.5)
bars <- barplot(as.matrix(plotdata1*100), ylim = c(range[1], range[2]), yaxt = "n", 
                xlab = "Weighting Model", ylab = "%MSE Reduction to SA", density = c(90, 45, 0), beside = TRUE)
axis(2, seq(range[1], range[2], by = stepsize), labels = paste(seq(range[1], range[2], by = stepsize), "%", sep = ""), las = 1)
segments(bars, as.vector(as.matrix(plotdata1))*100 - as.vector(as.matrix(plotdata2))*300, 
         bars, as.vector(as.matrix(plotdata1))*100 + as.vector(as.matrix(plotdata2))*300, lwd = 2)
arrows(bars, as.vector(as.matrix(plotdata1))*100 - as.vector(as.matrix(plotdata2))*300, 
       bars, as.vector(as.matrix(plotdata1))*100 + as.vector(as.matrix(plotdata2))*300, lwd = 2, angle = 90, code = 3, length = 0.05)
legend("bottomright", c("WA", "Test", "CV"), density = c(90, 45, 0), bg = "white", bty = "n")
dev.off()

#### MSE Reduction: aggregating for each sample size ####
samplesize <- c(16, 32, 64, 128, 256)
sim_runs <- min(data.frame(table(dat$samplesize))$Freq)
PMSEdata.indw <- c()
PMSEdata.test <- c()
PMSEdata.cv   <- c()
crit <- 0.05
for(i in 1:length(samplesize)){
  subdat <- dat[which(dat$samplesize == samplesize[i]), ]
  ## load the data 
  MSEmean.i <- subdat[ , 19:27]
  PMSEtoSA.i <- (MSEmean.i[ , 9] - MSEmean.i)/MSEmean.i[ , 9]
  cverror.i <- subdat[ , 11:18]
  pvalue.i <- subdat[ , 3:10] 
  ## find the MSE for test and cv
  model_k <- c(1:8)
  PMSEdata.test.i <- matrix(NA, nrow = nrow(PMSEtoSA.i), ncol = length(model_k))
  PMSEdata.cv.i   <- matrix(NA, nrow = nrow(PMSEtoSA.i), ncol = length(model_k))
  for(k in model_k){
    test_select_wa_index <- which(pvalue.i[ , k] <= crit)
    test_select_sa_index <- which(pvalue.i[ , k] > crit)
    cv_select_wa_index   <- which(cverror.i[ , k] <= 0)
    cv_select_sa_index   <- which(cverror.i[ , k] > 0)
    PMSEdata.test.i[test_select_wa_index, k] <- PMSEtoSA.i[test_select_wa_index, k]
    PMSEdata.test.i[test_select_sa_index, k] <- PMSEtoSA.i[test_select_sa_index, 9]
    PMSEdata.cv.i[cv_select_wa_index, k]     <- PMSEtoSA.i[cv_select_wa_index, k]
    PMSEdata.cv.i[cv_select_sa_index, k]     <- PMSEtoSA.i[cv_select_sa_index, 9]
  }
  PMSEtoSA.i <- PMSEtoSA.i[ , -c(7, 9)]
  PMSEdata.test.i <- PMSEdata.test.i[ , -7]
  PMSEdata.cv.i <- PMSEdata.cv.i[ , -7]
  PMSEdata.indw <- cbind(PMSEdata.indw, as.vector(as.matrix(PMSEtoSA.i[1:sim_runs, ])))
  PMSEdata.test <- cbind(PMSEdata.test, as.vector(PMSEdata.test.i[1:sim_runs, ]))
  PMSEdata.cv <- cbind(PMSEdata.cv, as.vector(PMSEdata.cv.i[1:sim_runs, ]))
}
mean.indw <- apply(PMSEdata.indw, 2, function(x) mean(na.omit(x)))
se.indw <- apply(PMSEdata.indw, 2, function(x) sd(na.omit(x))/sqrt(length(na.omit(x))))
mean.test <- apply(PMSEdata.test, 2, function(x) mean(na.omit(x)))
se.test <- apply(PMSEdata.test, 2, function(x) sd(na.omit(x))/sqrt(length(na.omit(x))))
mean.cv <- apply(PMSEdata.cv, 2, function(x) mean(na.omit(x)))
se.cv <- apply(PMSEdata.cv, 2, function(x) sd(na.omit(x))/sqrt(length(na.omit(x))))
# Visualization: including individual weighting methods   
plotdata1 <- data.frame(rbind(mean.indw, mean.test, mean.cv))
plotdata2 <- data.frame(rbind(se.indw, se.test, se.cv))
colnames(plotdata1) <- as.character(samplesize)
plotdata1
plotdata2
print(min(plotdata1))
print(max(plotdata1))
range <- c(-18, 0) # for unemployment rate data
stepsize <- 2
# plot
setEPS()
postscript("PMSE_size_equal.eps", width = 4, height = 4)
bars <- barplot(as.matrix(plotdata1*100), ylim = c(range[1], range[2]), yaxt = "n", 
                xlab = "Sample Size", ylab = "%MSE Reduction to SA", density = c(90, 45, 0), beside = TRUE)
axis(2, seq(range[1], range[2], by = stepsize), labels = paste(seq(range[1], range[2], by = stepsize), "%", sep = ""), las = 1)
segments(bars, as.vector(as.matrix(plotdata1))*100 - as.vector(as.matrix(plotdata2))*300, 
         bars, as.vector(as.matrix(plotdata1))*100 + as.vector(as.matrix(plotdata2))*300, lwd = 2)
arrows(bars, as.vector(as.matrix(plotdata1))*100 - as.vector(as.matrix(plotdata2))*300, 
       bars, as.vector(as.matrix(plotdata1))*100 + as.vector(as.matrix(plotdata2))*300, lwd = 2, angle = 90, code = 3, length = 0.05)
legend("bottomright", c("WA", "Test", "CV"), density = c(90, 45, 0), bg = "white", bty = "n")
dev.off()

#### FAR/HR table ####
if (!require('data.table')) install.packages('data.table'); library('data.table')
if (!require('mltools')) install.packages('mltools'); library('mltools')
#### FAR/HR table: for each weighting model and each sample size ####
n <- c(16, 32, 64, 128, 256)
weight.name <- c("OW", "TOP3", "RP", "CWM", "SSIN", "SSDE", "REGCOV", "LAS")
crit.outtest <- 0.001
crit <- 0.05
FAR <- array(NA, dim = c(length(n), length(weight.name), 2))
HR <- array(NA, dim = c(length(n), length(weight.name), 2))
CaseCount <- array(NA, dim = c(length(n), length(weight.name), 3))
for(i in 1:length(n)){
  subdat <- dat[which(dat$samplesize == n[i]), ]
  pvalue.i <- subdat[ , 3:10]
  cverror.i <- subdat[ , 11:18]
  MSEout.i <- subdat[ , 19:27]
  MSEout.ttest.i <- subdat[ , 28:35]
  MSEerror.i <- MSEout.i - MSEout.i[ , 9] # error(wa) - error(sa)
  MSEerror.i <- MSEerror.i[ , -9]
  MSEsign.i <- MSEout.ttest.i <= crit.outtest
  for(k in 1:length(weight.name)){
    pvalue.vec <- as.vector(pvalue.i[ ,k])
    cv.vec <- as.vector(cverror.i[ ,k]) # wa wins when cv is negative 
    sign.index <- which(MSEsign.i[ ,k] == TRUE)
    MSE.vec <- rep(0, nrow(MSEerror.i))
    MSE.vec[sign.index] <- as.vector(MSEerror.i[sign.index, k]) # wa wins when MSE is negative 
    tempdata <- cbind(pvalue.vec, cv.vec, MSE.vec)
    remove.index <- which(pvalue.vec == "No need to compare!")
    if(length(remove.index) > 0){
      tempdata <- tempdata[-remove.index, ]
    }
    # find the FAR and HR
    CaseCount[i, k, ] <- c(length(which(tempdata[ , 3] == 0)), # inconclusice
                           length(which(tempdata[ , 3] > 0)), # SA wins significantly
                           length(which(tempdata[ , 3] < 0))) # WA wins significantly
    NullNull_1 <- length(which(tempdata[ , 1] >= crit & tempdata[ , 3] > 0))
    NullNull_2 <- length(which(tempdata[ , 2] >= 0 & tempdata[ , 3] > 0))
    NullReject_1 <- length(which(tempdata[ , 1] < crit & tempdata[ , 3] > 0))
    NullReject_2 <- length(which(tempdata[ , 2] < 0 & tempdata[ , 3] > 0))
    FAR[i, k, ] <- c(NullReject_1/(NullReject_1 + NullNull_1), 
                     NullReject_2/(NullReject_2 + NullNull_2))
    AltRetain_1 <- length(which(tempdata[ , 1] >= crit & tempdata[ , 3] < 0))
    AltRetain_2 <- length(which(tempdata[ , 2] >= 0 & tempdata[ , 3] < 0))
    AltAlt_1 <- length(which(tempdata[ , 1] < crit & tempdata[ , 3] < 0))
    AltAlt_2 <- length(which(tempdata[ , 2] < 0 & tempdata[ , 3] < 0))
    HR[i, k, ] <- c(AltAlt_1/(AltAlt_1 + AltRetain_1), 
                    AltAlt_2/(AltAlt_2 + AltRetain_2))
  }
}
size_index <- 5
CaseCount[ size_index, , ]

t(FAR[ , , 1]*CaseCount[ , , 2])
t(FAR[ , , 2]*CaseCount[ , , 2])
t((1-HR[ , , 1])*CaseCount[ , , 3])
t((1-HR[ , , 2])*CaseCount[ , , 3])

##### FAR/HR table: aggregating for each weighting model ##### 
n <- c(16, 32, 64, 128, 256)
## loading all related data
pvalue <- dat[ , 3:10]
cverror <- dat[ , 11:18]
MSEout <- dat[ , 19:27]
MSEout.ttest <- dat[ , 28:35]
MSEerror <- MSEout - MSEout[ , 9] # error(wa) - error(sa)
MSEerror <- MSEerror[ , -9]
MSEsign <- MSEout.ttest <= 0.001
## find FAR, HR and AUC
weight.name <- c("OW", "TOP3", "RP", "CWM", "SSIN", "SSDE", "REGCOV", "LAS")
crit <- 0.05
FAR <- matrix(NA, nrow = 2, ncol = length(weight.name))
HR <- matrix(NA, nrow = 2, ncol = length(weight.name))
# results tables 
for(k in 1:ncol(pvalue)){
  pvalue.vec <- as.vector(pvalue[ ,k])
  cv.vec <- as.vector(cverror[ ,k]) # wa wins when cv is negative 
  sign.index <- which(MSEsign[ ,k] == TRUE)
  MSE.vec <- rep(0, nrow(MSEerror))
  MSE.vec[sign.index] <- as.vector(MSEerror[sign.index, k]) # wa wins when MSE is negative 
  tempdata <- cbind(pvalue.vec, cv.vec, MSE.vec)
  remove.index <- which(pvalue.vec == "No need to compare!")
  if(length(remove.index) > 0){
    tempdata <- tempdata[-remove.index, ]
  }
  # find the FAR and HR
  NullNull_1 <- length(which(tempdata[ , 1] >= crit & tempdata[ , 3] > 0))
  NullNull_2 <- length(which(tempdata[ , 2] >= 0 & tempdata[ , 3] > 0))
  NullReject_1 <- length(which(tempdata[ , 1] < crit & tempdata[ , 3] > 0))
  NullReject_2 <- length(which(tempdata[ , 2] < 0 & tempdata[ , 3] > 0))
  FAR[ , k] <- c(NullReject_1/(NullReject_1 + NullNull_1), 
                 NullReject_2/(NullReject_2 + NullNull_2))
  AltRetain_1 <- length(which(tempdata[ , 1] >= crit & tempdata[ , 3] < 0))
  AltRetain_2 <- length(which(tempdata[ , 2] >= 0 & tempdata[ , 3] < 0))
  AltAlt_1 <- length(which(tempdata[ , 1] < crit & tempdata[ , 3] < 0))
  AltAlt_2 <- length(which(tempdata[ , 2] < 0 & tempdata[ , 3] < 0))
  HR[ , k] <- c(AltAlt_1/(AltAlt_1 + AltRetain_1), 
                AltAlt_2/(AltAlt_2 + AltRetain_2))
}
# print out the pooled results: 
round(FAR, 3)
round(HR, 3)

###### FAR/HR table: aggregating for each sample size ########
crit.outtest <- 0.001
n.vec <- c(16, 32, 64, 128, 256)
crit <- 0.05
for(j in 1:length(n.vec)){
  n <- n.vec[j]
  subdat <- dat[which(dat$samplesize == n), ]
  pvalue <- subdat[ , 3:10]
  cverror <- subdat[ , 11:18]
  MSEout <- subdat[ , 19:27]
  MSEout.ttest <- subdat[ , 28:35]
  MSEerror <- MSEout - MSEout[ , 9] # error(wa) - error(sa)
  MSEerror <- MSEerror[ , -9]
  MSEsign <- MSEout.ttest <= crit.outtest
  ## find FAR, HR and AUC
  FAR <- c()
  HR <- c()
  pvalue.vec <- as.vector(as.matrix(pvalue[ , -7]))
  cv.vec <- as.vector(as.matrix(cverror[ , -7]))
  sign.index <- which(as.vector(as.matrix(MSEsign[ , -7])) == TRUE)
  MSE.vec <- rep(0, length(pvalue.vec))
  MSEerror.vec <- as.vector(as.matrix(MSEerror[ , -7]))
  MSE.vec[sign.index] <- MSEerror.vec[sign.index] # wa wins when MSE is negative 
  tempdata <- cbind(pvalue.vec, cv.vec, MSE.vec)
  remove.index <- which(as.character(pvalue.vec) == "No need to compare!")
  if(length(remove.index) > 0){
    tempdata <- tempdata[-remove.index, ]
  }
  # find the FAR
  NullNull_1 <- length(which(tempdata[ , 1] > crit & tempdata[ , 3] > 0))
  NullNull_2 <- length(which(tempdata[ , 2] > 0 & tempdata[ , 3] > 0))
  NullReject_1 <- length(which(tempdata[ , 1] < crit & tempdata[ , 3] > 0))
  NullReject_2 <- length(which(tempdata[ , 2] < 0 & tempdata[ , 3] > 0))
  FAR <- c(NullReject_1/(NullReject_1 + NullNull_1),
           NullReject_2/(NullReject_2 + NullNull_2))
  AltRetain_1 <- length(which(tempdata[ , 1] > crit & tempdata[ , 3] < 0))
  AltRetain_2 <- length(which(tempdata[ , 2] > 0 & tempdata[ , 3] < 0))
  AltAlt_1 <- length(which(tempdata[ , 1] < crit & tempdata[ , 3] < 0))
  AltAlt_2 <- length(which(tempdata[ , 2] < 0 & tempdata[ , 3] < 0))
  HR <- c(AltAlt_1/(AltAlt_1 + AltRetain_1), 
          AltAlt_2/(AltAlt_2 + AltRetain_2))
  print(round(FAR, 3))
  print(round(HR, 3))
}
#########################################################
################ Simulation: SPF DATA ###################
#########################################################
setwd("insert the path to the related results data. ")
library(ggplot2)
library(dplyr)
library(forcats)
library(viridis)
## Percentage of MSE reduction 
modelnames <- c("OW", "TOP3", "RP", "CWM", "SSIN", "SSDE", "REGCOV", "LAS", "SA")
samplesize <- c(16, 32, 64, 128, 256)
dat <- read.csv("InfStates.csv", header = T, stringsAsFactors = F)
#dat <- read.csv("UnempStates.csv", header = T, stringsAsFactors = F)
colnames(dat) <- c("samplesize", "round", "p1", "p2", "p3", "p4", "p5", "p6", "p7", "p8",
                   "cv1", "cv2", "cv3", "cv4", "cv5", "cv6", "cv7", "cv8",
                   "MSEout1", "MSEout2", "MSEout3", "MSEout4", "MSEout5", "MSEout6", "MSEout7", "MSEout8", "MSEout9",
                   "MSEttest1", "MSEttest2", "MSEttest3", "MSEttest4", "MSEttest5", "MSEttest6", "MSEttest7", "MSEttest8")
head(dat)
# revised optimal weighting method 
dat_ow <- read.csv("InfStates_OW.csv", header = T, stringsAsFactors = F)
#dat_ow <- read.csv("UnempStates_OW.csv", header = T, stringsAsFactors = F)
colnames(dat_ow) <- c("samplesize", "round", "p1", "p2", "cv1", "cv2", "MSEout1", "MSEout2", "MSEout3", "MSEttest1", "MSEttest2")
head(dat_ow)
#### MSE Reduction: aggregating for each weighting method ####
pvalue <- dat[ , 3:10]
cverror <- dat[ , 11:18]
MSEmean <- dat[ , 19:27]
PMSEtoSA <- (MSEmean[ , 9] - MSEmean)/MSEmean[ , 9]
PMSEdata.test <- matrix(NA, nrow = nrow(PMSEtoSA), ncol = ncol(PMSEtoSA) - 1)
PMSEdata.cv   <- matrix(NA, nrow = nrow(PMSEtoSA), ncol = ncol(PMSEtoSA) - 1)
crit <- 0.05
for(k in 1:(length(modelnames) - 1)){
  test_select_wa_index <- which(pvalue[ , k] <= crit)
  test_select_sa_index <- which(pvalue[ , k] > crit)
  cv_select_wa_index   <- which(cverror[ , k] <= 0)
  cv_select_sa_index   <- which(cverror[ , k] > 0)
  PMSEdata.test[test_select_wa_index, k] <- PMSEtoSA[test_select_wa_index, k]
  PMSEdata.test[test_select_sa_index, k] <- PMSEtoSA[test_select_sa_index, 9]
  PMSEdata.cv[cv_select_wa_index, k]     <- PMSEtoSA[cv_select_wa_index, k]
  PMSEdata.cv[cv_select_sa_index, k]     <- PMSEtoSA[cv_select_sa_index, 9]
}
# remove the outliers for RegCov weighting method 
PMSEtoSA <- PMSEtoSA[ , -9]
mean.indw <- apply(PMSEtoSA[ , -7], 2, function(x) mean(na.omit(x)))
mean.test <- apply(PMSEdata.test[ , -7], 2, function(x) mean(na.omit(x)))
mean.cv <- apply(PMSEdata.cv[ , -7], 2, function(x) mean(na.omit(x)))
se.indw <- apply(PMSEtoSA[ , -7], 2, function(x) sd(na.omit(x))/sqrt(length(na.omit(x))))
se.test <- apply(PMSEdata.test[ , -7], 2, function(x) sd(na.omit(x))/sqrt(length(na.omit(x))))
se.cv <- apply(PMSEdata.cv[ , -7], 2, function(x) sd(na.omit(x))/sqrt(length(na.omit(x))))
#Using more accurate optimal weights without the non-negative constraint using optw_sum1 function 
pvalue_ow <- dat_ow[ , 3:4]
cverror_ow <- dat_ow[ , 5:6]
MSEmean_ow <- dat_ow[ , 7:9]
PMSEtoSA_ow <- (MSEmean_ow[ , 3] - MSEmean_ow)/MSEmean_ow[ , 3]
PMSEdata.test_ow <- matrix(NA, nrow = nrow(PMSEtoSA_ow), ncol = ncol(PMSEtoSA_ow) - 1)
PMSEdata.cv_ow   <- matrix(NA, nrow = nrow(PMSEtoSA_ow), ncol = ncol(PMSEtoSA_ow) - 1)
crit <- 0.05
for(k in 1:2){
  test_select_wa_index <- which(pvalue_ow[ , k] <= crit)
  test_select_sa_index <- which(pvalue_ow[ , k] > crit)
  cv_select_wa_index   <- which(cverror_ow[ , k] <= 0)
  cv_select_sa_index   <- which(cverror_ow[ , k] > 0)
  PMSEdata.test_ow[test_select_wa_index, k] <- PMSEtoSA_ow[test_select_wa_index, k]
  PMSEdata.test_ow[test_select_sa_index, k] <- PMSEtoSA_ow[test_select_sa_index, 3]
  PMSEdata.cv_ow[cv_select_wa_index, k]     <- PMSEtoSA_ow[cv_select_wa_index, k]
  PMSEdata.cv_ow[cv_select_sa_index, k]     <- PMSEtoSA_ow[cv_select_sa_index, 3]
}
# remove the outliers for RegCov weighting method 
PMSEtoSA_ow <- PMSEtoSA_ow[ , -3]
mean.indw_ow <- apply(PMSEtoSA_ow, 2, function(x) mean(na.omit(x)))
mean.test_ow <- apply(PMSEdata.test_ow, 2, function(x) mean(na.omit(x)))
mean.cv_ow <- apply(PMSEdata.cv_ow, 2, function(x) mean(na.omit(x)))
se.indw_ow <- apply(PMSEtoSA_ow, 2, function(x) sd(na.omit(x))/sqrt(length(na.omit(x))))
se.test_ow <- apply(PMSEdata.test_ow, 2, function(x) sd(na.omit(x))/sqrt(length(na.omit(x))))
se.cv_ow <- apply(PMSEdata.cv_ow, 2, function(x) sd(na.omit(x))/sqrt(length(na.omit(x))))
# Visualization: including individual weighting methods 
plotdata1 <- data.frame(rbind(mean.indw, mean.test, mean.cv))
plotdata2 <- data.frame(rbind(se.indw, se.test, se.cv))
colnames(plotdata1) <- modelnames[c(1:6, 8)]
plotdata1 
plotdata2
plotdata1[ , 1] <- c(mean.indw_ow[1], mean.test_ow[1], mean.cv_ow[1])
plotdata2[ , 1] <- c(se.indw_ow[1], se.test_ow[1], se.cv_ow[1])
print(min(plotdata1, na.rm = T))
print(max(plotdata1, na.rm = T))
range <- c(0, 10) # for unemployment rate data
stepsize <- 2
#plotdata1[1, 7] <- min(range)*0.01

# plot
setEPS()
postscript("PMSE_model_inf.eps", width = 7, height = 3.5)
bars <- barplot(as.matrix(plotdata1*100), ylim = c(range[1], range[2]), yaxt = "n", 
                xlab = "Weighting Model", ylab = "%MSE Reduction to SA", density = c(90, 45, 0), beside = TRUE)
axis(2, seq(range[1], range[2], by = stepsize), labels = paste(seq(range[1], range[2], by = stepsize), "%", sep = ""), las = 1)
segments(bars, as.vector(as.matrix(plotdata1))*100 - as.vector(as.matrix(plotdata2))*300, 
         bars, as.vector(as.matrix(plotdata1))*100 + as.vector(as.matrix(plotdata2))*300, lwd = 2)
arrows(bars, as.vector(as.matrix(plotdata1))*100 - as.vector(as.matrix(plotdata2))*300, 
       bars, as.vector(as.matrix(plotdata1))*100 + as.vector(as.matrix(plotdata2))*300, lwd = 2, angle = 90, code = 3, length = 0.05)
legend("topleft", c("WA", "Test", "CV"), density = c(90, 45, 0), bg = "white", bty = "n")
dev.off()

#### MSE Reduction: aggregating for each sample size ####
table(dat$samplesize)
table(dat_ow$samplesize)
samplesize <- c(16, 32, 64, 128, 256)
sim_runs <- 925
PMSEdata.indw <- c()
PMSEdata.test <- c()
PMSEdata.cv   <- c()
crit <- 0.05
for(i in 1:length(samplesize)){
  subdat <- dat[which(dat$samplesize == samplesize[i]), ]
  ## load the data 
  MSEmean.i <- subdat[ , 19:27]
  PMSEtoSA.i <- (MSEmean.i[ , 9] - MSEmean.i)/MSEmean.i[ , 9]
  cverror.i <- subdat[ , 11:18]
  pvalue.i <- subdat[ , 3:10] 
  ## find the MSE for test and cv
  model_k <- 1:8
  PMSEdata.test.i <- matrix(NA, nrow = nrow(PMSEtoSA.i), ncol = length(model_k))
  PMSEdata.cv.i   <- matrix(NA, nrow = nrow(PMSEtoSA.i), ncol = length(model_k))
  for(k in model_k){
    test_select_wa_index <- which(pvalue.i[ , k] <= crit)
    test_select_sa_index <- which(pvalue.i[ , k] > crit)
    cv_select_wa_index   <- which(cverror.i[ , k] <= 0)
    cv_select_sa_index   <- which(cverror.i[ , k] > 0)
    PMSEdata.test.i[test_select_wa_index, k] <- PMSEtoSA.i[test_select_wa_index, k]
    PMSEdata.test.i[test_select_sa_index, k] <- PMSEtoSA.i[test_select_sa_index, 9]
    PMSEdata.cv.i[cv_select_wa_index, k]     <- PMSEtoSA.i[cv_select_wa_index, k]
    PMSEdata.cv.i[cv_select_sa_index, k]     <- PMSEtoSA.i[cv_select_sa_index, 9]
  }
  PMSEtoSA.i <- PMSEtoSA.i[ , -9]
  PMSEdata.test.i <- PMSEdata.test.i[ , -c(1, 7)]
  PMSEdata.cv.i <- PMSEdata.test.i[ , -c(1, 7)]
  PMSEtoSA.i <- PMSEtoSA.i[ , -c(1, 7)]
  
  subdat_ow <- dat_ow[which(dat_ow$samplesize == samplesize[i]), ]
  ## load the data for revised optimal weights 
  MSEmean.i_ow <- subdat_ow[ , 7:9]
  PMSEtoSA.i_ow <- (MSEmean.i_ow[ , 3] - MSEmean.i_ow)/MSEmean.i_ow[ , 3]
  cverror.i_ow <- subdat_ow[ , 5:6]
  pvalue.i_ow <- subdat_ow[ , 3:4] 
  ## find the MSE for test and cv
  model_k <- 1:2
  PMSEdata.test.i_ow <- matrix(NA, nrow = nrow(PMSEtoSA.i_ow), ncol = length(model_k))
  PMSEdata.cv.i_ow   <- matrix(NA, nrow = nrow(PMSEtoSA.i_ow), ncol = length(model_k))
  for(k in model_k){
    test_select_wa_index <- which(pvalue.i_ow[ , k] <= crit)
    test_select_sa_index <- which(pvalue.i_ow[ , k] > crit)
    cv_select_wa_index   <- which(cverror.i_ow[ , k] <= 0)
    cv_select_sa_index   <- which(cverror.i_ow[ , k] > 0)
    PMSEdata.test.i_ow[test_select_wa_index, k] <- PMSEtoSA.i_ow[test_select_wa_index, k]
    PMSEdata.test.i_ow[test_select_sa_index, k] <- PMSEtoSA.i_ow[test_select_sa_index, 3]
    PMSEdata.cv.i_ow[cv_select_wa_index, k]     <- PMSEtoSA.i_ow[cv_select_wa_index, k]
    PMSEdata.cv.i_ow[cv_select_sa_index, k]     <- PMSEtoSA.i_ow[cv_select_sa_index, 3]
  }
  PMSEtoSA.i_ow <- PMSEtoSA.i_ow[ , -3]
  # combine all weights 
  PMSEtoSA.i <- cbind(PMSEtoSA.i_ow[1:nrow(PMSEtoSA.i), 1], PMSEtoSA.i)
  PMSEdata.test.i <- cbind(PMSEdata.test.i_ow[1:nrow(PMSEdata.test.i), 1], PMSEdata.test.i)
  PMSEdata.cv.i <- cbind(PMSEdata.cv.i_ow[1:nrow(PMSEdata.cv.i), 1], PMSEdata.cv.i)
  # 
  PMSEdata.indw <- cbind(PMSEdata.indw, as.vector(as.matrix(PMSEtoSA.i[1:sim_runs, ])))
  PMSEdata.test <- cbind(PMSEdata.test, as.vector(PMSEdata.test.i[1:sim_runs, ]))
  PMSEdata.cv <- cbind(PMSEdata.cv, as.vector(PMSEdata.cv.i[1:sim_runs, ]))
}
mean.indw <- apply(PMSEdata.indw, 2, function(x) mean(na.omit(x)))
se.indw <- apply(PMSEdata.indw, 2, function(x) sd(na.omit(x))/sqrt(length(na.omit(x))))
mean.test <- apply(PMSEdata.test, 2, function(x) mean(na.omit(x)))
se.test <- apply(PMSEdata.test, 2, function(x) sd(na.omit(x))/sqrt(length(na.omit(x))))
mean.cv <- apply(PMSEdata.cv, 2, function(x) mean(na.omit(x)))
se.cv <- apply(PMSEdata.cv, 2, function(x) sd(na.omit(x))/sqrt(length(na.omit(x))))
# Visualization: including individual weighting methods   
plotdata1 <- data.frame(rbind(mean.indw, mean.test, mean.cv))
plotdata2 <- data.frame(rbind(se.indw, se.test, se.cv))
plotdata1
plotdata2
colnames(plotdata1) <- as.character(samplesize)
print(min(plotdata1))
print(max(plotdata1))
#range <- c(-20, 5) # for inflation rate data
range <- c(-4, 8) # for unemployment rate data
stepsize <- 2
# plot
setEPS()
postscript("PMSE_size_inf.eps", width = 4, height = 4)
bars <- barplot(as.matrix(plotdata1*100), ylim = c(range[1], range[2]), yaxt = "n", 
                xlab = "Sample Size", ylab = "%MSE Reduction to SA", density = c(90, 45, 0), beside = TRUE)
axis(2, seq(range[1], range[2], by = stepsize), labels = paste(seq(range[1], range[2], by = stepsize), "%", sep = ""), las = 1)
segments(bars, as.vector(as.matrix(plotdata1))*100 - as.vector(as.matrix(plotdata2))*300, 
         bars, as.vector(as.matrix(plotdata1))*100 + as.vector(as.matrix(plotdata2))*300, lwd = 2)
arrows(bars, as.vector(as.matrix(plotdata1))*100 - as.vector(as.matrix(plotdata2))*300, 
       bars, as.vector(as.matrix(plotdata1))*100 + as.vector(as.matrix(plotdata2))*300, lwd = 2, angle = 90, code = 3, length = 0.05)
legend("topleft", c("WA", "Test", "CV"), density = c(90, 45, 0), bg = "white", bty = "n")
dev.off()

#### FAR/HR table ####
if (!require('data.table')) install.packages('data.table'); library('data.table')
if (!require('mltools')) install.packages('mltools'); library('mltools')
#### FAR/HR table: for each weighting model and each sample size ####
n <- c(16, 32, 64, 128, 256)
weight.name <- c("OW", "TOP3", "RP", "CWM", "SSIN", "SSDE", "REGCOV", "LAS")
crit.outtest <- 0.001
crit <- 0.05
FAR <- array(NA, dim = c(length(n), length(weight.name), 2))
HR <- array(NA, dim = c(length(n), length(weight.name), 2))
CaseCount <- array(NA, dim = c(length(n), length(weight.name), 3))
for(i in 1:length(n)){
  subdat <- dat[which(dat$samplesize == n[i]), ]
  pvalue.i <- subdat[ , 3:10]
  cverror.i <- subdat[ , 11:18]
  MSEout.i <- subdat[ , 19:27]
  MSEout.ttest.i <- subdat[ , 28:35]
  MSEerror.i <- MSEout.i - MSEout.i[ , 9] # error(wa) - error(sa)
  MSEerror.i <- MSEerror.i[ , -9]
  MSEsign.i <- MSEout.ttest.i <= crit.outtest
  for(k in 1:length(weight.name)){
    pvalue.vec <- as.vector(pvalue.i[ ,k])
    cv.vec <- as.vector(cverror.i[ ,k]) # wa wins when cv is negative 
    sign.index <- which(MSEsign.i[ ,k] == TRUE)
    MSE.vec <- rep(0, nrow(MSEerror.i))
    MSE.vec[sign.index] <- as.vector(MSEerror.i[sign.index, k]) # wa wins when MSE is negative 
    tempdata <- cbind(pvalue.vec, cv.vec, MSE.vec)
    remove.index <- which(pvalue.vec == "No need to compare!")
    if(length(remove.index) > 0){
      tempdata <- tempdata[-remove.index, ]
    }
    # find the FAR and HR
    CaseCount[i, k, ] <- c(length(which(tempdata[ , 3] == 0)), # inconclusice
                           length(which(tempdata[ , 3] > 0)), # SA wins significantly
                           length(which(tempdata[ , 3] < 0))) # WA wins significantly
    NullNull_1 <- length(which(tempdata[ , 1] >= crit & tempdata[ , 3] > 0))
    NullNull_2 <- length(which(tempdata[ , 2] >= 0 & tempdata[ , 3] > 0))
    NullReject_1 <- length(which(tempdata[ , 1] < crit & tempdata[ , 3] > 0))
    NullReject_2 <- length(which(tempdata[ , 2] < 0 & tempdata[ , 3] > 0))
    FAR[i, k, ] <- c(NullReject_1/(NullReject_1 + NullNull_1), 
                     NullReject_2/(NullReject_2 + NullNull_2))
    AltRetain_1 <- length(which(tempdata[ , 1] >= crit & tempdata[ , 3] < 0))
    AltRetain_2 <- length(which(tempdata[ , 2] >= 0 & tempdata[ , 3] < 0))
    AltAlt_1 <- length(which(tempdata[ , 1] < crit & tempdata[ , 3] < 0))
    AltAlt_2 <- length(which(tempdata[ , 2] < 0 & tempdata[ , 3] < 0))
    HR[i, k, ] <- c(AltAlt_1/(AltAlt_1 + AltRetain_1), 
                    AltAlt_2/(AltAlt_2 + AltRetain_2))
  }
}
size_index <- 5
CaseCount[ size_index, , ]

t(FAR[ , , 1]*CaseCount[ , , 2])
t(FAR[ , , 2]*CaseCount[ , , 2])
t((1-HR[ , , 1])*CaseCount[ , , 3])
t((1-HR[ , , 2])*CaseCount[ , , 3])

# for revised optimal weights 
weight.name <- c("OW0", "OW1")
crit.outtest <- 0.001
crit <- 0.05
FAR_ow <- array(NA, dim = c(length(n), length(weight.name), 2))
HR_ow <- array(NA, dim = c(length(n), length(weight.name), 2))
CaseCount_ow <- array(NA, dim = c(length(n), length(weight.name), 3))
for(i in 1:length(n)){
  subdat <- dat_ow[which(dat_ow$samplesize == n[i]), ]
  pvalue.i <- subdat[ , 3:4]
  cverror.i <- subdat[ , 5:6]
  MSEout.i <- subdat[ , 7:9]
  MSEout.ttest.i <- subdat[ , 10:11]
  MSEerror.i <- MSEout.i - MSEout.i[ , 3] # error(wa) - error(sa)
  MSEerror.i <- MSEerror.i[ , -3]
  MSEsign.i <- MSEout.ttest.i <= crit.outtest
  for(k in 1:length(weight.name)){
    pvalue.vec <- as.vector(pvalue.i[ ,k])
    cv.vec <- as.vector(cverror.i[ ,k]) # wa wins when cv is negative 
    sign.index <- which(MSEsign.i[ ,k] == TRUE)
    MSE.vec <- rep(0, nrow(MSEerror.i))
    MSE.vec[sign.index] <- as.vector(MSEerror.i[sign.index, k]) # wa wins when MSE is negative 
    tempdata <- cbind(pvalue.vec, cv.vec, MSE.vec)
    remove.index <- which(pvalue.vec == "No need to compare!")
    if(length(remove.index) > 0){
      tempdata <- tempdata[-remove.index, ]
    }
    # find the FAR and HR
    CaseCount_ow[i, k, ] <- c(length(which(tempdata[ , 3] == 0)), # inconclusice
                           length(which(tempdata[ , 3] > 0)), # SA wins significantly
                           length(which(tempdata[ , 3] < 0))) # WA wins significantly
    NullNull_1 <- length(which(tempdata[ , 1] >= crit & tempdata[ , 3] > 0))
    NullNull_2 <- length(which(tempdata[ , 2] >= 0 & tempdata[ , 3] > 0))
    NullReject_1 <- length(which(tempdata[ , 1] < crit & tempdata[ , 3] > 0))
    NullReject_2 <- length(which(tempdata[ , 2] < 0 & tempdata[ , 3] > 0))
    FAR_ow[i, k, ] <- c(NullReject_1/(NullReject_1 + NullNull_1), 
                     NullReject_2/(NullReject_2 + NullNull_2))
    AltRetain_1 <- length(which(tempdata[ , 1] >= crit & tempdata[ , 3] < 0))
    AltRetain_2 <- length(which(tempdata[ , 2] >= 0 & tempdata[ , 3] < 0))
    AltAlt_1 <- length(which(tempdata[ , 1] < crit & tempdata[ , 3] < 0))
    AltAlt_2 <- length(which(tempdata[ , 2] < 0 & tempdata[ , 3] < 0))
    HR_ow[i, k, ] <- c(AltAlt_1/(AltAlt_1 + AltRetain_1), 
                    AltAlt_2/(AltAlt_2 + AltRetain_2))
  }
}
size_index <- 5
CaseCount_ow[ size_index, , ]

t(FAR_ow[ , , 1]*CaseCount_ow[ , , 2])
t(FAR_ow[ , , 2]*CaseCount_ow[ , , 2])
t((1-HR_ow[ , , 1])*CaseCount_ow[ , , 3])
t((1-HR_ow[ , , 2])*CaseCount_ow[ , , 3])

##### FAR/HR table: aggregating for each weighting model ##### 
n <- c(16, 32, 64, 128, 256)
## loading all related data
pvalue <- dat[ , 3:10]
cverror <- dat[ , 11:18]
MSEout <- dat[ , 19:27]
MSEout.ttest <- dat[ , 28:35]
MSEerror <- MSEout - MSEout[ , 9] # error(wa) - error(sa)
MSEerror <- MSEerror[ , -9]
MSEsign <- MSEout.ttest <= 0.001
## find FAR, HR and AUC
weight.name <- c("OW", "TOP3", "RP", "CWM", "SSIN", "SSDE", "REGCOV", "LAS")
crit <- 0.05
FAR <- matrix(NA, nrow = 2, ncol = length(weight.name))
HR <- matrix(NA, nrow = 2, ncol = length(weight.name))
# results tables 
for(k in 1:ncol(pvalue)){
  pvalue.vec <- as.vector(pvalue[ ,k])
  cv.vec <- as.vector(cverror[ ,k]) # wa wins when cv is negative 
  sign.index <- which(MSEsign[ ,k] == TRUE)
  MSE.vec <- rep(0, nrow(MSEerror))
  MSE.vec[sign.index] <- as.vector(MSEerror[sign.index, k]) # wa wins when MSE is negative 
  tempdata <- cbind(pvalue.vec, cv.vec, MSE.vec)
  remove.index <- which(pvalue.vec == "No need to compare!")
  if(length(remove.index) > 0){
    tempdata <- tempdata[-remove.index, ]
  }
  # find the FAR and HR
  NullNull_1 <- length(which(tempdata[ , 1] >= crit & tempdata[ , 3] > 0))
  NullNull_2 <- length(which(tempdata[ , 2] >= 0 & tempdata[ , 3] > 0))
  NullReject_1 <- length(which(tempdata[ , 1] < crit & tempdata[ , 3] > 0))
  NullReject_2 <- length(which(tempdata[ , 2] < 0 & tempdata[ , 3] > 0))
  FAR[ , k] <- c(NullReject_1/(NullReject_1 + NullNull_1), 
                 NullReject_2/(NullReject_2 + NullNull_2))
  AltRetain_1 <- length(which(tempdata[ , 1] >= crit & tempdata[ , 3] < 0))
  AltRetain_2 <- length(which(tempdata[ , 2] >= 0 & tempdata[ , 3] < 0))
  AltAlt_1 <- length(which(tempdata[ , 1] < crit & tempdata[ , 3] < 0))
  AltAlt_2 <- length(which(tempdata[ , 2] < 0 & tempdata[ , 3] < 0))
  HR[ , k] <- c(AltAlt_1/(AltAlt_1 + AltRetain_1), 
                              AltAlt_2/(AltAlt_2 + AltRetain_2))
}
# print out the pooled results: 
round(FAR, 3)
round(HR, 3)
# for revised optimal weights 
## loading all related data
pvalue <- dat_ow[ , 3:4]
cverror <- dat_ow[ , 5:6]
MSEout <- dat_ow[ , 7:9]
MSEout.ttest <- dat_ow[ , 10:11]
MSEerror <- MSEout - MSEout[ , 3] # error(wa) - error(sa)
MSEerror <- MSEerror[ , -3]
MSEsign <- MSEout.ttest <= 0.001
## find FAR, HR and AUC
weight.name <- c("OW0", "OW1")
crit <- 0.05
FAR_ow <- matrix(NA, nrow = 2, ncol = length(weight.name))
HR_ow <- matrix(NA, nrow = 2, ncol = length(weight.name))
# results tables 
for(k in 1:ncol(pvalue)){
  pvalue.vec <- as.vector(pvalue[ ,k])
  cv.vec <- as.vector(cverror[ ,k]) # wa wins when cv is negative 
  sign.index <- which(MSEsign[ ,k] == TRUE)
  MSE.vec <- rep(0, nrow(MSEerror))
  MSE.vec[sign.index] <- as.vector(MSEerror[sign.index, k]) # wa wins when MSE is negative 
  tempdata <- cbind(pvalue.vec, cv.vec, MSE.vec)
  remove.index <- which(pvalue.vec == "No need to compare!")
  if(length(remove.index) > 0){
    tempdata <- tempdata[-remove.index, ]
  }
  # find the FAR and HR
  NullNull_1 <- length(which(tempdata[ , 1] >= crit & tempdata[ , 3] > 0))
  NullNull_2 <- length(which(tempdata[ , 2] >= 0 & tempdata[ , 3] > 0))
  NullReject_1 <- length(which(tempdata[ , 1] < crit & tempdata[ , 3] > 0))
  NullReject_2 <- length(which(tempdata[ , 2] < 0 & tempdata[ , 3] > 0))
  FAR_ow[ , k] <- c(NullReject_1/(NullReject_1 + NullNull_1), 
                 NullReject_2/(NullReject_2 + NullNull_2))
  AltRetain_1 <- length(which(tempdata[ , 1] >= crit & tempdata[ , 3] < 0))
  AltRetain_2 <- length(which(tempdata[ , 2] >= 0 & tempdata[ , 3] < 0))
  AltAlt_1 <- length(which(tempdata[ , 1] < crit & tempdata[ , 3] < 0))
  AltAlt_2 <- length(which(tempdata[ , 2] < 0 & tempdata[ , 3] < 0))
  HR_ow[ , k] <- c(AltAlt_1/(AltAlt_1 + AltRetain_1), 
                AltAlt_2/(AltAlt_2 + AltRetain_2))
}
# print out the pooled results: 
round(FAR_ow, 3)
round(HR_ow, 3)
###### FAR/HR table: aggregating for each sample size ########
crit.outtest <- 0.001
n.vec <- c(16, 32, 64, 128, 256)
crit <- 0.05
for(j in 1:length(n.vec)){
  n <- n.vec[j]
  subdat <- dat[which(dat$samplesize == n), ]
  subdat_ow <- dat_ow[which(dat_ow$samplesize == n), ]
  
  pvalue <- subdat[ , c(4:8, 10)]
  cverror <- subdat[ , c(12:16, 18)]
  MSEout <- subdat[ , c(20:24, 26:27)]
  MSEout.ttest <- subdat[ , c(29:33, 35)]
  MSEerror <- MSEout - MSEout[ , 7] # error(wa) - error(sa)
  MSEerror <- MSEerror[ , -7]
  MSEsign <- MSEout.ttest <= crit.outtest
  
  pvalue_ow <- subdat_ow[ , 3]
  cverror_ow <- subdat_ow[ , 5]
  MSEout_ow <- subdat_ow[ , c(7, 9)]
  MSEout.ttest_ow <- subdat_ow[ , 10]
  MSEerror_ow <- MSEout_ow - MSEout_ow[ , 2] # error(wa) - error(sa)
  MSEerror_ow <- MSEerror_ow[ , -2]
  MSEsign_ow <- MSEout.ttest_ow <= crit.outtest
  
  ## find FAR, HR and AUC
  FAR <- c()
  HR <- c()
  # for first 4 weighted averages
  # pvalue.vec <- c(pvalue[,1], pvalue[,2], pvalue[,6], pvalue[,7]) 
  # cv.vec <- c(cverror[,1],cverror[,2],cverror[,6],cverror[,7])
  pvalue.vec <- c(as.vector(as.matrix(pvalue)), as.vector(pvalue_ow))
  cv.vec <- c(as.vector(as.matrix(cverror)), as.vector(cverror_ow))
  sign.index <- which(c(as.vector(as.matrix(MSEsign)), as.vector(MSEsign_ow)) == TRUE)
  MSE.vec <- rep(0, length(pvalue.vec))
  # MSEerror.vec <- c(MSEerror[,1],MSEerror[,2],MSEerror[,6],MSEerror[,7])
  MSEerror.vec <- c(as.vector(as.matrix(MSEerror)), as.vector(MSEerror_ow))
  MSE.vec[sign.index] <- MSEerror.vec[sign.index] # wa wins when MSE is negative 
  tempdata <- cbind(pvalue.vec, cv.vec, MSE.vec)
  remove.index <- which(as.character(pvalue.vec) == "No need to compare!")
  if(length(remove.index) > 0){
    tempdata <- tempdata[-remove.index, ]
  }
  # find the FAR
  NullNull_1 <- length(which(tempdata[ , 1] > crit & tempdata[ , 3] > 0))
  NullNull_2 <- length(which(tempdata[ , 2] > 0 & tempdata[ , 3] > 0))
  NullReject_1 <- length(which(tempdata[ , 1] < crit & tempdata[ , 3] > 0))
  NullReject_2 <- length(which(tempdata[ , 2] < 0 & tempdata[ , 3] > 0))
  FAR <- c(NullReject_1/(NullReject_1 + NullNull_1),
           NullReject_2/(NullReject_2 + NullNull_2))
  AltRetain_1 <- length(which(tempdata[ , 1] > crit & tempdata[ , 3] < 0))
  AltRetain_2 <- length(which(tempdata[ , 2] > 0 & tempdata[ , 3] < 0))
  AltAlt_1 <- length(which(tempdata[ , 1] < crit & tempdata[ , 3] < 0))
  AltAlt_2 <- length(which(tempdata[ , 2] < 0 & tempdata[ , 3] < 0))
  HR <- c(AltAlt_1/(AltAlt_1 + AltRetain_1), 
          AltAlt_2/(AltAlt_2 + AltRetain_2))
  print(round(FAR, 3))
  print(round(HR, 3))
}

####################################################################
################ Simulation: COVID19 forecasting ###################
####################################################################
setwd("insert the path to the related results data. ")
library(ggplot2)
library(dplyr)
library(forcats)
library(viridis)
## Percentage of MSE reduction 
modelnames <- c("OW", "TOP3", "RP", "CWM", "SSIN", "SSDE", "REGCOV", "LAS", "SA")
samplesize <- c(16, 32, 64, 128, 256)
dat <- read.csv("CovidStates.csv", header = T, stringsAsFactors = F)
colnames(dat) <- c("target", "fips", "samplesize", "p1", "p2", "p3", "p4", "p5", "p6", "p7", "p8",
                   "cv1", "cv2", "cv3", "cv4", "cv5", "cv6", "cv7", "cv8",
                   "MSEout1", "MSEout2", "MSEout3", "MSEout4", "MSEout5", "MSEout6", "MSEout7", "MSEout8", "MSEout9",
                   "MSEttest1", "MSEttest2", "MSEttest3", "MSEttest4", "MSEttest5", "MSEttest6", "MSEttest7", "MSEttest8")
head(dat)
dat <- dat[ , -2]
#### MSE Reduction: aggregating for each weighting method ####
pvalue <- dat[ , 3:10]
cverror <- dat[ , 11:18]
MSEmean <- dat[ , 19:27]
PMSEtoSA <- (MSEmean[ , 9] - MSEmean)/MSEmean[ , 9]
PMSEdata.test <- matrix(NA, nrow = nrow(PMSEtoSA), ncol = ncol(PMSEtoSA) - 1)
PMSEdata.cv   <- matrix(NA, nrow = nrow(PMSEtoSA), ncol = ncol(PMSEtoSA) - 1)
crit <- 0.05
for(k in 1:(length(modelnames) - 1)){
  test_select_wa_index <- which(pvalue[ , k] <= crit)
  test_select_sa_index <- which(pvalue[ , k] > crit)
  cv_select_wa_index   <- which(cverror[ , k] <= 0)
  cv_select_sa_index   <- which(cverror[ , k] > 0)
  PMSEdata.test[test_select_wa_index, k] <- PMSEtoSA[test_select_wa_index, k]
  PMSEdata.test[test_select_sa_index, k] <- PMSEtoSA[test_select_sa_index, 9]
  PMSEdata.cv[cv_select_wa_index, k]     <- PMSEtoSA[cv_select_wa_index, k]
  PMSEdata.cv[cv_select_sa_index, k]     <- PMSEtoSA[cv_select_sa_index, 9]
}
# remove the outliers for RegCov weighting method 
PMSEtoSA <- PMSEtoSA[ , -9]
mean.indw <- apply(PMSEtoSA[ , -7], 2, function(x) mean(na.omit(x)))
mean.test <- apply(PMSEdata.test[ , -7], 2, function(x) mean(na.omit(x)))
mean.cv <- apply(PMSEdata.cv[ , -7], 2, function(x) mean(na.omit(x)))
se.indw <- apply(PMSEtoSA[ , -7], 2, function(x) sd(na.omit(x))/sqrt(length(na.omit(x))))
se.test <- apply(PMSEdata.test[ , -7], 2, function(x) sd(na.omit(x))/sqrt(length(na.omit(x))))
se.cv <- apply(PMSEdata.cv[ , -7], 2, function(x) sd(na.omit(x))/sqrt(length(na.omit(x))))
# Visualization: also including individual weighting methods 
plotdata1 <- data.frame(rbind(mean.indw, mean.test, mean.cv))
plotdata2 <- data.frame(rbind(se.indw, se.test, se.cv))
colnames(plotdata1) <- modelnames[c(1:6, 8)]
plotdata1 
plotdata2
print(min(plotdata1, na.rm = T))
print(max(plotdata1, na.rm = T))
range <- c(0, 50) # for covid 19 forecasting data 
stepsize <- 5

# plot
setEPS()
postscript("PMSE_model_covid.eps", width = 7, height = 3.5)
bars <- barplot(as.matrix(plotdata1*100), ylim = c(range[1], range[2]), yaxt = "n", 
                xlab = "Weighting Model", ylab = "%MSE Reduction to SA", density = c(90, 45, 0), beside = TRUE)
axis(2, seq(range[1], range[2], by = stepsize), labels = paste(seq(range[1], range[2], by = stepsize), "%", sep = ""), las = 1)
segments(bars, as.vector(as.matrix(plotdata1))*100 - as.vector(as.matrix(plotdata2))*300, 
         bars, as.vector(as.matrix(plotdata1))*100 + as.vector(as.matrix(plotdata2))*300, lwd = 2)
arrows(bars, as.vector(as.matrix(plotdata1))*100 - as.vector(as.matrix(plotdata2))*300, 
       bars, as.vector(as.matrix(plotdata1))*100 + as.vector(as.matrix(plotdata2))*300, lwd = 2, angle = 90, code = 3, length = 0.05)
legend("top", c("WA", "Test", "CV"), density = c(90, 45, 0), bg = "white", bty = "n")
dev.off()

#### MSE Reduction: aggregating for each sample size ####
samplesize <- c(16, 32, 64, 128, 256)
sim_runs <- min(data.frame(table(dat$samplesize))$Freq)
PMSEdata.indw <- c()
PMSEdata.test <- c()
PMSEdata.cv   <- c()
crit <- 0.05
for(i in 1:length(samplesize)){
  subdat <- dat[which(dat$samplesize == samplesize[i]), ]
  ## load the data 
  MSEmean.i <- subdat[ , 19:27]
  PMSEtoSA.i <- (MSEmean.i[ , 9] - MSEmean.i)/MSEmean.i[ , 9]
  cverror.i <- subdat[ , 11:18]
  pvalue.i <- subdat[ , 3:10] 
  ## find the MSE for test and cv
  model_k <- 1:8
  PMSEdata.test.i <- matrix(NA, nrow = nrow(PMSEtoSA.i), ncol = length(model_k))
  PMSEdata.cv.i   <- matrix(NA, nrow = nrow(PMSEtoSA.i), ncol = length(model_k))
  for(k in model_k){
    test_select_wa_index <- which(pvalue.i[ , k] <= crit)
    test_select_sa_index <- which(pvalue.i[ , k] > crit)
    cv_select_wa_index   <- which(cverror.i[ , k] <= 0)
    cv_select_sa_index   <- which(cverror.i[ , k] > 0)
    PMSEdata.test.i[test_select_wa_index, k] <- PMSEtoSA.i[test_select_wa_index, k]
    PMSEdata.test.i[test_select_sa_index, k] <- PMSEtoSA.i[test_select_sa_index, 9]
    PMSEdata.cv.i[cv_select_wa_index, k]     <- PMSEtoSA.i[cv_select_wa_index, k]
    PMSEdata.cv.i[cv_select_sa_index, k]     <- PMSEtoSA.i[cv_select_sa_index, 9]
  }
  PMSEtoSA.i <- PMSEtoSA.i[ , -9]
  PMSEtoSA.i <- PMSEtoSA.i[ , -7]
  PMSEdata.test.i <- PMSEdata.test.i[ , -7]
  PMSEdata.cv.i <- PMSEdata.cv.i[ , -7]
  PMSEdata.indw <- cbind(PMSEdata.indw, as.vector(as.matrix(PMSEtoSA.i[1:sim_runs, ])))
  PMSEdata.test <- cbind(PMSEdata.test, as.vector(PMSEdata.test.i[1:sim_runs, ]))
  PMSEdata.cv <- cbind(PMSEdata.cv, as.vector(PMSEdata.cv.i[1:sim_runs, ]))
}
mean.indw <- apply(PMSEdata.indw, 2, function(x) mean(na.omit(x)))
se.indw <- apply(PMSEdata.indw, 2, function(x) sd(na.omit(x))/sqrt(length(na.omit(x))))
mean.test <- apply(PMSEdata.test, 2, function(x) mean(na.omit(x)))
se.test <- apply(PMSEdata.test, 2, function(x) sd(na.omit(x))/sqrt(length(na.omit(x))))
mean.cv <- apply(PMSEdata.cv, 2, function(x) mean(na.omit(x)))
se.cv <- apply(PMSEdata.cv, 2, function(x) sd(na.omit(x))/sqrt(length(na.omit(x))))
# Visualization: including all individual weighting methods 
plotdata1 <- data.frame(rbind(mean.indw, mean.test, mean.cv))
plotdata2 <- data.frame(rbind(se.indw, se.test, se.cv))
colnames(plotdata1) <- as.character(samplesize)
plotdata1
plotdata2
print(min(plotdata1))
print(max(plotdata1))
range <- c(0, 35) # for inflation rate data
stepsize <- 5
# plot
setEPS()
postscript("PMSE_size_covid.eps", width = 4, height = 4)
bars <- barplot(as.matrix(plotdata1*100), ylim = c(range[1], range[2]), yaxt = "n", 
                xlab = "Sample Size", ylab = "%MSE Reduction to SA", density = c(90, 45, 0), beside = TRUE)
axis(2, seq(range[1], range[2], by = stepsize), labels = paste(seq(range[1], range[2], by = stepsize), "%", sep = ""), las = 1)
segments(bars, as.vector(as.matrix(plotdata1))*100 - as.vector(as.matrix(plotdata2))*300, 
         bars, as.vector(as.matrix(plotdata1))*100 + as.vector(as.matrix(plotdata2))*300, lwd = 2)
arrows(bars, as.vector(as.matrix(plotdata1))*100 - as.vector(as.matrix(plotdata2))*300, 
       bars, as.vector(as.matrix(plotdata1))*100 + as.vector(as.matrix(plotdata2))*300, lwd = 2, angle = 90, code = 3, length = 0.05)
legend("topleft", c("WA", "Test", "CV"), density = c(90, 45, 0), bg = "white", bty = "n")
dev.off()

#### FAR/HR table ####
if (!require('data.table')) install.packages('data.table'); library('data.table')
if (!require('mltools')) install.packages('mltools'); library('mltools')
#### FAR/HR table: for each weighting model and each sample size ####
n <- c(16, 32, 64, 128, 256)
weight.name <- c("OW", "TOP3", "RP", "CWM", "SSIN", "SSDE", "REGCOV", "LAS")
crit.outtest <- 0.001
crit <- 0.05
FAR <- array(NA, dim = c(length(n), length(weight.name), 2))
HR <- array(NA, dim = c(length(n), length(weight.name), 2))
CaseCount <- array(NA, dim = c(length(n), length(weight.name), 3))
for(i in 1:length(n)){
  subdat <- dat[which(dat$samplesize == n[i]), ]
  pvalue.i <- subdat[ , 3:10]
  cverror.i <- subdat[ , 11:18]
  MSEout.i <- subdat[ , 19:27]
  MSEout.ttest.i <- subdat[ , 28:35]
  MSEerror.i <- MSEout.i - MSEout.i[ , 9] # error(wa) - error(sa)
  MSEerror.i <- MSEerror.i[ , -9]
  MSEsign.i <- MSEout.ttest.i <= crit.outtest
  for(k in 1:length(weight.name)){
    pvalue.vec <- as.vector(pvalue.i[ ,k])
    cv.vec <- as.vector(cverror.i[ ,k]) # wa wins when cv is negative 
    sign.index <- which(MSEsign.i[ ,k] == TRUE)
    MSE.vec <- rep(0, nrow(MSEerror.i))
    MSE.vec[sign.index] <- as.vector(MSEerror.i[sign.index, k]) # wa wins when MSE is negative 
    tempdata <- cbind(pvalue.vec, cv.vec, MSE.vec)
    remove.index <- which(pvalue.vec == "No need to compare!")
    if(length(remove.index) > 0){
      tempdata <- tempdata[-remove.index, ]
    }
    # find the FAR and HR
    CaseCount[i, k, ] <- c(length(which(tempdata[ , 3] == 0)), # inconclusice
                           length(which(tempdata[ , 3] > 0)), # SA wins significantly
                           length(which(tempdata[ , 3] < 0))) # WA wins significantly
    NullNull_1 <- length(which(tempdata[ , 1] >= crit & tempdata[ , 3] > 0))
    NullNull_2 <- length(which(tempdata[ , 2] >= 0 & tempdata[ , 3] > 0))
    NullReject_1 <- length(which(tempdata[ , 1] < crit & tempdata[ , 3] > 0))
    NullReject_2 <- length(which(tempdata[ , 2] < 0 & tempdata[ , 3] > 0))
    FAR[i, k, ] <- c(NullReject_1/(NullReject_1 + NullNull_1), 
                     NullReject_2/(NullReject_2 + NullNull_2))
    AltRetain_1 <- length(which(tempdata[ , 1] >= crit & tempdata[ , 3] < 0))
    AltRetain_2 <- length(which(tempdata[ , 2] >= 0 & tempdata[ , 3] < 0))
    AltAlt_1 <- length(which(tempdata[ , 1] < crit & tempdata[ , 3] < 0))
    AltAlt_2 <- length(which(tempdata[ , 2] < 0 & tempdata[ , 3] < 0))
    HR[i, k, ] <- c(AltAlt_1/(AltAlt_1 + AltRetain_1), 
                    AltAlt_2/(AltAlt_2 + AltRetain_2))
  }
}
size_index <- 5
CaseCount[ size_index, , ]

t(FAR[ , , 1]*CaseCount[ , , 2])
t(FAR[ , , 2]*CaseCount[ , , 2])
t((1-HR[ , , 1])*CaseCount[ , , 3])
t((1-HR[ , , 2])*CaseCount[ , , 3])

##### FAR/HR table: aggregating for each weighting model ##### 
n <- c(16, 32, 64, 128, 256)
## loading all related data
pvalue <- dat[ , 3:10]
cverror <- dat[ , 11:18]
MSEout <- dat[ , 19:27]
MSEout.ttest <- dat[ , 28:35]
MSEerror <- MSEout - MSEout[ , 9] # error(wa) - error(sa)
MSEerror <- MSEerror[ , -9]
MSEsign <- MSEout.ttest <= 0.001
## find FAR, HR and AUC
weight.name <- c("OW", "TOP3", "RP", "CWM", "SSIN", "SSDE", "REGCOV", "LAS")
crit <- 0.05
FAR <- matrix(NA, nrow = 2, ncol = length(weight.name))
HR <- matrix(NA, nrow = 2, ncol = length(weight.name))
# results tables 
for(k in 1:ncol(pvalue)){
  pvalue.vec <- as.vector(pvalue[ ,k])
  cv.vec <- as.vector(cverror[ ,k]) # wa wins when cv is negative 
  sign.index <- which(MSEsign[ ,k] == TRUE)
  MSE.vec <- rep(0, nrow(MSEerror))
  MSE.vec[sign.index] <- as.vector(MSEerror[sign.index, k]) # wa wins when MSE is negative 
  tempdata <- cbind(pvalue.vec, cv.vec, MSE.vec)
  remove.index <- which(pvalue.vec == "No need to compare!")
  if(length(remove.index) > 0){
    tempdata <- tempdata[-remove.index, ]
  }
  # find the FAR and HR
  NullNull_1 <- length(which(tempdata[ , 1] >= crit & tempdata[ , 3] > 0))
  NullNull_2 <- length(which(tempdata[ , 2] >= 0 & tempdata[ , 3] > 0))
  NullReject_1 <- length(which(tempdata[ , 1] < crit & tempdata[ , 3] > 0))
  NullReject_2 <- length(which(tempdata[ , 2] < 0 & tempdata[ , 3] > 0))
  FAR[ , k] <- c(NullReject_1/(NullReject_1 + NullNull_1), 
                 NullReject_2/(NullReject_2 + NullNull_2))
  AltRetain_1 <- length(which(tempdata[ , 1] >= crit & tempdata[ , 3] < 0))
  AltRetain_2 <- length(which(tempdata[ , 2] >= 0 & tempdata[ , 3] < 0))
  AltAlt_1 <- length(which(tempdata[ , 1] < crit & tempdata[ , 3] < 0))
  AltAlt_2 <- length(which(tempdata[ , 2] < 0 & tempdata[ , 3] < 0))
  HR[ , k] <- c(AltAlt_1/(AltAlt_1 + AltRetain_1), 
                AltAlt_2/(AltAlt_2 + AltRetain_2))
}
# print out the pooled results: 
round(FAR, 3)
round(HR, 3)

###### FAR/HR table: aggregating for each sample size ########
crit.outtest <- 0.001
n.vec <- c(16, 32, 64, 128, 256)
crit <- 0.05
for(j in 1:length(n.vec)){
  n <- n.vec[j]
  subdat <- dat[which(dat$samplesize == n), ]
  pvalue <- subdat[ , 3:10]
  cverror <- subdat[ , 11:18]
  MSEout <- subdat[ , 19:27]
  MSEout.ttest <- subdat[ , 28:35]
  MSEerror <- MSEout - MSEout[ , 9] # error(wa) - error(sa)
  MSEerror <- MSEerror[ , -9]
  MSEsign <- MSEout.ttest <= crit.outtest
  ## find FAR, HR and AUC
  FAR <- c()
  HR <- c()
  # for first 4 weighted averages
  # pvalue.vec <- c(pvalue[,1], pvalue[,2], pvalue[,6], pvalue[,7]) 
  # cv.vec <- c(cverror[,1],cverror[,2],cverror[,6],cverror[,7])
  pvalue.vec <- as.vector(as.matrix(pvalue[ , -7]))
  cv.vec <- as.vector(as.matrix(cverror[ , -7]))
  sign.index <- which(as.vector(as.matrix(MSEsign[ , -7])) == TRUE)
  MSE.vec <- rep(0, length(pvalue.vec))
  # MSEerror.vec <- c(MSEerror[,1],MSEerror[,2],MSEerror[,6],MSEerror[,7])
  MSEerror.vec <- as.vector(as.matrix(MSEerror[ , -7]))
  MSE.vec[sign.index] <- MSEerror.vec[sign.index] # wa wins when MSE is negative 
  tempdata <- cbind(pvalue.vec, cv.vec, MSE.vec)
  remove.index <- which(as.character(pvalue.vec) == "No need to compare!")
  if(length(remove.index) > 0){
    tempdata <- tempdata[-remove.index, ]
  }
  # find the FAR
  NullNull_1 <- length(which(tempdata[ , 1] > crit & tempdata[ , 3] > 0))
  NullNull_2 <- length(which(tempdata[ , 2] > 0 & tempdata[ , 3] > 0))
  NullReject_1 <- length(which(tempdata[ , 1] < crit & tempdata[ , 3] > 0))
  NullReject_2 <- length(which(tempdata[ , 2] < 0 & tempdata[ , 3] > 0))
  FAR <- c(NullReject_1/(NullReject_1 + NullNull_1),
           NullReject_2/(NullReject_2 + NullNull_2))
  AltRetain_1 <- length(which(tempdata[ , 1] > crit & tempdata[ , 3] < 0))
  AltRetain_2 <- length(which(tempdata[ , 2] > 0 & tempdata[ , 3] < 0))
  AltAlt_1 <- length(which(tempdata[ , 1] < crit & tempdata[ , 3] < 0))
  AltAlt_2 <- length(which(tempdata[ , 2] < 0 & tempdata[ , 3] < 0))
  HR <- c(AltAlt_1/(AltAlt_1 + AltRetain_1), 
          AltAlt_2/(AltAlt_2 + AltRetain_2))
  print(round(FAR, 3))
  print(round(HR, 3))
}

########################################################################
################ ROC for simulation  ###################################
########################################################################
#### ROC curves for different weighting methods ####
setwd("insert the path to the related results data. ")
if (!require('data.table')) install.packages('data.table'); library('data.table')
if (!require('mltools')) install.packages('mltools'); library('mltools')
n <- c(16, 32, 64, 128, 256)
## loading all simulation data and standardize the cross validation error 
dat1 <- read.csv("InfStates.csv", header = T, stringsAsFactors = F)
dat1_ow <- read.csv("InfStates_OW.csv", header = T, stringsAsFactors = F)
dat2 <- read.csv("UnempStates.csv", header = T, stringsAsFactors = F)
dat2_ow <- read.csv("UnempStates_OW.csv", header = T, stringsAsFactors = F)
dat3 <- read.csv("CovidStates.csv", header = T, stringsAsFactors = F)
colnames(dat1) <- c("samplesize", "round", "p1", "p2", "p3", "p4", "p5", "p6", "p7", "p8",
                    "cv1", "cv2", "cv3", "cv4", "cv5", "cv6", "cv7", "cv8",
                    "MSEout1", "MSEout2", "MSEout3", "MSEout4", "MSEout5", "MSEout6", "MSEout7", "MSEout8", "MSEout9",
                    "MSEttest1", "MSEttest2", "MSEttest3", "MSEttest4", "MSEttest5", "MSEttest6", "MSEttest7", "MSEttest8")
colnames(dat1_ow) <- c("samplesize", "round", "p1", "p2", "cv1", "cv2", "MSEout1", "MSEout2", "MSEout3", "MSEttest1", "MSEttest2")
colnames(dat2) <- c("samplesize", "round", "p1", "p2", "p3", "p4", "p5", "p6", "p7", "p8",
                    "cv1", "cv2", "cv3", "cv4", "cv5", "cv6", "cv7", "cv8",
                    "MSEout1", "MSEout2", "MSEout3", "MSEout4", "MSEout5", "MSEout6", "MSEout7", "MSEout8", "MSEout9",
                    "MSEttest1", "MSEttest2", "MSEttest3", "MSEttest4", "MSEttest5", "MSEttest6", "MSEttest7", "MSEttest8")
colnames(dat2_ow) <- c("samplesize", "round", "p1", "p2", "cv1", "cv2", "MSEout1", "MSEout2", "MSEout3", "MSEttest1", "MSEttest2")
colnames(dat3) <- c("target", "fips", "samplesize", "p1", "p2", "p3", "p4", "p5", "p6", "p7", "p8",
                    "cv1", "cv2", "cv3", "cv4", "cv5", "cv6", "cv7", "cv8",
                    "MSEout1", "MSEout2", "MSEout3", "MSEout4", "MSEout5", "MSEout6", "MSEout7", "MSEout8", "MSEout9",
                    "MSEttest1", "MSEttest2", "MSEttest3", "MSEttest4", "MSEttest5", "MSEttest6", "MSEttest7", "MSEttest8")
weight.name <- c("OW", "TOP3", "RP", "CWM", "SSIN", "SSDE", "REGCOV", "LAS")
crit.test <- seq(0, 1, by = 0.001)
crit.outtest <- 0.001
#### for optimal weighting method (OW) ####
FAR.test <- matrix(NA, nrow = length(crit.test), ncol = length(weight.name))
HR.test <- matrix(NA, nrow = length(crit.test), ncol = length(weight.name))
FAR.cv <- matrix(NA, nrow = 1000, ncol = length(weight.name))
HR.cv <- matrix(NA, nrow = 1000, ncol = length(weight.name))
AUC <- matrix(NA, nrow = length(weight.name), ncol = 2)

weight_index <- 1
pvalue.vec <- c(dat1_ow$p1, dat2_ow$p1, dat3$p1)
cv.vec <- c(dat1_ow$cv1/sd(dat1_ow$cv1), dat2_ow$cv1/sd(dat2_ow$cv1), dat3$cv1/sd(dat3$cv1))
MSEsign <- c(dat1_ow$MSEttest1<=crit.outtest, dat2_ow$MSEttest1<=crit.outtest, dat3$MSEttest1<=crit.outtest)
sign.index <- which(MSEsign == TRUE)
MSEerror <- c(dat1_ow$MSEout1-dat1_ow$MSEout3, dat2_ow$MSEout1-dat2_ow$MSEout3, dat3$MSEout1-dat3$MSEout9)
MSE.vec <- rep(0, length(MSEerror))
MSE.vec[sign.index] <- MSEerror[sign.index]
tempdat <- cbind(pvalue.vec, cv.vec, MSE.vec)
remove.index <- which(pvalue.vec == "No need to compare!")
if(length(remove.index) > 0){
  tempdat <- tempdat[-remove.index, ]
}
tempdat <- matrix(as.numeric(tempdat), ncol = 3, nrow = nrow(tempdat))
# decide the criterion for cross validation 
sd.cverror <- sd(cv.vec)
crit.cv <- sd.cverror*seq(-3, 3, length.out = 1000)
# find the FAR and HR for our test algorithm 
for(j in 1:length(crit.test)){
  NullNull_1 <- length(which(tempdat[ , 1] >= crit.test[j] & tempdat[ , 3] > 0))
  NullReject_1 <- length(which(tempdat[ , 1] < crit.test[j] & tempdat[ , 3] > 0))
  FAR.test[j, weight_index] <- NullReject_1/(NullReject_1 + NullNull_1)
  AltRetain_1 <- length(which(tempdat[ , 1] >= crit.test[j] & tempdat[ , 3] < 0))
  AltAlt_1 <- length(which(tempdat[ , 1] < crit.test[j] & tempdat[ , 3] < 0))
  HR.test[j, weight_index] <- AltAlt_1/(AltAlt_1 + AltRetain_1)
}
# find the FAR and HR for cross valitempdation 
for(j in 1:length(crit.cv)){
  NullNull_2 <- length(which(tempdat[ , 2] >= -crit.cv[j] & tempdat[ , 3] > 0))
  NullReject_2 <- length(which(tempdat[ , 2] < -crit.cv[j] & tempdat[ , 3] > 0))
  FAR.cv[j, weight_index] <- NullReject_2/(NullReject_2 + NullNull_2)
  AltRetain_2 <- length(which(tempdat[ , 2] >= -crit.cv[j] & tempdat[ , 3] < 0))
  AltAlt_2 <- length(which(tempdat[ , 2] < -crit.cv[j] & tempdat[ , 3] < 0))
  HR.cv[j, weight_index] <- AltAlt_2/(AltAlt_2 + AltRetain_2)
}
# find the AUC
MSE.vec <- rep(NA, length(MSEerror))
MSE.vec[sign.index] <- MSEerror[sign.index] # wa wins when MSE is negative 
MSE.vec[which(MSE.vec < 0)] <- 0
MSE.vec[which(MSE.vec > 0)] <- 1
# prediction from cross validation 
pred.cv <- sapply(1:length(cv.vec), function(x) length(which(-crit.cv < cv.vec[x]))/length(crit.cv))
pred.test <- pvalue.vec
temp <- cbind(pred.test, pred.cv, MSE.vec)
remove.index <- which(pvalue.vec == "No need to compare!")
if(length(remove.index) > 0){
  temp <- temp[-remove.index, ]
}
temp <- na.omit(temp)
AUC[weight_index, 1] <- auc_roc(as.numeric(temp[ , 1]), as.numeric(temp[ , 3])) 
AUC[weight_index, 2] <- auc_roc(as.numeric(temp[ , 2]), as.numeric(temp[ , 3]))
#### for other weighting methods ####
weight_index <- c(2:6, 8)
for(k in weight_index){
  pvalue.vec <- c(dat1[ , k+2], dat2[ , k+2], dat3[ , k+3])
  cv.vec <- c(dat1[ , k+10]/sd(dat1[ , k+10]), dat2[ , k+10]/sd(dat2[ , k+10]), dat3[ , k+11]/sd(dat3[ , k+11]))
  MSEsign <- c(dat1[ , k+27]<=crit.outtest, dat2[ , k+27]<=crit.outtest, dat3[ , k+28]<=crit.outtest)
  sign.index <- which(MSEsign == TRUE)
  MSEerror <- c(dat1[ , k+18]-dat1[ , 27], dat2[ , k+18]-dat2[ , 27], dat3[ , k+19]-dat3[ , 28])
  MSE.vec <- rep(0, length(MSEerror))
  MSE.vec[sign.index] <- MSEerror[sign.index]
  tempdat <- cbind(pvalue.vec, cv.vec, MSE.vec)
  remove.index <- which(pvalue.vec == "No need to compare!")
  if(length(remove.index) > 0){
    tempdat <- tempdat[-remove.index, ]
  }
  tempdat <- matrix(as.numeric(tempdat), ncol = 3, nrow = nrow(tempdat))
  # decide the criterion for cross validation 
  sd.cverror <- sd(cv.vec)
  crit.cv <- sd.cverror*seq(-3, 3, length.out = 1000)
  # find the FAR and HR for our test algorithm 
  for(j in 1:length(crit.test)){
    NullNull_1 <- length(which(tempdat[ , 1] >= crit.test[j] & tempdat[ , 3] > 0))
    NullReject_1 <- length(which(tempdat[ , 1] < crit.test[j] & tempdat[ , 3] > 0))
    FAR.test[j, k] <- NullReject_1/(NullReject_1 + NullNull_1)
    AltRetain_1 <- length(which(tempdat[ , 1] >= crit.test[j] & tempdat[ , 3] < 0))
    AltAlt_1 <- length(which(tempdat[ , 1] < crit.test[j] & tempdat[ , 3] < 0))
    HR.test[j, k] <- AltAlt_1/(AltAlt_1 + AltRetain_1)
  }
  # find the FAR and HR for cross valitempdation 
  for(j in 1:length(crit.cv)){
    NullNull_2 <- length(which(tempdat[ , 2] >= -crit.cv[j] & tempdat[ , 3] > 0))
    NullReject_2 <- length(which(tempdat[ , 2] < -crit.cv[j] & tempdat[ , 3] > 0))
    FAR.cv[j, k] <- NullReject_2/(NullReject_2 + NullNull_2)
    AltRetain_2 <- length(which(tempdat[ , 2] >= -crit.cv[j] & tempdat[ , 3] < 0))
    AltAlt_2 <- length(which(tempdat[ , 2] < -crit.cv[j] & tempdat[ , 3] < 0))
    HR.cv[j, k] <- AltAlt_2/(AltAlt_2 + AltRetain_2)
  }
  # find the AUC
  MSE.vec <- rep(NA, length(MSEerror))
  MSE.vec[sign.index] <- MSEerror[sign.index] # wa wins when MSE is negative 
  MSE.vec[which(MSE.vec < 0)] <- 0
  MSE.vec[which(MSE.vec > 0)] <- 1
  # prediction from cross validation 
  pred.cv <- sapply(1:length(cv.vec), function(x) length(which(-crit.cv < cv.vec[x]))/length(crit.cv))
  pred.test <- pvalue.vec
  temp <- cbind(pred.test, pred.cv, MSE.vec)
  remove.index <- which(pvalue.vec == "No need to compare!")
  if(length(remove.index) > 0){
    temp <- temp[-remove.index, ]
  }
  temp <- na.omit(temp)
  AUC[k, 1] <- auc_roc(as.numeric(temp[ , 1]), as.numeric(temp[ , 3])) 
  AUC[k, 2] <- auc_roc(as.numeric(temp[ , 2]), as.numeric(temp[ , 3]))
  #### for other weighting methods ####
}

#### plotting the ROC curve ####
weight.name <- c("OW", "TOP3", "RP", "CWM", "SSIN", "SSDE", "REGCOV", "LAS")
for(k in c(1:6, 8)){
  setEPS()
  postscript(paste("ROC_", weight.name[k], ".eps", sep = ""), width = 5, height = 5)
  plot(FAR.test[ , k], HR.test[ , k], ylim = c(0, 1), xlim = c(0, 1),
       type = "l", lwd = 3, lty = 1,
       xlab = "False Alarm Rate", ylab = "Hit Rate", main = "")
  lines(c(0, FAR.test[ , k], 1), c(0, HR.test[ , k], 1), type = "l", lwd = 3)
  segments(0, 0, 1, 1, lty = 3)
  dat.cv <- cbind(FAR.cv[ , k], HR.cv[ , k])
  dat.cv <- rbind(dat.cv, c(0, 0), c(1, 1))
  dat.cv <- dat.cv[order(dat.cv[ , 1]), ]
  lines(dat.cv[ , 1], dat.cv[ , 2], lwd = 3, lty = 4, type = "l")
  legend("bottomright", legend = c("Hypothesis test", "Cross validation"), lwd = 2, lty = c(1, 4))
  text(0.7, 0.6, "AUC of Hypothesis Test:")
  text(0.7, 0.5, round(AUC[k, 1], 3), cex = 1.5)
  text(0.7, 0.4, "AUC of Cross Validation:")
  text(0.7, 0.3, round(AUC[k, 2], 3), cex = 1.5)
  #text(0.7, 0.3, "0.730", cex = 1.5) # for k = 8
  dev.off()
}
##### ROC curves for different sample sizes #####
FAR.test <- matrix(NA, nrow = length(crit.test), ncol = length(n))
HR.test <- matrix(NA, nrow = length(crit.test), ncol = length(n))
FAR.cv <- matrix(NA, nrow = 1000, ncol = length(n))
HR.cv <- matrix(NA, nrow = 1000, ncol = length(n))
AUC <- matrix(NA, nrow = length(n), ncol = 2)

for(n_index in 1:length(n)){
  subdat1 <- dat1[which(dat1$samplesize == n[n_index]), ]
  subdat1_ow <- dat1_ow[which(dat1_ow$samplesize == n[n_index]), ]
  subdat2 <- dat2[which(dat2$samplesize == n[n_index]), ]
  subdat2_ow <- dat2_ow[which(dat2_ow$samplesize == n[n_index]), ]
  subdat3 <- dat3[which(dat3$samplesize == n[n_index]), ]
  # organize the data 
  pvalue.vec <- c(subdat1_ow$p1, as.vector(as.matrix(subdat1[ , c(4:8, 10)])), 
                  subdat2_ow$p1, as.vector(as.matrix(subdat2[ , c(4:8, 10)])), 
                  as.vector(as.matrix(subdat3[ , c(4:9, 11)])))
  cv.vec <- c(subdat1_ow$cv1/sd(subdat1_ow$cv1), as.vector(as.matrix(subdat1[ , c(12:16, 18)]))/sd(as.vector(as.matrix(subdat1[ , c(12:16, 18)]))), 
              subdat2_ow$cv1/sd(subdat2_ow$cv1), as.vector(as.matrix(subdat2[ , c(12:16, 18)]))/sd(as.vector(as.matrix(subdat2[ , c(12:16, 18)]))),
              as.vector(as.matrix(subdat3[ , c(12:17, 19)]))/sd(as.vector(as.matrix(subdat3[ , c(12:17, 19)]))))
  MSEsign <- c(subdat1_ow$MSEttest1<=crit.outtest, as.vector(as.matrix(subdat1[ , c(29:33, 35)]))<=crit.outtest, 
               subdat2_ow$MSEttest1<=crit.outtest, as.vector(as.matrix(subdat2[ , c(29:33, 35)]))<=crit.outtest, 
               as.vector(as.matrix(subdat3[ , c(29:34, 36)]))<=crit.outtest)
  sign.index <- which(MSEsign == TRUE)
  MSEerror <- c(subdat1_ow$MSEout1-subdat1_ow$MSEout3, as.vector(as.matrix(subdat1[ , c(20:24, 26)]-subdat1[ , 27])),
                subdat2_ow$MSEout1-subdat2_ow$MSEout3, as.vector(as.matrix(subdat2[ , c(20:24, 26)]-subdat2[ , 27])),
                as.vector(as.matrix(subdat3[ , c(20:25, 27)]-subdat3[ , 28])))
  MSE.vec <- rep(0, length(MSEerror))
  MSE.vec[sign.index] <- MSEerror[sign.index]
  tempdat <- cbind(pvalue.vec, cv.vec, MSE.vec)
  remove.index <- which(pvalue.vec == "No need to compare!")
  if(length(remove.index) > 0){
    tempdat <- tempdat[-remove.index, ]
  }
  tempdat <- matrix(as.numeric(tempdat), ncol = 3, nrow = nrow(tempdat))
  # decide the criterion for cross validation 
  sd.cverror <- sd(cv.vec)
  crit.cv <- sd.cverror*seq(-3, 3, length.out = 1000)
  # find the FAR and HR for our test algorithm 
  for(j in 1:length(crit.test)){
    NullNull_1 <- length(which(tempdat[ , 1] >= crit.test[j] & tempdat[ , 3] > 0))
    NullReject_1 <- length(which(tempdat[ , 1] < crit.test[j] & tempdat[ , 3] > 0))
    FAR.test[j, n_index] <- NullReject_1/(NullReject_1 + NullNull_1)
    AltRetain_1 <- length(which(tempdat[ , 1] >= crit.test[j] & tempdat[ , 3] < 0))
    AltAlt_1 <- length(which(tempdat[ , 1] < crit.test[j] & tempdat[ , 3] < 0))
    HR.test[j, n_index] <- AltAlt_1/(AltAlt_1 + AltRetain_1)
  }
  # find the FAR and HR for cross valitempdation 
  for(j in 1:length(crit.cv)){
    NullNull_2 <- length(which(tempdat[ , 2] >= -crit.cv[j] & tempdat[ , 3] > 0))
    NullReject_2 <- length(which(tempdat[ , 2] < -crit.cv[j] & tempdat[ , 3] > 0))
    FAR.cv[j, n_index] <- NullReject_2/(NullReject_2 + NullNull_2)
    AltRetain_2 <- length(which(tempdat[ , 2] >= -crit.cv[j] & tempdat[ , 3] < 0))
    AltAlt_2 <- length(which(tempdat[ , 2] < -crit.cv[j] & tempdat[ , 3] < 0))
    HR.cv[j, n_index] <- AltAlt_2/(AltAlt_2 + AltRetain_2)
  }
  # find the AUC
  MSE.vec <- rep(NA, length(MSEerror))
  MSE.vec[sign.index] <- MSEerror[sign.index] # wa wins when MSE is negative 
  MSE.vec[which(MSE.vec < 0)] <- 0
  MSE.vec[which(MSE.vec > 0)] <- 1
  # prediction from cross validation 
  pred.cv <- sapply(1:length(cv.vec), function(x) length(which(-crit.cv < cv.vec[x]))/length(crit.cv))
  pred.test <- pvalue.vec
  temp <- cbind(pred.test, pred.cv, MSE.vec)
  remove.index <- which(pvalue.vec == "No need to compare!")
  if(length(remove.index) > 0){
    temp <- temp[-remove.index, ]
  }
  temp <- na.omit(temp)
  AUC[n_index, 1] <- auc_roc(as.numeric(temp[ , 1]), as.numeric(temp[ , 3])) 
  AUC[n_index, 2] <- auc_roc(as.numeric(temp[ , 2]), as.numeric(temp[ , 3]))
}
#### visualization: plotting ROC curves by aggregating all weighting methods ####
for(n_index in 1:length(n)){
  setEPS()
  postscript(paste("ROC_n", n[n_index], ".eps", sep = ""), width = 5, height = 5)
  plotdata.test <- cbind(FAR.test[ , n_index], HR.test[ , n_index])
  plotdata.test <- plotdata.test[order(plotdata.test[ , 1]), ]
  plotdata.cv <- cbind(FAR.cv[ , n_index], HR.cv[ , n_index])
  plotdata.cv <- plotdata.cv[order(plotdata.cv[ , 1]), ]
  plotdata.test <- rbind(c(0, 0), plotdata.test, c(1, 1))
  plotdata.cv <- rbind(c(0, 0), plotdata.cv, c(1, 1))
  
  plot(plotdata.test[ , 1], plotdata.test[ , 2], ylim = c(0, 1), xlim = c(0, 1),
       type = "l", lwd = 3, lty = 1,
       xlab = "False Alarm Rate", ylab = "Hit Rate")
  lines(plotdata.cv[ , 1], plotdata.cv[ , 2], lwd = 3, lty = 4)
  segments(0, 0, 1, 1, lty = 3)
  legend("bottomright", legend = c("Hypothesis test", "Cross validation"),lwd = 2, lty = c(1, 4))
  text(0.75, 0.5, "AUC of Hypothesis Test:")
  text(0.75, 0.4, round(AUC[n_index, 1], 3), cex = 1.5)
  text(0.75, 0.3, "AUC of Cross Validation:")
  text(0.75, 0.2, round(AUC[n_index, 2], 3), cex = 1.5)
  dev.off()
}
#
#par(pty = "s")


##### for covid-19 forecasting data #####
dat <- read.csv("CovidStates.csv", header = T, stringsAsFactors = F)
dat <- dat[ , -2]
colnames(dat) <- c("target", "samplesize", "p1", "p2", "p3", "p4", "p5", "p6", "p7", "p8",
                   "cv1", "cv2", "cv3", "cv4", "cv5", "cv6", "cv7", "cv8",
                   "MSEout1", "MSEout2", "MSEout3", "MSEout4", "MSEout5", "MSEout6", "MSEout7", "MSEout8", "MSEout9",
                   "MSEttest1", "MSEttest2", "MSEttest3", "MSEttest4", "MSEttest5", "MSEttest6", "MSEttest7", "MSEttest8")

crit.outtest <- 0.001
n <- 32 # change the sample sizes
#### plot the ROC curve according to the sample size ####
subdat <- dat[which(dat$samplesize == n), ]
## reading all related data
pvalue <- subdat[ , 3:10]
cverror <- subdat[ , 11:18]
MSEmean <- subdat[ , 19:27]
MSEerror <- MSEmean - MSEmean[ , 9]
MSEerror <- MSEerror[ , -9]
MSEout.ttest <- subdat[ , 28:35]
MSEsign <- MSEout.ttest <= crit.outtest

pvalue <- pvalue[ , weights.num]
cverror <- cverror[, weights.num]
MSEerror <- MSEerror[ , weights.num]
MSEsign <- MSEsign[ , weights.num]
MSEout.ttest <- MSEout.ttest[ , weights.num]

crit.test <- seq(0, 1, by = 0.001)
sd.cverror <- sd(cverror)
crit.cv <- sd.cverror*seq(-1, 1, length.out = 1000)
FAR.test <- c()
HR.test <- c()
FAR.cv <- c()
HR.cv <- c()
AUC <- c()
pvalue.vec <- as.vector(pvalue)
remove.index <- which(as.character(pvalue.vec) == "No need to compare!")
cv.vec <- as.vector(cverror) # wa wins when cv is negative 
MSEsign.vec <- as.vector(MSEsign)
sign.index <- which(MSEsign.vec == TRUE)
MSE.vec <- rep(0, length(pvalue.vec))
MSE.vec[sign.index] <- as.vector(MSEerror)[sign.index] # wa wins when MSE is negative 
tempdata <- cbind(pvalue.vec, cv.vec, MSE.vec)
if(length(remove.index) > 0){
  tempdata <- tempdata[-remove.index, ]
}
# find the FAR and HR
tempdata <- matrix(as.numeric(tempdata), ncol = 3, nrow = nrow(tempdata))
for(j in 1:length(crit.test)){
  NullNull_1 <- length(which(tempdata[ , 1] >= crit.test[j] & tempdata[ , 3] > 0))
  NullReject_1 <- length(which(tempdata[ , 1] < crit.test[j] & tempdata[ , 3] > 0))
  FAR.test[j] <- NullReject_1/(NullReject_1 + NullNull_1)
  AltRetain_1 <- length(which(tempdata[ , 1] >= crit.test[j] & tempdata[ , 3] < 0))
  AltAlt_1 <- length(which(tempdata[ , 1] < crit.test[j] & tempdata[ , 3] < 0))
  HR.test[j] <- AltAlt_1/(AltAlt_1 + AltRetain_1)
}
for(j in 1:length(crit.cv)){
  NullNull_2 <- length(which(tempdata[ , 2] >= -crit.cv[j] & tempdata[ , 3] > 0))
  NullReject_2 <- length(which(tempdata[ , 2] < -crit.cv[j] & tempdata[ , 3] > 0))
  FAR.cv[j] <- NullReject_2/(NullReject_2 + NullNull_2)
  AltRetain_2 <- length(which(tempdata[ , 2] >= -crit.cv[j] & tempdata[ , 3] < 0))
  AltAlt_2 <- length(which(tempdata[ , 2] < -crit.cv[j] & tempdata[ , 3] < 0))
  HR.cv[j] <- AltAlt_2/(AltAlt_2 + AltRetain_2)
}
# find the AUC
MSE.vec.auc <- rep(NA, length(MSE.vec))
MSE.vec.auc[which(MSE.vec < 0)] <- 0
MSE.vec.auc[which(MSE.vec > 0)] <- 1
pred.cv <- sapply(1:length(cv.vec), function(x) length(which(-crit.cv < cv.vec[x]))/length(crit.cv))
pred.test <- pvalue.vec
tempdata.auc <- cbind(pred.test, pred.cv, MSE.vec.auc)
# remove all inconclusive cases
remove.index <- which(as.character(pvalue.vec) == "No need to compare!")
if(length(remove.index) > 0){
  tempdata.auc <- tempdata.auc[-remove.index, ]
}
tempdata.auc <- na.omit(tempdata.auc)
AUC[1] <- auc_roc(as.numeric(tempdata.auc[ , 1]), as.numeric(tempdata.auc[ , 3])) # for our hypothesis test algorithm
AUC[2] <- auc_roc(as.numeric(tempdata.auc[ , 2]), as.numeric(tempdata.auc[ , 3])) # for cross validation

## visualization: plotting ROC curves for optimally weighted average
plotdata.test <- cbind(FAR.test, HR.test)
plotdata.test <- plotdata.test[order(plotdata.test[ , 1]), ]
plotdata.cv <- cbind(FAR.cv, HR.cv)
plotdata.cv <- plotdata.cv[order(plotdata.cv[ , 1]), ]
plotdata.test <- rbind(c(0, 0), plotdata.test, c(1, 1))
plotdata.cv <- rbind(c(0, 0), plotdata.cv, c(1, 1))
#
setEPS()
postscript(paste("ROC_OW_n", n, "_covid.eps", sep = ""), width = 5, height = 5)
#par(pty = "s")
plot(plotdata.test[ , 1], plotdata.test[ , 2], ylim = c(0, 1), xlim = c(0, 1),
     type = "l", lwd = 3, lty = 1,
     xlab = "False Alarm Rate", ylab = "Hit Rate")
lines(plotdata.cv[ , 1], plotdata.cv[ , 2], lwd = 3, lty = 4)
segments(0, 0, 1, 1, lty = 3)
legend("bottomright", legend = c("Hypothesis test", "Cross validation"),lwd = 2, lty = c(1, 4))
text(0.75, 0.5, "AUC of Hypothesis Test:")
text(0.75, 0.4, round(AUC[1], 3), cex = 1.5)
text(0.75, 0.3, "AUC of Cross Validation:")
text(0.75, 0.2, round(AUC[2], 3), cex = 1.5)
dev.off()

##################################################################
############ Empirical Aalysis: SPF data #########################
##################################################################
setwd("insert the path to the related results data. ")
dat <- read.csv("res_inf.csv", header = F, stringsAsFactors = F)
MSE_SA <- read.csv("MSEout_inf.csv", header = F, stringsAsFactors = F)
#dat <- read.csv("res_unemp.csv", header = F, stringsAsFactors = F)
#MSE_SA <- read.csv("MSEout_unemp.csv", header = F, stringsAsFactors = F)

colnames(dat) <- c("datasplit", "samplesize", "round", "mse_ow0", "mse_ow1", "mse_cwm", "p_ow0", "p_ow1", "p_cwm")
mse_sa_split <- MSE_SA[ , c(1, 5)]
colnames(mse_sa_split) <- c("datasplit", "mse_sa")
##### seperate the data split to count the rejection proportion #####
data_split <- sort(unique(dat$datasplit))
samplesize <- sort(unique(dat$samplesize))
crit <- 0.05
rates <- array(NA, dim = c(length(data_split), length(samplesize), 3, 4)) # 3 weighting methods and 4 rates 
for(i in 1:length(data_split)){
  mse_sa <- mse_sa_split$mse_sa[which(mse_sa_split$datasplit == data_split[i])]
  for(j in 1:length(samplesize)){
    subdat <- dat[which(dat$datasplit == data_split[i] & dat$samplesize == samplesize[j]), ]
    rates[i, j, 1, 1] <- length(which(subdat$p_ow0 <= crit & subdat$mse_ow0 <= mse_sa))/nrow(subdat) # correct reject
    rates[i, j, 1, 2] <- length(which(subdat$p_ow0 <= crit & subdat$mse_ow0 >  mse_sa))/nrow(subdat) # incorrect reject
    rates[i, j, 1, 3] <- length(which(subdat$p_ow0 >  crit & subdat$mse_ow0 <  mse_sa))/nrow(subdat) # incorrect retain
    rates[i, j, 1, 4] <- length(which(subdat$p_ow0 >  crit & subdat$mse_ow0 >= mse_sa))/nrow(subdat) # correct retain
    rates[i, j, 2, 1] <- length(which(subdat$p_ow1 <= crit & subdat$mse_ow1 <= mse_sa))/nrow(subdat) # correct reject
    rates[i, j, 2, 2] <- length(which(subdat$p_ow1 <= crit & subdat$mse_ow1 >  mse_sa))/nrow(subdat) # incorrect reject
    rates[i, j, 2, 3] <- length(which(subdat$p_ow1 >  crit & subdat$mse_ow1 <  mse_sa))/nrow(subdat) # incorrect retain
    rates[i, j, 2, 4] <- length(which(subdat$p_ow1 >  crit & subdat$mse_ow1 >= mse_sa))/nrow(subdat) # correct retain
    rates[i, j, 3, 1] <- length(which(subdat$p_cwm <= crit & subdat$mse_cwm <= mse_sa))/nrow(subdat) # correct reject
    rates[i, j, 3, 2] <- length(which(subdat$p_cwm <= crit & subdat$mse_cwm >  mse_sa))/nrow(subdat) # incorrect reject
    rates[i, j, 3, 3] <- length(which(subdat$p_cwm >  crit & subdat$mse_cwm <  mse_sa))/nrow(subdat) # incorrect retain
    rates[i, j, 3, 4] <- length(which(subdat$p_cwm >  crit & subdat$mse_cwm >= mse_sa))/nrow(subdat) # correct retain
  }
}
ave.rates <- apply(rates, c(2, 3, 4), mean)
values.OW0 <- t(ave.rates[ , 1, 1:2])
values.OW1 <- t(ave.rates[ , 2, 1:2])
values.CWM <- t(ave.rates[ , 3, 1:2])
##### combine all data splits to count the rejection proportion #####
samplesize <- sort(unique(dat$samplesize))
crit <- 0.05
rates <- array(NA, dim = c(length(samplesize), 3, 4)) # 3 weighting methods and 4 rates 
mse_sa_vec <- unlist(sapply(dat$datasplit, function(x) mse_sa_split$mse_sa[which(mse_sa_split$datasplit == x)]))
dat$mse_sa <- mse_sa_vec
for(j in 1:length(samplesize)){
  subdat <- dat[which(dat$samplesize == samplesize[j]), ]
  rates[j, 1, 1] <- length(which(subdat$p_ow0 <= crit & subdat$mse_ow0 <= subdat$mse_sa))/nrow(subdat) # correct reject
  rates[j, 1, 2] <- length(which(subdat$p_ow0 <= crit & subdat$mse_ow0 >  subdat$mse_sa))/nrow(subdat) # incorrect reject
  rates[j, 1, 3] <- length(which(subdat$p_ow0 >  crit & subdat$mse_ow0 <  subdat$mse_sa))/nrow(subdat) # incorrect retain
  rates[j, 1, 4] <- length(which(subdat$p_ow0 >  crit & subdat$mse_ow0 >= subdat$mse_sa))/nrow(subdat) # correct retain
  rates[j, 2, 1] <- length(which(subdat$p_ow1 <= crit & subdat$mse_ow1 <= subdat$mse_sa))/nrow(subdat) # correct reject
  rates[j, 2, 2] <- length(which(subdat$p_ow1 <= crit & subdat$mse_ow1 >  subdat$mse_sa))/nrow(subdat) # incorrect reject
  rates[j, 2, 3] <- length(which(subdat$p_ow1 >  crit & subdat$mse_ow1 <  subdat$mse_sa))/nrow(subdat) # incorrect retain
  rates[j, 2, 4] <- length(which(subdat$p_ow1 >  crit & subdat$mse_ow1 >= subdat$mse_sa))/nrow(subdat) # correct retain
  rates[j, 3, 1] <- length(which(subdat$p_cwm <= crit & subdat$mse_cwm <= subdat$mse_sa))/nrow(subdat) # correct reject
  rates[j, 3, 2] <- length(which(subdat$p_cwm <= crit & subdat$mse_cwm >  subdat$mse_sa))/nrow(subdat) # incorrect reject
  rates[j, 3, 3] <- length(which(subdat$p_cwm >  crit & subdat$mse_cwm <  subdat$mse_sa))/nrow(subdat) # incorrect retain
  rates[j, 3, 4] <- length(which(subdat$p_cwm >  crit & subdat$mse_cwm >= subdat$mse_sa))/nrow(subdat) # correct retain
}
values.OW0 <- t(rates[ , 1, 1:2])
values.OW1 <- t(rates[ , 2, 1:2])
values.CWM <- t(rates[ , 3, 1:2])
##### visualization #####
xaxisname <- as.character(samplesize)

setEPS()
postscript("barchart_inf.eps", width = 10, height = 4.5)
par(mfrow = c(1, 2))
# plotdata0 <- values.OW0
# barplot(plotdata0, names.arg = xaxisname,
#         xlab = "Sample Size", ylab = "Proportion", main = "OW vs. SA", 
#         density = c(15, 0), ylim = c(0, 1))
# items <- c("Total rejections", "Correct rejections (OW0)")
# legend("topleft", items, density = c(0, 15), bg = "white")
plotdata1 <- values.OW1
barplot(plotdata1, names.arg = xaxisname,
         xlab = "Sample Size", ylab = "Proportion", main = "OW vs. SA", 
         density = c(15, 0), ylim = c(0, 1))
items <- c("Total rejections", "Correct rejections (OW)")
legend("bottomright", items, density = c(0, 15), bg = "white")
plotdata2 <- values.CWM
barplot(plotdata2, names.arg = xaxisname,
        xlab = "Sample Size", ylab = "Proportion", main = "CWM vs. SA", 
        density = c(75, 0), ylim = c(0, 1))
items <- c("Total rejections", "Correct rejections (CWM)")
legend("bottomright", items, density = c(0, 75), bg = "white")
dev.off()

##################################################################
############ Empirical Aalysis: Covid-19 cases ###################
##################################################################
setwd("insert the path to the related results data. ")
dat <- read.csv("res_covid_cases.csv", header = T, stringsAsFactors = F)
colnames(dat) <- c("target", "fips", "samplesize", "mse_ow0", "mse_ow1", "mse_cwm", "mse_sa", "p_ow0", "p_ow1", "p_cwm")
table(dat$samplesize)
target <- sort(unique(dat$target))
samplesize <- sort(unique(dat$samplesize))
#### seperate each target ####
crit <- 0.05
rates <- array(NA, dim = c(length(target), length(samplesize), 3, 4)) # 3 weighting mehtods, 4 rates
for(i in 1:length(target)){
  subdat <- dat[which(dat$target == target[i]), ]
  for(j in 1:length(samplesize)){
    subsubdat <- subdat[which(subdat$samplesize == samplesize[j]), ]
    rates[i, j, 1, 1] <- length(which(subsubdat$p_ow0 <= crit & subsubdat$mse_ow0 <= subsubdat$mse_sa))/nrow(subsubdat) # correct reject
    rates[i, j, 1, 2] <- length(which(subsubdat$p_ow0 <= crit & subsubdat$mse_ow0 >  subsubdat$mse_sa))/nrow(subsubdat) # incorrect reject
    rates[i, j, 1, 3] <- length(which(subsubdat$p_ow0 >  crit & subsubdat$mse_ow0 <  subsubdat$mse_sa))/nrow(subsubdat) # incorrect retain
    rates[i, j, 1, 4] <- length(which(subsubdat$p_ow0 >  crit & subsubdat$mse_ow0 >= subsubdat$mse_sa))/nrow(subsubdat) # correct retain
    rates[i, j, 2, 1] <- length(which(subsubdat$p_ow1 <= crit & subsubdat$mse_ow1 <= subsubdat$mse_sa))/nrow(subsubdat) # correct reject
    rates[i, j, 2, 2] <- length(which(subsubdat$p_ow1 <= crit & subsubdat$mse_ow1 >  subsubdat$mse_sa))/nrow(subsubdat) # incorrect reject
    rates[i, j, 2, 3] <- length(which(subsubdat$p_ow1 >  crit & subsubdat$mse_ow1 <  subsubdat$mse_sa))/nrow(subsubdat) # incorrect retain
    rates[i, j, 2, 4] <- length(which(subsubdat$p_ow1 >  crit & subsubdat$mse_ow1 >= subsubdat$mse_sa))/nrow(subsubdat) # correct retain
    rates[i, j, 3, 1] <- length(which(subsubdat$p_cwm <= crit & subsubdat$mse_cwm <= subsubdat$mse_sa))/nrow(subsubdat) # correct reject
    rates[i, j, 3, 2] <- length(which(subsubdat$p_cwm <= crit & subsubdat$mse_cwm >  subsubdat$mse_sa))/nrow(subsubdat) # incorrect reject
    rates[i, j, 3, 3] <- length(which(subsubdat$p_cwm >  crit & subsubdat$mse_cwm <  subsubdat$mse_sa))/nrow(subsubdat) # incorrect retain
    rates[i, j, 3, 4] <- length(which(subsubdat$p_cwm >  crit & subsubdat$mse_cwm >= subsubdat$mse_sa))/nrow(subsubdat) # correct retain
  }
}
# if aggregating all targets 
ave.rates <- apply(rates, c(2, 3, 4), mean)
values.OW0 <- t(ave.rates[ , 1, 1:2])
values.OW1 <- t(ave.rates[ , 2, 1:2])
values.CWM <- t(ave.rates[ , 3, 1:2])
# if seperate all targets
target_index <- 1
values.OW0 <- t(rates[target_index, , 1, 1:2])
values.OW1 <- t(rates[target_index, , 2, 1:2])
values.CWM <- t(rates[target_index, , 3, 1:2])
#### Alternative summary method: combining all four targets ####
rates <- array(NA, dim = c(length(samplesize), 3, 4)) # 3 weighting mehtods, 4 rates
crit <- 0.05
for(j in 1:length(samplesize)){
  subdat <- dat[which(dat$samplesize == samplesize[j]), ]
  rates[j, 1, 1] <- length(which(subdat$p_ow0 <= crit & subdat$mse_ow0 <= subdat$mse_sa))/nrow(subdat) # correct reject
  rates[j, 1, 2] <- length(which(subdat$p_ow0 <= crit & subdat$mse_ow0 >  subdat$mse_sa))/nrow(subdat) # incorrect reject
  rates[j, 1, 3] <- length(which(subdat$p_ow0 >  crit & subdat$mse_ow0 <  subdat$mse_sa))/nrow(subdat) # incorrect retain
  rates[j, 1, 4] <- length(which(subdat$p_ow0 >  crit & subdat$mse_ow0 >= subdat$mse_sa))/nrow(subdat) # correct retain
  rates[j, 2, 1] <- length(which(subdat$p_ow1 <= crit & subdat$mse_ow1 <= subdat$mse_sa))/nrow(subdat) # correct reject
  rates[j, 2, 2] <- length(which(subdat$p_ow1 <= crit & subdat$mse_ow1 >  subdat$mse_sa))/nrow(subdat) # incorrect reject
  rates[j, 2, 3] <- length(which(subdat$p_ow1 >  crit & subdat$mse_ow1 <  subdat$mse_sa))/nrow(subdat) # incorrect retain
  rates[j, 2, 4] <- length(which(subdat$p_ow1 >  crit & subdat$mse_ow1 >= subdat$mse_sa))/nrow(subdat) # correct retain
  rates[j, 3, 1] <- length(which(subdat$p_cwm <= crit & subdat$mse_cwm <= subdat$mse_sa))/nrow(subdat) # correct reject
  rates[j, 3, 2] <- length(which(subdat$p_cwm <= crit & subdat$mse_cwm >  subdat$mse_sa))/nrow(subdat) # incorrect reject
  rates[j, 3, 3] <- length(which(subdat$p_cwm >  crit & subdat$mse_cwm <  subdat$mse_sa))/nrow(subdat) # incorrect retain
  rates[j, 3, 4] <- length(which(subdat$p_cwm >  crit & subdat$mse_cwm >= subdat$mse_sa))/nrow(subdat) # correct retain
}
values.OW0 <- t(rates[ , 1, 1:2])
values.OW1 <- t(rates[ , 2, 1:2])
values.CWM <- t(rates[ , 3, 1:2])
#### visualization ####
plotdata.OW0 <- colSums(values.OW0)
plotdata.OW1 <- colSums(values.OW1)
plotdata.CWM <- colSums(values.CWM)
xaxisname <- as.character(samplesize)

setEPS()
postscript("barchart_covid_case.eps", width = 6, height = 4.5)
par(mfrow = c(1, 2))
barplot(plotdata.OW1, names.arg = xaxisname, col = "white", 
        xlab = "Sample Size", ylab = "Proportion", main = "OW vs. SA", ylim = c(0, 1))
items <- c("Total rejections")
legend("topleft", items, density = c(0), bg = "white")
barplot(plotdata.CWM, names.arg = xaxisname, col = "white", 
        xlab = "Sample Size", ylab = "Proportion", main = "CWM vs. SA", ylim = c(0, 1))
items <- c("Total rejections")
legend("topleft", items, density = c(0), bg = "white")
dev.off()
##################################################################
############ Empirical Aalysis: Covid-19 deaths ##################
##################################################################
setwd("insert the path to the related results data. ")
dat <- read.csv("res_covid_deaths.csv", header = T, stringsAsFactors = F)
colnames(dat) <- c("target", "fips", "samplesize", "mse_ow0", "mse_ow1", "mse_cwm", "mse_sa", "p_ow0", "p_ow1", "p_cwm")
target <- sort(unique(dat$target))
samplesize <- sort(unique(dat$samplesize))
#### seperate four targets ####
rates <- array(NA, dim = c(length(target), length(samplesize), 3, 4)) # 3 weighting mehtods, 4 rates
crit <- 0.05
for(i in 1:length(target)){
  subdat <- dat[which(dat$target == target[i]), ]
  for(j in 1:length(samplesize)){
    subsubdat <- subdat[which(subdat$samplesize == samplesize[j]), ]
    rates[i, j, 1, 1] <- length(which(subsubdat$p_ow0 <= crit & subsubdat$mse_ow0 <= subsubdat$mse_sa))/nrow(subsubdat) # correct reject
    rates[i, j, 1, 2] <- length(which(subsubdat$p_ow0 <= crit & subsubdat$mse_ow0 >  subsubdat$mse_sa))/nrow(subsubdat) # incorrect reject
    rates[i, j, 1, 3] <- length(which(subsubdat$p_ow0 >  crit & subsubdat$mse_ow0 <  subsubdat$mse_sa))/nrow(subsubdat) # incorrect retain
    rates[i, j, 1, 4] <- length(which(subsubdat$p_ow0 >  crit & subsubdat$mse_ow0 >= subsubdat$mse_sa))/nrow(subsubdat) # correct retain
    rates[i, j, 2, 1] <- length(which(subsubdat$p_ow1 <= crit & subsubdat$mse_ow1 <= subsubdat$mse_sa))/nrow(subsubdat) # correct reject
    rates[i, j, 2, 2] <- length(which(subsubdat$p_ow1 <= crit & subsubdat$mse_ow1 >  subsubdat$mse_sa))/nrow(subsubdat) # incorrect reject
    rates[i, j, 2, 3] <- length(which(subsubdat$p_ow1 >  crit & subsubdat$mse_ow1 <  subsubdat$mse_sa))/nrow(subsubdat) # incorrect retain
    rates[i, j, 2, 4] <- length(which(subsubdat$p_ow1 >  crit & subsubdat$mse_ow1 >= subsubdat$mse_sa))/nrow(subsubdat) # correct retain
    rates[i, j, 3, 1] <- length(which(subsubdat$p_cwm <= crit & subsubdat$mse_cwm <= subsubdat$mse_sa))/nrow(subsubdat) # correct reject
    rates[i, j, 3, 2] <- length(which(subsubdat$p_cwm <= crit & subsubdat$mse_cwm >  subsubdat$mse_sa))/nrow(subsubdat) # incorrect reject
    rates[i, j, 3, 3] <- length(which(subsubdat$p_cwm >  crit & subsubdat$mse_cwm <  subsubdat$mse_sa))/nrow(subsubdat) # incorrect retain
    rates[i, j, 3, 4] <- length(which(subsubdat$p_cwm >  crit & subsubdat$mse_cwm >= subsubdat$mse_sa))/nrow(subsubdat) # correct retain
  }
}
# if aggregating all targets 
ave.rates <- apply(rates, c(2, 3, 4), mean)
values.OW0 <- t(ave.rates[ , 1, 1:2])
values.OW1 <- t(ave.rates[ , 2, 1:2])
values.CWM <- t(ave.rates[ , 3, 1:2])
# if seperate all targets
target_index <- 1
values.OW0 <- t(rates[target_index, , 1, 1:2])
values.OW1 <- t(rates[target_index, , 2, 1:2])
values.CWM <- t(rates[target_index, , 3, 1:2])
#### Alterntive way to combine results: combining all targets ####
rates <- array(NA, dim = c(length(samplesize), 3, 4)) # 2 weighting mehtods, 4 rates
crit <- 0.05
for(j in 1:length(samplesize)){
  subdat <- dat[which(dat$samplesize == samplesize[j]), ]
  rates[j, 1, 1] <- length(which(subdat$p_ow0 <= crit & subdat$mse_ow0 <= subdat$mse_sa))/nrow(subdat) # correct reject
  rates[j, 1, 2] <- length(which(subdat$p_ow0 <= crit & subdat$mse_ow0 >  subdat$mse_sa))/nrow(subdat) # incorrect reject
  rates[j, 1, 3] <- length(which(subdat$p_ow0 >  crit & subdat$mse_ow0 <  subdat$mse_sa))/nrow(subdat) # incorrect retain
  rates[j, 1, 4] <- length(which(subdat$p_ow0 >  crit & subdat$mse_ow0 >= subdat$mse_sa))/nrow(subdat) # correct retain
  rates[j, 2, 1] <- length(which(subdat$p_ow1 <= crit & subdat$mse_ow1 <= subdat$mse_sa))/nrow(subdat) # correct reject
  rates[j, 2, 2] <- length(which(subdat$p_ow1 <= crit & subdat$mse_ow1 >  subdat$mse_sa))/nrow(subdat) # incorrect reject
  rates[j, 2, 3] <- length(which(subdat$p_ow1 >  crit & subdat$mse_ow1 <  subdat$mse_sa))/nrow(subdat) # incorrect retain
  rates[j, 2, 4] <- length(which(subdat$p_ow1 >  crit & subdat$mse_ow1 >= subdat$mse_sa))/nrow(subdat) # correct retain
  rates[j, 3, 1] <- length(which(subdat$p_cwm <= crit & subdat$mse_cwm <= subdat$mse_sa))/nrow(subdat) # correct reject
  rates[j, 3, 2] <- length(which(subdat$p_cwm <= crit & subdat$mse_cwm >  subdat$mse_sa))/nrow(subdat) # incorrect reject
  rates[j, 3, 3] <- length(which(subdat$p_cwm >  crit & subdat$mse_cwm <  subdat$mse_sa))/nrow(subdat) # incorrect retain
  rates[j, 3, 4] <- length(which(subdat$p_cwm >  crit & subdat$mse_cwm >= subdat$mse_sa))/nrow(subdat) # correct retain
}
values.OW0 <- t(rates[ , 1, 1:2])
values.OW1 <- t(rates[ , 2, 1:2])
values.CWM <- t(rates[ , 3, 1:2])
#### visualization ####
plotdata.OW0 <- colSums(values.OW0)
plotdata.OW1 <- colSums(values.OW1)
plotdata.CWM <- colSums(values.CWM)
xaxisname <- as.character(samplesize)

setEPS()
postscript("barchart_covid_deaths.eps", width = 6, height = 4.5)
par(mfrow = c(1, 2))
barplot(plotdata.OW1, names.arg = xaxisname, col = "white", 
        xlab = "Sample Size", ylab = "Proportion", main = "OW vs. SA", ylim = c(0, 1))
items <- c("Total rejections")
legend("bottomright", items, density = c(0), bg = "white")
barplot(plotdata.CWM, names.arg = xaxisname, col = "white", 
        xlab = "Sample Size", ylab = "Proportion", main = "CWM vs. SA", ylim = c(0, 1))
items <- c("Total rejections")
legend("bottomright", items, density = c(0), bg = "white")
dev.off()
##################################################################
############ Empirical Aalysis: Keck's study #####################
##################################################################
setwd("insert the path to the related results data. ")
dat <- read.csv("res_keck.csv", header = T, stringsAsFactors = F)
colnames(dat) <- c("condition", "draw_people", "draw_judgments", "samplesize", "mse_ow0", "mse_ow1", "mse_cwm", "mse_sa", "p_ow0", "p_ow1", "p_cwm")
table(dat$samplesize)
samplesize <- sort(unique(dat$samplesize))
conditions <- sort(unique(dat$condition))
#### seperate conditions and combining all random draws ####
statdat <- unique(dat[ , c(1, 2, 4)])
num_runs <- c()
rates <- matrix(NA, nrow = nrow(statdat), ncol = 3*4) # 3 weighting methods * 4 rates 
crit <- 0.05
for(i in 1:nrow(statdat)){
  subdat <- dat[which(dat$condition == statdat$condition[i] & dat$draw_people == statdat$draw_people[i] & dat$samplesize == statdat$samplesize[i]), ]
  num_runs[i] <- nrow(subdat)
  rates[i, 1] <- length(which(subdat$p_ow0 <= crit & subdat$mse_ow0 <= subdat$mse_sa))/nrow(subdat) # correct reject
  rates[i, 2] <- length(which(subdat$p_ow0 <= crit & subdat$mse_ow0 >  subdat$mse_sa))/nrow(subdat) # incorrect reject
  rates[i, 3] <- length(which(subdat$p_ow0 >  crit & subdat$mse_ow0 <  subdat$mse_sa))/nrow(subdat) # incorrect retain
  rates[i, 4] <- length(which(subdat$p_ow0 >  crit & subdat$mse_ow0 >= subdat$mse_sa))/nrow(subdat) # correct retain
  rates[i, 5] <- length(which(subdat$p_ow1 <= crit & subdat$mse_ow1 <= subdat$mse_sa))/nrow(subdat) # correct reject
  rates[i, 6] <- length(which(subdat$p_ow1 <= crit & subdat$mse_ow1 >  subdat$mse_sa))/nrow(subdat) # incorrect reject
  rates[i, 7] <- length(which(subdat$p_ow1 >  crit & subdat$mse_ow1 <  subdat$mse_sa))/nrow(subdat) # incorrect retain
  rates[i, 8] <- length(which(subdat$p_ow1 >  crit & subdat$mse_ow1 >= subdat$mse_sa))/nrow(subdat) # correct retain
  rates[i, 9] <- length(which(subdat$p_cwm <= crit & subdat$mse_cwm <= subdat$mse_sa))/nrow(subdat) # correct reject
  rates[i,10] <- length(which(subdat$p_cwm <= crit & subdat$mse_cwm >  subdat$mse_sa))/nrow(subdat) # incorrect reject
  rates[i,11] <- length(which(subdat$p_cwm >  crit & subdat$mse_cwm <  subdat$mse_sa))/nrow(subdat) # incorrect retain
  rates[i,12] <- length(which(subdat$p_cwm >  crit & subdat$mse_cwm >= subdat$mse_sa))/nrow(subdat) # correct retain
}
statdat <- cbind(statdat, num_runs, rates)
plotdata.OW0 <- mean(statdat[ , 5]) + mean(statdat[ , 6])
plotdata.OW1 <- mean(statdat[ , 9]) + mean(statdat[ , 10])
plotdata.CWM <- mean(statdat[ , 13]) + mean(statdat[ , 14])
#### combining all conditions and random draws of people ####
values.OW0 <- c()
values.OW1 <- c()
values.CWM <- c()
for(i in 1:length(samplesize)){
  tempdat <- statdat[which(statdat$samplesize == samplesize[i]), ]
  values.OW0 <- rbind(values.OW0, colMeans(tempdat[ , 5:6]))
  values.OW1 <- rbind(values.OW1, colMeans(tempdat[ , 9:10]))
  values.CWM <- rbind(values.CWM, colMeans(tempdat[ , 13:14]))
}
plotdata.OW0 <- rowSums(values.OW0)
plotdata.OW1 <- rowSums(values.OW1)
plotdata.CWM <- rowSums(values.CWM)
## visualization 
xaxisname <- as.character(samplesize)

setEPS()
postscript("barchart_keck.eps", width = 4, height = 4.5)
par(mfrow = c(1, 2))
barplot(plotdata.OW1, names.arg = xaxisname, col = "white", 
        xlab = "Sample Size", ylab = "Proportion", main = "OW vs. SA", ylim = c(0, 1))
items <- c("Total rejections")
legend("bottomright", items, density = c(0), bg = "white")
barplot(plotdata.CWM, names.arg = xaxisname, col = "white", 
        xlab = "Sample Size", ylab = "Proportion", main = "CWM vs. SA", ylim = c(0, 1))
items <- c("Total rejections")
legend("bottomright", items, density = c(0), bg = "white")dev.off()


#################################################
############ Power Analysis #####################
#################################################
setwd("insert the path to the related results data. ")
dat <- read.csv("power_inf.csv", header = F, stringsAsFactors = F)
#dat <- read.csv("power_unemp.csv", header = F, stringsAsFactors = F)
colnames(dat) <- c("samplesize", "round", paste("p", 1:7, sep = ""))

samplesize <- sort(unique(dat$samplesize))
table(dat$samplesize)
crit.test <- c(0.001, 0.01, 0.05, 0.1)
power_value <- array(NA, dim = c(length(samplesize), 7, length(crit.test)))
for(i in 1:length(samplesize)){
  subdat <- dat[which(dat$samplesize == samplesize[i]), ]
  for(j in 1:length(crit.test)){
    for(k in 1:7){
      tempvalue <- subdat[ , k+2]
      tempvalue[which(tempvalue == "No need to compare!")] <- NA
      tempvalue <- as.numeric(tempvalue)
      tempvalue <- na.omit(tempvalue)
      power_value[i, k, j] <- length(which(tempvalue < crit.test[j]))/length(tempvalue)
    }
  }
}
# record all results 
p.crit <- 0.001
write.csv(power_value[ , , which(crit.test == p.crit)], "res_inf_001.csv", row.names = F, col.names = F)
# get the target sample size 
p.crit <- 0.001
power.threshold <- 0.9
for(k in 1:7){
  samplesize_index <- which(power_value[ , k, which(crit.test == p.crit)] > power.threshold)
  print(samplesize[samplesize_index[1:5]])
}
# plot the power curves
#samplesize <- c(seq(16, 26, by = 2), seq(28, 68, by = 1), seq(70, 80, by = 2))
plotdata <- read.csv("res_inf_001.csv", header = T, stringsAsFactors = F)
#samplesize <- c(seq(16, 98, by = 2), seq(100, 128, by = 4))
samplesize <- c(16, 18, 20, seq(24, 98, by = 2), seq(100, 128, by = 4))

setEPS()
postscript("power_inf_001.eps", width = 5, height = 5)
plot(samplesize, plotdata[ , 2], type = "l", lwd = 2, col = "black", ylim = c(0, 1), 
     xlab = "Sample Size", ylab = "Power")
lines(samplesize, plotdata[ , 4], type = "l", lwd = 2, col = "purple")
lines(samplesize, plotdata[ , 6], type = "l", lwd = 2, col = "blue")
lines(samplesize, plotdata[ , 8], type = "l", lwd = 2, col = "green")
lines(samplesize, plotdata[ , 10], type = "l", lwd = 2, col = "orange")
lines(samplesize, plotdata[ , 12], type = "l", lwd = 2, col = "yellow")
lines(samplesize, plotdata[ , 14], type = "l", lwd = 2, col = "red")
abline(h = 0.8, lwd = 1, lty = 2)
legend("bottomright", c("OW", "TOP3", "RP", "CWM", "SSIN", "SSDE", "LAS"), col = c("black", "purple", "blue", "green", "orange", "yellow", "red"), lwd = 2)
dev.off()
