# load the path
# update the dictionary path to test algorithm/ weighting functions/ cross validation
setwd("add the dictonary path to the file here")
# load all functions 
source("0algorithm/test algorithm.R")
source("0algorithm/weighting functions.R")

#####################################
######### Data Analysis #############
#####################################
infdata <- read.csv("emp_filtered.csv", header = T, sep = ",")
infdata <- infdata[ , -1]

M <- 8 
X <- as.matrix(infdata[ , 4:11]) # this is for inflation rate data
Y <- infdata[ , 3]
train.num <- round(nrow(X) * 0.6, 0) # 80% data as training set and 20% as testing set 

round1 <- 5 # this is for 80% data selection
round2 <- 100 # this is for sampel size n's selection 
# sample size to be changed 
n.vector <- seq(16, 128, by = 16)
# settings for hypothesis test algorithm
reps1 <- 500
depth <- 20
reps2.fast <- 10000
reps2.slow <- 10000

for(i in 1:round1){
  set.seed(i*284)
  # select the training and testing data in each round 1
  train.index <- sample(1:nrow(X), size = train.num, replace = FALSE)
  train.index <- sort(train.index)
  X.train <- X[train.index, ]
  X.test <- X[-train.index, ]
  Y.train <- Y[train.index]
  Y.test <- Y[-train.index]
  
  # calculate the MSE of simple average method in the testing data set 
  Optweight0 <- optw(X.train, Y.train)
  Optweight1 <- optw.sum1pos(X.train, Y.train)
  CWMweight <- CWMall(X.train, Y.train)
  # get the MSE
  OW0.test <- X.test %*% Optweight0
  OW1.test <- X.test %*% Optweight1
  CWM.test <- X.test %*% CWMweight
  SA.test  <- X.test %*% rep(1/M, M)
  MSE <- c(mean((OW0.test - Y.test)^2), mean((OW1.test - Y.test)^2), 
           mean((CWM.test - Y.test)^2), mean((SA.test - Y.test)^2))
  write.table(t(c(i, MSE, Optweight0, Optweight1, CWMweight)), "MSEout_unemp.csv", sep = ",", append = TRUE, row.names = FALSE, col.names = FALSE)
  
  for(j in 1:round2){
    for(l in 1:length(n.vector)){
      n <- n.vector[l]
      tryCatch({
        set.seed(i*k*j*195)
        # randomly select the training data according to the sample size 
        sampleindex <- sample(1:nrow(X.train), size = n, replace = FALSE)
        x.train <- X.train[sampleindex, ]
        y.train <- Y.train[sampleindex]
        xy_all.train <- cbind(x.train, y.train)
        ##### step1: observe bias and covariance matrix #####
        optweight0 <- optw(x.train, y.train)
        optweight1 <- optw.sum1pos(x.train, y.train)
        cwmweight <- CWMall(x.train, y.train)
        # get the MSE
        ow0.test <- X.test %*% optweight0
        ow1.test <- X.test %*% optweight1
        cwm.test <- X.test %*% cwmweight
        mse <- c(mean((ow0.test - Y.test)^2), mean((ow1.test - Y.test)^2), mean((cwm.test - Y.test)^2))
        # run the test algorithm 
        mu_xy.obs <- colMeans(xy_all.train)
        Sigma_xy.obs <- cov(xy_all.train)
        weights <- cbind(optweight0, optweight1, cwmweight)
        #### calculate pvalue ####
        pvalue <- c()
        weights.remove.index <- which(colSums(abs(weights - rep(0.2, M))) == 0)
        if(length(weights.remove.index) > 0){
          pvalue[weights.remove.index] <- "No need to compare!"
          weights.retain.index <- c(1:ncol(weights))[-weights.remove.index]
          weights.retain <- weights[ , weights.retain.index]
        } else {
          weights.retain.index <- c(1:ncol(weights))
          weights.retain <- weights
        }
        equalweight <- rep(1/M, M)
        cons <- ese.oneset(mu_xy.obs, Sigma_xy.obs, cbind(weights.retain, equalweight))
        cons.wa <- cons[1:ncol(weights.retain)]
        cons.sa <- cons[length(cons)]
        weights.retain.fast.index <- weights.retain.index[which(cons.wa >= cons.sa)]
        weights.retain.slow.index <- weights.retain.index[which(cons.wa <  cons.sa)]
        weights.retain.fast <- weights[ , weights.retain.fast.index]
        weights.retain.slow <- weights[ , weights.retain.slow.index]
        if(length(weights.retain.fast.index) > 0){
          # for weights.retain.fast
          mu_xy.star <- mu_xy.obs
          Sigma_xy.star <- Sigma_xy.obs
          density.obs <- compdensity.oneset(n, M+1, mu_xy.star, Sigma_xy.star, mu_xy.obs, Sigma_xy.obs) # return the log objective function
          pvalue[weights.retain.fast.index] <- calPvalue(reps2 = reps2.fast, mu.star = mu_xy.star, sigma.star = Sigma_xy.star, density.obs = density.obs, n = n)
        }
        if(length(weights.retain.slow.index) == 1){
          mu_xy.star <- rep(mu_xy.obs[M+1], M+1)
          Sigma_xy.star <- getEqualCov(Sigma_xy.obs)
          density.obs.old <- compdensity.oneset(n, M+1, mu_xy.star, Sigma_xy.star, mu_xy.obs, Sigma_xy.obs)
          solvestar <- findStar(reps1 = reps1, mu.obs = mu_xy.obs, sigma.obs = Sigma_xy.obs, obj.star = density.obs.old, n = n, 
                                weight = cbind(weights.retain.slow, weights.retain.slow), depth = depth)
          mu_xy.star <- solvestar$mustar
          Sigma_xy.star <- solvestar$sigmastar
          density.obs <- compdensity.forstar(n, M+1, t(mu_xy.star), Sigma_xy.star, mu_xy.obs, Sigma_xy.obs)
          pv <- calPvalue(reps2 = reps2.slow, mu.star = mu_xy.star[ , 1], sigma.star = Sigma_xy.star[ , , 1], density.obs = density.obs[1], n = n)
          pvalue[weights.retain.slow.index] <- pv
        }
        if(length(weights.retain.slow.index) > 1){
          # for weights.retain.slow 
          ## find the unique weights vectors 
          dup.judge <- duplicated(t(weights.retain.slow))
          unique.index <- which(dup.judge == FALSE)
          dup.index <- which(dup.judge == TRUE)
          target.index <- 1:ncol(weights.retain.slow)
          for(k in 1:length(unique.index)){
            temp <- which(colSums(abs(weights.retain.slow - weights.retain.slow[ , unique.index[k]])) == 0)
            if(length(temp) > 1){
              target.index[temp] <- unique.index[k]
            } else {
              target.index[unique.index[k]] <- unique.index[k]
            }
          }
          weights.retain.slow.unique <- t(unique(t(weights.retain.slow)))
          ## get the p-value for unique weights 
          if(length(unique.index) == 1){
            mu_xy.star <- rep(mu_xy.obs[M+1], M+1)
            Sigma_xy.star <- getEqualCov(Sigma_xy.obs)
            density.obs.old <- compdensity.oneset(n, M+1, mu_xy.star, Sigma_xy.star, mu_xy.obs, Sigma_xy.obs)
            solvestar <- findStar(reps1 = reps1, mu.obs = mu_xy.obs, sigma.obs = Sigma_xy.obs, obj.star = density.obs.old, n = n, 
                                  weight = cbind(weights.retain.slow.unique, weights.retain.slow.unique), depth = depth)
            mu_xy.star <- solvestar$mustar
            Sigma_xy.star <- solvestar$sigmastar
            density.obs <- compdensity.forstar(n, M+1, t(mu_xy.star), Sigma_xy.star, mu_xy.obs, Sigma_xy.obs)
            pv <- calPvalue(reps2 = reps2.slow, mu.star = mu_xy.star[ , 1], sigma.star = Sigma_xy.star[ , , 1], density.obs = density.obs[1], n = n)
            pvalue[weights.retain.slow.index] <- pv
          } else {
            mu_xy.star <- rep(mu_xy.obs[M+1], M+1)
            Sigma_xy.star <- getEqualCov(Sigma_xy.obs)
            density.obs.old <- compdensity.oneset(n, M+1, mu_xy.star, Sigma_xy.star, mu_xy.obs, Sigma_xy.obs)
            solvestar <- findStar(reps1 = reps1, mu.obs = mu_xy.obs, sigma.obs = Sigma_xy.obs, obj.star = density.obs.old, n = n, 
                                  weight = weights.retain.slow.unique, depth = depth)
            mu_xy.star <- solvestar$mustar
            Sigma_xy.star <- solvestar$sigmastar
            density.obs <- compdensity.forstar(n, M+1, t(mu_xy.star), Sigma_xy.star, mu_xy.obs, Sigma_xy.obs)
            pv <- sapply(1:ncol(mu_xy.star), function(x) calPvalue(reps2 = reps2.slow, 
                                                                   mu.star = mu_xy.star[ , x], 
                                                                   sigma.star = Sigma_xy.star[ , , x],
                                                                   density.obs = density.obs[x], n = n))
            pv.full <- c()
            for(k in 1:length(unique.index)){
              pv.full[which(target.index == unique.index[k])] <- pv[k]
            }
            pvalue[weights.retain.slow.index] <- pv.full
          }
        }
        #### output: pvalue ####
        write.table(t(c(i, n, j, mse, pvalue)), file = "res_unemp.csv", sep = ",", append = TRUE, row.names = FALSE, col.names = FALSE)
        print(paste(i, "-", n, "-", j, sep = ""))
      }, error = function(e){
        cat(paste("One Error For", i, "-", n, "-", j, sep = ""))
      }
      )
    }
  }
}
