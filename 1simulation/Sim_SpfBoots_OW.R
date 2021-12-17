# update the dictionary path to test algorithm/ weighting functions/ cross validation
setwd("add the dictonary path to the file here")

# load all functions 
source("0algorithm/test algorithm.R")
source("0algorithm/weighting functions.R")
source("0algorithm/cross validation.R")

#####################################################
######### Data Analysis: Inflation Rate #############
#####################################################
## Prepare the data 
infdata <- read.csv("inf_filtered.csv", header = T, sep = ",")
infdata <- infdata[ , -1]
M <- 5
X <- as.matrix(infdata[ , 4:8]) # this is for inflation rate data
Y <- infdata[ , 3]
## Prepare the data: unemployment rate data 
# unempdata <- read.csv("unemp_filtered.csv", header = T, sep = ",")
# unempdata <- unempdata[ , -1]
# M <- 5
# X <- as.matrix(unempdata[ , 4:8]) # this is for inflation rate data
# Y <- unempdata[ , 3]

## Estimate the bias and covariance matrix
XY_all <- cbind(X, Y)
mu_xy.emp <- colMeans(XY_all)
Sigma_xy.emp <- cov(XY_all)

# parameter setting for hypothesis testing algorithm 
reps1 <- 500
depth <- 20
reps2.fast <- 10000
reps2.slow <- 10000

## parameters in the simulation 
# sample size to be changed 
n.vector <- 2^seq(4, 8, by = 1)
n.test <- 1000
draws <- 1000
for(j in 1:draws){
  for(i in 1:length(n.vector)){
    tryCatch({
      n <- n.vector[i]
      set.seed(i*j*785)
      # generate the true state
      XY_all.true <- rmvn(100, mu_xy.emp, Sigma_xy.emp)
      mu_xy.true <- colMeans(XY_all.true)
      Sigma_xy.true <- cov(XY_all.true)
      # generate a set of judgments 
      XY_all.train <- rmvn(n, mu_xy.true, Sigma_xy.true)
      XY_all.test  <- rmvn(n.test, mu_xy.true, Sigma_xy.true)
      x.train <- XY_all.train[ , 1:M]
      y.train <- XY_all.train[ , (M+1)]
      x.test  <- XY_all.test[ , 1:M]
      y.test  <- XY_all.test[ , (M+1)]
      # estimated parameters 
      mu_xy.obs <- colMeans(XY_all.train)
      Sigma_xy.obs <- cov(XY_all.train)
      # estimating all weights
      optweight1 <- optw.sum1(x.train, y.train) #1
      optweight2 <- optw.sum1pos(x.train, y.train) #1
      weights <- cbind(optweight1, optweight2)
      # Run hypothesis test for p-value 
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
      # cross validation 
      cv <- CrossValidation_subset(x.train, y.train, k.fold = 0.8, round.num = 10) # return a 9-dimension vector
      cv.error <- cv - cv[3]
      cv.wa.win <- cv.error[1:2]
      # really out-of-sample test 
      AllWeights <- cbind(weights, equalweight)
      mean.mse.error <- colMeans((x.test %*% AllWeights - y.test)^2)
      mean.mse.error <- mean.mse.error - mean.mse.error[3]
      mean.mse.error <- mean.mse.error[1:2]
      MSE.out <- colMeans((x.test %*% AllWeights - y.test)^2)
      
      mse.all <- (x.test %*% AllWeights - y.test)^2
      MSE.out.ttest.pvalue <- c()
      for(k in which(mean.mse.error != 0)){
        ttest.mse <- t.test(mse.all[ , k], mse.all[ , 3], paired = TRUE)
        MSE.out.ttest.pvalue[k] <- ttest.mse$p.value
      }
      write.table(t(c(n, j, pvalue, cv.wa.win, MSE.out, MSE.out.ttest.pvalue)), "InfStates_OW.csv", sep = ",", append = TRUE, row.names = FALSE, col.names = FALSE)
      #write.table(t(c(n, j, pvalue, cv.wa.win, MSE.out, MSE.out.ttest.pvalue)), "UnempStates_OW.csv", sep = ",", append = TRUE, row.names = FALSE, col.names = FALSE)
      print(paste(n, "-", j, sep = ""))  
    }, error = function(e){
      cat("One Error Found!")
    })
  }
}

