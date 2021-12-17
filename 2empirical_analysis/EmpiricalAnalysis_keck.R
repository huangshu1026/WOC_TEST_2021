# load the path
# update the dictionary path to test algorithm/ weighting functions/ cross validation
setwd("add the dictonary path to the file here")
# load all functions 
source("0algorithm/test algorithm.R")
source("0algorithm/weighting functions.R")

# settings for hypothesis test algorithm
reps1 <- 500
depth <- 20
reps2.fast <- 10000
reps2.slow <- 10000

# load the data 
#install.packages("readxl")
library(readxl)
exp_dat <- read_excel("KeckDataStudy1.xlsx")

# Start the anlaysis 
M <- 5
n <- 16 # can change the sample size c(16, 20, 24, 28, 32)
draws_people <- 10
draws_judgments <- 100
condition <- unique(exp_dat$con) # 3 conditions
for(i in 1:length(condition)){ # index of condition
  for(j in 1:draws_people){ # index of the random draw
    # prepare the subset of data for analysis 
    set.seed(i*j*123)
    subdata <- exp_dat[which(exp_dat$con == condition[i]), ]
    subdata <- na.omit(subdata)
    subdata <- subdata[sample(1:nrow(subdata), M, replace = FALSE), ]
    x.all <- t(subdata[ , 2:41])
    y.all <- t(subdata[1, 42:81])
    for(l in 1:draws_judgments){
      # split the training data and testing data 
      set.seed(l*522)
      index_train <- sample(1:nrow(x.all), n, replace = FALSE)
      index_test <- c(1:nrow(x.all))[-index_train]
      x.train <- as.matrix(x.all[index_train, ])
      y.train <- as.vector(y.all[index_train])
      x.test <- as.matrix(x.all[index_test, ])
      y.test <- as.vector(y.all[index_test])
      xy_all.train <- cbind(x.train, y.train)
      tryCatch({
        # calculate the p-value for two weighted averages 
        optweight0 <- optw(x.train, y.train)
        optweight1 <- optw.sum1pos(x.train, y.train) # update weighted average here
        cwmweight <- CWMall(x.train, y.train)
        # get the MSE
        ow0.test <- x.test %*% optweight0
        ow1.test <- x.test %*% optweight1
        cwm.test <- x.test %*% cwmweight
        sa.test  <- x.test %*% rep(1/M, M)
        mse <- c(mean((ow0.test - y.test)^2), mean((ow1.test - y.test)^2), 
                 mean((cwm.test - y.test)^2), mean((sa.test - y.test)^2))
        
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
        
        #### output: pvalue  ####
        # record all the results
        write.table(t(c(i, j, l, n, mse, pvalue)), 
                    file = paste("res_keck.csv", sep = ""), 
                    sep = ",", append = TRUE, row.names = FALSE, col.names = FALSE)
        print(paste(i, "-", j, "-", l, sep = ""))
      }, error = function(e){
        cat(paste("Error for", i, "-", j, "-", l, sep = ""))
      })
    }
  }
}
