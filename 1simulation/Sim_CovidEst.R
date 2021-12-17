# load the path
# update the dictionary path to test algorithm/ weighting functions/ cross validation
setwd("add the dictonary path to the file here")
# load all functions 
source("0algorithm/test algorithm.R")
source("0algorithm/weighting functions.R")
source("0algorithm/cross validation.R")

# parameter setting for hypothesis testing algorithm 
reps1 <- 500
depth <- 20
reps2.fast <- 10000
reps2.slow <- 10000

##################################################################
######### Data Analysis: Covid New Cases Forecasting #############
##################################################################
## Prepare the data 
data_case <- read.csv("all_data_case_county.csv", header = T, stringsAsFactors = F)
data_case <- data_case[ , -1] #1:7, 66:70 are indicators
targets <- unique(data_case$target)
fips <- unique(data_case$fips)
na_threshold <- 0.1
n.vector <- 2^seq(4, 8, by = 1)
n.test <- 1000
for(i in 1:4){
  for(j in 1:500){
    dat <- data_case[which(data_case$target == targets[i] & data_case$fips == fips[j]), ]
    na_ratio <- apply(dat, 2, function(x) length(which(is.na(x)))/nrow(dat))
    dat <- dat[ , -which(na_ratio >= na_threshold)]
    dat <- na.omit(dat)
    print(paste("#Observation: ", nrow(dat), sep = ""))
    print(paste("#Forecasters: ", ncol(dat)-11, sep = ""))
    if(nrow(dat) >= 32 & ncol(dat) >= 14){
      M <- ncol(dat)-11
      X <- as.matrix(dat[ , c(8:(ncol(dat)-4))]) # this is for inflation rate data
      Y <- dat$inc_case.y
      ## Estimate the bias and covariance matrix
      XY_all <- cbind(X, Y)
      mu_xy.emp <- colMeans(XY_all)
      Sigma_xy.emp <- cov(XY_all)
      for(r in 1:length(n.vector)){
        tryCatch({
          n <- n.vector[r]
          set.seed(i*j*r*234)
          # generate a set of judgments with sample size n 
          XY_all.train <- rmvn(n, mu_xy.emp, Sigma_xy.emp)
          XY_all.test  <- rmvn(n.test, mu_xy.emp, Sigma_xy.emp)
          x.train <- XY_all.train[ , 1:M]
          y.train <- XY_all.train[ , (M+1)]
          x.test  <- XY_all.test[ , 1:M]
          y.test  <- XY_all.test[ , (M+1)]
          # estimated parameters 
          mu_xy.obs <- colMeans(XY_all.train)
          Sigma_xy.obs <- cov(XY_all.train)
          # estimating all weights
          optweight <- optw(x.train, y.train) #1
          top3equalweight <- topNew(x.train, y.train, N = 3) #2
          rankweight <- RankPerformance(x.train, y.train) #3
          cwmweight <- CWMall(x.train, y.train) #4
          SSweights <- SeqSearch(x.train, y.train)
          ssinweight <- SSweights$weight.SSin #5
          ssdeweight <- SSweights$weight.SSde #6
          regcovweight <- RegCovOptw(x.train, y.train) #7
          reglassoweight <- ConstrENetZERO(x.train, y.train, addl1 = TRUE, addl2 = FALSE, l1.seq = 2^(seq(0, 10, by = 0.2)), l2.seq = NULL) #8
          weights <- cbind(optweight, top3equalweight, rankweight, cwmweight, 
                           ssinweight, ssdeweight, regcovweight, reglassoweight)
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
          cv <- CrossValidation(x.train, y.train, k.fold = 0.8, round.num = 10) # return a 9-dimension vector
          cv.error <- cv - cv[9]
          cv.wa.win <- cv.error[1:8]
          # really out-of-sample test 
          AllWeights <- cbind(weights, equalweight)
          mean.mse.error <- colMeans((x.test %*% AllWeights - y.test)^2)
          mean.mse.error <- mean.mse.error - mean.mse.error[9]
          mean.mse.error <- mean.mse.error[1:8]
          MSE.out <- colMeans((x.test %*% AllWeights - y.test)^2)
          
          mse.all <- (x.test %*% AllWeights - y.test)^2
          MSE.out.ttest.pvalue <- c()
          for(k in which(mean.mse.error != 0)){
            ttest.mse <- t.test(mse.all[ , k], mse.all[ , 9], paired = TRUE)
            MSE.out.ttest.pvalue[k] <- ttest.mse$p.value
          }
          write.table(t(c(i, j, n, pvalue, cv.wa.win, MSE.out, MSE.out.ttest.pvalue)), "CovidStates.csv", sep = ",", append = TRUE, row.names = FALSE, col.names = FALSE)
          print(paste(i, "-", j, "-", n, sep = ""))
        }, error = function(e){
          cat("One Error Found!")
        })
      }
    }
  }
}