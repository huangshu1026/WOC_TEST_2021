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


############################################################
######### Data Analysis: COVID-19 forecasting  #############
############################################################
data_case <- read.csv("all_data_case_county.csv", header = T, stringsAsFactors = F)
data_case <- data_case[ , -1] #1:7, 66:70 are indicators
na_ratio <- apply(data_case, 2, function(x) length(which(is.na(x)))/nrow(data_case))
# grouped data 
targets <- unique(data_case$target)
fips <- unique(data_case$fips)

na_threshold <- 0.1 # removing missing data
n_vec <- c(12, 16, 20, 24, 28) # n = c(12, 16, 20, 24, 28)
for(i in 1:4){
  for(j in 501:1000){
    subdata <- data_case[which(data_case$target == targets[i] & data_case$fips == fips[j]), ]
    na_ratio <- apply(subdata, 2, function(x) length(which(is.na(x)))/nrow(subdata))
    subdata <- subdata[ , -which(na_ratio >= na_threshold)]
    subdata <- na.omit(subdata)
    print(paste("#Observation: ", nrow(subdata), sep = ""))
    print(paste("#Forecasters: ", ncol(subdata)-11, sep = ""))
    for(r in 1:length(n_vec)){
      n <- n_vec[r]
      set.seed(r*i*j*324)
      tryCatch({
        if(nrow(subdata) >= 32 & ncol(subdata) >= 14){
          data_train <- subdata[1:n, ]
          data_test  <- subdata[(n+1):nrow(subdata), ]
          x.train <- as.matrix(data_train[ , c(8:(ncol(data_train)-4))])
          y.train <- as.vector(data_train$inc_case.y)
          x.test <- as.matrix(data_test[ , c(8:(ncol(data_test)-4))])
          y.test <- as.vector(data_test$inc_case.y)
          M <- ncol(x.train)
          xy_all.train <- cbind(x.train, y.train)
          # calcualte two weights
          optweight0 <- optw(x.train, y.train)
          optweight1 <- optw.sum1pos(x.train, y.train) # update weighted average here
          cwmweight <- CWMall(x.train, y.train)
          # get the MSE
          ow0.test <- x.test %*% optweight0
          ow1.test <- x.test %*% optweight1
          cwm.test <- x.test %*% cwmweight
          sa.test <- x.test %*% rep(1/M, M)
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
          #### output: pvalue ####
          write.table(t(c(i, j, n, mse, pvalue)), 
                      file = "res_covid_cases.csv", 
                      sep = ",", append = TRUE, row.names = FALSE, col.names = FALSE)
          print(paste(i, "-", j, "-", n, sep = ""))
        }
      }, error = function(e){ 
        cat(paste("Error for", i, "-", j, sep = ""))
      })
    }
  }
}

          