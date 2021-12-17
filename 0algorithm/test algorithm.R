##########################################################################
########## This file calcualtes p-value for given weights ################
##########################################################################
##### load the library #####
if (!require('mgcv')) install.packages('mgcv'); library('mgcv')
if (!require('matrixcalc')) install.packages('matrixcalc'); library('matrixcalc')
if (!require('parallel')) install.packages('parallel'); library('parallel')
if (!require('MASS')) install.packages('MASS'); library('MASS')
if (!require('MBESS')) install.packages('MBESS'); library('MBESS')
if (!require('foreach')) install.packages('foreach'); library('foreach')
if (!require('doParallel')) install.packages('doParallel'); library('doParallel')
if (!require('CholWishart')) install.packages('CholWishart'); library('CholWishart')

##### marginal functions for calculating p-values #####
getEqualCov <- function(Sigma_xy){
  M <- nrow(Sigma_xy) - 1
  variance <- diag(Sigma_xy)
  variance_new <- c(rep(mean(variance[1:M]), M), variance[M+1])
  Corr_xy <- cov2cor(Sigma_xy)
  Corr_xx <- Corr_xy[1:M, 1:M]
  cor1 <- mean(Corr_xx[upper.tri(Corr_xx)])
  cor2 <- mean(Corr_xy[1:M, M+1])
  Corr_xx.new <- matrix(rep(cor1, M*M), nrow = M, ncol = M)
  diag(Corr_xx.new) <- 1
  Corr_xy.new <- as.matrix(cbind(rbind(Corr_xx.new, rep(cor2, M)), c(rep(cor2, M), 1)))
  if(is.positive.definite(Corr_xy.new) == FALSE){
    Corr_xy.new <- nearPD(Corr_xy.new)$mat
  }
  Sigma_xy.new <- cor2cov(Corr_xy.new, sqrt(variance_new))
  return(Sigma_xy.new)
}
ese.oneset <- function(mu_xy, Sigma_xy, weights){
  M <- length(mu_xy)-1
  mu_x <- mu_xy[1:M]
  mu_y <- mu_xy[(M+1)]
  sigma_xx <- Sigma_xy[1:M, 1:M]
  sigma_xy <- Sigma_xy[1:M, (M+1)]
  sigma_y2 <- Sigma_xy[(M+1), (M+1)]
  term1 <- (t(mu_x)%*%weights - mu_y)^2 
  term2 <- diag(t(weights)%*%sigma_xx%*%weights)
  term3 <- 2*t(weights)%*%sigma_xy
  term4 <- sigma_y2
  res <- as.vector(term1) + as.vector(term2) - as.vector(term3) + term4
  return(res)
}
ese <- function(mu_xy, Sigma_xy, weights){
  # mu = d * M, sigma = M * M * d, weights = M * k
  # return = d * k
  # this function calculates the expected squared errors of simple average
  # given statistics (bias, covariance matrix)
  # mu and sigma should be the true one instead of observed ones 
  # observed bias and covariance matrix are used to calculate optimal weights (w)
  d <- dim(Sigma_xy)[3]
  res <- matrix(unlist(sapply(1:d, function(x) ese.oneset(mu_xy[x, ], Sigma_xy[ , , x], weights))), nrow = d, ncol = ncol(weights), byrow = TRUE)
  return(res)
}

genMu <- function(mu.old, mu.obs, reps1){
  # this function generates the random bias vector between old one to the observed one
  M <- length(mu.old)
  w <- array(runif(M*reps1, 0, 1), dim = c(M, 1, reps1))
  mu.new <- w * mu.old + (1-w) * mu.obs
  return(mu.new)
}
genMu.once <- function(mu.old, mu.obs){
  # mu.old = d * M
  # this function generates the random bias vector between old one to the observed one
  M <- length(mu.obs)
  d <- nrow(mu.old)
  w <- matrix(runif(d*M, 0, 1), nrow = d, ncol = M)
  mu.new <- w * mu.old + (1-w) * matrix(rep(mu.obs, d), nrow = d, ncol = M, byrow = TRUE)
  return(mu.new)
}
genMu.random <- function(mu.old, dist){
  # this function generates the random bias around the bias of interest 
  M <- length(mu.old)
  mu.new <- mu.old +  runif(M, -dist, dist)
  return(mu.new)
}
genSigma <- function(sigma.old, sigma.obs){
  # this function generates random covariance matrix: convex combination of sigma.old and sigma.obs
  M <- nrow(sigma.old)
  w <- matrix(runif(M*M, 0, 1), M, M)
  w[upper.tri(w)] <- t(w)[upper.tri(t(w))]
  sigma.new <- w*sigma.old + (1-w)*sigma.obs
  return(sigma.new)
}
genSigma.random <- function(sigma.old, dist){
  M <- nrow(sigma.old)
  mat <- matrix(runif(M*M, -dist, dist), M, M)
  mat[upper.tri(mat)] <- t(mat)[upper.tri(t(mat))]
  sigma.new <- sigma.old + mat
  return(sigma.new)
}
compdensity.oneset <- function(n, M, mu.star, sigma.star, mu, sigma){
  # if(is.positive.definite(sigma.star) == FALSE){
  #   sigma.star <- nearPD(sigma.star)$mat
  # }
  # this function calculates the densities to compare 
  term1 <- ((n-M-2)/2)*log(det(n*sigma))
  term2 <- -(1/2)*sum(diag(solve(sigma.star)%*%(n*sigma)))
  term3 <- -(1/2)*t(mu - mu.star)%*%(solve(sigma.star/n))%*%(mu - mu.star)
  return(as.vector(term1 + term2 + term3))
}
compdensity <- function(n, M, mu.star, sigma.star, mu, sigma){
  # mu = d * M, sigma = M * M * d
  # this function calculates the densities to compare 
  d <- dim(sigma)[3]
  res <- unlist(sapply(1:d, function(x) ((n-M-2)/2)*log(det(n*sigma[,,x])) -(1/2)*sum(diag(solve(sigma.star)%*%(n*sigma[,,x]))) -(1/2)*t(mu[x,] - mu.star)%*%(solve(sigma.star/n))%*%(mu[x,] - mu.star)))
  #term1 <- ((n-M-2)/2)*log(det(n*sigma))
  #term2 <- -(1/2)*sum(diag(solve(sigma.star)%*%(n*sigma)))
  #term3 <- -(1/2)*t(mu - mu.star)%*%(solve(sigma.star/n))%*%(mu - mu.star)
  return(as.vector(res))
}
compdensity.forstar <- function(n, M, mu.star, sigma.star, mu, sigma){
  # mu.star = d * M, sigma.star = M * M * d
  # mu = mu.obs, sigma = sigma.obs
  # this function calculates the densities to compare 
  d <- dim(sigma.star)[3]
  res <- unlist(sapply(1:d, function(x) ((n-M-2)/2)*log(det(n*sigma)) -(1/2)*sum(diag(solve(sigma.star[,,x])%*%(n*sigma))) -(1/2)*t(mu - mu.star[x,])%*%(solve(sigma.star[,,x]/n))%*%(mu - mu.star[x,])))
  return(as.vector(res))
}
findStar <- function(reps1, mu.obs = mu.obs, sigma.obs = sigma.obs, obj.star = density.obs.old, n = n, weight = weights, depth = depth){
  M <- nrow(sigma.obs)
  obj.old <- obj.star
  
  # strategy 1: average value as mu.old and sigma.old
  mu.old1 <- rep(mu.obs[M], M)
  sigma.old1 <- getEqualCov(sigma.obs)
  #mean.offdiag1 <- (sum(sigma.obs) - sum(diag(sigma.obs)))/(M*M - M)
  #sigma.old1 <- matrix(rep(mean.offdiag1, M*M), M, M)
  #diag(sigma.old1) <- mean(diag(sigma.obs))
  
  mu.rand.start <- genMu(mu.old1, mu.obs, reps1)
  mu.rand.start <- matrix(mu.rand.start, nrow = reps1, ncol = M, byrow = TRUE)
  sigma.rand.start <- sapply(1:reps1, function(x) genSigma(sigma.old1, sigma.obs))
  sigma.rand.start <- round(sigma.rand.start, 6)
  sigma.rand.start <- array(sigma.rand.start, dim = c(M, M, reps1))
  obj.rand.start <- compdensity.forstar(n, M, 
                                        mu.star = mu.rand.start, 
                                        sigma.star = sigma.rand.start, 
                                        mu = mu.obs, 
                                        sigma = sigma.obs)
  # select part of starting point 
  cond1 <- sapply(1:length(obj.rand.start), function(x) is.positive.definite(sigma.rand.start[,,x]))
  cond2 <- obj.rand.start > obj.old
  equalweight <- rep(1/nrow(weight), nrow(weight))
  ese.forcond3 <- ese(mu_xy = mu.rand.start, Sigma_xy = sigma.rand.start, weights = cbind(weight, equalweight))
  cond3 <- (ese.forcond3[ , 1:ncol(weight)] - ese.forcond3[ , ncol(weight)+1]) > 0
  
  mu.star.allweights <- matrix(NA, nrow = M, ncol = ncol(weight))
  sigma.star.allweights <- array(NA, dim = c(M, M, ncol(weight)))
  # Here below are for each weight
  for(k in 1:ncol(weight)){
    index <- which(cond1 == TRUE & cond2 == TRUE & cond3[ , k] == TRUE)
    if(length(index) == 0){
      crit <- 1
      while (length(index) == 0 & crit <= 10) {
        set.seed(crit*round(runif(1, 694, 2286), 0))
        mu.rand <- genMu(mu.old1, mu.obs, reps1)
        mu.rand <- matrix(mu.rand, nrow = reps1, ncol = M, byrow = TRUE)
        sigma.rand <- sapply(1:reps1, function(x) genSigma(sigma.old1, sigma.obs))
        sigma.rand <- round(sigma.rand, 6)
        sigma.rand <- array(sigma.rand, dim = c(M, M, reps1))
        obj.rand <- compdensity.forstar(n, M, 
                                        mu.star = mu.rand, 
                                        sigma.star = sigma.rand, 
                                        mu = mu.obs, 
                                        sigma = sigma.obs)
        
        # select part of starting point 
        cond1.k <- sapply(1:length(obj.rand), function(x) is.positive.definite(sigma.rand[,,x]))
        cond2.k <- obj.rand > rep(obj.old, length(obj.rand))
        ese.forcond3.k <- ese(mu_xy = mu.rand, Sigma_xy = sigma.rand, weights = cbind(weight[, k], equalweight))
        cond3.k <- (ese.forcond3.k[ , 1] - ese.forcond3.k[ , 2]) > 0
        index <- which(cond1.k == TRUE & cond2.k == TRUE & cond3.k == TRUE)
        crit <- crit + 1
      }
      if(length(index) == 0){ 
        # if we can't find an appropriate starting point -> using identity matrix as a starting point 
        mu.old <- rep(mu.obs[M], M)
        sigma.old <- getEqualCov(sigma.obs)
        obj.old <- compdensity.oneset(n, M, mu.old, sigma.old, mu.obs, sigma.obs)
      } else {
        mu.old <- mu.rand[index, ]
        sigma.old <- sigma.rand[ , , index]
        obj.old <- obj.rand[index]
      }
    } else {
      mu.old <- mu.rand.start[index, ]
      sigma.old <- sigma.rand.start[ , , index]
      obj.old <- obj.rand.start[index]
    } # output the mu.old, sigma.old, obj.old 
    
    # from the starting point to search more points 
    if(is.matrix(mu.old) == TRUE){
      # input: mu.old, mu.obs, sigma.old, sigma.obs, n, M, obj.old
      # output: updated, mu.old, sigma.old, obj.old
      for(j in 1:depth){
        # generate new random mu and sigma according to 2 rules: linear combination and randomization given a tolerant distance 
        mu.rand1    <- genMu.once(mu.old, mu.obs)
        sigma.rand1 <- sapply(1:length(index), function(x) genSigma(matrix(sigma.old[, , x], M, M), sigma.obs))
        sigma.rand1 <- round(sigma.rand1, digits = 6)
        sigma.rand1 <- array(sigma.rand1, dim = c(M, M, length(index)))
        obj.rand1   <- compdensity.forstar(n, M, 
                                           mu.star = mu.rand1, 
                                           sigma.star = sigma.rand1, 
                                           mu = mu.obs, 
                                           sigma = sigma.obs)
        # judge whether to replace 
        cond1.rand1 <- sapply(1:length(index), function(x) is.positive.definite(sigma.rand1[,,x]))
        cond2.rand1 <- obj.rand1 > obj.old
        ese.forcond3.rand1 <- ese(mu_xy = mu.rand1, Sigma_xy = sigma.rand1, weights = cbind(weight[, k], equalweight))
        cond3.rand1 <- (ese.forcond3.rand1[ , 1] - ese.forcond3.rand1[ , 2]) > 0
        
        replace.index <- which(cond1.rand1 == TRUE & cond2.rand1 == TRUE & cond3.rand1 == TRUE)
        if(length(replace.index) > 0){
          mu.old[replace.index, ] <- mu.rand1[replace.index, ]
          sigma.old[ , , replace.index] <- sigma.rand1[ , , replace.index]
          obj.old[replace.index] <- obj.rand1[replace.index]
        }
        # second rule 
        # dist1 <- max(abs(mu.old)) * 0.5
        # dist2 <- max(abs(sigma.old)) * 0.5
        # mu.rand2 <- sapply(1:length(index), function(x) genMu.random(mu.old[x, ], dist1))
        # mu.rand2 <- t(mu.rand2)
        # sigma.rand2 <- sapply(1:length(index), function(x) genSigma.random(sigma.old[ , , x], dist2))
        # sigma.rand2 <- round(sigma.rand2, digits = 6)
        # sigma.rand2 <- array(sigma.rand2, dim = c(M, M, length(index)))
        # obj.rand2 <- compdensity.forstar(n, M, 
        #                                  mu.star = mu.rand2, 
        #                                  sigma.star = sigma.rand2, 
        #                                  mu = mu.obs, 
        #                                  sigma = sigma.obs)
        
        # judge whether to replace 
        # cond1.rand2 <- sapply(1:length(index), function(x) is.positive.definite(sigma.rand2[,,x]))
        # cond2.rand2 <- obj.rand2 > obj.old
        # ese.forcond3.rand2 <- ese(mu_xy = mu.rand2, Sigma_xy = sigma.rand2, weights = cbind(weight[, k], equalweight))
        # cond3.rand2 <- (ese.forcond3.rand2[ , 1] - ese.forcond3.rand2[ , 2]) > 0
        # replace.index <- which(cond1.rand2 == TRUE & cond2.rand2 == TRUE & cond3.rand2 == TRUE)
        # if(length(replace.index) > 0){
        #   mu.old[replace.index, ] <- mu.rand2[replace.index, ]
        #   sigma.old[ , , replace.index] <- sigma.rand2[ , , replace.index]
        #   obj.old[replace.index] <- obj.rand2[replace.index]
        # }
      }
      star.index <- which.max(obj.old)
      mu.star <- mu.old[star.index, ]
      sigma.star <- sigma.old[ , , star.index]
    } else {
      # initial mu.star and sigma.star 
      mu.star <- mu.old
      sigma.star <- matrix(sigma.old, M, M)
      obj.star <- compdensity.oneset(n, M, mu.star, sigma.star, mu.obs, sigma.obs)
      for(j in 1:depth){
        mu.rand <- genMu(mu.star, mu.obs, reps1)
        mu.rand <- matrix(mu.rand, nrow = reps1, ncol = M, byrow = TRUE)
        sigma.rand <- sapply(1:reps1, function(x) genSigma(sigma.star, sigma.obs))
        sigma.rand <- round(sigma.rand, 6)
        sigma.rand <- array(sigma.rand, dim = c(M, M, reps1))
        obj.rand <- compdensity.forstar(n, M, 
                                        mu.star = mu.rand, 
                                        sigma.star = sigma.rand, 
                                        mu = mu.obs, 
                                        sigma = sigma.obs)
        # judge whether to replace 
        cond1.once <- sapply(1:length(obj.rand), function(x) is.positive.definite(sigma.rand[,,x]))
        cond2.once <- obj.rand > obj.star
        ese.forcond3.once <- ese(mu_xy = mu.rand, Sigma_xy = sigma.rand, weights = cbind(weight[ , k], equalweight))
        cond3.once <- ese.forcond3.once[ ,1] - ese.forcond3.once[ , 2] > 0
        index <- which(cond1.once == TRUE & cond2.once == TRUE & cond3.once == TRUE)
        if(length(index) == 1){
          mu.star <- mu.rand[index, ]
          sigma.star <- sigma.rand[ , , index]
          obj.star <- obj.rand[index]
        }
        if(length(index) > 1){
          mu.candidate <- mu.rand[index, ]
          sigma.candidate <- sigma.rand[ , , index]
          obj.candidate <- obj.rand[index]
          star.index <- which.max(obj.candidate)
          mu.star <- mu.candidate[star.index, ]
          sigma.star <- sigma.candidate[ , , star.index]
          obj.star <- obj.candidate[star.index]
        }  
      }
    }
    mu.star.allweights[ , k] <- mu.star
    sigma.star.allweights[ , , k] <- sigma.star
  }
  return(list(mustar = mu.star.allweights, sigmastar = sigma.star.allweights))
}
calPvalue <- function(reps2, mu.star = mu.star, sigma.star = sigma.star, density.obs = density.obs, n = n){
  
  M <- nrow(sigma.star)
  
  mu.tilde <- rmvn(reps2, mu.star, sigma.star/n)
  sigma.tilde <- rWishart(reps2, n-1, sigma.star)/n
  
  density.tilde <- compdensity(n, M, mu.star, sigma.star, mu.tilde, sigma.tilde)
  ratio <- length(which(density.tilde <= density.obs))/reps2
  
  return(ratio)
}
validation <- function(y, n, mu.true, sigma.true){
  ###### Validation: validate the algorithm by using new judgments ######
  x.train <- y + rmvn(n, mu.true, sigma.true)
  x.test <- y + rmvn(n, mu.true, sigma.true)
  
  ##### calculate weights based on x.train
  # calculate the estimated weights
  mu.obs <- colMeans(x.train - y)
  # estimated covariance matrix 
  x.bar <- colMeans(x.train - y)
  xbar.matrix <- matrix(rep(x.bar, n), ncol = n)
  xbar.matrix <- t(xbar.matrix)
  sigma.obs <- (t(x.train - y - xbar.matrix) %*% (x.train - y - xbar.matrix)) / n
  ##### step2: find the candidates of null true bias and covariance matrix #####
  optweight <- optw(mu.obs, sigma.obs)
  
  ##### work on x.test
  SA <- rowMeans(x.test)
  OW <- x.test %*% optweight
  judge <- mean((OW - y)^2) < mean((SA - y)^2)
  return(judge)
}



