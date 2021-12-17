##### load the library #####
if (!require('mgcv')) install.packages('mgcv'); library('mgcv')
if (!require('matrixcalc')) install.packages('matrixcalc'); library('matrixcalc')
if (!require('parallel')) install.packages('parallel'); library('parallel')
if (!require('MASS')) install.packages('MASS'); library('MASS')
if (!require('foreach')) install.packages('foreach'); library('foreach')
if (!require('doParallel')) install.packages('doParallel'); library('doParallel')
if (!require('CholWishart')) install.packages('CholWishart'); library('CholWishart')
if (!require('Matrix')) install.packages('Matrix'); library('Matrix')
if (!require('quadprog')) install.packages('quadprog'); library('quadprog')

##### functions #####
## different weighting functions
# optiomal weighting 
optw <- function(x, y){
  n <- nrow(x)
  # estimated bias 
  mu.obs <- colMeans(x - y)
  # estimated covariance matrix 
  sigma.obs <- cov(x)
  
  # this function calculates the observed optimal weights 
  # by using estimated bias and covariance matrix 
  mat <- sigma.obs + mu.obs%*%t(mu.obs)
  sigma.xy <- cov(x, y)
  A <- t(rep(1, ncol(x)))
  QPmodel <- solve.QP(Dmat = nearPD(mat, doSym = TRUE)$mat, dvec = sigma.xy, Amat = t(A), bvec = c(1), meq = 1) 
  w <- QPmodel$solution
  return(w)
}
optw.sum1 <- function(x, y){
  M <- ncol(x)
  mat <- t(x) %*% x
  if(!is.positive.definite(mat)){
    mat <- nearPD(mat)$mat
  } # change the optimal weights
  Rinv <- solve(chol(mat))
  C <- cbind(rep(1, M))
  b <- c(1)
  d <- t(y) %*% x
  weights <- rep(NA, M)
  tryCatch({
    qp.model <- solve.QP(Dmat = Rinv, factorized = TRUE, dvec = d, Amat = C, bvec = b, meq = 1)
    weights <- qp.model$solution
  }, error = function(e){
    cat("No solution for OPTW.SUM1!")
  })
  return(weights)
}
optw.sum1pos <- function(x, y){
  M <- ncol(x)
  mat <- t(x) %*% x
  if(!is.positive.definite(mat)){
    mat <- nearPD(mat)$mat
  } # change the optimal weights
  Rinv <- solve(chol(mat))
  C <- cbind(rep(1, M), diag(M))
  b <- c(1, rep(0, M))
  d <- t(y) %*% x
  weights <- rep(NA, M)
  tryCatch({
    qp.model <- solve.QP(Dmat = Rinv, factorized = TRUE, dvec = d, Amat = C, bvec = b, meq = 1)
    weights <- qp.model$solution
  }, error = function(e){
    cat("No solution for OPTW.SUM1POS!")
  })
  return(weights)
}
# top five judges and equal weighting them 
topNew <- function(x, y, N){
  # this function calcualtes the individual MSE of all judges and select the top 5 judges 
  # returning the weights where top 5 persons have equal weights 
  IndMse <- colMeans((x - y)^2)
  weight <- rep(0, ncol(x))
  weight[which(rank(IndMse) <= N)] <- 1/N
  return(weight)
}
# rank performance 
RankPerformance <- function(x, y){
  # set-up
  M <- ncol(x)
  MSE.person <- colMeans((x - y)^2)
  rank.person <- rank(MSE.person, ties.method = "random") # from lowest mse to the highest mse
  # rank performance model
  MSE.RANK <- sapply(2:M, function(i) mean((rowMeans(x[ , which(rank.person <= i)]) - y)^2))
  select.num <- which.min(MSE.RANK)
  RANK.select <- which(rank.person <= select.num + 1)
  weight <- rep(0, M)
  weight[RANK.select] <- 1/length(RANK.select) 
  return(weight)
}
# CWM all models
score <- function(x, y){
  outcome <- y
  predict <- rowMeans(x)
  s <- sum((outcome - predict)^2)
  return(s)
}
contribution <- function(x, y){
  s <- score(x, y)
  sj <- c()
  for(j in 1:ncol(x)){
    sj[j] <- score(x[ , -j], y)
  }
  c <- -(s - sj)/nrow(x)
  return(c)
}
CWMall <- function(x, y){
  cont <- contribution(x, y)
  ## CEWM: equal weighting people with positive contribution
  CEWM.select <- which(cont >= 0)
  weight <- rep(0, ncol(x))
  weight[CEWM.select] <- 1/length(CEWM.select) 
  return(weight)
}
# sequencial searching 
SeqSearch <- function(x, y){
  # set-up
  M <- ncol(x)
  MSE.person <- colMeans((x - y)^2) 
  rank.person <- rank(MSE.person, ties.method = "random") # from lowest mse to the highest mse
  ## increasing way 
  SSin.select <- c(which(rank.person == 1))
  mse.train.ssin <- c(mean((x[ , SSin.select] - y)^2))
  for(i in 2:M){
    candidate <- c(1:M)[-SSin.select]
    new.mse <- sapply(1:length(candidate), function(j) mean((rowMeans(x[ , c(SSin.select , candidate[j])]) - y)^2))
    SSin.select[i] <- candidate[which.min(new.mse)]
    mse.train.ssin[i] <- min(new.mse)
  }
  SSin.select <- SSin.select[1:which.min(mse.train.ssin)]
  weight1 <- rep(0, ncol(x))
  weight1[SSin.select] <- 1/length(SSin.select)
  
  ## decreasing way 
  SSde.select <- c(1:M)
  mse.train.ssde <- c()
  remove <- c()
  for(i in 1:(M-2)){
    new.mse <- sapply(1:length(SSde.select), function(j) mean((rowMeans(x[ , SSde.select[-which(SSde.select == SSde.select[j])]]) - y)^2))
    remove[i] <- SSde.select[which.min(new.mse)]
    mse.train.ssde[i] <- mean((rowMeans(x[ , -remove]) - y)^2)
    SSde.select <- c(1:M)[-remove]
  }
  new.mse1 <- mean((x[ , SSde.select[1]] - y)^2) 
  new.mse2 <- mean((x[ , SSde.select[2]] - y)^2) 
  remove[M-1] <- SSde.select[3 - which.min(c(new.mse1, new.mse2))]
  mse.train.ssde[M-1] <- min(c(new.mse1, new.mse2))
  remove <- remove[1:which.min(mse.train.ssde)]
  SSde.select <- c(1:M)[-remove]
  weight2 <- rep(0, ncol(x))
  weight2[SSde.select] <- 1/length(SSde.select)
  
  return(list(weight.SSin = weight1, weight.SSde = weight2))
}
# regularized weighting method -> regularized covariance matrix  
#install.packages("CVTuningCov")
library(CVTuningCov)
RegCovOptw <- function(x, y){
  M <- ncol(x)
  # estimated bias 
  mu.obs <- colMeans(x - y)
  # estimated covariance matrix 
  CV.F.fit <- regular.CV(x-y, k.grid = seq(0, 100, by = 1), method = "Tapering", fold = 5, norm = "F")
  sigma.obs <- tapering(ncol(x), k = CV.F.fit$CV.k)
  # this function calculates the observed optimal weights 
  # by using estimated bias and covariance matrix 
  mat <- sigma.obs + mu.obs%*%t(mu.obs)
  sigma.xy <- cov(x, y)
  A <- t(rep(1, M))
  QPmodel <- solve.QP(Dmat = nearPD(mat, doSym = TRUE)$mat, dvec = sigma.xy, Amat = t(A), bvec = c(1), meq = 0) 
  w <- QPmodel$solution
  return(w)
}
# regularized weighting method -> LASSO 
CalGCV.lasso <- function(x, y, lambda1){
  ## gamma = 1
  M <- ncol(x)
  n <- nrow(x)
  xTx <- t(x) %*% x
  xTY <- t(x) %*% y
  ## step 1: calculate beta hat
  D <- rbind(cbind(xTx, -xTx), cbind(-xTx, xTx))
  if(!is.positive.definite(D)){
    D <- nearPD(D)$mat
  }
  A <- cbind(c(rep(1, M), rep(-1, M)), diag(2*M))
  b0 <- c(1, rep(0, 2*M))
  d <- c(xTY, -xTY) - lambda1*rep(1, 2*M)
  result <- solve.QP(D, d, A, b0, meq = 1)$solution
  weights <- result[1:M] - result[(M+1):(2*M)]
  ## step 2: p(lambda1)
  # remove zero weights
  remove.index <- which(abs(weights) <= 10^(-8))
  n0 <- length(remove.index)
  if(n0 == 0){
    invW <- diag(1/(2*abs(weights)))
    B <- t(x)%*%x + lambda1*invW
    BB <- x %*% solve(B) %*% t(x)
    plambda <- sum(diag(BB)) - n0
  }
  else if(n0 >= M-1){
    plambda <- M - length(remove.index)
  }
  else if(n0 > 0 && n0 < M-1){
    invW <- diag(1/(2*abs(weights))[-remove.index])
    B <- t(x[ , -remove.index])%*%x[ , -remove.index] + lambda1*invW
    BB <- x[ , -remove.index] %*% solve(B) %*% t(x[ , -remove.index])
    plambda <- sum(diag(BB)) - n0
  }
  ## step 3: GCV
  GCV <- (t(y-x%*%weights) %*% (y-x%*%weights))/(n * (1 - plambda/n)^2)
  return(GCV)
}
CalGCV.ridge <- function(x, y, lambda2){
  ## gamma = 2
  M <- ncol(x)
  n <- nrow(x)
  xTx <- t(x) %*% x
  xTY <- t(x) %*% y
  ## step 1: calculate beta hat
  mat <- cbind(xTx + lambda2*diag(M), rep(1, M))
  mat <- rbind(mat, t(c(rep(1, M), 0)))
  result <- solve(mat) %*% c(xTY, 1)
  weights <- result[1:M]
  ## step 2: p(lambda2)
  B <- x %*% solve(t(x)%*%x + lambda2*diag(M)) %*% t(x)
  plambda <- sum(diag(B))
  ## step 3: GCV
  GCV <- (t(y-x%*%weights) %*% (y-x%*%weights))/(n * (1 - plambda/n)^2)
  return(GCV)
}
ConstrENetZERO <- function(x, y, addl1, addl2, l1.seq, l2.seq){
  bestl1 <- NA
  bestl2 <- NA
  M <- ncol(x)
  n <- nrow(x)
  xTx <- t(x) %*% x
  xTY <- t(x) %*% y
  if (addl1 == TRUE && addl2 == FALSE){ 
    gcv <- sapply(1:length(l1.seq), function(i) CalGCV.lasso(x, y, l1.seq[i]))
    bestl1 <- l1.seq[which.min(gcv)]
    D <- rbind(cbind(xTx, -xTx), cbind(-xTx, xTx))
    if(!is.positive.definite(D)){
      D <- nearPD(D)$mat
    }
    A <- cbind(c(rep(1, M), rep(-1, M)), diag(2*M))
    b0 <- c(1, rep(0, 2*M))
    d <- c(xTY, -xTY) - bestl1*rep(1, 2*M)
    result <- solve.QP(D, d, A, b0, meq = 1)$solution
    weights <- result[1:M] - result[(M+1):(2*M)]
  }### Lasso regression
  else if (addl1 == FALSE && addl2 == TRUE){
    gcv <- sapply(1:length(l2.seq), function(i) CalGCV.ridge(x, y, l2.seq[i]))
    bestl2 <- l2.seq[which.min(gcv)]
    mat <- cbind(xTx + bestl2*diag(M), rep(1, M))
    mat <- rbind(mat, t(c(rep(1, M), 0)))
    result <- solve(mat) %*% c(xTY, 1)
    weights <- result[1:M]
  } 
  return(weights)
}



