##### load the library #####
if (!require('mgcv')) install.packages('mgcv'); library('mgcv')
if (!require('matrixcalc')) install.packages('matrixcalc'); library('matrixcalc')
if (!require('parallel')) install.packages('parallel'); library('parallel')
if (!require('MASS')) install.packages('MASS'); library('MASS')
if (!require('foreach')) install.packages('foreach'); library('foreach')
if (!require('doParallel')) install.packages('doParallel'); library('doParallel')
if (!require('CholWishart')) install.packages('CholWishart'); library('CholWishart')

## cross validation functions 
CrossValidation <- function(x, y, k.fold, round.num){
  # k.fold is a number between 0 to 1, representing the percent of training data 
  # round.num is the total number of rounds for cross validation 
  outMSE <- c()
  for(i in 1:round.num){
    # prepare the training data and testing data
    train.index <- sort(sample(1:nrow(x), floor(nrow(x)*k.fold)))
    test.index <- c(1:nrow(x))[-train.index]
    train.x <- x[train.index, ]
    train.y <- y[train.index]
    test.x <- x[test.index, ]
    test.y <- y[test.index]
    # obtain all weights based on test.index
    optweight <- optw(train.x, train.y) # 1
    top3equalweight <- topNew(train.x, train.y, N = 3) # 2
    rankweight <- RankPerformance(train.x, train.y) # 3
    cwmweight <- CWMall(train.x, train.y) #4
    SSweights <- SeqSearch(train.x, train.y)
    ssinweight <- SSweights$weight.SSin #5
    ssdeweight <- SSweights$weight.SSde #6
    regcovweight <- RegCovOptw(train.x, train.y) #7
    reglassoweight <- ConstrENetZERO(train.x, train.y, addl1 = TRUE, addl2 = FALSE, l1.seq = 2^(seq(0, 10, by = 0.2)), l2.seq = NULL) #8
    equalweight <- 1/ncol(train.x)
    weights <- cbind(optweight, top3equalweight, rankweight, cwmweight, 
                     ssinweight, ssdeweight, regcovweight, reglassoweight, 
                     equalweight)
    # calculate mse
    mse.test <- colMeans((test.x %*% weights - test.y)^2)
    outMSE <- rbind(outMSE, mse.test)
  }
  results <- colMeans(outMSE)
  return(results)  
}

CrossValidation_subset <- function(x, y, k.fold, round.num){
  # k.fold is a number between 0 to 1, representing the percent of training data 
  # round.num is the total number of rounds for cross validation 
  outMSE <- c()
  for(i in 1:round.num){
    # prepare the training data and testing data
    train.index <- sort(sample(1:nrow(x), floor(nrow(x)*k.fold)))
    test.index <- c(1:nrow(x))[-train.index]
    train.x <- x[train.index, ]
    train.y <- y[train.index]
    test.x <- x[test.index, ]
    test.y <- y[test.index]
    # obtain all weights based on test.index
    optweight1 <- optw.sum1(train.x, train.y) # 1
    optweight2 <- optw.sum1pos(train.x, train.y) # 1
    equalweight <- 1/ncol(train.x)
    weights <- cbind(optweight1, optweight2, equalweight)
    # calculate mse
    mse.test <- colMeans((test.x %*% weights - test.y)^2)
    outMSE <- rbind(outMSE, mse.test)
  }
  results <- colMeans(outMSE)
  return(results)  
}
