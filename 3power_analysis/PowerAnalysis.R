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
reps2.slow.power <- 20000

infdata <- read.csv("inf_filtered.csv", header = T, sep = ",") # change empirical data
infdata <- infdata[ , -1]
M <- 9 # this is for inflation rate data
X <- as.matrix(infdata[ , 4:12]) # this is for inflation rate data
Y <- infdata[ , 3]

# unempdata <- read.csv("unemp_filtered.csv", header = T, sep = ",") # change empirical data
# unempdata <- unempdata[ , -1]
# M <- 8 # this is for inflation rate data
# X <- as.matrix(unempdata[ , 4:11]) # this is for inflation rate data
# Y <- unempdata[ , 3]

## Step 1: estimate the empirical mean and covariance matrix 
XY_all <- cbind(X, Y)
mu_xy.emp <- colMeans(XY_all)
Sigma_xy.emp <- cov(XY_all)

## Step 2: find the Sigma^* and mu^* for each weighting method by using the whole set of data 
mu_xy.obs <- colMeans(XY_all)
Sigma_xy.obs <- cov(XY_all)
# get all 8 different weights
optweight <- optw.sum1(X, Y) #1
top3equalweight <- topNew(X, Y, N = 3) #2
rankweight <- RankPerformance(X, Y) #3
cwmweight <- CWMall(X, Y) #4
SSweights <- SeqSearch(X, Y)
ssinweight <- SSweights$weight.SSin #5
ssdeweight <- SSweights$weight.SSde #6
#regcovweight <- RegCovOptw(X, Y) #7
reglassoweight <- ConstrENetZERO(X, Y, addl1 = TRUE, addl2 = FALSE, l1.seq = 2^(seq(0, 10, by = 0.2)), l2.seq = NULL) #8
weights <- cbind(optweight, top3equalweight, rankweight, cwmweight, 
                 ssinweight, ssdeweight, reglassoweight)
mu_xy.star <- rep(mu_xy.obs[M+1], M+1)
Sigma_xy.star <- getEqualCov(Sigma_xy.obs)
density.obs.old <- compdensity.oneset(nrow(X), M+1, mu_xy.star, Sigma_xy.star, mu_xy.obs, Sigma_xy.obs)
set.seed(23467)
solvestar <- findStar(reps1 = reps1, mu.obs = mu_xy.obs, sigma.obs = Sigma_xy.obs, obj.star = density.obs.old, n = nrow(X), 
                      weight = weights, depth = depth)
mu_xy.star <- solvestar$mustar
Sigma_xy.star <- solvestar$sigmastar

## Step 3: Generate Sigma_hat and mu_hat to calculate the p-value 
n.vec <- c(seq(66, 96, by = 2), seq(100, 128, by = 4))
draws <- 100
for(i in 1:length(n.vec)){
  n <- n.vec[i]
  for(j in 1:draws){
    set.seed(j*i*234)
    # random bootstrap of sets of judgments 
    xy_all <- rmvn(n, mu_xy.emp, Sigma_xy.emp)
    # estimated parameters: mu_hat and Sigma_hat 
    mu_xy.obs <- colMeans(xy_all)
    Sigma_xy.obs <- cov(xy_all)
    # calculate p-value for each weighting method
    density.obs <- compdensity.forstar(n, M+1, t(mu_xy.star), Sigma_xy.star, mu_xy.obs, Sigma_xy.obs)
    pv <- sapply(1:ncol(mu_xy.star), function(x) calPvalue(reps2 = reps2.slow.power, 
                                                           mu.star = mu_xy.star[ , x], 
                                                           sigma.star = Sigma_xy.star[ , , x],
                                                           density.obs = density.obs[x], n = n))
    write.table(t(c(n, j, pv)), "power_inf.csv", sep = ",", append = TRUE, row.names = FALSE, col.names = FALSE)
    #write.table(t(c(n, j, pv)), "power_unemp.csv", sep = ",", append = TRUE, row.names = FALSE, col.names = FALSE)
    print(paste(n, "-", j, sep = ""))
  }
}
