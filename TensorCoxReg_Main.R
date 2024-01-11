##############################################################################

rm(list=ls())

##---------------
## Download files of R code from 
## https://github.com/ChinChunChen0921/Tensor-Cox-Regression
## and store the under "YourDirectory"
##---------------

set.seed(123)
library(survival)

# codedir <- "YourDirectory/"
source(paste0(codedir,"TensorCoxReg_Function.r"))
`%hp%` <- function (X, B) sapply(1:dim(X)[3], function(i) sum(X[, , i] * B))

# B with T pattern
B_T <- matrix(0,6,6)
B_T[1,] <- rnorm(6, mean = 0.5, sd = 0.25)
B_T[,3] <- rnorm(6, mean = 0.5, sd = 0.25)

# Simulation data
n <- 500 # sample size
cen.theta <- 1 # theta to control the censoring rate
lambda <- 1 # baseline hazard
p <- 0.5 # parameter for binary covariate
z1 <- as.matrix(rbinom(n, 1, p), nrow = n) # binary covariate
z2 <- as.matrix(rnorm(n), nrow = n) # normal covariate
z <- cbind(z1, z2)
betas <- c(0.2,0.5) # effects related to z
B_shape <- B_T
n_P <- nrow(B_shape)
n_G <- ncol(B_shape)
# tensor covariate
X <- array(rnorm(n*n_P*n_G, mean = 0, sd = 1), dim=c(n_P, n_G, n)) 
# censoring time
cen.t <- runif(n, 0, cen.theta) 
# actual survival time
surv.t <- as.numeric(-log(runif(n,0,1))/(lambda*exp(z%*%betas + X%hp%B_shape))) 
# status of the occurrence of event
status <- ifelse(surv.t < cen.t, 1, 0)
surv.dat <- list(n = n, z = z, surv.t = surv.t, status = status, X = X, 
                 cen.theta = cen.theta)

# Main code
max_ite <- 200 # The maximum number of iterations
tol <- 10^-8 # Stopping criterions.
# n_R is the rank for rank-R decomposition
res_nR1 <- TensorCox(DATA = surv.dat, n_R = 1, max_ite = max_ite, tol = tol)
res_nR2 <- TensorCox(DATA = surv.dat, n_R = 2, max_ite = max_ite, tol = tol)
res_nR3 <- TensorCox(DATA = surv.dat, n_R = 3, max_ite = max_ite, tol = tol)

# Visualize the estimates of B
par(mfrow = c(1,4))
image(B_T)
image(res_nR1$B_EST)
image(res_nR2$B_EST)
image(res_nR3$B_EST)

# Check the AIC
res_nR1$IC$AIC; res_nR2$IC$AIC; res_nR3$IC$AIC

