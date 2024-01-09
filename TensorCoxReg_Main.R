##############################################################################

rm(list=ls())

##---------------
## 1. Download files of R code and example data from 
## https://github.com/ChinChunChen0921/Tensor-Cox-Regression
##    and store the under "YourDirectory"
##---------------

library(survival)

codedir <- "YourDirectory/"
source(paste0(codedir,"TensorCoxReg_Function.r"))

# Main function
TensorCox <- function(DATA, n_R = n_R, opt = opt, 
                      max_ite = max_ite, tol = tol){
  n_d <- ncol(DATA$z)
  n_P <- dim(DATA$X)[1]
  n_G <- dim(DATA$X)[2]
  n <- dim(DATA$X)[3]
  df <- n_d + (n_G + n_P - n_R) * n_R
  
  Estimate <- CoxTensor_Est(DATA = DATA, n_R = n_R, opt = opt, 
                            max_ite = max_ite, tol = tol)
  Variance <- CoxTensor_Test(DATA = DATA, n_R = n_R, beta = Estimate$b_Est,
                             B1 = Estimate$B1_Est, B2 = Estimate$B2_Est)
  IC <- Calculate_IC(DATA = DATA, n_R = n_R, df = df, 
                     B = Estimate$B_Est, beta = Estimate$b_Est)
  sigular <- Variance$singular_text
  if (inherits(Variance, "try-error")) {
    df <- prod(dim(B))
    Std_B <- matrix(NA, df, df)
    Std_b <- matrix(NA, n_d, n_d)
    # print('Std_B NA')
  }else{
    V_B <- Variance$V_B
    V_b <- Variance$V_b
    Std_B <- sqrt(t(matrix(diag(Variance$V_B), n_G, n_P)))
    Std_b <- sqrt(as.matrix(diag(Variance$V_b)))
  }
  if (is.na(Std_B[1])) {
    B_PV <- Std_B
    b_PV <- Std_b
  }else{
    B_PV <- pnorm(-abs(Estimate$B_Est/Std_B)) * 2
    b_PV <- pnorm(-abs(Estimate$b_Est/Std_b)) * 2
  }
  result <- list(ite = Estimate$ite, 
                 b_EST = Estimate$b_Est, b_SD = Std_b,
                 B_EST = Estimate$B_Est, B_SD = Std_B,
                 IC = IC, df = df)
  return(result)
}

# B with T pattern

B_T <- matrix(0,6,6)
B_T[1,] <- rnorm(6, mean = 0.5, sd = 0.25)
B_T[,3] <- rnorm(6, mean = 0.5, sd = 0.25)

# Simulation data
n <- 500
p <- 0.5
cen.theta <- 1
lambda <- 1
betas <- c(0.2,0.5)
z1 <- as.matrix(rbinom(n, 1, p), nrow = n) # binary covariate
z2 <- as.matrix(rnorm(n), nrow = n) # normal covariate
z <- cbind(z1, z2)
B_shape <- B_T
n_P <- nrow(B_shape)
n_G <- ncol(B_shape)
X <- array(rnorm(n*n_P*n_G, mean = 0, sd = 1), dim=c(n_P, n_G, n))
cen.t <- runif(n, 0, cen.theta) # censoring time
surv.t <- as.numeric(- log(runif(n,0,1)) / (lambda*exp(z %*% betas+ X %hp% B_shape)))
status <- ifelse(surv.t < cen.t, 1, 0)
surv.dat <- list(n = n, z = z, surv.t = surv.t, status = status, X = X, 
                 cen.theta = cen.theta)

# Main code

opt <- 1 # 'opt' can be set to 1 or 2 to represent different optimization 
max_ite <- 200 # The maximum number of iterations
tol <- 10^-6 # Stopping criterion
# n_R is the rank for rank-R decomposition
res_nR1 <- TensorCox(DATA = surv.dat, n_R = 1, opt = opt, max_ite = max_ite, tol)
res_nR2 <- TensorCox(DATA = surv.dat, n_R = 2, opt = opt, max_ite = max_ite, tol)
res_nR3 <- TensorCox(DATA = surv.dat, n_R = 3, opt = opt, max_ite = max_ite, tol)

# Visualize the estimates of B
par(mfrow = c(1,4))
image(B_T)
image(res_nR1$B_EST)
image(res_nR2$B_EST)
image(res_nR3$B_EST)

# Check the BIC (consider the number of events in the penalty term)
res_nR1$IC$BIC; res_nR2$IC$BIC; res_nR3$IC$BIC


