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


# load the CRC data
load(paste0(codedir,'CRC_Data.RData'))

## The tensor covariates of 555 subject, each subject has a 3-by-13 matrix 
## represents 9 genes across 3 omic platforms.
# crc_data$X; dim(crc_data$X)

## The progression free survial (PFS) time
# crc_data$surv.t

## Status of PFS
# crc_data$status

## Three clinical covariates: age, gender, stage
# crc_data$z
TensorTest2D::

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

# The input data shoud include four items X, surv.t, status and z
# n_R is the rank for rank-R decomposition
# 'opt' can be set to 1 or 2 to represent different optimization 
# stopping criterions.


# tensor Cox regression using rank-1 decomposition
res1 <- TensorCox(DATA = crc_data, n_R = 1, opt = 1, max_ite = 1000, tol = 1e-8)
# tensor Cox regression using rank-2 decomposition
res2 <- TensorCox(DATA = crc_data, n_R = 2, opt = 1, max_ite = 1000, tol = 1e-8)
# tensor Cox regression using rank-3 decomposition
res3 <- TensorCox(DATA = crc_data, n_R = 3, opt = 1, max_ite = 1000, tol = 1e-8)

rbind(res1$IC, res2$IC, res3$IC)

# beta and B estimates
res1$b_EST; res1$B_EST

# visualize B estimates
image(res1$B_EST)
