
rm(list=ls())

library(survival)

####### Import data ###########################################################

# set your working directory 
setwd('YourDirectory/')
source("TensorCoxReg_Function.r")
clinical_data <- read.csv('clinical_data.csv')
tensor_data <- read.csv('tensor_data.csv')

####### Preprocessing #########################################################

# survival time
time <- clinical_data$time
# status
status <- clinical_data$status
# clinical covariates
z <- as.matrix(clinical_data[,!colnames(clinical_data)%in%c('time', 'status')])
# transform multi-omic data into tensor covariates
n_id <- length(unique(tensor_data$id))
n_platform <- length(unique(tensor_data$platform))
n_gene <- length(unique(tensor_data$gene))
X <- array(NA, dim = c(n_platform, n_gene, n_id))
for (i in 1:n_id){
  tmp <- tensor_data[tensor_data$id == i, 2:4]
  X[,,i] <- as.matrix(reshape(tmp, idvar = c('platform'), 
                              timevar =  'gene', 
                              direction = 'wide')[,-1])
}

####### Main Function #########################################################

# model applying rank-1 decomposition
res_nR1 <- TensorCox(time = time, status = status, X = X, z = z, n_R = 1)
# model applying rank-2 decomposition
res_nR2 <- TensorCox(time = time, status = status, X = X, z = z, n_R = 2)
# model applying rank-3 decomposition
res_nR3 <- TensorCox(time = time, status = status, X = X, z = z, n_R = 3)

###### Results ################################################################

# Visualize the estimates of B
par(mfrow = c(1,3))
image(t(res_nR1$B_EST[6:1,]))
image(t(res_nR2$B_EST[6:1,]))
image(t(res_nR3$B_EST[6:1,]))

# Check the AIC
res_nR1$IC$AIC; res_nR2$IC$AIC; res_nR3$IC$AIC

