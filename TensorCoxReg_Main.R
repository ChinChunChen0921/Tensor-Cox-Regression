
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
z <- as.matrix(clinical_data[,!colnames(clinical_data)%in%c('id', 'time', 'status')])
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
res <- TensorCox(time = time, status = status, X = X, z = z, n_R = 1)

res$ite
res$IC
res$df
res$b_EST
res$b_SE
res$b_PV
res$B_EST
res$B_SE
res$B_PV
