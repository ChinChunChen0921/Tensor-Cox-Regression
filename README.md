# Multi-Omics Integration Analysis Using Tensor-Based Covariates in Cox Regression Modeling

## Installation and File Download

1. Install [R](https://www.r-project.org/).
2. Download the two R files (TensorCoxReg_Main.R and TensorCoxReg_Function.R) and the demo CRC data (crc_data.RData) into the same directory.

## Execution

1. In TensorCoxReg_Main.R, change your working directory to source TensorCoxReg_Function.R and load crc_data.RData.

## About the CRC Data

1. The dataset includes information on 555 subjects.
2. It focuses on progression-free survival time and status as the model's outcomes.
3. Each subject has three clinical covariates (age, gender, and stage) and second-order tensor covariates for 9 genes across 3 omic platforms.
4. The data within each platform are standardized, having a mean of 0 and a standard deviation of 1 across samples.

## TensorCox() Function Instructions

### Description

### Usage
    TensorCox(DATA, n_R = n_R, opt = opt, max_ite = max_ite, tol = tol)
### Arguments
* `DATA`: A list containing **DATA$surv.t** (survival time vector), **DATA$status** (status vector), **DATA$z** (matrix of classical covariates), and **DATA$X** (three-dimensional array of tensor covariates).

* `n_R`: A numerical constant specifying the rank for the rank-R decomposition.

* `opt`: Options for the optimization stopping criterion (**opt=1** or **opt=2**).

* `max_ite`: Maximum number of iterations for the algorithm.

* `tol`: Tolerance level for optimization.

### Value

* `ite`: Number of iterations completed by the algorithm.

* `b_EST`: Estimated coefficients for classical covariates.

* `b_SE`: Estimated standard deviation of coefficients for classical covariates.

* `B_EST`: Estimated coefficient matrix for tensor covariates.

* `B_SE`: Estimated standard deviation of coefficients for tensor covariates.

* `IC`: Information criteria, including AIC and BIC with the number of events in the penalty term.

* `df`: Degrees of freedom.
