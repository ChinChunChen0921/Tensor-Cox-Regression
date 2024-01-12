# Multi-Omics Integration Analysis Using Tensor-Based Covariates in Cox Regression Modeling

## Installation and File Download

1. Install [R](https://www.r-project.org/).
2. Download the two R files (TensorCoxReg_Main.R and TensorCoxReg_Function.R).

## Execution

1. In TensorCoxReg_Main.R, change your working directory to source TensorCoxReg_Function.R.
2. In TensorCoxReg_Main.R, we provide a simplec code using simulated data for illustration purposes.
 
## TensorCox() Function Instructions


### Usage
    TensorCox(DATA, n_R = n_R, max_ite = max_ite, tol = tol)
### Arguments
* `DATA`: A list containing **DATA$surv.t** (survival time vector), **DATA$status** (status vector), **DATA$z** (matrix of classical covariates), and **DATA$X** (three-dimensional array of tensor covariates).

* `n_R`: A numerical constant specifying the rank for the rank-R decomposition.

* `max_ite`: Maximum number of iterations for the algorithm.

* `tol`: Tolerance level for optimization.

### Value

* `ite`: The number of iterations completed by the algorithm.

* `b_EST`: The estimated coefficients for classical covariates.

* `b_SE`: The estimated standard deviation of coefficients for classical covariates.

* `b_PV`: The p-value for classical covariates.

* `B_EST`: The estimated coefficient matrix for tensor covariates.

* `B_SE`: The estimated standard deviation of coefficients for tensor covariates.

* `B_PV`: The p-value for tensor covariates.

* `IC`: Information criteria, including AIC and BIC with the number of events in the penalty term.

* `df`: Degrees of freedom.
