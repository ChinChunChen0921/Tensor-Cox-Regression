# Multi-Omics Integration Analysis Using Tensor-Based Covariates in Cox Regression Modeling

## Installation and File Download

1. Install [R](https://www.r-project.org/).
2. Download the two R files (TensorCoxReg_Main.R and TensorCoxReg_Function.R).

## Execution

1. In TensorCoxReg_Main.R, change your working directory to source TensorCoxReg_Function.R.

## TensorCox() Function Instructions


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
