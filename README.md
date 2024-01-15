# Multi-Omics Integration Analysis Using Tensor-Based Covariates in Cox Regression Modeling

## Installation and File Download

1. Install [R](https://www.r-project.org/).
2. Download the two R files (TensorCoxReg_Main.R and TensorCoxReg_Function.R).

## Execution

1. In TensorCoxReg_Main.R, change your working directory to source TensorCoxReg_Function.R.
2. In TensorCoxReg_Main.R, we provide a simplec code using simulated data for illustration purposes.
 
## TensorCox() Function Instructions


### Usage
    TensorCox(time, status, X, z, n_R, max_ite = 200, tol = 10^-6)
### Arguments
* `time`: survival time vector.

* `status`: status vector.

* `X`: three-dimensional array of tensor covariates.

* `z`: matrix of classical covariates.

* `n_R`: A numerical constant specifying the rank for the rank-R decomposition.

* `max_ite`: The maximum number of iterations for the algorithm. The default is set to 200 iterations.

* `tol`: Tolerance level for the difference in partial likelihood values between two iterations. The default value is set to $10^{-6}$.

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
