Sur_Coef <- function (X, surv.t, status, offset_vec) {
  coefficients(coxph(Surv(surv.t, status) ~ X + offset(offset_vec)))
}
CoxTensor_Est <- function(DATA, n_R, opt, max_ite, tol){
  
  `%w%` <- function(X, B) sapply(1:dim(X)[3], function(i) X[, , i] %*% B)
  `%wt%` <- function(X, B) sapply(1:dim(X)[3], function(i) t(X[, , i]) %*% B)
  `%hp%` <- function (X, B) sapply(1:dim(X)[3], function(i) sum(X[, , i] * B))
  
  X <- DATA$X; y <- DATA$surv.t; status <- DATA$status; W <- as.matrix(DATA$z)
  n_vec <- dim(X); n <- n_vec[3]; n_P <- n_vec[1]; n_G <- n_vec[2]
  n_d <- ncol(W)
  
  Z <- matrix(0, n, n_P * n_G)
  for (i in 1:n) Z[i, ] <- c(X[, , i])
  Q <- cbind(W, Z)
  B_1_temp <- matrix(rnorm(n_P * n_R), n_P, n_R)
  offset_vec <- rep(0, n)
  # initial: beta
  beta <- Sur_Coef(W, y, status, offset_vec)
  # initial: B1
  B_1 <- matrix(B_1_temp[, 1:n_R], n_P, n_R)
  H <- matrix(B_1[1:n_R, 1:n_R], n_R, n_R)
  C <- matrix(B_1[-c(1:n_R), 1:n_R], n_P - n_R, n_R)
  # initial: B2
  B_2 <- matrix(0, n_G, n_R)
  B <- B_1 %*% t(B_2)
  
  test <- 1
  for (ite_index in 1:max_ite) {
    offset_vec <- W %*% beta
    # Update B2
    B_2_new <- matrix(Sur_Coef(t(X %wt% B_1), y, status,
                               offset_vec), n_G, n_R)
    # Update B1
    B_1_new <- matrix(Sur_Coef(t(X %w% B_2_new), y, status, 
                               offset_vec), n_P, n_R)
    G <- matrix(B_1_new[n_R, n_R], n_R, n_R)
    B_new <- B_1_new %*% t(B_2_new)
    offset_vec <- X %hp% B_new
    # Update beta
    beta_new <- Sur_Coef(W, y, status, offset_vec)
    # two options for 
    if (opt == 1) {
      test <- max(max(abs(B_new - B)), max(abs(beta_new - beta)))
    }
    if (opt == 2) {
      test <- max(c(abs(B_1 - B_1_new)), c(abs(B_2 - B_2_new)), abs(beta_new - beta))
    }
    B_1 <- B_1_new
    B_2 <- B_2_new
    beta <- beta_new
    B <- B_new
    if (test < tol) 
      break
  }
  return(list(b_Est = beta, B_Est = B, 
              B1_Est = B_1, B2_Est = B_2, 
              ite = ite_index))
}

CoxTensor_Test <- function(DATA, n_R, B1, B2, beta) {
  
  `%hp%` <- function (X, B) sapply(1:dim(X)[3], function(i) sum(X[, , i] * B))
  # box product
  `%b%` <- function(A, B) {
    n_A <- ncol(A)
    n_B <- ncol(B)
    O <- matrix(0, nrow(A) * nrow(B), n_A * n_B)
    for (i in 1:n_B) for (j in 1:n_A) O[, (i - 1) *
                                          n_A + j] <- A[, j] %x% B[, i]
    return(O)
  }
  surv.t <- DATA$surv.t
  y <- DATA$surv.t[order(surv.t)]
  X <- DATA$X[,,order(surv.t)]
  status <- DATA$status[order(surv.t)]
  Z <- DATA$z[order(surv.t),]
  n_vec <- dim(X)
  n <- n_vec[3]
  n_P <- n_vec[1]
  n_G <- n_vec[2]
  n_d <- ncol(Z)
  
  B <- B1 %*% t(B2)
  B_12 <- matrix(B1[-c(1:n_R), 1:n_R], n_P - n_R, n_R)
  C <- matrix(B1[1:n_R, 1:n_R], n_R, n_R)
  df <- (n_P - n_R + n_G) * n_R + n_d
  vec_c <- c(C)
  b12 <- c(B_12)
  b2 <- c(B2)
  b1 <- c(B1)
  w <- exp(X %hp% B + Z %*% beta)
  X1 <- lapply(1:dim(X)[3], function(i) matrix(X[1:n_R, 1:n_G, i], ncol = n_G))
  X2 <- lapply(1:dim(X)[3], function(i) matrix(X[(n_R+1):n_P, 1:n_G, i], ncol = n_G))
  K1 <- lapply(X1, function(x)  matrix(diag(n_R) %x% x, ncol = n_R*n_G))
  K2 <- lapply(X2, function(x) matrix(diag(n_R) %x% x, ncol = n_R*n_G))
  W1 <- lapply(K2, function(k2) k2%*%b2)
  W2 <- mapply(function(k1,k2) t(k1)%*%vec_c + t(k2)%*%b12, k1 = K1, k2 = K2, SIMPLIFY = F)
  
  Hessian_beta <- matrix(0, n_d, n_d)
  Hessian_B12 <- matrix(0, n_R*(n_P-n_R), n_R*(n_P-n_R))
  Hessian_B2 <- matrix(0, n_R*n_G, n_R*n_G)
  Hessian_beta_B12 <- matrix(0, n_d, n_R*(n_P-n_R))
  Hessian_beta_B2 <- matrix(0, n_d, n_R*n_G)
  Hessian_B12_B2 <- matrix(0, n_R*(n_P-n_R), n_R*n_G)
  for (i in 1:n){
    w_risk <- w[i:n]
    wij <- w_risk/sum(w_risk)
    V1 <- matrix(Z[i:n, ], ncol = n_d)
    V2 <- W1[i:n]
    V3 <- W2[i:n]
    Zwi <- t(wij) %*% V1
    Awi <- Reduce('+',mapply('*', wij, V2, SIMPLIFY = F))
    Bwi <- Reduce('+',mapply('*', wij, V3, SIMPLIFY = F))
    # diagonal
    beta_H_i <- Reduce('+',lapply(1:(n-i+1), function(j) wij[j]*t(V1[j, ] - Zwi) %*% (V1[j, ] - Zwi)))
    B12_H_i <- Reduce('+',lapply(1:(n-i+1), function(j) wij[j]*(V2[[j]] - Awi) %*% t(V2[[j]] - Awi)))
    B2_H_i <- Reduce('+',lapply(1:(n-i+1), function(j) wij[j]*(V3[[j]] - Bwi) %*% t(V3[[j]] - Bwi)))
    # off-diagonal
    beta_B12_H_i <- Reduce('+',lapply(1:(n-i+1), function(j) wij[j]*t(V1[j, ] - Zwi) %*% t(V2[[j]] - Awi)))
    beta_B2_H_i <- Reduce('+',lapply(1:(n-i+1), function(j) wij[j]*t(V1[j, ] - Zwi) %*% t(V3[[j]] - Bwi)))
    B12_B2_H_i <- Reduce('+',lapply(1:(n-i+1), function(j) wij[j]*(V2[[j]] - Awi) %*% t(V3[[j]] - Bwi)))
    
    Hessian_beta <- Hessian_beta-status[i]*beta_H_i
    Hessian_B12 <- Hessian_B12-status[i]*B12_H_i
    Hessian_B2 <- Hessian_B2-status[i]*B2_H_i
    Hessian_beta_B12 <- Hessian_beta_B12-status[i]*beta_B12_H_i
    Hessian_beta_B2 <- Hessian_beta_B2-status[i]*beta_B2_H_i
    Hessian_B12_B2 <- Hessian_B12_B2-status[i]*B12_B2_H_i
  }
  R1 <- cbind(Hessian_beta, Hessian_beta_B12, Hessian_beta_B2)
  R2 <- cbind(t(Hessian_beta_B12), Hessian_B12, Hessian_B12_B2)
  R3 <- cbind(t(Hessian_beta_B2), t(Hessian_B12_B2), Hessian_B2)
  Hessian <- rbind(R1, R2, R3)
  I <- -Hessian
  Iinv <- try(solve(I), T)
  if (inherits(Iinv,"try-error")){
    Iinv <- Ginv(I)
    singular_text <- 'Ginv'
  }else{
    singular_text <- 'solve'
  }
  V_beta <- as.matrix(Iinv[1:n_d, 1:n_d])
  Iinv <- Iinv[-c(1:n_d), -c(1:n_d)]
  if (!is.null(attr(Iinv, "class"))) {
    df <- prod(dim(B))
    V_B <- matrix(NA, df, df)
  }else{
    A1 <- cbind(matrix(0, n_G * n_R, (n_P - n_R) * n_R),
                diag(n_G * n_R))
    A2 <- cbind(diag(n_P - n_R) %b% B2, B_12 %x% diag(n_G))
    A3 <- B_12 %x% diag(n_G)
    A4 <- C %x% diag(n_G)
    Vff <- A1 %*% Iinv %*% t(A1)
    V11 <- A4 %*% Vff %*% t(A4)
    V12 <- A4 %*% Vff %*% t(A3)
    V22 <- A2 %*% Iinv %*% t(A2)
    V_B <- rbind(cbind(V11, V12), cbind(t(V12), V22))
  }
  return(list(V_B = V_B, V_b = V_beta, singular_text = singular_text))
}
Calculate_IC <- function(DATA, n_R, B, beta, df){
  `%hp%` <- function (X, B) sapply(1:dim(X)[3], function(i) sum(X[, , i] * B))
  surv.t <- DATA$surv.t
  X <- DATA$X[,,order(surv.t)]
  status <- DATA$status[order(surv.t)]
  W <- DATA$z[order(surv.t),]
  n_vec <- dim(X)
  n <- n_vec[3]
  S <- X %hp% B + W %*% beta
  exp_S <- exp(S)
  mlpl <- 0
  for (i in 1:n){
    mlpl <- mlpl + status[i]*(S[i]-log(sum(exp_S[i:n])))
  }
  AIC <- -2*mlpl+2*df
  # BIC use the number of events in the penalty term
  BIC_NoE<- -2*mlpl+log(sum(status))*df
  return(list(AIC = AIC, BIC_NoE = BIC_NoE))
}
