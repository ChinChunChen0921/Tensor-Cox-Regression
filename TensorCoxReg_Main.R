CoxTensor_Main <- function(DATA, n_R = n_R, opt = opt, max_ite = max_ite, tol = tol){
  n_d <- ncol(DATA$z)
  n_P <- dim(DATA$X)[1]
  n_G <- dim(DATA$X)[2]
  n <- dim(DATA$X)[3]
  df <- n_d + (n_G + n_P - n_R) * n_R
  
  Estimate <- CoxTensor_Est(DATA = DATA, n_R = n_R, opt = opt, max_ite = max_ite, tol = tol)
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
  result <- list(ite = Estimate$ite, b_EST = Estimate$b_Est, b_SD = Std_b,
                 b_PV = b_PV, B_EST = Estimate$B_Est, c_index = Estimate$c_index,
                 B_SD = Std_B, B_PV = B_PV,
                 IC = IC, df = df, sigular = sigular)
  return(result)
}