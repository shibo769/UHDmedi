library(glmnet)

safe_solve <- function(A, b = NULL, ridge = 1e-8){
  A <- as.matrix(A)
  if(is.null(b)){
    out <- tryCatch(solve(A), error = function(e) solve(A + ridge * diag(nrow(A))))
  } else {
    out <- tryCatch(solve(A, b), error = function(e) solve(A + ridge * diag(nrow(A)), b))
  }
  return(out)
}

deSCAD <- function(z, lamb, a = 3.7){
  return(1 * (z <= lamb) + pmax((a * lamb - z), 0) / ((a - 1) * lamb) * (lamb < z))
}

LLA_h1 <- function(X, Y, M, lamb, n, p, q, n_imp = 0, S = NULL, penalize_S = TRUE){
  if(is.null(S) || length(S) == 0){
    s <- 0
    V <- X
  } else {
    s <- ncol(S)
    V <- cbind(X, S)
  }
  
  D <- cbind(M, V)
  
  w1 <- rep(0, p + q + s)
  if(p - n_imp > 0) w1[1:(p - n_imp)] <- 1
  if(s > 0 && penalize_S) w1[(p + q + 1):(p + q + s)] <- 1
  
  alpha_int <- as.vector(coef(glmnet(D, Y, family = "gaussian", alpha = 1,
                                     lambda = lamb, penalty.factor = w1,
                                     intercept = FALSE))[-1])
  
  w2 <- rep(0, p + q + s)
  if(p - n_imp > 0){
    for(j in 1:(p - n_imp)) w2[j] <- deSCAD(abs(alpha_int[j]), lamb)
  }
  if(s > 0 && penalize_S){
    for(j in (p + q + 1):(p + q + s)) w2[j] <- deSCAD(abs(alpha_int[j]), lamb)
  }
  
  alpha <- as.vector(coef(glmnet(D, Y, family = "gaussian", alpha = 1,
                                 lambda = lamb, penalty.factor = w2,
                                 intercept = FALSE))[-1])
  return(alpha)
}

HBIC_calc <- function(lamb, xx, yy, mm, S = NULL, n_imp = 0, penalize_S = TRUE){
  n <- nrow(xx)
  p <- ncol(mm)
  q <- ncol(xx)
  
  if(is.null(S) || length(S) == 0){
    s <- 0
    result <- LLA_h1(xx, yy, mm, lamb, n, p, q, n_imp = n_imp, S = NULL, penalize_S = penalize_S)
    alpha0 <- result[1:p]
    alpha1 <- result[(p + 1):(p + q)]
    alpha2 <- NULL
    tmp <- yy - mm %*% alpha0 - xx %*% alpha1
    df <- sum(alpha0 != 0) + q
  } else {
    s <- ncol(S)
    result <- LLA_h1(xx, yy, mm, lamb, n, p, q, n_imp = n_imp, S = S, penalize_S = penalize_S)
    alpha0 <- result[1:p]
    alpha1 <- result[(p + 1):(p + q)]
    alpha2 <- result[(p + q + 1):(p + q + s)]
    tmp <- yy - mm %*% alpha0 - xx %*% alpha1 - S %*% alpha2
    df <- sum(alpha0 != 0) + q + ifelse(penalize_S, sum(alpha2 != 0), s)
  }
  
  sigma_hat <- as.numeric(t(tmp) %*% tmp / n)
  BIC <- log(sigma_hat + 1e-12) + df * log(log(n)) * log(p + q + s) / n
  
  return(list(BIC = BIC, alpha0 = alpha0, alpha1 = alpha1,
              alpha2 = alpha2, sigma1_hat = sigma_hat, df = df))
}

mediationInference <- function(X, Y, M = NULL, S = NULL){
  X <- as.matrix(X)
  Y <- as.vector(Y)
  n <- nrow(X)
  q <- ncol(X)
  
  if(is.null(M) || length(M) == 0) M <- matrix(numeric(0), nrow = n, ncol = 0)
  if(is.null(S) || length(S) == 0) S <- matrix(numeric(0), nrow = n, ncol = 0)
  
  M <- as.matrix(M)
  S <- as.matrix(S)
  p <- ncol(M)
  s <- ncol(S)
  
  Z <- cbind(M, X, S)
  MS <- cbind(M, S)
  V <- cbind(X, S)
  
  if(ncol(MS) == 0){
    RSS02 <- as.numeric(t(Y) %*% Y)
  } else {
    alpha0_tld <- safe_solve(t(MS) %*% MS, t(MS) %*% Y)
    RSS02 <- as.numeric(t(Y - MS %*% alpha0_tld) %*% (Y - MS %*% alpha0_tld))
  }
  
  alpha_rf <- safe_solve(t(Z) %*% Z, t(Z) %*% Y)
  
  alpha0_hat <- if(p > 0) alpha_rf[1:p] else numeric(0)
  alpha1_hat <- alpha_rf[(p + 1):(p + q)]
  alpha2_hat <- if(s > 0) alpha_rf[(p + q + 1):(p + q + s)] else numeric(0)
  
  res <- Y - Z %*% alpha_rf
  RSS12 <- as.numeric(t(res) %*% res)
  df <- p + q + s
  sigma1_hat <- RSS12 / max(n - df, 1)
  
  gamma_hat <- safe_solve(t(V) %*% V, t(V) %*% Y)
  sigmaT_hat <- as.numeric(t(Y - V %*% gamma_hat) %*% (Y - V %*% gamma_hat) / max(n - q - s, 1))
  beta_hat <- gamma_hat[1:q] - alpha1_hat
  sigma2_hat <- pmax(0, sigmaT_hat - sigma1_hat)
  
  Sigma_VV <- t(V) %*% V / n
  invVV <- safe_solve(Sigma_VV)
  
  if(p > 0){
    Sigma_MM <- t(M) %*% M / n
    Sigma_MV <- t(M) %*% V / n
    Sigma_VM <- t(Sigma_MV)
    Schur <- Sigma_MM - Sigma_MV %*% invVV %*% Sigma_VM
    Bmat <- invVV %*% Sigma_VM %*% safe_solve(Schur) %*% Sigma_MV %*% invVV
  } else {
    Bmat <- matrix(0, nrow = q + s, ncol = q + s)
  }
  
  var_alpha1_hat <- sigma1_hat * (invVV + Bmat)[1:q, 1:q, drop = FALSE]
  cov_beta_hat <- sigma2_hat * invVV + sigma1_hat * Bmat
  cov_beta_hat <- cov_beta_hat[1:q, 1:q, drop = FALSE]
  
  Sn <- as.numeric(n * t(beta_hat) %*% safe_solve(cov_beta_hat) %*% beta_hat)
  Tn <- as.numeric((n - q - s) * (RSS02 - RSS12) / RSS12)
  
  p_beta <- 1 - pchisq(Sn, q)
  p_alpha1 <- 1 - pchisq(Tn, q)
  
  std_alpha1 <- sqrt(diag(var_alpha1_hat / n))
  std_beta <- sqrt(diag(cov_beta_hat / n))
  
  result.df <- data.frame(
    Coefficient = c("alpha1", "beta"),
    Estimate = c(as.numeric(alpha1_hat), as.numeric(beta_hat)),
    SE = c(std_alpha1, std_beta),
    Test_statistics = c(Tn, Sn),
    p_value = c(p_alpha1, p_beta)
  )
  
  return(list(Sn = Sn, Tn = Tn, sigma1_hat = sigma1_hat,
              beta_hat = beta_hat, alpha0_hat = alpha0_hat,
              alpha1_hat = alpha1_hat, alpha2_hat = alpha2_hat,
              B = Bmat, var_beta = cov_beta_hat,
              var_alpha1_hat = var_alpha1_hat,
              p_beta = p_beta, p_alpha1 = p_alpha1,
              summary_result = result.df))
}