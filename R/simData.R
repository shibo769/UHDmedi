 
#' Simulate High-Dimensional Mediation Data
#'
#' This function simulates data for estimating mediation effects under ultra-high dimensional settings.
#' It generates a design matrix X, binary treatment assignment A, mediators M, and outcome Y.
#' In addition to continuous predictors (generated from a normal distribution), the function can add
#' extra categorical predictors (with levels 1, 2, 3) and binary predictors (0/1) if the argument
#' \code{categorical = TRUE}. These extra predictors are appended to the continuous ones.
#'
#' @param nn Integer. Number of samples (default = 100).
#' @param p Integer. Number of continuous predictors in X (default = 200).
#' @param q Integer. Number of mediators (default = 200).
#' @param variance Numeric. Variance of the noise (default = 1).
#' @param k1 Integer. Sparsity level for beta coefficients (default = 3).
#' @param k2 Integer. Sparsity level for gamma coefficients (default = 3).
#' @param s_X Integer. Number of active indices in each row of the loading matrices B0 and B1 (default = 3).
#' @param categorical Logical. If TRUE, extra categorical and binary predictors are added to X (default = FALSE).
#'
#' @return A list containing:
#' \describe{
#'   \item{X}{The design matrix of dimension nn x p_total, where p_total = p + p_cat + p_bin if categorical = TRUE, otherwise p_total = p.}
#'   \item{A}{Treatment assignment vector (binary) of length nn.}
#'   \item{Y}{Outcome vector of length nn, assembled from treatment and control outcomes.}
#'   \item{M}{Mediator matrix of dimension nn x q, assembled from treatment and control mediators.}
#' }
#'
#' @examples
#' # Generate simulated data with default parameters (only continuous predictors)
#' sim_data <- simData(nn = 30, p = 50, q = 20, variance = 0.1, k1 = 3, k2 = 6, s_X = 3, categorical = TRUE)
#' str(sim_data)
#'
#' # Generate simulated data with extra categorical and binary predictors
#' sim_data2 <- simData(categorical = TRUE)
#' str(sim_data2)
#'
#' @export
simData <- function(nn = 100, p = 200, q = 200, variance = 1, k1 = 3, k2 = 3, s_X = 3, categorical = FALSE) {
  # -------------------------
  # Generate design matrix X
  # -------------------------
  # Continuous predictors from N(0,1)
  X_cont <- matrix(rnorm(nn * p), nrow = nn, ncol = p)
  
  if (categorical) {
    # Set numbers for extra categorical and binary predictors (10% each of continuous predictors)
    p_cat <- round(0.1 * p)
    p_bin <- round(0.1 * p)
    
    # Categorical predictors: generate integer levels 1,2,3 (multinomial)
    X_cat <- matrix(sample(1:3, nn * p_cat, replace = TRUE), nrow = nn, ncol = p_cat)
    
    # Binary predictors: generate from Bernoulli with probability 0.5
    X_bin <- matrix(rbinom(nn * p_bin, 1, 0.5), nrow = nn, ncol = p_bin)
    
    # Combine all predictors
    X <- cbind(X_cont, X_cat, X_bin)
  } else {
    X <- X_cont
  }
  
  # Update total number of predictors
  p_total <- ncol(X)
  
  # -------------------------
  # Generate true beta coefficients for outcomes (length = p_total) using sparsity level k1
  # -------------------------
  beta_0 <- rep(0, p_total)
  beta_0_index <- sample(1:p_total, k1, replace = FALSE)
  beta_0[beta_0_index] <- runif(k1, 1, 2)
  
  beta_1 <- rep(0, p_total)
  beta_1_index <- sample(1:p_total, k1, replace = FALSE)
  beta_1[beta_1_index] <- runif(k1, 1, 2)
  
  # -------------------------
  # Generate true gamma coefficients for mediators (length = q) using sparsity level k2
  # -------------------------
  gamma_0 <- rep(0, q)
  gamma_0_index <- sample(1:q, k2, replace = FALSE)
  gamma_0[gamma_0_index] <- runif(k2, 1.5, 2)
  
  gamma_1 <- rep(0, q)
  gamma_1_index <- sample(1:q, k2, replace = FALSE)
  gamma_1[gamma_1_index] <- runif(k2, 1.5, 2)
  
  # -------------------------
  # Generate loading matrices B0 and B1 for mediators (dimension: q x p_total)
  # -------------------------
  B0 <- matrix(0, nrow = q, ncol = p_total)
  B1 <- matrix(0, nrow = q, ncol = p_total)
  
  for (i in 1:q) {
    # For B0, randomly select s_X active indices and assign random values
    active_indices <- sample(1:p_total, s_X, replace = FALSE)
    B0[i, active_indices] <- runif(s_X)
    
    # For B1, randomly select s_X active indices and assign random values
    active_indices <- sample(1:p_total, s_X, replace = FALSE)
    B1[i, active_indices] <- runif(s_X)
  }
  
  # -------------------------
  # Generate treatment assignment A using a logistic model based on X
  # -------------------------
  alpha <- rep(0, p_total)
  alpha_index <- sample(1:p_total, k1, replace = FALSE)
  alpha[alpha_index] <- runif(k1, 0.5, 1)
  logit_vec <- 1 / (1 + exp(-X %*% alpha))
  A <- rbinom(nn, 1, logit_vec)
  
  # -------------------------
  # Identify treatment and control indices
  # -------------------------
  treat_index <- which(A == 1)
  control_index <- which(A == 0)
  nt <- length(treat_index)
  nc <- length(control_index)
  
  # -------------------------
  # Generate noise terms
  # -------------------------
  noise_Yt <- rnorm(nt, 0, variance)
  noise_Yc <- rnorm(nc, 0, variance)
  noise_M0 <- matrix(rnorm(nc * q, 0, variance), nrow = nc, ncol = q)
  noise_M1 <- matrix(rnorm(nt * q, 0, variance), nrow = nt, ncol = q)
  
  # -------------------------
  # Generate mediator matrices for control and treatment groups
  # -------------------------
  M0 <- X[control_index, , drop = FALSE] %*% t(B0) + noise_M0
  M1 <- X[treat_index, , drop = FALSE] %*% t(B1) + noise_M1
  
  # # --- Transform mediators to be between 0 and 1 using logistic transformation ---
  # M0 <- plogis(M0)
  # M1 <- plogis(M1)
  
  # -------------------------
  # Generate outcome variables for control and treatment groups
  # -------------------------
  Yc <- X[control_index, , drop = FALSE] %*% beta_0 + M0 %*% gamma_0 + noise_Yc
  Yt <- X[treat_index, , drop = FALSE] %*% beta_1 + M1 %*% gamma_1 + noise_Yt
  
  # -------------------------
  # Assemble outcome vector Y and mediator matrix M in the original order
  # -------------------------
  Y <- numeric(nn)
  M <- matrix(NA, nrow = nn, ncol = q)
  
  Y[treat_index] <- Yt
  Y[control_index] <- Yc
  
  M[treat_index, ] <- M1
  M[control_index, ] <- M0
  
  # -------------------------
  # Return only X, A, Y, and M
  # -------------------------
  cat("The dim of X" , dim(X), "\n")
  cat("The dim of M" , dim(M), "\n")
  cat("The dim of A" , dim(as.matrix(A)), "\n")
  cat("The dim of Y" , dim(as.matrix(Y)), "\n")
  cat("Dim of (X,M):", dim(X)[2] + dim(M)[2])
  
  return(list(
    X = X,
    A = A,
    Y = Y,
    M = M
  ))
}

