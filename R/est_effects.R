 
#' Estimate Mediation Effects (ATE, NIE, NDE)
#'
#' This function estimates the total effect (TE), natural indirect effect (NIE), and natural direct effect (NDE)
#' using a debiased approach based on high-dimensional mediation analysis. It fits lasso models for the treatment
#' and control groups and then uses residual balancing for debiasing. Optionally, the outcome Y can be log-transformed.
#'
#' @param X Design matrix (n x p) for continuous predictors.
#' @param Y Outcome vector (length n).
#' @param M Mediator matrix (n x q).
#' @param A Treatment assignment vector (binary, length n).
#' @param logY Logical. If TRUE, outcome Y is log-transformed. Default is FALSE.
#' @param parallel Logical. If TRUE, parallel computing is used. Default is TRUE.
#' @param nfold1 Integer. Number of folds in cross-validation for outcome models. Default is 5.
#' @param nfold2 Integer. Number of folds in cross-validation for mediator models. Default is 3.
#' @param net_alpha Numeric. Mixing parameter for elastic net in cv.glmnet (default is 1 for lasso).
#' @param zeta Numeric. Tuning parameter for residual balancing. Default is 0.1.
#' @param optimizer Character. Optimizer method for residual balancing. Default is "pogs".
#' @param intercept Logical. If TRUE, models include an intercept. Default is FALSE.
#'
#' @return A list containing:
#' \describe{
#'   \item{TE}{Total effect estimate.}
#'   \item{NIE}{Natural indirect effect estimate.}
#'   \item{NDE}{Natural direct effect estimate.}
#'   \item{Y11}{Debiased outcome estimate for the treatment group under treatment.}
#'   \item{Y10}{Debiased outcome estimate for the control group under the treatment counterfactual.}
#'   \item{Y00}{Debiased outcome estimate for the control group.}
#' }
#'
#' @examples
#' \dontrun{
#'   effects_result <- est_effects(X = sim_data$X,
#'                                 Y = sim_data$Y,
#'                                 M = sim_data$M,
#'                                 A = sim_data$A,
#'                                  logY = FALSE,
#'                                  parallel = TRUE,
#'                                  nfold1 = 5,
#'                                  nfold2 = 5,
#'                                  net_alpha = 1,
#'                                  zeta = 0.5,
#'                                  optimizer = "pogs",
#'                                  intercept = FALSE,
#'                                  verbose = TRUE)
#'    
#'    print(effects_result)
#' }
#'
#' @export
est_effects <- function(X, Y, M, A, logY = FALSE, parallel = TRUE,
                        nfold1 = 5, nfold2 = 3, net_alpha = 1,
                        zeta = 0.1, optimizer = "pogs", intercept = FALSE, verbose = T) {
  require(doParallel)
  require(foreach)
  #suppressPackageStartupMessages(library(pogs))
  
  if (!requireNamespace("knitr", quietly = TRUE)) {
    install.packages("knitr")
  }
  library(knitr)
  
  if(logY) {
    Y <- log(Y)
  }
  
  # Remove zero variance columns from X
  if ((is.matrix(X) || is.data.frame(X)) && ncol(X) > 0) {
    var_X <- apply(X, 2, var, na.rm = TRUE)
    if(any(var_X == 0)) {
      if(verbose) message("Warning: Found predictors with zero variance in X; removing them.")
      X <- X[, var_X > 0, drop = FALSE]
      p <- ncol(X)  # update p
    }
  }
  
  p = dim(X)[2]
  q = dim(M)[2]
  
  treat_index <- which(A == 1)
  control_index <- which(A == 0)
  
  Xt <- data.matrix(X[treat_index, , drop = FALSE])
  Xc <- data.matrix(X[control_index, , drop = FALSE])
  M1 <- M[treat_index, , drop = FALSE]
  M0 <- M[control_index, , drop = FALSE]
  Yt <- Y[treat_index]
  Yc <- Y[control_index]
  
  X_bar <- colMeans(X)
  bold_Wt <- data.matrix(cbind(Xt, M1))
  bold_Wc <- data.matrix(cbind(Xc, M0))
  
  # Setup parallel backend if requested
  if(parallel) {
    num_cores <- parallel::detectCores() - 1
    cl <- parallel::makeCluster(num_cores)
    doParallel::registerDoParallel(cl)
  }
  
  # Lasso for treatment group
  fit_cv_1 <- glmnet::cv.glmnet(bold_Wt, Yt, nfolds = nfold1, intercept = intercept,
                                alpha = net_alpha, parallel = parallel)
  phi_1_hat <- as.vector(coef(fit_cv_1, s = "lambda.min"))[-1]
  beta_1_hat <- phi_1_hat[1:p]
  gamma_1_hat <- phi_1_hat[(p + 1):(p + q)]
  if(all(gamma_1_hat == 0)) {
    message("Warning: In treatment group, gamma1_hat is all zeros. Skipping this iteration.")
    return(NULL)
  }
  # Lasso for control group
  fit_cv_2 <- glmnet::cv.glmnet(bold_Wc, Yc, nfolds = nfold1, intercept = intercept,
                                alpha = net_alpha, parallel = parallel)
  phi_0_hat <- as.vector(coef(fit_cv_2, s = "lambda.min"))[-1]
  beta_0_hat <- phi_0_hat[1:p]
  gamma_0_hat <- phi_0_hat[(p + 1):(p + q)]
  if(all(gamma_0_hat == 0)) {
    message("Warning: In control group, gamma0_hat is all zeros. Skipping this bootstrap iteration.")
    return(NULL)
  }
  ### print ### 
  if (verbose == T){
    cat("-------------------------------------------------------------------------------\n")
    cat("Treatment group - Nonzero coefficients: total:", length(which(phi_1_hat != 0)),
        " beta:", length(which(beta_1_hat != 0)),
        " gamma:", length(which(gamma_1_hat != 0)), "\n")
    cat("L1 norms: beta:", sum(abs(beta_1_hat)), " gamma:", sum(abs(gamma_1_hat)), "\n")
    cat("-------------------------------------------------------------------------------\n")
    
    cat("-------------------------------------------------------------------------------\n")
    cat("Control group - Nonzero coefficients: total:", length(which(phi_0_hat != 0)),
        " beta:", length(which(beta_0_hat != 0)),
        " gamma:", length(which(gamma_0_hat != 0)), "\n")
    cat("L1 norms: beta:", sum(abs(beta_0_hat)), " gamma:", sum(abs(gamma_0_hat)), "\n")
    cat("-------------------------------------------------------------------------------\n")
  }
  ########################### Second Part of Algorithm ###########################
  # Estimate B0_hat via parallel foreach over mediators in control group
  B0_hat <- foreach::foreach(i = 1:ncol(M0), .combine = rbind, .packages = "glmnet") %dopar% {
    Mi <- M0[, i]
    cv_model <- glmnet::cv.glmnet(data.matrix(Xc), Mi, alpha = net_alpha, intercept = intercept, nfolds = nfold2)
    coef_values <- as.vector(coef(cv_model, s = "lambda.min"))[-1]
    return(coef_values)
  }
  
  contrast_Y10_1 <- c(X_bar, B0_hat %*% X_bar)
  outputY10_1 <- balanceHD::residualBalance.mean(bold_Wt, Yt,
                                                   balance.target = contrast_Y10_1,
                                                   zeta = zeta,
                                                   fit.method = "elnet",
                                                   alpha = net_alpha,
                                                   optimizer = optimizer)
                                                   # intercept = intercept)
  M0_gamma_1 <- M0 %*% gamma_1_hat
  contrast_Y10_2 <- c(X_bar)
  outputY10_2 <- balanceHD::residualBalance.mean(data.matrix(Xc),
                                                   M0_gamma_1,
                                                   balance.target = contrast_Y10_2,
                                                   zeta = zeta,
                                                   fit.method = "elnet",
                                                   alpha = net_alpha,
                                                   optimizer = optimizer)
                                                   # intercept = intercept)
  Y10_db <- outputY10_1[1] + outputY10_2[1] - t(B0_hat %*% X_bar) %*% gamma_1_hat
  
  # Estimate B1_hat via parallel foreach over mediators in treatment group
  B1_hat <- foreach::foreach(i = 1:ncol(M1), .combine = rbind, .packages = "glmnet") %dopar% {
    Mi <- M1[, i]
    cv_model <- glmnet::cv.glmnet(data.matrix(Xt), Mi, alpha = net_alpha, intercept = intercept, nfolds = nfold2)
    coef_values <- as.vector(coef(cv_model, s = "lambda.min"))[-1]
    return(coef_values)
  }
  
  contrast_Y11_1 <- c(X_bar, B1_hat %*% X_bar)
  outputY11_1 <- balanceHD::residualBalance.mean(bold_Wt, Yt,
                                                   balance.target = contrast_Y11_1,
                                                   zeta = zeta,
                                                   fit.method = "elnet",
                                                   alpha = net_alpha,
                                                   optimizer = optimizer)
                                                   # intercept = intercept)
  M1_gamma_1 <- M1 %*% gamma_1_hat
  contrast_Y11_2 <- c(X_bar)
  outputY11_2 <- balanceHD::residualBalance.mean(data.matrix(Xt),
                                                   M1_gamma_1,
                                                   balance.target = contrast_Y11_2,
                                                   zeta = zeta,
                                                   fit.method = "elnet",
                                                   alpha = net_alpha,
                                                   optimizer = optimizer)
                                                   # intercept = intercept)
  Y11_db <- outputY11_1[1] + outputY11_2[1] - t(B1_hat %*% X_bar) %*% gamma_1_hat
  
  # For control group: Y00
  contrast_Y00_1 <- c(X_bar, B0_hat %*% X_bar)
  outputY00_1 <- balanceHD::residualBalance.mean(data.matrix(bold_Wc), Yc,
                                                   balance.target = contrast_Y00_1,
                                                   zeta = zeta,
                                                   fit.method = "elnet",
                                                   alpha = net_alpha,
                                                   optimizer = optimizer)
                                                   # intercept = intercept)
  M0_gamma_0 <- M0 %*% gamma_0_hat
  contrast_Y00_2 <- c(X_bar)
  outputY00_2 <- balanceHD::residualBalance.mean(data.matrix(Xc), M0_gamma_0,
                                                   balance.target = contrast_Y00_2,
                                                   zeta = zeta,
                                                   fit.method = "elnet",
                                                   alpha = net_alpha,
                                                   optimizer = optimizer)
                                                   # intercept = intercept)
  Y00_db <- outputY00_1[1] + outputY00_2[1] - t(B0_hat %*% X_bar) %*% gamma_0_hat
  
  ########################### Results ###########################
  NIE_debiasing <- Y11_db - Y10_db
  NDE_debiasing <- Y10_db - Y00_db
  TE_de <- NIE_debiasing + NDE_debiasing
  
  if (verbose) {
    if (!requireNamespace("knitr", quietly = TRUE)) {
      install.packages("knitr")
    }
    library(knitr)
    
    # 1) Y Estimates
    Y_table <- data.frame(
      Estimate = c(Y11 = Y11_db, Y10 = Y10_db, Y00 = Y00_db)
    )
    cat("\n----------------------------------------------------------\n")
    cat("               Debiased Y(.) Estimates\n")
    cat("----------------------------------------------------------\n")
    Y_kable <- kable(
      Y_table, 
      caption = "Debiased Y(.) Estimates", 
      align = "c", 
      format = "pipe"
    )
    cat(Y_kable, sep = "\n")  
    
    cat("\n\n")  #
    
    # 2) Effect Estimates
    Effect_table <- data.frame(
      Estimate = c(NIE = NIE_debiasing, NDE = NDE_debiasing, TE = TE_de)
    )
    cat("\n----------------------------------------------------------\n")
    cat("               Effects Estimates\n")
    cat("----------------------------------------------------------\n")
    Effect_kable <- kable(
      Effect_table, 
      caption = "Effects Estimates", 
      align = "c",
      format = "pipe"
    )
    cat(Effect_kable, sep = "\n")
    cat("\n\n")  
  }
  
  
  
  # Stop parallel cluster if used
  if(parallel) {
    parallel::stopCluster(cl)
  }
  
  return(list(
    TE = TE_de,
    NIE = NIE_debiasing,
    NDE = NDE_debiasing,
    Y11 = Y11_db,
    Y10 = Y10_db,
    Y00 = Y00_db
  ))
}

