 
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
                        zeta = 0.1, optimizer = "pogs.new", intercept = FALSE, verbose = T, K = 2.75) {
  require(doParallel)
  require(foreach)
  # require(pogs)
  # source(utilis)
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
  outputY10_1 <- residualBalance.mean(bold_Wt, Yt,
                                                   balance.target = contrast_Y10_1,
                                                   zeta = zeta,
                                                   fit.method = "elnet",
                                                   alpha = net_alpha,
                                                   optimizer = optimizer,
                                                   K = 2.75)
  M0_gamma_1 <- M0 %*% gamma_1_hat
  contrast_Y10_2 <- c(X_bar)
  outputY10_2 <- residualBalance.mean(data.matrix(Xc),
                                                   M0_gamma_1,
                                                   balance.target = contrast_Y10_2,
                                                   zeta = zeta,
                                                   fit.method = "elnet",
                                                   alpha = net_alpha,
                                                   optimizer = optimizer,
                                                   K = 2.75)
  Y10_db <- outputY10_1[1] + outputY10_2[1] - t(B0_hat %*% X_bar) %*% gamma_1_hat
  
  # Estimate B1_hat via parallel foreach over mediators in treatment group
  B1_hat <- foreach::foreach(i = 1:ncol(M1), .combine = rbind, .packages = "glmnet") %dopar% {
    Mi <- M1[, i]
    cv_model <- glmnet::cv.glmnet(data.matrix(Xt), Mi, alpha = net_alpha, intercept = intercept, nfolds = nfold2)
    coef_values <- as.vector(coef(cv_model, s = "lambda.min"))[-1]
    return(coef_values)
  }
  
  contrast_Y11_1 <- c(X_bar, B1_hat %*% X_bar)
  outputY11_1 <- residualBalance.mean(bold_Wt, Yt,
                                                   balance.target = contrast_Y11_1,
                                                   zeta = zeta,
                                                   fit.method = "elnet",
                                                   alpha = net_alpha,
                                                   optimizer = optimizer,
                                                   K = 2.75)
  M1_gamma_1 <- M1 %*% gamma_1_hat
  contrast_Y11_2 <- c(X_bar)
  outputY11_2 <- residualBalance.mean(data.matrix(Xt),
                                                   M1_gamma_1,
                                                   balance.target = contrast_Y11_2,
                                                   zeta = zeta,
                                                   fit.method = "elnet",
                                                   alpha = net_alpha,
                                                   optimizer = optimizer,
                                                   K = 2.75)
  Y11_db <- outputY11_1[1] + outputY11_2[1] - t(B1_hat %*% X_bar) %*% gamma_1_hat
  
  # For control group: Y00
  contrast_Y00_1 <- c(X_bar, B0_hat %*% X_bar)
  outputY00_1 <- residualBalance.mean(data.matrix(bold_Wc), Yc,
                                                   balance.target = contrast_Y00_1,
                                                   zeta = zeta,
                                                   fit.method = "elnet",
                                                   alpha = net_alpha,
                                                   optimizer = optimizer,
                                                   K = 2.75)
  M0_gamma_0 <- M0 %*% gamma_0_hat
  contrast_Y00_2 <- c(X_bar)
  outputY00_2 <- residualBalance.mean(data.matrix(Xc), M0_gamma_0,
                                                   balance.target = contrast_Y00_2,
                                                   zeta = zeta,
                                                   fit.method = "elnet",
                                                   alpha = net_alpha,
                                                   optimizer = optimizer,
                                                   K = 2.75)
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


############# Utilis Functions ##############
#' Estimate mean outcome at balance.target via residual balancing
#'
#' @param XW the input features for the sub-population of interest
#' @param YW the observed responses for the sub-population of interest
#' @param balance.target the desired center of the dataset
#' @param allow.negative.weights whether negative gammas are allowed for balancing
#' @param zeta tuning parameter for selecting approximately balancing weights
#' @param fit.method the method used to fit mu(x) = E[YW | XW = x]
#' @param alpha tuning paramter for glmnet
#' @param optimizer which optimizer to use for approximate balancing
#' @param bound.gamma whether upper bound on gamma should be imposed
#' @param verbose whether the optimizer should print progress information
#' @param K hyperparameter for optimization
#'
#' @return Estimate for E[YW | XW = balance.target], along with variance estimate
#'
#' @export residualBalance.mean
residualBalance.mean = function(XW, YW,
                                balance.target,
                                allow.negative.weights = FALSE,
                                zeta,
                                K = 2.75,
                                fit.method = c("elnet", "none"),
                                alpha,
                                optimizer = c("mosek", "pogs", "pogs.dual", "quadprog", "pogs.new"),
                                bound.gamma = TRUE,
                                verbose = FALSE) {
  
  fit.method = match.arg(fit.method)
  optimizer = match.arg(optimizer)
  
  gamma = approx.balance(XW, balance.target, zeta = zeta, K, allow.negative.weights = allow.negative.weights, optimizer = optimizer, bound.gamma=bound.gamma, verbose=verbose)
  
  if (fit.method == "elnet") {
    
    lasso.fit = glmnet::cv.glmnet(XW, YW, alpha = alpha, nfold = 3)
    mu.lasso = predict(lasso.fit, newx = matrix(balance.target, 1, length(balance.target)))
    
    residuals = YW - predict(lasso.fit, newx = XW)
    mu.residual = sum(gamma * residuals)
    
    var.hat = sum(gamma^2 * residuals^2) *
      # degrees of freedom correction
      length(gamma) / max(1, length(gamma) - sum(coef(lasso.fit) != 0))
    
  } else if (fit.method == "none") {
    
    mu.lasso = 0
    mu.residual = sum(gamma * YW)
    
    var.hat = NA
    
  } else {
    
    stop("Invalid choice of fitting method.")
    
  }
  
  mu.hat = mu.lasso + mu.residual
  c(mu.hat, var.hat)
}

#' Compute approximately balancing weights
#'
#' Returns the minimizer of:
#'   (1 - zeta) ||gamma||^2 + zeta ||M'gamma - balance.target||_infty^2 (*)
#'
#' @param M the feature matrix, see (*)
#' @param balance.target the target solution, see (*)
#' @param zeta tuning parameter, see (*)
#' @param allow.negative.weights are the gammas allowed to be negative?
#' @param optimizer Which optimizer to use? Mosek is a commercial solver, but free
#'                  academic licenses are available. Needs to be installed separately.
#'                  Pogs runs ADMM and may be useful for large problems, and
#'                  must be installed separately. Quadprog is the default
#'                  R solver.
#' @param bound.gamma whether upper bound on gamma should be imposed
#' @param gamma.max specific upper bound for gamma (ignored if bound.gamma = FALSE)
#' @param verbose whether the optimizer should print progress information
#'
#' @return gamma, the minimizer of (*)
#'
#' @export approx.balance
approx.balance = function(M,
                          balance.target,
                          zeta = 0.5,
                          K = 1,
                          cv.K = FALSE,  # NEW: 是否对 K 进行交叉验证
                          allow.negative.weights = FALSE,
                          optimizer = c("mosek", "pogs", "pogs.dual", "quadprog", "pogs.new"),
                          bound.gamma = FALSE,
                          gamma.max = 1/nrow(M)^(2/3),
                          verbose = FALSE) {
  
  if (zeta <= 0 || zeta >= 1) {
    stop("approx.balance: zeta must be between 0 and 1")
  }
  
  optimizer = match.arg(optimizer)
  
  # 当 optimizer = "pogs.new" 且 cv.K = TRUE 时，进行搜索选择最优 K
  if (cv.K && optimizer == "pogs.new") {
    K_values <- seq(1, 3, length.out = 15)  # 生成 K 取值范围
    errors <- numeric(length(K_values))  # 存储 balance error
    
    for (i in seq_along(K_values)) {
      gamma_cv <- approx.balance.pogs.new(M, balance.target, K = K_values[i],
                                          allow.negative.weights = allow.negative.weights, 
                                          bound.gamma = bound.gamma, gamma.max = gamma.max)
      
      balance_error <- max(abs(t(M) %*% gamma_cv - balance.target))  # 计算 balance error
      errors[i] <- balance_error
    }
    
    # 选择 balance error 最小的 K
    K <- K_values[which.min(errors)]
    if (verbose) cat("Selected optimal K:", K, "\n")
  }
  
  if (optimizer == "mosek") {
    if (suppressWarnings(require("Rmosek", quietly = TRUE))) {
      gamma = approx.balance.mosek.dual(M, balance.target, zeta, allow.negative.weights, bound.gamma, gamma.max, verbose)
    } else {
      if (suppressWarnings(require("pogs", quietly = TRUE))) {
        warning("The mosek optimizer is not installed. Using pogs instead.")
        optimizer = "pogs"
      } else {
        warning("Neither mosek nor pogs optimizers are installed. Using quadprog instead.")
        optimizer = "quadprog"
      }
    }
  }
  
  if (optimizer %in% c("pogs", "pogs.dual", "pogs.new")) {
    if (suppressWarnings(require("pogs", quietly = TRUE))) {
      if (optimizer == "pogs") {
        gamma = approx.balance.pogs(M, balance.target, zeta, allow.negative.weights, bound.gamma, gamma.max, verbose)
      } else if (optimizer == "pogs.new"){
        gamma = approx.balance.pogs.new(M, balance.target, K, allow.negative.weights, bound.gamma, gamma.max)
      } else {
        if (bound.gamma) {warning("bound.gamma = TRUE not implemented for this optimizer")}
        gamma = approx.balance.pogs.dual(M, balance.target, zeta, allow.negative.weights, verbose)
      }
    } else {
      warning("The POGS optimizer is not installed. Using quadprog instead.")
      optimizer = "quadprog"
    }
  }
  
  if (optimizer == "quadprog") {
    if (bound.gamma) {warning("bound.gamma = TRUE not implemented for this optimizer")}
    gamma = approx.balance.quadprog(M, balance.target, zeta, allow.negative.weights)
  }
  
  gamma
}

#################################################################################################################
approx.balance.pogs.new = function(M,                # Design matrix X_c
                                   balance.target,   # Target vector 
                                   K,            # Trade-off parameter for balance constraints
                                   allow.negative.weights = FALSE,  # Whether gamma can be negative
                                   bound.gamma = FALSE,  # Whether to impose max gamma constraint
                                   verbose = FALSE) {
  
  # Ensure balance.target is a column vector
  balance.target = matrix(balance.target, ncol = 1) # xi
  
  # Compute the bound for the ||ξ - X_c^T γ||_∞ constraint
  nn = nrow(M)
  balance.bound = K * sqrt(log(ncol(M)) / (nn))
  
  # Compute the maximum allowed γ_i values, ensuring it's not too small
  gamma.max = max(nn^(-2/3))  # the gamma constraint
  
  # Define the quadratic loss function ||γ||_2^2
  g = list(h = kSquare())
  
  # Define constraints:
  # (1) Infinity norm constraint: ||ξ - X_c^T γ||_∞ ≤ balance.bound
  # (2) Sum constraint: sum(gamma_i) = 1
  f = list(h = c(kIndLe0(2 * ncol(M)), kIndEq0(2)),  # (1) 2p inequality constraints, (2) sum constraint
           b = c(rep(balance.bound, 2 * ncol(M)), 1, 1))  # Right-hand side values for the constraints
  
  # Construct the constraint matrix A
  # The first two blocks enforce ||M^T gamma - balance.target||_∞ ≤ balance.bound
  # The third block enforces sum(gamma) = 1
  A = rbind(
    cbind(t(M), -balance.target, -balance.bound),  # M^T γ - balance.target ≤ balance.bound
    cbind(-t(M), balance.target, -balance.bound),  # -M^T γ + balance.target ≤ balance.bound
    c(rep(1, nrow(M)), 0, 0),  # sum(gamma) = 1
    c(rep(0, nrow(M)), 1, 0)   # Auxiliary variable Z constraint
  )
  
  # If we do not allow negative weights, we enforce gamma ≥ 0
  if (!allow.negative.weights) {
    f$h = c(f$h, kIndLe0(nrow(M)))  # Add inequality constraint gamma ≥ 0
    f$b = c(f$b, rep(0, nrow(M)))  # Right-hand side for the new constraints
    A = rbind(A, cbind(diag(-1, nrow(M)), 0, 0))  # Add negative identity to enforce non-negativity
  }
  
  # If we bound gamma, we enforce gamma_i ≤ gamma.max
  if (bound.gamma) {
    f$h = c(f$h, kIndLe0(nrow(M)))  # Add upper bound constraint
    f$b = c(f$b, rep(0, nrow(M)))  # Right-hand side for the constraints
    A = rbind(A, cbind(diag(1, nrow(M)), -gamma.max, 0))  # Enforce gamma ≤ gamma.max
  }
  
  # Solve using POGS optimizer
  pogs.solution = tryCatch({
    pogs(A, f, g, params = list(rel_tol=1e-4, abs_tol=1e-5, verbose=2*as.numeric(verbose)))
  }, error = function(e) {
    cat("Warning: POGS failed, returning zero vector.\n")
    return(list(x = rep(0, nrow(M))))
  })
  
  # Extract the gamma solution and normalize it to sum to 1
  gamma = pogs.solution$x[1:nrow(M)]
  gamma / sum(gamma)
}



