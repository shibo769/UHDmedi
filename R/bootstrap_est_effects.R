 
#' Bootstrap Estimation of Mediation Effects
#'
#' This function performs a non-parametric bootstrap to estimate the total effect (TE), natural indirect effect (NIE),
#' and natural direct effect (NDE), as well as intermediate debiased outcome estimates (Y11, Y10, Y00).
#' For each bootstrap sample, the function resamples the observations (with replacement) and calls \code{est_effects()}.
#' It then summarizes the bootstrap estimates by computing the mean, standard deviation, and 95\% confidence intervals
#' (using the 2.5\% and 97.5\% quantiles) for each effect.
#'
#' @param X Design matrix of dimension n x p.
#' @param Y Outcome vector of length n.
#' @param M Mediator matrix of dimension n x q.
#' @param A Treatment assignment vector (binary) of length n.
#' @param B Integer. Number of bootstrap replications (default = 200).
#' @param logY Logical. If TRUE, outcome Y is log-transformed (passed to \code{est_effects()}). Default is FALSE.
#' @param parallel Logical. If TRUE, parallel computing is used in \code{est_effects()} (and in each bootstrap iteration).
#'        Default is TRUE.
#' @param nfold1 Integer. Number of folds for cross-validation in outcome models. Default is 5.
#' @param nfold2 Integer. Number of folds for cross-validation in mediator models. Default is 3.
#' @param net_alpha Numeric. Mixing parameter for elastic net in \code{cv.glmnet} (default is 1 for lasso).
#' @param zeta Numeric. Tuning parameter for residual balancing. Default is 0.1.
#' @param optimizer Character. Optimizer method for residual balancing. Default is \code{"pogs"}.
#' @param intercept Logical. If TRUE, models include an intercept. Default is FALSE.
#'
#' @return A list containing:
#' \describe{
#'   \item{summary}{A data.frame summarizing the mean, standard deviation, lower and upper 95\% CI for each effect (TE, NIE, NDE, Y11, Y10, Y00).}
#'   \item{boot_samples}{A data.frame containing the bootstrap estimates for each replication. Columns include TE, NIE, NDE, Y11, Y10, and Y00.}
#' }
#'
#' @examples
#' \dontrun{
#'   # Suppose sim_data is generated from \code{simData()}
#'   sim_data <- simData()
#'   boot_results <- bootstrap_est_effects(X = sim_data$X,
#'                                        Y = sim_data$Y,
#'                                       M = sim_data$M,
#'                                       A = sim_data$A,
#'                                        B = 100,
#'                                        logY = FALSE,
#'                                        parallel = TRUE,
#'                                        nfold1 = 5,
#'                                        nfold2 = 3,
#'                                        net_alpha = 1,
#'                                        zeta = 0.5,
#'                                        optimizer = "pogs",
#'                                        intercept = FALSE,
#'                                        verbose = FALSE)
#'  
# cat("\nBootstrap summary:\n")
# print(boot_results$summary)
# cat("\nFirst few rows of bootstrap samples:\n")
# head(boot_results$boot_samples)
#'}
#' @export
bootstrap_est_effects <- function(X, Y, M, A, B = 200, logY = FALSE, parallel = TRUE,
                                  nfold1 = 5, nfold2 = 3, net_alpha = 1, zeta = 0.1,
                                  optimizer = "pogs", intercept = FALSE, verbose = FALSE,
                                  max_attempts = 5 * B) {
  n <- nrow(X)
  # Initialize matrix to store bootstrap estimates for each replication
  boot_mat <- matrix(NA, nrow = B, ncol = 6)
  colnames(boot_mat) <- c("TE", "NIE", "NDE", "Y11", "Y10", "Y00")
  
  # If parallel is TRUE, create the cluster once outside the loop.
  if (parallel) {
    num_cores <- parallel::detectCores() - 1
    cl <- parallel::makeCluster(num_cores)
    doParallel::registerDoParallel(cl)
    # Override the parallel parameter for est_effects to avoid re-creating the cluster in each loop.
    est_parallel <- FALSE
  } else {
    est_parallel <- FALSE
  }
  
  # Record start time for bootstrap iterations
  start_time <- Sys.time()
  b <- 1
  attempts <- 0
  
  while(b <= B && attempts < max_attempts) {
    attempts <- attempts + 1
    current_time <- Sys.time()
    elapsed <- as.numeric(difftime(current_time, start_time, units = "secs"))
    avg_time <- if(b > 0) elapsed / b else 0
    remaining_time <- avg_time * (B - b)
    
    message(sprintf("Bootstrap iteration %d of %d. Estimated remaining time: %.1f seconds.", 
                    b, B, remaining_time))
    if(verbose) {
      message(sprintf("Iteration %d: Resampling data and estimating effects...", b))
    }
    
    sample_indices <- sample(1:n, size = n, replace = TRUE)
    X_b <- X[sample_indices, , drop = FALSE]
    Y_b <- Y[sample_indices]
    M_b <- M[sample_indices, , drop = FALSE]
    A_b <- A[sample_indices]
    
    est_out <- est_effects(X = X_b, Y = Y_b, M = M_b, A = A_b,
                           logY = logY, parallel = est_parallel,
                           nfold1 = nfold1, nfold2 = nfold2, net_alpha = net_alpha,
                           zeta = zeta, optimizer = optimizer, intercept = intercept, verbose = verbose)
    
    if(is.null(est_out)) {
      message(sprintf("Iteration %d skipped due to gamma estimation issue.", b))
      next  # skip this iteraction
    }
    
    boot_mat[b, ] <- c(est_out$TE, est_out$NIE, est_out$NDE, est_out$Y11, est_out$Y10, est_out$Y00)
    b <- b + 1
  }
  
  if (parallel) {
    parallel::stopCluster(cl)
  }
  
  if(attempts >= max_attempts && b <= B) {
    message("Warning: Maximum attempts reached; bootstrap sample size may be smaller than intended.")
    boot_mat <- boot_mat[1:(b - 1), , drop = FALSE]
  }
  
  # Compute summary statistics for each effect
  summary_df <- data.frame(
    Effect = colnames(boot_mat),
    Mean = apply(boot_mat, 2, mean, na.rm = TRUE),
    SD = apply(boot_mat, 2, sd, na.rm = TRUE),
    Lower_95CI = apply(boot_mat, 2, quantile, probs = 0.025, na.rm = TRUE),
    Upper_95CI = apply(boot_mat, 2, quantile, probs = 0.975, na.rm = TRUE)
  )
  
  boot_samples_df <- as.data.frame(boot_mat)
  
  # Print a summary table using knitr::kable
  if (!requireNamespace("knitr", quietly = TRUE)) {
    install.packages("knitr")
  }
  cat("\n----------------------------------------------------------\n")
  cat("               Bootstrap Summary Statistics\n")
  cat("----------------------------------------------------------\n")
  print(knitr::kable(summary_df, caption = "Bootstrap Summary Statistics", align = "c"))
  
  return(list(
    summary = summary_df,
    boot_samples = boot_samples_df
  ))
}

