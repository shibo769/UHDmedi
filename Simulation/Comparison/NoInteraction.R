############################################################
# High-dimensional mediation: no interaction
# Structure follows original code EXACTLY
############################################################

rm(list = ls())
source("/usr3/graduate/shibo/UHDmediation/utilis.R")

library(doRNG)
library(pracma)
library(MASS)

######################### CSV SAVE PATH ###################################

get_script_dir <- function() {
  
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  
  if (length(file_arg) > 0) {
    script_path <- sub("^--file=", "", file_arg[1])
    return(dirname(normalizePath(script_path)))
  }
  
  if (requireNamespace("rstudioapi", quietly = TRUE)) {
    if (rstudioapi::isAvailable()) {
      script_path <- rstudioapi::getActiveDocumentContext()$path
      if (!is.null(script_path) && nzchar(script_path)) {
        return(dirname(normalizePath(script_path)))
      }
    }
  }
  
  return(getwd())
}

script_dir <- get_script_dir()
cat("CSV files will be saved to:", script_dir, "\n")

######################### PARALLEL COMPUTING ###################################

list.of.packages <- c(
  "foreach",
  "doParallel",
  "glmnet",
  "tidyverse",
  "kableExtra"
)

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages) > 0){
  install.packages(new.packages, dep = TRUE)
}

for(package.i in list.of.packages){
  suppressPackageStartupMessages(library(package.i, character.only = TRUE))
}

n.cores <- 31
my.cluster <- parallel::makeCluster(n.cores, type = "PSOCK")
doParallel::registerDoParallel(cl = my.cluster)
registerDoRNG(2000)

################### Parameters ###################

Iter      <- 500
pq_vec   = c(750)# c(50, 400, 750, 1250)
n_vec    = c(1250)# c(250, 500, 750, 1000, 1250)

# s         <- 3
# s_B       <- 3

s_map <- list(
  "50"   = list(s = 5,  s_B = 5),
  "400"  = list(s = 10, s_B = 10),
  "750"  = list(s = 10, s_B = 10),
  "1250" = list(s = 20, s_B = 20)
)

sigma_eps <- 1
mu        <- 0.5
K1        <- 2.75
K2        <- 2.75

res <- list()

################### Simulation ###################
for (pq in pq_vec){
  
  p <- pq
  q <- pq
  
  # ---- set sparsity according to dimension ----
  s   <- s_map[[as.character(pq)]]$s
  s_B <- s_map[[as.character(pq)]]$s_B
  
  #### true beta ####
  beta <- rep(0, p)
  beta[sample(1:p, s, FALSE)] <- runif(s, 0.05, 0.2)
  
  #### gamma & delta1 with partial overlap ####
  idx_gamma  <- sample(1:q, s, FALSE)
  idx_shared <- sample(idx_gamma, floor(s/2), FALSE)
  idx_delta  <- c(
    idx_shared,
    sample(setdiff(1:q, idx_gamma), s - length(idx_shared), FALSE)
  )
  
  gamma <- rep(0, q)
  gamma[idx_gamma] <- runif(s, 0.05, 0.2)
  
  delta1 <- rep(0, q)
  delta1[idx_delta] <- runif(s, 0.05, 0.2)
  
  alpha1 <- 0.15
  
  TE_true  <- alpha1 + sum(delta1 * gamma)
  NDE_true <- alpha1
  NIE_true <- sum(delta1 * gamma)
  
  #### B matrix ####
  B <- matrix(0, q, p)
  
  for (iB in 1:q) {
    B[iB, sample(1:p, s_B, FALSE)] <- runif(s_B, 0.05, 0.2)
  }
  
  for (nn in n_vec){
    
    ### X ###
    mu_vec <- rep(mu, p)
    
    U_mat <- randortho(p)
    
    Sigma_X <- U_mat %*%
      diag(runif(p, 0.1, 0.5)) %*%
      t(U_mat)
    
    X <- mvrnorm(nn, mu_vec, Sigma_X)
    
    X_bar <- colMeans(X)
    
    ### A ###
    alpha_A <- rep(0, p)
    
    alpha_A[sample(1:p, s, FALSE)] <- runif(s, 0, 2)
    
    A <- rbinom(
      nn,
      1,
      1 / (1 + exp(-X %*% alpha_A))
    )
    
    ### Monte Carlo ###
    err_vec <- foreach(
      i = 1:Iter,
      .combine = "rbind"
    ) %dopar% {
      
      require(glmnet)
      require(pracma)
      require(MASS)
      
      ############################################################
      ### Generate Data
      ############################################################
      
      ### M ###
      lambda_U <- runif(q, 0.1, 0.5)
      
      R_U <- randortho(q)
      
      Sigma_U <- R_U %*%
        diag(lambda_U) %*%
        t(R_U)
      
      U <- mvrnorm(nn, rep(0, q), Sigma_U)
      
      M <- A %*% t(delta1) +
        X %*% t(B) +
        U
      
      ### Y ###
      Y <- alpha1 * A +
        X %*% beta +
        M %*% gamma +
        rnorm(nn, 0, sigma_eps)
      
      ############################################################
      ### Original Non-Cross-Fitted Estimator
      ############################################################
      
      #### Step 1: Total Effect ####
      W1 <- cbind(A, X)
      
      contrast_a <- c(1, rep(0, p))
      
      out_TE <- residualBalance.mean(
        W1,
        Y,
        balance.target = contrast_a,
        zeta = 0.5,
        K = K1,
        fit.method = "elnet",
        alpha = 1,
        optimizer = "pogs.new"
      )
      
      TE_hat <- out_TE[1]
      
      #### Step 2: Natural Direct Effect ####
      W2 <- cbind(A, X, M)
      
      contrast_b <- c(1, rep(0, p + q))
      
      out_NDE <- residualBalance.mean(
        W2,
        Y,
        balance.target = contrast_b,
        zeta = 0.5,
        K = K2,
        fit.method = "elnet",
        alpha = 1,
        optimizer = "pogs.new"
      )
      
      NDE_hat <- out_NDE[1]
      
      #### Step 3: NIE ####
      NIE_hat_naive <- TE_hat - NDE_hat
      
      ############################################################
      ### Cross Fitting Split
      ############################################################
      
      all_index <- sample(1:nn, nn, replace = FALSE)
      
      D1 <- all_index[1:floor(nn / 2)]
      
      D2 <- all_index[(floor(nn / 2) + 1):nn]
      
      ############################################################
      ### First Direction
      ### D1 estimates TE
      ### D2 estimates NDE
      ############################################################
      
      ### TE using D1 ###
      A_D1 <- A[D1]
      X_D1 <- X[D1, , drop = FALSE]
      Y_D1 <- Y[D1]
      
      W1_D1 <- cbind(A_D1, X_D1)
      
      contrast_a_1 <- c(1, rep(0, p))
      
      out_TE_1 <- residualBalance.mean(
        W1_D1,
        Y_D1,
        balance.target = contrast_a_1,
        zeta = 0.5,
        K = K1,
        fit.method = "elnet",
        alpha = 1,
        optimizer = "pogs.new"
      )
      
      TE_hat_1 <- out_TE_1[1]
      
      ### NDE using D2 ###
      A_D2 <- A[D2]
      X_D2 <- X[D2, , drop = FALSE]
      M_D2 <- M[D2, , drop = FALSE]
      Y_D2 <- Y[D2]
      
      W2_D2 <- cbind(A_D2, X_D2, M_D2)
      
      contrast_b_1 <- c(1, rep(0, p + q))
      
      out_NDE_1 <- residualBalance.mean(
        W2_D2,
        Y_D2,
        balance.target = contrast_b_1,
        zeta = 0.5,
        K = K2,
        fit.method = "elnet",
        alpha = 1,
        optimizer = "pogs.new"
      )
      
      NDE_hat_1 <- out_NDE_1[1]
      
      NIE_hat_1 <- TE_hat_1 - NDE_hat_1
      
      ############################################################
      ### Second Direction
      ### D2 estimates TE
      ### D1 estimates NDE
      ############################################################
      
      ### TE using D2 ###
      W1_D2 <- cbind(A_D2, X_D2)
      
      contrast_a_2 <- c(1, rep(0, p))
      
      out_TE_2 <- residualBalance.mean(
        W1_D2,
        Y_D2,
        balance.target = contrast_a_2,
        zeta = 0.5,
        K = K1,
        fit.method = "elnet",
        alpha = 1,
        optimizer = "pogs.new"
      )
      
      TE_hat_2 <- out_TE_2[1]
      
      ### NDE using D1 ###
      M_D1 <- M[D1, , drop = FALSE]
      
      W2_D1 <- cbind(A_D1, X_D1, M_D1)
      
      contrast_b_2 <- c(1, rep(0, p + q))
      
      out_NDE_2 <- residualBalance.mean(
        W2_D1,
        Y_D1,
        balance.target = contrast_b_2,
        zeta = 0.5,
        K = K2,
        fit.method = "elnet",
        alpha = 1,
        optimizer = "pogs.new"
      )
      
      NDE_hat_2 <- out_NDE_2[1]
      
      NIE_hat_2 <- TE_hat_2 - NDE_hat_2
      
      ############################################################
      ### Final Cross-Fitted Estimator
      ############################################################
      
      NIE_hat_cf <- 0.5 * (NIE_hat_1 + NIE_hat_2)
      
      ############################################################
      ### Bias
      ############################################################
      
      bias_NIE_cf <- NIE_hat_cf - NIE_true
      
      bias_NIE_naive <- NIE_hat_naive - NIE_true
      
      ############################################################
      ### Return
      ############################################################
      
      c(
        NIE_hat_cf,
        NIE_hat_naive,
        NIE_true,
        bias_NIE_cf,
        bias_NIE_naive
      )
    }
    
    ### Save every Monte Carlo iteration result
    err_df <- as.data.frame(err_vec)
    
    colnames(err_df) <- c(
      "NIE_hat_cf",
      "NIE_hat_naive",
      "NIE_true",
      "bias_NIE_cf",
      "bias_NIE_naive"
    )
    
    err_df$iter <- 1:nrow(err_df)
    err_df$pq <- pq
    err_df$n <- nn
    
    err_df <- err_df[, c(
      "pq",
      "n",
      "iter",
      "NIE_true",
      "NIE_hat_cf",
      "NIE_hat_naive",
      "bias_NIE_cf",
      "bias_NIE_naive"
    )]
    
    write.csv(
      err_df,
      file = file.path(script_dir, paste0("mc_results_pq_", pq, "_n_", nn, ".csv")),
      row.names = FALSE
    )
    
    cat(
      "Finished |",
      "pq =", pq,
      "| n =", nn,
      "| Mean Bias CF =", round(mean(err_vec[,4], na.rm = TRUE), 4),
      "| RMSE CF =", round(sqrt(mean(err_vec[,4]^2, na.rm = TRUE)), 4),
      "| SD CF =", round(sd(err_vec[,1], na.rm = TRUE), 4),
      "| Mean Bias Naive =", round(mean(err_vec[,5], na.rm = TRUE), 4),
      "| RMSE Naive =", round(sqrt(mean(err_vec[,5]^2, na.rm = TRUE)), 4),
      "\n"
    )
    
    res[[paste0("pq_", pq, "_n_", nn)]] <- err_vec
  }
}

################### Summary ###################

results_df <- data.frame(
  n = integer(),
  p = integer(),
  NIE_true = numeric(),
  bias_mean = numeric(),
  rmse = numeric(),
  sd = numeric(),
  ci_length = numeric(),
  coverage = numeric()
)

for (pq in pq_vec){
  for (nn in n_vec){
    
    key <- paste0("pq_", pq, "_n_", nn)
    
    if (is.null(res[[key]])) {
      cat("Missing:", key, "\n")
      next
    }
    
    seg <- res[[key]]
    
    if (nrow(seg) == 0) {
      cat("Empty:", key, "\n")
      next
    }
    
    NIE_hat  <- seg[,1]
    NIE_true <- seg[1,3]   # 修正这里
    
    bias <- NIE_hat - NIE_true
    
    bias_mean <- mean(bias, na.rm=TRUE)
    rmse <- sqrt(mean(bias^2, na.rm=TRUE))
    sd_hat <- sd(bias, na.rm=TRUE)
    
    ci_lower <- NIE_hat - 1.96 * sd_hat
    ci_upper <- NIE_hat + 1.96 * sd_hat
    
    ci_length <- mean(ci_upper - ci_lower, na.rm=TRUE)
    
    coverage <- mean(
      ci_lower <= NIE_true &
        NIE_true <= ci_upper,
      na.rm=TRUE
    )
    
    results_df <- rbind(
      results_df,
      data.frame(
        n = nn,
        p = pq,
        NIE_true = NIE_true,
        bias_mean = round(bias_mean,3),
        rmse = round(rmse,3),
        sd = round(sd_hat,3),
        ci_length = round(ci_length,3),
        coverage = round(coverage,3)
      )
    )
  }
}

print(results_df)

write.csv(
  results_df,
  file = file.path(script_dir, "summary_results_all_settings.csv"),
  row.names = FALSE
)

################### Plots ###################

# plot_df <- data.frame(
#   error = sqrt(n_vec[1]) * res[[1]][,4]
# )
# 
# p1 <- ggplot(plot_df, aes(x = error)) +
#   geom_histogram(bins = 30, fill = "cornflowerblue", alpha = 0.8) +
#   theme_minimal() +
#   ggtitle(expression(sqrt(n)*(hat(NIE)-NIE)))
# 
# p2 <- ggplot(plot_df, aes(sample = error)) +
#   stat_qq() +
#   stat_qq_line() +
#   theme_minimal() +
#   ggtitle("QQ plot")
# 
# grid.arrange(p1, p2, ncol = 2)

############################################################