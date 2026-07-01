############################################################
# NIE via TE - NDE using SCAD (ncvreg) + hand-written HBIC
# Sparsity varies with dimension p = q
# Save every Monte Carlo iteration result to CSV
############################################################

rm(list = ls())

library(doRNG)
library(pracma)
library(MASS)
source("/usr3/graduate/shibo/UHDmediation/Simulation/Revision/utils_mediation.R")

######################### Get R file directory #################################

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
  "ncvreg",
  "glmnet",
  "tidyverse",
  "kableExtra",
  "gridExtra"
)

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
if (length(new.packages) > 0) {
  install.packages(new.packages, dep = TRUE)
}

for (package.i in list.of.packages) {
  suppressPackageStartupMessages(
    library(package.i, character.only = TRUE)
  )
}

n.cores <- 31
my.cluster <- parallel::makeCluster(n.cores, type = "PSOCK")
doParallel::registerDoParallel(cl = my.cluster)
registerDoRNG(1998)

################### Parameters ################### 

Iter <- 500
pq_vec <- c(50, 400, 750, 1250)
n_vec <- c(250, 500, 750, 1000, 1250)

res <- list()

sigma_eps <- 1
mu <- 0.5
SCAD_a <- 3.7

s_map <- list(
  "50"   = list(s = 5,  s_B = 5),
  "400"  = list(s = 10, s_B = 10),
  "750"  = list(s = 10, s_B = 10),
  "1250" = list(s = 20, s_B = 20)
)

################### Simulation ###################

for (q in pq_vec) {
  
  p <- q
  
  s   <- s_map[[as.character(q)]]$s
  s_B <- s_map[[as.character(q)]]$s_B
  
  beta <- rep(0, p)
  beta[sample(1:p, s, FALSE)] <- runif(s, 0.05, 0.2)
  
  idx_gamma  <- sample(1:q, s, FALSE)
  idx_shared <- sample(idx_gamma, floor(s / 2), FALSE)
  idx_delta  <- c(
    idx_shared,
    sample(setdiff(1:q, idx_gamma), s - length(idx_shared), FALSE)
  )
  
  gamma <- rep(0, q)
  gamma[idx_gamma] <- runif(s, 0.05, 0.2)
  
  delta1 <- rep(0, q)
  delta1[idx_delta] <- runif(s, 0.05, 0.2)
  
  alpha1_true <- 0.15
  
  TE_true  <- alpha1_true + sum(delta1 * gamma)
  NDE_true <- alpha1_true
  NIE_true <- sum(delta1 * gamma)
  
  B <- matrix(0, q, p)
  for (iB in 1:q) {
    B[iB, sample(1:p, s_B, FALSE)] <- runif(s_B, 0.05, 0.2)
  }
  
  for (nn in n_vec) {
    
    mu_vec <- rep(mu, p)
    U_mat <- randortho(p)
    Sigma_X <- U_mat %*% diag(runif(p, 0.1, 0.5)) %*% t(U_mat)
    X <- mvrnorm(nn, mu_vec, Sigma_X)
    
    alpha_A <- rep(0, p)
    alpha_A[sample(1:p, s, FALSE)] <- runif(s, 0, 2)
    
    A <- rbinom(nn, 1, 1 / (1 + exp(-X %*% alpha_A)))
    
    err_vec <- foreach(
      i = 1:Iter,
      .combine = "rbind",
      .packages = c("MASS", "pracma", "glmnet", "ncvreg")
    ) %dopar% {
      
      lambda_U <- runif(q, 0.1, 0.5)
      R_U <- randortho(q)
      Sigma_U <- R_U %*% diag(lambda_U) %*% t(R_U)
      U <- mvrnorm(nn, rep(0, q), Sigma_U)
      
      M <- A %*% t(delta1) + X %*% t(B) + U
      Y <- alpha1_true * A + X %*% beta + M %*% gamma + rnorm(nn, 0, sigma_eps)
      
      lamb_grid <- seq(0.01, 0.25, length.out = 20)
      
      hbic_fit <- lapply(
        lamb_grid,
        HBIC_calc,
        xx = as.matrix(A),
        yy = Y,
        mm = M,
        S = X,
        n_imp = 0,
        penalize_S = TRUE
      )
      
      hbic_vals <- sapply(hbic_fit, function(z) z$BIC)
      id_best <- tail(which(hbic_vals == min(hbic_vals)), 1)
      result_best <- hbic_fit[[id_best]]
      
      active_M <- which(result_best$alpha0 != 0)
      active_S <- which(result_best$alpha2 != 0)
      
      M_refit <- if (length(active_M) == 0) {
        matrix(numeric(0), nrow = nn, ncol = 0)
      } else {
        M[, active_M, drop = FALSE]
      }
      
      S_refit <- if (length(active_S) == 0) {
        matrix(numeric(0), nrow = nn, ncol = 0)
      } else {
        X[, active_S, drop = FALSE]
      }
      
      infer <- mediationInference(
        X = as.matrix(A),
        Y = Y,
        M = M_refit,
        S = S_refit
      )
      
      NIE_hat <- as.numeric(infer$beta_hat)
      bias_NIE <- NIE_hat - NIE_true
      
      c(
        iter = i,
        NIE_hat = NIE_hat,
        NIE_true = NIE_true,
        bias_NIE = bias_NIE
      )
    }
    
    err_df <- as.data.frame(err_vec)
    
    err_df$q <- q
    err_df$p <- p
    err_df$n <- nn
    err_df$s <- s
    err_df$s_B <- s_B
    err_df$Var <- sigma_eps
    
    err_df <- err_df[, c(
      "iter",
      "n",
      "p",
      "q",
      "s",
      "s_B",
      "Var",
      "NIE_true",
      "NIE_hat",
      "bias_NIE"
    )]
    
    csv_name <- paste0(
      "MC_iter_results_SCAD_q_", q,
      "_n_", nn,
      ".csv"
    )
    # 
    # write.csv(
    #   err_df,
    #   file = file.path(script_dir, csv_name),
    #   row.names = FALSE
    # )
    
    res[[paste0("q_", q, "_n_", nn)]] <- err_df
    
    cat(
      "Finished q =", q,
      "n =", nn,
      "| CSV saved:",
      file.path(script_dir, csv_name),
      "\n"
    )
  }
}

parallel::stopCluster(my.cluster)

################### Summary ###################

results_df <- data.frame(
  n = integer(),
  p = integer(),
  q = integer(),
  s = integer(),
  s_B = integer(),
  Var = numeric(),
  
  NIE_true = numeric(),
  
  bias_mean_NIE = numeric(),
  rmse_NIE = numeric(),
  sd_NIE = numeric(),
  ci_length_NIE = numeric(),
  coverage_NIE = numeric()
)

for (q in pq_vec) {
  
  p <- q
  
  s   <- s_map[[as.character(q)]]$s
  s_B <- s_map[[as.character(q)]]$s_B
  
  for (nn in n_vec) {
    
    key <- paste0("q_", q, "_n_", nn)
    
    if (!is.null(res[[key]])) {
      
      segment <- res[[key]]
      
      NIE_hat  <- segment$NIE_hat
      NIE_true <- segment$NIE_true[1]
      
      bias_NIE <- NIE_hat - NIE_true
      
      bias_mean_NIE <- mean(bias_NIE, na.rm = TRUE)
      rmse_NIE <- sqrt(mean(bias_NIE^2, na.rm = TRUE))
      sd_NIE <- sd(bias_NIE, na.rm = TRUE)
      
      ci_lower_NIE <- NIE_hat - 1.96 * sd_NIE
      ci_upper_NIE <- NIE_hat + 1.96 * sd_NIE
      
      coverage_NIE <- mean(
        ci_lower_NIE <= NIE_true &
          NIE_true <= ci_upper_NIE,
        na.rm = TRUE
      )
      
      ci_length_NIE <- mean(
        ci_upper_NIE - ci_lower_NIE,
        na.rm = TRUE
      )
      
      results_df <- rbind(
        results_df,
        data.frame(
          n = nn,
          p = p,
          q = q,
          s = s,
          s_B = s_B,
          Var = sigma_eps,
          
          NIE_true = NIE_true,
          
          bias_mean_NIE = round(bias_mean_NIE, 3),
          rmse_NIE = round(rmse_NIE, 3),
          sd_NIE = round(sd_NIE, 3),
          ci_length_NIE = round(ci_length_NIE, 3),
          coverage_NIE = round(coverage_NIE, 3)
        )
      )
    }
  }
}

print(results_df)

summary_csv_name <- "Summary_results_SCAD_all_settings.csv"

write.csv(
  results_df,
  file = file.path(script_dir, summary_csv_name),
  row.names = FALSE
)

cat(
  "Summary CSV saved:",
  file.path(script_dir, summary_csv_name),
  "\n"
)

######################## Plot ########################

# res_plot_df <- data.frame()
# 
# for (q in pq_vec) {
#   for (nn in n_vec) {
#     
#     key <- paste0("q_", q, "_n_", nn)
#     
#     if (!is.null(res[[key]])) {
#       
#       segment <- res[[key]]
#       
#       temp_df <- data.frame(
#         pq = q,
#         n = nn,
#         bias_scaled = sqrt(nn) * segment$bias_NIE,
#         type = "NIE"
#       )
#       
#       res_plot_df <- rbind(res_plot_df, temp_df)
#     }
#   }
# }
# 
# res_plot_df$pq <- factor(res_plot_df$pq, levels = pq_vec)
# res_plot_df$n  <- factor(res_plot_df$n, levels = n_vec)
# 
# hist_plot <- ggplot(res_plot_df, aes(x = bias_scaled, fill = type)) +
#   geom_histogram(
#     bins = 35,
#     color = "#e9ecef",
#     alpha = 0.8,
#     position = "identity"
#   ) +
#   scale_fill_manual(values = c("cornflowerblue")) +
#   theme_minimal() +
#   ggtitle(expression(paste("Histogram of ", sqrt(n) * (hat(NIE) - NIE)))) +
#   xlab(expression(paste(sqrt(n) * (hat(NIE) - NIE)))) +
#   ylab("Frequency") +
#   theme(
#     plot.title = element_text(size = 15),
#     strip.text = element_text(size = 10)
#   )
# 
# qq_plot <- ggplot(res_plot_df, aes(sample = bias_scaled, color = type)) +
#   stat_qq(alpha = 0.8) +
#   stat_qq_line(color = "black", linetype = "dashed") +
#   scale_color_manual(values = c("cornflowerblue")) +
#   theme_minimal() +
#   ggtitle("QQ Plot with 45 Degree Line") +
#   xlab("Theoretical Normal Quantiles") +
#   ylab(expression(paste(sqrt(n) * (hat(NIE) - NIE)))) +
#   theme(
#     plot.title = element_text(size = 15),
#     strip.text = element_text(size = 10)
#   )
# 
# grid.arrange(hist_plot, qq_plot, ncol = 2)
