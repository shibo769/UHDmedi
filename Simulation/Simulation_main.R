rm(list = ls())
source("/usr3/graduate/shibo/UHDmediation/utilis.R")
library(doRNG)
library(pracma)
library(MASS)
######################### PARALLEL COMPUTING ###################################

##### Preparing cluster ##########
list.of.packages <- c(
  "foreach",
  "doParallel",
  "ranger",
  "palmerpenguins",
  "tidyverse",
  "kableExtra"
)

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if(length(new.packages) > 0){
  install.packages(new.packages, dep=TRUE)
}

#loading packages
for(package.i in list.of.packages){
  suppressPackageStartupMessages(
    library(
      package.i, 
      character.only = TRUE
    )
  )
}

n.cores = 30

my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)

doParallel::registerDoParallel(cl = my.cluster)
registerDoRNG()

################### Choose Parameters ################### 
Iter     = 500
pq_vec   = c(50, 400, 750, 1250)
n_vec    = c(250, 500, 750, 1000, 1250)
res      = list()
s        = 3
s_B      = 3
sigma_eps = 1
mu = 0.5
K1 = 2.75
K2 = 2.75

for (pq in pq_vec){
  
  p = pq
  q = pq
  
  ### true beta_0 and beta_1 ###
  beta_0 = rep(0, p)
  
  beta_0_index = sample(1:p, s, F)
  beta_0[beta_0_index] = runif(s, 0.5, 1.5)
  
  beta_1 = rep(0, p)
  beta_1_index = sample(1:p, s, F)
  beta_1[beta_1_index] = runif(s, 0.5, 1.5)
  
  ### true gamma_0 and gamma_1 ###
  gamma_0 = rep(0, q)
  gamma_0_index = sample(1:q, s, F)
  gamma_0[gamma_0_index] = runif(s, 0.5, 1.5)
  
  gamma_1 = rep(0, q)
  gamma_1_index = sample(1:q, s, F)
  gamma_1[gamma_1_index] = runif(s, 0.5, 1.5)
  
  ### B0 and B1 ###
  B0 <- matrix(0, q, p)
  B1 <- matrix(0, q, p)
  for (i in 1:q) {
    B0[i, sample(1:p, s_B, replace = FALSE)] <- runif(s_B, 0.5, 1.5)
    B1[i, sample(1:p, s_B, replace = FALSE)] <- runif(s_B, 0.5, 1.5)
  }
  
  ### Loop over different sample sizes ###
  for (nn in n_vec){
    
    ### X generation ###
    mu_vec <- rep(mu, p)
    eigen_vec = diag(sqrt(runif(p, 1, 2)))
    U_mat = randortho(p)
    Sigma_X <- U_mat %*% eigen_vec %*% t(eigen_vec) %*% t(U_mat)
    X <- mvrnorm(nn, mu_vec, Sigma_X)
    
    ### A generation ###
    alpha = rep(0, p)
    alpha_index = sample(1:p, s, F)
    alpha[alpha_index] = runif(s, 0, 2)
    
    logit_vec <- 1 / (1 + exp(-X %*% alpha))
    A <- as.numeric(rbinom(n = nn, 1, logit_vec))
    
    treat_index   <- which(A == 1)
    control_index <- which(A == 0)
    nt <- length(treat_index)
    nc <- length(control_index)
    
    Xc <- X[control_index, ]
    Xt <- X[treat_index, ]
    
    ### true value calculation ###
    X_bar = colMeans(X)
    theta_0_true  = t(X_bar) %*% (beta_1 + (t(B0) %*% gamma_1))
    
    ### Monte Carlo Simulation ###
    err_vec <- foreach(
      i = 1:Iter, 
      .combine = 'rbind'
    ) %dopar% {
      require(MASS)
      require(glmnet)
      require(pracma)
      
      noise_Yt = rnorm(nt, 0, sigma_eps)
      noise_Yc = rnorm(nc, 0, sigma_eps)
      
      lambda_U  <- runif(q, 1, 2)
      D_half    <- diag(sqrt(lambda_U))
      R_U       <- randortho(q)
      Sigma_U  <- R_U %*% D_half %*% D_half %*% t(R_U)
      noise_M0 <- mvrnorm(nc, mu = rep(0, q), Sigma = Sigma_U)
      noise_M1 <- mvrnorm(nt, mu = rep(0, q), Sigma = Sigma_U)
      
      M0 <-  Xc %*% t(B0) + noise_M0
      M1 <-  Xt %*% t(B1) + noise_M1
      
      Wt = cbind(Xt, M1)
      
      Yt <-  Xt %*% beta_1 + M1 %*% gamma_1 + noise_Yt
      Yc <-  Xc %*% beta_0 + M0 %*% gamma_0 + noise_Yc
      
      ############################################################
      ### Build full sample M and Y for cross fitting
      ############################################################
      M <- matrix(0, nrow = nn, ncol = q)
      M[control_index, ] <- M0
      M[treat_index, ]   <- M1
      
      Y <- rep(0, nn)
      Y[control_index] <- as.vector(Yc)
      Y[treat_index]   <- as.vector(Yt)
      
      ############################################################
      ### Original naive estimator, kept for comparison
      ############################################################
      fit_cv = cv.glmnet(Wt, Yt, nfolds = 3, intercept = FALSE)
      coeff = as.vector(coef(fit_cv, s = "lambda.1se"))[-1]
      beta_1_hat = coeff[1:p]
      gamma_1_hat = as.matrix(coeff[(p+1):(p+q)])
      
      ### Parallelized regression for B0_hat ###
      B0_hat <- matrix(0, nrow = ncol(M0), ncol = ncol(Xc))
      
      for (j in 1:ncol(M0)){
        Mi <- M0[, j]
        cv_model <- cv.glmnet(Xc, Mi, alpha = 1, intercept = FALSE, nfolds = 3)
        B0_hat[j, ] <- as.vector(coef(cv_model, s = "lambda.1se")[-1])
      }
      
      theta_0_hat_naive = t(X_bar) %*% (beta_1_hat + (t(B0_hat) %*% gamma_1_hat))
      
      ############################################################
      ### Cross fitting split
      ############################################################
      all_index = sample(1:nn, nn, replace = FALSE)
      D1 = all_index[1:floor(nn / 2)]
      D2 = all_index[(floor(nn / 2) + 1):nn]
      
      ############################################################
      ### First direction: D1 estimates B0, D2 estimates rest
      ############################################################
      
      ### Step 1: use D1 controls to estimate B0
      D1_control = intersect(D1, control_index)
      Xc_D1 = X[D1_control, , drop = FALSE]
      Mc_D1 = M[D1_control, , drop = FALSE]
      
      B0_hat_1 <- matrix(0, nrow = q, ncol = p)
      
      for (j in 1:q){
        Mi <- Mc_D1[, j]
        cv_model <- cv.glmnet(Xc_D1, Mi, alpha = 1, intercept = FALSE, nfolds = 3)
        B0_hat_1[j, ] <- as.vector(coef(cv_model, s = "lambda.1se")[-1])
      }
      
      ### Step 2: use D2 treated to estimate phi_1
      D2_treat = intersect(D2, treat_index)
      Xt_D2 = X[D2_treat, , drop = FALSE]
      Mt_D2 = M[D2_treat, , drop = FALSE]
      Yt_D2 = Y[D2_treat]
      Wt_D2 = cbind(Xt_D2, Mt_D2)
      
      fit_cv_1 = cv.glmnet(Wt_D2, Yt_D2, nfolds = 3, intercept = FALSE)
      coeff_1 = as.vector(coef(fit_cv_1, s = "lambda.1se"))[-1]
      beta_1_hat_1 = coeff_1[1:p]
      gamma_1_hat_1 = as.matrix(coeff_1[(p+1):(p+q)])
      
      ### Step 3: compute tau_1 and theta_{0,1}
      contrast_1 = c(X_bar, B0_hat_1 %*% X_bar)
      
      output1_1 = residualBalance.mean(Wt_D2,
                                       Yt_D2,
                                       balance.target = contrast_1,
                                       zeta = 0.5, 
                                       K = K1,
                                       fit.method = "elnet",
                                       alpha = 1,
                                       optimizer = "pogs.new")
      
      theta_0_1_hat_1 = output1_1[1]
      
      ### Step 4: use D2 controls to estimate b = B0^T gamma_1
      D2_control = intersect(D2, control_index)
      Xc_D2 = X[D2_control, , drop = FALSE]
      Mc_D2 = M[D2_control, , drop = FALSE]
      
      M0_gamma_1 = Mc_D2 %*% gamma_1_hat_1
      
      fit_b_1 = cv.glmnet(Xc_D2, M0_gamma_1, alpha = 1, intercept = FALSE, nfolds = 3)
      b_hat_1 = as.vector(coef(fit_b_1, s = "lambda.1se")[-1])
      
      ### Step 5: compute tau_2 and theta_{0,2}
      output2_1 = residualBalance.mean(Xc_D2,
                                       M0_gamma_1,
                                       balance.target = X_bar,
                                       zeta = 0.5,
                                       K = K2,
                                       fit.method = "elnet",
                                       alpha = 1,
                                       optimizer = 'pogs.new')
      
      theta_0_2_hat_1 = output2_1[1]
      
      ### Step 6: first cross-fitted estimator
      theta_0_hat_1 = theta_0_1_hat_1 + theta_0_2_hat_1 -
        t(B0_hat_1 %*% X_bar) %*% gamma_1_hat_1
      
      ############################################################
      ### Second direction: D2 estimates B0, D1 estimates rest
      ############################################################
      
      ### Step 1: use D2 controls to estimate B0
      D2_control = intersect(D2, control_index)
      Xc_D2 = X[D2_control, , drop = FALSE]
      Mc_D2 = M[D2_control, , drop = FALSE]
      
      B0_hat_2 <- matrix(0, nrow = q, ncol = p)
      
      for (j in 1:q){
        Mi <- Mc_D2[, j]
        cv_model <- cv.glmnet(Xc_D2, Mi, alpha = 1, intercept = FALSE, nfolds = 3)
        B0_hat_2[j, ] <- as.vector(coef(cv_model, s = "lambda.1se")[-1])
      }
      
      ### Step 2: use D1 treated to estimate phi_1
      D1_treat = intersect(D1, treat_index)
      Xt_D1 = X[D1_treat, , drop = FALSE]
      Mt_D1 = M[D1_treat, , drop = FALSE]
      Yt_D1 = Y[D1_treat]
      Wt_D1 = cbind(Xt_D1, Mt_D1)
      
      fit_cv_2 = cv.glmnet(Wt_D1, Yt_D1, nfolds = 3, intercept = FALSE)
      coeff_2 = as.vector(coef(fit_cv_2, s = "lambda.1se"))[-1]
      beta_1_hat_2 = coeff_2[1:p]
      gamma_1_hat_2 = as.matrix(coeff_2[(p+1):(p+q)])
      
      ### Step 3: compute tau_1 and theta_{0,1}
      contrast_2 = c(X_bar, B0_hat_2 %*% X_bar)
      
      output1_2 = residualBalance.mean(Wt_D1,
                                       Yt_D1,
                                       balance.target = contrast_2,
                                       zeta = 0.5, 
                                       K = K1,
                                       fit.method = "elnet",
                                       alpha = 1,
                                       optimizer = "pogs.new")
      
      theta_0_1_hat_2 = output1_2[1]
      
      ### Step 4: use D1 controls to estimate b = B0^T gamma_1
      D1_control = intersect(D1, control_index)
      Xc_D1 = X[D1_control, , drop = FALSE]
      Mc_D1 = M[D1_control, , drop = FALSE]
      
      M0_gamma_2 = Mc_D1 %*% gamma_1_hat_2
      
      fit_b_2 = cv.glmnet(Xc_D1, M0_gamma_2, alpha = 1, intercept = FALSE, nfolds = 3)
      b_hat_2 = as.vector(coef(fit_b_2, s = "lambda.1se")[-1])
      
      ### Step 5: compute tau_2 and theta_{0,2}
      output2_2 = residualBalance.mean(Xc_D1,
                                       M0_gamma_2,
                                       balance.target = X_bar,
                                       zeta = 0.5,
                                       K = K2,
                                       fit.method = "elnet",
                                       alpha = 1,
                                       optimizer = 'pogs.new')
      
      theta_0_2_hat_2 = output2_2[1]
      
      ### Step 6: second cross-fitted estimator
      theta_0_hat_2 = theta_0_1_hat_2 + theta_0_2_hat_2 -
        t(B0_hat_2 %*% X_bar) %*% gamma_1_hat_2
      
      ############################################################
      ### Final cross-fitted estimator
      ############################################################
      theta_0_hat_debiasing = 0.5 * (theta_0_hat_1 + theta_0_hat_2)
      
      bias_debiasing = theta_0_hat_debiasing - theta_0_true
      bias_naive = theta_0_hat_naive - theta_0_true
      
      c(theta_0_hat_debiasing, theta_0_hat_naive, theta_0_true, bias_debiasing, bias_naive)
    }
    
    cat(
      "Loop Finished |",
      "pq =", pq,
      "| n =", nn,
      "| sigma =", sigma_eps,
      "| Mean Bias CF =", round(mean(err_vec[,4], na.rm = TRUE), 4),
      "| RMSE CF =", round(sqrt(mean(err_vec[,4]^2, na.rm = TRUE)), 4),
      "| SD CF =", round(sd(err_vec[,1], na.rm = TRUE), 4),
      "| Valid CF =", sum(!is.na(err_vec[,1])),
      "/", Iter,
      "\n"
    )
    res[[paste0("pq_", pq, "_n_", nn)]] <- err_vec
  }
}

# Initialize the data frame
results_df <- data.frame(
  # Iter = integer(), 
  # S = integer(), 
  n = integer(), 
  p = integer(), 
  q = integer(), 
  Var = numeric(),
  theta_true = numeric(),
  
  # Debiased
  bias_mean = numeric(),
  rmse = numeric(),
  sd = numeric(),
  ci_lower = numeric(),
  ci_upper = numeric(),
  ci_length = numeric(),
  coverage = numeric(),
  
  # Naive
  bias_mean_naive = numeric(),
  rmse_naive = numeric(),
  sd_naive = numeric(),
  ci_length_naive = numeric(),
  coverage_naive = numeric()
)

# Iterate over pq_vec and n_vec
for (pq in pq_vec) {
  for (nn in n_vec) {
    
    key <- paste0("pq_", pq, "_n_", nn)
    
    if (!is.null(res[[key]])) {
      
      segment <- res[[key]]
      
      # Extract estimates
      theta_hat_db <- segment[, 1]
      theta_hat_nv <- segment[, 2]
      theta_true   <- segment[, 3][1]
      
      ###################### Debiased ######################
      bias_mean <- mean(segment[,4]) # mean(theta_hat_db - theta_true)
      rmse      <- sqrt(mean(segment[,4]^2)) # sqrt(mean((theta_hat_db - theta_true)^2))
      sd_hat    <- sqrt(rmse^2 - bias_mean^2) 
      
      ci_lower  <- theta_hat_db - 1.96 * sd_hat
      ci_upper  <- theta_hat_db + 1.96 * sd_hat
      coverage  <- mean(ci_lower <= theta_true & theta_true <= ci_upper)
      ci_length <- mean(ci_upper - ci_lower)
      
      ###################### Naive ######################
      bias_mean_naive <- mean(segment[,5]) # mean(theta_hat_nv - theta_true)
      rmse_naive      <- sqrt(mean(segment[,5]^2)) # sqrt(mean((theta_hat_nv - theta_true)^2))
      sd_naive        <- sqrt(rmse_naive^2 - bias_mean_naive^2) # sd(theta_hat_nv)
      
      ci_lower_naive  <- theta_hat_nv - 1.96 * sd_naive
      ci_upper_naive  <- theta_hat_nv + 1.96 * sd_naive
      coverage_naive  <- mean(ci_lower_naive <= theta_true & theta_true <= ci_upper_naive)
      ci_length_naive <- mean(ci_upper_naive - ci_lower_naive)
      
      # Store results
      results_df <- rbind(
        results_df,
        data.frame(
          # Iter = Iter,
          # S = s,
          n = nn,
          p = pq,
          q = pq,
          Var = sigma_eps,
          theta_true = theta_true,
          
          bias_mean = round(bias_mean, 3),
          rmse = round(rmse, 3),
          sd = round(sd_hat, 3),
          ci_lower = round(mean(ci_lower), 3),
          ci_upper = round(mean(ci_upper), 3),
          ci_length = round(ci_length, 3),
          coverage = round(coverage, 3),
          
          bias_mean_naive = round(bias_mean_naive, 3),
          rmse_naive = round(rmse_naive, 3),
          sd_naive = round(sd_naive, 3),
          ci_length_naive = round(ci_length_naive, 3),
          coverage_naive = round(coverage_naive, 3)
        )
      )
    }
  }
}

# Print the results table
print(results_df)



# library(ggplot2)
library(hrbrthemes)
library(gridExtra)

# Prepare Data for Plotting
res_plot_df <- data.frame()

for (pq in pq_vec) {
  for (nn in n_vec) {
    key <- paste0("pq_", pq, "_n_", nn)
    if (!is.null(res[[key]])) {
      segment <- res[[key]]

      temp_df <- data.frame(
        pq = pq,
        n = nn,
        bias_scaled = sqrt(nn) * segment[, 4],  # Debiasing Bias
        type = "Debiasing"
      )

      temp_df_naive <- data.frame(
        pq = pq,
        n = nn,
        bias_scaled = sqrt(nn) * segment[, 5],  # Naive Bias
        type = "Naive"
      )

      res_plot_df <- rbind(res_plot_df, temp_df, temp_df_naive)
    }
  }
}

########################Single Graph without Grid#################################

# Extract data for pq = 400, n = 1000
res_debiasing <- sqrt(1000) * res$pq_400_n_1000[, 4]  # Debiasing errors
res_naive <- sqrt(1000) * res$pq_400_n_1000[, 5]  # Naive errors

# Create combined data frame
res_df <- data.frame(
  x = c(res_debiasing, res_naive),
  type = rep(c("Debiasing", "Naive"), each = length(res_debiasing))
)

# 增加 x 轴刻度的版本
hist_plot <- ggplot(res_df, aes(x = x, fill = type)) +
  geom_histogram(bins = 35, color = "#e9ecef", alpha = 0.8, position = "identity") +
  scale_fill_manual(values = c("cornflowerblue", "tomato")) +
  scale_x_continuous(
    breaks = seq(floor(min(res_df$x)), ceiling(max(res_df$x)), by = 8)
  ) +
  theme_minimal() +
  ggtitle(expression(paste("Histogram of ", sqrt(n) * (hat(theta) - theta)))) +
  xlab(expression(paste(sqrt(n) * (hat(theta) - theta)))) +
  ylab("Frequency") +
  theme(
    plot.title = element_text(size = 15),
    strip.text = element_text(size = 10)
  )

# Compute theoretical normal quantiles based on Debiasing
qq_debiasing <- qqnorm(res_debiasing, plot.it = FALSE)

res$pq_400_n_1000[,1]

# Create data frame for plotting
qq_df <- data.frame(
  theoretical_quantiles = qq_debiasing$x,
  sample_quantiles = qq_debiasing$y
)

# # Compute best-fit line for Debiasing
# debiasing_slope <- sd(res_debiasing)
# debiasing_intercept <- mean(res_debiasing)

# QQ Plot with properly aligned reference line
qq_plot <- ggplot(qq_df, aes(x = theoretical_quantiles, y = sample_quantiles)) +
  geom_point(color = "cornflowerblue", alpha = 0.8) +  # Scatter plot for QQ points
  # geom_abline(slope = debiasing_slope, intercept = debiasing_intercept, color = "black", linetype = "dashed") +  # Reference line based on Debiasing
  theme_minimal() +
  ggtitle("QQ Plot with 45 Degree Line") +
  xlab("Theoretical Normal Quantiles") +
  ylab(expression(paste(sqrt(n) * (hat(theta) - theta)))) +
  theme(
    plot.title = element_text(size = 15),
    strip.text = element_text(size = 10)
  )
grid.arrange(hist_plot, qq_plot, ncol = 2)

########################################################
# Compute theoretical normal quantiles based on Debiasing
qq_debiasing <- qqnorm(res_debiasing, plot.it = FALSE)
qq_naive <- qqnorm(res_naive, plot.it = FALSE)

# Create data frame for plotting
qq_df <- data.frame(
  theoretical_quantiles = c(qq_debiasing$x, qq_naive$x),
  sample_quantiles = c(qq_debiasing$y, qq_naive$y),
  type = rep(c("Debiasing", "Naive"), each = length(res_debiasing))
)

# Compute best-fit line for Debiasing
debiasing_slope <- sd(res_debiasing)
debiasing_intercept <- mean(res_debiasing)

# QQ Plot with properly aligned reference line
qq_plot <- ggplot(qq_df, aes(x = theoretical_quantiles, y = sample_quantiles, color = type)) +
  geom_point(alpha = 0.8) +  # Scatter plot for QQ points
  geom_abline(slope = debiasing_slope, intercept = debiasing_intercept, color = "black", linetype = "dashed") +  # Reference line based on Debiasing
  scale_color_manual(values = c("cornflowerblue", "tomato")) +
  theme_minimal() +
  ggtitle("Improved QQ Plot with 45 Degree Line") +
  xlab("Theoretical Normal Quantiles") +
  ylab(expression(paste(sqrt(n) * (hat(theta) - theta)))) +
  theme(
    plot.title = element_text(size = 15),
    strip.text = element_text(size = 10)
  )

# Display plots
grid.arrange(hist_plot, qq_plot, ncol = 2)
