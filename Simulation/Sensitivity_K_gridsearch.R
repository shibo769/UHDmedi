rm(list = ls())
source("/usr3/graduate/shibo/utilis.R")
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

n.cores = 29

my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)

doParallel::registerDoParallel(cl = my.cluster)
registerDoRNG()

################### Choose Parameters ###################
# Parameters
Iter     = 300
res      = list()
s        = 5
s_B      = 5
sigma_eps = 0.5
nn       = 1250
pq       = 750
K_grid   = expand.grid(K1 = seq(1, 3, 0.25), K2 = seq(1, 3, 0.25))
mu = 0.5

# Set p and q
p = pq
q = pq

### Generate True Parameters ###
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

### Generate B0 and B1 Matrices ###
B0 <- matrix(0, q, p)
B1 <- matrix(0, q, p)
for (i in 1:q) {
  B0[i, sample(1:p, s_B, replace = FALSE)] <- runif(s_B, 0.5, 1)
  B1[i, sample(1:p, s_B, replace = FALSE)] <- runif(s_B, 0.5, 1)
}

### Generate X ###
mu_vec <- rep(mu, p)
eigen_vec = diag(sqrt(runif(p, 1, 2)))
U_mat = randortho(p)
Sigma_X <- U_mat %*% eigen_vec %*% t(eigen_vec) %*% t(U_mat)
X <- mvrnorm(nn, mu_vec, Sigma_X)

### Generate Treatment Assignment A ###
alpha = rep(0, p)
alpha[sample(1:p, s, F)] = runif(s, 0, 2)

logit_vec <- 1 / (1 + exp(-X %*% alpha))
A <- as.numeric(rbinom(n = nn, 1, logit_vec))

treat_index   <- which(A == 1)
control_index <- which(A == 0)
nt <- length(treat_index)
nc <- length(control_index)

Xc <- X[control_index, ]
Xt <- X[treat_index, ]

### Compute True Value ###
X_bar = colMeans(X)
theta_0_true = t(X_bar) %*% (beta_1 + (t(B0) %*% gamma_1))

for (row in 1:nrow(K_grid)) {
  K1_use <- K_grid$K1[row]
  K2_use <- K_grid$K2[row]
  
  err_vec <- foreach(i = 1:Iter, .combine = 'rbind') %dopar% {
    require(MASS)
    require(glmnet)
    
    noise_Yt = rnorm(nt, 0, sigma_eps)
    noise_Yc = rnorm(nc, 0, sigma_eps)
    noise_M0 <- mvrnorm(nc, mu = rep(0, q), Sigma = diag(sigma_eps, q, q))
    noise_M1 <- mvrnorm(nt, mu = rep(0, q), Sigma = diag(sigma_eps, q, q))
    
    M0 <-  Xc %*% t(B0) + noise_M0
    M1 <-  Xt %*% t(B1) + noise_M1
    
    Wt = cbind(Xt, M1)
    
    Yt <-  Xt %*% beta_1 + M1 %*% gamma_1 + noise_Yt
    Yc <-  Xc %*% beta_0 + M0 %*% gamma_0 + noise_Yc
    
    ############################################################
    ### Build full sample M and Y for cross fitting
    ############################################################
    M_full <- matrix(0, nrow = nn, ncol = q)
    M_full[control_index, ] <- M0
    M_full[treat_index, ]   <- M1
    
    Y_full <- rep(0, nn)
    Y_full[control_index] <- as.vector(Yc)
    Y_full[treat_index]   <- as.vector(Yt)
    
    ############################################################
    ### Two-fold stratified cross-fitting split
    ############################################################
    treat_perm <- sample(treat_index, length(treat_index), replace = FALSE)
    control_perm <- sample(control_index, length(control_index), replace = FALSE)
    
    D1 <- c(
      treat_perm[seq_len(floor(length(treat_perm) / 2))],
      control_perm[seq_len(floor(length(control_perm) / 2))]
    )
    
    D2 <- setdiff(seq_len(nn), D1)
    
    ############################################################
    ### Helper function: one CF direction
    ############################################################
    estimate_one_direction <- function(fold_B, fold_eval) {
      
      ############################
      ### Step 1: estimate B0 on fold_B controls
      ############################
      fold_B_control <- intersect(fold_B, control_index)
      
      Xc_B <- X[fold_B_control, , drop = FALSE]
      M0_B <- M_full[fold_B_control, , drop = FALSE]
      
      B0_hat <- matrix(0, nrow = q, ncol = p)
      
      for (jj in 1:q) {
        Mi <- M0_B[, jj]
        
        cv_model <- cv.glmnet(
          Xc_B,
          Mi,
          alpha = 1,
          intercept = FALSE,
          nfolds = 3
        )
        
        B0_hat[jj, ] <- as.vector(coef(cv_model, s = "lambda.1se")[-1])
      }
      
      ############################
      ### Step 2: estimate gamma_1 on fold_eval treated
      ############################
      fold_eval_treat <- intersect(fold_eval, treat_index)
      
      Xt_eval <- X[fold_eval_treat, , drop = FALSE]
      M1_eval <- M_full[fold_eval_treat, , drop = FALSE]
      Yt_eval <- Y_full[fold_eval_treat]
      
      Wt_eval <- cbind(Xt_eval, M1_eval)
      
      fit_cv <- cv.glmnet(
        Wt_eval,
        Yt_eval,
        nfolds = 3,
        intercept = FALSE
      )
      
      coeff <- as.vector(coef(fit_cv, s = "lambda.1se"))[-1]
      beta_1_hat <- coeff[1:p]
      gamma_1_hat <- as.matrix(coeff[(p + 1):(p + q)])
      
      ############################
      ### Step 3: first residual balancing part
      ############################
      B0_Xbar <- B0_hat %*% X_bar
      contrast <- c(X_bar, B0_Xbar)
      
      output1 <- residualBalance.mean(
        Wt_eval,
        Yt_eval,
        balance.target = contrast,
        zeta = 0.5,
        K = K1_use,
        fit.method = "elnet",
        alpha = 1,
        optimizer = "pogs.new"
      )
      
      ############################
      ### Step 4: second residual balancing part
      ############################
      fold_eval_control <- intersect(fold_eval, control_index)
      
      Xc_eval <- X[fold_eval_control, , drop = FALSE]
      M0_eval <- M_full[fold_eval_control, , drop = FALSE]
      
      M0_gamma <- M0_eval %*% gamma_1_hat
      
      output2 <- residualBalance.mean(
        Xc_eval,
        M0_gamma,
        balance.target = X_bar,
        zeta = 0.5,
        K = K2_use,
        fit.method = "elnet",
        alpha = 1,
        optimizer = "pogs.new"
      )
      
      ############################
      ### Step 5: one-direction CF estimator
      ############################
      theta_hat <- output1[1] + output2[1] -
        t(B0_Xbar) %*% gamma_1_hat
      
      return(as.numeric(theta_hat))
    }
    
    ############################################################
    ### Two directions of cross fitting
    ############################################################
    theta_hat_1 <- estimate_one_direction(fold_B = D1, fold_eval = D2)
    theta_hat_2 <- estimate_one_direction(fold_B = D2, fold_eval = D1)
    
    ############################################################
    ### Final cross-fitted estimator
    ############################################################
    theta_0_hat_debiasing <- 0.5 * (theta_hat_1 + theta_hat_2)
    bias_debiasing <- theta_0_hat_debiasing - theta_0_true
    
    c(theta_0_hat_debiasing, theta_0_true, bias_debiasing)
  }
  
  cat("Loop Finished for K1 =", K1_use, "K2 =", K2_use, "\n")
  key <- paste0("K1_", K1_use, "_K2_", K2_use)
  res[[key]] <- err_vec
}

# Initialize results data frame
results_df <- data.frame(Iter = integer(), S = integer(), n = integer(), K1 = numeric(), K2 = numeric(), Var = numeric(), theta_bias_mean = numeric(), theta_bias_rmse = numeric(), theta_bias_sd = numeric())

# Compute Bias, RMSE, and SD
for (row in 1:nrow(K_grid)) {
  K1_use <- K_grid$K1[row]
  K2_use <- K_grid$K2[row]
  key <- paste0("K1_", K1_use, "_K2_", K2_use)
  
  if (!is.null(res[[key]])) {
    segment <- res[[key]]
    segment_bias_mean <- mean(segment[,3])
    segment_bias_rmse <- sqrt(mean(segment[,3]^2)) 
    segment_bias_sd <- sqrt(segment_bias_rmse^2 - segment_bias_mean^2)
    
    results_df <- rbind(results_df, data.frame(Iter = Iter, S = s, n = nn, K1 = K1_use, K2 = K2_use, Var = sigma_eps, theta_bias_mean = segment_bias_mean, theta_bias_rmse = round(segment_bias_rmse, 4), theta_bias_sd = round(segment_bias_sd, 4)))
  }
}

print(results_df)

sort(results_df$theta_bias_rmse)
