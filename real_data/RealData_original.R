################################################################################
### Real-data example (lung-cancer DNA-methylation mediation) using the
### ORIGINAL debiased algorithm (explicit B0 estimation), matching the user's
### RealDataNoBoot analysis. B0 is estimated by BIC-lasso (b0_bic_matrix),
### which we verified is equivalent to per-column cv.glmnet but ~5.6x faster.
###
### theta_10 = outY10_1 + outY10_2 - (B0_hat Xbar)^T gamma1_hat,  etc.
### NIE = theta_11 - theta_10, NDE = theta_10 - theta_00, TE = NIE + NDE.
### Inference: nonparametric bootstrap; report point estimate, bootstrap mean,
### and 2.5/97.5 quantile CI, plus histograms.
################################################################################

suppressPackageStartupMessages({
  library(glmnet); library(pogs)
  library(foreach); library(doParallel); library(doRNG)
})

repo <- dirname(SHARED)   # kept for any legacy refs; deps now via SHARED
SHARED <- file.path((function(){a<-commandArgs(FALSE);f<-sub("^--file=","",a[grep("^--file=",a)]);if(length(f)) dirname(normalizePath(f)) else normalizePath(getwd())})(), "..", "shared")   # bundled deps (was original_v1/*)
source(file.path(SHARED, "mc_common.R"))                 # env_*, cluster, b0_bic_matrix, safe_cv_glmnet
source(file.path(SHARED, "utilis.R"))       # residualBalance.mean / approx.balance
set_single_thread_blas()

################################################################################
### Configuration
################################################################################
n_boot  <- env_int("MC_BOOT", 300L)
seed    <- env_int("MC_SEED", 20260719L)
nfolds  <- env_int("MC_NFOLDS", 3L)
K1 <- env_num("MC_K1", 2.75); K2 <- env_num("MC_K2", 2.75)
zeta <- 0.5; optimizer <- "pogs.new"
lam  <- env_chr("MC_LAMBDA", "lambda.min")                 # matches the old code default
data_dir <- file.path(repo, "modified", "real_data")
out_dir  <- file.path(repo, "modified", "real_data_outputs")
ensure_dir(out_dir)

################################################################################
### One full estimate on (X, M, Y, A) using the original algorithm.
################################################################################
estimate_effects <- function(X, M, Y, A) {
  p <- ncol(X); q <- ncol(M)
  ti <- which(A == 1); ci <- which(A == 0)
  Xt <- X[ti, , drop = FALSE]; Xc <- X[ci, , drop = FALSE]
  M1 <- M[ti, , drop = FALSE]; M0 <- M[ci, , drop = FALSE]
  Yt <- Y[ti]; Yc <- Y[ci]
  X_bar <- colMeans(X)
  Wt <- cbind(Xt, M1); Wc <- cbind(Xc, M0)

  ## outcome lasso (treated / control) -> gamma1, gamma0
  f1 <- safe_cv_glmnet(Wt, Yt, alpha = 1, intercept = FALSE, nfolds = nfolds)
  if (is.null(f1)) stop("cv1 failed")
  phi1 <- as.vector(coef(f1, s = lam))[-1]; gamma1 <- phi1[(p + 1):(p + q)]
  f0 <- safe_cv_glmnet(Wc, Yc, alpha = 1, intercept = FALSE, nfolds = nfolds)
  if (is.null(f0)) stop("cv0 failed")
  phi0 <- as.vector(coef(f0, s = lam))[-1]; gamma0 <- phi0[(p + 1):(p + q)]

  ## explicit B0 / B1 via BIC-lasso (control / treated arms)
  B0 <- b0_bic_matrix(Xc, M0); B0_Xbar <- as.vector(B0 %*% X_bar)
  B1 <- b0_bic_matrix(Xt, M1); B1_Xbar <- as.vector(B1 %*% X_bar)

  ## ---- theta_10 ----
  o10a <- residualBalance.mean(Wt, Yt, balance.target = c(X_bar, B0_Xbar),
                               zeta = zeta, fit.method = "elnet", alpha = 1, K = K1, optimizer = optimizer)
  o10b <- residualBalance.mean(Xc, as.vector(M0 %*% gamma1), balance.target = X_bar,
                               zeta = zeta, fit.method = "elnet", alpha = 1, K = K2, optimizer = optimizer)
  theta10 <- o10a[1] + o10b[1] - sum(B0_Xbar * gamma1)

  ## ---- theta_11 ----
  o11a <- residualBalance.mean(Wt, Yt, balance.target = c(X_bar, B1_Xbar),
                               zeta = zeta, fit.method = "elnet", alpha = 1, K = K1, optimizer = optimizer)
  o11b <- residualBalance.mean(Xt, as.vector(M1 %*% gamma1), balance.target = X_bar,
                               zeta = zeta, fit.method = "elnet", alpha = 1, K = K2, optimizer = optimizer)
  theta11 <- o11a[1] + o11b[1] - sum(B1_Xbar * gamma1)

  ## ---- theta_00 ----
  o00a <- residualBalance.mean(Wc, Yc, balance.target = c(X_bar, B0_Xbar),
                               zeta = zeta, fit.method = "elnet", alpha = 1, K = K1, optimizer = optimizer)
  o00b <- residualBalance.mean(Xc, as.vector(M0 %*% gamma0), balance.target = X_bar,
                               zeta = zeta, fit.method = "elnet", alpha = 1, K = K2, optimizer = optimizer)
  theta00 <- o00a[1] + o00b[1] - sum(B0_Xbar * gamma0)

  NIE <- theta11 - theta10; NDE <- theta10 - theta00; TE <- NIE + NDE
  c(NIE = NIE, NDE = NDE, TE = TE, theta11 = theta11, theta10 = theta10, theta00 = theta00)
}

################################################################################
### Load data
################################################################################
cat("Loading data ...\n")
A_Smoking <- readRDS(file.path(data_dir, "A_Smoking.rds"))
X_all     <- data.matrix(readRDS(file.path(data_dir, "X_fourthorder_MICE.rds")))
Y_all     <- log(readRDS(file.path(data_dir, "Y_SurvivalTime.rds")))
M_raw     <- readRDS(file.path(data_dir, "M_raw_cleaned.rds"))
M_all     <- M_raw[, apply(M_raw, 2, sd) > 0.1]; rm(M_raw); gc()
nn <- length(A_Smoking)
cat(sprintf("n = %d | p = %d | q = %d | treated = %d | control = %d\n",
            nn, ncol(X_all), ncol(M_all), sum(A_Smoking == 1), sum(A_Smoking == 0)))

################################################################################
### Point estimate
################################################################################
set.seed(seed)
point_est <- estimate_effects(X_all, M_all, Y_all, A_Smoking)
cat("\nPoint estimate (original algorithm):\n"); print(round(point_est, 5))

################################################################################
### Bootstrap
################################################################################
if (n_boot > 0) {
  cat(sprintf("\nRunning %d bootstrap replicates ...\n", n_boot))
  cl <- make_mc_cluster(n_boot, seed)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  boot_mat <- foreach::foreach(bidx = 1:n_boot, .combine = rbind, .errorhandling = "remove") %dopar% {
    worker_single_thread_blas()
    idx <- sample(1:nn, size = nn, replace = TRUE)
    tryCatch(estimate_effects(X_all[idx, , drop = FALSE], M_all[idx, , drop = FALSE],
                              Y_all[idx], A_Smoking[idx]), error = function(e) NULL)
  }
  boot_mat <- as.matrix(boot_mat)
  cat(sprintf("Bootstrap done: %d / %d succeeded\n", nrow(boot_mat), n_boot))

  saveRDS(list(point = point_est, boot = boot_mat, n_boot = n_boot, K = c(K1, K2),
               lam = lam, seed = seed),
          file.path(out_dir, "real_data_original_results.rds"))

  summ <- function(nm, label) {
    v <- boot_mat[, nm]; qs <- quantile(v, c(0.025, 0.975), na.rm = TRUE)
    data.frame(Effect = label, Estimate = round(point_est[nm], 4),
               BootMean = round(mean(v, na.rm = TRUE), 4),
               CI_low = round(qs[1], 4), CI_high = round(qs[2], 4), row.names = NULL)
  }
  tab <- rbind(summ("NIE", "Indirect Effect"), summ("NDE", "Direct Effect"), summ("TE", "Total Effect"))
  cat("\n===== SUMMARY =====\n"); print(tab, row.names = FALSE)
  write.csv(tab, file.path(out_dir, "real_data_original_summary.csv"), row.names = FALSE)

  for (dev_fun in list(
    function() pdf(file.path(out_dir, "real_data_original_bootstrap_hist.pdf"), width = 11, height = 3.6),
    function() png(file.path(out_dir, "real_data_original_bootstrap_hist.png"), width = 1650, height = 540, res = 150))) {
    dev_fun(); par(mfrow = c(1, 3), mar = c(4, 4, 2, 1))
    for (nm in c("NIE", "NDE", "TE")) {
      v <- boot_mat[, nm]
      hist(v, breaks = 30, col = "#4a6fa5", border = "white", main = nm,
           xlab = expression(hat(theta)), ylab = "Count")
      abline(v = point_est[nm], col = "red", lwd = 2)
      abline(v = quantile(v, c(0.025, 0.975), na.rm = TRUE), col = "black", lty = 2)
    }
    dev.off()
  }
  cat("\nOutputs written to:", out_dir, "\n")
}
