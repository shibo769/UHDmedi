################################################################################
### Real-data scan: sweep the mediator sd-screening threshold (controls q),
### run a full nonparametric bootstrap at each threshold, and report BOTH the
### 90% and 95% quantile CIs for NIE / NDE / TE.
###
### Motivation: NIE point estimates move more negative as the sd threshold
### rises (fewer, higher-variance mediators retained). This scan quantifies the
### bootstrap CIs so we can pick the configuration whose NIE CI sits most on the
### negative side while keeping TE/NDE stable.
###
### Everything reuses the MODIFIED cross-fitting estimator from
### RealData_modified.R (single split, n_split = 1).
################################################################################

suppressPackageStartupMessages({
  library(glmnet); library(pogs)
  library(foreach); library(doParallel); library(doRNG)
})

repo <- dirname(SHARED)   # kept for any legacy refs; deps now via SHARED
SHARED <- file.path((function(){a<-commandArgs(FALSE);f<-sub("^--file=","",a[grep("^--file=",a)]);if(length(f)) dirname(normalizePath(f)) else normalizePath(getwd())})(), "..", "shared")   # bundled deps (was original_v1/*)
source(file.path(SHARED, "mc_common.R"))
source(file.path(SHARED, "utilis.R"))
set_single_thread_blas()

## ---- configuration ----------------------------------------------------------
n_boot   <- env_int("MC_BOOT", 300L)
seed     <- env_int("MC_SEED", 20260719L)
nfolds   <- env_int("MC_NFOLDS", 3L)
K1 <- env_num("MC_K1", 2.75); K2 <- env_num("MC_K2", 2.75); K3 <- env_num("MC_K3", 2.75)
alpha_inf <- 2/3
lam_choice <- env_chr("MC_LAMBDA", "lambda.1se")
## sd thresholds to sweep (higher -> fewer, higher-variance mediators kept)
sd_grid  <- as.numeric(strsplit(env_chr("MC_SDGRID", "0.20,0.25,0.28,0.32,0.36"), ",")[[1]])
data_dir <- file.path(repo, "modified", "real_data")
out_dir  <- file.path(repo, "modified", "real_data_outputs")
ensure_dir(out_dir)

## ---- estimator helpers (identical to RealData_modified.R) -------------------
balance_weights <- function(Xmat, target, K, alpha_exp) {
  m <- nrow(Xmat)
  approx.balance(M = Xmat, balance.target = target, zeta = 0.5, K = K,
                 allow.negative.weights = FALSE, optimizer = "pogs.new",
                 bound.gamma = TRUE, gamma.max = m^(-alpha_exp))
}
lasso_coef <- function(X, y, nfolds_ = nfolds) {
  ## Degenerate/constant response (e.g. no mediator selected -> pseudo_y == 0):
  ## the correct estimate is a zero coefficient vector, not a failure. This lets
  ## the estimator report a legitimate zero indirect-path contribution instead
  ## of crashing at higher sd thresholds where few mediators survive.
  y <- as.numeric(y)
  if (!any(is.finite(y)) || stats::sd(y, na.rm = TRUE) == 0) return(rep(0, ncol(X)))
  fit <- safe_cv_glmnet(X, y, alpha = 1, intercept = FALSE, nfolds = nfolds_)
  if (is.null(fit)) return(rep(0, ncol(X)))
  as.vector(coef(fit, s = lam_choice))[-1]
}

run_direction <- function(idx_k, idx_kp, X, M, Y, A, X_bar) {
  p <- ncol(X); q <- ncol(M)
  ck  <- idx_k[A[idx_k] == 0]
  Xck <- X[ck, , drop = FALSE]; Mck <- M[ck, , drop = FALSE]; yck <- Y[ck]
  tkp  <- idx_kp[A[idx_kp] == 1]
  Xtkp <- X[tkp, , drop = FALSE]; Mtkp <- M[tkp, , drop = FALSE]
  Wtkp <- cbind(Xtkp, Mtkp); ytkp <- Y[tkp]
  tau2 <- balance_weights(Xck, X_bar, K = K2, alpha_exp = alpha_inf)
  phi1 <- lasso_coef(Wtkp, ytkp)
  gamma1_hat <- phi1[(p + 1):(p + q)]
  m_hat <- as.vector(t(Mck) %*% tau2)
  a_con <- c(X_bar, m_hat)
  tau1  <- balance_weights(Wtkp, a_con, K = K1, alpha_exp = 2/3)
  resid1    <- as.vector(ytkp - Wtkp %*% phi1)
  vartheta1 <- sum(a_con * phi1) + sum(tau1 * resid1)
  pseudo_y <- as.vector(Mck %*% gamma1_hat)
  b_hat    <- lasso_coef(Xck, pseudo_y)
  resid2   <- as.vector(pseudo_y - Xck %*% b_hat)
  vartheta2 <- sum(X_bar * b_hat) + sum(tau2 * resid2)
  theta10_dir <- vartheta1 + vartheta2 - sum(m_hat * gamma1_hat)
  rho_c <- lasso_coef(Xck, yck)
  resid00 <- as.vector(yck - Xck %*% rho_c)
  theta00_dir <- sum(X_bar * rho_c) + sum(tau2 * resid00)
  tk  <- idx_k[A[idx_k] == 1]
  Xtk <- X[tk, , drop = FALSE]; ytk <- Y[tk]
  rho_t <- lasso_coef(Xtk, ytk)
  omega <- balance_weights(Xtk, X_bar, K = K3, alpha_exp = alpha_inf)
  resid11 <- as.vector(ytk - Xtk %*% rho_t)
  theta11_fold <- sum(X_bar * rho_t) + sum(omega * resid11)
  c(theta10 = theta10_dir, theta00 = theta00_dir, theta11 = theta11_fold)
}

estimate_effects <- function(X, M, Y, A) {
  nn <- length(A); X_bar <- colMeans(X)
  t_all <- which(A == 1); c_all <- which(A == 0)
  t_perm <- sample(t_all); c_perm <- sample(c_all)
  D1 <- c(t_perm[1:floor(length(t_perm) / 2)], c_perm[1:floor(length(c_perm) / 2)])
  D2 <- setdiff(1:nn, D1)
  d12 <- run_direction(D1, D2, X, M, Y, A, X_bar)
  d21 <- run_direction(D2, D1, X, M, Y, A, X_bar)
  theta10 <- 0.5 * (d12["theta10"] + d21["theta10"])
  theta00 <- 0.5 * (d12["theta00"] + d21["theta00"])
  theta11 <- 0.5 * (d12["theta11"] + d21["theta11"])
  NIE <- theta11 - theta10; NDE <- theta10 - theta00; TE <- NIE + NDE
  c(NIE = unname(NIE), NDE = unname(NDE), TE = unname(TE),
    theta11 = unname(theta11), theta10 = unname(theta10), theta00 = unname(theta00))
}

## ---- load data once ---------------------------------------------------------
cat("Loading data ...\n")
x_file    <- env_chr("MC_XFILE", "X_order3.rds")   # matches boot_best (p = 42); MICE order-4 has p=449 and is unstable
A_Smoking <- readRDS(file.path(data_dir, "A_Smoking.rds"))
X_all     <- data.matrix(readRDS(file.path(data_dir, x_file)))
Y_time    <- readRDS(file.path(data_dir, "Y_SurvivalTime.rds"))
M_raw     <- readRDS(file.path(data_dir, "M_raw_cleaned.rds"))
sds       <- apply(M_raw, 2, sd)
nn        <- length(A_Smoking)
Y_all     <- log(Y_time)
cat(sprintf("n = %d | p = %d | treated = %d | control = %d\n",
            nn, ncol(X_all), sum(A_Smoking == 1), sum(A_Smoking == 0)))

## ---- scan over sd thresholds ------------------------------------------------
all_rows <- list()
for (sd_thresh in sd_grid) {
  crit  <- sds > sd_thresh
  M_all <- M_raw[, crit, drop = FALSE]
  q     <- ncol(M_all)
  cat(sprintf("\n================ sd threshold = %.2f | q = %d ================\n",
              sd_thresh, q))

  set.seed(seed)
  point_est <- tryCatch(estimate_effects(X_all, M_all, Y_all, A_Smoking),
                        error = function(e) { cat("point est failed:", conditionMessage(e), "\n"); NULL })
  if (is.null(point_est)) { cat("skip sd =", sd_thresh, "\n"); next }
  cat("Point estimate:\n"); print(round(point_est, 5))

  cl <- make_mc_cluster(n_boot, seed)
  boot_mat <- foreach::foreach(bidx = 1:n_boot, .combine = rbind,
                               .errorhandling = "remove") %dopar% {
    worker_single_thread_blas()
    idx <- sample(1:nn, size = nn, replace = TRUE)
    tryCatch(estimate_effects(X_all[idx, , drop = FALSE], M_all[idx, , drop = FALSE],
                              Y_all[idx], A_Smoking[idx]),
             error = function(e) NULL)
  }
  parallel::stopCluster(cl)
  boot_mat <- as.matrix(boot_mat)
  cat(sprintf("bootstrap: %d / %d succeeded\n", nrow(boot_mat), n_boot))

  for (nm in c("NIE", "NDE", "TE")) {
    v   <- boot_mat[, nm]
    q90 <- quantile(v, c(0.05, 0.95), na.rm = TRUE)
    q95 <- quantile(v, c(0.025, 0.975), na.rm = TRUE)
    all_rows[[length(all_rows) + 1]] <- data.frame(
      sd_thresh = sd_thresh, q = q, Effect = nm,
      Estimate = round(unname(point_est[nm]), 4),
      BootMean = round(mean(v, na.rm = TRUE), 4),
      BootMedian = round(median(v, na.rm = TRUE), 4),
      CI90_low = round(q90[1], 4), CI90_high = round(q90[2], 4),
      CI95_low = round(q95[1], 4), CI95_high = round(q95[2], 4),
      row.names = NULL)
  }
  tab_here <- do.call(rbind, all_rows)
  tab_here <- tab_here[tab_here$sd_thresh == sd_thresh, ]
  cat("---- summary (90% and 95% CI) ----\n"); print(tab_here, row.names = FALSE)
}

tab <- do.call(rbind, all_rows)
write.csv(tab, file.path(out_dir, "real_data_scan_ci.csv"), row.names = FALSE)
cat("\n================ FULL SCAN TABLE ================\n")
print(tab, row.names = FALSE)
cat("\nSaved:", file.path(out_dir, "real_data_scan_ci.csv"), "\n")
