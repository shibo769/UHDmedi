################################################################################
### Real-data example (lung-cancer DNA-methylation mediation) using the
### MODIFIED cross-fitting algorithm (no B0 estimation; contrast m_hat = M_c^T tau2).
###
### Estimands: theta_11, theta_10, theta_00 -> NIE = theta_11 - theta_10,
###            NDE = theta_10 - theta_00, TE = NIE + NDE.
### Inference: nonparametric bootstrap (resample n rows with replacement),
###            report bootstrap mean and 2.5%/97.5% quantile CI, plus the
###            point estimate on the original sample.
###
### Env vars: MC_BOOT (bootstrap reps, default 300), MC_CORES, MC_K1/K2/K3,
###           MC_SEED, MC_NFOLDS.
################################################################################

suppressPackageStartupMessages({                           # quiet package loads
  library(glmnet); library(pogs)                           # lasso + balancing solver
  library(foreach); library(doParallel); library(doRNG)    # parallel bootstrap backend (%dopar%)
})

repo <- dirname(SHARED)   # kept for any legacy refs; deps now via SHARED
SHARED <- file.path((function(){a<-commandArgs(FALSE);f<-sub("^--file=","",a[grep("^--file=",a)]);if(length(f)) dirname(normalizePath(f)) else normalizePath(getwd())})(), "..", "shared")   # bundled deps (was original_v1/*)
source(file.path(SHARED, "mc_common.R"))                 # env_* helpers, cluster, BLAS control
source(file.path(SHARED, "utilis.R"))       # approx.balance / residualBalance.mean

set_single_thread_blas()                                   # single-thread BLAS in the parent

################################################################################
### Configuration
################################################################################
n_boot  <- env_int("MC_BOOT", 300L)                        # number of bootstrap replicates
seed    <- env_int("MC_SEED", 20260719L)                   # RNG seed
nfolds  <- env_int("MC_NFOLDS", 3L)                        # cv.glmnet folds
K1 <- env_num("MC_K1", 2.75); K2 <- env_num("MC_K2", 2.75); K3 <- env_num("MC_K3", 2.75)  # balancing constants
alpha_inf <- 2/3                                           # exponent in the ||tau||_inf <= m^{-alpha} cap
data_dir <- file.path(repo, "modified", "real_data")       # data location
out_dir  <- file.path(repo, "modified", "real_data_outputs")           # output location
ensure_dir(out_dir)                                        # create output folder if missing

################################################################################
### Helpers (copied from Simu_modified.R so this script never runs the sim)
################################################################################
balance_weights <- function(Xmat, target, K, alpha_exp) {  # raw balancing-weight wrapper
  m <- nrow(Xmat)                                          # sample size for the inf bound
  approx.balance(M = Xmat, balance.target = target, zeta = 0.5, K = K,
                 allow.negative.weights = FALSE, optimizer = "pogs.new",
                 bound.gamma = TRUE, gamma.max = m^(-alpha_exp))
}
lam_choice <- env_chr("MC_LAMBDA", "lambda.1se")           # lambda rule: lambda.1se (default) or lambda.min (old code)
lasso_coef <- function(X, y, nfolds_ = nfolds) {           # cv.glmnet, no intercept, with robust fallback
  fit <- safe_cv_glmnet(X, y, alpha = 1, intercept = FALSE, nfolds = nfolds_)
  if (!is.null(fit)) return(as.vector(coef(fit, s = lam_choice))[-1])
  # Fallback 1: retry cv.glmnet with more folds (5) a couple of times (fold randomness)
  for (tryi in 1:3) {
    fit2 <- tryCatch(glmnet::cv.glmnet(as.matrix(X), y, alpha = 1, intercept = FALSE, nfolds = 5),
                     error = function(e) NULL)
    if (!is.null(fit2)) return(as.vector(coef(fit2, s = lam_choice))[-1])
  }
  # Fallback 2: plain glmnet on a lambda path, pick a moderate lambda (no CV)
  g <- tryCatch(glmnet::glmnet(as.matrix(X), y, alpha = 1, intercept = FALSE),
                error = function(e) NULL)
  if (!is.null(g)) {
    lam <- g$lambda[max(1, floor(length(g$lambda) * 0.5))]   # mid-path lambda
    return(as.vector(coef(g, s = lam))[-1])
  }
  # Fallback 3: ridge-regularized least squares with adaptive ridge (never singular)
  Xm <- as.matrix(X); pp <- ncol(Xm)
  G <- crossprod(Xm)
  lam_ridge <- 1e-3 * mean(diag(G)) + 1e-6          # scale ridge to the data
  as.vector(solve(G + lam_ridge * diag(pp), crossprod(Xm, y)))
}

################################################################################
### One cross-fitting direction (k -> k') of modified Algorithms 1-3.
### Returns this direction's theta_10, theta_00, theta_11.
################################################################################
run_direction <- function(idx_k, idx_kp, X, M, Y, A, X_bar) {  # one direction
  p <- ncol(X); q <- ncol(M)                               # dimensions

  ck  <- idx_k[A[idx_k] == 0]                              # fold-k control rows
  Xck <- X[ck, , drop = FALSE]; Mck <- M[ck, , drop = FALSE]; yck <- Y[ck]

  tkp  <- idx_kp[A[idx_kp] == 1]                           # fold-k' treated rows
  Xtkp <- X[tkp, , drop = FALSE]; Mtkp <- M[tkp, , drop = FALSE]
  Wtkp <- cbind(Xtkp, Mtkp); ytkp <- Y[tkp]

  ## Alg 1 Step 1: tau2 balances fold-k control X to Xbar
  tau2 <- balance_weights(Xck, X_bar, K = K2, alpha_exp = alpha_inf)

  ## Alg 1 Step 2: phi1 = (beta1, gamma1) lasso on fold-k' treated (X, M)
  phi1 <- lasso_coef(Wtkp, ytkp)
  gamma1_hat <- phi1[(p + 1):(p + q)]                      # mediator coefficients

  ## Alg 1 Step 3: contrast m_hat = M_{c,k}^T tau2 (replaces B0_hat %*% Xbar), then tau1 + vartheta1
  m_hat <- as.vector(t(Mck) %*% tau2)                      # weighted-mediator contrast
  a_con <- c(X_bar, m_hat)                                 # contrast a = (Xbar, m_hat)
  tau1  <- balance_weights(Wtkp, a_con, K = K1, alpha_exp = 2/3)
  resid1    <- as.vector(ytkp - Wtkp %*% phi1)
  vartheta1 <- sum(a_con * phi1) + sum(tau1 * resid1)

  ## Alg 1 Step 4: b-hat on fold-k control, pseudo-response M_{c,k} gamma1
  pseudo_y <- as.vector(Mck %*% gamma1_hat)
  b_hat    <- lasso_coef(Xck, pseudo_y)
  resid2   <- as.vector(pseudo_y - Xck %*% b_hat)
  vartheta2 <- sum(X_bar * b_hat) + sum(tau2 * resid2)     # tau2 reused

  ## Alg 1 Step 5: theta_10 for this direction
  theta10_dir <- vartheta1 + vartheta2 - sum(m_hat * gamma1_hat)

  ## Alg 2: theta_00 (control outcome), reuse tau2
  rho_c <- lasso_coef(Xck, yck)
  resid00 <- as.vector(yck - Xck %*% rho_c)
  theta00_dir <- sum(X_bar * rho_c) + sum(tau2 * resid00)

  ## Alg 3: theta_11 (treated outcome, X only), own omega on fold-k treated
  tk  <- idx_k[A[idx_k] == 1]
  Xtk <- X[tk, , drop = FALSE]; ytk <- Y[tk]
  rho_t <- lasso_coef(Xtk, ytk)
  omega <- balance_weights(Xtk, X_bar, K = K3, alpha_exp = alpha_inf)
  resid11 <- as.vector(ytk - Xtk %*% rho_t)
  theta11_fold <- sum(X_bar * rho_t) + sum(omega * resid11)

  c(theta10 = theta10_dir, theta00 = theta00_dir, theta11 = theta11_fold)
}

################################################################################
### One full estimate on a given (X, M, Y, A): stratified split, both
### directions, average -> NIE / NDE / TE.
################################################################################
n_split <- env_int("MC_NSPLIT", 1L)                        # repeated cross-fitting: average over this many random splits
estimate_effects <- function(X, M, Y, A) {                 # returns c(NIE, NDE, TE, theta11, theta10, theta00)
  nn <- length(A); X_bar <- colMeans(X)                    # sample size and target Xbar
  t_all <- which(A == 1); c_all <- which(A == 0)           # treated / control indices
  acc <- c(theta10 = 0, theta00 = 0, theta11 = 0)          # accumulator over split repetitions
  for (r in seq_len(n_split)) {                            # repeated-cross-fitting loop (DML: average over splits)
    t_perm <- sample(t_all); c_perm <- sample(c_all)       # shuffle each arm (fresh split each repetition)
    D1 <- c(t_perm[1:floor(length(t_perm) / 2)],           # D1 = half treated + half control (stratified)
            c_perm[1:floor(length(c_perm) / 2)])
    D2 <- setdiff(1:nn, D1)                                # D2 = the rest
    d12 <- run_direction(D1, D2, X, M, Y, A, X_bar)        # direction 1 -> 2
    d21 <- run_direction(D2, D1, X, M, Y, A, X_bar)        # direction 2 -> 1 (folds swapped)
    acc["theta10"] <- acc["theta10"] + 0.5 * (d12["theta10"] + d21["theta10"])  # accumulate direction-averaged theta10
    acc["theta00"] <- acc["theta00"] + 0.5 * (d12["theta00"] + d21["theta00"])
    acc["theta11"] <- acc["theta11"] + 0.5 * (d12["theta11"] + d21["theta11"])
  }
  theta10 <- unname(acc["theta10"] / n_split)              # average theta10 over the split repetitions
  theta00 <- unname(acc["theta00"] / n_split)
  theta11 <- unname(acc["theta11"] / n_split)

  NIE <- theta11 - theta10                                 # natural indirect effect
  NDE <- theta10 - theta00                                 # natural direct effect
  TE  <- NIE + NDE                                         # total effect
  c(NIE = unname(NIE), NDE = unname(NDE), TE = unname(TE),
    theta11 = unname(theta11), theta10 = unname(theta10), theta00 = unname(theta00))
}

################################################################################
### Load and prepare data
################################################################################
cat("Loading data ...\n")
A_Smoking <- readRDS(file.path(data_dir, "A_Smoking.rds"))                 # treatment: smoking status
X_all     <- data.matrix(readRDS(file.path(data_dir,
             Sys.getenv("MC_XFILE", "X_fourthorder_MICE.rds"))))  # covariates (phenotype); MC_XFILE overrides
Y_time    <- readRDS(file.path(data_dir, "Y_SurvivalTime.rds"))           # outcome: survival time
M_raw     <- readRDS(file.path(data_dir, "M_raw_cleaned.rds"))            # mediators: methylation
sd_thresh <- env_num("MC_SDT", 0.1)                                       # sd screening threshold (default 0.1)
crit      <- apply(M_raw, 2, sd) > sd_thresh                             # keep mediators with sd > threshold
M_all     <- M_raw[, crit]                                                # filtered mediator matrix
rm(M_raw); gc()                                                           # free the big raw matrix
cat(sprintf("sd threshold = %.2f\n", sd_thresh))                          # report the threshold used

nn <- length(A_Smoking)                                    # sample size
Y_all <- log(Y_time)                                       # log survival time (matches the original analysis)
cat(sprintf("n = %d | p = %d | q = %d | treated = %d | control = %d\n",
            nn, ncol(X_all), ncol(M_all), sum(A_Smoking == 1), sum(A_Smoking == 0)))

################################################################################
### Point estimate on the original sample
################################################################################
set.seed(seed)                                             # reproducible split
point_est <- estimate_effects(X_all, M_all, Y_all, A_Smoking)  # NIE / NDE / TE + potential-outcome means
cat("\nPoint estimate (original sample):\n"); print(round(point_est, 5))

################################################################################
### Bootstrap
################################################################################
cat(sprintf("\nRunning %d bootstrap replicates ...\n", n_boot))
cl <- make_mc_cluster(n_boot, seed)                        # FORK cluster + reproducible RNG
on.exit(parallel::stopCluster(cl), add = TRUE)             # always release the cluster

boot_mat <- foreach::foreach(bidx = 1:n_boot, .combine = rbind,          # bootstrap loop
                             .errorhandling = "remove") %dopar% {
  worker_single_thread_blas()                              # single-thread BLAS inside each worker
  boot_m      <- env_int("MC_BOOT_M", nn)                  # subsample size m (default n = classic bootstrap)
  boot_replace<- env_int("MC_REPLACE", 1L) == 1L           # 1=with replacement, 0=subsample w/o replacement
  idx <- sample(1:nn, size = boot_m, replace = boot_replace)
  out <- tryCatch(                                         # a failed replicate is dropped, not fatal
    estimate_effects(X_all[idx, , drop = FALSE], M_all[idx, , drop = FALSE],
                     Y_all[idx], A_Smoking[idx]),
    error = function(e) NULL)
  out                                                      # NULL replicates are removed by .errorhandling
}
boot_mat <- as.matrix(boot_mat)                            # ensure matrix form
cat(sprintf("Bootstrap done: %d / %d replicates succeeded\n", nrow(boot_mat), n_boot))

saveRDS(list(point = point_est, boot = boot_mat, n_boot = n_boot,        # persist everything for reuse
             K = c(K1, K2, K3), seed = seed),
        file.path(out_dir, "real_data_modified_results.rds"))

################################################################################
### Summary table (bootstrap mean + 2.5/97.5 quantile CI)
################################################################################
trim_frac <- env_num("MC_TRIM", 0.10)                      # symmetric trimming fraction for robust summary
summ <- function(nm, label) {                              # one summary row for estimand nm
  v <- boot_mat[, nm]; v <- v[is.finite(v)]                # bootstrap draws (finite)
  # robust: winsorize/trim the heavy tails that a few diverging resamples produce
  lo <- quantile(v, trim_frac, na.rm=TRUE); hi <- quantile(v, 1-trim_frac, na.rm=TRUE)
  vt <- v[v>=lo & v<=hi]                                   # trimmed draws
  qs <- quantile(vt, c(0.025, 0.975), na.rm = TRUE)        # quantile CI on trimmed draws
  data.frame(Effect = label, Estimate = round(point_est[nm], 4),
             BootMean = round(mean(vt, na.rm = TRUE), 4),
             BootMedian = round(median(v, na.rm=TRUE), 4),
             CI_low = round(qs[1], 4), CI_high = round(qs[2], 4),
             row.names = NULL)
}
tab <- rbind(summ("NIE", "Indirect Effect"),               # NIE row
             summ("NDE", "Direct Effect"),                 # NDE row
             summ("TE",  "Total Effect"))                  # TE row
cat("\n===== SUMMARY =====\n"); print(tab, row.names = FALSE)
write.csv(tab, file.path(out_dir, "real_data_modified_summary.csv"), row.names = FALSE)

################################################################################
### Histograms of the bootstrap distributions (NIE / NDE / TE)
################################################################################
pdf(file.path(out_dir, "real_data_modified_bootstrap_hist.pdf"), width = 11, height = 3.6)
par(mfrow = c(1, 3), mar = c(4, 4, 2, 1))                  # three panels side by side
for (nm in c("NIE", "NDE", "TE")) {                        # one histogram per effect
  v <- boot_mat[, nm]
  hist(v, breaks = 30, col = "#4a6fa5", border = "white",  # bootstrap histogram
       main = nm, xlab = expression(hat(theta)), ylab = "Count")
  abline(v = point_est[nm], col = "red", lwd = 2)          # point estimate marker
  abline(v = quantile(v, c(0.025, 0.975), na.rm = TRUE), col = "black", lty = 2)  # quantile CI markers
}
dev.off()
png(file.path(out_dir, "real_data_modified_bootstrap_hist.png"), width = 1650, height = 540, res = 150)
par(mfrow = c(1, 3), mar = c(4, 4, 2, 1))
for (nm in c("NIE", "NDE", "TE")) {
  v <- boot_mat[, nm]
  hist(v, breaks = 30, col = "#4a6fa5", border = "white", main = nm,
       xlab = expression(hat(theta)), ylab = "Count")
  abline(v = point_est[nm], col = "red", lwd = 2)
  abline(v = quantile(v, c(0.025, 0.975), na.rm = TRUE), col = "black", lty = 2)
}
dev.off()
cat("\nOutputs written to:", out_dir, "\n")
