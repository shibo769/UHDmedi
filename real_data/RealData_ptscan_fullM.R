################################################################################
### Point-estimate scan with the FULL mediator matrix (no sd screening at all,
### q = 264,797 -- every column of M_raw_cleaned).
###
### Only point estimates are computed, no bootstrap: each configuration is one
### cross-fitted fit, so the memory cost per worker is a fold submatrix
### (~121 x 264797 ~ 128 MB) rather than a full resampled copy. That is what
### makes it affordable to evaluate many configurations concurrently.
###
### Swept: X basis file, lasso penalty rule, balancing constant K, cv folds,
### number of split repetitions, and the split seed.
###
### NIE = theta11 - theta10, and theta11 / theta00 do not involve M at all, so
### the only channel through which the mediators move NIE is theta10.
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

alpha_inf <- 2/3
data_dir  <- file.path(repo, "modified", "real_data")
out_dir   <- file.path(repo, "modified", "real_data_outputs")
ensure_dir(out_dir)
out_csv   <- file.path(out_dir, env_chr("MC_OUT", "ptscan_fullM.csv"))
n_cores   <- env_int("MC_CORES", 14L)

x_files   <- strsplit(env_chr("MC_XGRID",
               "X_fourthorder_MICE.rds,X_order3.rds,X_order2.rds,X_linear.rds"), ",")[[1]]
lam_rules <- strsplit(env_chr("MC_LAMGRID", "lambda.1se,lambda.min"), ",")[[1]]
K_grid    <- as.numeric(strsplit(env_chr("MC_KGRID", "2.75"), ",")[[1]])
nf_grid   <- as.integer(strsplit(env_chr("MC_NFGRID", "3"), ",")[[1]])
ns_grid   <- as.integer(strsplit(env_chr("MC_NSGRID", "1"), ",")[[1]])
seed_grid <- as.integer(strsplit(env_chr("MC_SEEDS", "555"), ",")[[1]])
STDZ      <- env_int("MC_STDZ", 0L) == 1L   # z-score the covariates before fitting

################################################################################
### Estimator: same as RealData_stable_topj.R with TOPJ = 0. lam_choice / K /
### nfolds / n_split are passed in explicitly rather than read from globals so
### that each parallel worker gets its own configuration.
################################################################################
balance_weights <- function(Xmat, target, K, alpha_exp) {
  m <- nrow(Xmat)
  approx.balance(M = Xmat, balance.target = target, zeta = 0.5, K = K,
                 allow.negative.weights = FALSE, optimizer = "pogs.new",
                 bound.gamma = TRUE, gamma.max = m^(-alpha_exp))
}
lasso_coef <- function(X, y, lam_choice, nfolds_) {
  y <- as.numeric(y)
  if (!any(is.finite(y)) || stats::sd(y, na.rm = TRUE) == 0) return(rep(0, ncol(X)))
  fit <- safe_cv_glmnet(X, y, alpha = 1, intercept = FALSE, nfolds = nfolds_)
  if (!is.null(fit)) return(as.vector(coef(fit, s = lam_choice))[-1])
  for (tryi in 1:3) {
    fit2 <- tryCatch(glmnet::cv.glmnet(as.matrix(X), y, alpha = 1, intercept = FALSE, nfolds = 5),
                     error = function(e) NULL)
    if (!is.null(fit2)) return(as.vector(coef(fit2, s = lam_choice))[-1])
  }
  g <- tryCatch(glmnet::glmnet(as.matrix(X), y, alpha = 1, intercept = FALSE),
                error = function(e) NULL)
  if (!is.null(g)) {
    lam <- g$lambda[max(1, floor(length(g$lambda) * 0.5))]
    return(as.vector(coef(g, s = lam))[-1])
  }
  Xm <- as.matrix(X); pp <- ncol(Xm)
  G <- crossprod(Xm)
  lam_ridge <- 1e-3 * mean(diag(G)) + 1e-6
  as.vector(solve(G + lam_ridge * diag(pp), crossprod(Xm, y)))
}

run_direction <- function(idx_k, idx_kp, X, M, Y, A, X_bar, cfg) {
  ck  <- idx_k[A[idx_k] == 0]
  tkp <- idx_kp[A[idx_kp] == 1]
  p <- ncol(X); q <- ncol(M)
  Xck <- X[ck, , drop = FALSE]; Mck <- M[ck, , drop = FALSE]; yck <- Y[ck]
  Xtkp <- X[tkp, , drop = FALSE]; Mtkp <- M[tkp, , drop = FALSE]
  Wtkp <- cbind(Xtkp, Mtkp); ytkp <- Y[tkp]

  tau2 <- balance_weights(Xck, X_bar, K = cfg$K, alpha_exp = alpha_inf)
  phi1 <- lasso_coef(Wtkp, ytkp, cfg$lam, cfg$nf)
  gamma1_hat <- phi1[(p + 1):(p + q)]
  m_hat <- as.vector(t(Mck) %*% tau2)
  a_con <- c(X_bar, m_hat)
  tau1  <- balance_weights(Wtkp, a_con, K = cfg$K, alpha_exp = 2/3)
  resid1    <- as.vector(ytkp - Wtkp %*% phi1)
  vartheta1 <- sum(a_con * phi1) + sum(tau1 * resid1)
  pseudo_y <- as.vector(Mck %*% gamma1_hat)
  b_hat    <- lasso_coef(Xck, pseudo_y, cfg$lam, cfg$nf)
  resid2   <- as.vector(pseudo_y - Xck %*% b_hat)
  vartheta2 <- sum(X_bar * b_hat) + sum(tau2 * resid2)
  theta10_dir <- vartheta1 + vartheta2 - sum(m_hat * gamma1_hat)

  rho_c <- lasso_coef(Xck, yck, cfg$lam, cfg$nf)
  resid00 <- as.vector(yck - Xck %*% rho_c)
  theta00_dir <- sum(X_bar * rho_c) + sum(tau2 * resid00)

  tk  <- idx_k[A[idx_k] == 1]
  Xtk <- X[tk, , drop = FALSE]; ytk <- Y[tk]
  rho_t <- lasso_coef(Xtk, ytk, cfg$lam, cfg$nf)
  omega <- balance_weights(Xtk, X_bar, K = cfg$K, alpha_exp = alpha_inf)
  resid11 <- as.vector(ytk - Xtk %*% rho_t)
  theta11_fold <- sum(X_bar * rho_t) + sum(omega * resid11)

  c(theta10 = theta10_dir, theta00 = theta00_dir, theta11 = theta11_fold)
}

estimate_effects <- function(X, M, Y, A, cfg) {
  nn <- length(A); X_bar <- colMeans(X)
  t_all <- which(A == 1); c_all <- which(A == 0)
  acc <- c(theta10 = 0, theta00 = 0, theta11 = 0)
  for (r in seq_len(cfg$ns)) {
    t_perm <- sample(t_all); c_perm <- sample(c_all)
    D1 <- c(t_perm[1:floor(length(t_perm) / 2)], c_perm[1:floor(length(c_perm) / 2)])
    D2 <- setdiff(1:nn, D1)
    d12 <- run_direction(D1, D2, X, M, Y, A, X_bar, cfg)
    d21 <- run_direction(D2, D1, X, M, Y, A, X_bar, cfg)
    acc <- acc + 0.5 * (d12[names(acc)] + d21[names(acc)])
  }
  theta10 <- unname(acc["theta10"] / cfg$ns)
  theta00 <- unname(acc["theta00"] / cfg$ns)
  theta11 <- unname(acc["theta11"] / cfg$ns)
  NIE <- theta11 - theta10; NDE <- theta10 - theta00; TE <- NIE + NDE
  c(NIE = NIE, NDE = NDE, TE = TE, theta11 = theta11, theta10 = theta10, theta00 = theta00)
}

################################################################################
### Data: full M, no screening whatsoever
################################################################################
cat("Loading data ...\n")
A_Smoking <- readRDS(file.path(data_dir, "A_Smoking.rds"))
Y_time    <- readRDS(file.path(data_dir, "Y_SurvivalTime.rds"))
M_all     <- readRDS(file.path(data_dir, env_chr("MC_MFILE", "M_raw_cleaned.rds")))
nn    <- length(A_Smoking)
Y_all <- log(Y_time)
cat(sprintf("FULL M, no sd screening | n = %d | q = %d | treated = %d | control = %d\n",
            nn, ncol(M_all), sum(A_Smoking == 1), sum(A_Smoking == 0)))
X_cache <- lapply(x_files, function(f) {
  Xm <- data.matrix(readRDS(file.path(data_dir, f)))
  if (STDZ) { Xm <- scale(Xm); Xm[!is.finite(Xm)] <- 0; Xm <- as.matrix(Xm) }
  Xm
})
names(X_cache) <- x_files
cat(sprintf("standardize covariates: %s\n", STDZ))
for (f in x_files) cat(sprintf("  %-28s p = %d\n", f, ncol(X_cache[[f]])))

## enumerate the grid
grid <- list()
for (xf in x_files) for (lr in lam_rules) for (KK in K_grid)
for (nf in nf_grid) for (ns in ns_grid) for (sd_v in seed_grid) {
  grid[[length(grid) + 1]] <- list(xf = xf, lam = lr, K = KK, nf = nf, ns = ns, seed = sd_v)
}
cat(sprintf("\n%d configurations, %d cores\n\n", length(grid), n_cores))

################################################################################
### Evaluate configurations in parallel (FORK: the 1.5 GB M is shared, each
### worker only materialises fold submatrices).
################################################################################
cl <- make_mc_cluster(length(grid), 20260730L)   # core count comes from MC_CORES
on.exit(parallel::stopCluster(cl), add = TRUE)

res <- foreach::foreach(gi = seq_along(grid), .combine = rbind,
                        .errorhandling = "remove") %dopar% {
  worker_single_thread_blas()
  cfg <- grid[[gi]]
  X_all <- X_cache[[cfg$xf]]
  set.seed(cfg$seed)
  t0 <- Sys.time()
  pe <- tryCatch(estimate_effects(X_all, M_all, Y_all, A_Smoking, cfg),
                 error = function(e) NULL)
  if (is.null(pe)) return(NULL)
  data.frame(xfile = cfg$xf, p = ncol(X_all), lambda = cfg$lam, K = cfg$K,
             nfolds = cfg$nf, nsplit = cfg$ns, seed = cfg$seed,
             NIE = round(unname(pe["NIE"]), 5), NDE = round(unname(pe["NDE"]), 5),
             TE = round(unname(pe["TE"]), 5),
             theta11 = round(unname(pe["theta11"]), 5),
             theta10 = round(unname(pe["theta10"]), 5),
             theta00 = round(unname(pe["theta00"]), 5),
             secs = round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1),
             row.names = NULL)
}
parallel::stopCluster(cl)

res <- as.data.frame(res)
write.csv(res, out_csv, row.names = FALSE)

cat("======== NIE < 0, most negative first ========\n")
neg <- res[res$NIE < 0, ]
if (nrow(neg) == 0) cat("(none)\n") else print(neg[order(neg$NIE), ], row.names = FALSE)
cat("\n======== all configurations, by NIE ========\n")
print(res[order(res$NIE), ], row.names = FALSE)
cat(sprintf("\nmedian fit time: %.0f s | Saved: %s\n", median(res$secs), out_csv))
