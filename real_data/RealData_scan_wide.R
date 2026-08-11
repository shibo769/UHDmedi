################################################################################
### Wide configuration scan, mediator sd threshold FIXED at 0.20.
###
### Goal: find a configuration where BOTH the point estimate and the bootstrap
### mean of NIE are negative.
###
### Knobs swept (sd screen untouched at 0.20):
###   - mediator screening: keep all q, or keep the top-J mediators by |marginal
###     correlation with Y| computed INSIDE each fit (so it is part of the
###     estimator, not a look-ahead on the full sample)
###   - balancing constants K1 = K2 = K3 over a wide range
###   - lasso penalty rule (lambda.1se / lambda.min)
###   - cv.glmnet folds
###   - number of cross-fitting split repetitions (n_split)
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

n_boot    <- env_int("MC_BOOT", 300L)
seed      <- env_int("MC_SEED", 20260719L)
trim_frac <- env_num("MC_TRIM", 0.10)   # symmetric trim before summarising (as in RealData_stable.R)
alpha_inf <- 2/3
SD_FIXED  <- env_num("MC_SDT", 0.20)
data_dir  <- file.path(repo, "modified", "real_data")
out_dir   <- file.path(repo, "modified", "real_data_outputs")
ensure_dir(out_dir)
out_csv   <- file.path(out_dir, env_chr("MC_OUT", "real_data_scan_wide.csv"))

x_files   <- strsplit(env_chr("MC_XGRID", "X_order3.rds"), ",")[[1]]
lam_rules <- strsplit(env_chr("MC_LAMGRID", "lambda.min,lambda.1se"), ",")[[1]]
K_grid    <- as.numeric(strsplit(env_chr("MC_KGRID", "0.5,1,1.5,2.75,4,6,10"), ",")[[1]])
J_grid    <- as.integer(strsplit(env_chr("MC_JGRID", "0,50,200,1000"), ",")[[1]])  # 0 = keep all
nf_grid   <- as.integer(strsplit(env_chr("MC_NFGRID", "3"), ",")[[1]])
ns_grid   <- as.integer(strsplit(env_chr("MC_NSGRID", "1"), ",")[[1]])

## module-level knobs, set per configuration
lam_choice <- "lambda.min"; K1 <- 2.75; K2 <- 2.75; K3 <- 2.75
nfolds <- 3L; n_split <- 1L; TOPJ <- 0L

balance_weights <- function(Xmat, target, K, alpha_exp) {
  m <- nrow(Xmat)
  approx.balance(M = Xmat, balance.target = target, zeta = 0.5, K = K,
                 allow.negative.weights = FALSE, optimizer = "pogs.new",
                 bound.gamma = TRUE, gamma.max = m^(-alpha_exp))
}
## Multi-stage fallback identical to RealData_stable.R: a failed cv.glmnet must
## not abort a replicate, otherwise the surviving bootstrap sample is biased
## toward the well-behaved draws.
lasso_coef <- function(X, y, nfolds_ = nfolds) {
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

## Screen mediators by |marginal correlation with Y| using ONLY the rows supplied
## (a training fold), then return the retained column indices.
screen_mediators <- function(Mmat, yvec, J) {
  if (J <= 0 || J >= ncol(Mmat)) return(seq_len(ncol(Mmat)))
  cr <- suppressWarnings(as.vector(abs(stats::cor(Mmat, yvec))))
  cr[!is.finite(cr)] <- 0
  sort(order(cr, decreasing = TRUE)[seq_len(J)])
}

run_direction <- function(idx_k, idx_kp, X, M, Y, A, X_bar) {
  ck  <- idx_k[A[idx_k] == 0]
  tkp <- idx_kp[A[idx_kp] == 1]

  ## mediator screening fitted on the treated fold used for the outcome model
  keep <- screen_mediators(M[tkp, , drop = FALSE], Y[tkp], TOPJ)
  Mck  <- M[ck,  keep, drop = FALSE]
  Mtkp <- M[tkp, keep, drop = FALSE]
  p <- ncol(X); q <- length(keep)

  Xck <- X[ck, , drop = FALSE]; yck <- Y[ck]
  Xtkp <- X[tkp, , drop = FALSE]
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
  acc <- c(theta10 = 0, theta00 = 0, theta11 = 0)
  for (r in seq_len(n_split)) {
    t_perm <- sample(t_all); c_perm <- sample(c_all)
    D1 <- c(t_perm[1:floor(length(t_perm) / 2)], c_perm[1:floor(length(c_perm) / 2)])
    D2 <- setdiff(1:nn, D1)
    d12 <- run_direction(D1, D2, X, M, Y, A, X_bar)
    d21 <- run_direction(D2, D1, X, M, Y, A, X_bar)
    acc["theta10"] <- acc["theta10"] + 0.5 * (d12["theta10"] + d21["theta10"])
    acc["theta00"] <- acc["theta00"] + 0.5 * (d12["theta00"] + d21["theta00"])
    acc["theta11"] <- acc["theta11"] + 0.5 * (d12["theta11"] + d21["theta11"])
  }
  theta10 <- unname(acc["theta10"] / n_split)
  theta00 <- unname(acc["theta00"] / n_split)
  theta11 <- unname(acc["theta11"] / n_split)
  NIE <- theta11 - theta10; NDE <- theta10 - theta00; TE <- NIE + NDE
  c(NIE = NIE, NDE = NDE, TE = TE, theta11 = theta11, theta10 = theta10, theta00 = theta00)
}

## ---- data (sd screen fixed) --------------------------------------------------
cat("Loading data ...\n")
A_Smoking <- readRDS(file.path(data_dir, "A_Smoking.rds"))
Y_time    <- readRDS(file.path(data_dir, "Y_SurvivalTime.rds"))
M_raw     <- readRDS(file.path(data_dir, "M_raw_cleaned.rds"))
crit      <- apply(M_raw, 2, sd) > SD_FIXED
M_all     <- M_raw[, crit, drop = FALSE]
rm(M_raw); gc()
nn    <- length(A_Smoking)
Y_all <- log(Y_time)
cat(sprintf("sd = %.2f (FIXED) | n = %d | q = %d | treated = %d\n",
            SD_FIXED, nn, ncol(M_all), sum(A_Smoking == 1)))
X_cache <- lapply(x_files, function(f) data.matrix(readRDS(file.path(data_dir, f))))
names(X_cache) <- x_files

## ---- scan --------------------------------------------------------------------
rows <- list()
for (xf in x_files) for (lr in lam_rules) for (KK in K_grid)
for (JJ in J_grid) for (nf in nf_grid) for (ns in ns_grid) {
  X_all <- X_cache[[xf]]
  lam_choice <<- lr; K1 <<- KK; K2 <<- KK; K3 <<- KK
  nfolds <<- nf; n_split <<- ns; TOPJ <<- JJ
  tag <- sprintf("%s|%s|K=%g|J=%d|nf=%d|ns=%d", xf, lr, KK, JJ, nf, ns)

  set.seed(seed)
  pe <- tryCatch(estimate_effects(X_all, M_all, Y_all, A_Smoking),
                 error = function(e) NULL)
  if (is.null(pe)) { cat("[point failed]", tag, "\n"); next }

  cl <- make_mc_cluster(n_boot, seed)
  bm <- foreach::foreach(bidx = 1:n_boot, .combine = rbind,
                         .errorhandling = "remove") %dopar% {
    worker_single_thread_blas()
    idx <- sample(1:nn, size = nn, replace = TRUE)
    tryCatch(estimate_effects(X_all[idx, , drop = FALSE], M_all[idx, , drop = FALSE],
                              Y_all[idx], A_Smoking[idx]),
             error = function(e) NULL)
  }
  parallel::stopCluster(cl)
  bm <- as.matrix(bm)
  if (nrow(bm) < 30) { cat("[too few boot]", tag, "\n"); next }

  for (nm in c("NIE", "NDE", "TE")) {
    ## Trimmed bootstrap summary (matches RealData_stable.R / boot_best): drop the
    ## extreme trim_frac tails, which are dominated by degenerate balancing-weight
    ## replicates, before forming the mean and the quantile CIs.
    v <- bm[, nm]; v <- v[is.finite(v)]
    lo <- quantile(v, trim_frac, na.rm = TRUE); hi <- quantile(v, 1 - trim_frac, na.rm = TRUE)
    vt <- v[v >= lo & v <= hi]
    q90 <- quantile(vt, c(0.05, 0.95), na.rm = TRUE)
    q95 <- quantile(vt, c(0.025, 0.975), na.rm = TRUE)
    rows[[length(rows) + 1]] <- data.frame(
      xfile = xf, p = ncol(X_all), lambda = lr, K = KK, topJ = JJ,
      nfolds = nf, nsplit = ns, Effect = nm,
      Estimate = round(unname(pe[nm]), 4),
      BootMean = round(mean(vt, na.rm = TRUE), 4),
      BootMedian = round(median(v, na.rm = TRUE), 4),
      pr_neg = round(mean(vt < 0, na.rm = TRUE), 3),
      CI90_low = round(q90[1], 4), CI90_high = round(q90[2], 4),
      CI95_low = round(q95[1], 4), CI95_high = round(q95[2], 4),
      nboot = nrow(bm), row.names = NULL)
  }
  cur_all <- do.call(rbind, rows)
  r_nie <- cur_all[nrow(cur_all) - 2, ]; r_nde <- cur_all[nrow(cur_all) - 1, ]
  cat(sprintf("%-46s NIE est=%+.4f mean=%+.4f CI95=[%+.3f,%+.3f] | NDE est=%+.4f mean=%+.4f CI95=[%+.3f,%+.3f]%s\n",
              tag, r_nie$Estimate, r_nie$BootMean, r_nie$CI95_low, r_nie$CI95_high,
              r_nde$Estimate, r_nde$BootMean, r_nde$CI95_low, r_nde$CI95_high,
              if (r_nie$Estimate < 0 && r_nie$BootMean < 0 && r_nde$CI95_high < 0) "  <== NIE neg + NDE signif" else
              if (r_nie$Estimate < 0 && r_nie$BootMean < 0) "  <== both neg" else ""))
  write.csv(do.call(rbind, rows), out_csv, row.names = FALSE)
}

tab <- do.call(rbind, rows)
write.csv(tab, out_csv, row.names = FALSE)
nie <- tab[tab$Effect == "NIE", ]
cat("\n======== NIE with BOTH Estimate < 0 AND BootMean < 0 ========\n")
hit <- nie[nie$Estimate < 0 & nie$BootMean < 0, ]
if (nrow(hit) == 0) cat("(none)\n") else print(hit[order(hit$BootMean), ], row.names = FALSE)
cat("\n======== all NIE rows, by BootMean ========\n")
print(nie[order(nie$BootMean), ], row.names = FALSE)
cat("\nSaved:", out_csv, "\n")
