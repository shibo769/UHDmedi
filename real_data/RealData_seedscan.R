################################################################################
### Multi-seed bootstrap at ONE fixed configuration.
###
### Purpose: the bootstrap mean at the paper's configuration sits far from the
### point estimate (NIE: est -0.256 vs trimmed mean +0.014). A single seed cannot
### tell us whether that gap is a property of the estimator or an artefact of one
### unlucky split, so we re-run the WHOLE analysis (point estimate + bootstrap)
### under several seeds and report the gap for each.
###
### The mediator matrix is loaded and sd-screened ONCE (~2 min, ~40 GB peak) and
### reused across seeds; restarting Rscript per seed would pay that cost every
### time.
###
### Config comes from the same env vars as RealData_stable_topj.R
### (MC_SDT / MC_XFILE / MC_LAMBDA / MC_K* / MC_NFOLDS / MC_TOPJ / MC_BOOT /
### MC_CORES / MC_TRIM); MC_SEEDS is the comma-separated seed list.
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
nfolds    <- env_int("MC_NFOLDS", 3L)
K1 <- env_num("MC_K1", 2.75); K2 <- env_num("MC_K2", 2.75); K3 <- env_num("MC_K3", 2.75)
alpha_inf <- 2/3
trim_frac <- env_num("MC_TRIM", 0.10)
TOPJ      <- env_int("MC_TOPJ", 0L)
n_split   <- env_int("MC_NSPLIT", 1L)
data_dir  <- file.path(repo, "modified", "real_data")
out_dir   <- file.path(repo, "modified", "real_data_outputs")
ensure_dir(out_dir)
out_csv   <- file.path(out_dir, env_chr("MC_OUT", "real_data_seedscan.csv"))
seeds     <- as.integer(strsplit(env_chr("MC_SEEDS", "555"), ",")[[1]])

################################################################################
### Estimator: identical to RealData_stable_topj.R (same fallback chain, same
### in-fold screening, same balancing calls) so the numbers are comparable.
################################################################################
balance_weights <- function(Xmat, target, K, alpha_exp) {
  m <- nrow(Xmat)
  approx.balance(M = Xmat, balance.target = target, zeta = 0.5, K = K,
                 allow.negative.weights = FALSE, optimizer = "pogs.new",
                 bound.gamma = TRUE, gamma.max = m^(-alpha_exp))
}
lam_choice <- env_chr("MC_LAMBDA", "lambda.1se")
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

run_direction <- function(idx_k, idx_kp, X, M, Y, A, X_bar) {
  ck  <- idx_k[A[idx_k] == 0]
  tkp <- idx_kp[A[idx_kp] == 1]
  keep <- seq_len(ncol(M))
  if (TOPJ > 0 && TOPJ < ncol(M)) {
    cr <- suppressWarnings(as.vector(abs(stats::cor(M[tkp, , drop = FALSE], Y[tkp]))))
    cr[!is.finite(cr)] <- 0
    keep <- sort(order(cr, decreasing = TRUE)[seq_len(TOPJ)])
  }
  p <- ncol(X); q <- length(keep)
  Xck <- X[ck, , drop = FALSE]; Mck <- M[ck, keep, drop = FALSE]; yck <- Y[ck]
  Xtkp <- X[tkp, , drop = FALSE]; Mtkp <- M[tkp, keep, drop = FALSE]
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

################################################################################
### Load data once, reuse across seeds
################################################################################
cat("Loading data ...\n")
A_Smoking <- readRDS(file.path(data_dir, "A_Smoking.rds"))
X_all     <- data.matrix(readRDS(file.path(data_dir, env_chr("MC_XFILE", "X_order3.rds"))))
Y_time    <- readRDS(file.path(data_dir, "Y_SurvivalTime.rds"))
M_raw     <- readRDS(file.path(data_dir, "M_raw_cleaned.rds"))
sd_thresh <- env_num("MC_SDT", 0.10)
crit      <- apply(M_raw, 2, sd) > sd_thresh
M_all     <- M_raw[, crit, drop = FALSE]
rm(M_raw); gc()
nn    <- length(A_Smoking)
Y_all <- log(Y_time)
cat(sprintf("sd = %.2f | n = %d | p = %d | q = %d | treated = %d | control = %d\n",
            sd_thresh, nn, ncol(X_all), ncol(M_all),
            sum(A_Smoking == 1), sum(A_Smoking == 0)))
cat(sprintf("lambda = %s | K = %g | nfolds = %d | topJ = %d | B = %d | trim = %.2f\n",
            lam_choice, K1, nfolds, TOPJ, n_boot, trim_frac))
cat("seeds:", paste(seeds, collapse = ", "), "\n\n")

################################################################################
### Loop over seeds; append each result to the csv as soon as it is done so a
### crash or a kill never loses the seeds already finished.
################################################################################
rows <- list()
for (sd_i in seq_along(seeds)) {
  sdv <- seeds[sd_i]
  set.seed(sdv)
  pe <- tryCatch(estimate_effects(X_all, M_all, Y_all, A_Smoking), error = function(e) NULL)
  if (is.null(pe)) { cat(sprintf("[seed %d] point estimate failed\n", sdv)); next }

  cl <- make_mc_cluster(n_boot, sdv)
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
  if (nrow(bm) < 30) { cat(sprintf("[seed %d] only %d replicates, skip\n", sdv, nrow(bm))); next }
  saveRDS(list(point = pe, boot = bm, seed = sdv, n_boot = n_boot),
          file.path(out_dir, sprintf("seedscan_%d.rds", sdv)))

  for (nm in c("NIE", "NDE", "TE")) {
    v <- bm[, nm]; v <- v[is.finite(v)]
    lo <- quantile(v, trim_frac, na.rm = TRUE); hi <- quantile(v, 1 - trim_frac, na.rm = TRUE)
    vt <- v[v >= lo & v <= hi]                     # trimmed core, as in RealData_stable.R
    q90 <- quantile(vt, c(0.05, 0.95), na.rm = TRUE)
    q95 <- quantile(vt, c(0.025, 0.975), na.rm = TRUE)
    rows[[length(rows) + 1]] <- data.frame(
      seed = sdv, Effect = nm,
      Estimate = round(unname(pe[nm]), 4),
      BootMean = round(mean(vt, na.rm = TRUE), 4),
      BootMedian = round(median(vt, na.rm = TRUE), 4),
      gap = round(mean(vt, na.rm = TRUE) - unname(pe[nm]), 4),   # bootstrap mean - point estimate
      BootSD = round(sd(vt, na.rm = TRUE), 4),
      RawMean = round(mean(v, na.rm = TRUE), 4),
      RawSD = round(sd(v, na.rm = TRUE), 4),
      pr_neg = round(mean(vt < 0, na.rm = TRUE), 3),
      CI90_low = round(q90[1], 4), CI90_high = round(q90[2], 4),
      CI95_low = round(q95[1], 4), CI95_high = round(q95[2], 4),
      nboot = nrow(bm), row.names = NULL)
  }
  cur <- do.call(rbind, rows)
  write.csv(cur, out_csv, row.names = FALSE)          # checkpoint after every seed
  r <- cur[cur$seed == sdv, ]
  cat(sprintf("[seed %8d] NIE est=%+.4f mean=%+.4f gap=%+.4f | NDE est=%+.4f mean=%+.4f gap=%+.4f | TE est=%+.4f mean=%+.4f gap=%+.4f\n",
              sdv,
              r$Estimate[1], r$BootMean[1], r$gap[1],
              r$Estimate[2], r$BootMean[2], r$gap[2],
              r$Estimate[3], r$BootMean[3], r$gap[3]))
  flush.console()
}

tab <- do.call(rbind, rows)
write.csv(tab, out_csv, row.names = FALSE)
cat("\n======== NIE rows, smallest |gap| first ========\n")
nie <- tab[tab$Effect == "NIE", ]
print(nie[order(abs(nie$gap)), ], row.names = FALSE)
cat("\nSaved:", out_csv, "\n")
