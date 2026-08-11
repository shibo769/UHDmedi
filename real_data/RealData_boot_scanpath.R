################################################################################
### Bootstrap for a configuration whose point estimate came from the ptscan
### foreach loop.
###
### Why this script exists: RealData_stable_topj.R seeds the split with
### set.seed() in the MASTER process, whereas RealData_ptscan_fullM.R evaluates
### each configuration inside a doRNG-managed foreach WORKER. The same nominal
### seed therefore selects a different cross-fitting split in the two scripts
### (verified on seed 111: -0.0520 in the scan vs -0.0205 in the master path).
### To bootstrap the configuration the scan reported, we replay the scan's RNG
### path: same cluster size, same grid, refit only the target index.
###
### Trimming: MC_TRIMGRID reports several symmetric trimming fractions side by
### side. The primary number stays MC_TRIM (default 0.10, as elsewhere in the
### pipeline); the extra rows expose sensitivity rather than licensing a fraction
### chosen after seeing which one clears zero.
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

n_boot      <- env_int("MC_BOOT", 300L)
alpha_inf   <- 2/3
data_dir    <- file.path(repo, "modified", "real_data")
out_dir     <- file.path(repo, "modified", "real_data_outputs")
ensure_dir(out_dir)
tag         <- env_chr("MC_TAG", "p449_s333")
xfile       <- env_chr("MC_XFILE", "X_fourthorder_MICE.rds")
lam_choice  <- env_chr("MC_LAMBDA", "lambda.1se")
KK          <- env_num("MC_K", 2.75)
nfolds      <- env_int("MC_NFOLDS", 3L)
n_split     <- env_int("MC_NSPLIT", 1L)
target_seed <- env_int("MC_TARGET_SEED", 333L)
trim_grid   <- as.numeric(strsplit(env_chr("MC_TRIMGRID", "0,0.05,0.10,0.15,0.20,0.25"), ",")[[1]])
trim_main   <- env_num("MC_TRIM", 0.10)

balance_weights <- function(Xmat, target, K, alpha_exp) {
  m <- nrow(Xmat)
  approx.balance(M = Xmat, balance.target = target, zeta = 0.5, K = K,
                 allow.negative.weights = FALSE, optimizer = "pogs.new",
                 bound.gamma = TRUE, gamma.max = m^(-alpha_exp))
}
lasso_coef <- function(X, y) {
  y <- as.numeric(y)
  if (!any(is.finite(y)) || stats::sd(y, na.rm = TRUE) == 0) return(rep(0, ncol(X)))
  fit <- safe_cv_glmnet(X, y, alpha = 1, intercept = FALSE, nfolds = nfolds)
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
  as.vector(solve(G + (1e-3 * mean(diag(G)) + 1e-6) * diag(pp), crossprod(Xm, y)))
}

run_direction <- function(idx_k, idx_kp, X, M, Y, A, X_bar) {
  ck  <- idx_k[A[idx_k] == 0]
  tkp <- idx_kp[A[idx_kp] == 1]
  p <- ncol(X); q <- ncol(M)
  Xck <- X[ck, , drop = FALSE]; Mck <- M[ck, , drop = FALSE]; yck <- Y[ck]
  Xtkp <- X[tkp, , drop = FALSE]; Mtkp <- M[tkp, , drop = FALSE]
  Wtkp <- cbind(Xtkp, Mtkp); ytkp <- Y[tkp]

  tau2 <- balance_weights(Xck, X_bar, K = KK, alpha_exp = alpha_inf)
  phi1 <- lasso_coef(Wtkp, ytkp)
  gamma1_hat <- phi1[(p + 1):(p + q)]
  m_hat <- as.vector(t(Mck) %*% tau2)
  a_con <- c(X_bar, m_hat)
  tau1  <- balance_weights(Wtkp, a_con, K = KK, alpha_exp = 2/3)
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
  omega <- balance_weights(Xtk, X_bar, K = KK, alpha_exp = alpha_inf)
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
    acc <- acc + 0.5 * (d12[names(acc)] + d21[names(acc)])
  }
  theta10 <- unname(acc["theta10"] / n_split); theta00 <- unname(acc["theta00"] / n_split)
  theta11 <- unname(acc["theta11"] / n_split)
  NIE <- theta11 - theta10; NDE <- theta10 - theta00
  c(NIE = NIE, NDE = NDE, TE = NIE + NDE,
    theta11 = theta11, theta10 = theta10, theta00 = theta00)
}

cat("Loading data ...\n")
A_Smoking <- readRDS(file.path(data_dir, "A_Smoking.rds"))
Y_time    <- readRDS(file.path(data_dir, "Y_SurvivalTime.rds"))
M_all     <- readRDS(file.path(data_dir, "M_raw_cleaned.rds"))
X_all     <- data.matrix(readRDS(file.path(data_dir, xfile)))
nn        <- length(A_Smoking)
Y_all     <- log(Y_time)
cat(sprintf("FULL M | n = %d | p = %d | q = %d | X = %s | lambda = %s | K = %g\n",
            nn, ncol(X_all), ncol(M_all), xfile, lam_choice, KK))

scan_xgrid <- strsplit(env_chr("MC_SCAN_XGRID",
                "X_fourthorder_MICE.rds,X_order3.rds,X_order2.rds,X_linear.rds"), ",")[[1]]
scan_lam   <- strsplit(env_chr("MC_SCAN_LAMGRID", "lambda.1se,lambda.min"), ",")[[1]]
scan_seeds <- as.integer(strsplit(env_chr("MC_SCAN_SEEDS", "555,111,333"), ",")[[1]])
grid <- list()
for (xf in scan_xgrid) for (lr in scan_lam) for (sv in scan_seeds)
  grid[[length(grid) + 1]] <- list(xf = xf, lam = lr, seed = sv)
gi_target <- which(vapply(grid, function(g)
  g$xf == xfile && g$lam == lam_choice && g$seed == target_seed, logical(1)))[1]
cat(sprintf("target config = grid index %d of %d\n", gi_target, length(grid)))

## Resume: if a previous run of this tag already produced the point estimate and
## some bootstrap draws, reuse them instead of paying for step 1 again.
RESUME    <- env_int("MC_RESUME", 0L) == 1L
resume_rds <- file.path(out_dir, sprintf("boot_%s.rds", tag))
prev <- if (RESUME && file.exists(resume_rds)) readRDS(resume_rds) else NULL
if (!is.null(prev)) {
  point_est <- prev$point
  cat(sprintf("Resuming: reusing point estimate and %d existing replicates from %s\n",
              nrow(prev$boot), basename(resume_rds)))
  print(round(point_est, 5))
} else {
cat("\nStep 1: replaying scan RNG path for the point estimate ...\n")
Sys.setenv(MC_CORES = env_chr("MC_SCAN_CORES", "16"))
cl <- make_mc_cluster(length(grid), 20260730L)
pt_all <- foreach::foreach(gi = seq_along(grid), .combine = rbind,
                           .errorhandling = "remove") %dopar% {
  worker_single_thread_blas()
  if (gi != gi_target) return(NULL)
  set.seed(grid[[gi]]$seed)
  pe <- tryCatch(estimate_effects(X_all, M_all, Y_all, A_Smoking), error = function(e) NULL)
  if (is.null(pe)) return(NULL)
  as.data.frame(as.list(pe))
}
parallel::stopCluster(cl)
point_est <- unlist(as.data.frame(pt_all)[1, ])
cat("Point estimate (scan RNG path):\n"); print(round(point_est, 5))
saveRDS(point_est, file.path(out_dir, sprintf("point_%s.rds", tag)))
}

################################################################################
### Trimming table for a set of bootstrap draws. Factored out so it can be
### applied to a partial batch mid-run as well as to the full set at the end.
################################################################################
make_trim_table <- function(bm, point_est, trim_grid) {
  rows <- list()
  for (tf in trim_grid) for (nm in c("NIE", "NDE", "TE")) {
    v <- bm[, nm]; v <- v[is.finite(v)]
    if (tf > 0) {
      lo <- quantile(v, tf, na.rm = TRUE); hi <- quantile(v, 1 - tf, na.rm = TRUE)
      vt <- v[v >= lo & v <= hi]
    } else vt <- v
    q90 <- quantile(vt, c(0.05, 0.95), na.rm = TRUE)
    q95 <- quantile(vt, c(0.025, 0.975), na.rm = TRUE)
    rows[[length(rows) + 1]] <- data.frame(
      trim = tf, Effect = nm, kept = length(vt),
      Estimate = round(unname(point_est[nm]), 4),
      BootMean = round(mean(vt, na.rm = TRUE), 4),
      BootMedian = round(median(vt, na.rm = TRUE), 4),
      gap = round(mean(vt, na.rm = TRUE) - unname(point_est[nm]), 4),
      BootSD = round(sd(vt, na.rm = TRUE), 4),
      pr_neg = round(mean(vt < 0, na.rm = TRUE), 3),
      CI90_low = round(q90[1], 4), CI90_high = round(q90[2], 4),
      CI95_low = round(q95[1], 4), CI95_high = round(q95[2], 4),
      excl0_90 = (q90[1] < 0) == (q90[2] < 0),
      excl0_95 = (q95[1] < 0) == (q95[2] < 0),
      row.names = NULL)
  }
  do.call(rbind, rows)
}

cat(sprintf("\nStep 2: %d bootstrap replicates ...\n", n_boot))
Sys.setenv(MC_CORES = env_chr("MC_BOOT_CORES", "20"))
chunk_sz <- env_int("MC_CHUNK", 50L)          # report a full table every chunk
boot_mat <- if (!is.null(prev)) prev$boot else NULL
n_done  <- if (!is.null(prev)) nrow(prev$boot) else 0L   # count successes, not attempts
chunk_i <- 0L
max_chunks <- env_int("MC_MAXCHUNKS", 40L)               # backstop against endless retries
while (n_done < n_boot && chunk_i < max_chunks) {
  chunk_i <- chunk_i + 1L
  this_n  <- min(chunk_sz, n_boot - n_done)
  ## A fresh cluster per chunk, seeded distinctly, so the draws stay independent
  ## across chunks while each chunk's results reach the master immediately.
  ## The whole chunk is guarded: a cluster-level failure (e.g. a worker dying and
  ## breaking the socket) discards that chunk's draws and the loop simply tries
  ## another one, rather than killing the run and losing the draws already banked.
  bm_chunk <- tryCatch({
    cl <- make_mc_cluster(this_n, target_seed + 1000L * chunk_i)
    on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
    out <- foreach::foreach(bidx = 1:this_n, .combine = rbind,
                            .errorhandling = "remove") %dopar% {
      worker_single_thread_blas()
      idx <- sample(1:nn, size = nn, replace = TRUE)
      tryCatch(estimate_effects(X_all[idx, , drop = FALSE], M_all[idx, , drop = FALSE],
                                Y_all[idx], A_Smoking[idx]),
               error = function(e) NULL)
    }
    try(parallel::stopCluster(cl), silent = TRUE)
    out
  }, error = function(e) {
    cat(sprintf("[chunk %d] discarded: %s\n", chunk_i, conditionMessage(e)))
    NULL
  })
  if (is.null(bm_chunk) || !length(bm_chunk)) { gc(); next }
  bm_chunk <- as.matrix(bm_chunk)
  boot_mat <- if (is.null(boot_mat)) bm_chunk else rbind(boot_mat, bm_chunk)
  n_done   <- nrow(boot_mat)
  cat(sprintf("\n########## chunk %d done: %d successful replicates banked (target %d) ##########\n",
              chunk_i, n_done, n_boot))
  gc()
  saveRDS(list(point = point_est, boot = boot_mat, n_done = n_done, xfile = xfile,
               lambda = lam_choice, K = KK, target_seed = target_seed),
          file.path(out_dir, sprintf("boot_%s.rds", tag)))
  if (nrow(boot_mat) >= 10) {
    tab_i <- make_trim_table(boot_mat, point_est, trim_grid)
    write.csv(tab_i, file.path(out_dir, sprintf("boot_%s_trimgrid.csv", tag)), row.names = FALSE)
    print(tab_i, row.names = FALSE)
    cat(sprintf("raw (untrimmed) sd: NIE %.3f | NDE %.3f | TE %.3f\n",
        sd(boot_mat[,"NIE"], na.rm=TRUE), sd(boot_mat[,"NDE"], na.rm=TRUE),
        sd(boot_mat[,"TE"], na.rm=TRUE)))
  }
  flush.console()
}
cat(sprintf("\nBootstrap done: %d / %d replicates succeeded\n", nrow(boot_mat), n_boot))

tab <- make_trim_table(boot_mat, point_est, trim_grid)
write.csv(tab, file.path(out_dir, sprintf("boot_%s_trimgrid.csv", tag)), row.names = FALSE)
cat("\n===== ALL TRIMMING LEVELS =====\n"); print(tab, row.names = FALSE)
cat(sprintf("\n===== PRIMARY (trim = %.2f) =====\n", trim_main))
print(tab[abs(tab$trim - trim_main) < 1e-9, ], row.names = FALSE)
cat(sprintf("\nraw (untrimmed) sd: NIE %.3f | NDE %.3f | TE %.3f\n",
    sd(boot_mat[,"NIE"], na.rm=TRUE), sd(boot_mat[,"NDE"], na.rm=TRUE),
    sd(boot_mat[,"TE"], na.rm=TRUE)))
cat("Saved:", file.path(out_dir, sprintf("boot_%s_trimgrid.csv", tag)), "\n")
