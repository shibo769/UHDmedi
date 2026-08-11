############################################################
# High-dimensional mediation: no interaction
# Our debiased method (TE/NDE approach)
#
# DGP options (env var DGP_NAME):
#   base | increasing_sparsity | minimal_signal
#
# Output (per DGP, per sigma):
#   1. per-iteration CSV for each (p=q, n)
#   2. full RDS object
#   3. MC CI summary (CSV + RDS)
#   4. analytical-variance CI summary (CSV + RDS)
#   5. histogram + QQ plots (PNG + PDF)
#
# Env-controlled. Examples:
#   MC_ITER=2 MC_PQ=50 MC_N=250 MC_CORES=4 DGP_NAME=base Rscript NoInteraction.R
############################################################

rm(list = ls())

################################################################################
### Thread guard: single-thread BLAS/OpenMP BEFORE loading glmnet / linalg
################################################################################

Sys.setenv(
  OMP_NUM_THREADS        = "1",
  OPENBLAS_NUM_THREADS   = "1",
  MKL_NUM_THREADS        = "1",
  BLIS_NUM_THREADS       = "1",
  VECLIB_MAXIMUM_THREADS = "1"
)

################################################################################
### Paths (portable)
################################################################################

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  }
  normalizePath(getwd())
}

base_dir <- normalizePath(Sys.getenv("BASE_DIR", unset = get_script_dir()), mustWork = TRUE)

## Shared plumbing at repo root (parent of Comparison/).
SHARED <- file.path(base_dir, "..", "shared")   # bundled deps
source(file.path(SHARED, "mc_common.R"))
## Residual-balancing estimator functions live in Main_simu/utilis.R.
source(file.path(SHARED, "utilis.R"))

################################################################################
### Package handling (check only)
################################################################################

required_pkgs <- c(
  "foreach", "doParallel", "doRNG", "pracma", "MASS",
  "glmnet", "quadprog", "ggplot2", "gridExtra", "dplyr"
)
require_packages(required_pkgs)

suppressPackageStartupMessages({
  library(foreach)
  library(doParallel)
  library(doRNG)
  library(pracma)
  library(MASS)
  library(glmnet)
  library(quadprog)
  library(ggplot2)
  library(gridExtra)
  library(dplyr)
})

has_pogs <- requireNamespace("pogs", quietly = TRUE)
allow_quadprog_fallback <- env_lgl("MC_ALLOW_QUADPROG_FALLBACK", FALSE)

if (!has_pogs) {
  msg <- paste0(
    "The 'pogs' package is not installed, so optimizer = 'pogs.new' will fall ",
    "back to quadprog, which is a DIFFERENT estimator.\n",
    "Install it with: Rscript setup_packages.R"
  )
  if (allow_quadprog_fallback) {
    warning("*** ", msg, "\n*** Proceeding with quadprog fallback because ",
            "MC_ALLOW_QUADPROG_FALLBACK=TRUE (intended for quick tests only).",
            immediate. = TRUE)
  } else {
    stop(msg, "\nTo run a quick smoke test anyway, set ",
         "MC_ALLOW_QUADPROG_FALLBACK=TRUE.", call. = FALSE)
  }
}

################################################################################
### Parameters (env-controlled)
################################################################################

Iter      <- env_int("MC_ITER", 500L)
pq_vec    <- env_int_vec("MC_PQ", "50,400,750,1250")
n_vec     <- env_int_vec("MC_N", "250,500,750,1000,1250")
sigma_eps <- env_num("MC_SIGMA", 1)
DGP_NAME  <- env_chr("DGP_NAME", "base")
seed      <- env_int("MC_SEED", 1998L)
nfolds    <- env_int("MC_NFOLDS", 3L)
run_plots <- env_lgl("RUN_PLOTS", TRUE)

options(MC_NFOLDS = nfolds)

s_base   <- 3
s_B_base <- 3
mu       <- 0.5
K1       <- 2.75
K2       <- 2.75

if (!(DGP_NAME %in% c("base", "increasing_sparsity", "minimal_signal"))) {
  stop("DGP_NAME must be one of: base, increasing_sparsity, minimal_signal")
}

set.seed(seed)

################################################################################
### Output folder (no-interaction debiasing outputs live under Comparison/)
################################################################################

sigma_tag <- make_sigma_tag(sigma_eps)
file_tag  <- paste0("debiasing_nointeraction_", DGP_NAME, "_", sigma_tag)

output_root <- file.path(base_dir, "simulation_outputs_debiasing_no_interaction")
output_dir  <- file.path(output_root, DGP_NAME, sigma_tag)
ensure_dir(output_dir)

cat("No-interaction debiasing DGP version:", DGP_NAME, "\n")
cat("Iter:", Iter, "| pq:", paste(pq_vec, collapse = ","),
    "| n:", paste(n_vec, collapse = ","), "\n")
cat("Sigma:", sigma_eps, "| seed:", seed, "| nfolds:", nfolds, "\n")
cat("Outputs will be saved to:", output_dir, "\n")

################################################################################
### Cluster
################################################################################

my.cluster <- make_mc_cluster(Iter, seed)
on.exit({
  try(parallel::stopCluster(my.cluster), silent = TRUE)
}, add = TRUE)

################################################################################
### DGP helpers
################################################################################

s_map <- list(
  "50"   = list(s = 5,  s_B = 5),
  "400"  = list(s = 10, s_B = 10),
  "750"  = list(s = 10, s_B = 10),
  "1250" = list(s = 20, s_B = 20)
)

get_s_values <- function(pq, dgp_name) {

  if (dgp_name == "base") {
    return(list(s = s_base, s_B = s_B_base))
  }

  return(s_map[[as.character(pq)]])
}

draw_signal_nointeraction <- function(m, dgp_name, ref_n = max(n_vec)) {

  if (dgp_name == "minimal_signal") {

    m_weak   <- max(1, floor(m / 2))
    m_strong <- m - m_weak

    return(c(
      runif(m_strong, 0.05, 0.2),
      runif(m_weak, 0.8, 1.2) / sqrt(ref_n)
    ))

  } else {

    return(runif(m, 0.05, 0.2))

  }
}

################### Small helpers for residualBalance.mean output ###################

get_mu <- function(out) {
  if ("mu.hat" %in% names(out)) {
    as.numeric(out["mu.hat"])
  } else {
    as.numeric(out[1])
  }
}

get_var <- function(out) {
  if ("var.hat" %in% names(out)) {
    as.numeric(out["var.hat"])
  } else if (length(out) >= 2) {
    as.numeric(out[2])
  } else {
    NA_real_
  }
}

get_cov_rb <- function(out1, out2) {

  gamma1 <- attr(out1, "gamma")
  gamma2 <- attr(out2, "gamma")
  res1   <- attr(out1, "residuals")
  res2   <- attr(out2, "residuals")

  if (is.null(gamma1) || is.null(gamma2) || is.null(res1) || is.null(res2)) {
    return(0)
  }

  if (length(gamma1) != length(gamma2) ||
      length(res1) != length(res2) ||
      length(gamma1) != length(res1)) {
    return(0)
  }

  out <- sum(as.numeric(gamma1) * as.numeric(gamma2) * as.numeric(res1) * as.numeric(res2))

  if (!is.finite(out)) out <- 0

  out
}

################################################################################
### Result column template (keeps failed replications rbind-compatible)
################################################################################

result_cols <- c(
  "iter", "pq", "n",
  "TE_true", "NDE_true", "NIE_true",
  "NIE_hat_cf", "NIE_hat_noncf",
  "bias_NIE_cf", "bias_NIE_noncf",
  "se_NIE_cf", "se_NIE_cf_nocov", "se_NIE_noncf",
  "var_NIE_cf", "var_NIE_cf_nocov", "var_NIE_noncf",
  "var_TE_1", "var_NDE_1", "var_TE_2", "var_NDE_2",
  "cov_D1", "cov_D2",
  "failed"
)

################### Simulation ###################

res <- list()

for (pq in pq_vec){

  p <- pq
  q <- pq

  s_info <- get_s_values(pq, DGP_NAME)
  s   <- s_info$s
  s_B <- s_info$s_B

  cat("\n============================\n")
  cat("Starting pq =", pq, "| s =", s, "| s_B =", s_B, "\n")
  cat("============================\n")

  ############################################################
  ### True beta
  ############################################################

  beta <- rep(0, p)
  beta_index <- sample(1:p, s, FALSE)
  beta[beta_index] <- draw_signal_nointeraction(s, DGP_NAME)

  ############################################################
  ### True gamma and delta1 with partial overlap
  ############################################################

  idx_gamma  <- sample(1:q, s, FALSE)
  idx_shared <- sample(idx_gamma, floor(s/2), FALSE)
  idx_delta  <- c(
    idx_shared,
    sample(setdiff(1:q, idx_gamma), s - length(idx_shared), FALSE)
  )

  gamma <- rep(0, q)
  gamma[idx_gamma] <- draw_signal_nointeraction(s, DGP_NAME)

  delta1 <- rep(0, q)
  delta1[idx_delta] <- draw_signal_nointeraction(s, DGP_NAME)

  alpha1 <- 0.15

  TE_true  <- alpha1 + sum(delta1 * gamma)
  NDE_true <- alpha1
  NIE_true <- sum(delta1 * gamma)

  ############################################################
  ### True B matrix
  ############################################################

  B <- matrix(0, q, p)

  for (iB in 1:q) {
    B[iB, sample(1:p, s_B, FALSE)] <- runif(s_B, 0.05, 0.2)
  }

  for (nn in n_vec){

    ############################################################
    ### X generation
    ############################################################

    mu_vec <- rep(mu, p)

    U_mat <- randortho(p)

    Sigma_X <- U_mat %*%
      diag(runif(p, 0.1, 0.5)) %*%
      t(U_mat)

    X <- mvrnorm(nn, mu_vec, Sigma_X)

    ############################################################
    ### A generation
    ############################################################

    alpha_A <- rep(0, p)
    alpha_A[sample(1:p, s, FALSE)] <- runif(s, 0, 2)

    linpred <- as.vector(X %*% alpha_A)
    linpred <- linpred - mean(linpred)
    sd_linpred <- sd(linpred)

    if (is.finite(sd_linpred) && sd_linpred > 1e-8) {
      linpred <- linpred / sd_linpred
    }

    prob_A <- plogis(linpred)
    A <- as.numeric(rbinom(nn, 1, prob_A))

    if (length(unique(A)) < 2) {
      warning("Skipping pq = ", pq, ", n = ", nn, " because A has no variation.")
      next
    }

    ############################################################
    ### Monte Carlo
    ############################################################

    err_vec <- foreach(
      i = 1:Iter,
      .combine = "rbind",
      .packages = c("glmnet", "pracma", "MASS", "quadprog"),
      .export = c(
        "residualBalance.mean",
        "approx.balance",
        "approx.balance.pogs.new",
        "approx.balance.pogs",
        "approx.balance.pogs.dual",
        "approx.balance.quadprog",
        "approx.balance.mosek.dual",
        "get_mu",
        "get_var",
        "get_cov_rb",
        "safe_cv_glmnet",
        "cv_nfolds",
        "na_result_row",
        "result_cols",
        "nfolds",
        "worker_single_thread_blas"
      ),
      .options.multicore = mc_foreach_options_multicore,
      .options.snow = mc_foreach_options_snow
    ) %dopar% {

      require(glmnet)
      require(pracma)
      require(MASS)
      worker_single_thread_blas()

      one_rep <- tryCatch({

        ############################################################
        ### Generate data
        ############################################################

        lambda_U <- runif(q, 0.1, 0.5)

        R_U <- randortho(q)

        Sigma_U <- R_U %*%
          diag(lambda_U) %*%
          t(R_U)

        U <- mvrnorm(nn, rep(0, q), Sigma_U)

        A_mat <- matrix(A, ncol = 1)

        M <- A_mat %*% t(delta1) +
          X %*% t(B) +
          U

        Y <- as.vector(
          alpha1 * A +
            X %*% beta +
            M %*% gamma +
            rnorm(nn, 0, sigma_eps)
        )

        ############################################################
        ### Original non-cross-fitted estimator
        ############################################################

        W1 <- cbind(A, X)
        contrast_a <- c(1, rep(0, p))

        out_TE <- residualBalance.mean(
          W1,
          Y,
          balance.target = contrast_a,
          zeta = 0.5,
          K = K1,
          fit.method = "elnet",
          alpha = 1,
          optimizer = "pogs.new"
        )

        TE_hat <- get_mu(out_TE)
        var_TE <- get_var(out_TE)

        W2 <- cbind(A, X, M)
        contrast_b <- c(1, rep(0, p + q))

        out_NDE <- residualBalance.mean(
          W2,
          Y,
          balance.target = contrast_b,
          zeta = 0.5,
          K = K2,
          fit.method = "elnet",
          alpha = 1,
          optimizer = "pogs.new"
        )

        NDE_hat <- get_mu(out_NDE)
        var_NDE <- get_var(out_NDE)

        cov_TE_NDE <- get_cov_rb(out_TE, out_NDE)
        var_NIE_noncf <- var_TE + var_NDE - 2 * cov_TE_NDE
        var_NIE_noncf <- pmax(var_NIE_noncf, 0)
        se_NIE_noncf <- sqrt(var_NIE_noncf)

        NIE_hat_noncf <- TE_hat - NDE_hat

        ############################################################
        ### Cross fitting split
        ############################################################

        all_index <- sample(1:nn, nn, replace = FALSE)

        D1 <- all_index[1:floor(nn / 2)]
        D2 <- all_index[(floor(nn / 2) + 1):nn]

        ############################################################
        ### First direction:
        ### TE using D1, NDE using D2
        ############################################################

        A_D1 <- A[D1]
        X_D1 <- X[D1, , drop = FALSE]
        M_D1 <- M[D1, , drop = FALSE]
        Y_D1 <- Y[D1]

        W1_D1 <- cbind(A_D1, X_D1)
        W2_D1 <- cbind(A_D1, X_D1, M_D1)

        A_D2 <- A[D2]
        X_D2 <- X[D2, , drop = FALSE]
        M_D2 <- M[D2, , drop = FALSE]
        Y_D2 <- Y[D2]

        W1_D2 <- cbind(A_D2, X_D2)
        W2_D2 <- cbind(A_D2, X_D2, M_D2)

        contrast_a_1 <- c(1, rep(0, p))
        contrast_b_1 <- c(1, rep(0, p + q))

        out_TE_1 <- residualBalance.mean(
          W1_D1,
          Y_D1,
          balance.target = contrast_a_1,
          zeta = 0.5,
          K = K1,
          fit.method = "elnet",
          alpha = 1,
          optimizer = "pogs.new"
        )

        out_NDE_1 <- residualBalance.mean(
          W2_D2,
          Y_D2,
          balance.target = contrast_b_1,
          zeta = 0.5,
          K = K2,
          fit.method = "elnet",
          alpha = 1,
          optimizer = "pogs.new"
        )

        TE_hat_1  <- get_mu(out_TE_1)
        NDE_hat_1 <- get_mu(out_NDE_1)
        NIE_hat_1 <- TE_hat_1 - NDE_hat_1

        ############################################################
        ### Second direction:
        ### TE using D2, NDE using D1
        ############################################################

        contrast_a_2 <- c(1, rep(0, p))
        contrast_b_2 <- c(1, rep(0, p + q))

        out_TE_2 <- residualBalance.mean(
          W1_D2,
          Y_D2,
          balance.target = contrast_a_2,
          zeta = 0.5,
          K = K1,
          fit.method = "elnet",
          alpha = 1,
          optimizer = "pogs.new"
        )

        out_NDE_2 <- residualBalance.mean(
          W2_D1,
          Y_D1,
          balance.target = contrast_b_2,
          zeta = 0.5,
          K = K2,
          fit.method = "elnet",
          alpha = 1,
          optimizer = "pogs.new"
        )

        TE_hat_2  <- get_mu(out_TE_2)
        NDE_hat_2 <- get_mu(out_NDE_2)
        NIE_hat_2 <- TE_hat_2 - NDE_hat_2

        ############################################################
        ### Final cross-fitted estimator
        ############################################################

        NIE_hat_cf <- 0.5 * (NIE_hat_1 + NIE_hat_2)

        ############################################################
        ### Analytical variance for cross-fitted estimator
        ############################################################

        var_TE_1  <- get_var(out_TE_1)
        var_NDE_2 <- get_var(out_NDE_2)
        cov_D1    <- get_cov_rb(out_TE_1, out_NDE_2)

        var_TE_2  <- get_var(out_TE_2)
        var_NDE_1 <- get_var(out_NDE_1)
        cov_D2    <- get_cov_rb(out_TE_2, out_NDE_1)

        ## TE and NDE are estimated on COMPLEMENTARY cross-fitting folds
        ## (TE on D_k, NDE on D_{k'}), so they are approximately independent and
        ## the cross-fold covariance is ~0. The empirical get_cov_rb() estimate
        ## multiplies weights/residuals from two DIFFERENT regressions (cbind(A,X)
        ## vs cbind(A,X,M)) and does not form a valid covariance; subtracting it
        ## systematically deflates the variance (SE underestimated ~2x). We drop
        ## it and keep cov_D1/cov_D2 only as diagnostics.
        var_diff_D1 <- var_TE_1 + var_NDE_2
        var_diff_D2 <- var_TE_2 + var_NDE_1

        var_diff_D1 <- pmax(var_diff_D1, 0)
        var_diff_D2 <- pmax(var_diff_D2, 0)

        var_NIE_cf <- (var_diff_D1 + var_diff_D2) / 4
        se_NIE_cf  <- sqrt(var_NIE_cf)

        ### Also save a no-covariance version for diagnostics
        var_NIE_cf_nocov <- (var_TE_1 + var_NDE_2 + var_TE_2 + var_NDE_1) / 4
        se_NIE_cf_nocov  <- sqrt(pmax(var_NIE_cf_nocov, 0))

        ############################################################
        ### Bias
        ############################################################

        bias_NIE_cf    <- NIE_hat_cf - NIE_true
        bias_NIE_noncf <- NIE_hat_noncf - NIE_true

        ############################################################
        ### Return
        ############################################################

        c(
          iter = i,
          pq = pq,
          n = nn,

          TE_true = TE_true,
          NDE_true = NDE_true,
          NIE_true = NIE_true,

          NIE_hat_cf = as.numeric(NIE_hat_cf),
          NIE_hat_noncf = as.numeric(NIE_hat_noncf),

          bias_NIE_cf = as.numeric(bias_NIE_cf),
          bias_NIE_noncf = as.numeric(bias_NIE_noncf),

          se_NIE_cf = as.numeric(se_NIE_cf),
          se_NIE_cf_nocov = as.numeric(se_NIE_cf_nocov),
          se_NIE_noncf = as.numeric(se_NIE_noncf),

          var_NIE_cf = as.numeric(var_NIE_cf),
          var_NIE_cf_nocov = as.numeric(var_NIE_cf_nocov),
          var_NIE_noncf = as.numeric(var_NIE_noncf),

          var_TE_1 = as.numeric(var_TE_1),
          var_NDE_1 = as.numeric(var_NDE_1),
          var_TE_2 = as.numeric(var_TE_2),
          var_NDE_2 = as.numeric(var_NDE_2),
          cov_D1 = as.numeric(cov_D1),
          cov_D2 = as.numeric(cov_D2),
          failed = 0
        )
      }, error = function(e) {
        ## Failed replication: keep identifying fields + true values, NA the rest.
        row <- na_result_row(result_cols)
        row["iter"] <- i
        row["pq"] <- pq
        row["n"] <- nn
        row["TE_true"] <- TE_true
        row["NDE_true"] <- NDE_true
        row["NIE_true"] <- NIE_true
        row
      })

      out <- na_result_row(result_cols)
      out[names(one_rep)] <- one_rep[names(one_rep)]
      out
    }

    if (is.null(dim(err_vec))) {
      err_vec <- matrix(err_vec, nrow = 1, dimnames = list(NULL, names(err_vec)))
    }

    ############################################################
    ### Collect this setting's iterations (kept in memory + full RDS).
    ### No per-(pq,n) CSV is written; only integrated summaries are saved.
    ############################################################

    err_df <- as.data.frame(err_vec)
    err_df$DGP_NAME <- DGP_NAME
    err_df$seed <- seed
    err_df$Iter <- Iter
    err_df$p <- pq
    err_df$q <- pq
    err_df$sigma <- sigma_eps
    err_df$s <- s
    err_df$s_B <- s_B

    cat(
      "Finished |",
      "pq =", pq,
      "| n =", nn,
      "| Mean Bias CF =", round(mean(err_df$bias_NIE_cf, na.rm = TRUE), 4),
      "| RMSE CF =", round(sqrt(mean(err_df$bias_NIE_cf^2, na.rm = TRUE)), 4),
      "| SD CF =", round(sd(err_df$NIE_hat_cf, na.rm = TRUE), 4),
      "| Mean Analytical SE CF =", round(mean(err_df$se_NIE_cf, na.rm = TRUE), 4),
      "| Mean Bias NonCF =", round(mean(err_df$bias_NIE_noncf, na.rm = TRUE), 4),
      "| RMSE NonCF =", round(sqrt(mean(err_df$bias_NIE_noncf^2, na.rm = TRUE)), 4),
      "| Failed =", sum(err_df$failed, na.rm = TRUE),
      "\n"
    )

    res[[paste0("pq_", pq, "_n_", nn)]] <- err_df
  }
}

############################################################
### Merge with any previous stage's RDS (two-stage low/high-pq runs)
### then save the combined object. This lets pq=1250 be run separately
### (fewer cores, to avoid OOM) and merged with the low-pq results.
############################################################

rds_path <- file.path(output_dir, paste0("iter_outputs_", file_tag, "_res.rds"))
if (file.exists(rds_path)) {
  prev <- tryCatch(readRDS(rds_path), error = function(e) NULL)
  if (is.list(prev)) {
    for (k in names(prev)) if (is.null(res[[k]])) res[[k]] <- prev[[k]]
  }
}
saveRDS(res, file = rds_path)

## Rebuild summaries over ALL pq present in the merged object.
pq_vec <- sort(unique(as.integer(sub("^pq_([0-9]+)_n_.*$", "\\1", names(res)))))

############################################################
### Summary 1: Monte Carlo CI table
############################################################

results_mc_df <- data.frame()

for (pq in pq_vec){

  s_info <- get_s_values(pq, DGP_NAME)
  s   <- s_info$s
  s_B <- s_info$s_B

  for (nn in n_vec){

    key <- paste0("pq_", pq, "_n_", nn)

    if (is.null(res[[key]])) {
      cat("Missing:", key, "\n")
      next
    }

    seg <- res[[key]]

    if (nrow(seg) == 0) {
      cat("Empty:", key, "\n")
      next
    }

    NIE_hat  <- seg$NIE_hat_cf
    NIE_true_setting <- seg$NIE_true[1]

    bias <- NIE_hat - NIE_true_setting

    bias_mean <- mean(bias, na.rm = TRUE)
    rmse <- sqrt(mean(bias^2, na.rm = TRUE))
    sd_hat <- sd(NIE_hat, na.rm = TRUE)

    ci_lower <- NIE_hat - 1.96 * sd_hat
    ci_upper <- NIE_hat + 1.96 * sd_hat

    ci_length <- mean(ci_upper - ci_lower, na.rm = TRUE)

    coverage <- mean(
      ci_lower <= NIE_true_setting &
        NIE_true_setting <= ci_upper,
      na.rm = TRUE
    )

    results_mc_df <- rbind(
      results_mc_df,
      data.frame(
        DGP_NAME = DGP_NAME,
        seed = seed,
        Iter = Iter,
        n = nn,
        p = pq,
        q = pq,
        s = s,
        s_B = s_B,
        Var = sigma_eps,
        NIE_true = NIE_true_setting,
        bias_mean_NIE = round(bias_mean, 3),
        rmse_NIE = round(rmse, 3),
        sd_NIE = round(sd_hat, 3),
        ci_length_NIE = round(ci_length, 3),
        coverage_NIE = round(coverage, 3),
        n_valid = sum(!is.na(NIE_hat)),
        n_failed = sum(seg$failed, na.rm = TRUE)
      )
    )
  }
}

cat("\n================ Debiasing Monte Carlo CI Table ================\n")
print(results_mc_df)

write.csv(
  results_mc_df,
  file = file.path(output_dir, paste0("summary_", file_tag, "_MC.csv")),
  row.names = FALSE
)

saveRDS(
  results_mc_df,
  file = file.path(output_dir, paste0("summary_", file_tag, "_MC.rds"))
)

############################################################
### Summary 2: analytical-variance CI table
############################################################

results_var_df <- data.frame()

for (pq in pq_vec){

  s_info <- get_s_values(pq, DGP_NAME)
  s   <- s_info$s
  s_B <- s_info$s_B

  for (nn in n_vec){

    key <- paste0("pq_", pq, "_n_", nn)

    if (is.null(res[[key]])) {
      cat("Missing:", key, "\n")
      next
    }

    seg <- res[[key]]

    if (nrow(seg) == 0) {
      cat("Empty:", key, "\n")
      next
    }

    NIE_hat  <- seg$NIE_hat_cf
    NIE_true_setting <- seg$NIE_true[1]
    se_hat <- seg$se_NIE_cf

    bias <- NIE_hat - NIE_true_setting

    bias_mean <- mean(bias, na.rm = TRUE)
    rmse <- sqrt(mean(bias^2, na.rm = TRUE))
    sd_hat <- sd(NIE_hat, na.rm = TRUE)

    ci_lower <- NIE_hat - 1.96 * se_hat
    ci_upper <- NIE_hat + 1.96 * se_hat

    ci_length <- mean(ci_upper - ci_lower, na.rm = TRUE)

    coverage <- mean(
      ci_lower <= NIE_true_setting &
        NIE_true_setting <= ci_upper,
      na.rm = TRUE
    )

    results_var_df <- rbind(
      results_var_df,
      data.frame(
        DGP_NAME = DGP_NAME,
        seed = seed,
        Iter = Iter,
        n = nn,
        p = pq,
        q = pq,
        s = s,
        s_B = s_B,
        Var = sigma_eps,
        NIE_true = NIE_true_setting,
        bias_mean_NIE = round(bias_mean, 3),
        rmse_NIE = round(rmse, 3),
        sd_NIE = round(sd_hat, 3),
        ci_length_NIE = round(ci_length, 3),
        coverage_NIE = round(coverage, 3),
        mean_se_NIE = round(mean(se_hat, na.rm = TRUE), 4),
        mean_se_NIE_nocov = round(mean(seg$se_NIE_cf_nocov, na.rm = TRUE), 4),
        n_valid = sum(!is.na(NIE_hat)),
        n_failed = sum(seg$failed, na.rm = TRUE)
      )
    )
  }
}

cat("\n================ Debiasing Analytical Variance CI Table ================\n")
print(results_var_df)

write.csv(
  results_var_df,
  file = file.path(output_dir, paste0("summary_", file_tag, "_Analytical.csv")),
  row.names = FALSE
)

saveRDS(
  results_var_df,
  file = file.path(output_dir, paste0("summary_", file_tag, "_Analytical.rds"))
)

############################################################
### Plots: best setting
############################################################

if (run_plots && nrow(results_var_df) > 0) {

plot_choice <- results_var_df %>%
  filter(is.finite(coverage_NIE), is.finite(rmse_NIE)) %>%
  mutate(score = abs(coverage_NIE - 0.95) + 0.2 * rmse_NIE) %>%
  arrange(score) %>%
  slice(1)

if (nrow(plot_choice) > 0) {

  plot_pq <- plot_choice$p[1]
  plot_n  <- plot_choice$n[1]
  plot_key <- paste0("pq_", plot_pq, "_n_", plot_n)

  segment_plot <- res[[plot_key]]

  res_cf    <- sqrt(plot_n) * segment_plot$bias_NIE_cf
  res_noncf <- sqrt(plot_n) * segment_plot$bias_NIE_noncf

  res_df <- data.frame(
    x = c(res_cf, res_noncf),
    type = rep(c("Debiasing_CF", "NonCF"), each = length(res_cf))
  )

  x_lim_best <- c(min(min(res_df$x, na.rm = TRUE), 0), max(max(res_df$x, na.rm = TRUE), 0))
  if (!is.finite(diff(x_lim_best)) || diff(x_lim_best) == 0) {
    x_lim_best <- x_lim_best + c(-1, 1)
  }
  x_breaks_best <- sort(unique(c(pretty(x_lim_best, n = 8), 0)))

  hist_plot <- ggplot(res_df, aes(x = x, fill = type)) +
    geom_histogram(bins = 35, color = "#e9ecef", alpha = 0.8, position = "identity") +
    scale_fill_manual(values = c("Debiasing_CF" = "cornflowerblue", "NonCF" = "tomato")) +
    scale_x_continuous(limits = x_lim_best, breaks = x_breaks_best) +
    theme_minimal() +
    ggtitle(expression(paste("Histogram of ", sqrt(n) * (hat(NIE) - NIE)))) +
    xlab(expression(paste(sqrt(n) * (hat(NIE) - NIE)))) +
    ylab("Frequency")

  z_cf <- segment_plot$bias_NIE_cf / segment_plot$se_NIE_cf
  z_cf <- z_cf[is.finite(z_cf)]

  if (length(z_cf) >= 2) {
    qq_cf <- qqnorm(z_cf, plot.it = FALSE)
    qq_df <- data.frame(
      theoretical_quantiles = qq_cf$x,
      sample_quantiles = qq_cf$y
    )
  } else {
    qq_df <- data.frame(theoretical_quantiles = numeric(0), sample_quantiles = numeric(0))
  }

  qq_plot <- ggplot(qq_df, aes(x = theoretical_quantiles, y = sample_quantiles)) +
    geom_point(color = "cornflowerblue", alpha = 0.8) +
    theme_minimal() +
    ggtitle("QQ Plot") +
    xlab("Theoretical Normal Quantiles") +
    ylab(expression(paste((hat(NIE) - NIE) / hat(se))))

  combined_plot <- arrangeGrob(hist_plot, qq_plot, ncol = 2)

  ggsave(
    filename = file.path(output_dir, paste0("hist_qq_best_", file_tag, "_pq_", plot_pq, "_n_", plot_n, ".png")),
    plot = combined_plot,
    width = 12,
    height = 5,
    dpi = 300
  )

  ggsave(
    filename = file.path(output_dir, paste0("hist_qq_best_", file_tag, "_pq_", plot_pq, "_n_", plot_n, ".pdf")),
    plot = combined_plot,
    width = 12,
    height = 5
  )
}

}  ### end if (run_plots)

cat("\nAll debiasing no-interaction outputs saved to:", output_dir, "\n")

try(parallel::stopCluster(my.cluster), silent = TRUE)
