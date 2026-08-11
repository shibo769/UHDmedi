################################################################################
### Simu_modified.R  (PROTOTYPE -- new cross-fitting algorithm)
###
### Implements the MODIFIED estimation procedure (Algorithms 1-3):
###   Alg 1: mediation functional  theta_10 = E[Y(1, M(0))] at Xbar
###   Alg 2: theta_00 (debiased control-outcome mean at Xbar) -> NDE
###   Alg 3: theta_11 (debiased treated-outcome mean at Xbar) -> NIE
###
### Key differences from the original algorithm (original_v1):
###   * NO estimation of the full B0 matrix (no q-column lasso loop).
###     The contrast component B0_hat %*% Xbar is replaced by
###     m_hat = M_c^T tau2 (balancing-weighted observed mediators).
###   * tau2 is solved on fold k and REUSED twice: to build the contrast
###     and to debias vartheta_2.
###   * The b-hat regression also lives on fold k (same fold as tau2).
###   * Fold division of labour (direction k -> k'): fold k handles the
###     control side, fold k' handles the treated side.
###
### Deb's implementation requirements:
###   1. ONE data split shared by all three algorithms; swap folds and rerun;
###      average the two directions.
###   2. STRATIFIED split: treated cut in half + control cut in half,
###      combined into D1/D2.
###   3. Report only theta_10 for now, but SAVE everything.
###   4. Analytical variance keeps the old sandwich form (Deb will re-check).
###   5. mu_X must keep mu^T Sigma^{-1} mu bounded (controlled by
###      MC_MUX_SPARSE).
###
### Quick test:
###   MC_ITER=2 MC_PQ=50 MC_N=250 MC_CORES=4 Rscript Simu_modified.R
################################################################################

rm(list = ls())                                            # clear workspace to avoid stale objects

################################################################################
### Thread guard: force single-threaded BLAS/OpenMP BEFORE loading glmnet
################################################################################

Sys.setenv(                                                # set env vars (must precede library loads)
  OMP_NUM_THREADS        = "1",                            # OpenMP single thread
  OPENBLAS_NUM_THREADS   = "1",                            # OpenBLAS single thread
  MKL_NUM_THREADS        = "1",                            # MKL single thread (if present)
  BLIS_NUM_THREADS       = "1",                            # BLIS single thread (if present)
  VECLIB_MAXIMUM_THREADS = "1"                             # macOS vecLib single thread (compat)
)

################################################################################
### Paths: reuse shared plumbing + balancing solver from original_v1
################################################################################

get_script_dir <- function() {                             # helper: locate this script's directory
  args <- commandArgs(trailingOnly = FALSE)                # full command line of Rscript
  file_arg <- grep("^--file=", args, value = TRUE)         # find the --file=... entry
  if (length(file_arg) > 0) {                              # if run via Rscript
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))  # return the script directory
  }                                                        #
  normalizePath(getwd())                                   # fallback (interactive): working directory
}                                                          #

base_dir <- normalizePath(Sys.getenv("BASE_DIR", unset = get_script_dir()), mustWork = TRUE)  # this script's dir (BASE_DIR overrides)
repo_root <- dirname(base_dir)                             # repository root ~/paper_code/mediation_sparsity

SHARED <- file.path((function(){a<-commandArgs(FALSE);f<-sub("^--file=","",a[grep("^--file=",a)]);if(length(f)) dirname(normalizePath(f)) else normalizePath(getwd())})(), "..", "shared")   # bundled deps (was original_v1/*)
source(file.path(SHARED, "mc_common.R"))                 # shared utils (env parsing / cluster / safe_cv_glmnet)
source(file.path(SHARED, "utilis.R"))       # approx.balance & friends (balancing-weight QP solver)

################################################################################
### Package check (never auto-install)
################################################################################

required_pkgs <- c(                                        # packages this script needs
  "foreach", "doParallel", "doRNG", "pracma", "MASS",      # parallel + randortho + mvrnorm
  "glmnet", "quadprog", "ggplot2", "gridExtra", "dplyr"    # lasso + QP fallback + plotting + data wrangling
)                                                          #
require_packages(required_pkgs)                            # stop with instructions if anything missing

suppressPackageStartupMessages({                           # load quietly
  library(foreach); library(doParallel); library(doRNG)    # parallel framework + reproducible RNG
  library(pracma);  library(MASS);       library(glmnet)   # randortho + mvrnorm + lasso
  library(quadprog); library(ggplot2);   library(gridExtra); library(dplyr)  # remaining deps
})                                                         #

if (!requireNamespace("pogs", quietly = TRUE)) {           # check the pogs solver is available
  stop("The 'pogs' package is required (optimizer = 'pogs.new'). ",  # hard stop if missing
       "Install with: Rscript shared/setup_packages.R")         # with install instructions
}                                                          #

################################################################################
### Parameters (all controllable via environment variables)
################################################################################

Iter      <- env_int("MC_ITER", 500L)                      # number of Monte Carlo replications
pq_vec    <- env_int_vec("MC_PQ", "50,400,750,1250")       # grid of p = q values
n_vec     <- env_int_vec("MC_N", "250,500,750,1000,1250")  # grid of sample sizes
sigma_eps <- env_num("MC_SIGMA", 1)                        # response noise sd
DGP_NAME  <- env_chr("DGP_NAME", "base")                   # DGP name (prototype: base only)
seed      <- env_int("MC_SEED", 20260704L)                 # RNG seed
nfolds    <- env_int("MC_NFOLDS", 3L)                      # folds for cv.glmnet
run_plots <- env_lgl("RUN_PLOTS", TRUE)                    # plotting switch (unused in prototype)
skip_naive <- env_lgl("SKIP_NAIVE", FALSE)                 # TRUE = skip the naive estimator (fast test runs)

alpha_inf <- env_num("MC_ALPHA", 2/3)                      # exponent in ||tau||_inf <= n^{-alpha} (Alg1 Step1 / Alg3)
K1 <- env_num("MC_K1", 2.75)                               # balancing constant K1 (env-tunable)
K2 <- env_num("MC_K2", 2.75)                               # balancing constant K2 (Deb: try reducing in high sparsity)
K3 <- env_num("MC_K3", 2.75)                               # balancing constant K3
normalize_B0 <- env_lgl("MC_NORM_B0", FALSE)               # TRUE = normalize rows of B0/B1 to l1 norm 1 (Deb suggestion)

mux_mode <- env_chr("MC_MUX_MODE", "scaled")               # mu_X mode: "scaled" (c/sqrt(p) * 1, bounded), "dense" (c * 1, old), "sparse" (first k coords)
mux_c    <- env_num("MC_MUX", 0.5 * sqrt(50))              # constant c; default 0.5*sqrt(50) so p = 50 matches the old 0.5 mean
mux_sparse <- env_int("MC_MUX_SPARSE", 5L)                 # nonzero coords when mux_mode = "sparse"

s   <- env_int("MC_S", 3L)                                 # sparsity of beta/gamma (base/hetero DGP; env-tunable for sensitivity runs)
s_B <- env_int("MC_SB", 3L)                                # row sparsity of B0/B1 (base/hetero DGP; env-tunable for sensitivity runs)

if (!(DGP_NAME %in% c("base", "increasing_sparsity", "increasing_sparsity_B", "minimal_signal", "hetero"))) {  # supported DGPs
  stop("DGP_NAME must be one of: base, increasing_sparsity, increasing_sparsity_B, minimal_signal, hetero")  # guard
}                                                          #
is_hetero <- (DGP_NAME == "hetero")                        # heteroscedastic-noise flag

## DGP helpers (same conventions as original_v1/Main_simu/Simu_1_normal.R):
s_map_interaction <- c("50" = 5L, "400" = 10L, "750" = 15L, "1250" = 20L)  # sparsity grows with dimension

draw_signal <- function(m, dgp_name, ref_n = max(n_vec)) { # signal magnitudes per DGP
  if (dgp_name == "minimal_signal") {                      # minimal-signal regime
    m_weak   <- max(1, floor(m / 2))                       # about half the coords are weak
    m_strong <- m - m_weak                                 # the rest are moderate
    return(c(runif(m_strong, 0.5, 1.5),                    # moderate signals
             runif(m_weak, 0.8, 1.2) / sqrt(ref_n)))       # weak signals ~ n_max^{-1/2}
  }                                                        #
  runif(m, 0.5, 1.5)                                       # base / sparsity / hetero: moderate signals
}                                                          #

set.seed(seed)                                             # fix seed (true-parameter draws reproducible)
options(MC_NFOLDS = nfolds)                                # default folds for safe_cv_glmnet

################################################################################
### Output folder
################################################################################

sigma_tag <- make_sigma_tag(sigma_eps)                     # e.g. sigma_1 / sigma_0p5
file_tag  <- paste0("modified_", DGP_NAME, "_", sigma_tag) # file tag: modified_base_sigma_1

output_dir <- file.path(base_dir, "simulation_outputs_modified", DGP_NAME, sigma_tag)  # output path
ensure_dir(output_dir)                                     # create if missing

cat("MODIFIED algorithm | DGP:", DGP_NAME, "| sigma:", sigma_eps, "\n")            # print run config
cat("Iter:", Iter, "| pq:", paste(pq_vec, collapse = ","),                         # print iterations and pq grid
    "| n:", paste(n_vec, collapse = ","), "\n")                                    # print n grid
cat("alpha (tau_inf exponent):", alpha_inf,                                        # print alpha
    "| mu_X mode:", mux_mode, "| c:", round(mux_c, 4), "\n")                       # print mu_X setting
cat("Outputs:", output_dir, "\n")                                                  # print output location

################################################################################
### Parallel cluster
################################################################################

my.cluster <- make_mc_cluster(Iter, seed)                  # FORK cluster + doRNG registration (reproducible)
on.exit(try(parallel::stopCluster(my.cluster), silent = TRUE), add = TRUE)  # stop cluster on exit

################################################################################
### Helper: raw balancing weights, thin wrapper around approx.balance
###
### Solves:  min 1/2 ||tau||_2^2
###          s.t. || target - Xmat^T tau ||_inf <= K * sqrt(log(d)/m)
###               || tau ||_inf <= m^{-alpha}
### where d = ncol(Xmat), m = nrow(Xmat).
###
### NOTE (check!): approx.balance additionally imposes sum(tau) = 1 and
### tau >= 0 (allow.negative.weights = FALSE) and returns tau normalised to
### sum one. The written algorithm does not state these constraints -- they
### are inherited from the original residual-balancing implementation.
################################################################################

balance_weights <- function(Xmat, target, K, alpha_exp) {  # balancing-weight wrapper
  m <- nrow(Xmat)                                          # sample size m (for the inf bound)
  approx.balance(                                          # call the solver from utilis.R
    M = Xmat,                                              # design matrix
    balance.target = target,                               # balance target vector
    zeta = 0.5,                                            # unused by the pogs.new path (placeholder)
    K = K,                                                 # balancing bound constant K
    allow.negative.weights = FALSE,                        # nonnegative weights (inherited)
    optimizer = "pogs.new",                                # K-bound formulation via pogs
    bound.gamma = TRUE,                                    # enable the ||tau||_inf cap
    gamma.max = m^(-alpha_exp)                             # cap value m^{-alpha}
  )                                                        #
}                                                          #

################################################################################
### Helper: lasso coefficients (cv.glmnet + lambda.1se, no intercept)
################################################################################

lasso_coef <- function(X, y, nfolds_ = nfolds) {           # lasso wrapper returning coefficient vector
  fit <- safe_cv_glmnet(X, y, alpha = 1, intercept = FALSE, nfolds = nfolds_)  # CV lasso without intercept
  if (is.null(fit)) stop("cv.glmnet failed")               # throw on failure (caught by outer tryCatch)
  as.vector(coef(fit, s = "lambda.1se"))[-1]               # coefficients at lambda.1se, drop intercept slot
}                                                          #

################################################################################
### Helper: sandwich variance of a weighted-residual mean (same form as
### residualBalance.mean's var.hat). Deb: keep the old formula for now.
################################################################################

sandwich_var <- function(tau, resid, df_used) {            # variance of sum(tau_i * resid_i)
  m <- length(tau)                                         # number of weights (= sample size)
  sum(tau^2 * resid^2) * m / max(1, m - df_used)           # sum(tau_i^2 r_i^2) with df correction
}                                                          #

################################################################################
### One direction (k -> k') of the modified Algorithms 1-3.
###
### Inputs: index sets of fold k and fold k' (already stratified), full data,
###         Xbar, alpha exponent.
### Output: named list with this direction's theta_10 estimate, theta_00,
###         theta_11 and all intermediate quantities (everything is saved).
################################################################################

run_direction <- function(idx_k, idx_kp, X, M, Y, A, X_bar, alpha_exp) {  # one cross-fitting direction

  p <- ncol(X); q <- ncol(M)                               # dimensions p, q

  ## ----- fold k: control-side data -------------------------------------------
  ck  <- idx_k[A[idx_k] == 0]                              # control indices in fold k
  Xck <- X[ck, , drop = FALSE]                             # X_{c,k}: control X in fold k
  Mck <- M[ck, , drop = FALSE]                             # M_{c,k}: control M in fold k
  yck <- Y[ck]                                             # y_{c,k}: control Y in fold k (used by Alg 2)
  nck <- length(ck)                                        # n_{c,k}: control count in fold k

  ## ----- fold k': treated-side data -------------------------------------------
  tkp  <- idx_kp[A[idx_kp] == 1]                           # treated indices in fold k'
  Xtkp <- X[tkp, , drop = FALSE]                           # X_{t,k'}: treated X in fold k'
  Mtkp <- M[tkp, , drop = FALSE]                           # M_{t,k'}: treated M in fold k'
  Wtkp <- cbind(Xtkp, Mtkp)                                # W_{t,k'} = (X, M) for treated in fold k'
  ytkp <- Y[tkp]                                           # y_{t,k'}: treated Y in fold k'
  ntkp <- length(tkp)                                      # n_{t,k'}: treated count in fold k'

  ## ============ Alg 1, Step 1: tau2 on fold-k control X =======================
  tau2 <- balance_weights(Xck, X_bar, K = K2, alpha_exp = alpha_exp)  # tau2^{(k)}: target Xbar, bound K2*sqrt(log p/n_ck)

  ## ============ Alg 1, Step 2: phi1 lasso on fold-k' treated ==================
  phi1 <- lasso_coef(Wtkp, ytkp)                           # phi1^{(k')} = (beta1, gamma1) lasso estimate
  beta1_hat  <- phi1[1:p]                                  # first p entries: beta1_hat
  gamma1_hat <- phi1[(p + 1):(p + q)]                      # last q entries: gamma1_hat
  df_phi1 <- sum(phi1 != 0)                                # nonzero count of phi1 (df correction)

  ## ============ Alg 1, Step 3: contrast + tau1 + vartheta_1 ===================
  m_hat  <- as.vector(t(Mck) %*% tau2)                     # m_hat_{c,k} = M_{c,k}^T tau2 (replaces B0_hat %*% Xbar)
  a_con  <- c(X_bar, m_hat)                                # contrast a^{(k)} = (Xbar, m_hat)

  tau1 <- balance_weights(Wtkp, a_con, K = K1, alpha_exp = 2/3)  # tau1^{(k')}: target a, bound K1*sqrt(log(p+q)/n_tk'), inf cap fixed 2/3

  resid1     <- as.vector(ytkp - Wtkp %*% phi1)            # treated-outcome lasso residuals
  vartheta1  <- sum(a_con * phi1) + sum(tau1 * resid1)     # vartheta_1 = a^T phi1 + tau1^T residuals
  var_vth1   <- sandwich_var(tau1, resid1, df_phi1)        # sandwich variance of vartheta_1

  ## ============ Alg 1, Step 4: b-hat on fold-k control + vartheta_2 ===========
  pseudo_y <- as.vector(Mck %*% gamma1_hat)                # pseudo-response M_{c,k} %*% gamma1_hat^{(k')}
  b_hat    <- lasso_coef(Xck, pseudo_y)                    # b_hat^{(k)}: lasso of pseudo-response on X_{c,k}
  df_b     <- sum(b_hat != 0)                              # nonzero count of b_hat

  resid2     <- as.vector(pseudo_y - Xck %*% b_hat)        # pseudo-response lasso residuals
  vartheta2  <- sum(X_bar * b_hat) + sum(tau2 * resid2)    # vartheta_2 = Xbar^T b_hat + tau2^T residuals (tau2 REUSED)
  var_vth2   <- sandwich_var(tau2, resid2, df_b)           # sandwich variance of vartheta_2

  ## ============ Alg 1, Step 5: this direction's theta_10 ======================
  theta10_dir <- vartheta1 + vartheta2 - sum(m_hat * gamma1_hat)  # theta10 = v1 + v2 - m_hat^T gamma1 (remove overlap)

  ## ============ Alg 2 (theta_00): rho_c on fold-k control + reuse tau2 ========
  rho_c <- lasso_coef(Xck, yck)                            # rho_c^{(k)}: lasso of y_{c,k} on X_{c,k}
  df_rho_c <- sum(rho_c != 0)                              # nonzero count of rho_c
  resid00  <- as.vector(yck - Xck %*% rho_c)               # control-outcome residuals
  theta00_dir <- sum(X_bar * rho_c) + sum(tau2 * resid00)  # theta00 = Xbar^T rho_c + tau2^T residuals (tau2 reused again)
  var_th00    <- sandwich_var(tau2, resid00, df_rho_c)     # sandwich variance of theta00

  ## ============ Alg 3 (theta_11): fold-k treated, own omega ===================
  ## NOTE: Alg 3 works within a single fold (rho_t and omega both on fold-k
  ## treated); it has no k -> k' structure.
  tk  <- idx_k[A[idx_k] == 1]                              # treated indices in fold k
  Xtk <- X[tk, , drop = FALSE]                             # treated X in fold k
  ytk <- Y[tk]                                             # treated Y in fold k

  rho_t <- lasso_coef(Xtk, ytk)                            # rho_t^{(k)}: lasso of y_{t,k} on X_{t,k} (no M!)
  df_rho_t <- sum(rho_t != 0)                              # nonzero count of rho_t
  omega <- balance_weights(Xtk, X_bar, K = K3, alpha_exp = alpha_exp)  # omega^{(k)}: target Xbar, bound K3*sqrt(log p/n_tk)
  resid11 <- as.vector(ytk - Xtk %*% rho_t)                # treated-outcome (X-only) residuals
  theta11_fold <- sum(X_bar * rho_t) + sum(omega * resid11)  # theta11 = Xbar^T rho_t + omega^T residuals
  var_th11     <- sandwich_var(omega, resid11, df_rho_t)   # sandwich variance of theta11

  list(                                                    # bundle every output of this direction
    theta10_dir = theta10_dir,                             # this direction's theta_10
    vartheta1 = vartheta1, vartheta2 = vartheta2,          # the two components
    overlap_term = sum(m_hat * gamma1_hat),                # the subtracted overlap term
    var_vth1 = var_vth1, var_vth2 = var_vth2,              # component variances
    theta00_dir = theta00_dir, var_th00 = var_th00,        # Alg 2 outputs
    theta11_fold = theta11_fold, var_th11 = var_th11,      # Alg 3 outputs
    tau1_norm = sum(tau1^2), tau2_norm = sum(tau2^2),      # squared 2-norms of weights (diagnostics)
    omega_norm = sum(omega^2),                             # squared 2-norm of omega (diagnostics)
    tau1_max = max(abs(tau1)), tau2_max = max(abs(tau2)),  # inf-norms of weights (check the n^-alpha cap binds)
    omega_max = max(abs(omega)),                           # inf-norm of omega
    sigma1_hat = mean(resid1^2),                           # residual variance of treated outcome lasso (~ sigma^2)
    vU_hat = mean(resid2^2),                               # residual variance of pseudo-response lasso (~ gamma' Sigma_U gamma)
    sigma0_hat = mean(resid00^2),                          # residual variance of control outcome lasso (Alg 2)
    sigma11_hat = mean(resid11^2),                         # residual variance of treated X-only lasso (Alg 3)
    df_phi1 = df_phi1, df_b = df_b,                        # lasso support sizes (Alg 1)
    df_rho_c = df_rho_c, df_rho_t = df_rho_t,              # lasso support sizes (Alg 2 / Alg 3)
    m_hat_norm = sum(m_hat^2),                             # squared 2-norm of the weighted-mediator contrast
    n_ck = nck, n_tkp = ntkp, n_tk = length(tk)            # subsample sizes (diagnostics)
  )                                                        #
}                                                          #

################################################################################
### Result column template (failed replications return all-NA + failed = 1)
################################################################################

result_cols <- c(                                          # column names returned per replication
  "theta10_cf", "theta10_true", "bias10",                  # theta10 CF estimate / truth / bias
  "theta10_naive", "bias10_naive",                         # naive plug-in estimate / bias
  "se10_analytical",                                       # analytical SE of theta10
  "se00_analytical", "se11_analytical",                    # analytical SEs of theta00 / theta11
  "theta00_cf", "theta00_true",                            # theta00 estimate / truth
  "theta11_cf", "theta11_true",                            # theta11 estimate / truth
  "NDE_cf", "NDE_true", "NIE_cf", "NIE_true",              # NDE / NIE estimates and truths
  "theta10_d12", "theta10_d21",                            # per-direction theta10
  "theta00_d12", "theta00_d21",                            # per-direction theta00
  "theta11_f1", "theta11_f2",                              # per-fold theta11
  "vartheta1_d12", "vartheta2_d12", "overlap_d12",         # direction 1->2 components
  "vartheta1_d21", "vartheta2_d21", "overlap_d21",         # direction 2->1 components
  "var_vth1_d12", "var_vth2_d12", "var_vth1_d21", "var_vth2_d21",  # component variances
  "var_th00_d12", "var_th00_d21", "var_th11_f1", "var_th11_f2",    # theta00 / theta11 variances
  "tau1_norm_d12", "tau2_norm_d12", "omega_norm_d12",      # weight 2-norms, direction 1->2
  "tau1_norm_d21", "tau2_norm_d21", "omega_norm_d21",      # weight 2-norms, direction 2->1
  "tau1_max_d12", "tau2_max_d12", "omega_max_d12",         # weight inf-norms, direction 1->2
  "tau1_max_d21", "tau2_max_d21", "omega_max_d21",         # weight inf-norms, direction 2->1
  "sigma1_hat_d12", "sigma1_hat_d21",                      # treated-outcome residual variances
  "vU_hat_d12", "vU_hat_d21",                              # pseudo-response residual variances
  "sigma0_hat_d12", "sigma0_hat_d21",                      # control-outcome residual variances (Alg 2)
  "sigma11_hat_d12", "sigma11_hat_d21",                    # X-only treated residual variances (Alg 3)
  "df_phi1_d12", "df_b_d12", "df_rho_c_d12", "df_rho_t_d12",  # lasso support sizes, direction 1->2
  "df_phi1_d21", "df_b_d21", "df_rho_c_d21", "df_rho_t_d21",  # lasso support sizes, direction 2->1
  "m_hat_norm_d12", "m_hat_norm_d21",                      # contrast norms
  "n_ck_d12", "n_tkp_d12", "n_tk_d12",                     # subsample sizes, direction 1->2
  "n_ck_d21", "n_tkp_d21", "n_tk_d21",                     # subsample sizes, direction 2->1
  "failed"                                                 # failure flag
)                                                          #

################################################################################
### Main simulation loop
################################################################################

res <- list()                                              # container for all (pq, n) results

for (pq in pq_vec) {                                       # outer loop: dimension

  p <- pq; q <- pq                                         # p = q = pq

  cat("\n============================\n")                  # separator
  cat("Starting pq =", pq, "\n")                           # print current dimension
  cat("============================\n")                    # separator

  ### Per-DGP sparsity ###
  ## increasing_sparsity: both s(gamma/beta) and s_B grow with dimension (old setting)
  ## increasing_sparsity_B: only s_B grows, gamma/beta sparsity k stays fixed at s=3
  s_cur <- if (DGP_NAME == "increasing_sparsity") as.integer(s_map_interaction[[as.character(pq)]]) else s               # beta/gamma sparsity
  s_B_cur <- if (DGP_NAME %in% c("increasing_sparsity", "increasing_sparsity_B")) as.integer(s_map_interaction[[as.character(pq)]]) else s_B  # B row sparsity

  cat("s =", s_cur, "| s_B =", s_B_cur, "\n")              # print effective sparsity

  ### True coefficients (fixed across n and replications) ###
  beta_0 <- rep(0, p); beta_0[sample(1:p, s_cur)] <- draw_signal(s_cur, DGP_NAME)   # control-arm X coefficients
  beta_1 <- rep(0, p); beta_1[sample(1:p, s_cur)] <- draw_signal(s_cur, DGP_NAME)   # treated-arm X coefficients
  gamma_0 <- rep(0, q); gamma_0[sample(1:q, s_cur)] <- draw_signal(s_cur, DGP_NAME) # control-arm M coefficients
  gamma_1 <- rep(0, q); gamma_1[sample(1:q, s_cur)] <- draw_signal(s_cur, DGP_NAME) # treated-arm M coefficients

  B0 <- matrix(0, q, p); B1 <- matrix(0, q, p)             # the two mediator regression matrices
  for (i in 1:q) {                                         # fill nonzeros row by row
    B0[i, sample(1:p, s_B_cur)] <- runif(s_B_cur, 0.5, 1.5)  # s_B nonzeros in row i of B0
    B1[i, sample(1:p, s_B_cur)] <- runif(s_B_cur, 0.5, 1.5)  # s_B nonzeros in row i of B1
  }                                                        #
  if (normalize_B0) {                                      # Deb: normalize rows so l1 norm is 1
    B0 <- B0 / pmax(rowSums(abs(B0)), 1e-12)               # each row of B0 has ||row||_1 = 1
    B1 <- B1 / pmax(rowSums(abs(B1)), 1e-12)               # same for B1 (keep arms symmetric)
  }                                                        #

  ### Heteroscedastic-noise loadings (fixed per setting; only used if hetero) ###
  a_eps <- rep(0, p); a_eps[sample(1:p, 3)] <- runif(3, 0, 1)  # sparse X loadings for v(X,M)
  b_eps <- rep(0, q); b_eps[sample(1:q, 3)] <- runif(3, 0, 1)  # sparse M loadings for v(X,M)

  for (nn in n_vec) {                                      # inner loop: sample size

    ### X generation; mu_X can be sparse per Deb's boundedness requirement ###
    mu_vec <- switch(mux_mode,                             # build mu_X according to the chosen mode
      scaled = rep(mux_c / sqrt(p), p),                    # dense, dimension-rescaled: c/sqrt(p)*1 => ||mu||^2 = c^2 bounded
      dense  = rep(mux_c, p),                              # constant dense mean (old setting; mu' Sig^{-1} mu ~ p, unbounded)
      sparse = c(rep(mux_c, min(mux_sparse, p)), rep(0, max(0, p - mux_sparse))),  # first k coords nonzero
      stop("Unknown MC_MUX_MODE: ", mux_mode)              # guard against typos
    )                                                      #
    U_mat   <- randortho(p)                                # random orthonormal matrix
    Sigma_X <- U_mat %*% diag(runif(p, 1, 2)) %*% t(U_mat) # Sigma_X with eigenvalues in [1, 2]
    X <- mvrnorm(nn, mu_vec, Sigma_X)                      # draw X ~ N(mu, Sigma_X)

    ### A generation (same as original base DGP) ###
    alpha_A <- rep(0, p); alpha_A[sample(1:p, s)] <- runif(s, 0, 2)  # sparse propensity coefficients
    linpred <- as.vector(X %*% alpha_A)                    # linear predictor
    linpred <- linpred - mean(linpred)                     # centre
    if (is.finite(sd(linpred)) && sd(linpred) > 1e-8) linpred <- linpred / sd(linpred)  # scale (avoid all-treated/control)
    A <- as.numeric(rbinom(nn, 1, plogis(linpred)))        # draw treatment indicator via logistic

    treat_index   <- which(A == 1)                         # treated indices
    control_index <- which(A == 0)                         # control indices
    if (min(length(treat_index), length(control_index)) < 20) {  # skip if either arm too small
      warning("Skipping pq=", pq, " n=", nn, ": too few treated/control.")  # notify
      next                                                 # skip this cell
    }                                                      #

    ### True values (conditional on this X) ###
    X_bar <- colMeans(X)                                   # sample mean Xbar
    theta10_true <- as.numeric(t(X_bar) %*% (beta_1 + t(B0) %*% gamma_1))  # true theta_10
    theta00_true <- as.numeric(t(X_bar) %*% (beta_0 + t(B0) %*% gamma_0))  # true theta_00
    theta11_true <- as.numeric(t(X_bar) %*% (beta_1 + t(B1) %*% gamma_1))  # true theta_11
    NDE_true <- theta10_true - theta00_true                # true NDE
    NIE_true <- theta11_true - theta10_true                # true NIE

    ### Monte Carlo ###
    err_vec <- foreach(                                    # parallel replications
      i = 1:Iter,                                          # replication index
      .combine = "rbind",                                  # row-bind the results
      .packages = c("MASS", "glmnet", "pracma", "quadprog", "pogs"),  # packages needed by workers
      .options.multicore = mc_foreach_options_multicore,   # FORK: no prescheduling (load balancing)
      .options.snow = mc_foreach_options_snow              # PSOCK equivalent (placeholder)
    ) %dopar% {                                            # parallel body begins

      worker_single_thread_blas()                          # re-cap BLAS threads inside the worker

      one_rep <- tryCatch({                                # one replication (NA row on failure)

        ### Generate M and Y (same interaction DGP as original base) ###
        lambda_U <- runif(q, 1, 2)                         # eigenvalues of Sigma_U
        R_U <- randortho(q)                                # random orthonormal matrix
        Sigma_U <- R_U %*% diag(lambda_U) %*% t(R_U)       # mediator noise covariance
        noise_M0 <- mvrnorm(length(control_index), rep(0, q), Sigma_U)  # control-arm M noise
        noise_M1 <- mvrnorm(length(treat_index),  rep(0, q), Sigma_U)   # treated-arm M noise

        M <- matrix(0, nn, q)                              # initialise M
        M[control_index, ] <- X[control_index, , drop = FALSE] %*% t(B0) + noise_M0  # control: M = B0 X + U
        M[treat_index, ]   <- X[treat_index, , drop = FALSE]   %*% t(B1) + noise_M1  # treated: M = B1 X + U'

        ### Response noise: homoscedastic N(0, sigma^2), or heteroscedastic
        ### sd * sqrt(v(X,M)) with v normalised to mean one (hetero DGP). ###
        if (is_hetero) {                                   # heteroscedastic-noise regime
          vc_raw <- 1 + 0.5 * ((X[control_index, , drop = FALSE] %*% a_eps)^2 +  # v raw, control arm
                               (M[control_index, , drop = FALSE] %*% b_eps)^2)   #
          vt_raw <- 1 + 0.5 * ((X[treat_index, , drop = FALSE] %*% a_eps)^2 +    # v raw, treated arm
                               (M[treat_index, , drop = FALSE] %*% b_eps)^2)     #
          v_mean <- mean(c(vc_raw, vt_raw))                # normalise so E[v] = 1 (marginal level = sigma^2)
          noise_Yc <- rnorm(length(control_index), 0, sigma_eps) * sqrt(as.vector(vc_raw) / v_mean)  # control noise
          noise_Yt <- rnorm(length(treat_index),  0, sigma_eps) * sqrt(as.vector(vt_raw) / v_mean)   # treated noise
        } else {                                           # homoscedastic regimes
          noise_Yc <- rnorm(length(control_index), 0, sigma_eps)  # control noise
          noise_Yt <- rnorm(length(treat_index),  0, sigma_eps)   # treated noise
        }                                                  #

        Y <- rep(0, nn)                                    # initialise Y
        Y[control_index] <- as.vector(                     # control-arm Y
          X[control_index, , drop = FALSE] %*% beta_0 +    #   X beta_0
          M[control_index, , drop = FALSE] %*% gamma_0 +   # + M gamma_0
          noise_Yc)                                        # + noise
        Y[treat_index] <- as.vector(                       # treated-arm Y
          X[treat_index, , drop = FALSE] %*% beta_1 +      #   X beta_1
          M[treat_index, , drop = FALSE] %*% gamma_1 +     # + M gamma_1
          noise_Yt)                                        # + noise

        ### NAIVE estimator of theta_10 (plug-in, no debiasing): lasso of Y_t
        ### on W_t = (X_t, M_t) for phi, per-column BIC-lasso of M_c on X_c
        ### for B0, then theta_naive = Xbar' beta_hat + (B0_hat Xbar)' gamma_hat.
        ### Skipped when SKIP_NAIVE=TRUE (naive's B0 loop dominates runtime).
        if (skip_naive) {                                  # test runs: skip the expensive naive estimator
          theta10_naive <- NA_real_                        # store NA for naive
        } else {                                           # full runs: compute naive
          Wt_all <- cbind(X[treat_index, , drop = FALSE], M[treat_index, , drop = FALSE])  # full treated (X, M)
          phi_nv <- lasso_coef(Wt_all, Y[treat_index])     # naive lasso on ALL treated obs
          beta_nv  <- phi_nv[1:p]                          # naive beta_1 estimate
          gamma_nv <- phi_nv[(p + 1):(p + q)]              # naive gamma_1 estimate
          B0_nv <- b0_bic_matrix(X[control_index, , drop = FALSE],  # naive B0 (BIC lasso per column)
                                 M[control_index, , drop = FALSE]) #   on ALL control obs
          theta10_naive <- sum(X_bar * beta_nv) + sum((B0_nv %*% X_bar) * gamma_nv)  # naive plug-in theta_10
        }                                                  #

        ### Step 0: ONE stratified split shared by all three algorithms (Deb) ###
        t_perm <- sample(treat_index)                      # shuffle treated indices
        c_perm <- sample(control_index)                    # shuffle control indices
        D1 <- c(t_perm[1:floor(length(t_perm) / 2)],       # D1 = half of treated
                c_perm[1:floor(length(c_perm) / 2)])       #    + half of control
        D2 <- setdiff(1:nn, D1)                            # D2 = everything else

        ### Two directions (each runs all three algorithms internally) ###
        d12 <- run_direction(D1, D2, X, M, Y, A, X_bar, alpha_exp = alpha_inf)  # direction 1 -> 2
        d21 <- run_direction(D2, D1, X, M, Y, A, X_bar, alpha_exp = alpha_inf)  # direction 2 -> 1 (folds swapped)

        ### Cross-fitted combinations ###
        theta10_cf <- 0.5 * (d12$theta10_dir + d21$theta10_dir)   # theta10: average of the two directions
        theta00_cf <- 0.5 * (d12$theta00_dir + d21$theta00_dir)   # theta00: average of the two directions
        theta11_cf <- 0.5 * (d12$theta11_fold + d21$theta11_fold) # theta11: average of the two folds

        NDE_cf <- theta10_cf - theta00_cf                  # NDE = theta10 - theta00
        NIE_cf <- theta11_cf - theta10_cf                  # NIE = theta11 - theta10

        ### Analytical SE of theta10 (old combination form: average the four
        ### component variances; /4 comes from the 1/2 direction weights) ###
        se10 <- sqrt((d12$var_vth1 + d12$var_vth2 +        # direction 1->2 component variances
                      d21$var_vth1 + d21$var_vth2) / 4)    # + direction 2->1, divide by 4, sqrt
        se00 <- sqrt((d12$var_th00 + d21$var_th00) / 4)    # analytical SE of theta00 (two directions, 1/2 weights)
        se11 <- sqrt((d12$var_th11 + d21$var_th11) / 4)    # analytical SE of theta11 (two folds, 1/2 weights)

        c(                                                 # bundle all outputs of this replication
          theta10_cf = theta10_cf, theta10_true = theta10_true,       # theta10 estimate and truth
          bias10 = theta10_cf - theta10_true,                         # theta10 bias
          theta10_naive = theta10_naive,                              # naive plug-in estimate
          bias10_naive = theta10_naive - theta10_true,                # naive bias
          se10_analytical = se10,                                     # theta10 analytical SE
          se00_analytical = se00, se11_analytical = se11,             # theta00 / theta11 analytical SEs
          theta00_cf = theta00_cf, theta00_true = theta00_true,       # theta00
          theta11_cf = theta11_cf, theta11_true = theta11_true,       # theta11
          NDE_cf = NDE_cf, NDE_true = NDE_true,                       # NDE
          NIE_cf = NIE_cf, NIE_true = NIE_true,                       # NIE
          theta10_d12 = d12$theta10_dir, theta10_d21 = d21$theta10_dir,   # per-direction theta10
          theta00_d12 = d12$theta00_dir, theta00_d21 = d21$theta00_dir,   # per-direction theta00
          theta11_f1 = d12$theta11_fold, theta11_f2 = d21$theta11_fold,   # per-fold theta11
          vartheta1_d12 = d12$vartheta1, vartheta2_d12 = d12$vartheta2,   # direction 1->2 components
          overlap_d12 = d12$overlap_term,                                 # direction 1->2 overlap term
          vartheta1_d21 = d21$vartheta1, vartheta2_d21 = d21$vartheta2,   # direction 2->1 components
          overlap_d21 = d21$overlap_term,                                 # direction 2->1 overlap term
          var_vth1_d12 = d12$var_vth1, var_vth2_d12 = d12$var_vth2,       # component variances 1->2
          var_vth1_d21 = d21$var_vth1, var_vth2_d21 = d21$var_vth2,       # component variances 2->1
          var_th00_d12 = d12$var_th00, var_th00_d21 = d21$var_th00,       # theta00 variances
          var_th11_f1 = d12$var_th11, var_th11_f2 = d21$var_th11,         # theta11 variances
          tau1_norm_d12 = d12$tau1_norm, tau2_norm_d12 = d12$tau2_norm,   # weight 2-norms 1->2
          omega_norm_d12 = d12$omega_norm,                                # omega 2-norm 1->2
          tau1_norm_d21 = d21$tau1_norm, tau2_norm_d21 = d21$tau2_norm,   # weight 2-norms 2->1
          omega_norm_d21 = d21$omega_norm,                                # omega 2-norm 2->1
          tau1_max_d12 = d12$tau1_max, tau2_max_d12 = d12$tau2_max,       # weight inf-norms 1->2
          omega_max_d12 = d12$omega_max,                                  # omega inf-norm 1->2
          tau1_max_d21 = d21$tau1_max, tau2_max_d21 = d21$tau2_max,       # weight inf-norms 2->1
          omega_max_d21 = d21$omega_max,                                  # omega inf-norm 2->1
          sigma1_hat_d12 = d12$sigma1_hat, sigma1_hat_d21 = d21$sigma1_hat,     # treated-outcome resid vars
          vU_hat_d12 = d12$vU_hat, vU_hat_d21 = d21$vU_hat,                     # pseudo-response resid vars
          sigma0_hat_d12 = d12$sigma0_hat, sigma0_hat_d21 = d21$sigma0_hat,     # control-outcome resid vars
          sigma11_hat_d12 = d12$sigma11_hat, sigma11_hat_d21 = d21$sigma11_hat, # X-only treated resid vars
          df_phi1_d12 = d12$df_phi1, df_b_d12 = d12$df_b,                 # lasso support sizes 1->2
          df_rho_c_d12 = d12$df_rho_c, df_rho_t_d12 = d12$df_rho_t,       # lasso support sizes 1->2 (Alg2/3)
          df_phi1_d21 = d21$df_phi1, df_b_d21 = d21$df_b,                 # lasso support sizes 2->1
          df_rho_c_d21 = d21$df_rho_c, df_rho_t_d21 = d21$df_rho_t,       # lasso support sizes 2->1 (Alg2/3)
          m_hat_norm_d12 = d12$m_hat_norm, m_hat_norm_d21 = d21$m_hat_norm,  # contrast norms
          n_ck_d12 = d12$n_ck, n_tkp_d12 = d12$n_tkp, n_tk_d12 = d12$n_tk,   # subsample sizes 1->2
          n_ck_d21 = d21$n_ck, n_tkp_d21 = d21$n_tkp, n_tk_d21 = d21$n_tk,   # subsample sizes 2->1
          failed = 0                                                      # success flag
        )                                                  #
      }, error = function(e) {                             # on any failure inside the replication
        na_result_row(result_cols)                         # return all-NA row with failed = 1
      })                                                   #

      out <- na_result_row(result_cols)                    # start from the standard NA row (consistent columns)
      out[names(one_rep)] <- one_rep[names(one_rep)]       # overwrite with this replication's values
      out                                                  # return to foreach
    }                                                      # end parallel body

    if (is.null(dim(err_vec))) {                           # foreach returns a vector when Iter = 1
      err_vec <- matrix(err_vec, nrow = 1, dimnames = list(NULL, names(err_vec)))  # normalise to 1-row matrix
    }                                                      #

    cat("Finished | pq =", pq, "| n =", nn,                                        # per-cell summary line
        "| Bias10 =", round(mean(err_vec[, "bias10"], na.rm = TRUE), 4),           # mean bias
        "| RMSE10 =", round(sqrt(mean(err_vec[, "bias10"]^2, na.rm = TRUE)), 4),   # RMSE
        "| SD10 =", round(sd(err_vec[, "theta10_cf"], na.rm = TRUE), 4),           # empirical SD
        "| SE10 =", round(mean(err_vec[, "se10_analytical"], na.rm = TRUE), 4),    # mean analytical SE
        "| Valid =", sum(!is.na(err_vec[, "theta10_cf"])), "/", Iter,              # valid count
        "| Failed =", sum(err_vec[, "failed"], na.rm = TRUE), "\n")                # failed count

    res[[paste0("pq_", pq, "_n_", nn)]] <- err_vec         # store into the result list
  }                                                        # end n loop
}                                                          # end pq loop

################################################################################
### Save: full per-replication RDS + summary tables (theta_10 reported;
### everything else lives in the RDS). RDS ONLY -- no CSV.
################################################################################

saveRDS(res, file.path(output_dir, paste0("iter_outputs_", file_tag, "_res.rds")))  # save all raw replications

## Summary tables cover THREE estimands: theta_10 (mediation functional),
## NDE and NIE. For each estimand we compute MC-CI stats (SD-based) and
## Analytical-CI stats (per-replication analytical SE).
##
## Analytical SEs: theta_10 has its own (se10_analytical). For NDE/NIE the
## components are built on approximately independent folds/arms, so we combine
## variances WITHOUT covariance terms (same convention as before):
##   se_NDE = sqrt(se10^2 + se00^2),  se_NIE = sqrt(se11^2 + se10^2).

summ_mc <- data.frame(); summ_an <- data.frame()           # the two summary tables (MC / analytical)

for (key in names(res)) {                                  # loop over each (pq, n) cell
  seg <- res[[key]]                                        # this cell's replication matrix
  pq <- as.integer(sub("^pq_([0-9]+)_n_.*$", "\\1", key))  # parse pq from the key
  nn <- as.integer(sub("^.*_n_([0-9]+)$", "\\1", key))     # parse n from the key

  se_nde <- sqrt(seg[, "se10_analytical"]^2 + seg[, "se00_analytical"]^2)  # analytical SE for NDE (no cov term)
  se_nie <- sqrt(seg[, "se11_analytical"]^2 + seg[, "se10_analytical"]^2)  # analytical SE for NIE (no cov term)

  first_ok <- function(x) x[which(!is.na(x))[1]]           # truth from the first NON-FAILED row (failed rows store NA truth)
  ests <- list(                                            # the estimands to summarise
    theta10 = list(est = seg[, "theta10_cf"], tru = first_ok(seg[, "theta10_true"]), se = seg[, "se10_analytical"]),  # mediation functional
    NDE     = list(est = seg[, "NDE_cf"],     tru = first_ok(seg[, "NDE_true"]),     se = se_nde),                     # natural direct effect
    NIE     = list(est = seg[, "NIE_cf"],     tru = first_ok(seg[, "NIE_true"]),     se = se_nie),                     # natural indirect effect
    naive   = list(est = seg[, "theta10_naive"], tru = first_ok(seg[, "theta10_true"]),                               # naive plug-in for theta10
                   se = rep(sd(seg[, "theta10_naive"], na.rm = TRUE), nrow(seg)))                                      #   (no analytical SE; MC SD)
  )                                                        #

  for (nm in names(ests)) {                                # loop over the three estimands
    e  <- ests[[nm]]                                       # this estimand's series
    th <- e$est; tt <- e$tru; se <- e$se                   # unpack
    bias <- th - tt                                        # bias series
    sdh  <- sd(th, na.rm = TRUE)                           # empirical SD

    mk_row <- function(se_used) {                          # build one summary row (given the SE to use)
      lo <- th - 1.96 * se_used; hi <- th + 1.96 * se_used # CI bounds
      data.frame(                                          # the summary row
        DGP_NAME = DGP_NAME, seed = seed, Iter = Iter,     # identifier columns
        estimand = nm,                                     # which estimand this row describes
        n = nn, p = pq, q = pq, Var = sigma_eps, s = s, s_B = s_B,  # setting columns
        true_value = tt,                                   # truth
        bias_mean = round(mean(bias, na.rm = TRUE), 3),    # mean bias
        rmse = round(sqrt(mean(bias^2, na.rm = TRUE)), 3), # RMSE
        sd = round(sdh, 3),                                # empirical SD
        ci_length = round(mean(hi - lo, na.rm = TRUE), 3), # mean CI length
        coverage = round(mean(lo <= tt & tt <= hi, na.rm = TRUE), 3),  # empirical coverage
        mean_se_analytical = round(mean(se, na.rm = TRUE), 4),         # mean analytical SE
        n_valid = sum(!is.na(th)),                         # valid replication count
        n_failed = sum(seg[, "failed"], na.rm = TRUE)      # failed replication count
      )                                                    #
    }                                                      #

    summ_mc <- rbind(summ_mc, mk_row(rep(sdh, length(th))))  # MC version: SE = empirical SD (constant)
    summ_an <- rbind(summ_an, mk_row(se))                    # analytical version: per-replication SE
  }                                                        # end estimand loop
}                                                          #

cat("\n=========== MODIFIED: Monte Carlo CI table (theta10 / NDE / NIE) ===========\n")  # print MC table
print(summ_mc)                                             #
cat("\n=========== MODIFIED: Analytical CI table (theta10 / NDE / NIE) ===========\n")   # print analytical table
print(summ_an)                                             #

saveRDS(summ_mc, file.path(output_dir, paste0("table_", file_tag, "_monte_carlo_ci.rds")))          # save MC table (RDS only)
saveRDS(summ_an, file.path(output_dir, paste0("table_", file_tag, "_analytical_variance_ci.rds")))  # save analytical table (RDS only)

cat("\nAll modified-algorithm outputs saved to:", output_dir, "\n")  # done

try(parallel::stopCluster(my.cluster), silent = TRUE)      # stop cluster (on.exit is the backstop)
