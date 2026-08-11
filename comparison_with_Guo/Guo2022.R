############################################################
# NIE via TE - NDE using SCAD + hand-written HBIC
# No-interaction comparison with modified Guo et al. method
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
#   MC_ITER=2 MC_PQ=50 MC_N=250 MC_CORES=4 DGP_NAME=base Rscript Guo2022.R
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
## Guo helper functions (HBIC_calc, mediationInference, safe_solve, deSCAD, LLA_h1).
source(file.path(SHARED, "utils_mediation.R"))

################################################################################
### Package handling (check only)
################################################################################

required_pkgs <- c(
  "foreach", "doParallel", "doRNG", "pracma", "MASS",
  "glmnet", "ncvreg", "ggplot2", "gridExtra", "dplyr"
)
require_packages(required_pkgs)

suppressPackageStartupMessages({
  library(foreach)
  library(doParallel)
  library(doRNG)
  library(pracma)
  library(MASS)
  library(glmnet)
  library(ncvreg)
  library(ggplot2)
  library(gridExtra)
  library(dplyr)
})

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
SCAD_a   <- 3.7

if (!(DGP_NAME %in% c("base", "increasing_sparsity", "minimal_signal"))) {
  stop("DGP_NAME must be one of: base, increasing_sparsity, minimal_signal")
}

set.seed(seed)

################################################################################
### Output folder (Guo no-interaction outputs live under Comparison/)
################################################################################

sigma_tag <- make_sigma_tag(sigma_eps)
file_tag  <- paste0("guo_nointeraction_", DGP_NAME, "_", sigma_tag)

output_root <- file.path(base_dir, "simulation_outputs_guo_no_interaction")
output_dir  <- file.path(output_root, DGP_NAME, sigma_tag)
ensure_dir(output_dir)

cat("No-interaction Guo DGP version:", DGP_NAME, "\n")
cat("Iter:", Iter, "| pq:", paste(pq_vec, collapse = ","),
    "| n:", paste(n_vec, collapse = ","), "\n")
cat("Sigma:", sigma_eps, "| seed:", seed, "\n")
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

################################################################################
### Result column template (keeps failed replications rbind-compatible)
################################################################################

result_cols <- c(
  "iter", "n", "p", "q", "s", "s_B", "Var",
  "TE_true", "NDE_true", "NIE_true",
  "TE_hat", "NDE_hat", "NIE_hat", "bias_NIE",
  "se_NIE_hat", "var_NIE_hat",
  "lambda_best", "hbic_best", "active_M_size", "active_S_size",
  "failed"
)

################### Simulation ###################

res <- list()

for (q in pq_vec) {

  p <- q

  s_info <- get_s_values(q, DGP_NAME)
  s   <- s_info$s
  s_B <- s_info$s_B

  cat("\n============================\n")
  cat("Starting q =", q, "| s =", s, "| s_B =", s_B, "\n")
  cat("============================\n")

  ##########################################################################
  ### True parameters
  ##########################################################################

  beta <- rep(0, p)
  beta_index <- sample(1:p, s, FALSE)
  beta[beta_index] <- draw_signal_nointeraction(s, DGP_NAME)

  idx_gamma  <- sample(1:q, s, FALSE)
  idx_shared <- sample(idx_gamma, floor(s / 2), FALSE)
  idx_delta  <- c(
    idx_shared,
    sample(setdiff(1:q, idx_gamma), s - length(idx_shared), FALSE)
  )

  gamma <- rep(0, q)
  gamma[idx_gamma] <- draw_signal_nointeraction(s, DGP_NAME)

  delta1 <- rep(0, q)
  delta1[idx_delta] <- draw_signal_nointeraction(s, DGP_NAME)

  alpha1_true <- 0.15

  TE_true  <- alpha1_true + sum(delta1 * gamma)
  NDE_true <- alpha1_true
  NIE_true <- sum(delta1 * gamma)

  B <- matrix(0, q, p)
  for (iB in 1:q) {
    B[iB, sample(1:p, s_B, FALSE)] <- runif(s_B, 0.05, 0.2)
  }

  for (nn in n_vec) {

    ########################################################################
    ### Generate X and A once for this setting
    ########################################################################

    mu_vec <- rep(mu, p)
    U_mat <- randortho(p)
    Sigma_X <- U_mat %*% diag(runif(p, 0.1, 0.5)) %*% t(U_mat)
    X <- mvrnorm(nn, mu_vec, Sigma_X)

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
      warning("Skipping q = ", q, ", n = ", nn, " because A has no variation.")
      next
    }

    ########################################################################
    ### Monte Carlo loop
    ########################################################################

    err_vec <- foreach(
      i = 1:Iter,
      .combine = "rbind",
      .packages = c("MASS", "pracma", "glmnet", "ncvreg"),
      .export = c(
        "HBIC_calc", "mediationInference", "safe_solve", "deSCAD", "LLA_h1",
        "na_result_row", "result_cols", "worker_single_thread_blas"
      ),
      .options.multicore = mc_foreach_options_multicore,
      .options.snow = mc_foreach_options_snow
    ) %dopar% {

      worker_single_thread_blas()

      one_rep <- tryCatch({

        lambda_U <- runif(q, 0.1, 0.5)
        R_U <- randortho(q)
        Sigma_U <- R_U %*% diag(lambda_U) %*% t(R_U)
        U <- mvrnorm(nn, rep(0, q), Sigma_U)

        A_mat <- matrix(A, ncol = 1)

        M <- A_mat %*% t(delta1) + X %*% t(B) + U
        Y <- as.vector(alpha1_true * A + X %*% beta + M %*% gamma + rnorm(nn, 0, sigma_eps))

        ######################################################################
        ### Step 1: SCAD + HBIC selection
        ######################################################################

        lamb_grid <- seq(0.01, 0.25, length.out = 20)

        hbic_fit <- lapply(
          lamb_grid,
          HBIC_calc,
          xx = as.matrix(A),
          yy = Y,
          mm = M,
          S = X,
          n_imp = 0,
          penalize_S = TRUE
        )

        hbic_vals <- sapply(hbic_fit, function(z) z$BIC)
        id_best <- tail(which(hbic_vals == min(hbic_vals)), 1)
        result_best <- hbic_fit[[id_best]]

        active_M <- which(result_best$alpha0 != 0)
        active_S <- which(result_best$alpha2 != 0)

        M_refit <- if (length(active_M) == 0) {
          matrix(numeric(0), nrow = nn, ncol = 0)
        } else {
          M[, active_M, drop = FALSE]
        }

        S_refit <- if (length(active_S) == 0) {
          matrix(numeric(0), nrow = nn, ncol = 0)
        } else {
          X[, active_S, drop = FALSE]
        }

        ######################################################################
        ### Step 2: low-dimensional refit/inference after selection
        ######################################################################

        infer <- mediationInference(
          X = as.matrix(A),
          Y = Y,
          M = M_refit,
          S = S_refit
        )

        NIE_hat <- as.numeric(infer$beta_hat[1])
        NDE_hat <- as.numeric(infer$alpha1_hat[1])
        TE_hat <- NIE_hat + NDE_hat

        bias_NIE <- NIE_hat - NIE_true

        ######################################################################
        ### Analytical variance from Guo-style post-selection inference
        ######################################################################

        se_NIE_hat <- as.numeric(infer$summary_result$SE[infer$summary_result$Coefficient == "beta"][1])

        if (!is.finite(se_NIE_hat)) {
          se_NIE_hat <- sqrt(pmax(as.numeric(infer$var_beta[1, 1]) / nn, 0))
        }

        var_NIE_hat <- se_NIE_hat^2

        c(
          iter = i,
          n = nn,
          p = p,
          q = q,
          s = s,
          s_B = s_B,
          Var = sigma_eps,

          TE_true = TE_true,
          NDE_true = NDE_true,
          NIE_true = NIE_true,

          TE_hat = TE_hat,
          NDE_hat = NDE_hat,
          NIE_hat = NIE_hat,
          bias_NIE = bias_NIE,

          se_NIE_hat = se_NIE_hat,
          var_NIE_hat = var_NIE_hat,
          lambda_best = lamb_grid[id_best],
          hbic_best = hbic_vals[id_best],
          active_M_size = length(active_M),
          active_S_size = length(active_S),
          failed = 0
        )
      }, error = function(e) {
        row <- na_result_row(result_cols)
        row["iter"] <- i
        row["n"] <- nn
        row["p"] <- p
        row["q"] <- q
        row["s"] <- s
        row["s_B"] <- s_B
        row["Var"] <- sigma_eps
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

    ## Kept in memory + full RDS. No per-(q,n) CSV; only integrated summaries.
    err_df <- as.data.frame(err_vec)
    err_df$DGP_NAME <- DGP_NAME
    err_df$seed <- seed
    err_df$Iter <- Iter

    res[[paste0("q_", q, "_n_", nn)]] <- err_df

    cat(
      "Finished q =", q,
      "n =", nn,
      "| Mean Bias =", round(mean(err_df$bias_NIE, na.rm = TRUE), 4),
      "| RMSE =", round(sqrt(mean(err_df$bias_NIE^2, na.rm = TRUE)), 4),
      "| MC SD =", round(sd(err_df$NIE_hat, na.rm = TRUE), 4),
      "| Mean analytical SE =", round(mean(err_df$se_NIE_hat, na.rm = TRUE), 4),
      "| Failed =", sum(err_df$failed, na.rm = TRUE),
      "\n"
    )
  }
}

################################################################################
### Merge with any previous stage's RDS (two-stage low/high-pq runs), then save.
################################################################################

rds_path <- file.path(output_dir, paste0("iter_outputs_", file_tag, "_res.rds"))
if (file.exists(rds_path)) {
  prev <- tryCatch(readRDS(rds_path), error = function(e) NULL)
  if (is.list(prev)) {
    for (k in names(prev)) if (is.null(res[[k]])) res[[k]] <- prev[[k]]
  }
}
saveRDS(res, file = rds_path)

## Rebuild summaries over ALL q present in the merged object.
pq_vec <- sort(unique(as.integer(sub("^q_([0-9]+)_n_.*$", "\\1", names(res)))))

################### Summary Table 1: Monte Carlo CI ###################

results_mc_df <- data.frame()

for (q in pq_vec) {

  p <- q

  s_info <- get_s_values(q, DGP_NAME)
  s   <- s_info$s
  s_B <- s_info$s_B

  for (nn in n_vec) {

    key <- paste0("q_", q, "_n_", nn)

    if (!is.null(res[[key]])) {

      segment <- res[[key]]

      NIE_hat  <- segment$NIE_hat
      NIE_true_setting <- segment$NIE_true[1]

      bias_NIE <- NIE_hat - NIE_true_setting

      bias_mean_NIE <- mean(bias_NIE, na.rm = TRUE)
      rmse_NIE <- sqrt(mean(bias_NIE^2, na.rm = TRUE))
      sd_NIE <- sd(NIE_hat, na.rm = TRUE)

      ci_lower_NIE <- NIE_hat - 1.96 * sd_NIE
      ci_upper_NIE <- NIE_hat + 1.96 * sd_NIE

      coverage_NIE <- mean(
        ci_lower_NIE <= NIE_true_setting &
          NIE_true_setting <= ci_upper_NIE,
        na.rm = TRUE
      )

      ci_length_NIE <- mean(
        ci_upper_NIE - ci_lower_NIE,
        na.rm = TRUE
      )

      results_mc_df <- rbind(
        results_mc_df,
        data.frame(
          DGP_NAME = DGP_NAME,
          seed = seed,
          Iter = Iter,
          n = nn,
          p = p,
          q = q,
          s = s,
          s_B = s_B,
          Var = sigma_eps,

          NIE_true = NIE_true_setting,

          bias_mean_NIE = round(bias_mean_NIE, 3),
          rmse_NIE = round(rmse_NIE, 3),
          sd_NIE = round(sd_NIE, 3),
          ci_lower_NIE = round(mean(ci_lower_NIE, na.rm = TRUE), 3),
          ci_upper_NIE = round(mean(ci_upper_NIE, na.rm = TRUE), 3),
          ci_length_NIE = round(ci_length_NIE, 3),
          coverage_NIE = round(coverage_NIE, 3),

          active_M_mean = round(mean(segment$active_M_size, na.rm = TRUE), 2),
          active_S_mean = round(mean(segment$active_S_size, na.rm = TRUE), 2),
          n_valid = sum(!is.na(NIE_hat)),
          n_failed = sum(segment$failed, na.rm = TRUE)
        )
      )
    }
  }
}

cat("\n================ Guo Monte Carlo CI Table ================\n")
print(results_mc_df)

saveRDS(
  results_mc_df,
  file = file.path(output_dir, paste0("table_", file_tag, "_monte_carlo_ci.rds"))
)

write.csv(
  results_mc_df,
  file = file.path(output_dir, paste0("table_", file_tag, "_monte_carlo_ci.csv")),
  row.names = FALSE
)

################### Summary Table 2: Analytical Variance CI ###################

results_var_df <- data.frame()

for (q in pq_vec) {

  p <- q

  s_info <- get_s_values(q, DGP_NAME)
  s   <- s_info$s
  s_B <- s_info$s_B

  for (nn in n_vec) {

    key <- paste0("q_", q, "_n_", nn)

    if (!is.null(res[[key]])) {

      segment <- res[[key]]

      NIE_hat  <- segment$NIE_hat
      NIE_true_setting <- segment$NIE_true[1]
      se_NIE_hat <- segment$se_NIE_hat

      bias_NIE <- NIE_hat - NIE_true_setting

      bias_mean_NIE <- mean(bias_NIE, na.rm = TRUE)
      rmse_NIE <- sqrt(mean(bias_NIE^2, na.rm = TRUE))
      sd_NIE <- sd(NIE_hat, na.rm = TRUE)

      ci_lower_NIE <- NIE_hat - 1.96 * se_NIE_hat
      ci_upper_NIE <- NIE_hat + 1.96 * se_NIE_hat

      coverage_NIE <- mean(
        ci_lower_NIE <= NIE_true_setting &
          NIE_true_setting <= ci_upper_NIE,
        na.rm = TRUE
      )

      ci_length_NIE <- mean(
        ci_upper_NIE - ci_lower_NIE,
        na.rm = TRUE
      )

      results_var_df <- rbind(
        results_var_df,
        data.frame(
          DGP_NAME = DGP_NAME,
          seed = seed,
          Iter = Iter,
          n = nn,
          p = p,
          q = q,
          s = s,
          s_B = s_B,
          Var = sigma_eps,

          NIE_true = NIE_true_setting,

          bias_mean_NIE = round(bias_mean_NIE, 3),
          rmse_NIE = round(rmse_NIE, 3),
          sd_NIE = round(sd_NIE, 3),
          ci_lower_NIE = round(mean(ci_lower_NIE, na.rm = TRUE), 3),
          ci_upper_NIE = round(mean(ci_upper_NIE, na.rm = TRUE), 3),
          ci_length_NIE = round(ci_length_NIE, 3),
          coverage_NIE = round(coverage_NIE, 3),
          mean_se_NIE = round(mean(se_NIE_hat, na.rm = TRUE), 4),

          active_M_mean = round(mean(segment$active_M_size, na.rm = TRUE), 2),
          active_S_mean = round(mean(segment$active_S_size, na.rm = TRUE), 2),
          n_valid = sum(!is.na(NIE_hat)),
          n_failed = sum(segment$failed, na.rm = TRUE)
        )
      )
    }
  }
}

cat("\n================ Guo Analytical Variance CI Table ================\n")
print(results_var_df)

saveRDS(
  results_var_df,
  file = file.path(output_dir, paste0("table_", file_tag, "_analytical_variance_ci.rds"))
)

write.csv(
  results_var_df,
  file = file.path(output_dir, paste0("table_", file_tag, "_analytical_variance_ci.csv")),
  row.names = FALSE
)

######################## Plot: one best setting ########################

if (run_plots && nrow(results_var_df) > 0) {

plot_choice <- results_var_df %>%
  filter(is.finite(coverage_NIE), is.finite(rmse_NIE)) %>%
  mutate(score = abs(coverage_NIE - 0.95) + 0.2 * rmse_NIE) %>%
  arrange(score) %>%
  slice(1)

if (nrow(plot_choice) > 0) {

  plot_q  <- plot_choice$q[1]
  plot_n  <- plot_choice$n[1]
  plot_key <- paste0("q_", plot_q, "_n_", plot_n)

  cat("\nPlotting best combination based on analytical coverage and RMSE:\n")
  print(plot_choice)

  segment_plot <- res[[plot_key]]

  res_scaled <- sqrt(plot_n) * segment_plot$bias_NIE

  res_df <- data.frame(
    x = res_scaled,
    type = "Guo"
  )

  x_lim_best <- c(min(min(res_df$x, na.rm = TRUE), 0), max(max(res_df$x, na.rm = TRUE), 0))
  if (!is.finite(diff(x_lim_best)) || diff(x_lim_best) == 0) {
    x_lim_best <- x_lim_best + c(-1, 1)
  }
  x_breaks_best <- sort(unique(c(pretty(x_lim_best, n = 8), 0)))

  hist_plot <- ggplot(res_df, aes(x = x, fill = type)) +
    geom_histogram(
      bins = 35,
      color = "#e9ecef",
      alpha = 0.8,
      position = "identity"
    ) +
    scale_fill_manual(values = c("Guo" = "tomato")) +
    scale_x_continuous(limits = x_lim_best, breaks = x_breaks_best) +
    theme_minimal() +
    ggtitle(expression(paste("Histogram of ", sqrt(n) * (hat(NIE) - NIE)))) +
    xlab(expression(paste(sqrt(n) * (hat(NIE) - NIE)))) +
    ylab("Frequency")

  z_NIE <- segment_plot$bias_NIE / segment_plot$se_NIE_hat
  z_NIE <- z_NIE[is.finite(z_NIE)]

  if (length(z_NIE) >= 2) {
    qq_Guo <- qqnorm(z_NIE, plot.it = FALSE)
    qq_df <- data.frame(
      theoretical_quantiles = qq_Guo$x,
      sample_quantiles = qq_Guo$y
    )
  } else {
    ## Not enough finite standardized errors to form a QQ plot (e.g. tiny Iter).
    qq_df <- data.frame(theoretical_quantiles = numeric(0), sample_quantiles = numeric(0))
  }

  qq_plot <- ggplot(qq_df, aes(x = theoretical_quantiles, y = sample_quantiles)) +
    geom_point(color = "tomato", alpha = 0.8) +
    theme_minimal() +
    ggtitle("QQ Plot") +
    xlab("Theoretical Normal Quantiles") +
    ylab(expression(paste((hat(NIE) - NIE) / hat(se))))

  combined_plot <- gridExtra::arrangeGrob(hist_plot, qq_plot, ncol = 2)

  ggsave(
    filename = file.path(output_dir, paste0("hist_qq_Guo_best_", file_tag, "_q_", plot_q, "_n_", plot_n, ".png")),
    plot = combined_plot,
    width = 12,
    height = 5,
    dpi = 300
  )

  ggsave(
    filename = file.path(output_dir, paste0("hist_qq_Guo_best_", file_tag, "_q_", plot_q, "_n_", plot_n, ".pdf")),
    plot = combined_plot,
    width = 12,
    height = 5
  )
}

}  ### end if (run_plots)

cat("\nAll Guo no-interaction simulation outputs saved to:", output_dir, "\n")

try(parallel::stopCluster(my.cluster), silent = TRUE)
