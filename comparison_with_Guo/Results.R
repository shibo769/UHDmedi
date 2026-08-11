############################################################
# Results.R
#
# Combine the no-interaction debiasing method (NoInteraction.R) and the
# Guo/SCAD method (Guo2022.R) into comparison tables and comparison plots,
# for ONE (DGP_NAME, sigma) at a time.
#
# Reads (per (pq, n)):
#   simulation_outputs_debiasing_no_interaction/<DGP>/<sigma_tag>/
#       mc_results_debiasing_nointeraction_<DGP>_<sigma_tag>_pq_<q>_n_<n>.csv
#   simulation_outputs_guo_no_interaction/<DGP>/<sigma_tag>/
#       mc_results_guo_nointeraction_<DGP>_<sigma_tag>_q_<q>_n_<n>.csv
#
# Writes to:
#   comparison_outputs_no_interaction/<DGP>/<sigma_tag>/
#
# Env-controlled: DGP_NAME, MC_SIGMA, MC_PQ, MC_N
############################################################

rm(list = ls())

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

SHARED <- file.path(base_dir, "..", "shared")   # bundled deps
source(file.path(SHARED, "mc_common.R"))

################################################################################
### Package handling (check only)
################################################################################

required_pkgs <- c("dplyr", "readr", "stringr", "ggplot2", "tidyr", "purrr", "gridExtra")
require_packages(required_pkgs)

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(tidyr)
  library(purrr)
  library(gridExtra)
})

# ============================================================
# 0. Settings (env-controlled). Do NOT mix DGPs or sigmas.
# ============================================================

DGP_NAME  <- env_chr("DGP_NAME", "base")
sigma_eps <- env_num("MC_SIGMA", 1)
pq_vec    <- env_int_vec("MC_PQ", "50,400,750,1250")
n_vec     <- env_int_vec("MC_N", "250,500,750,1000,1250")
run_plots <- env_lgl("RUN_PLOTS", TRUE)

if (!(DGP_NAME %in% c("base", "increasing_sparsity", "minimal_signal"))) {
  stop("DGP_NAME must be one of: base, increasing_sparsity, minimal_signal")
}

sigma_tag <- make_sigma_tag(sigma_eps)

debias_tag <- paste0("debiasing_nointeraction_", DGP_NAME, "_", sigma_tag)
guo_tag    <- paste0("guo_nointeraction_", DGP_NAME, "_", sigma_tag)

debias_dir <- file.path(base_dir, "simulation_outputs_debiasing_no_interaction", DGP_NAME, sigma_tag)
guo_dir    <- file.path(base_dir, "simulation_outputs_guo_no_interaction", DGP_NAME, sigma_tag)

compare_out_dir <- file.path(base_dir, "comparison_outputs_no_interaction", DGP_NAME, sigma_tag)
ensure_dir(compare_out_dir)

cat("Comparison for DGP:", DGP_NAME, "| sigma:", sigma_eps, "\n")
cat("  debiasing dir:", debias_dir, "\n")
cat("  guo dir      :", guo_dir, "\n")
cat("  output dir   :", compare_out_dir, "\n")

s_map <- list(
  "50"   = list(s = 5,  s_B = 5),
  "400"  = list(s = 10, s_B = 10),
  "750"  = list(s = 10, s_B = 10),
  "1250" = list(s = 20, s_B = 20)
)

s_base   <- 3
s_B_base <- 3

get_s_values <- function(pq, dgp_name) {
  if (dgp_name == "base") return(list(s = s_base, s_B = s_B_base))
  s_map[[as.character(pq)]]
}

# ============================================================
# 1. Helper functions
# ============================================================

find_existing_file <- function(candidates) {
  hit <- candidates[file.exists(candidates)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

find_col <- function(dat, candidates, file_name = "") {
  hit <- intersect(candidates, names(dat))
  if (length(hit) == 0) {
    stop(
      paste0(
        "Cannot find any of these columns in file: ", file_name, "\n",
        "Candidates: ", paste(candidates, collapse = ", "), "\n",
        "Available columns: ", paste(names(dat), collapse = ", ")
      )
    )
  }
  hit[1]
}

find_col_optional <- function(dat, candidates) {
  hit <- intersect(candidates, names(dat))
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

compute_one_summary <- function(dat, file_name, method, n, p, q, s, s_B, ci_type = c("MC", "Analytical")) {

  ci_type <- match.arg(ci_type)

  true_col <- find_col(
    dat,
    candidates = c("NIE_true", "true_NIE", "theta_true", "theta0"),
    file_name = file_name
  )

  if (method == "Debiasing") {
    est_col <- find_col(
      dat,
      candidates = c("NIE_hat_cf", "NIE_hat_debias", "NIE_hat_debiased",
                     "NIE_hat", "theta_hat"),
      file_name = file_name
    )
    se_col <- find_col_optional(
      dat,
      candidates = c("se_NIE_cf", "se_NIE_hat", "se_analytical", "NIE_se")
    )
  } else if (method == "Guo") {
    est_col <- find_col(
      dat,
      candidates = c("NIE_hat", "NIE_hat_Guo", "NIE_hat_SCAD", "theta_hat"),
      file_name = file_name
    )
    se_col <- find_col_optional(
      dat,
      candidates = c("se_NIE_hat", "se_NIE", "SE_NIE", "NIE_se")
    )
  } else {
    stop("Unknown method.")
  }

  NIE_hat  <- dat[[est_col]]
  NIE_true <- dat[[true_col]][1]

  bias_NIE <- NIE_hat - NIE_true

  bias_mean_NIE <- mean(bias_NIE, na.rm = TRUE)
  rmse_NIE      <- sqrt(mean(bias_NIE^2, na.rm = TRUE))
  sd_NIE        <- sd(NIE_hat, na.rm = TRUE)

  if (ci_type == "MC") {
    se_used <- rep(sd_NIE, length(NIE_hat))
  } else {
    if (is.na(se_col)) {
      warning("No analytical SE column found for ", method, " in file: ", file_name)
      se_used <- rep(NA_real_, length(NIE_hat))
    } else {
      se_used <- dat[[se_col]]
    }
  }

  ci_lower_NIE <- NIE_hat - 1.96 * se_used
  ci_upper_NIE <- NIE_hat + 1.96 * se_used

  coverage_NIE <- mean(
    ci_lower_NIE <= NIE_true & NIE_true <= ci_upper_NIE,
    na.rm = TRUE
  )

  ci_length_NIE <- mean(ci_upper_NIE - ci_lower_NIE, na.rm = TRUE)

  data.frame(
    DGP_NAME = DGP_NAME,
    method = method,
    ci_type = ci_type,
    n = n,
    p = p,
    q = q,
    s = s,
    s_B = s_B,
    Var = sigma_eps,
    NIE_true = NIE_true,
    bias_mean_NIE = round(bias_mean_NIE, 3),
    rmse_NIE = round(rmse_NIE, 3),
    sd_NIE = round(sd_NIE, 3),
    ci_length_NIE = round(ci_length_NIE, 3),
    coverage_NIE = round(coverage_NIE, 3),
    mean_se_NIE = round(mean(se_used, na.rm = TRUE), 4)
  )
}

make_hist_data <- function(dat, file_name, method, n, p, q) {

  true_col <- find_col(
    dat,
    candidates = c("NIE_true", "true_NIE", "theta_true", "theta0"),
    file_name = file_name
  )

  if (method == "Debiasing") {
    est_col <- find_col(
      dat,
      candidates = c("NIE_hat_cf", "NIE_hat_debias", "NIE_hat_debiased",
                     "NIE_hat", "theta_hat"),
      file_name = file_name
    )
    se_col <- find_col_optional(
      dat,
      candidates = c("se_NIE_cf", "se_NIE_hat", "se_analytical", "NIE_se")
    )
  } else if (method == "Guo") {
    est_col <- find_col(
      dat,
      candidates = c("NIE_hat", "NIE_hat_Guo", "NIE_hat_SCAD", "theta_hat"),
      file_name = file_name
    )
    se_col <- find_col_optional(
      dat,
      candidates = c("se_NIE_hat", "se_NIE", "SE_NIE", "NIE_se")
    )
  }

  NIE_hat  <- dat[[est_col]]
  NIE_true <- dat[[true_col]][1]

  z_value <- rep(NA_real_, length(NIE_hat))
  if (!is.na(se_col)) {
    z_value <- (NIE_hat - NIE_true) / dat[[se_col]]
  }

  data.frame(
    method = method,
    n = n,
    p = p,
    q = q,
    value = sqrt(n) * (NIE_hat - NIE_true),
    z_value = z_value
  )
}

# ============================================================
# 2. Load the full per-iteration RDS objects and compute result tables
#
# The simulation scripts no longer write per-(q,n) CSVs; the raw iterations
# live in the integrated RDS objects. We read each (q,n) segment from there.
#   debiasing RDS keys: "pq_<q>_n_<n>"
#   guo RDS keys:       "q_<q>_n_<n>"
# ============================================================

debias_rds <- file.path(debias_dir, paste0("iter_outputs_", debias_tag, "_res.rds"))
guo_rds    <- file.path(guo_dir,    paste0("iter_outputs_", guo_tag,    "_res.rds"))

debias_list <- if (file.exists(debias_rds)) readRDS(debias_rds) else NULL
guo_list    <- if (file.exists(guo_rds))    readRDS(guo_rds)    else NULL

if (is.null(debias_list)) warning("Missing debiasing RDS: ", debias_rds)
if (is.null(guo_list))    warning("Missing Guo RDS: ", guo_rds)

results_mc_all  <- data.frame()
results_var_all <- data.frame()
hist_all        <- data.frame()

for (q in pq_vec) {

  p <- q

  s_info <- get_s_values(q, DGP_NAME)
  s   <- s_info$s
  s_B <- s_info$s_B

  for (nn in n_vec) {

    dat_debias <- if (!is.null(debias_list)) debias_list[[paste0("pq_", q, "_n_", nn)]] else NULL
    dat_guo    <- if (!is.null(guo_list))    guo_list[[paste0("q_", q, "_n_", nn)]]     else NULL

    # ----------------------------
    # Debiasing / new method
    # ----------------------------
    if (!is.null(dat_debias)) {

      dbg <- paste0("debiasing RDS segment pq_", q, "_n_", nn)

      results_mc_all <- rbind(
        results_mc_all,
        compute_one_summary(
          dat = dat_debias, file_name = dbg, method = "Debiasing",
          n = nn, p = p, q = q, s = s, s_B = s_B, ci_type = "MC"
        )
      )

      results_var_all <- rbind(
        results_var_all,
        compute_one_summary(
          dat = dat_debias, file_name = dbg, method = "Debiasing",
          n = nn, p = p, q = q, s = s, s_B = s_B, ci_type = "Analytical"
        )
      )

      hist_all <- rbind(
        hist_all,
        make_hist_data(
          dat = dat_debias, file_name = dbg, method = "Debiasing",
          n = nn, p = p, q = q
        )
      )

    } else {
      warning("Missing debiasing segment for q = ", q, ", n = ", nn)
    }

    # ----------------------------
    # Guo / SCAD method
    # ----------------------------
    if (!is.null(dat_guo)) {

      gbg <- paste0("guo RDS segment q_", q, "_n_", nn)

      results_mc_all <- rbind(
        results_mc_all,
        compute_one_summary(
          dat = dat_guo, file_name = gbg, method = "Guo",
          n = nn, p = p, q = q, s = s, s_B = s_B, ci_type = "MC"
        )
      )

      results_var_all <- rbind(
        results_var_all,
        compute_one_summary(
          dat = dat_guo, file_name = gbg, method = "Guo",
          n = nn, p = p, q = q, s = s, s_B = s_B, ci_type = "Analytical"
        )
      )

      hist_all <- rbind(
        hist_all,
        make_hist_data(
          dat = dat_guo, file_name = gbg, method = "Guo",
          n = nn, p = p, q = q
        )
      )

    } else {
      warning("Missing Guo segment for q = ", q, ", n = ", nn)
    }
  }
}

# ============================================================
# 3. Print and save two comparison tables
# ============================================================

if (nrow(results_mc_all) > 0) {
  results_mc_all <- results_mc_all %>% arrange(q, n, method)
}
if (nrow(results_var_all) > 0) {
  results_var_all <- results_var_all %>% arrange(q, n, method)
}

cat("\n================ Comparison Table: Monte Carlo CI ================\n")
print(results_mc_all)

cat("\n================ Comparison Table: Analytical Variance CI ================\n")
print(results_var_all)

write_csv(results_mc_all,  file.path(compare_out_dir, "comparison_table_MC_CI.csv"))
write_csv(results_var_all, file.path(compare_out_dir, "comparison_table_Analytical_CI.csv"))
saveRDS(results_mc_all,   file.path(compare_out_dir, "comparison_table_MC_CI.rds"))
saveRDS(results_var_all,  file.path(compare_out_dir, "comparison_table_Analytical_CI.rds"))

# ============================================================
# 4. Histogram: one figure for each (q, n)
#    x-axis: sqrt(n)(NIE_hat - NIE_true)
# ============================================================

if (run_plots && nrow(hist_all) > 0) {

hist_all <- hist_all %>%
  filter(!is.na(value), is.finite(value)) %>%
  mutate(
    method = factor(method, levels = c("Debiasing", "Guo"))
  )

for (q_now in pq_vec) {

  for (n_now in n_vec) {

    plot_dat <- hist_all %>%
      filter(q == q_now, n == n_now)

    if (nrow(plot_dat) == 0) {
      warning("No histogram data for q = ", q_now, ", n = ", n_now)
      next
    }

    p_hist <- ggplot(
      plot_dat,
      aes(x = value, fill = method, color = method)
    ) +
      geom_histogram(
        aes(y = after_stat(density)),
        bins = 80,
        alpha = 0.35,
        position = "identity"
      ) +
      scale_fill_manual(
        values = c("Debiasing" = "blue", "Guo" = "red"),
        labels = c("Debiasing" = "Ours", "Guo" = "Guo et al.")
      ) +
      scale_color_manual(
        values = c("Debiasing" = "blue", "Guo" = "red"),
        labels = c("Debiasing" = "Ours", "Guo" = "Guo et al.")
      ) +
      labs(
        title = expression("Histogram of " * sqrt(n) * "(" * widehat(NIE) - NIE * ")"),
        subtitle = paste0(
          "DGP: ", DGP_NAME,
          "; Dimension: p = ", q_now,
          ", q = ", q_now,
          "; Sample size: n = ", n_now
        ),
        x = expression(sqrt(n) * "(" * widehat(NIE) - NIE * ")"),
        y = "Density",
        fill = "Method",
        color = "Method"
      ) +
      theme_bw(base_size = 13) +
      theme(
        panel.grid.major = element_line(color = "grey80", linewidth = 0.35),
        panel.grid.minor = element_line(color = "grey90", linewidth = 0.25),
        legend.position = "top",
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 11)
      )

    ggsave(
      filename = file.path(compare_out_dir, paste0("hist_compare_q_", q_now, "_n_", n_now, ".png")),
      plot = p_hist,
      width = 8,
      height = 5,
      dpi = 300
    )
  }
}

# ============================================================
# 5. Analytical standardized QQ plots:
#    z = (NIE_hat - NIE_true) / se_hat
# ============================================================

z_all <- hist_all %>%
  filter(!is.na(z_value), is.finite(z_value)) %>%
  mutate(
    method = factor(method, levels = c("Debiasing", "Guo"))
  )

for (q_now in pq_vec) {

  for (n_now in n_vec) {

    plot_dat <- z_all %>%
      filter(q == q_now, n == n_now)

    if (nrow(plot_dat) == 0) {
      warning("No analytical QQ data for q = ", q_now, ", n = ", n_now)
      next
    }

    p_qq <- ggplot(
      plot_dat,
      aes(sample = z_value, color = method)
    ) +
      stat_qq(alpha = 0.8) +
      stat_qq_line(color = "black", linetype = "dashed") +
      scale_color_manual(
        values = c("Debiasing" = "blue", "Guo" = "red"),
        labels = c("Debiasing" = "Ours", "Guo" = "Guo et al.")
      ) +
      labs(
        title = "QQ Plot of Analytical Standardized Errors",
        subtitle = paste0(
          "DGP: ", DGP_NAME,
          "; Dimension: p = ", q_now,
          ", q = ", q_now,
          "; Sample size: n = ", n_now
        ),
        x = "Theoretical Normal Quantiles",
        y = expression((widehat(NIE) - NIE) / widehat(se)),
        color = "Method"
      ) +
      theme_bw(base_size = 13) +
      theme(
        panel.grid.major = element_line(color = "grey80", linewidth = 0.35),
        panel.grid.minor = element_line(color = "grey90", linewidth = 0.25),
        legend.position = "top",
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 11)
      )

    ggsave(
      filename = file.path(compare_out_dir, paste0("qq_compare_analytic_q_", q_now, "_n_", n_now, ".png")),
      plot = p_qq,
      width = 8,
      height = 5,
      dpi = 300
    )
  }
}

# ============================================================
# 6. One combined figure for the best setting
# ============================================================

best_setting <- results_var_all %>%
  filter(method == "Debiasing", is.finite(coverage_NIE), is.finite(rmse_NIE)) %>%
  mutate(score = abs(coverage_NIE - 0.95) + 0.2 * rmse_NIE) %>%
  arrange(score) %>%
  slice(1)

if (nrow(best_setting) > 0) {

  best_q <- best_setting$q[1]
  best_n <- best_setting$n[1]

  best_hist <- hist_all %>% filter(q == best_q, n == best_n)
  best_z    <- z_all    %>% filter(q == best_q, n == best_n)

  p_best_hist <- ggplot(best_hist, aes(x = value, fill = method, color = method)) +
    geom_histogram(
      aes(y = after_stat(density)),
      bins = 80,
      alpha = 0.35,
      position = "identity"
    ) +
    scale_fill_manual(values = c("Debiasing" = "blue", "Guo" = "red"),
                      labels = c("Debiasing" = "Ours", "Guo" = "Guo et al.")) +
    scale_color_manual(values = c("Debiasing" = "blue", "Guo" = "red"),
                       labels = c("Debiasing" = "Ours", "Guo" = "Guo et al.")) +
    labs(
      title = expression("Histogram of " * sqrt(n) * "(" * widehat(NIE) - NIE * ")"),
      subtitle = paste0("DGP: ", DGP_NAME, "; p = q = ", best_q, ", n = ", best_n),
      x = expression(sqrt(n) * "(" * widehat(NIE) - NIE * ")"),
      y = "Density",
      fill = "Method",
      color = "Method"
    ) +
    theme_bw(base_size = 13) +
    theme(legend.position = "top")

  p_best_qq <- ggplot(best_z, aes(sample = z_value, color = method)) +
    stat_qq(alpha = 0.8) +
    stat_qq_line(color = "black", linetype = "dashed") +
    scale_color_manual(values = c("Debiasing" = "blue", "Guo" = "red"),
                       labels = c("Debiasing" = "Ours", "Guo" = "Guo et al.")) +
    labs(
      title = "QQ Plot of Analytical Standardized Errors",
      subtitle = paste0("DGP: ", DGP_NAME, "; p = q = ", best_q, ", n = ", best_n),
      x = "Theoretical Normal Quantiles",
      y = expression((widehat(NIE) - NIE) / widehat(se)),
      color = "Method"
    ) +
    theme_bw(base_size = 13) +
    theme(legend.position = "top")

  combined_best <- gridExtra::arrangeGrob(p_best_hist, p_best_qq, ncol = 2)

  ggsave(
    filename = file.path(compare_out_dir, paste0("best_hist_qq_compare_q_", best_q, "_n_", best_n, ".png")),
    plot = combined_best,
    width = 13,
    height = 5,
    dpi = 300
  )
}

}  ### end if (run_plots)

cat("\nComparison outputs saved to:", compare_out_dir, "\n")
