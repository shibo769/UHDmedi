library(dplyr)
library(readr)
library(stringr)
library(ggplot2)
library(tidyr)
library(purrr)

# ============================================================
# 0. User settings
# ============================================================

base_dir <- "/usr3/graduate/shibo/UHDmediation/Simulation/Revision"

pq_vec <- c(50, 400, 750, 1250)
n_vec  <- c(250, 500, 750, 1000, 1250)

sigma_eps <- 1

s_map <- list(
  "50"   = list(s = 5,  s_B = 5),
  "400"  = list(s = 10, s_B = 10),
  "750"  = list(s = 10, s_B = 10),
  "1250" = list(s = 20, s_B = 20)
)

# ============================================================
# 1. Helper functions
# ============================================================

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

compute_one_summary <- function(dat, file_name, method, n, p, q, s, s_B) {
  
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
  } else if (method == "Guo") {
    est_col <- find_col(
      dat,
      candidates = c("NIE_hat", "NIE_hat_Guo", "NIE_hat_SCAD", 
                     "NIE_hat_naive", "theta_hat"),
      file_name = file_name
    )
  } else {
    stop("Unknown method.")
  }
  
  NIE_hat  <- dat[[est_col]]
  NIE_true <- dat[[true_col]][1]
  
  bias_NIE <- NIE_hat - NIE_true
  
  bias_mean_NIE <- mean(bias_NIE, na.rm = TRUE)
  rmse_NIE      <- sqrt(mean(bias_NIE^2, na.rm = TRUE))
  sd_NIE        <- sd(bias_NIE, na.rm = TRUE)
  
  ci_lower_NIE <- NIE_hat - 1.96 * sd_NIE
  ci_upper_NIE <- NIE_hat + 1.96 * sd_NIE
  
  coverage_NIE <- mean(
    ci_lower_NIE <= NIE_true & NIE_true <= ci_upper_NIE,
    na.rm = TRUE
  )
  
  ci_length_NIE <- mean(ci_upper_NIE - ci_lower_NIE, na.rm = TRUE)
  
  data.frame(
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
    coverage_NIE = round(coverage_NIE, 3)
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
  } else if (method == "Guo") {
    est_col <- find_col(
      dat,
      candidates = c("NIE_hat", "NIE_hat_Guo", "NIE_hat_SCAD", 
                     "NIE_hat_naive", "theta_hat"),
      file_name = file_name
    )
  }
  
  NIE_hat  <- dat[[est_col]]
  NIE_true <- dat[[true_col]][1]
  
  data.frame(
    method = method,
    n = n,
    p = p,
    q = q,
    value = sqrt(n) * (NIE_hat - NIE_true)
  )
}

# ============================================================
# 2. Read all files, compute two result tables
# ============================================================

results_debiasing <- data.frame()
results_guo       <- data.frame()
hist_all          <- data.frame()

for (q in pq_vec) {
  
  p <- q
  
  s   <- s_map[[as.character(q)]]$s
  s_B <- s_map[[as.character(q)]]$s_B
  
  for (nn in n_vec) {
    
    debias_file <- file.path(
      base_dir,
      paste0("mc_results_pq_", q, "_n_", nn, ".csv")
    )
    
    guo_file <- file.path(
      base_dir,
      paste0("MC_iter_results_SCAD_q_", q, "_n_", nn, ".csv")
    )
    
    # ----------------------------
    # Debiasing / new method
    # ----------------------------
    if (file.exists(debias_file)) {
      
      dat_debias <- read_csv(debias_file, show_col_types = FALSE)
      
      results_debiasing <- rbind(
        results_debiasing,
        compute_one_summary(
          dat = dat_debias,
          file_name = debias_file,
          method = "Debiasing",
          n = nn,
          p = p,
          q = q,
          s = s,
          s_B = s_B
        )
      )
      
      hist_all <- rbind(
        hist_all,
        make_hist_data(
          dat = dat_debias,
          file_name = debias_file,
          method = "Debiasing",
          n = nn,
          p = p,
          q = q
        )
      )
      
    } else {
      warning("Missing debiasing file: ", debias_file)
    }
    
    # ----------------------------
    # Guo / SCAD method
    # ----------------------------
    if (file.exists(guo_file)) {
      
      dat_guo <- read_csv(guo_file, show_col_types = FALSE)
      
      results_guo <- rbind(
        results_guo,
        compute_one_summary(
          dat = dat_guo,
          file_name = guo_file,
          method = "Guo",
          n = nn,
          p = p,
          q = q,
          s = s,
          s_B = s_B
        )
      )
      
      hist_all <- rbind(
        hist_all,
        make_hist_data(
          dat = dat_guo,
          file_name = guo_file,
          method = "Guo",
          n = nn,
          p = p,
          q = q
        )
      )
      
    } else {
      warning("Missing Guo file: ", guo_file)
    }
  }
}

# ============================================================
# 3. Print two tables
# ============================================================

results_debiasing <- results_debiasing %>%
  arrange(q, n)

results_guo <- results_guo %>%
  arrange(q, n)

cat("\n================ Debiasing / New Method ================\n")
print(results_debiasing)

cat("\n================ Guo / SCAD Method ================\n")
print(results_guo)

# ============================================================
# 4. Histogram: print one figure for each (q, n)
#    x-axis: sqrt(n)(NIE_hat - NIE_true)
# ============================================================

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
    
    p_now <- ggplot(
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
        values = c("Debiasing" = "blue", "Guo" = "red")
      ) +
      scale_color_manual(
        values = c("Debiasing" = "blue", "Guo" = "red")
      ) +
      labs(
        title = expression("Histogram of " * sqrt(n) * "(" * widehat(NIE) - NIE * ")"),
        subtitle = paste0(
          "Dimension: p = ", q_now,
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
    
    print(p_now)
  }
}

