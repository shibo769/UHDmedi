#################################
##   UHDmedi Package Example   ##
#################################
library(UHDmedi)

# Generate simulated data
cat("----- 1) Generating Simulated Data -----\n")
sim_data <- simData(nn = 30, p = 50, q = 20, variance = 0.1, k1 = 3, k2 = 6, s_X = 3, categorical = TRUE)
str(sim_data)

# 3. Run basic EDA on the data
cat("\n----- 2) Running EDA -----\n")
eda_results <- eda_report(X = sim_data$X,
                          M = sim_data$M,
                          A = sim_data$A,
                          Y = sim_data$Y,
                          plot_corr = TRUE,
                          plot_Y = TRUE,
                          verbose = TRUE)

# View EDA output
print(eda_results$summary)
# To view plots in an interactive environment, run:
print(eda_results$plots$scatter_plot_minvar)
print(eda_results$plots$scatter_plot_maxvar)
print(eda_results$plots$corr_plot)
print(eda_results$plots$Y_density)
print(eda_results$plots$Y_boxplot)

# 4. Estimate effects (ATE, NIE, NDE)
cat("\n----- 3) Estimating Effects (ATE, NIE, NDE) -----\n")
effects_result <- est_effects(X = sim_data$X,
                              Y = sim_data$Y,
                              M = sim_data$M,
                              A = sim_data$A,
                              logY = FALSE,
                              parallel = TRUE,
                              nfold1 = 5,
                              nfold2 = 5,
                              net_alpha = 1,
                              zeta = 0.5,
                              optimizer = "pogs",
                              intercept = FALSE,
                              verbose = TRUE)

cat("\nEffects result:\n")
print(effects_result)

# 5. Perform non-parametric bootstrap estimation using bootstrap_est_effects
cat("\n----- 4) Bootstrap Estimation -----\n")
# For quick testing, we use B = 30; for more accurate results, increase B.
boot_results <- bootstrap_est_effects(X = sim_data$X,
                                      Y = sim_data$Y,
                                      M = sim_data$M,
                                      A = sim_data$A,
                                      B = 100,
                                      logY = FALSE,
                                      parallel = TRUE,
                                      nfold1 = 5,
                                      nfold2 = 3,
                                      net_alpha = 1,
                                      zeta = 0.5,
                                      optimizer = "pogs",
                                      intercept = FALSE,
                                      verbose = FALSE)

cat("\nBootstrap summary:\n")
print(boot_results$summary)
cat("\nFirst few rows of bootstrap samples:\n")
head(boot_results$boot_samples)

# 6. Plot checking for bootstrap samples
cat("\n----- 5) Plot Checking for Bootstrap Samples -----\n")
boot_plots <- plots_check(boot_results$boot_samples)
# In an interactive environment, view the plots:
print(boot_plots$histogram)
print(boot_plots$qqplot)

