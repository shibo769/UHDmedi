# UHDmedi

**UHDmedi** is an R package for estimating total, indirect, and direct effects under ultra-high dimensional mediators (UHD) via debiasing, including potential interaction effects. This package implements methods described in:

Bo, S., Ghassami, A., & Mukherjee, D. (2024). "A Debiased Estimator for the Mediation Functional in Ultra-High-Dimensional Setting in the Presence of Interaction Effects."[arXiv:2412.08827](https://arxiv.org/abs/2412.08827)

The package provides the following key functionalities:
- **simData()**: Generates synthetic data for high-dimensional mediation analysis.
- **real_data()**: Loads real lung cancer data (M, X, A, Y) from local files.
- **eda_report()**: Performs exploratory data analysis (EDA), including missing value checks, summary statistics, range and correlation analysis, and scatterplots.
- **est_effects()**: Estimates total effect (TE), natural indirect effect (NIE), and natural direct effect (NDE) via debiasing.
- **bootstrap_est_effects()**: Conducts non-parametric bootstrap to compute confidence intervals for effect estimates.
- **plots_check()**: Generates diagnostic plots (histograms and QQ plots) of bootstrap samples.

## Installation

You can install **UHDmedi** in two ways:

### 1. Install from GitHub

```r
# Using devtools:
library(devtools)
install_github("shibo769/UHDmedi")
```

### 2. Install from Local Tarball

If you have the source tarball (e.g., `UHDmedi_0.1.0.tar.gz`):

```r
install.packages("path/to/UHDmedi_0.1.0.tar.gz", repos = NULL, type = "source")
```

## Optimizers

This package works with two optimizers:

- **pogs**: A first-order optimizer based on ADMM. We generally achieve the best performance using pogs. *Note:* pogs needs to be installed separately. You can install one of the pre-compiled binaries available from the [pogs project repository](https://github.com/foges/pogs). See [this page](https://github.com/foges/pogs/blob/master/src/interface_r/README.md) for further instructions.
- **quadprog**: A standard R optimization library.

We recommend trying **pogs** first for optimal performance.

## Overview
We consider a causal mediation model where the treatment \( A \in \{0,1\} \) affects the outcome \( Y \) both directly and indirectly through a set of mediators \( \mathbf{M} \in \mathbb{R}^{q} \). The system is defined as:

\begin{equation}
Y =  (1-A)\left(\alpha_{0}+X^{\top} \beta_{0}+M^{\top} \gamma_{0}\right)+A\left(\alpha_{1}+X^{\top} \beta_{1}+M^{\top} \gamma_{1}\right)+\epsilon
\end{equation}

\[
M = (1-A)\left(\delta_{0}+\mathbf{B}_{0} X\right)+A\left(\delta_{1}+\mathbf{B}_{1} X\right)+  U.
\]

Here, \( X \in \mathbb{R}^{p} \) represents high-dimensional pre-treatment covariates, while \( \mathbf{M} \in \mathbb{R}^{q} \) are mediators that transmit part of the treatment effect. The parameters \( \beta_0, \beta_1 \in \mathbb{R}^{p} \) capture the direct effect of covariates on \( Y \), while \( \gamma_0, \gamma_1 \in \mathbb{R}^{q} \) quantify the influence of mediators. The matrices \( \mathbf{B}_0, \mathbf{B}_1 \in \mathbb{R}^{q \times p} \) determine how covariates influence the mediators.

This model accommodates **treatment-mediator interactions** (\(\gamma_1 \neq \gamma_0\)), meaning the effect of mediators on \( Y \) depends on treatment status, and **treatment-covariate interactions** (\(\beta_1 \neq \beta_0\)), allowing covariates to have different influences on the outcome under treatment and control. Such flexibility is essential for capturing heterogeneous treatment effects in high-dimensional settings.

For detailed esimation procedure of mediation functional and effects, please see Section 3 and Appendix in (<https://arxiv.org/pdf/2412.08827>).

## Example Usage

Below is an example script demonstrating the main functions of the package:

```r
#################################
##   UHDmedi Package Example   ##
#################################

# 1) Generate simulated data
cat("----- 1) Generating Simulated Data -----\n")
sim_data <- simData(nn = 30, p = 50, q = 20, variance = 0.1, k1 = 3, k2 = 6, s_X = 3, categorical = TRUE)
str(sim_data)

# 2) Run basic EDA on the data
cat("\n----- 2) Running EDA -----\n")
eda_results <- eda_report(X = sim_data$X,
                          M = sim_data$M,
                          A = sim_data$A,
                          Y = sim_data$Y,
                          plot_corr = TRUE,
                          plot_Y = TRUE,
                          verbose = TRUE)

# View EDA summary
print(eda_results$summary)
# To view additional plots in an interactive environment, run:
print(eda_results$plots$scatter_plot_minvar)
print(eda_results$plots$scatter_plot_maxvar)
print(eda_results$plots$corr_plot)
print(eda_results$plots$Y_density)
print(eda_results$plots$Y_boxplot)

# 3) Estimate effects (ATE, NIE, NDE)
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

# 4) Perform non-parametric bootstrap estimation
cat("\n----- 4) Bootstrap Estimation -----\n")
# For quick testing, we use B = 100; increase B for more accurate results.
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

# 5) Plot checking for bootstrap samples
cat("\n----- 5) Plot Checking for Bootstrap Samples -----\n")
boot_plots <- plots_check(boot_results$boot_samples)
# In an interactive environment, view the plots:
print(boot_plots$histogram)
print(boot_plots$qqplot)
```


# References
Bo, S., Ghassami, A., & Mukherjee, D. (2024). "A Debiased Estimator for the Mediation Functional in Ultra-High-Dimensional Setting in the Presence of Interaction Effects."

