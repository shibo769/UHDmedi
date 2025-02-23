# UHDmedi

**UHDmedi** is an R package for estimating total, indirect, and direct effects under ultra-high dimensional mediators (UHD), including potential interaction effects. This package implements methods described in the preprint:

Bo, S., Ghassami, A., & Mukherjee, D. (2024). "A Debiased Estimator for the Mediation Functional in Ultra-High-Dimensional Setting in the Presence of Interaction Effects." <arXiv:2412.08827>.

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
install_github("XXX/UHDmedi")
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

## Estimation Procedure (LaTeX)

Below is the high-level debiasing algorithm for mediation effect estimation as described in the preprint:

```latex
\begin{algorithm}[H]
\caption{Estimation procedure} 
\begin{algorithmic}\small
\State \textbf{Step 1:} 
Split the data into two parts, say $\cD_1$ and $\cD_2$. Use $\cD_1$ to run $q$ parallel regressions of $(\mathbf M_c)_{.,j}, j = 1, \ldots, q$ on $\mathbf X_c$, to obtain $\hat \bB_{0_{1,.}}, \ldots, \hat \bB_{0_{q,.}}$ where each $\hat \bB_{0_{i,.}} \in \reals^p$. Concatenate them to obtain the following estimator for $\mathbf B_0$: 
\[
\hat{\mathbf B}_0 = [\hat \bB_{0_{1,.}} \quad \hat \bB_{0_{2,.}} \quad \cdots \quad \hat \bB_{0_{q,.}}]^\top.
\]
Use $\cD_2$ to estimate ${\phi}_1 = [{\beta}_1^\top \quad {\gamma}_1^\top]^\top$ by regressing $Y_t$ on $\mathbf W_t=[\mathbf X_t \quad \mathbf M_t]$ with an $\ell_1$ penalty:
\[
\hat{\phi}_1 = [{\hat \beta}_1^\top \quad  {\hat \gamma}_1^\top]^\top = \arg\min_{\tilde \phi_1} \left\{ \sum_{\{i: A_i = 1\}} \left( Y_i - \mathbf W_{i,.} \cdot \tilde\phi_1 \right)^2 + \lambda_1 \|\tilde \phi_1\|_1 \right\}.
\]
\State \textbf{Step 2:} Compute bias correction weights $\tau_1$ with contrast 
\[
a = \left[ \bar X^\top \quad (\hat{\mathbf B}_0 \bar{X})^\top \right]^\top,
\]
and tuning parameter $K_1$ using $\cD_2$:
\begin{equation}
        \tau_1 =  \operatorname{argmin}_{\tilde{\tau}_1}\left\{ \|\tilde{\tau}_1\|_2^2  \quad \text{subject to} \quad \|a - \mathbf{W}_t^\top \tilde{\tau}_1\|_{\infty} \leq K_1 \sqrt{\frac{\log(p+q)}{n_t}}, \|\tilde \tau_{1}\|_\infty \leq n_t^{-2/3} \right\}.
\end{equation}
Set 
\[
\hat{\theta}_{0,1} = a^\top\hat{\phi}_1 + \sum_{\{i: A_i = 1\}} \tau_{1,i} \left( Y_i - \mathbf W_{i,.} \cdot \hat{\phi}_1 \right).
\]
\State \textbf{Step 3:} Estimate 
\[
b = \bB_0^\top \hat \gamma_1
\]
by regressing $M_i^\top \hat{\gamma}_1$ on $X_i$ (with $A_i = 0$) with an $\ell_1$ penalty using $\cD_2$:
\[
\hat{b} = \arg\min_{\tilde b} \left\{ \sum_{\{i: A_i = 0\}} \left( M_i^\top \hat{\gamma}_1 - X_i^\top\tilde b \right)^2 + \lambda_2  \|\tilde b \|_1 \right\}.
\]
\State \textbf{Step 4:} Compute bias correction weight $\tau_2$ with contrast $\bar{X}$ and tuning parameter $K_2$ using $\cD_2$:
\begin{equation}
        \tau_2 = \operatorname{argmin}_{\tilde{\tau}_2} \left\{  \|\tilde{\tau}_2\|_2^2 \quad \text{ subject to } \quad \|\bar X - \mathbf X_c^\top \tilde{\tau}_2\|_\infty \leq K_2\sqrt{\frac{\log(p)}{n_c}}, \|\tilde \tau_{2}\|_\infty \leq n_c^{-2/3}  \right\}.
\end{equation}
Set 
\[
\hat{\theta}_{0,2} = \bar{X}^\top\hat{b} + \sum_{\{i: A_i = 0\}} \tau_{2,i} \left( M_i^\top \hat{\gamma}_1 - X_i^\top \hat{b} \right).
\]
\State \textbf{Step 5:} Return the final estimator 
\[
\hat{\theta}_0 = \hat{\theta}_{0,1} + \hat{\theta}_{0,2} - (\hat{\mathbf B}_0 \bar{X})^\top \hat{\gamma}_1.
\]
\end{algorithmic}
\label{algomain}
\end{algorithm}
```

