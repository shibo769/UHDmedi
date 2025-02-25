

# UHDmedi

**UHDmedi** is an R package for estimating total, indirect, and direct effects under ultra-high dimensional setting (UHD) via debiasing, including potential treatment-covariates, treatment-mediators interaction effects. This package implements methods described in:

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

## Packages and Optimizers

This package works with two optimizers:

- **pogs**: A first-order optimizer based on ADMM. We generally achieve the best performance using pogs. *Note:* pogs needs to be installed separately. You can install one of the pre-compiled binaries available from the [pogs project repository](https://github.com/foges/pogs). See [this page](https://github.com/foges/pogs/blob/master/src/interface_r/README.md) for further instructions.
- **quadprog**: A standard R optimization library.

We recommend trying **pogs** first for optimal performance.

To run the debiasing in the estimation procedure, we are very grateful to [balancHD](https://github.com/swager) for providing the `balancHD` package. For details and installation instructions, please visit [https://github.com/swager/balanceHD](https://github.com/swager/balanceHD).

## Overview
We consider a causal mediation model where the treatment $ A \in \{0,1\} $ affects the outcome $ Y $ both directly and indirectly through a set of mediators $ \mathbf{M} \in \mathbb{R}^{q} $. The system is defined as:


$$Y =  (1-A)\left(\alpha_{0}+X^{\top} \beta_{0}+M^{\top} \gamma_{0}\right)+A\left(\alpha_{1}+X^{\top} \beta_{1}+M^{\top} \gamma_{1}\right)+\epsilon$$

$$M = (1-A)\left(\delta_{0}+B_{0} X\right)+A\left(\delta_{1}+B_{1} X\right) + U$$


Here, $ X \in \mathbb{R}^{p} $ represents high-dimensional pre-treatment covariates, while $ \mathbf{M} \in \mathbb{R}^{q} $ are mediators that transmit part of the treatment effect. The parameters $ \beta_0, \beta_1 \in \mathbb{R}^{p} $ capture the direct effect of covariates on $ Y $, while $ \gamma_0, \gamma_1 \in \mathbb{R}^{q} $ quantify the influence of mediators. The matrices $ \mathbf{B}_0, \mathbf{B}_1 \in \mathbb{R}^{q \times p} $ determine how covariates influence the mediators.

This model accommodates **treatment-mediator interactions** ($\gamma_1 \neq \gamma_0$), meaning the effect of mediators on $ Y $ depends on treatment status, and **treatment-covariate interactions** ($\beta_1 \neq \beta_0$), allowing covariates to have different influences on the outcome under treatment and control. Such flexibility is essential for capturing heterogeneous treatment effects in high-dimensional settings.

For detailed estimation procedure of mediation functional and effects, please see Section 3 and Appendix in [Arxiv](https://arxiv.org/pdf/2412.08827).


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
## Example Output

### **Simulated Data Output**
```r
> sim_data <- simData(nn = 30, p = 50, q = 20, variance = 0.1, k1 = 3, k2 = 6, s_X = 3, categorical = TRUE)
The dim of X 30 60 
The dim of M 30 20 
The dim of A 30 1 
The dim of Y 30 1 
Dim of (X,M): 80
> str(sim_data)
List of 4
 $ X: num [1:30, 1:60] 0.6201 2.1231 -0.4608 -0.2354 0.0602 ...
 $ A: int [1:30] 0 0 1 1 0 0 1 1 0 0 ...
 $ Y: num [1:30] 0.979 2.729 2.325 5.099 12.471 ...
 $ M: num [1:30, 1:20] 0.574 0.408 0.686 0.439 0.707 ...

eda_report()
==========================================================
                  EDA REPORT SUMMARY
==========================================================

----------------------------------------------------------
                     DIMENSIONS
----------------------------------------------------------
X dimension: [1] 30 60
M dimension: [1] 30 20
Length of Y: 30 
Length of A: 30 

----------------------------------------------------------
                    MISSING VALUES
----------------------------------------------------------
Missing values in X:
  No missing values in X.

Missing values in M:
  No missing values in M.

Missing values in Y: 0 
Missing values in A: 0 

----------------------------------------------------------
                      SUMMARY OF Y
----------------------------------------------------------
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
 -4.261   2.510   5.058   4.597   6.030  12.471 

----------------------------------------------------------
                    DISTRIBUTION OF A
----------------------------------------------------------
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
 0.0000  0.0000  1.0000  0.5667  1.0000  1.0000 

----------------------------------------------------------
                 CHARACTERISTICS FOR X
----------------------------------------------------------
Number of numeric columns in X: 50 
Number of categorical columns in X: 10 
Range of X: [ -3.298562 ,  3.174249 ]

----------------------------------------------------------
                 CHARACTERISTICS FOR M
----------------------------------------------------------
Range of M: [ 0.05860127 ,  0.987907 ]

----------------------------------------------------------
     TOP 10 CORRELATIONS (BY ABS. VALUE) WITH Y
----------------------------------------------------------
M with Y:
 [1]  0.3240070 -0.3183457  0.2830364 -0.2557979  0.2479194  0.2473952 -0.2398461  0.2188564  0.2102471  0.2078837

X with Y:
 [1]  0.5563897  0.4977343 -0.4336580  0.4038666 -0.3934824  0.3362304  0.3352937  0.3197300 -0.3153891 -0.2988701

----------------------------------------------------------
               VARIANCE INFORMATION FOR M
----------------------------------------------------------
Minimum variance in M:
  Column: 17 with variance: 0.01505739 
Maximum variance in M:
  Column: 2 with variance: 0.07436384 

Scatterplot for the M column with minimum variance vs Y is stored in plots$scatter_plot.

==========================================================
                  END OF EDA REPORT
==========================================================


effects_result()
----------------------------------------------------------
               Debiased Y(.) Estimates
----------------------------------------------------------
Table: Debiased Y(.) Estimates

|    | Estimate |
|:---|:--------:|
|Y11 | 5.440090 |
|Y10 | 5.222723 |
|Y00 | 5.102479 |

----------------------------------------------------------
               Effects Estimates
----------------------------------------------------------
Table: Effects Estimates

|    | Estimate  |
|:---|:---------:|
|NIE | 0.2173661 |
|NDE | 0.1202444 |
|TE  | 0.3376104 |

bootstrap_est_effects()
----------------------------------------------------------
               Bootstrap Summary Statistics
----------------------------------------------------------

Table: Bootstrap Summary Statistics

|    | Effect |    Mean    |    SD     | Lower_95CI | Upper_95CI |
|:---|:------:|:----------:|:---------:|:----------:|:----------:|
|TE  |   TE   | 0.7890942  | 1.3970916 | -1.840702  |  3.387049  |
|NIE |  NIE   | -0.1514042 | 0.7229927 | -1.465652  |  1.082800  |
|NDE |  NDE   | 0.9404984  | 1.6913323 | -2.364615  |  3.984233  |
|Y11 |  Y11   | 5.5248530  | 0.6630805 |  4.329847  |  6.902516  |
|Y10 |  Y10   | 5.6762572  | 1.0505229 |  3.862940  |  7.820340  |
|Y00 |  Y00   | 4.7357587  | 1.2804549 |  2.244986  |  6.964613  |

plots_check(boot_results$boot_samples)
```
<img src="https://github.com/user-attachments/assets/72452237-7876-421a-a607-eac2e530b6ae" width="500">
<img src="https://github.com/user-attachments/assets/d49901bf-f1cc-46d5-96b1-f8a4b4e107b3" width="500">

# References
Bo, S., Ghassami, A., & Mukherjee, D. (2024). "A Debiased Estimator for the Mediation Functional in Ultra-High-Dimensional Setting in the Presence of Interaction Effects." [arXiv:2412.08827](https://arxiv.org/abs/2412.08827)

Athey, S., Imbens, G. W., & Wager, S. (2018). Approximate residual balancing: debiased inference of average treatment effects in high dimensions. Journal of the Royal Statistical Society Series B: Statistical Methodology, 80(4), 597-623.
