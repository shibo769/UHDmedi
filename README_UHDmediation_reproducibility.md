# Reproducibility Materials for the UHD Mediation Manuscript

This repository contains the R code, processed data files, saved numerical outputs, and final tables/figures used to reproduce the numerical results reported in the revised manuscript on ultra-high-dimensional mediation analysis.

The repository is organized so that readers can either:

1. run the main simulation and real-data scripts from the raw/processed inputs, or  
2. inspect the saved outputs in `BootResults/` and `Tables and Figures/` that were used to prepare the manuscript tables and figures.

The scripts use fixed random seeds where applicable. Because the simulation studies and bootstrap analyses involve Monte Carlo replications, parallel computation, and numerical optimization, small numerical differences may occur across operating systems, BLAS/LAPACK implementations, and package versions. The saved outputs included in this repository correspond to the numerical results reported in the revised manuscript.

---

## Repository structure

```text
UHDmediation/
├── README.md
├── utilis.R
├── BootResults/
├── RealData/
├── Simulation/
└── Tables and Figures/
```

### Top-level files

- `README.md`  
  This file. It describes how to reproduce the numerical results in the manuscript.

- `utilis.R`  
  Common helper functions used by the simulation and real-data analysis scripts, including residual balancing and related numerical routines.

---

## Folder descriptions

### `Simulation/`

This folder contains the scripts for the simulation studies.

```text
Simulation/
├── Simulation_main.R
├── Sensitivity_K_gridsearch.R
└── Comparison/
    ├── Guo2022.R
    ├── NoInteraction.R
    ├── Results.R
    └── utils_mediation.R
```

- `Simulation_main.R`  
  Main simulation script for the proposed estimator. This script generates synthetic data under the simulation designs described in the manuscript and computes the simulation summaries used for the main simulation tables.

- `Sensitivity_K_gridsearch.R`  
  Sensitivity analysis over the tuning parameters `K1` and `K2` used in the residual balancing step. This script produces the grid-search results for choosing and checking the sensitivity of the tuning parameters.

- `Comparison/Guo2022.R`  
  Implementation of the comparison method based on Guo et al. (2022), as used in the simulation comparison study.

- `Comparison/NoInteraction.R`  
  Simulation script for the no-interaction setting used in the comparison study.

- `Comparison/Results.R`  
  Script used to summarize the comparison results and prepare the corresponding manuscript table entries.

- `Comparison/utils_mediation.R`  
  Helper functions used by the comparison scripts.

---

### `RealData/`

This folder contains the real-data analysis scripts and processed data files.

```text
RealData/
├── RealData_without_intercept_threshold_main.Rmd
├── RealDataLowD_without_intercept.Rmd
└── real_data/
    ├── A_Smoking.rds
    ├── gene_sites_for_M_raw_cleaned.rds
    ├── M_important_by_paper.rds
    ├── M_raw_cleaned.rds
    ├── X_fourthorder_MICE.rds
    ├── X_NoPolyInteraction.rds
    └── Y_SurvivalTime.rds
```

- `RealData_without_intercept_threshold_main.Rmd`  
  Main high-dimensional real-data analysis script. This file analyzes the lung cancer methylation data using the proposed high-dimensional mediation estimator.

- `RealDataLowD_without_intercept.Rmd`  
  Low-dimensional real-data analysis script using the selected mediator set.

- `real_data/A_Smoking.rds`  
  Treatment indicator used in the real-data analysis.

- `real_data/X_NoPolyInteraction.rds`  
  Covariate matrix used in the low-dimensional real-data analysis.

- `real_data/X_fourthorder_MICE.rds`  
  Covariate matrix used in the high-dimensional real-data analysis.

- `real_data/Y_SurvivalTime.rds`  
  Outcome variable used in the real-data analysis.

- `real_data/M_important_by_paper.rds`  
  Selected mediator matrix used in the low-dimensional real-data analysis.

- `real_data/M_raw_cleaned.rds`  
  Full cleaned mediator matrix used in the high-dimensional real-data analysis.

- `real_data/gene_sites_for_M_raw_cleaned.rds`  
  Annotation file for the mediator matrix.

The real-data scripts read these processed `.rds` files directly. If the repository is shared publicly, these files should only be included if permitted by the relevant data-use and privacy policies. If the original real data cannot be redistributed, see the section **Real-data availability note** below.

---

### `BootResults/`

This folder contains saved bootstrap outputs from the real-data analyses.

```text
BootResults/
├── bootresults_real_data_fixed_dim.xlsx
├── bootresults_real_data_high_dim.xlsx
├── NDE_histogram_95pct_quantile_CI.png
├── NIE_histogram_95pct_quantile_CI.png
└── TE_histogram_95pct_quantile_CI.png
```

- `bootresults_real_data_fixed_dim.xlsx`  
  Bootstrap results for the fixed-/low-dimensional real-data analysis.

- `bootresults_real_data_high_dim.xlsx`  
  Bootstrap results for the high-dimensional real-data analysis.

- `NIE_histogram_95pct_quantile_CI.png`  
  Histogram of bootstrap estimates for the natural indirect effect.

- `NDE_histogram_95pct_quantile_CI.png`  
  Histogram of bootstrap estimates for the natural direct effect.

- `TE_histogram_95pct_quantile_CI.png`  
  Histogram of bootstrap estimates for the total effect.

These files are included so that readers can reproduce the reported real-data point estimates and confidence intervals without rerunning the full bootstrap procedure.

---

### `Tables and Figures/`

This folder contains the final manuscript tables and figures.

```text
Tables and Figures/
├── Figure.1/
├── Table.1/
├── Table.2/
├── Table.3/
├── Table.4-7/
├── Table.8/
└── Table.9 and Figure 3/
```

- `Figure.1/`  
  Final figure files for Figure 1.

- `Table.1/`  
  Final LaTeX/table output for Table 1.

- `Table.2/`  
  Final LaTeX/table output for Table 2.

- `Table.3/`  
  Final LaTeX/table output for Table 3.

- `Table.4-7/`  
  Final LaTeX/table outputs for Tables 4--7.

- `Table.8/`  
  Final LaTeX/table output for Table 8.

- `Table.9 and Figure 3/`  
  Final table and figure outputs for Table 9 and Figure 3.

---

## Software requirements

The analyses were conducted in R. The main scripts require the following R packages:

```r
install.packages(c(
  "data.table",
  "MASS",
  "glmnet",
  "foreach",
  "doParallel",
  "doRNG",
  "boot",
  "tidyverse",
  "pracma",
  "caret",
  "kableExtra"
))
```

The following packages are also used by the debiased estimation and residual balancing scripts:

```r
# Install if not already available
install.packages("pogs")
```

The package `balanceHD` is required for residual balancing. Please install it following the package instructions before running the real-data and debiased simulation scripts.

---

## Recommended R version

The code was prepared for R 4.x. We recommend using a recent R 4.x release and running the scripts from the root directory of the repository:

```r
setwd("path/to/UHDmediation")
```

All paths in the scripts should be relative to the repository root.

---

## Reproducing the simulation results

From the repository root, run:

```r
source("Simulation/Simulation_main.R")
```

This script runs the main simulation study for the proposed estimator. It generates synthetic data, computes the proposed estimator and comparison summaries, and produces the numerical results used in the simulation tables.

The simulation can be computationally intensive because it uses Monte Carlo repetitions, high-dimensional regressions, residual balancing, and parallel computation.

---

## Reproducing the tuning-parameter sensitivity analysis

From the repository root, run:

```r
source("Simulation/Sensitivity_K_gridsearch.R")
```

This script performs a grid search over the residual balancing tuning parameters `K1` and `K2`. The output is used to assess the sensitivity of the estimator to the tuning parameter choice.

---

## Reproducing the comparison-method results

The comparison scripts are located in:

```text
Simulation/Comparison/
```

A typical running order is:

```r
source("Simulation/Comparison/NoInteraction.R")
source("Simulation/Comparison/Guo2022.R")
source("Simulation/Comparison/Results.R")
```

The file `Simulation/Comparison/utils_mediation.R` contains helper functions used by these scripts.

---

## Reproducing the real-data analysis

### High-dimensional real-data analysis

From the repository root, run:

```r
rmarkdown::render("RealData/RealData_without_intercept_threshold_main.Rmd")
```

This script reproduces the high-dimensional lung cancer methylation analysis.

### Low-dimensional real-data analysis

From the repository root, run:

```r
rmarkdown::render("RealData/RealDataLowD_without_intercept.Rmd")
```

This script reproduces the low-dimensional real-data analysis based on the selected mediator set.

### Bootstrap summaries

The full bootstrap procedures may take a long time. The saved bootstrap outputs are provided in:

```text
BootResults/
```

The reported point estimates and confidence intervals can be recomputed from these saved bootstrap files.

---

## Mapping between manuscript results and repository files

| Manuscript item | Main script/source | Saved output location |
|---|---|---|
| Figure 1 | `Simulation/Simulation_main.R` | `Tables and Figures/Figure.1/` |
| Table 1 | `Simulation/Simulation_main.R` | `Tables and Figures/Table.1/` |
| Table 2 | `Simulation/Simulation_main.R` | `Tables and Figures/Table.2/` |
| Table 3 | `Simulation/Comparison/` | `Tables and Figures/Table.3/` |
| Tables 4--7 | `Simulation/Simulation_main.R` and related simulation settings | `Tables and Figures/Table.4-7/` |
| Table 8 | `RealData/RealDataLowD_without_intercept.Rmd` | `Tables and Figures/Table.8/` |
| Table 9 | `RealData/RealData_without_intercept_threshold_main.Rmd` | `Tables and Figures/Table.9 and Figure 3/` |
| Figure 3 | `RealData/RealData_without_intercept_threshold_main.Rmd` | `Tables and Figures/Table.9 and Figure 3/` |
| Real-data bootstrap histograms | Real-data bootstrap scripts/Rmd files | `BootResults/` |

---

## Notes on random seeds and parallel computation

The simulation and bootstrap scripts use fixed random seeds where applicable. For parallel computation, the scripts use `foreach`, `doParallel`, and `doRNG` so that Monte Carlo repetitions are reproducible under the same computational setup.

Users may observe small numerical differences across machines due to:

- different R versions,
- different versions of numerical optimization packages,
- different BLAS/LAPACK libraries,
- different parallel backends,
- convergence tolerances in high-dimensional regressions or residual balancing.

The saved numerical outputs in `BootResults/` and `Tables and Figures/` correspond to the results used in the revised manuscript.

---

## Real-data availability note

If the processed real-data files in `RealData/real_data/` can be shared under the relevant data-use and privacy policies, the files in that folder allow the real-data analysis to be reproduced directly.

If the real-data files cannot be publicly redistributed, the analysis scripts are still provided to document the full workflow. In that case, users should replace the real-data files with data objects of the same names and dimensions, or use a simulated data-generating script with the same structure as the real-data analysis.

The required real-data objects are:

```text
A_Smoking.rds
X_NoPolyInteraction.rds
X_fourthorder_MICE.rds
Y_SurvivalTime.rds
M_important_by_paper.rds
M_raw_cleaned.rds
gene_sites_for_M_raw_cleaned.rds
```

---

## Full reproduction workflow

To reproduce the main numerical results from the repository root, use:

```r
# Main simulation
source("Simulation/Simulation_main.R")

# Sensitivity analysis over K1 and K2
source("Simulation/Sensitivity_K_gridsearch.R")

# Comparison methods
source("Simulation/Comparison/NoInteraction.R")
source("Simulation/Comparison/Guo2022.R")
source("Simulation/Comparison/Results.R")

# Real-data analyses
rmarkdown::render("RealData/RealDataLowD_without_intercept.Rmd")
rmarkdown::render("RealData/RealData_without_intercept_threshold_main.Rmd")
```

For a faster check of the manuscript results, inspect the saved outputs in:

```text
BootResults/
Tables and Figures/
```

---

## Suggested citation statement for the manuscript

The following sentence can be included in the main text or appendix of the revised manuscript:

```latex
The code and processed data required to reproduce the numerical results, simulation studies, real-data analyses, tables, and figures are available at \url{INSERT_GITHUB_REPOSITORY_LINK}. The repository contains the R scripts, processed data files, saved bootstrap outputs, final tables and figures, and documentation describing the correspondence between each script and the numerical results reported in the manuscript.
```

If the real data cannot be redistributed, use the following version instead:

```latex
The code required to reproduce the simulation studies, tables, and figures is available at \url{INSERT_GITHUB_REPOSITORY_LINK}. Due to data-use restrictions, the original real-data files cannot be redistributed. The repository provides the full analysis scripts and saved reproducibility outputs, and documents the required processed data structure for reproducing the real-data workflow.
```

---

## Contact

For questions about the reproducibility materials, please contact the corresponding author of the manuscript.
