# UHDmedi

**UHDmedi** provides R research code for estimating natural direct, natural indirect, and total effects when both covariates and mediators are ultra-high-dimensional. The method uses multi-stage debiasing and accommodates treatment–covariate and treatment–mediator interactions.

## Paper

Shi Bo, AmirEmad Ghassami, and Debarghya Mukherjee. **Debiased Causal Mediation Analysis in Ultra-High-Dimensional Settings in the Presence of Interaction Effects.** [arXiv:2412.08827](https://arxiv.org/abs/2412.08827).

## Repository structure

| Directory | Contents |
| --- | --- |
| [`simulation/`](simulation/) | Main simulation (`Simu_modified.R`) and sensitivity-analysis launcher |
| [`comparison_with_Guo/`](comparison_with_Guo/) | Guo et al. comparison, the no-interaction version of our method, and result summaries |
| [`real_data/`](real_data/) | Lung-cancer DNA methylation analyses and related analysis scripts |
| [`shared/`](shared/) | Common utilities, balancing solvers, and dependency installation |

This repository is organized as research scripts. Processed real-data inputs and experiment outputs are not included.

## Setup

Install R with `Rscript` available, then clone the repository and install the dependencies:

```bash
git clone https://github.com/shibo769/UHDmedi.git
cd UHDmedi
Rscript shared/setup_packages.R
```

The dependency list is maintained in [`shared/setup_packages.R`](shared/setup_packages.R). The debiased simulation requires the **pogs** solver. The installer attempts to obtain it from GitHub when it is not installed; the vendored source archive referenced in that script is not currently included in this repository.

## Run a simulation

From the repository root, the following Bash command selects a small configuration for checking the setup:

```bash
MC_ITER=2 MC_PQ=50 MC_N=250 MC_CORES=4 \
  Rscript simulation/Simu_modified.R
```

For the script's default simulation configuration:

```bash
Rscript simulation/Simu_modified.R
```

Simulation settings are controlled by environment variables in the script, including `MC_ITER` (replications), `MC_PQ` (dimension settings), `MC_N` (sample sizes), and `MC_CORES` (workers). Full experiment grids require substantially more computation than the small configuration above.

By default, results are written under `simulation/simulation_outputs_modified/`, organized by data-generating process and noise level. The script saves raw replication results and summaries as RDS files.

## Method comparison

The comparison scripts generate results for the no-interaction setting. Run both methods before aggregating their outputs:

```bash
Rscript comparison_with_Guo/Guo2022.R
Rscript comparison_with_Guo/NoInteraction.R
Rscript comparison_with_Guo/Results.R
```

Use matching data-generating-process and simulation settings across the two methods and the summary script. See the settings at the beginning of each file.

## Real-data analysis

The scripts in [`real_data/`](real_data/) require preprocessed `.rds` files containing treatment, covariate, outcome, and methylation data. These files are not distributed here. The paper describes the TCGA application and preprocessing.

Real-data scripts retain script-specific input and output paths that need to be configured before use. In `RealData_modified.R`, the `SHARED` path must also be initialized before `repo <- dirname(SHARED)` for standalone execution. The shell launchers contain environment-specific settings, including an R interpreter path in `simulation/run_sensitivity.sh`.

## License

The code is distributed under the [GNU General Public License v3.0](LICENSE).
