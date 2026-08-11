# High-dimensional Mediation — code

One method, three experiments. Code only (no data, no outputs).

```
mediation_code/
├── simulation/           # Simu_modified.R (main simulation) + run_sensitivity.sh
├── real_data/            # RealData_*.R (real-data analysis), check_M_raw.R, launch*/scan* .sh
├── comparison_with_Guo/  # Guo2022.R (Guo et al. 2022), NoInteraction.R (our method), Results.R
└── shared/               # dependencies sourced by all of the above:
                          #   mc_common.R, utilis.R (balancing QP),
                          #   utils_mediation.R, setup_packages.R
```

## Run
1. Install packages: `Rscript shared/setup_packages.R`
2. Each script auto-detects `shared/` relative to its own location — run from anywhere:
   `Rscript simulation/Simu_modified.R`, `Rscript comparison_with_Guo/Guo2022.R`, etc.
3. Real-data scripts read `.rds` inputs from `real_data/` (or set env `DATA_DIR`).
   The data files themselves are NOT included.

No absolute paths; all cross-file `source()` resolves via `shared/`.
