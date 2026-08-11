#!/bin/bash
cd "$(dirname "$0")"   # run from this script's folder
for sdt in 0.15 0.25; do
  MC_BOOT=200 MC_CORES=15 MC_SEED=333 MC_SDT=$sdt MC_NSPLIT=6 MC_XFILE=X_reconstructed_order2.rds \
    MC_BOOT_M=520 MC_REPLACE=0 MC_TRIM=0.05 \
    Rscript RealData_stable.R > real_data_outputs/boot_sd${sdt}.log 2>&1 &
done
wait
