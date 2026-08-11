#!/bin/bash
cd "$(dirname "$0")"   # run from this script's folder
for ns in 10 12; do
  MC_BOOT=300 MC_CORES=15 MC_SEED=333 MC_SDT=0.2 MC_NSPLIT=$ns MC_XFILE=X_reconstructed_order2.rds \
    MC_BOOT_M=520 MC_REPLACE=0 MC_TRIM=0.05 \
    Rscript RealData_stable.R > real_data_outputs/boot_ns${ns}.log 2>&1 &
done
wait
