#!/bin/bash
cd "$(dirname "$0")"   # run from this script's folder
for sdt in 0.20 0.22 0.25 0.28; do
 for sd in 333 444 555; do
  MC_BOOT=2 MC_CORES=2 MC_SEED=$sd MC_SDT=$sdt MC_NSPLIT=3 MC_XFILE=X_reconstructed_order2.rds \
    Rscript RealData_stable.R > real_data_outputs/scan_sd${sdt}_s${sd}.log 2>&1 &
 done
done
wait
