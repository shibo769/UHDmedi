#!/bin/bash
cd "$(dirname "$0")"   # run from this script's folder
for sd in 333 444 555; do
  MC_BOOT=2 MC_CORES=2 MC_SEED=$sd MC_SDT=0.2 MC_NSPLIT=3 MC_XFILE=X_reconstructed_order2.rds MC_LAMBDA=lambda.min \
    Rscript RealData_stable.R > real_data_outputs/lammin_s${sd}.log 2>&1 &
done
wait
