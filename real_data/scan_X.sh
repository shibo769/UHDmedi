#!/bin/bash
cd "$(dirname "$0")"   # run from this script's folder
for xf in X_linear X_order2 X_order3; do
 for sd in 333 444 555 666; do
  MC_BOOT=2 MC_CORES=2 MC_SEED=$sd MC_SDT=0.2 MC_NSPLIT=4 MC_XFILE=${xf}.rds MC_LAMBDA=lambda.min \
    Rscript RealData_stable.R > real_data_outputs/scanX_${xf}_s${sd}.log 2>&1 &
 done
done
wait
