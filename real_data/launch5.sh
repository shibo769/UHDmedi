#!/bin/bash
cd "$(dirname "$0")"   # run from this script's folder
for sd in 111 222 333 444 555; do
  MC_BOOT=2 MC_CORES=2 MC_SEED=$sd MC_SDT=0.2 MC_NSPLIT=1 Rscript RealData_stable.R > real_data_outputs/pt2_seed$sd.log 2>&1 &
done
wait
