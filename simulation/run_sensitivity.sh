#!/bin/bash
# Sensitivity analysis of (K1, K2) for the modified algorithm.
# Tables 8-10: diagonal K1 = K2 = K in {1, ..., 3}, pq in {50,400,750,1250}, n in {750,1000,1250}.
# Table 11:    full 9x9 (K1, K2) grid at pq = 750, n = 1250.
# Common: base DGP with s = s_B = 3, sigma^2 = 0.5, rep = 400, no naive.
set -u
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
RS=$HOME/micromamba/envs/mc/bin/Rscript
OUT=simulation_outputs_modified/base
LOG=/tmp/sensitivity_run.log

# Wait for any other Simu_modified run to release the cores first.
while pgrep -f "Rscript Simu_modified.R" > /dev/null; do sleep 60; done

echo "===== SENSITIVITY START $(date +%T) =====" >> "$LOG"

KS="1 1.25 1.5 1.75 2 2.25 2.5 2.75 3"

# --- Diagonal runs (tables 8-10) ---
for K in $KS; do
  tag="diag_K${K/./p}"
  [ -d "$OUT/sens_$tag" ] && { echo "skip $tag (exists)" >> "$LOG"; continue; }
  DGP_NAME=base MC_SIGMA=0.5 MC_ITER=400 MC_S=3 MC_SB=3 SKIP_NAIVE=TRUE RUN_PLOTS=FALSE \
    MC_K1=$K MC_K2=$K MC_PQ=50,400,750,1250 MC_N=750,1000,1250 \
    "$RS" Simu_modified.R >> "$LOG" 2>&1
  mv "$OUT/sigma_0p5" "$OUT/sens_$tag"
  echo "done $tag $(date +%T)" >> "$LOG"
done

# --- Full grid at pq = 750, n = 1250 (table 11) ---
for K1 in $KS; do
  for K2 in $KS; do
    tag="grid_K1_${K1/./p}_K2_${K2/./p}"
    [ -d "$OUT/sens_$tag" ] && { echo "skip $tag (exists)" >> "$LOG"; continue; }
    DGP_NAME=base MC_SIGMA=0.5 MC_ITER=400 MC_S=3 MC_SB=3 SKIP_NAIVE=TRUE RUN_PLOTS=FALSE \
      MC_K1=$K1 MC_K2=$K2 MC_PQ=750 MC_N=1250 \
      "$RS" Simu_modified.R >> "$LOG" 2>&1
    mv "$OUT/sigma_0p5" "$OUT/sens_$tag"
    echo "done $tag $(date +%T)" >> "$LOG"
  done
done

echo "===== SENSITIVITY DONE $(date +%T) =====" >> "$LOG"
