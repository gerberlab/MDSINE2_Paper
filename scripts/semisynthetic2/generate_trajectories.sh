#!/bin/bash


n_perts=$1
require_variable "n_perts" "${n_perts}"
base_seed=314159


for (( i = 0; i < ${N_TRAJ_SEEDS}; ++i)); do
  traj_seed=$((base_seed + i))  # increment for next iter
  python semisynthetic2/python_helpers/simulate_trajectories.py \
  -r "${REAL_DATA_PKL}" \
  -n "${MAX_N_MICE}" \
  -p "${n_perts} " \
  -s "${traj_seed}" \
  -o "${DATASET_DIR}/perts_${n_perts}/trajectory_replicate_${i}" \
  -g "${DATASET_DIR}/truth" \
  -dt "${SIM_DT}"
done
