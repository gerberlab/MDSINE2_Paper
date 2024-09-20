#!/bin/bash


n_perts=$1
require_variable "n_perts" "${n_perts}"


for (( i = 0; i < ${N_TRAJ_SEEDS}; ++i)); do
  python semisynthetic2/python_helpers/simulate_trajectories.py \
  -r "${REAL_DATA_PKL}" \
  -n "${MAX_N_MICE}" \
  -p "${n_perts} " \
  -i "${DATASET_DIR}/trajectory_replicate_${i}/initial_condition.npy" \
  -o "${DATASET_DIR}/trajectory_replicate_${i}/perts_${n_perts}" \
  -g "${DATASET_DIR}/truth" \
  -dt "${SIM_DT}" \
  -sm "${SIM_MAX}"
done
