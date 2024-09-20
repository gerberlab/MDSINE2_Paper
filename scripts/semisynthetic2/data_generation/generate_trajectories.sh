#!/bin/bash


for n_perts in 0 1 2 3; do
  for (( replicate = 0; replicate < ${N_TRAJ_SEEDS}; ++replicate)); do
    echo "[* generate_trajectories.sh] Generating forward simulation trajectories (n_perts=${n_perts} | replicate=${replicate})"
    python semisynthetic2/data_generation/python_helpers/simulate_trajectories.py \
    -r "${REAL_DATA_PKL}" \
    -n "${MAX_N_MICE}" \
    -p "${n_perts} " \
    -i "${DATASET_DIR}/trajectory_replicate_${replicate}/initial_condition.npy" \
    -o "${DATASET_DIR}/trajectory_replicate_${replicate}/perts_${n_perts}" \
    -g "${DATASET_DIR}/truth" \
    -dt "${SIM_DT}" \
    -sm "${SIM_MAX}"
  done
done
