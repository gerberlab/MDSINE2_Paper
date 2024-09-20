#!/bin/bash


base_seed=314159


for (( i = 0; i < ${N_TRAJ_SEEDS}; ++i)); do
  traj_seed=$((base_seed + i))  # increment for next iter
  out_dir=${DATASET_DIR}/trajectory_replicate_${i}
  mkdir -p "${out_dir}"

  echo "[* create_initial_condition] Handling replicate ${i} (seed=${traj_seed})"
  outpath="${out_dir}/initial_condition.npy"
  python semisynthetic2/data_generation/python_helpers/create_initial_condition.py \
  -r "${REAL_DATA_PKL}" \
  -n "${MAX_N_MICE}" \
  -s "${traj_seed}" \
  -o "${outpath}"
  echo "[* create_initial_condition] Saved initial conditions to ${outpath}."
done
