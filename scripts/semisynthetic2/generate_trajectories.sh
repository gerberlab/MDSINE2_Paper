#!/bin/bash


python semisynthetic2/python_helpers/simulate_trajectories.py \
-r "${real_data_pkl}" \
-n "${MAX_N_MICE}" \
-p "${n_perts} " \
-s "${traj_seed}" \
-o "${DATASET_DIR}/trajectories/seed_${traj_seed}" \
-g "${ground_truth_dir}" \
-dt "${sim_dt}"

