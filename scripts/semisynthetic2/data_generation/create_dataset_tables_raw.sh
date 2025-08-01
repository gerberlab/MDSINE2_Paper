#!/bin/bash
set -e
source semisynthetic2/settings.sh


base_seed=271828

for n_perts in 0 1 2 3; do
  for (( traj_repl = 0; traj_repl < ${N_TRAJ_SEEDS}; ++traj_repl)); do
    for (( data_repl = 0; data_repl < ${N_DATA_SEEDS}; ++data_repl )); do
      data_gen_seed="${base_seed}${data_repl}"  # increment for next iter

      echo "[* create_dataset_tables_raw.sh] Generating mdsine2-readable tables for (n_perts=${n_perts} | traj_repl=${traj_repl} | data_repl=${data_repl})"
      python semisynthetic2/data_generation/python_helpers/dataset_tables_from_fwsim.py \
      --input-study "${REAL_DATA_PKL}" \
      --sim-dir "${DATASET_DIR}/trajectory_replicate_${traj_repl}/perts_${n_perts}/trajectories" \
      --outdir "${DATASET_DIR}/trajectory_replicate_${traj_repl}/perts_${n_perts}/data_replicate_${data_repl}/tables_raw" \
      --read-depth "${READ_DEPTH}" \
      --calibration-replicates 6 \
      --a0 "${NEGBIN_A0}" \
      --a1 "${NEGBIN_A1}" \
      --qpcr-noise-scale "${QPCR_NOISE_SCALE}" \
      --seed "${data_gen_seed}"
    done
  done
done
