#!/bin/bash
set -e
source semisynthetic/settings.sh


for read_depth in "${READ_DEPTHS[@]}"; do
  for (( trial = 0; trial < ${NUM_SAMPLE_TRIALS}; trial++ )); do
    seed=${read_depth}${trial}

    echo "[*] Sampling dataset for read_depth=${read_depth} | trial=${trial}"
#    python semisynthetic/generate_dataset_fwsim.py \
#      --study ${STUDY_PKL} \
#      --fixed-module-pkl ${FIXED_MODULE_MCMC} \
#      --cache-dir ${DATASET_DIR}/_tmp \
#      --out ${DATASET_DIR}/read_depth_${read_depth}/trial_${trial}/dataset.pkl \
#      --read-depth ${read_depth} \
#      --alpha ${ALPHA_DIRICHLET} \
#      --qpcr-noise-scale ${QPCR_NOISE_SCALE} \
#      --seed ${seed} \
#      --sim-dt 0.01 \
#      --sim-max 1e20
    echo "using generate_dataset2 (median traj)"
    python semisynthetic/generate_dataset2.py \
      --study ${STUDY_PKL} \
      --fixed-module-pkl ${FIXED_MODULE_MCMC} \
      --out ${DATASET_DIR}/read_depth_${read_depth}/trial_${trial}/dataset.pkl \
      --replicate-out ${DATASET_DIR}/read_depth_${read_depth}/trial_${trial}/replicate.pkl \
      --read-depth ${read_depth} \
      --alpha ${ALPHA_DIRICHLET} \
      --qpcr-noise-scale ${QPCR_NOISE_SCALE} \
      --seed ${seed}
  done
done
