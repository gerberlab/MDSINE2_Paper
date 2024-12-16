#!/bin/bash
set -e
source semisynthetic2/settings.sh


REAL_POSTERIOR_DIR=/data/local/MDSINE2_files/merged_studies_2

read_depth=75000
for (( trial = 0; trial < ${NUM_SAMPLE_TRIALS}; trial++ )); do
  sample_seed=${trial}

  echo "[*] Sampling dataset for trial=${trial} (read depth=${read_depth})"
  echo "Using negbin noise model."
  raw_dataset=${DATASET_DIR}/trial_${trial}/raw/dataset.pkl
  raw_replicate=${DATASET_DIR}/trial_${trial}/raw/replicate.pkl

  python semisynthetic/generate_dataset_fwsim.py \
    --study ${STUDY_PKL} \
    --growths ${REAL_POSTERIOR_DIR}/growth.npy \
    --interactions ${REAL_POSTERIOR_DIR}/interactions.npy \
    --perts ${REAL_POSTERIOR_DIR}/perturbations.npz \
    --coclust ${REAL_POSTERIOR_DIR}/coclusters.npy \
    --truth-dir ${DATASET_DIR}/truth \
    --out ${raw_dataset} \
    --replicate-out ${raw_replicate} \
    --read-depth ${read_depth} \
    --a0 ${NEGBIN_A0} \
    --a1 ${NEGBIN_A1} \
    --qpcr-noise-scale ${QPCR_NOISE_SCALE} \
    --seed ${sample_seed} \
    --sim-max 1e20 \
    --sim-dt 0.01
#    --sim-max 1e20
done
