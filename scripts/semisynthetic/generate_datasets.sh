#!/bin/bash
set -e
source semisynthetic/settings.sh


for read_depth in "${READ_DEPTHS[@]}"; do
  for (( trial = 0; trial < ${NUM_SAMPLE_TRIALS}; trial++ )); do
    seed=${read_depth}${trial}
    python semisynthetic/generate_dataset.py \
      --study ${STUDY_PKL} \
      --fixed-module-pkl ${FIXED_MODULE_MCMC} \
      --read-depth ${read_depth} \
      --alpha ${ALPHA_DIRICHLET} \
      --qpcr-noise-scale ${QPCR_NOISE_SCALE} \
      --seed ${seed}
  done
done
