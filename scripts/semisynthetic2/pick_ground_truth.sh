#!/bin/bash
set -e
source semisynthetic2/settings.sh


REAL_POSTERIOR_DIR=/data/local/MDSINE2_files/merged_studies_2
echo "[!] Picking ground truth."

python semisynthetic2/python_helpers/extract_glv_model.py \
-s "${STUDY_PKL}" \
--growths ${REAL_POSTERIOR_DIR}/growth.npy \
--interactions ${REAL_POSTERIOR_DIR}/interactions.npy \
--perts ${REAL_POSTERIOR_DIR}/perturbations.npz \
--coclust ${REAL_POSTERIOR_DIR}/coclusters.npy \
-o "${DATASET_DIR}/truth"
