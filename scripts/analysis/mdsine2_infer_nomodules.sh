#!/bin/bash
set -e
source analysis/settings.sh

echo "Running MDSINE2 model"
echo "Writing files to ${MDSINE2_OUT_DIR}"

# Healthy cohort
# --------------

seed=0
study_name="healthy-seed${seed}-nomodule"

# Seed 0
echo "[*] Performing MDSINE2 inference on ${study_name}"
export MDSINE2_LOG_INI="${PROJECT_DIR}/scripts/analysis/logging_to_file.ini"
export LOG_FILEPATH="${MDSINE2_OUT_DIR}/mdsine2_inference_${study_name}_NOMODULE.log"
mkdir -p ${MDSINE2_OUT_DIR}
mdsine2 infer \
		--input $HEALTHY_DSET \
		--negbin $REPLICATE_MCMC \
		--seed $seed \
		--burnin $BURNIN \
		--n-samples $N_SAMPLES \
		--checkpoint $CHECKPOINT \
		--rename-study $study_name \
		--basepath $MDSINE2_OUT_DIR \
		--interaction-ind-prior $INTERACTION_IND_PRIOR \
		--perturbation-ind-prior $PERTURBATION_IND_PRIOR \
		--nomodules

echo "[*] Visualizing output of ${study_name}"
export LOG_FILEPATH="${MDSINE2_OUT_DIR}/mdsine2_visualization_${study_name}_NOMODULE.log"
mdsine2 visualize-posterior \
		--chain $MDSINE2_OUT_DIR/$study_name/mcmc.pkl \
		--output-basepath $MDSINE2_OUT_DIR/$study_name/posterior

echo "Finished inference."
