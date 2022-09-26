#!/bin/bash
set -e
source analysis/settings.sh

# Healthy cohort
# --------------

seed=$1
require_variable "seed" $seed


study_name="healthy-seed${seed}"

# Seed 0
echo "[*] Performing MDSINE2 inference on ${study_name}"
echo "[*] Writing files to ${MDSINE2_OUT_DIR}"

export MDSINE2_LOG_INI="${PROJECT_DIR}/scripts/analysis/logging_to_file.ini"
export LOG_FILEPATH="${MDSINE2_OUT_DIR}/mdsine2_inference_${study_name}.log"

touch $LOG_FILEPATH
mdsine2 infer \
		--input $HEALTHY_DSET \
		--negbin $REPLICATE_MCMC \
		--seed $seed \
		--burnin $BURNIN \
		--n-samples $N_SAMPLES \
		--checkpoint $CHECKPOINT \
		--multiprocessing $MULTIPROCESSING \
		--rename-study $study_name \
		--basepath $MDSINE2_OUT_DIR \
		--interaction-ind-prior $INTERACTION_IND_PRIOR \
		--perturbation-ind-prior $PERTURBATION_IND_PRIOR

echo "[*] Visualizing output of ${study_name}"
export LOG_FILEPATH="${MDSINE2_OUT_DIR}/mdsine2_visualization_${study_name}.log"
mdsine2 visualize-posterior \
		--chain $MDSINE2_OUT_DIR/$study_name/mcmc.pkl \
		--output-basepath $MDSINE2_OUT_DIR/$study_name/posterior

echo "Finished inference."
