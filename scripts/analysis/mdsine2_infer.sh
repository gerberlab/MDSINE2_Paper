#!/bin/bash
set -e
data_modality=$1
if [ "${data_modality}" == "healthy" ]; then
  source analysis/settings_healthy.sh
elif [ "${data_modality}" == "uc" ]; then
  source analysis/settings_uc.sh
else
  echo "data_modality argument is required and must be either 'healthy' or 'uc'. Exiting."
  exit 1
fi

# Healthy cohort
# --------------

seed=$2
require_variable "seed" $seed


study_name="${data_modality}-seed${seed}"

# Seed 0
echo "[*] Performing MDSINE2 inference on ${study_name}"
echo "[*] Writing files to ${MDSINE2_OUT_DIR}"

export MDSINE2_LOG_INI="${PROJECT_DIR}/scripts/analysis/logging_to_file.ini"
export LOG_FILEPATH="${MDSINE2_OUT_DIR}/mdsine2_inference_${study_name}.log"

mkdir -p "${MDSINE2_OUT_DIR}"
touch $LOG_FILEPATH
echo "LOGGING START/END TIMES to ${MDSINE2_OUT_DIR}"
start_time=$(date +%s.%N)

mdsine2 infer \
		--input $MAIN_DSET \
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

# record runtime for inference
end_time=$(date +%s.%N)
duration=$(echo "$end_time - $start_time" | bc)
echo "${duration}" >> ${MDSINE2_OUT_DIR}/runtime_inference.txt


echo "[*] Visualizing output of ${study_name}"
export LOG_FILEPATH="${MDSINE2_OUT_DIR}/mdsine2_visualization_${study_name}.log"
mdsine2 visualize-posterior \
		--chain $MDSINE2_OUT_DIR/$study_name/mcmc.pkl \
		--output-basepath $MDSINE2_OUT_DIR/$study_name/posterior

echo "Finished inference."
