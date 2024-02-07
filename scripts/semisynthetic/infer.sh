#!/bin/bash
set -e
source semisynthetic/settings.sh

# Healthy cohort
# --------------

read_depth=$1
trial=$2
require_variable "read_depth" $read_depth
require_variable "trial" $trial


echo "[*] Performing MDSINE2 inference for semisynthetic (read_depth=${read_depth} | trial=${trial})"

inference_dir=${DATASET_DIR}/read_depth_${read_depth}/trial_${trial}/inference
mkdir -p ${inference_dir}

breadcrumb=${inference_dir}/inference.DONE
if [ -f ${breadcrumb} ]; then
  echo "[**] Inference already done."
  exit 0
fi


export MDSINE2_LOG_INI="${PROJECT_DIR}/scripts/analysis/logging_to_file.ini"
export LOG_FILEPATH="${inference_dir}/inference.log"
touch $LOG_FILEPATH

# ======= Fit NegBin qPCR model
echo "[*] Fitting Negative binomial model."
replicate_study=${DATASET_DIR}/read_depth_${read_depth}/trial_${trial}/replicate.pkl
mdsine2 infer-negbin --input ${replicate_study} --seed ${NEGBIN_SEED} --burnin 2000 --n-samples 6000 --checkpoint 200 --basepath $inference_dir

# ======= Perform inference
echo "[*] Performing inference on dataset."
synth_study=${DATASET_DIR}/read_depth_${read_depth}/trial_${trial}/dataset.pkl
output_study_name="simulated"
mdsine2 infer \
		--input $synth_study \
		--negbin $inference_dir/synthetic-replicate/mcmc.pkl \
		--seed $INFERENCE_SEED \
		--burnin $BURNIN \
		--n-samples $N_SAMPLES \
		--checkpoint $CHECKPOINT \
		--multiprocessing $MULTIPROCESSING \
		--rename-study "${output_study_name}" \
		--basepath $inference_dir \
		--interaction-ind-prior $INTERACTION_IND_PRIOR \
		--perturbation-ind-prior $PERTURBATION_IND_PRIOR

echo "[*] Visualizing output of ${output_study_name}"
export LOG_FILEPATH="${inference_dir}/visualization.log"
mdsine2 visualize-posterior \
		--chain ${inference_dir}/${output_study_name}/mcmc.pkl \
		--output-basepath ${inference_dir}/${output_study_name}/posterior

echo "[*] Finished inference."
touch ${breadcrumb}
