#!/bin/bash
set -e
source settings.sh


require_program "mdsine2"
require_program "date"

n_taxa=$1
trial=$2
model_name=$3
regression_type=$4

require_variable "n_taxa" $n_taxa
require_variable "model_name" $model_name
require_variable "regression_type" $regression_type


input_dataset_dir="${DATASET_DIR}/trimmed_${n_taxa}"
trial_dir="${OUTPUT_DIR}/taxa_top_${n_taxa}/trial_${trial}"

output_dir="${trial_dir}/${model_name}/${regression_type}"
runtime_file="${output_dir}/runtime.txt"

mkdir -p $output_dir
echo "[*] Saving output to ${output_dir}"

echo "[*] Running model ${model_name} with ${regression_type} for top "\
"${n_taxa} OTUs"
start_time=$(date +%s%N)  # nanoseconds

python ${REGRESSION_CODE_DIR}/healthy_prediction_experiments.py \
    -m "${model_name}" \
    -r "${regression_type}" \
    -o "${output_dir}" \
    -i "${input_dataset_dir}"

end_time=$(date +%s%N)
elapsed_time=$(( $(($end_time-$start_time)) / 1000000 ))
echo "[*] Finished mdsine2 inference. (${elapsed_time} ms)"
echo "${elapsed_time}" > $runtime_file




