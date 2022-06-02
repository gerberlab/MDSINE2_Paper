#!/bin/bash
set -e
source settings.sh


require_program "mdsine2"
require_program "date"

n_taxa=$1
model_name=$2
regression_type=$3

require_variable "n_taxa" $n_taxa
require_variable "model_name" $model_name
require_variable "regression_type" $regression_type


input_dataset_dir="${OTHER_METHODS_DATA_DIR}/trimmed_${n_taxa}"
output_dir="${OTHER_METHODS_OUTPUT_DIR}/${model_name}/taxa_top${n_taxa}"
inference_out_dir="${output_dir}/inference"
runtime_file="${inference_out_dir}/runtime.txt"

mkdir -p $inference_out_dir
echo "[*] Saving output to ${inference_out_dir}"

echo "[*] Running model ${model_name} with ${regression_type} for top "\
"${n_taxa} OTUs"
start_time=$(date +%s%N)  # nanoseconds

python ${OTHER_METHODS_CODE_DIR}/healthy_prediction_experiments.py \
    -m "${model_name}" \
    -r "${regression_type}" \
    -o "${inference_out_dir}" \
    -i "${input_dataset_dir}" \
    -t "${OTHER_METHODS_RUNTIMES}"

end_time=$(date +%s%N)
elapsed_time=$(( $(($end_time-$start_time)) / 1000000 ))
echo "[*] Finished mdsine2 inference. (${elapsed_time} ms)"
echo "${elapsed_time}" > $runtime_file




