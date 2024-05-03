#!/bin/bash
set -e
source semisynthetic/settings.sh


read_depth=$1
trial=$2
require_variable "read_depth" $read_depth
require_variable "trial" $trial



instance_dir=${DATASET_DIR}/read_depth_${read_depth}/trial_${trial}
dataset_file=${instance_dir}/dataset_filtered.pkl
regression_dir=${instance_dir}/regression
regression_input_dir=${regression_dir}/inputs
mkdir -p $regression_input_dir

echo "[*] Creating CLV files in ${regression_dir}"
python semisynthetic/regression/create_clv_inputs.py -i ${dataset_file} -o ${regression_input_dir}


run_regression() {
	model=$1
	regression_type=$2
	echo "[*] Running ${model}, ${regression_type} (reads=${read_depth} | trial=${trial})"
	out_dir=${regression_dir}/output/$model/$regression_type
	mkdir -p $out_dir
	python ${CLV_DIR}/healthy_prediction_experiments.py \
	  -m $model \
	  -r $regression_type \
	  -o $out_dir \
	  -i ${regression_input_dir} \
	  --limit_of_detection 1.0
}

# ======= Run regression methods
#run_regression "clv" "elastic_net"
#run_regression "lra" "elastic_net"
run_regression "glv" "elastic_net"
run_regression "glv" "ridge"
#run_regression "glv-ra" "elastic_net"
#run_regression "glv-ra" "ridge"
