#!/bin/bash
set -e
source synthetic/settings.sh


run_regression() {
	model=$1
	regression_type=$2
	echo "[*] Running ${model}, ${regression_type} (reads=${read_depth}, trial=${trial}, noise level=${noise_level})"
	out_dir=$intermediate_dir/$model/$regression_type
	mkdir -p $out_dir
	python ${CLV_DIR}/healthy_prediction_experiments.py -m $model -r $regression_type -o $out_dir -i $intermediate_dir --limit_of_detection 1.0
}


for read_depth in 1000 25000; do
	for (( trial = 0; trial < ${NUM_SAMPLE_TRIALS}; trial++ )); do
		for noise_level in "low" "medium" "high"; do
			dataset_file=${DATASET_DIR}/data/reads_${read_depth}/trial_${trial}/subjset_${noise_level}.pkl
			intermediate_dir=${OUTPUT_DIR}/reads_${read_depth}/trial_${trial}/${noise_level}_noise
			mkdir -p $intermediate_dir

			echo "[*] Creating CLV files in ${intermediate_dir}"
			python synthetic/helpers/create_clv_inputs.py -i $dataset_file -o $intermediate_dir

			# ======= Run regression methods
			run_regression "clv" "elastic_net"
			run_regression "lra" "elastic_net"
			run_regression "glv" "elastic_net"
			run_regression "glv" "ridge"
			run_regression "glv-ra" "elastic_net"
			run_regression "glv-ra" "ridge"
		done
	done
done
