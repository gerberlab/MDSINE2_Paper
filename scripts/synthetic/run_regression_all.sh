#!/bin/bash
set -e
source synthetic/settings.sh


for read_depth in 1000 25000; do
	for (( trial = 0; trial < ${NUM_SAMPLE_TRIALS}; trial++ )); do
		for noise_level in "low" "medium" "high"; do
			dataset_file=${DATASET_DIR}/data/reads_${read_depth}/trial_${trial}/subjset_${noise_level}.pkl
			inference_out_dir=${OUTPUT_DIR}/reads_${read_depth}/trial_${trial}/${noise_level}_noise/${model_name}/${regression_type}
			mkdir -p $inference_out_dir

			echo "[*] Saving output to ${inference_out_dir}"
			python synthetic/helpers/create_clv_inputs.py -i $dataset_file -o $inference_out_dir

			# ======= Run regression methods
			echo "[*] Running CLV, elastic-net (reads=${read_depth}, trial=${trial}, noise level=${noise_level})"
			python ${CLV_DIR}/healthy_prediction_experiments.py -m "clv" -r "elastic-net" -o $inference_out_dir -i $inference_out_dir

			echo "[*] Running LRA, elastic-net (reads=${read_depth}, trial=${trial}, noise level=${noise_level})"
			python ${CLV_DIR}/healthy_prediction_experiments.py -m "lra" -r "elastic-net" -o $inference_out_dir -i $inference_out_dir

			echo "[*] Running GLV, elastic-net (reads=${read_depth}, trial=${trial}, noise level=${noise_level})"
			python ${CLV_DIR}/healthy_prediction_experiments.py -m "glv" -r "elastic-net" -o $inference_out_dir -i $inference_out_dir

			echo "[*] Running GLV, ridge (reads=${read_depth}, trial=${trial}, noise level=${noise_level})"
			python ${CLV_DIR}/healthy_prediction_experiments.py -m "glv" -r "ridge" -o $inference_out_dir -i $inference_out_dir

			echo "[*] Running GLV-ra, elastic-net (reads=${read_depth}, trial=${trial}, noise level=${noise_level})"
			python ${CLV_DIR}/healthy_prediction_experiments.py -m "glv-ra" -r "elastic-net" -o $inference_out_dir -i $inference_out_dir

			echo "[*] Running GLV-ra, ridge (reads=${read_depth}, trial=${trial}, noise level=${noise_level})"
			python ${CLV_DIR}/healthy_prediction_experiments.py -m "glv-ra" -r "ridge" -o $inference_out_dir -i $inference_out_dir
		done
	done
done
