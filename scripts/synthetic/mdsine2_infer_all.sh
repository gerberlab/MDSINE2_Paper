#!/bin/bash
set -e
source synthetic/settings.sh

for read_depth in 1000 25000; do
	for (( trial = 0; trial < ${NUM_SAMPLE_TRIALS}; trial++ )); do
		for noise_level in "low" "medium" "high"; do
			dataset=${DATASET_DIR}/data/reads_${read_depth}/trial_${trial}/subjset_${noise_level}.pkl
			trial_output_dir=${OUTPUT_DIR}/reads_${read_depth}/trial_${trial}/${noise_level}_noise
			inference_out_dir=${trial_output_dir}/inference
			mkdir -p $inference_out_dir

			# ======= Run inference (and compute runtime)
			echo "[*] Running mdsine2 inference (trial=${trial}, noise level=${noise_level})"
			python synthetic/helpers/inference.py \
					--input $dataset \
					--negbin ${NEGBIN_A0} ${NEGBIN_A1} \
					--seed 0 \
					--burnin 5000 \
					--n-samples 15000 \
					--checkpoint 1000 \
					--multiprocessing 0 \
					--basepath $inference_out_dir \
					--interaction-ind-prior "strong-sparse" \
					--perturbation-ind-prior "weak-agnostic"
			echo "[*] Finished mdsine2 inference."

			# ======= Draw visualizations.
			echo "[*] Drawing mdsine2 posterior visualization."
			mkdir -p $inference_out_dir/posterior
			mdsine2 visualize-posterior \
					--chain $inference_out_dir/mcmc.pkl \
					--output-basepath $inference_out_dir/posterior
			echo "[*] Finished mdsine2 posterior visualization."
		done
	done
done
