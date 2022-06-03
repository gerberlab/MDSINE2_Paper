#!/bin/bash
set -e
source runtime_benchmark/settings.sh


for n_taxa in 10 25 50 100; do
	for (( trial = 1; trial < ${MDSINE2_NUM_TRIALS}+1; trial++ )); do
		dataset=${DATASET_DIR}/data/trial_${trial}/subjset_${noise_level}.pkl
		trial_dir=${OUTPUT_DIR}/trial_${trial}/${noise_level}_noise
		inference_out_dir=${trial_dir}/inference
		mkdir -p $inference_out_dir

		# ======= Run inference (and compute runtime)
		echo "[*] Running mdsine2 inference (trial=${trial}, noise level=${noise_level})"
		mdsine2 infer \
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
		mdsine2 visualize-posterior \
				--chain $inference_out_dir/mcmc.pkl \
				--output-basepath $out_dir/posterior
		echo "[*] Finished mdsine2 posterior visualization."
	done
done
