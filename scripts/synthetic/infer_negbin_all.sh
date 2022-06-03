#!/bin/bash
set -e
source synthetic/settings.sh


require_program "mdsine2"

for (( trial = 0; trial < ${NUM_SAMPLE_TRIALS}; trial++ )); do
	for noise_level in "low" "medium" "high"; do
		seed=0

		dataset=${DATASET_DIR}/trial_${trial}/subjset_${noise_level}.pkl
		trial_dir=${OUTPUT_DIR}/trial_${trial}/${noise_level}_noise
		negbin_out_dir=${trial_dir}/negbin
		mkdir -p $negbin_out_dir

		echo "[*] Running negative-binomial inference (trial=${trial}, noise level=${noise_level})"
		mdsine2 infer-negbin --input $dataset --seed 0 --burnin 2000 --n-samples 6000 --checkpoint 200 --basepath $negbin_out_dir
		echo "[*] Finished negative-binomial inference."

		echo "[*] Drawing negative-binomial visualization."
		mdsine2 visualize-negbin \
				--chain "${negbin_out_dir}/replicates/mcmc.pkl" \
				--output-basepath "${negbin_out_dir}/replicates/posterior"
		echo "[*] Finished negative-binomial visualization."
	done
done
