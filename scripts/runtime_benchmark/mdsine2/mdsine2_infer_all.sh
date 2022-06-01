#!/bin/bash
set -e
source runtime_benchmark/settings.sh


require_program "mdsine2"
require_program "date"

for n_taxa in 10 25 50 100; do
	for (( trial = 1; trial < ${MDSINE2_NUM_TRIALS}+1; trial++ )); do
		dataset=${DATASET_DIR}/trimmed_${n_taxa}/dataset.pkl

		trial_dir=${MDSINE2_OUTPUT_DIR}/taxa_top_${n_taxa}/trial_${trial}
		negbin=${trial_dir}/negbin/replicates/mcmc.pkl

		inference_out_dir=${trial_dir}/inference
		runtime_file=${inference_out_dir}/inference_runtime.txt

		# ======= Generate and report seed
		seed=$((n_taxa * trial + 1))
		seed_file=${inference_out_dir}/seed.txt
		echo "${seed}" > $seed_file

		# ======= Run inference (and compute runtime)
		echo "[*] Running mdsine2 inference (n_taxa=${n_taxa}, trial=${trial}; using seed=${seed})"
		start_time=$(date +%s%N)  # nanoseconds

		mdsine2 infer \
				--input $dataset \
				--negbin $negbin \
				--seed 0 \
				--burnin 5000 \
				--n-samples 15000 \
				--checkpoint 1000 \
				--multiprocessing 0 \
				--basepath $out_dir \
				--interaction-ind-prior "strong-sparse" \
				--perturbation-ind-prior "weak-agnostic"

		end_time=$(date +%s%N)
		elapsed_time=$(( $(($end_time-$start_time)) / 1000000 ))
		echo "[*] Finished mdsine2 inference. (${elapsed_time} ms)"
		echo "${elapsed_time}" > $runtime_file

		# ======= Draw visualizations.
		echo "[*] Drawing mdsine2 posterior visualization."
		mdsine2 visualize-posterior \
				--chain $out_dir/mcmc.pkl \
				--output-basepath $out_dir/posterior
		echo "[*] Finished mdsine2 posterior visualization."
done
