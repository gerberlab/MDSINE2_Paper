#!/bin/bash
set -e
source runtime_benchmark/settings.sh


require_program "mdsine2"
require_program "date"

for n_taxa in 10 25 50 100 125 150 175 200; do
	for (( trial = 1; trial < ${NUM_TRIALS}+1; trial++ )); do
		seed=$((n_taxa * trial))
		dataset=${DATASET_DIR}/trimmed_${n_taxa}/replicates.pkl
		negbin_out_dir=${OUTPUT_DIR}/taxa_top_${n_taxa}/trial_${trial}/mdsine2_negbin
		mkdir -p $negbin_out_dir

		runtime_file=${negbin_out_dir}/runtime.txt
		seed_file=${negbin_out_dir}/seed.txt
		echo "${seed}" > $seed_file

		echo "[*] Running negative-binomial inference (n_taxa=${n_taxa}, trial=${trial}; using seed=${seed})"
		start_time=$(date +%s%N)  # nanoseconds

		mdsine2 infer-negbin --input $dataset --seed ${seed} --burnin 2000 --n-samples 6000 --checkpoint 200 --basepath $negbin_out_dir

		end_time=$(date +%s%N)
		elapsed_time=$(( $(($end_time-$start_time)) / 1000000 ))
		echo "[*] Finished negative-binomial inference. (${elapsed_time} ms)"
		echo "${elapsed_time}" > $runtime_file

		echo "[*] Drawing negative-binomial visualization."
		mdsine2 visualize-negbin \
				--chain "${negbin_out_dir}/replicates/mcmc.pkl" \
				--output-basepath "${negbin_out_dir}/replicates/posterior"
		echo "[*] Finished negative-binomial visualization."
	done
done
