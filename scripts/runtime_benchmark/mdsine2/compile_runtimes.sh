#!/bin/bash
set -e
source runtime_benchmark/settings.sh

echo "[*] Compiling runtimes."

echo "MethodName,NumTaxa,Trial,Seed,Runtime" > $csv_path  # Header line

for n_taxa in 10 25 50 100; do
	for (( trial = 1; trial < ${N_TRIALS}+1; trial++ )); do
		trial_dir=${MDSINE2_OUTPUT_DIR}/taxa_top_${n_taxa}/trial_${trial}
		elapsed_time=$(cat ${trial_dir}/inference/inference_runtime.txt)
		seed=$(cat ${trial_dir}/seed.txt)

		echo "n_taxa=${n_taxa}, trial=${trial} -> ${elapsed_time} ms"
		echo "MDSINE2,${n_taxa},${trial},${seed},${elapsed_time}" >> $csv_path
done
