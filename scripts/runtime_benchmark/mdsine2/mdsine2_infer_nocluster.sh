#!/bin/bash
set -e
source runtime_benchmark/settings.sh


require_program "mdsine2"
require_program "date"

n_taxa=$1
trial=$2

require_variable "n_taxa" $n_taxa
require_variable "trial" $trial


dataset=${DATASET_DIR}/trimmed_${n_taxa}/dataset.pkl

trial_dir=${OUTPUT_DIR}/taxa_top_${n_taxa}/trial_${trial}
negbin=${trial_dir}/mdsine2_negbin/replicates/mcmc.pkl
inference_out_dir=${trial_dir}/mdsine2_nocluster
mkdir -p $inference_out_dir

runtime_file=${inference_out_dir}/runtime.txt

# ======= Generate and report seed
seed=$((n_taxa * trial + 1))
seed_file=${inference_out_dir}/seed.txt
echo "${seed}" > $seed_file

# ======= Run inference (and compute runtime)
echo "[*] Running mdsine2 inference (n_taxa=${n_taxa}, trial=${trial}; using seed=${seed})"
start_time=$(date +%s%N)  # nanoseconds

python runtime_benchmark/helpers/mdsine2_nocluster.py \
		--input $dataset \
		--negbin $negbin \
		--seed ${seed} \
		--burnin 5000 \
		--n-samples 15000 \
		--checkpoint 1000 \
		--basepath $inference_out_dir \
		--interaction-ind-prior "strong-sparse" \
		--perturbation-ind-prior "weak-agnostic"

end_time=$(date +%s%N)
elapsed_time=$(( $(($end_time-$start_time)) / 1000000 ))
echo "[*] Finished mdsine2 (NON-CLUSTERED) inference. (${elapsed_time} ms)"
echo "${elapsed_time}" > $runtime_file
