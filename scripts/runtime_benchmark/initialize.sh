#!/bin/bash
set -e
source runtime_benchmark/settings.sh


require_program python

cd runtime_benchmark/helpers
study=${DATASET_DIR}/preprocessed/gibson_healthy_agg.pkl
replicate_study=${DATASET_DIR}/preprocessed/gibson_replicates_agg_taxa.pkl

for n_taxa in 10 25 50 100; do
	echo "[*] Preparing dataset of ${n_taxa} taxa..."
	target_dataset=${DATASET_DIR}/trimmed_${n_taxa}/dataset.pkl
	target_replicate=${DATASET_DIR}/trimmed_${n_taxa}/replicates.pkl

	python select_topk_bugs.py -f ${study} -n $n_taxa -o $target_dataset
	python filter_replicates_like_other_dataset.py \
			--replicate-dataset $replicate_study \
			--like-other $target_dataset \
			--out-path $target_replicate
done
echo "[*] Done."
