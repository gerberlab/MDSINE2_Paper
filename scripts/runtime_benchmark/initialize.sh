#!/bin/bash
set -e
source runtime_benchmark/settings.sh


require_program python

cd runtime_benchmark/helpers
study=${DATASET_DIR}/preprocessed/gibson_healthy_agg_taxa_filtered.pkl
replicate_study=${DATASET_DIR}/preprocessed/gibson_replicates_agg_taxa.pkl

for n_taxa in 10 25 50 100; do
	target_dir=${DATASET_DIR}/trimmed_${n_taxa}

	echo "[*] Preparing dataset of ${n_taxa} taxa... (target dir = ${target_dir})"
	mkdir -p ${target_dir}

	target_dataset=${target_dir}/dataset.pkl
	target_replicate=${target_dir}/replicates.pkl

	python select_topk_bugs.py -f ${study} -n $n_taxa -o $target_dataset
	python filter_replicates_like_other_dataset.py \
			--replicate-dataset $replicate_study \
			--like-other $target_dataset \
			--out-path $target_replicate
	python create_biomass_table.py \
	       -s ${DATASET_DIR}/trimmed_${n_taxa}/dataset \
	       -o "${target_dir}" \
	       -on "abundance" \
	       -n $n_taxa
	python format_data.py \
	       -a "${target_dir}/abundance.txt"\
	       -m "${REGRESSION_DATA_DIR}/basic_info/metadata.txt"\
	       -o "${target_dir}"
done
echo "[*] Done."
