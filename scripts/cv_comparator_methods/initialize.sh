#!/bin/bash
# creates the files necessary to run inference on gibson dataset using comparator methods
set -e 
source cv_comparator_methods/settings.sh 

echo "[*] Preparing dataset for healthy data and saving to ${INPUT_DATASET_DIR}"

python cv_comparator_methods/helpers/create_biomass_table.py \
	-s "${PREPROCESS_DIR}/gibson_healthy_agg_filtered.pkl" \
	-o "${INPUT_DATASET_DIR}" \
	-on "abundance"

python cv_comparator_methods/helpers/format_data.py \
	-a "${INPUT_DATASET_DIR}/abundance.txt"\
	-m "${CLV_DIR}/data/gibson/basic_info/metadata.txt"\
	-o "${INPUT_DATASET_DIR}"

echo "[*] Done."
