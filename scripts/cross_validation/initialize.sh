#!/bin/bash
set -e
source cross_validation/settings.sh


echo "[*] Preparing dataset for data (target: ${REGRESSION_DATASET_DIR})"
mkdir -p $REGRESSION_DATASET_DIR

python cv_comparator_methods/helpers/create_biomass_table.py \
	-s "${PREPROCESS_DIR}/gibson_healthy_agg_filtered.pkl" \
	-o "${REGRESSION_DATASET_DIR}" \
	-on "abundance"

python cv_comparator_methods/helpers/format_data.py \
	-a "${REGRESSION_DATASET_DIR}/abundance.txt"\
	-m "${CLV_DIR}/data/gibson/basic_info/metadata.txt"\
	-o "${REGRESSION_DATASET_DIR}"

echo "[*] Done."
