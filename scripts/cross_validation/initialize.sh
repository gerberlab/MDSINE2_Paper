#!/bin/bash
set -e
source cross_validation/settings.sh


echo "[*] Preparing dataset for data (target: ${REGRESSION_DATASET_DIR})"
mkdir -p $REGRESSION_DATASET_DIR

python cross_validation/helpers/create_biomass_table.py \
	-s ${DATASET_PKL} \
	-o ${REGRESSION_DATASET_DIR} \
	-on "abundance"

python cross_validation/helpers/format_data.py \
	-a ${REGRESSION_DATASET_DIR}/abundance.txt \
	-m ${CLV_DIR}/data/gibson/basic_info/metadata.txt \
	-o ${REGRESSION_DATASET_DIR}

echo "[*] Done."
