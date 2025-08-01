#!/bin/bash
set -e
source cross_validation/settings.sh


echo "[*] Formatting data for regression code."
mkdir -p $REGRESSION_DATASET_DIR
env PYTHONPATH=../submodules/clv_fork \
  python cross_validation/helpers/format_data.py \
	-i ${DATASET_PKL} \
	-o ${REGRESSION_DATASET_DIR}

echo "[*] Done."
