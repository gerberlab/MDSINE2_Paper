#!/bin/bash

set -e
source cv_comparator_methods/settings.sh

model=$1
reg=$2

echo "Running ${model}, regression ${reg}. Output dir=${CV_OUTPUT_DIR}"
python ${CLV_DIR}/healthy_prediction_experiments_cv.py \
		-m "${model}" \
		-r "${reg}" \
		-o "${OUTPUT_DIR}/regression_${model}_${reg}" \
		-i "${REGRESSION_DATASET_DIR}"
