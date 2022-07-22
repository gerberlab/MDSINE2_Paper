#!/bin/bash
set -e 
source cv_comparator_methods/settings.sh 

reg=$1

echo "Running gLV ${reg} regression and saving output to ${CV_OUTPUT_DIR}"
python ${CLV_DIR}/healthy_prediction_experiments_cv.py -m "glv" -r "${reg}" -o "${CV_OUTPUT_DIR}/results_abs_${reg}" -i "${INPUT_DATASET_DIR}"


