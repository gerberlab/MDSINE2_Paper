#!/bin/bash
set -e
source semisynthetic2/settings.sh


out_dir=$1
input_dir=$2
model=$3
regularization=$4
require_variable "out_dir" "${out_dir}"
require_variable "input_dir" "${input_dir}"
require_variable "model" "${model}"
require_variable "regularization" "${regularization}"


set -e
echo "[* run_regression.sh] Running regression: ${model}, ${regularization}"
mkdir -p "$out_dir"
python ${CLV_DIR}/healthy_prediction_experiments.py \
  -m "$model" \
  -r "$regularization" \
  -o "$out_dir" \
  -i "${input_dir}" \
  --limit_of_detection 1.0
