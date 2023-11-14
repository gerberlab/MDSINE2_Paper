#!/bin/bash
set -e
source cross_validation/settings.sh


evaluate_subdir()
{
  out_subdir=$1
  python cross_validation/helpers/prior_comparison_errors.py \
    --study ${DATASET_PKL} \
    --mdsine2_outdir ${OUTPUT_DIR}/${out_subdir} \
    --out_dir ${OUTPUT_DIR}/${out_subdir} \
    --subsample_every 100
}


evaluate_subdir "mdsine2-modules-interaction-density-dense"
evaluate_subdir "mdsine2-modules-interaction-density-sparse"
evaluate_subdir "mdsine2-modules-interaction-density-medium"
evaluate_subdir "mdsine2-modules-interaction-strength-negative"
evaluate_subdir "mdsine2-modules-interaction-strength-positive"
evaluate_subdir "mdsine2-modules-interaction-strength-zero"
