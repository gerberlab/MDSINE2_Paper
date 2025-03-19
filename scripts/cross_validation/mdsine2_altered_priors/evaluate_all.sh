#!/bin/bash
set -e
source cross_validation/settings.sh


evaluate_subdir()
{
  out_subdir=$1
  python cross_validation/helpers/prior_comparison_errors.py \
    --study ${DATASET_PKL} \
    --mdsine_outdir ${OUTPUT_DIR}/${out_subdir} \
    --out_dir ${OUTPUT_DIR}/${out_subdir} \
    --subsample_every 100
}


evaluate_regression_subdir()
{
  python cross_validation/helpers/regression_errors.py \
  --study ${DATASET_PKL} \
  --regression_inputs_dir ${REGRESSION_DATASET_DIR} \
  --regression_outdir /data/local/MDSINE2_files_other/regression/regression_glv_elastic-net \
  --out_dir /data/local/MDSINE2_files_other/regression/regression_glv_elastic-net \
  --model_name 'glv-elastic-net'
}


evaluate_subdir "default_new"
evaluate_subdir "edge_density_1"
evaluate_subdir "edge_density_2"
evaluate_subdir "edge_density_3"
evaluate_subdir "edge_density_4"
evaluate_subdir "edge_density_5"
evaluate_subdir "growth_si_var_1"
evaluate_subdir "growth_si_var_2"
evaluate_subdir "growth_si_var_3"
evaluate_subdir "growth_si_var_4"
evaluate_subdir "interaction_var_1"
evaluate_subdir "interaction_var_2"
evaluate_subdir "interaction_var_3"
evaluate_subdir "interaction_var_4"
evaluate_subdir "pert_var_1"
evaluate_subdir "pert_var_2"
evaluate_subdir "pert_var_3"
evaluate_subdir "pert_var_4"
evaluate_regression_subdir


