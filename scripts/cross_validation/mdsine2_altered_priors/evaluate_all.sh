#!/bin/bash
set -e
source cross_validation/settings.sh


evaluate_mdsine2_subdir()
{
  out_subdir=$1
  python cross_validation/helpers/mdsine2_errors.py \
    --study ${DATASET_PKL} \
    --mdsine_outdir ${OUTPUT_DIR}/other_priors/${out_subdir} \
    --out_dir ${OUTPUT_DIR}/other_priors/${out_subdir} \
    --subsample_every 100
}


evaluate_mdsine2_subdir "default_new"
evaluate_mdsine2_subdir "edge_density_1"
evaluate_mdsine2_subdir "edge_density_2"
evaluate_mdsine2_subdir "edge_density_3"
evaluate_mdsine2_subdir "edge_density_4"
evaluate_mdsine2_subdir "edge_density_5"
evaluate_mdsine2_subdir "growth_si_var_1"
evaluate_mdsine2_subdir "growth_si_var_2"
evaluate_mdsine2_subdir "growth_si_var_3"
evaluate_mdsine2_subdir "growth_si_var_4"
evaluate_mdsine2_subdir "interaction_var_1"
evaluate_mdsine2_subdir "interaction_var_2"
evaluate_mdsine2_subdir "interaction_var_3"
evaluate_mdsine2_subdir "interaction_var_4"
evaluate_mdsine2_subdir "pert_var_1"
evaluate_mdsine2_subdir "pert_var_2"
evaluate_mdsine2_subdir "pert_var_3"
evaluate_mdsine2_subdir "pert_var_4"


