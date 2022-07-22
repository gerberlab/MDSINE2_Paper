#!/bin/bash

set -e
source figures/settings.sh
echo "Saving figure to ${FIGURES_OUTPUT_DIR}"

python figures/figure3.py \
    --mdsine_path "/Users/microbiome/Desktop/figures_code/cv_mixed/forward_sims" \
    --non_glv_elas_path "${FW_SIM_OUTPUT_DIR}/results_rel_elastic-net"\
    --non_glv_ridge_path "${FW_SIM_OUTPUT_DIR}/results_rel_ridge"\
    --glv_elas_path "${FW_SIM_OUTPUT_DIR}/results_abs_elastic-net"\
    --glv_ridge_path "${FW_SIM_OUTPUT_DIR}/results_abs_ridge"\
    --study_file "${PREPROCESS_ALL_DIR}/gibson_healthy_agg_taxa.pkl"\
    --output_path "${FIGURES_OUTPUT_DIR}"\
    --output_name "figure3"\
    --rel_abund_lim 1e-6

