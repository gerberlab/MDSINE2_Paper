#!/bin/bash
set -e
source figures/settings.sh

python figures/functional_enrichment_plot.py \
    --base_loc "${ENRICHMENT_OUTPUT_DIR}" \
    --module_result_pkl df_kegg_modules \
    --pathways_result_pkl df_kegg_pathways \
    --cazyme_result_pkl df_cazymes \
    --output_loc "${FIGURES_OUTPUT_DIR}"\
    --mcmc_file_loc "${OUTPUT_DIR}/mdsine2/inference/healthy-seed0-fixed-cluster/mcmc.pkl"\

