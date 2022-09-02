#!/bin/bash
set -e
source figures/settings.sh

#Plot figure4
echo "Preprocess dir: ${PREPROCESS_DIR}"
python figures/figure4.py \
    --chain_healthy "${OUTPUT_DIR}/inference/healthy-seed0-fixed-cluster/mcmc.pkl"\
    --tree_fname "${DATASET_DIR}/metadata_OTU/newick_tree_query_reads.nhx"\
    --study_healthy "${PREPROCESS_DIR}/gibson_healthy_agg.pkl"\
    --study_inoc "${PREPROCESS_DIR}/gibson_inoculum_agg.pkl"\
    --detected_study_healthy "${PREPROCESSED_ALL_DIR}/gibson_healthy_agg_filtered3.pkl"\
    --output_loc "${FIGURES_OUTPUT_DIR}"
