#!/bin/bash

set -e
source figures/settings.sh

python figures/enrichment_analysis.py \
    -loc1 "${OUTPUT_DIR}/inference/healthy-seed0-fixed-cluster/mcmc.pkl" \
    -o_loc "${FIGURES_OUTPUT_DIR}"

