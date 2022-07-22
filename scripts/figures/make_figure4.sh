#!/bin/bash
set -e
source figures/settings.sh
: '
echo "Generating Study object to represent info about taxa that are consistently detected"
mdsine2 filter \
    --dataset ${PREPROCESSED_ALL_DIR}/gibson_healthy_agg_taxa.pkl \
    --outfile ${PREPROCESSED_ALL_DIR}/gibson_healthy_agg_taxa_filtered3.pkl \
    --dtype rel \
    --threshold 0.0001 \
    --min-num-consecutive 3 \
    --min-num-subjects 1 \
    --colonization-time 0
'

#Plot figure4
python figures/figure4.py \
    --chain_healthy "/Users/microbiome/Desktop/figures_code/healthy-seed0-mixed/mcmc.pkl" \
    --tree_fname "/Users/microbiome/Desktop/figures_code/phylogenetic_placement_OTUs/phylogenetic_tree_only_query.nhx" \
    --study_healthy "${PREPROCESS_DIR}/gibson_healthy_agg_taxa.pkl"  \
    --study_inoc "${PREPROCESS_DIR}/gibson_inoculum_agg_taxa.pkl"\
    --detected_study_healthy "${PREPROCESSED_ALL_DIR}/gibson_healthy_agg_taxa_filtered3.pkl" \
    --keystoneness_file "/Users/microbiome/Desktop/figures_code/keystoneness/healthy_fwsim_day20.h5"\
    --output_loc "${FIGURES_OUTPUT_DIR}"
