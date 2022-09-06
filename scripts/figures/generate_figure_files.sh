#!/bin/bash

set -e
source figures/settings.sh

# directory where the preprocessed pkl files containing time 0 and 0.5 are saved
echo "Preprocessing raw data files and generating pkl files containing time points 0 and 0.5. The pkl files are used to make figure 2"
echo ""
echo "Making directory to keep Study pkl files containing time-points 0 and 0.5"
mkdir -p ${PREPROCESSED_ALL_DIR}

dataset="healthy"

python preprocess/helpers/preprocess.py \
    --hamming-distance 3 \
    --rename-prefix OTU \
    --sequences ${DATASET_DIR}/metadata_ASV/prefiltered_asvs.fa \
    --max-n-species 2 \
    --dataset-dir ${DATASET_DIR}/raw_tables \
    --dataset-name ${dataset} \
    --output-prefix "${PREPROCESSED_ALL_DIR}/gibson_${dataset}_agg" \
    --trim-option ALL_GAPS

rdp_file="${DATASET_DIR}/metadata_OTU/taxonomy_RDP.txt"
python preprocess/helpers/assign_taxonomy_for_consensus_seqs.py \
    --rdp-table ${rdp_file} \
    --confidence-threshold 50 \
    --output-study ${PREPROCESSED_ALL_DIR}/gibson_${dataset}_agg.pkl \
    --input-study ${PREPROCESSED_ALL_DIR}/gibson_${dataset}_agg.pkl
    
    
echo "Generating Study objects containing info about taxa that are consistently detected. Used to make figure 4"
mdsine2 filter \
    --dataset ${PREPROCESSED_ALL_DIR}/gibson_healthy_agg.pkl \
    --outfile ${PREPROCESSED_ALL_DIR}/gibson_healthy_agg_filtered3.pkl \
    --dtype rel \
    --threshold 0.0001 \
    --min-num-consecutive 3 \
    --min-num-subjects 1 \
    --colonization-time 0

echo "Generated relevant files. Files saved to ${PREPROCESSED_ALL_DIR}"

