#!/bin/bash
set -e
source figures/settings.sh

# directory where the preprocessed pkl files containing time 0 and 0.5 are saved 
echo "Performing preprocessing on raw data files and generate pkl files containing time points 0 and 0.5"

# generate the preprocessed pickle files containing timepoints 0 and 0.5. 
# If files in $PREPROCESS_DIR contains timepoints 0 and 0.5, can ignore this code 
#python preprocess/helpers/preprocess.py \
#    --hamming-distance 2 \
#    --rename-prefix OTU \
#    --sequences ${PREPROCESS_DIR}/gibson_16S_rRNA_v4_ASV_seqs_aligned_filtered.fa \
#    --output-basepath ${PREPROCESSED_ALL_DIR} \
#    --max-n-species 2 \
#    --dataset_dir ${DATASET_DIR}

#python preprocess/helpers/assign_taxonomy_for_consensus_seqs.py \
#    --rdp-table ${PREPROCESS_DIR}/taxonomy_RDP.txt \
#    --confidence-threshold 50 \
#    --output-basepath ${PREPROCESSED_ALL_DIR}

echo "Generated relevant files. Now generating the figure"

python figures/figure2.py \
       -f1 "${PREPROCESSED_ALL_DIR}/gibson_healthy_agg_taxa.pkl" \
       -f2 "${PREPROCESSED_ALL_DIR}/gibson_inoculum_agg_taxa.pkl" \
       -f_loc "figures/figure2_files" \
       -de_names "healthy_order1" "healthy_order2" "healthy_order3" \
       -n 100 \
       -m "2" \
       -t 10 \
       -o "${FIGURES_OUTPUT_DIR}"

echo "Generating the deseq heatmaps"

python figures/deseq_heatmap_order.py \
    -loc "figures/figure2_files" \
    -abund "high" \
    -txt "abundant_species" \
    -taxo "order" \
    -o "mat_order_high" \
    -o_loc "${FIGURES_OUTPUT_DIR}"


python figures/deseq_heatmap_order.py \
    -loc "figures/figure2_files" \
    -abund "low" \
    -txt "abundant_species" \
    -taxo "order" \
    -o "mat_order_low" \
    -o_loc "${FIGURES_OUTPUT_DIR}"
