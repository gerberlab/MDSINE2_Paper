#!/bin/bash
set -e
source figures/settings.sh
echo "Saving figure to ${FIGURES_OUTPUT_DIR}"

#bash figures/generate_figure_files.sh

python figures/figure2.py \
       -f1 "${PREPROCESSED_ALL_DIR}/gibson_healthy_agg.pkl" \
       -f2 "${PREPROCESS_DIR}/gibson_inoculum_agg.pkl" \
       -f_loc "${DATASET_DIR}/differential_abundance" \
       -de_names "healthy_of_window_1" "healthy_of_window_2" "healthy_of_window_3" \
       -n 100 \
       -m "2" \
       -t 10 \
       -o "${FIGURES_OUTPUT_DIR}"

echo "Generating the deseq heatmaps"

python figures/deseq_heatmap_order.py \
    -loc "${DATASET_DIR}/differential_abundance" \
    -abund "high" \
    -txt "abundant_species" \
    -de_name "healthy_of_window_"\
    -taxo "order" \
    -o "mat_order_high" \
    -o_loc "${FIGURES_OUTPUT_DIR}"


python figures/deseq_heatmap_order.py \
    -loc "${DATASET_DIR}/differential_abundance" \
    -abund "low" \
    -txt "abundant_species" \
    -de_name "healthy_of_window_"\
    -taxo "order" \
    -o "mat_order_low" \
    -o_loc "${FIGURES_OUTPUT_DIR}"

