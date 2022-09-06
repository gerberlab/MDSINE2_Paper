#!/bin/bash
set -e
source figures/settings.sh
echo "Saving figure to ${FIGURES_OUTPUT_DIR}"

bash figures/generate_figure_files.sh

python figures/figure2.py \
       -f1 "${PREPROCESSED_ALL_DIR}/gibson_healthy_agg.pkl" \
       -f2 "${PREPROCESS_DIR}/gibson_inoculum_agg.pkl" \
       -f_loc "${DATASET_DIR}/differential_abundance" \
       -de_name "healthy_of_window_"\
       -n 100 \
       -m "2" \
       -t 10 \
       -o "${FIGURES_OUTPUT_DIR}"

