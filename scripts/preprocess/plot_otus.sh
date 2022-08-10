#!/bin/bash
set -e
source preprocess/settings.sh


require_program python
top_n="200"

echo "Plotting the Aggregation of the OTUs."
echo "Output dir: ${PLOTS_OUT_DIR}."
echo "Only rendering the top ${top_n} OTUs."

# Plot the OTU aggregates
python preprocess/helpers/plot_otus.py \
    --study ${PREPROCESS_DIR}/gibson_healthy_agg.pkl \
    --outdir "${DATASET_DIR}/metadata_OTU/plots" \
    --top $top_n

echo "Done."
