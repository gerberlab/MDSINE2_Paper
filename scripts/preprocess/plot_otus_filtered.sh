#!/bin/bash
set -e
data_modality=$1
if [ "${data_modality}" == "healthy" ]; then
  source preprocess/settings_healthy.sh
elif [ "${data_modality}" == "uc" ]; then
  source preprocess/settings_uc.sh
else
  echo "data_modality argument is required and must be either 'healthy' or 'uc'. Exiting."
  exit 1
fi


require_program python
top_n="200"

plots_out_dir="${DATASET_DIR}/plots"
echo "Plotting the Aggregation of the OTUs."
echo "Output dir: ${plots_out_dir}."
echo "Only rendering the top ${top_n} OTUs."

# Plot the OTU aggregates
python preprocess/helpers/plot_otus.py \
    --study ${PREPROCESS_DIR}/gibson_healthy_agg_filtered.pkl \
    --outdir "${plots_out_dir}" \
    --top $top_n

echo "Done."
