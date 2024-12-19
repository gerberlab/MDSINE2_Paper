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


# A wrapper around the DADA2 R script.

require_program Rscript
require_dir ${READS_DIR}

echo "Running DADA2 on the raw reads."
Rscript --vanilla preprocess/helpers/dada2.R \
-r ${READS_DIR} \
-o ${DATASET_DIR}/raw_tables

echo "Done."
