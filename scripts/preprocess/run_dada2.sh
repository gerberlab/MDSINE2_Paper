#!/bin/bash
set -e
source preprocess/settings.sh


# A wrapper around the DADA2 R script.

require_program Rscript
require_dir ${READS_DIR}

echo "Running DADA2 on the raw reads."
Rscript --vanilla preprocess/helpers/dada2.R \
-r ${READS_DIR} \
-o ${DATASET_DIR}/raw_tables

echo "Done."
