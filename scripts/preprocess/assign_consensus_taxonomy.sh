#!/bin/bash
set -e
source settings.sh

require_program python
echo "Assigning Taxonomy for OTUs..."

# Assign taxonomy for OTUs
python preprocess/helpers/assign_taxonomy_for_consensus_seqs.py \
    --rdp-table ${PREPROCESS_DIR}/taxonomy_RDP.txt \
    --confidence-threshold 50 \
    --output-basepath ${PREPROCESS_DIR}

echo "Done."
