#!/bin/bash
set -e
source preprocess/settings.sh

require_program python
echo "Assigning Taxonomy for OTUs..."


rdp_file=${DATASET_DIR}/metadata_ASV/taxonomy_RDP.txt
if ![ -f $rdp_file ]; then
	echo "Couldn't locate taxonomy file ${rdp_file}"
	exit 1
fi


for dataset in healthy replicates inoculum; do
	python preprocess/helpers/assign_taxonomy_for_consensus_seqs.py \
	--rdp-table $rdp_file \
	--confidence-threshold 50 \
	-i ${PREPROCESS_DIR}/gibson_${dataset}_agg.pkl \
	-o ${PREPROCESS_DIR}/gibson_${dataset}_agg.pkl \

echo "Done."
