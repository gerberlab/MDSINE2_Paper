#!/bin/bash
set -e
source preprocess/settings.sh

require_program python
echo "Assigning Taxonomy for OTUs..."


rdp_file=${DATASET_DIR}/raw_tables/rdp_species.tsv
if ! [ -f $rdp_file ]; then
	echo "Couldn't locate taxonomy file ${rdp_file}"
	exit 1
fi
echo "[*] This script assumes each OTU is just an ASV. (Todo: deprecate OTU agglomeration.)"


for dataset in healthy replicates inoculum; do
	python preprocess/helpers/assign_taxonomy_for_consensus_seqs.py \
	--rdp-table $rdp_file \
	--confidence-threshold 50 \
	-i ${PREPROCESS_DIR}/gibson_${dataset}_agg.pkl \
	-o ${PREPROCESS_DIR}/gibson_${dataset}_agg.pkl
done

echo "Done."
