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
echo "Assigning Taxonomy for OTUs..."


rdp_file=${DATASET_DIR}/raw_tables/rdp_species.tsv
if ! [ -f $rdp_file ]; then
	echo "Couldn't locate taxonomy file ${rdp_file}"
	exit 1
fi
echo "[*] This script assumes each OTU is just an ASV. (Todo: deprecate OTU agglomeration.)"


for dataset in ${data_modality} replicates inoculum; do
	python preprocess/helpers/assign_taxonomy_for_consensus_seqs.py \
	--rdp-table $rdp_file \
	--confidence-threshold 50 \
	-i ${PREPROCESS_DIR}/gibson_${dataset}_agg.pkl \
	-o ${PREPROCESS_DIR}/gibson_${dataset}_agg.pkl
done

echo "Done."
