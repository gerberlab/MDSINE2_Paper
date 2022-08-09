#!/bin/bash
set -e
source preprocess/settings.sh

require_program python
echo "Agglomerating ASVs into OTUs"
echo "Using dataset found in ${DATASET_DIR}."
echo "Writing files into ${PREPROCESS_DIR}."

# Agglomerate ASVs into OTUs
for dataset in healthy replicate inoculum; do
	echo "[*] Extracting dataset: ${dataset}"
	python preprocess/helpers/preprocess.py \
			--hamming-distance 2 \
			--rename-prefix OTU \
			--sequences ${PREPROCESS_DIR}/gibson_16S_rRNA_v4_ASV_seqs_aligned_filtered.fa \
			--output-basepath ${PREPROCESS_DIR} \
			--remove-timepoints 0 0.5 \
			--max-n-species 2 \
			--dataset-name ${dataset} \
			--dataset-dir ${DATASET_DIR}
done

echo "Done."
