#!/bin/bash
set -e
source preprocess/settings.sh

require_program python
echo "Agglomerating ASVs into OTUs"
echo "Using dataset found in ${DATASET_DIR}."
echo "Writing files into ${PREPROCESS_DIR}."


asv_aln_fasta=${DATASET_DIR}/metadata_ASV/aligned_asvs.fa
prefilt_aln_fasta=${DATASET_DIR}/metadata_ASV/prefiltered_asvs.fa


if ![ -f asv_aln_fasta ]; then
	echo "Couldn't locate pre-aligned ASV sequence file: ${asv_aln_fasta}"
	exit 1
fi


echo "[*] Pre-filtering ASVs by length ${PREFILT_LEN}"
python preprocess/helpers/prefilter_asvs.py \
	-i $asv_aln_fasta \
	-o $prefilt_aln_fasta \
	-t $PREFILT_LEN


echo "[*] Agglomerating ASVs into OTUs."
for dataset in healthy replicates inoculum; do
	echo "[*] Extracting dataset: ${dataset}"
	python preprocess/helpers/preprocess.py \
			--hamming-distance 2 \
			--rename-prefix OTU \
			--sequences $prefilt_aln_fasta \
			--remove-timepoints 0 0.5 \
			--max-n-species 2 \
			--dataset-name ${dataset} \
			--dataset-dir ${DATASET_DIR}/raw_tables \
			--output-prefix "${PREPROCESS_DIR}/gibson_${dataset}_agg"
done

echo "Done."
