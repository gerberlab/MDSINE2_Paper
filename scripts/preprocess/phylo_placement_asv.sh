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

require_dir ${REFERENCE_RDP_DIR}
require_file ${ASV_FASTA}

outdir=${DATASET_DIR}/metadata_ASV

echo "[*] Aligning ASV sequences and performing phylogenetic placement."
mkdir -p ${outdir}
bash preprocess/place_seqs.sh \
  ${data_modality} \
  ${REFERENCE_RDP_DIR} \
  "RDP-11-5_TS_Processed" \
  ${outdir} \
  ${ASV_FASTA} \
  1045 1374


echo "[*] Converting .sto to .fasta format."
python preprocess/helpers/phylo_placement/sto_to_fasta.py \
  -i ${outdir}/aligned_sequences.sto \
  -o ${outdir}/aligned_asvs.fa
