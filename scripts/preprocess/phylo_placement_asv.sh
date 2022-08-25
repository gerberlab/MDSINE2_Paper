#!/bin/bash
set -e
source preprocess/settings.sh

require_dir ${REFERENCE_RDP_DIR}
require_file ${ASV_FASTA}

outdir=${DATASET_DIR}/metadata_ASV

echo "[*] Aligning ASV sequences and performing phylogenetic placement."
python preprocess/place_seqs.sh \
${REFERENCE_RDP_DIR} \
"RDP-11-5_TS_Processed" \
${outdir} \
${ASV_FASTA} \
1045 1374


echo "[*] Converting .sto to .fasta format."
python preprocess/helpers/sto_to_fasta.py \
-i ${outdir}/aligned_sequences.sto \
-o ${outdir}/aligned_asvs.fa
