#!/bin/bash
set -e
source preprocess/settings.sh

require_dir ${REFERENCE_RDP_DIR}
require_file ${ASV_FASTA}

echo "[*] Aligning ASV sequences and performing phylogenetic placement."
python preprocess/helpers/place_seqs.py \
    --v4-region-start 1045 \
    --v4-region-end 1374 \
    --refpkg ${REFERENCE_RDP_DIR}\
    --query-reads ${ASV_FASTA} \
    --output-folder ${DATASET_DIR}/metadata_ASV \
    --temp-folder ${DATASET_DIR}/metadata_ASV/_tmp

echo "[*] Converting .sto to .fasta format."
python preprocess/helpers/sto_to_fasta.py \
-i output_ASVs/placed_sequences_on_v4_region.sto \
-o ${DATASET_DIR}/metadata_asv/aligned_asvs.fa
