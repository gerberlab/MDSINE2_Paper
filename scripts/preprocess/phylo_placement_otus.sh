#!/bin/bash
set -e
source preprocess/settings.sh

require_dir ${REFERENCE_RDP_DIR}
require_file ${OTU_FASTA}

echo "[*] Performing phylogenetic placement of OTUs."
outdir=${DATASET_DIR}/metadata_OTU

bash preprocess/place_seqs.sh \
${REFERENCE_RDP_DIR} \
"RDP-11-5_TS_Processed" \
${outdir} \
${OTU_FASTA} \
1045 1374
