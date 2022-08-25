#!/bin/bash
set -e
source preprocess/settings.sh

require_dir ${REFERENCE_RDP_DIR}
require_file ${OTU_FASTA}

echo "[*] Performing phylogenetic placement of OTUs."
python scripts/place_seqs.py \
    --v4-region-start 1045 \
    --v4-region-end 1374 \
    --refpkg ${REFERENCE_RDP_DIR} \
    --query-reads ${OTU_FASTA} \
    --output-folder ${DATASET_DIR}/metadata_OTU \
    --temp-folder ${DATASET_DIR}/metadata_OTU/_tmp
