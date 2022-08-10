#!/bin/bash
set -e
source preprocess/settings.sh

phylo_out_dir=${OUTPUT_DIR}/phylogeny

require_program mdsine2
echo "Making the phylogenetic subtrees."
echo "Writing the files to ${phylo_out_dir}"

# Plot the phylogenetic subtrees for each OTU
mdsine2 render-phylogeny \
    --study ${PREPROCESS_DIR}/gibson_healthy_agg.pkl \
    --output-basepath ${phylo_out_dir} \
    --tree ${DATASET_DIR}/metadata_OTU/phylogenetic_tree_full_taxid.nhx \
    --seq-info ${DATASET_DIR}/metadata_OTU/RDP-11-5_BA_TS_info.tsv \
    --family-radius-factor 1.5 \
    --top 200

echo "Done."
