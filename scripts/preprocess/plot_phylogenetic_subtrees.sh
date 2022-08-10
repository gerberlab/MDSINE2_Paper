#!/bin/bash
set -e
source preprocess/settings.sh


require_program mdsine2
echo "Making the phylogenetic subtrees."
echo "Writing the files to ${PHYLOGENY_OUT_DIR}"

# Plot the phylogenetic subtrees for each OTU
mdsine2 render-phylogeny \
    --study ${PREPROCESS_DIR}/gibson_healthy_agg.pkl \
    --output-basepath ${PHYLOGENY_OUT_DIR} \
    --tree files/phylogenetic_placement_OTUs/phylogenetic_tree_full_taxid.nhx \
    --seq-info files/subtrees/RDP-11-5_BA_TS_info.tsv \
    --family-radius-factor 1.5 \
    --top 200

echo "Done."
