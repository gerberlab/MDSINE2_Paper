#!/bin/bash
set -e
source settings.sh


require_program python
require_program mdsine2
echo "Performing consistency filtering over Study objects."

# Filter the OTUs using consistency filtering

# ==== Healthy
mdsine2 filter \
    --dataset ${PREPROCESS_DIR}/gibson_healthy_agg_taxa.pkl \
    --outfile ${PREPROCESS_DIR}/gibson_healthy_agg_taxa_filtered.pkl \
    --dtype rel \
    --threshold 0.0001 \
    --min-num-consecutive 7 \
    --min-num-subjects 2 \
    --colonization-time 5

python preprocess/helpers/filter_replicates_like_other_dataset.py \
    --replicate-dataset ${PREPROCESS_DIR}/gibson_replicates_agg_taxa.pkl \
    --like-other ${PREPROCESS_DIR}/gibson_healthy_agg_taxa_filtered.pkl \
    --output-basepath ${PREPROCESS_DIR}/gibson_replicates_agg_taxa_filtered.pkl

echo "Done."
