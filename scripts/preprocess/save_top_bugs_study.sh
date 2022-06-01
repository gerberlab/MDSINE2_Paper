#!/bin/bash

python scripts/preprocess/helpers/select_top_bugs.py \
    -f "datasets/gibson/preprocessed/gibson_healthy_agg_taxa_filtered.pkl"\
    -o "output/top_OTUs" \
    -n 10

python scripts/preprocess/helpers/select_top_bugs.py \
    -f "datasets/gibson/preprocessed/gibson_healthy_agg_taxa_filtered.pkl"\
    -o "output/top_OTUs" \
    -n 25

python scripts/preprocess/helpers/select_top_bugs.py \
    -f "datasets/gibson/preprocessed/gibson_healthy_agg_taxa_filtered.pkl"\
    -o "output/top_OTUs" \
    -n 50

  python scripts/preprocess/helpers/select_top_bugs.py \
    -f "datasets/gibson/preprocessed/gibson_healthy_agg_taxa_filtered.pkl"\
    -o "output/top_OTUs" \
    -n 100
