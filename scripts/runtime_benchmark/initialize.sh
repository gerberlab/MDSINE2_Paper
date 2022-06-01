#!/bin/bash
set -e
source runtime_benchmark/settings.sh


require_program python

cd runtime_benchmark/helpers
study=${DATASET_DIR}/preprocessed/gibson_healthy_agg.pkl

echo "[*] Preparing dataset of 10 taxa..."
python select_topk_bugs.py -f ${study} -n 10 -o ${DATASET_DIR}/trimmed_10/healthy.pkl

echo "[*] Preparing dataset of 25 taxa..."
python select_topk_bugs.py -f ${study} -n 25 -o ${DATASET_DIR}/trimmed_25/healthy.pkl

echo "[*] Preparing dataset of 50 taxa..."
python select_topk_bugs.py -f ${study} -n 50 -o ${DATASET_DIR}/trimmed_50/healthy.pkl

echo "[*] Preparing dataset of 100 taxa..."
python select_topk_bugs.py -f ${study} -n 100 -o ${DATASET_DIR}/trimmed_100/healthy.pkl

echo "[*] Done."
