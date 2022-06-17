#!/bin/bash
set -e
source synthetic/settings.sh


export PYTHONPATH=${CLV_DIR}
echo "[*] Running evaluation script."
python synthetic/helpers/evaluate.py \
-r ${OUTPUT_DIR} \
-o ${OUTPUT_DIR} \
-g ${DATASET_DIR}/glv.npz \
-m 3723.06368786 \
-s 6567.91725508
