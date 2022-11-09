#!/bin/bash
set -e
source synthetic/small/settings.sh


export PYTHONPATH=${CLV_DIR}
echo "[*] Running evaluation script."
python synthetic/helpers/evaluate.py \
-r ${OUTPUT_DIR} \
-o ${OUTPUT_DIR} \
-g ${DATASET_DIR}/glv.npz \
-m 7.09060361e+08 \
-s 4.25878564e+09
