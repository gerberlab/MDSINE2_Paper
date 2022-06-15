#!/bin/bash
set -e
source synthetic/settings.sh


export PYTHONPATH=${CLV_DIR}
bash synthetic/helpers/evaluate.py -g ${DATASET_DIR}/glv.npz
