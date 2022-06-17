#!/bin/bash
set -e
source synthetic/settings.sh


echo "[*] Running plotting script."
python synthetic/helpers/render_evaluation.py \
-o ${OUTPUT_DIR} \
-f pdf
