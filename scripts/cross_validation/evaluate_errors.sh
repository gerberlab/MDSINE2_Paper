#!/bin/bash
set -e
source cross_validation/settings.sh


echo "[*] Evaluating cross-validation runs."
python cross_validation/helpers/evaluate_errors.py \
--study ${} \
--mdsine_outdir ${OUTPUT_DIR}/mdsine2 \
--clv_elastic_outdir ${} \
--glv_elastic_outdir ${} \
--glv_ra_elastic_outdir ${} \
--glv_ra_ridge_outdir ${} \
--glv_ridge_outdir ${} \
--lra_elastic_outdir ${}