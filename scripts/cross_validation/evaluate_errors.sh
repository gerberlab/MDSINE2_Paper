#!/bin/bash
set -e
source cross_validation/settings.sh


echo "[*] Evaluating cross-validation runs."
python cross_validation/helpers/evaluate_errors.py \
--study ${DATASET_PKL} \
--regression_inputs_dir ${REGRESSION_DATASET_DIR} \
--mdsine_outdir ${OUTPUT_DIR}/mdsine2 \
--clv_elastic_outdir ${OUTPUT_DIR}/regression_clv_elastic-net \
--glv_elastic_outdir ${OUTPUT_DIR}/regression_glv_elastic-net \
--glv_ra_elastic_outdir ${OUTPUT_DIR}/regression_glv-ra_elastic-net \
--glv_ra_ridge_outdir ${OUTPUT_DIR}/regression_glv-ra_ridge \
--glv_ridge_outdir ${OUTPUT_DIR}/regression_glv_ridge \
--lra_elastic_outdir ${OUTPUT_DIR}/regression_lra_elastic-net
