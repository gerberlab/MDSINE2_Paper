#!/bin/bash
set -e
source cross_validation/settings.sh


export PYTHONPATH="${PYTHONPATH}:${CLV_DIR}"
echo "[*] Evaluating cross-validation runs."
python cross_validation/helpers/evaluate_errors.py \
--study ${DATASET_PKL} \
--regression_inputs_dir ${REGRESSION_DATASET_DIR} \
--mdsine_outdir ${OUTPUT_DIR}/mdsine2-modules \
--mdsine_nomodule_outdir ${OUTPUT_DIR}/mdsine2-nomodules \
--mdsine_weak_prior_outdir ${OUTPUT_DIR}/mdsine2-modules-weak-interaction-prior \
--clv_elastic_outdir ${OUTPUT_DIR}/regression_clv_elastic-net \
--glv_elastic_outdir ${OUTPUT_DIR}/regression_glv_elastic-net \
--glv_ra_elastic_outdir ${OUTPUT_DIR}/regression_glv-ra_elastic-net \
--glv_ra_ridge_outdir ${OUTPUT_DIR}/regression_glv-ra_ridge \
--glv_ridge_outdir ${OUTPUT_DIR}/regression_glv_ridge \
--plot_dir ${OUTPUT_DIR} \
--lra_elastic_outdir ${OUTPUT_DIR}/regression_lra_elastic-net \
--subsample_every 1
#--recompute_cache \
