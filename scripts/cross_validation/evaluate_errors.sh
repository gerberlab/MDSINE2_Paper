#!/bin/bash
set -e
#source cross_validation/settings.sh


# MODIFY THESE BASED ON LOCAL ENVIRONMENT.
export DATA_MODALITY="uc"
export DATASET_NAME="gibson/${DATA_MODALITY}"
source ./settings.sh  # use parent settings (scripts/settings.sh)

_THIS_PATH="${PROJECT_DIR}/scripts/cross_validation/settings.sh"
echo "[*] Using additional settings from ${_THIS_PATH}"


export OUTPUT_DIR="${DATASET_DIR}/cross_validation"
export DATASET_PKL=${DATASET_DIR}/preprocessed/gibson_${DATA_MODALITY}_agg_filtered.pkl
export REPLICATE_MCMC=${DATASET_DIR}/output/mdsine2/negbin/replicates/mcmc.pkl

export REGRESSION_DATASET_DIR="${DATASET_DIR}/regression_files"

export MDSINE2_LOG_INI="${PROJECT_DIR}/scripts/cross_validation/logging.ini"
export MDSINE2_BURNIN=5000
export MDSINE2_SAMPLES=15000
export MDSINE2_SAVE_EVERY=1000
export MDSINE2_INTERACTION_INDICATOR_PRIOR="strong-sparse"
export MDSINE2_PERTURBATION_INDICATOR_PRIOR="weak-agnostic"



evaluate_mdsine2_subdir()
{
  out_subdir=$1
  python cross_validation/helpers/mdsine2_errors.py \
    --study "${DATASET_PKL}" \
    --mdsine_outdir "${OUTPUT_DIR}/${out_subdir}" \
    --out_dir "${OUTPUT_DIR}/${out_subdir}" \
    --subsample_every 100 \
    --cohort "${DATA_MODALITY}" \
    --recompute_cache
}


evaluate_regression_subdir()
{
  out_subdir=$1
  model_name=$2

  export PYTHONPATH="${PYTHONPATH}:${CLV_DIR}"
  python cross_validation/helpers/regression_errors.py \
  --study "${DATASET_PKL}" \
  --regression_inputs_dir "${REGRESSION_DATASET_DIR}" \
  --regression_outdir "${OUTPUT_DIR}/${out_subdir}" \
  --out_dir "${OUTPUT_DIR}/${out_subdir}" \
  --model_name "${model_name}"
}


#evaluate_mdsine2_subdir "mdsine2-modules"
#evaluate_mdsine2_subdir "mdsine2-nomodules"
#evaluate_mdsine2_subdir "mdsine2-noqpcr-modules"
#evaluate_mdsine2_subdir "mdsine2-noqpcr-nomodules"
evaluate_mdsine2_subdir "mdsine2-modules-seed10"
#evaluate_regression_subdir "regression_clv_elastic_net" "clv"
#evaluate_regression_subdir "regression_glv_elastic_net" "glv-elastic-net"
#evaluate_regression_subdir "regression_glv_ridge" "glv-ridge"
#evaluate_regression_subdir "regression_glv-ra_elastic_net" "glv-ra-elastic-net"
#evaluate_regression_subdir "regression_glv-ra_ridge" "glv-ra-ridge"
#evaluate_regression_subdir "regression_lra_elastic_net" "lra"


#export PYTHONPATH="${PYTHONPATH}:${CLV_DIR}"
#echo "[*] Evaluating cross-validation runs."
#python cross_validation/helpers/evaluate_errors.py \
#--study ${DATASET_PKL} \
#--regression_inputs_dir ${REGRESSION_DATASET_DIR} \
#--mdsine_outdir ${OUTPUT_DIR}/mdsine2-modules \
#--mdsine_nomodule_outdir ${OUTPUT_DIR}/mdsine2-nomodules \
#--clv_elastic_outdir ${OUTPUT_DIR}/regression_clv_elastic_net \
#--glv_elastic_outdir ${OUTPUT_DIR}/regression_glv_elastic_net \
#--glv_ra_elastic_outdir ${OUTPUT_DIR}/regression_glv-ra_elastic_net \
#--glv_ra_ridge_outdir ${OUTPUT_DIR}/regression_glv-ra_ridge \
#--glv_ridge_outdir ${OUTPUT_DIR}/regression_glv_ridge \
#--out_dir ${OUTPUT_DIR} \
#--lra_elastic_outdir ${OUTPUT_DIR}/regression_lra_elastic_net \
#--subsample_every 100 \
#--healthy_or_uc "${DATA_MODALITY}" \
#--recompute_cache \
