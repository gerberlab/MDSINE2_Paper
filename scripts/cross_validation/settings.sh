# MODIFY THESE BASED ON LOCAL ENVIRONMENT.
export DATASET_NAME="gibson"
source ./settings.sh  # use parent settings (scripts/settings.sh)

_THIS_PATH="${PROJECT_DIR}/scripts/cross_validation/settings.sh"
echo "[*] Using additional settings from ${_THIS_PATH}"


export OUTPUT_DIR="${DATASET_DIR}/cross_validation"
export DATASET_PKL=${DATASET_DIR}/preprocessed/gibson_healthy_agg_filtered.pkl
export REPLICATE_MCMC=${DATASET_DIR}/output/mdsine2/negbin/replicates/mcmc.pkl

export REGRESSION_DATASET_DIR="${DATASET_DIR}/regression_files"

export MDSINE2_LOG_INI="${PROJECT_DIR}/scripts/cross_validation/logging.ini"
export MDSINE2_BURNIN=5000
export MDSINE2_SAMPLES=15000
export MDSINE2_SAVE_EVERY=1000
export MDSINE2_INTERACTION_INDICATOR_PRIOR="strong-sparse"
export MDSINE2_PERTURBATION_INDICATOR_PRIOR="weak-agnostic"
