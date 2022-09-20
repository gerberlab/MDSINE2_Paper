# MODIFY THESE BASED ON LOCAL ENVIRONMENT.
export DATASET_NAME="gibson"
source ./settings.sh  # use parent settings (scripts/settings.sh)

_THIS_PATH="${PROJECT_DIR}/scripts/runtime_benchmark/settings.sh"
echo "[*] Using additional settings from ${_THIS_PATH}"

export MDSINE2_LOG_INI="${PROJECT_DIR}/scripts/analysis/logging.ini"
export NEGBIN_OUT_DIR=${OUTPUT_DIR}/mdsine2/negbin
export MDSINE2_OUT_DIR=${OUTPUT_DIR}/mdsine2/inference

export HEALTHY_DSET=${DATASET_DIR}/preprocessed/gibson_healthy_agg_filtered.pkl
export REPLICATE_DSET=${DATASET_DIR}/preprocessed/gibson_replicates_agg_filtered.pkl
export REPLICATE_MCMC=${NEGBIN_OUT_DIR}/replicates/mcmc.pkl
export REPLICATE_PLOTS=${NEGBIN_OUT_DIR}/replicates/posterior

export BURNIN=50
export N_SAMPLES=150
export CHECKPOINT=10
export MULTIPROCESSING=0
export INTERACTION_IND_PRIOR="strong-sparse"
export PERTURBATION_IND_PRIOR="weak-agnostic"
