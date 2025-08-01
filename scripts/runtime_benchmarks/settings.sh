# MODIFY THESE BASED ON LOCAL ENVIRONMENT.
export DATASET_NAME="gibson/healthy"
source settings.sh  # use parent settings (scripts/settings.sh)

_THIS_PATH="${PROJECT_DIR}/scripts/runtime_benchmark/settings.sh"
echo "[*] Using additional settings from ${_THIS_PATH}"

export MDSINE2_LOG_INI="${PROJECT_DIR}/scripts/analysis/logging.ini"
export FILTER_PICKLE_DIR=${OUTPUT_DIR}/mdsine2_runtime/filtered
export NEGBIN_OUT_DIR=${OUTPUT_DIR}/mdsine2_runtime/negbin
export MDSINE2_OUT_DIR=${OUTPUT_DIR}/mdsine2_runtime/inference

export SOURCE_FILTERED_DSET=${DATASET_DIR}/preprocessed/gibson_healthy_agg_filtered.pkl
export SOURCE_REPLICATE_DSET=${DATASET_DIR}/preprocessed/gibson_replicates_agg_filtered.pkl
#export REPLICATE_MCMC=${NEGBIN_OUT_DIR}/replicates/mcmc.pkl
#export REPLICATE_PLOTS=${NEGBIN_OUT_DIR}/replicates/posterior

export BURNIN=5000
export N_SAMPLES=15000
export CHECKPOINT=100
export MULTIPROCESSING=0
export INTERACTION_IND_PRIOR="strong-sparse"
export PERTURBATION_IND_PRIOR="weak-agnostic"
