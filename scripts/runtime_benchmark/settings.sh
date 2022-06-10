# MODIFY THESE BASED ON LOCAL ENVIRONMENT.
source ./settings.sh  # use parent settings (scripts/settings.sh)

_THIS_PATH="${PROJECT_DIR}/scripts/runtime_benchmark/settings.sh"
echo "[*] Using additional settings from ${_THIS_PATH}"

export OUTPUT_DIR=${OUTPUT_DIR}/runtime_benchmark
export NUM_TRIALS=5

export MDSINE2_LOG_INI="${PROJECT_DIR}/scripts/runtime_benchmark/logging.ini"
export REGRESSION_CODE_DIR="${PROJECT_DIR}/submodules/clv_fork"
export REGRESSION_DATA_DIR="${REGRESSION_CODE_DIR}/data/gibson"
