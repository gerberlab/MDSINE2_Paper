# MODIFY THESE BASED ON LOCAL ENVIRONMENT.
source ./settings.sh  # use parent settings (scripts/settings.sh)

_THIS_PATH="${PROJECT_DIR}/scripts/runtime_benchmark/settings.sh"
echo "[*] Using additional settings from ${_THIS_PATH}"

export OUTPUT_DIR=${OUTPUT_DIR}/runtime_benchmark

export MDSINE2_LOG_INI="${PROJECT_DIR}/scripts/runtime_benchmark/logging.ini"
export MDSINE2_NUM_TRIALS=5
export MDSINE2_OUTPUT_DIR=${OUTPUT_DIR}/mdsine2
export MDSINE2_RUNTIMES=${MDSINE2_OUTPUT_DIR}/runtimes.txt

export OTHER_METHODS_OUTPUT_DIR="${OUTPUT_DIR}/other_methods"
export OTHER_METHODS_RUNTIMES="${OTHER_METHODS_OUTPUT_DIR}/other_methods_runtimes.txt"
export OTHER_METHODS_CODE_DIR="${PROJECT_DIR}/scripts/runtime_benchmark/other_methods/code_repo"
export OTHER_METHODS_DATA_DIR="${OTHER_METHODS_CODE_DIR}/data/gibson"


