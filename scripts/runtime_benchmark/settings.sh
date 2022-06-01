# MODIFY THESE BASED ON LOCAL ENVIRONMENT.
source ./settings.sh  # use parent settings (scripts/settings.sh)

_THIS_PATH="${PROJECT_DIR}/scripts/runtime_benchmark/settings.sh"
echo "[*] Using additional settings from ${_THIS_PATH}"

export OUTPUT_DIR=${OUTPUT_DIR}/runtime_benchmark

export MDSINE2_LOG_INI="${PROJECT_DIR}/scripts/runtime_benchmark/logging.ini"
export MDSINE2_NUM_TRIALS=5
export MDSINE2_OUTPUT_DIR=${OUTPUT_DIR}/mdsine2
export MDSINE2_RUNTIMES=${MDSINE2_OUTPUT_DIR}/runtimes.txt
