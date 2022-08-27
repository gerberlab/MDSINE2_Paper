# MODIFY THESE BASED ON LOCAL ENVIRONMENT.
export DATASET_NAME="gibson"
source ./settings.sh  # use parent settings (scripts/settings.sh)

_THIS_PATH="${PROJECT_DIR}/scripts/synthetic/settings.sh"
echo "[*] Using additional settings from ${_THIS_PATH}"

export PREPROCESS_DIR="${DATASET_DIR}/preprocessed"
export INPUT_DATASET_DIR="${DATASET_DIR}/cv_comparator_methods"
export CV_OUTPUT_DIR="${INPUT_DATASET_DIR}/cv_inference_output"
export FW_SIM_OUTPUT_DIR="${INPUT_DATASET_DIR}/forward_sim_lv"
