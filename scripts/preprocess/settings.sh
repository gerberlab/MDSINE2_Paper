# MODIFY THESE BASED ON LOCAL ENVIRONMENT.
export DATASET_NAME="gibson"
source ./settings.sh  # use parent settings (scripts/settings.sh)
export MDSINE2_LOG_INI="${PROJECT_DIR}/scripts/preprocess/logging.ini"

export PREPROCESS_DIR="${DATASET_DIR}/preprocessed"
export PHYLOGENY_OUT_DIR="${OUTPUT_DIR}/phylogeny"
export PLOTS_OUT_DIR="${OUTPUT_DIR}/plots"
