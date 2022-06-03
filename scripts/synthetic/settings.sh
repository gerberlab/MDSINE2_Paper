# MODIFY THESE BASED ON LOCAL ENVIRONMENT.
export DATASET_NAME="synthetic"
source ./settings.sh  # use parent settings (scripts/settings.sh)

_THIS_PATH="${PROJECT_DIR}/scripts/synthetic/settings.sh"
echo "[*] Using additional settings from ${_THIS_PATH}"

export MDSINE2_LOG_INI="${PROJECT_DIR}/scripts/synthetic/logging.ini"

export GLV_PARAMS=${DATASET_DIR}/glv.npz
export TIME_POINTS=${DATASET_DIR}/time_points.txt

export NUM_SEEDS=5
export COHORT_SIZE=5
export PROCESS_VAR=0.01
export SIMULATION_DT=0.01
export NEGBIN_A0=1e-10
export NEGBIN_A1=0.05
export READ_DEPTH=50000
export LOW_NOISE_SCALE=0.01
export MEDIUM_NOISE_SCALE=0.1
export HIGH_NOISE_SCALE=0.2
