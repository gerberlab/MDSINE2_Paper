# MODIFY THESE BASED ON LOCAL ENVIRONMENT.
export DATASET_NAME="synthetic"
source ./settings.sh  # use parent settings (scripts/settings.sh)

_THIS_PATH="${PROJECT_DIR}/scripts/synthetic/settings.sh"
echo "[*] Using additional settings from ${_THIS_PATH}"

export MDSINE2_LOG_INI="${PROJECT_DIR}/scripts/synthetic/logging.ini"

export GLV_PARAMS=${DATASET_DIR}/glv.npz
export TIME_POINTS=${DATASET_DIR}/time_points.txt

export NUM_SAMPLE_TRIALS=10
export COHORT_SIZE=10
export PROCESS_VAR=0.00
export SIMULATION_DT=0.01
export LOW_NOISE_SCALE=0.01
export MEDIUM_NOISE_SCALE=0.15
export HIGH_NOISE_SCALE=0.3

export QPCR_TRIPLICATES=6  # The number of physical+technical replicate combinations to create.

export OUTPUT_DIR="${DATASET_DIR}/output"
export NUM_CORES=1
export MDSINE1_DIR="/home/younhun/mdsine"

export MATLAB_DIR="/mnt/c/Program\ Files/MATLAB/R2015b/bin"
export PATH=${PATH}:${MATLAB_DIR}
export MATLAB="matlab.exe"
