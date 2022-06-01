# MODIFY THESE BASED ON LOCAL ENVIRONMENT.
export PROJECT_DIR="/mnt/d/Projects/MDSINE2_figures"

# ========= Don't modify below
cd ${PROJECT_DIR}/scripts
source ${PROJECT_DIR}/scripts/runtime_benchmark/settings.sh
cd -

export LSF_GROUP=LSF_JOB_GROUP
export LSF_DIR="./lsf_files"
