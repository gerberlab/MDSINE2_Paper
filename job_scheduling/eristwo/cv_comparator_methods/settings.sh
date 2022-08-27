# MODIFY THESE BASED ON LOCAL ENVIRONMENT.
_cwd=$(pwd)
cd ../../..
export PROJECT_DIR="$(pwd)"
cd -

_this_path="${_cwd}/settings.sh"
echo "[*] Using settings from ${_this_path}"

# ========= Don't modify below
export LSF_GROUP=LSF_JOB_GROUP
export LSF_DIR="${PROJECT_DIR}/job_scheduling/eristwo/cv_comparator_methods/lsf_files"
export LSF_QUEUE="normal"
